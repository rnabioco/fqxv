//! Clean-room LZMA-class codec for the global reference bases.
//!
//! The order-k context model ([`fqxv_seq`], the shipped reference coder) and the
//! block-sorting coder ([`super::refbwt`]) both leave the reference's *long-range*
//! redundancy — near-duplicate contigs scattered across the whole reference — on
//! the table, because neither *copies* a repeat: the context model re-codes it and
//! BWT only models it statistically. The way to capture it is explicit LZ77
//! matching with strong entropy coding, which is what xz/LZMA do (they reach
//! ~1.79 b/base here vs the order-k model's 1.98). This module is a clean-room
//! LZMA-class coder built only on the [`fqxv_range`] range coder — it is *not* a
//! port of liblzma/xz; the bit-models, the match finder, the parse, and the
//! stream layout are all local, following the well-known LZMA *design* (adaptive
//! bit-models, a 12-state machine, context literals with matched-byte prediction,
//! a length coder, position-slot + aligned distance coding, and rep0–3 short
//! codes for repeated distances).
//!
//! ## Window and determinism
//!
//! Unlike the per-read blocks (which must decode in parallel), the reference is
//! coded and decoded exactly once per file, so this uses a single whole-reference
//! window (a serial LZ decode) to give the match finder the full ~89 MB of reach
//! — the whole point, since the repeats are far apart. Encoding is deterministic
//! (a fixed greedy/lazy parse over a fixed hash chain), so the output is identical
//! regardless of thread count. The container gates it never-worse against the
//! order-k path, so it can only shrink the archive.

use fqxv_bytes::{read_varint, write_varint};
use fqxv_range::{Decoder, Encoder};

use crate::{Error, Result};

// --- adaptive binary model over the range coder ------------------------------

/// Probability scale for the adaptive bit-models (LZMA uses 11 bits). Must stay
/// below the range coder's `BOT` (2^16) so `range / total` is exact.
const PROB_BITS: u32 = 11;
const PROB_TOTAL: u32 = 1 << PROB_BITS;
const PROB_INIT: u16 = (PROB_TOTAL / 2) as u16;
/// Adaptation rate (LZMA `kNumMoveBits`). The shift-toward update never drives the
/// probability to 0 or `PROB_TOTAL`, so both intervals stay non-empty.
const MOVE_BITS: u32 = 5;

/// One adaptive bit (probability that the next bit is 0).
#[derive(Clone)]
struct Bit(u16);

impl Bit {
    #[inline]
    fn new() -> Self {
        Bit(PROB_INIT)
    }

    #[inline]
    fn encode(&mut self, enc: &mut Encoder, bit: u32) {
        let p0 = u32::from(self.0);
        if bit == 0 {
            enc.encode(0, p0, PROB_TOTAL);
            self.0 += ((PROB_TOTAL - p0) >> MOVE_BITS) as u16;
        } else {
            enc.encode(p0, PROB_TOTAL - p0, PROB_TOTAL);
            self.0 -= (p0 >> MOVE_BITS) as u16;
        }
    }

    #[inline]
    fn decode(&mut self, dec: &mut Decoder<'_>) -> u32 {
        let p0 = u32::from(self.0);
        let t = dec.freq(PROB_TOTAL);
        if t < p0 {
            dec.decode(0, p0);
            self.0 += ((PROB_TOTAL - p0) >> MOVE_BITS) as u16;
            0
        } else {
            dec.decode(p0, PROB_TOTAL - p0);
            self.0 -= (p0 >> MOVE_BITS) as u16;
            1
        }
    }
}

/// Encode `num_bits` equiprobable ("direct") bits, MSB first.
#[inline]
fn encode_direct(enc: &mut Encoder, value: u32, num_bits: u32) {
    for i in (0..num_bits).rev() {
        enc.encode((value >> i) & 1, 1, 2);
    }
}

#[inline]
fn decode_direct(dec: &mut Decoder<'_>, num_bits: u32) -> u32 {
    let mut v = 0u32;
    for _ in 0..num_bits {
        let b = dec.freq(2).min(1);
        dec.decode(b, 1);
        v = (v << 1) | b;
    }
    v
}

/// A bit-tree of `2^n - 1` adaptive bits coding an `n`-bit value MSB-first — the
/// LZMA primitive for length slots, position slots, and literals.
#[derive(Clone)]
struct Tree {
    bits: Vec<Bit>,
    n: u32,
}

impl Tree {
    fn new(n: u32) -> Self {
        Tree {
            bits: vec![Bit::new(); 1usize << n],
            n,
        }
    }

    #[inline]
    fn encode(&mut self, enc: &mut Encoder, value: u32) {
        let mut m = 1usize;
        for i in (0..self.n).rev() {
            let bit = (value >> i) & 1;
            self.bits[m].encode(enc, bit);
            m = (m << 1) | bit as usize;
        }
    }

    #[inline]
    fn decode(&mut self, dec: &mut Decoder<'_>) -> u32 {
        let mut m = 1usize;
        for _ in 0..self.n {
            let bit = self.bits[m].decode(dec);
            m = (m << 1) | bit as usize;
        }
        (m - (1 << self.n)) as u32
    }

    /// Reverse (LSB-first) bit-tree, used for the aligned low bits of a distance.
    #[inline]
    fn encode_reverse(&mut self, enc: &mut Encoder, mut value: u32) {
        let mut m = 1usize;
        for _ in 0..self.n {
            let bit = value & 1;
            value >>= 1;
            self.bits[m].encode(enc, bit);
            m = (m << 1) | bit as usize;
        }
    }

    #[inline]
    fn decode_reverse(&mut self, dec: &mut Decoder<'_>) -> u32 {
        let mut m = 1usize;
        let mut value = 0u32;
        for i in 0..self.n {
            let bit = self.bits[m].decode(dec);
            m = (m << 1) | bit as usize;
            value |= bit << i;
        }
        value
    }
}

// --- LZMA structural constants -----------------------------------------------

const NUM_STATES: usize = 12;
const MIN_MATCH: usize = 2;
/// Length coder split points (LZMA): 8 low, 8 mid, 256 high.
const LEN_LOW: u32 = 8;
const LEN_MID: u32 = 8;
const LEN_HIGH: u32 = 256;
const MAX_MATCH: usize = MIN_MATCH + (LEN_LOW + LEN_MID + LEN_HIGH) as usize - 1; // 273

/// Distance-slot coding constants (LZMA).
const NUM_LEN_TO_POS: usize = 4; // length states that select a pos-slot tree
const NUM_POS_SLOT_BITS: u32 = 6;
const END_POS_MODEL: u32 = 14; // slots < this use SpecPos reverse trees
const ALIGN_BITS: u32 = 4;

// --- length coder ------------------------------------------------------------

/// LZMA length coder: a `choice`/`choice2` gate selecting a low (3-bit), mid
/// (3-bit), or high (8-bit) tree. Codes a match length in `MIN_MATCH..=MAX_MATCH`.
#[derive(Clone)]
struct LenCoder {
    choice: Bit,
    choice2: Bit,
    low: Tree,
    mid: Tree,
    high: Tree,
}

impl LenCoder {
    fn new() -> Self {
        LenCoder {
            choice: Bit::new(),
            choice2: Bit::new(),
            low: Tree::new(3),
            mid: Tree::new(3),
            high: Tree::new(8),
        }
    }

    #[inline]
    fn encode(&mut self, enc: &mut Encoder, len: usize) {
        let mut l = (len - MIN_MATCH) as u32;
        if l < LEN_LOW {
            self.choice.encode(enc, 0);
            self.low.encode(enc, l);
        } else {
            self.choice.encode(enc, 1);
            l -= LEN_LOW;
            if l < LEN_MID {
                self.choice2.encode(enc, 0);
                self.mid.encode(enc, l);
            } else {
                self.choice2.encode(enc, 1);
                self.high.encode(enc, l - LEN_MID);
            }
        }
    }

    #[inline]
    fn decode(&mut self, dec: &mut Decoder<'_>) -> usize {
        let l = if self.choice.decode(dec) == 0 {
            self.low.decode(dec)
        } else if self.choice2.decode(dec) == 0 {
            LEN_LOW + self.mid.decode(dec)
        } else {
            LEN_LOW + LEN_MID + self.high.decode(dec)
        };
        MIN_MATCH + l as usize
    }
}

// --- the full model set ------------------------------------------------------

/// Number of `Bit` probabilities per `lc`-context literal model (LZMA `0x300`):
/// `0x100` for the non-matched region plus two `0x100` regions for matched-bit
/// prediction.
const LIT_SIZE: usize = 0x300;

/// All adaptive models for one coding stream. `lc` literal-context bits select a
/// literal sub-model by the high bits of the previous byte.
struct Models {
    is_match: [Bit; NUM_STATES],
    is_rep: [Bit; NUM_STATES],
    is_rep_g0: [Bit; NUM_STATES],
    is_rep0_long: [Bit; NUM_STATES],
    is_rep_g1: [Bit; NUM_STATES],
    is_rep_g2: [Bit; NUM_STATES],
    pos_slot: [Tree; NUM_LEN_TO_POS],
    /// One reverse bit-tree per pos-slot in `4..END_POS_MODEL` (indexed by
    /// `slot - 4`); each has `(slot >> 1) - 1` direct bits.
    spec: Vec<Tree>,
    align: Tree,
    len: LenCoder,
    rep_len: LenCoder,
    /// Literal models, one `LIT_SIZE`-probability set per `lc`-bit context.
    lit: Vec<Vec<Bit>>,
    lc: u32,
}

impl Models {
    fn new(lc: u32) -> Self {
        let spec = (4..END_POS_MODEL)
            .map(|slot| Tree::new((slot >> 1) - 1))
            .collect();
        Models {
            is_match: std::array::from_fn(|_| Bit::new()),
            is_rep: std::array::from_fn(|_| Bit::new()),
            is_rep_g0: std::array::from_fn(|_| Bit::new()),
            is_rep0_long: std::array::from_fn(|_| Bit::new()),
            is_rep_g1: std::array::from_fn(|_| Bit::new()),
            is_rep_g2: std::array::from_fn(|_| Bit::new()),
            pos_slot: std::array::from_fn(|_| Tree::new(NUM_POS_SLOT_BITS)),
            spec,
            align: Tree::new(ALIGN_BITS),
            len: LenCoder::new(),
            rep_len: LenCoder::new(),
            lit: vec![vec![Bit::new(); LIT_SIZE]; 1usize << lc],
            lc,
        }
    }
}

/// LZMA state transitions.
#[inline]
fn state_after_lit(s: usize) -> usize {
    if s < 4 {
        0
    } else if s < 10 {
        s - 3
    } else {
        s - 6
    }
}
#[inline]
fn state_after_match(s: usize) -> usize {
    if s < 7 {
        7
    } else {
        10
    }
}
#[inline]
fn state_after_rep(s: usize) -> usize {
    if s < 7 {
        8
    } else {
        11
    }
}
#[inline]
fn state_after_shortrep(s: usize) -> usize {
    if s < 7 {
        9
    } else {
        11
    }
}

/// LZMA pos-slot for a 0-based distance.
#[inline]
fn pos_slot(dist: u32) -> u32 {
    if dist < 4 {
        dist
    } else {
        let n = 31 - dist.leading_zeros();
        (n << 1) | ((dist >> (n - 1)) & 1)
    }
}

#[inline]
fn len_to_pos_state(len: usize) -> usize {
    (len - MIN_MATCH).min(NUM_LEN_TO_POS - 1)
}

// --- literal coding ----------------------------------------------------------

#[inline]
fn lit_state(lc: u32, prev: u8) -> usize {
    (u32::from(prev) >> (8 - lc)) as usize
}

/// Encode a literal byte into a `LIT_SIZE` probability set. When the previous op
/// was a match, the byte is coded against `match_byte` (the byte at distance
/// `rep0 + 1`) until they diverge — LZMA's matched-literal prediction.
#[inline]
fn encode_literal(enc: &mut Encoder, probs: &mut [Bit], byte: u8, match_byte: u8, matched: bool) {
    let sym = u32::from(byte);
    let mut ctx = 1u32;
    let mut still = matched;
    for i in (0..8).rev() {
        let bit = (sym >> i) & 1;
        if still {
            let mbit = (u32::from(match_byte) >> i) & 1;
            let idx = (((1 + mbit) << 8) | ctx) as usize;
            probs[idx].encode(enc, bit);
            if mbit != bit {
                still = false;
            }
        } else {
            probs[ctx as usize].encode(enc, bit);
        }
        ctx = (ctx << 1) | bit;
    }
}

#[inline]
fn decode_literal(dec: &mut Decoder<'_>, probs: &mut [Bit], match_byte: u8, matched: bool) -> u8 {
    let mut ctx = 1u32;
    let mut still = matched;
    for i in (0..8).rev() {
        let bit = if still {
            let mbit = (u32::from(match_byte) >> i) & 1;
            let idx = (((1 + mbit) << 8) | ctx) as usize;
            let b = probs[idx].decode(dec);
            if mbit != b {
                still = false;
            }
            b
        } else {
            probs[ctx as usize].decode(dec)
        };
        ctx = (ctx << 1) | bit;
    }
    (ctx & 0xff) as u8
}

// --- distance coding ---------------------------------------------------------

/// Encode a 0-based distance `dist` for a match of length `len` (LZMA pos-slot +
/// SpecPos/align coding).
#[inline]
fn encode_distance(m: &mut Models, enc: &mut Encoder, dist: u32, len: usize) {
    let len_state = len_to_pos_state(len);
    let slot = pos_slot(dist);
    m.pos_slot[len_state].encode(enc, slot);
    if slot >= 4 {
        let num_direct = (slot >> 1) - 1;
        let base = (2 | (slot & 1)) << num_direct;
        let rest = dist - base;
        if slot < END_POS_MODEL {
            m.spec[(slot - 4) as usize].encode_reverse(enc, rest);
        } else {
            encode_direct(enc, rest >> ALIGN_BITS, num_direct - ALIGN_BITS);
            m.align.encode_reverse(enc, rest & ((1 << ALIGN_BITS) - 1));
        }
    }
}

#[inline]
fn decode_distance(m: &mut Models, dec: &mut Decoder<'_>, len: usize) -> Result<u32> {
    let len_state = len_to_pos_state(len);
    let slot = m.pos_slot[len_state].decode(dec);
    if slot < 4 {
        return Ok(slot);
    }
    let num_direct = (slot >> 1) - 1;
    let base = (2 | (slot & 1)) << num_direct;
    let rest = if slot < END_POS_MODEL {
        let idx = (slot - 4) as usize;
        if idx >= m.spec.len() {
            return Err(Error::Malformed("reflzma: bad pos-slot"));
        }
        m.spec[idx].decode_reverse(dec)
    } else {
        let hi = decode_direct(dec, num_direct - ALIGN_BITS);
        let lo = m.align.decode_reverse(dec);
        (hi << ALIGN_BITS) | lo
    };
    Ok(base + rest)
}

// --- match finder (hash-chain) -----------------------------------------------

const NIL: u32 = u32::MAX;
/// Bytes hashed to seed the chain; see [`hash_at`]. A 4-byte seed keys a 16-base
/// match on the 2-bit-packed reference and finds more (shorter) matches than a
/// long seed.
const HASH_LEN: usize = 4;
const HASH_BITS: u32 = 22;
/// Stop searching once a match at least this long is found.
const NICE_LEN: usize = 273;
/// Max hash-chain nodes to visit per position.
const MAX_CHAIN: usize = 1024;

#[inline]
fn hash_at(s: &[u8], pos: usize) -> usize {
    // Hash `HASH_LEN` bytes -> bucket. Any deterministic mix works (affects only
    // which matches are found, never correctness). A short seed finds more (and
    // shorter) matches; on the 2-bit-packed reference each byte is 4 bases, so a
    // 4-byte seed already keys a 16-base match.
    let v = u32::from_le_bytes(s[pos..pos + 4].try_into().unwrap());
    (u64::from(v).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - HASH_BITS)) as usize
}

/// Longest match for the suffix at `pos` among earlier positions on its chain.
/// Returns `(len, dist)` with `dist` 1-based (`pos - candidate`), or `(0, 0)`.
#[inline]
fn find_match(s: &[u8], pos: usize, head: &[u32], prev: &[u32], mask: usize) -> (usize, usize) {
    let n = s.len();
    if pos + HASH_LEN > n {
        return (0, 0);
    }
    let max_len = (n - pos).min(MAX_MATCH);
    let h = hash_at(s, pos) & mask;
    let mut cand = head[h];
    let mut best_len = 0usize;
    let mut best_dist = 0usize;
    let mut chain = 0usize;
    while cand != NIL && chain < MAX_CHAIN {
        let c = cand as usize;
        if best_len == 0 || s[c + best_len] == s[pos + best_len] {
            let mut l = 0usize;
            while l < max_len && s[c + l] == s[pos + l] {
                l += 1;
            }
            if l > best_len {
                best_len = l;
                best_dist = pos - c;
                // Stop when maximal — both to avoid wasted work and, critically, so
                // the quick-reject `s[pos + best_len]` below never reads past `n`.
                if best_len >= max_len || best_len >= NICE_LEN {
                    break;
                }
            }
        }
        cand = prev[c];
        chain += 1;
    }
    (best_len, best_dist)
}

#[inline]
fn insert(s: &[u8], pos: usize, head: &mut [u32], prev: &mut [u32], mask: usize) {
    if pos + HASH_LEN <= s.len() {
        let h = hash_at(s, pos) & mask;
        prev[pos] = head[h];
        head[h] = pos as u32;
    }
}

/// Match length of `s[pos..]` against the earlier copy at 0-based distance `d`
/// (i.e. starting `d + 1` bytes back), capped at `MAX_MATCH`.
#[inline]
fn rep_len_at(s: &[u8], pos: usize, d: u32) -> usize {
    let back = d as usize + 1;
    if back > pos {
        return 0;
    }
    let start = pos - back;
    let n = s.len();
    let max_len = (n - pos).min(MAX_MATCH);
    let mut l = 0usize;
    while l < max_len && s[start + l] == s[pos + l] {
        l += 1;
    }
    l
}

// --- main encode / decode ----------------------------------------------------

/// Literal-context bits: full previous byte (DNA's tiny alphabet makes this
/// cheap and gives order-1-byte literal context).
const LC: u32 = 8;
/// Shortest normal (non-rep) match worth coding at all; shorter repeats are only
/// taken as reps (cheap distance). Distance-scaled below.
const NORMAL_MIN: usize = 4;

/// LZMA-encode the whole byte slice `s` into one range-coded stream. Exposed to
/// sibling modules (e.g. [`super::refpack`]) that want the raw byte-domain LZMA
/// core without the reference `(lens, seq)` framing.
pub(crate) fn lzma_encode(s: &[u8]) -> Vec<u8> {
    let n = s.len();
    let mut m = Models::new(LC);
    let mut enc = Encoder::new();
    let mut state = 0usize;
    let mut reps = [0u32; 4];

    let mask = (1usize << HASH_BITS) - 1;
    let mut head = vec![NIL; 1usize << HASH_BITS];
    let mut prev = vec![NIL; n.max(1)];

    let mut pos = 0usize;
    while pos < n {
        let prev_byte = if pos > 0 { s[pos - 1] } else { 0 };
        // Candidate normal match at `pos` (searches earlier positions only), then
        // insert `pos` so the lazy peek at `pos + 1` can see it.
        let (mut mlen, mdist) = find_match(s, pos, &head, &prev, mask);
        insert(s, pos, &mut head, &mut prev, mask);
        let mut best_rep_len = 0usize;
        let mut best_rep_idx = 0usize;
        for (i, &d) in reps.iter().enumerate() {
            let l = rep_len_at(s, pos, d);
            if l > best_rep_len {
                best_rep_len = l;
                best_rep_idx = i;
            }
        }
        // Distance-scaled acceptance for normal matches (a far match must be
        // longer to beat coding it as literals).
        let normal_min = NORMAL_MIN
            + if mdist >= (1 << 16) {
                2
            } else if mdist >= 512 {
                1
            } else {
                0
            };
        if mlen < normal_min {
            mlen = 0;
        }

        // Decide op: prefer a rep when it is competitive (its distance is nearly
        // free), else a good normal match, else a literal.
        let take_rep = best_rep_len >= MIN_MATCH && best_rep_len + 1 >= mlen;

        // Lazy matching: if we'd take a NORMAL match but an even longer match
        // starts at `pos + 1`, emit a literal now and let `pos + 1` take the
        // better match. Classic ~10-15% ratio win over pure greedy. (Reps are
        // near-free, so never defer them.)
        if !take_rep && mlen >= MIN_MATCH && pos + 1 < n {
            let (nlen, _) = find_match(s, pos + 1, &head, &prev, mask);
            if nlen > mlen {
                m.is_match[state].encode(&mut enc, 0);
                let match_byte = {
                    let back = reps[0] as usize + 1;
                    if back <= pos {
                        s[pos - back]
                    } else {
                        0
                    }
                };
                let ls = lit_state(m.lc, prev_byte);
                let matched = state >= 7;
                encode_literal(&mut enc, &mut m.lit[ls], s[pos], match_byte, matched);
                state = state_after_lit(state);
                pos += 1; // `pos` already inserted above
                continue;
            }
        }

        if take_rep {
            let len = best_rep_len;
            m.is_match[state].encode(&mut enc, 1);
            m.is_rep[state].encode(&mut enc, 1);
            // Rep-index selector; `is_rep0_long` is coded ONLY for rep0 (matching
            // the decoder, which reads it only when idx == 0).
            if best_rep_idx == 0 {
                m.is_rep_g0[state].encode(&mut enc, 0);
                m.is_rep0_long[state].encode(&mut enc, 1); // len >= 2, not a short rep
            } else {
                m.is_rep_g0[state].encode(&mut enc, 1);
                if best_rep_idx == 1 {
                    m.is_rep_g1[state].encode(&mut enc, 0);
                } else {
                    m.is_rep_g1[state].encode(&mut enc, 1);
                    m.is_rep_g2[state].encode(&mut enc, (best_rep_idx - 2) as u32);
                }
            }
            m.rep_len.encode(&mut enc, len);
            // Move the used rep to front.
            let d = reps[best_rep_idx];
            for j in (1..=best_rep_idx).rev() {
                reps[j] = reps[j - 1];
            }
            reps[0] = d;
            state = state_after_rep(state);
            // `pos` already inserted; insert the rest of the covered positions.
            for k in 1..len {
                insert(s, pos + k, &mut head, &mut prev, mask);
            }
            pos += len;
        } else if mlen >= MIN_MATCH {
            let len = mlen;
            let dist0 = (mdist - 1) as u32;
            m.is_match[state].encode(&mut enc, 1);
            m.is_rep[state].encode(&mut enc, 0);
            m.len.encode(&mut enc, len);
            encode_distance(&mut m, &mut enc, dist0, len);
            reps = [dist0, reps[0], reps[1], reps[2]];
            state = state_after_match(state);
            for k in 1..len {
                insert(s, pos + k, &mut head, &mut prev, mask);
            }
            pos += len;
        } else {
            // Literal (`pos` already inserted at the top of the loop).
            m.is_match[state].encode(&mut enc, 0);
            let match_byte = {
                let back = reps[0] as usize + 1;
                if back <= pos {
                    s[pos - back]
                } else {
                    0
                }
            };
            let ls = lit_state(m.lc, prev_byte);
            let matched = state >= 7;
            encode_literal(&mut enc, &mut m.lit[ls], s[pos], match_byte, matched);
            state = state_after_lit(state);
            pos += 1;
        }
    }
    enc.finish()
}

/// LZMA-decode `n` bytes from `buf`. Exposed to sibling modules.
pub(crate) fn lzma_decode(buf: &[u8], n: usize) -> Result<Vec<u8>> {
    let mut m = Models::new(LC);
    let mut dec = Decoder::new(buf);
    let mut state = 0usize;
    let mut reps = [0u32; 4];
    // `n` is an untrusted decoded-length header; reserve fallibly so a corrupt
    // value errors rather than aborting on a huge infallible allocation.
    let mut out: Vec<u8> = Vec::new();
    out.try_reserve(n)
        .map_err(|_| Error::Malformed("reflzma: output too large to allocate"))?;

    while out.len() < n {
        let pos = out.len();
        let prev_byte = if pos > 0 { out[pos - 1] } else { 0 };
        if m.is_match[state].decode(&mut dec) == 0 {
            // Literal.
            let match_byte = {
                let back = reps[0] as usize + 1;
                if back <= pos {
                    out[pos - back]
                } else {
                    0
                }
            };
            let ls = lit_state(m.lc, prev_byte);
            let matched = state >= 7;
            let b = decode_literal(&mut dec, &mut m.lit[ls], match_byte, matched);
            out.push(b);
            state = state_after_lit(state);
            continue;
        }
        // Match or rep.
        let (len, dist0) = if m.is_rep[state].decode(&mut dec) == 1 {
            // Rep match: select which rep.
            let idx = decode_rep_index(&mut m, &mut dec, state);
            if idx == 0 && m.is_rep0_long[state].decode(&mut dec) == 0 {
                // Short rep: length 1 at rep0.
                let d = reps[0];
                state = state_after_shortrep(state);
                (1usize, d)
            } else {
                let d = reps[idx];
                // Move used rep to front.
                for j in (1..=idx).rev() {
                    reps[j] = reps[j - 1];
                }
                reps[0] = d;
                let len = m.rep_len.decode(&mut dec);
                state = state_after_rep(state);
                (len, d)
            }
        } else {
            // Normal match.
            let len = m.len.decode(&mut dec);
            let d = decode_distance(&mut m, &mut dec, len)?;
            reps = [d, reps[0], reps[1], reps[2]];
            state = state_after_match(state);
            (len, d)
        };
        // Copy the match from the already-decoded output.
        let back = dist0 as usize + 1;
        if back > out.len() {
            return Err(Error::Malformed("reflzma: distance past output"));
        }
        if out.len() + len > n {
            return Err(Error::Malformed("reflzma: match past end"));
        }
        let start = out.len() - back;
        for i in 0..len {
            let b = out[start + i];
            out.push(b);
        }
    }
    Ok(out)
}

#[inline]
fn decode_rep_index(m: &mut Models, dec: &mut Decoder<'_>, state: usize) -> usize {
    if m.is_rep_g0[state].decode(dec) == 0 {
        0
    } else if m.is_rep_g1[state].decode(dec) == 0 {
        1
    } else if m.is_rep_g2[state].decode(dec) == 0 {
        2
    } else {
        3
    }
}

// --- whole-reference framing -------------------------------------------------

/// LZMA-code a whole reference given per-contig `lens` and concatenated `seq`.
/// One window over the whole reference (the reference is coded once, so a serial
/// decode is fine and the full reach is what captures the distant repeats).
/// Frame: `[varint n_contigs][lens...][varint n_bases][lzma payload]`.
pub(crate) fn encode(lens: &[u32], seq: &[u8]) -> Result<Vec<u8>> {
    let sum: usize = lens.iter().map(|&l| l as usize).sum();
    if sum != seq.len() {
        return Err(Error::Malformed("reflzma: length/seq disagreement"));
    }
    let payload = lzma_encode(seq);
    let mut out = Vec::with_capacity(payload.len() + lens.len() * 2 + 16);
    write_varint(&mut out, lens.len() as u64);
    for &l in lens {
        write_varint(&mut out, u64::from(l));
    }
    write_varint(&mut out, seq.len() as u64);
    out.extend_from_slice(&payload);
    Ok(out)
}

/// Inverse of [`encode`]: returns `(lens, seq)`.
pub(crate) fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut p = 0usize;
    let nc =
        read_varint(src, &mut p).ok_or(Error::Malformed("reflzma: bad contig count"))? as usize;
    if nc > src.len() + 1 {
        return Err(Error::Malformed("reflzma: implausible contig count"));
    }
    let mut lens = Vec::with_capacity(nc.min(1 << 20));
    for _ in 0..nc {
        lens.push(
            read_varint(src, &mut p).ok_or(Error::Malformed("reflzma: bad contig len"))? as u32,
        );
    }
    let nb = read_varint(src, &mut p).ok_or(Error::Malformed("reflzma: bad n_bases"))? as usize;
    // Bomb guard: the range coder emits at least ~1 bit per op and each op yields
    // >=1 byte, so a payload cannot legitimately expand beyond this factor.
    let max_plausible = src.len().saturating_mul(1 << 18).saturating_add(1 << 18);
    if nb > max_plausible {
        return Err(Error::Malformed(
            "reflzma: declared length exceeds capacity",
        ));
    }
    let seq = lzma_decode(&src[p..], nb)?;
    let sum: usize = lens.iter().map(|&l| l as usize).sum();
    if sum != seq.len() {
        return Err(Error::Malformed("reflzma: block length/seq disagreement"));
    }
    Ok((lens, seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(lens: &[u32], seq: &[u8]) {
        let coded = encode(lens, seq).expect("encode");
        let (out_lens, out_seq) = decode(&coded).expect("decode");
        assert_eq!(out_lens, lens, "lens mismatch");
        assert_eq!(out_seq, seq, "seq mismatch");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[], b"");
    }

    #[test]
    fn roundtrip_short() {
        roundtrip(&[4], b"ACGT");
        roundtrip(&[12], b"ACGTNNRYACGT");
    }

    #[test]
    fn roundtrip_repeats_compress() {
        let unit: &[u8] = b"ACGTACGGTATTCCGGAACCTTGGACGTACGGTATTCCGGAACCTTGG";
        let mut seq = Vec::new();
        let mut lens = Vec::new();
        for _ in 0..400 {
            seq.extend_from_slice(unit);
            lens.push(unit.len() as u32);
        }
        roundtrip(&lens, &seq);
        let coded = encode(&lens, &seq).unwrap();
        assert!(
            coded.len() < seq.len() / 20,
            "LZMA should crush exact repeats: {} -> {}",
            seq.len(),
            coded.len()
        );
    }

    #[test]
    fn decode_rejects_garbage_without_panic() {
        for seed in 0u8..64 {
            let junk: Vec<u8> = (0..300u16)
                .map(|i| (i as u8).wrapping_mul(31) ^ seed)
                .collect();
            let _ = decode(&junk);
        }
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(b"ACGTNacgtRYKM".to_vec()), 0..90),
                0..70),
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            let coded = encode(&lens, &seq).expect("encode");
            let (out_lens, out_seq) = decode(&coded).expect("decode");
            proptest::prop_assert_eq!(out_lens, lens);
            proptest::prop_assert_eq!(out_seq, seq);
        }

        #[test]
        fn roundtrip_arbitrary_bytes(seq in proptest::collection::vec(0u8..=255, 0..400)) {
            let coded = encode(&[seq.len() as u32], &seq).expect("encode");
            let (_lens, out_seq) = decode(&coded).expect("decode");
            proptest::prop_assert_eq!(out_seq, seq);
        }

        #[test]
        fn roundtrip_repetitive(
            unit in proptest::collection::vec(proptest::sample::select(b"ACGT".to_vec()), 4..40),
            reps in 1usize..300,
        ) {
            let seq: Vec<u8> = unit.iter().cycle().take(unit.len() * reps).copied().collect();
            let coded = encode(&[seq.len() as u32], &seq).expect("encode");
            let (_l, out_seq) = decode(&coded).expect("decode");
            proptest::prop_assert_eq!(out_seq, seq);
        }

        #[test]
        fn decode_never_panics(bytes in proptest::collection::vec(0u8..=255, 0..500)) {
            let _ = decode(&bytes);
        }
    }

    #[test]
    fn lzma_decode_rejects_huge_len_without_aborting() {
        // A corrupt decoded-length header must fail the reservation, not abort
        // the process on a huge infallible allocation.
        assert!(lzma_decode(&[], usize::MAX).is_err());
    }
}
