//! Nucleotide sequence coding via an order-k adaptive context model.
//!
//! Each base is one of A/C/G/T (2 bits) plus a fifth symbol for `N`; the model
//! conditions on the previous `k` ACGT bases and is range-coded ([`fqxv_range`]).
//! `N` — overwhelmingly the most common non-ACGT byte, and clustered in runs at
//! read ends — is coded directly by the model, so an `N` run costs almost
//! nothing and needs no side data. Rarer non-ACGT bytes (IUPAC codes, lowercase)
//! are coded as the same fifth symbol and their true byte is restored from a
//! small exception list, so the codec stays byte-exact. Non-ACGT symbols do
//! *not* advance the context: they neither pollute the model with a spurious
//! base nor cost one. Context carries across reads within a block (blocks stay
//! independent for parallelism).
//!
//! While an order-`k` context is still cold (few observations), the coder
//! escapes to a warm order-`k/2` model instead — most high-order contexts are
//! seen only a handful of times, so this cuts the cold-start tax. Both models
//! observe every symbol, so the decoder replays the same escape decisions.
//!
//! This is the sequence path for reads that are *not* reordered; the reordered
//! path lives in `fqxv-reorder`. Because the model is adaptive it uses the range
//! coder, not rANS (whose reverse encode can't carry adaptive state).
//!
//! ```
//! use fqxv_seq::{encode, decode};
//! let lens = [4u32, 5];
//! let seq = b"ACGTACGTN"; // note the trailing N
//! let enc = encode(&lens, seq, 6).unwrap();
//! let (out_lens, out_seq) = decode(&enc).unwrap();
//! assert_eq!(out_lens, lens);
//! assert_eq!(out_seq, seq);
//! ```

use fqxv_bytes::{read_lens, write_lens, write_varint, ReaderError};
use fqxv_range::{Decoder, Encoder};
use thiserror::Error;

/// Errors returned by the sequence codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed sequence stream: {0}")]
    Malformed(&'static str),
    /// The provided lengths do not sum to the sequence-buffer size.
    #[error("read lengths ({lens}) do not match sequence bytes ({seq})")]
    LengthMismatch {
        /// Sum of the provided read lengths.
        lens: usize,
        /// Number of sequence bytes provided.
        seq: usize,
    },
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

const FORMAT_VERSION: u8 = 1;

/// Loose upper bound on bases the range coder can emit per compressed byte, used
/// as a decompression-bomb guard in [`decode`]. The adaptive model caps its
/// total frequency at `1<<13`, so with 4 symbols the most-skewed base still
/// costs ~0.0005 bits (~15k bases/byte is the true ceiling); this leaves a wide
/// margin so legitimate streams never trip it.
const MAX_BASES_PER_BYTE: usize = 1 << 18;
/// Largest context order (4^11 contexts ≈ 4.2M models).
const MAX_ORDER: usize = 11;
/// Largest hashed high-order tier order. `2 * MAX_HASH_ORDER` bits must fit the
/// `u64` rolling context used to index the hash table.
const MAX_HASH_ORDER: usize = 24;
/// Largest hashed-tier table size (`1 << MAX_HASH_BITS` slots). Caps a hostile
/// header from requesting an enormous allocation before any decode happens.
const MAX_HASH_BITS: u32 = 26;
/// Odd multipliers for the hashed tier's index and collision-check (Fibonacci
/// hashing; the check uses an independent constant so index and check don't
/// correlate). See [`HashSlot`].
const HASH_MUL: u64 = 0x9E37_79B9_7F4A_7C15;
const CHECK_MUL: u64 = 0xD6E8_FEB8_6659_FD93;
/// Escape threshold: use the order-k model once its context total reaches this,
/// else fall back to the warm low-order model. A `NucModel` starts at `tot == 5`
/// and gains `STEP == 16` per observation, so this fires after ~4 observations —
/// below that the high-order context is too cold and the low-order model
/// predicts better. 69 was the sweep optimum across MiSeq/GAIIx/NovaSeq heads
/// (helps RNA-seq's many cold high-order contexts without over-escaping on
/// ultra-deep data). Tunable ratio knob.
const SEQ_ESCAPE_TOT: u32 = 69;

/// Order of the low-order fallback model, derived from the primary order so the
/// decoder reconstructs it from the stored `k` (no extra header byte).
#[inline]
fn lo_order(k: usize) -> usize {
    (k / 2).max(1)
}

/// byte -> 2-bit symbol, 255 for non-ACGT (coded as [`NSYM`]).
const BASE_LUT: [u8; 256] = base_lut();
const SYM2BASE: [u8; 4] = *b"ACGT";
/// The fifth model symbol: `N` (or any non-ACGT byte, restored via exceptions).
const NSYM: usize = 4;

const fn base_lut() -> [u8; 256] {
    let mut t = [255u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
}

/// Compact adaptive 5-symbol frequency model driving the range coder.
///
/// The alphabet is A/C/G/T plus a fifth "N/other" symbol ([`NSYM`]). Same
/// increment and halving cap as `fqxv_range::SimpleModel::<5>`, so it emits the
/// same `(cum, freq, tot)` intervals and the stream stays byte-exact — but it
/// stores *only* the five frequencies (10 bytes) instead of also caching the
/// total. The order-k context index walks this array with a near-random access
/// pattern over up to 4^11 ≈ 4.2M entries, which the profile shows dominates
/// both encode and decode; a small entry keeps the memory traffic down. The
/// total is a 5-add recompute — far cheaper than the miss it rides alongside.
#[derive(Clone)]
struct NucModel {
    freq: [u16; 5],
}

impl NucModel {
    /// Total frequency cap; matches `SimpleModel::MAX_TOT`.
    const MAX_TOT: u32 = 1 << 13;
    /// Frequency increment applied to the observed symbol.
    const STEP: u16 = 16;

    #[inline]
    fn new() -> Self {
        NucModel { freq: [1; 5] }
    }

    #[inline]
    fn tot(&self) -> u32 {
        let f = &self.freq;
        u32::from(f[0]) + u32::from(f[1]) + u32::from(f[2]) + u32::from(f[3]) + u32::from(f[4])
    }

    #[inline]
    fn encode(&mut self, enc: &mut Encoder, sym: usize) {
        let f = &self.freq;
        let mut cum = 0u32;
        for s in 0..sym {
            cum += u32::from(f[s]);
        }
        enc.encode(cum, u32::from(f[sym]), self.tot());
        self.update(sym);
    }

    #[inline]
    fn decode(&mut self, dec: &mut Decoder<'_>) -> usize {
        let target = dec.freq(self.tot());
        let f = &self.freq;
        let mut cum = 0u32;
        let mut sym = 0usize;
        while sym < 4 {
            let next = cum + u32::from(f[sym]);
            if target < next {
                break;
            }
            cum = next;
            sym += 1;
        }
        dec.decode(cum, u32::from(f[sym]));
        self.update(sym);
        sym
    }

    #[inline]
    fn update(&mut self, sym: usize) {
        self.freq[sym] += Self::STEP;
        if self.tot() > Self::MAX_TOT {
            for x in &mut self.freq {
                *x = (*x + 1) >> 1; // halve, keep >= 1
            }
        }
    }
}

/// One slot of the hashed high-order tier: a [`NucModel`] plus a `check` tag of
/// the context that owns it. On access, a slot whose `check` disagrees with the
/// current context is evicted (reset) — so a hash collision just costs recall (the
/// coder escapes to the dense order-k model), never correctness: encoder and
/// decoder derive `check` identically and evict in lock-step, so the stream stays
/// byte-exact. `check == 0` marks an empty slot ([`hash_check`] never returns 0).
#[derive(Clone)]
struct HashSlot {
    check: u32,
    m: NucModel,
}

impl HashSlot {
    #[inline]
    fn empty() -> Self {
        HashSlot {
            check: 0,
            m: NucModel::new(),
        }
    }
}

/// Table slot index for `ctx` (Fibonacci hash, top `bits` of the product).
#[inline]
fn hash_index(ctx: u64, bits: u32) -> usize {
    (ctx.wrapping_mul(HASH_MUL) >> (64 - bits)) as usize
}

/// Nonzero collision-check tag for `ctx`, independent of [`hash_index`].
#[inline]
fn hash_check(ctx: u64) -> u32 {
    (ctx.wrapping_mul(CHECK_MUL) as u32) | 1
}

/// Encode per-read sequences with an order-`order` context model.
pub fn encode(lens: &[u32], seq: &[u8], order: usize) -> Result<Vec<u8>> {
    encode_impl(lens, seq, order, 0, 0)
}

/// Encode with an extra **hashed** high-order tier above the dense order-`order`
/// model: an order-`hash_order` context indexed into a `1 << hash_bits`-slot
/// table (see [`HashSlot`]). Per base the warmest of `hashed → order-k → order-k/2`
/// codes the symbol; collisions in the table cost only recall (escape to the dense
/// tier), never correctness. `hash_order <= order` or `hash_bits == 0` disables it,
/// making this identical to [`encode`]. Byte-exactly reversible by [`decode`]
/// (the tier params are self-describing in the header).
pub fn encode_hashed(
    lens: &[u32],
    seq: &[u8],
    order: usize,
    hash_order: usize,
    hash_bits: u32,
) -> Result<Vec<u8>> {
    encode_impl(lens, seq, order, hash_order, hash_bits)
}

fn encode_impl(
    lens: &[u32],
    seq: &[u8],
    order: usize,
    hash_order: usize,
    hash_bits: u32,
) -> Result<Vec<u8>> {
    let total: usize = lens.iter().map(|&l| l as usize).sum();
    if total != seq.len() {
        return Err(Error::LengthMismatch {
            lens: total,
            seq: seq.len(),
        });
    }
    let k = order.clamp(1, MAX_ORDER);
    let klo = lo_order(k);
    let ctx_mask = (1usize << (2 * k)) - 1;
    let lo_mask = (1usize << (2 * klo)) - 1;
    // Optional hashed high-order tier. Only meaningful strictly above the dense
    // order-k model; `hash_bits` is clamped so the table stays bounded.
    let h = hash_order.min(MAX_HASH_ORDER);
    let hb = hash_bits.min(MAX_HASH_BITS);
    let use_hash = h > k && hb > 0;
    let h_mask: u64 = if use_hash { (1u64 << (2 * h)) - 1 } else { 0 };
    // The rolling context must retain the WIDEST tier's history (the hashed order
    // when enabled); each tier masks it down to its own order at lookup. Masking
    // it to `ctx_mask` here would strip the extra bases the hashed tier needs,
    // silently degrading that tier to an order-k copy.
    let top_mask: usize = if use_hash {
        (1usize << (2 * h)) - 1
    } else {
        ctx_mask
    };

    let mut hi = vec![NucModel::new(); ctx_mask + 1];
    let mut lo = vec![NucModel::new(); lo_mask + 1];
    let mut hh: Vec<HashSlot> = if use_hash {
        vec![HashSlot::empty(); 1usize << hb]
    } else {
        Vec::new()
    };
    let mut enc = Encoder::new();
    let mut exceptions: Vec<(usize, u8)> = Vec::new();
    let mut idx = 0usize;
    // Context carries across reads within this block (blocks stay independent).
    let mut ctx = 0usize;
    // Cursor over `seq` so the per-base read carries no bounds check; `idx` is
    // kept only to tag exception positions. `total == seq.len()` (checked
    // above), so every `split_at` is in range.
    let mut rest: &[u8] = seq;
    for &l in lens {
        let (read, tail) = rest.split_at(l as usize);
        rest = tail;
        for &byte in read {
            let raw = BASE_LUT[byte as usize];
            let sym = if raw == 255 { NSYM } else { raw as usize };
            let lc = ctx & lo_mask;
            let hi_ctx = ctx & ctx_mask;
            debug_assert!(hi_ctx <= ctx_mask && lc <= lo_mask && ctx <= top_mask);
            // Hashed slot for this context, evicting a colliding occupant first so
            // its total reads cold. `hslot` indexes `hh` for the rest of the base.
            let hslot = if use_hash {
                let hc = (ctx as u64) & h_mask;
                let idxh = hash_index(hc, hb);
                let ck = hash_check(hc);
                debug_assert!(idxh < hh.len());
                // SAFETY: `idxh < 1 << hb == hh.len()`.
                let s = unsafe { hh.get_unchecked_mut(idxh) };
                if s.check != ck {
                    s.check = ck;
                    s.m = NucModel::new();
                }
                idxh
            } else {
                0
            };
            // Code with the warmest tier (hashed → hi → lo); every tier observes
            // the symbol so encoder and decoder stay in sync.
            // SAFETY: `ctx <= ctx_mask == hi.len()-1`, `lc <= lo_mask == lo.len()-1`,
            // and `hslot < hh.len()` (set above; unread when `!use_hash`).
            unsafe {
                let hashed_warm = use_hash && hh.get_unchecked(hslot).m.tot() >= SEQ_ESCAPE_TOT;
                if hashed_warm {
                    hh.get_unchecked_mut(hslot).m.encode(&mut enc, sym);
                    hi.get_unchecked_mut(hi_ctx).update(sym);
                    lo.get_unchecked_mut(lc).update(sym);
                } else if hi.get_unchecked(hi_ctx).tot() >= SEQ_ESCAPE_TOT {
                    hi.get_unchecked_mut(hi_ctx).encode(&mut enc, sym);
                    if use_hash {
                        hh.get_unchecked_mut(hslot).m.update(sym);
                    }
                    lo.get_unchecked_mut(lc).update(sym);
                } else {
                    lo.get_unchecked_mut(lc).encode(&mut enc, sym);
                    if use_hash {
                        hh.get_unchecked_mut(hslot).m.update(sym);
                    }
                    hi.get_unchecked_mut(hi_ctx).update(sym);
                }
            }
            if raw == 255 {
                // Non-ACGT: `N` needs no side data; rarer bytes are recorded for
                // verbatim restore. Context does not advance — N is transparent.
                if byte != b'N' {
                    exceptions.push((idx, byte));
                }
            } else {
                ctx = ((ctx << 2) | sym) & top_mask;
            }
            idx += 1;
        }
    }
    let payload = enc.finish();

    let mut out = Vec::with_capacity(16 + lens.len() + exceptions.len() * 2 + payload.len());
    out.push(FORMAT_VERSION);
    out.push(k as u8);
    // Hashed-tier params (0 = none), self-describing so decode rebuilds the ladder.
    out.push(if use_hash { h as u8 } else { 0 });
    out.push(if use_hash { hb as u8 } else { 0 });
    write_lens(&mut out, lens);
    write_exceptions(&mut out, &exceptions);
    out.extend_from_slice(&payload);
    Ok(out)
}

/// Decode a stream produced by [`encode`], returning `(lengths, sequence)`.
pub fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut r = ByteReader::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let k = r.u8()? as usize;
    if !(1..=MAX_ORDER).contains(&k) {
        return Err(Error::Malformed("order out of range"));
    }
    // Hashed-tier params (0 = none). Validate before any shift/allocation, since
    // they come from an untrusted header.
    let h = r.u8()? as usize;
    let hb = u32::from(r.u8()?);
    let use_hash = h > k && hb > 0;
    if use_hash && (h > MAX_HASH_ORDER || hb > MAX_HASH_BITS) {
        return Err(Error::Malformed("hashed-tier params out of range"));
    }
    let klo = lo_order(k);
    let ctx_mask = (1usize << (2 * k)) - 1;
    let lo_mask = (1usize << (2 * klo)) - 1;
    let h_mask: u64 = if use_hash { (1u64 << (2 * h)) - 1 } else { 0 };
    let top_mask: usize = if use_hash {
        (1usize << (2 * h)) - 1
    } else {
        ctx_mask
    };
    let lens = read_lens(&mut r)?;
    let exceptions = read_exceptions(&mut r)?;

    let mut hi = vec![NucModel::new(); ctx_mask + 1];
    let mut lo = vec![NucModel::new(); lo_mask + 1];
    let mut hh: Vec<HashSlot> = Vec::new();
    if use_hash {
        // Fallible: the (capped) table must error, not abort, if it won't allocate.
        hh.try_reserve_exact(1usize << hb)
            .map_err(|_| Error::Malformed("hashed-tier table too large to allocate"))?;
        hh.resize(1usize << hb, HashSlot::empty());
    }
    let payload_len = r.rest().len();
    // Checked sum: a malformed stream can declare lengths whose total wraps
    // `usize`, which would under-allocate `seq` and then over-push.
    let total = lens
        .iter()
        .try_fold(0usize, |acc, &l| acc.checked_add(l as usize))
        .ok_or(Error::Malformed("total length overflows usize"))?;
    // Decompression-bomb guard. The range model caps its total frequency at
    // `1<<13`, so with 4 nucleotide symbols the coder emits at most ~15k bases
    // per compressed byte. A header that declares far more output than the
    // payload could possibly encode is malformed — reject before allocating or
    // looping so a tiny stream can't request a terabyte-scale decode.
    let max_plausible = payload_len
        .saturating_mul(MAX_BASES_PER_BYTE)
        .saturating_add(MAX_BASES_PER_BYTE);
    if total > max_plausible {
        return Err(Error::Malformed("declared length exceeds payload capacity"));
    }
    let mut dec = Decoder::new(r.rest());
    let mut seq = Vec::new();
    seq.try_reserve(total)
        .map_err(|_| Error::Malformed("declared total length too large to allocate"))?;
    let mut ctx = 0usize;
    for &l in &lens {
        for _ in 0..l {
            let lc = ctx & lo_mask;
            let hi_ctx = ctx & ctx_mask;
            debug_assert!(hi_ctx <= ctx_mask && lc <= lo_mask && ctx <= top_mask);
            // Hashed slot, evicting a collision first (mirrors the encoder).
            let hslot = if use_hash {
                let hc = (ctx as u64) & h_mask;
                let idxh = hash_index(hc, hb);
                let ck = hash_check(hc);
                debug_assert!(idxh < hh.len());
                // SAFETY: `idxh < 1 << hb == hh.len()`.
                let s = unsafe { hh.get_unchecked_mut(idxh) };
                if s.check != ck {
                    s.check = ck;
                    s.m = NucModel::new();
                }
                idxh
            } else {
                0
            };
            // SAFETY: `ctx <= ctx_mask == hi.len()-1`, `lc <= lo_mask == lo.len()-1`,
            // `hslot < hh.len()` (unread when `!use_hash`). Mirror the encoder's
            // tier choice exactly.
            let sym = unsafe {
                let hashed_warm = use_hash && hh.get_unchecked(hslot).m.tot() >= SEQ_ESCAPE_TOT;
                if hashed_warm {
                    let s = hh.get_unchecked_mut(hslot).m.decode(&mut dec);
                    hi.get_unchecked_mut(hi_ctx).update(s);
                    lo.get_unchecked_mut(lc).update(s);
                    s
                } else if hi.get_unchecked(hi_ctx).tot() >= SEQ_ESCAPE_TOT {
                    let s = hi.get_unchecked_mut(hi_ctx).decode(&mut dec);
                    if use_hash {
                        hh.get_unchecked_mut(hslot).m.update(s);
                    }
                    lo.get_unchecked_mut(lc).update(s);
                    s
                } else {
                    let s = lo.get_unchecked_mut(lc).decode(&mut dec);
                    if use_hash {
                        hh.get_unchecked_mut(hslot).m.update(s);
                    }
                    hi.get_unchecked_mut(hi_ctx).update(s);
                    s
                }
            };
            if sym == NSYM {
                // Default the N/other symbol to 'N'; the exception pass below
                // overwrites the rarer bytes. Context does not advance.
                seq.push(b'N');
            } else {
                seq.push(SYM2BASE[sym]);
                ctx = ((ctx << 2) | sym) & top_mask;
            }
        }
    }
    // Restore rarer non-ACGT bytes verbatim (N was already emitted).
    for (pos, byte) in exceptions {
        *seq.get_mut(pos)
            .ok_or(Error::Malformed("exception out of range"))? = byte;
    }
    Ok((lens, seq))
}

// --- length + exception side streams (LEB128 varints) ------------------------

fn write_exceptions(out: &mut Vec<u8>, exceptions: &[(usize, u8)]) {
    write_varint(out, exceptions.len() as u64);
    let mut prev = 0usize;
    for &(pos, byte) in exceptions {
        write_varint(out, (pos - prev) as u64);
        out.push(byte);
        prev = pos;
    }
}

fn read_exceptions(r: &mut ByteReader<'_>) -> Result<Vec<(usize, u8)>> {
    let n = r.varint()? as usize;
    let mut v = Vec::with_capacity(n.min(1 << 20));
    let mut pos = 0usize;
    for _ in 0..n {
        pos += r.varint()? as usize;
        let byte = r.u8()?;
        v.push((pos, byte));
    }
    Ok(v)
}

/// Shared byte cursor specialized to this crate's [`Error`].
type ByteReader<'a> = fqxv_bytes::Reader<'a, Error>;

impl ReaderError for Error {
    fn truncated() -> Self {
        Error::Malformed("truncated stream")
    }
    fn bad_varint() -> Self {
        Error::Malformed("varint too long")
    }
    fn oversized() -> Self {
        Error::Malformed("length count too large to allocate")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(lens: &[u32], seq: &[u8], order: usize) {
        let enc = encode(lens, seq, order).expect("encode");
        let (out_lens, out_seq) = decode(&enc).expect("decode");
        assert_eq!(out_lens, lens, "lengths mismatch");
        assert_eq!(out_seq, seq, "sequence mismatch (order {order})");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[], b"", 4);
    }

    #[test]
    fn roundtrip_basic() {
        roundtrip(&[4, 5], b"ACGTACGTN", 6);
    }

    #[test]
    fn roundtrip_with_exceptions() {
        roundtrip(&[10], b"ACGTNNRYAC", 4);
    }

    #[test]
    fn roundtrip_orders() {
        let seq: Vec<u8> = (0..1000).map(|i| SYM2BASE[(i * 7 + i / 3) % 4]).collect();
        for order in [1usize, 2, 4, 8, 11] {
            roundtrip(&[100; 10], &seq, order);
        }
    }

    #[test]
    fn n_runs_are_cheap() {
        // N is coded by the model (not one exception per base) and is
        // transparent to the context, so a poly-N tail must compress hard —
        // the whole point of the fifth symbol. 100 reads of ACGT(50) + N(50).
        let mut seq = Vec::new();
        let mut lens = Vec::new();
        for i in 0..100u32 {
            for j in 0..50 {
                seq.push(SYM2BASE[((i + j) % 4) as usize]);
            }
            seq.extend(std::iter::repeat_n(b'N', 50));
            lens.push(100);
        }
        let enc = encode(&lens, &seq, 8).expect("encode");
        // 5000 N bases must not cost anywhere near a byte each.
        assert!(
            enc.len() < 2_000,
            "poly-N should be nearly free, got {} bytes for {} bases",
            enc.len(),
            seq.len()
        );
        let (_, out) = decode(&enc).expect("decode");
        assert_eq!(out, seq);
    }

    #[test]
    fn repetitive_compresses() {
        // A repeated motif should compress well at a modest order.
        let unit = b"ACGTACGGTA";
        let seq: Vec<u8> = unit.iter().cycle().take(50_000).copied().collect();
        let enc = encode(&[50_000], &seq, 8).expect("encode");
        assert!(
            enc.len() < seq.len() / 4,
            "expected strong compression on repeats, got {} -> {}",
            seq.len(),
            enc.len()
        );
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(b"ACGTNacgtRYKM".to_vec()), 0..60),
                0..40),
            order in 1usize..=11,
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            roundtrip(&lens, &seq, order);
        }

        // Same, with the hashed tier on and a deliberately tiny table (hb=6) so
        // collisions and evictions are frequent — the eviction path must stay
        // byte-exact for arbitrary inputs.
        #[test]
        fn hashed_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(b"ACGTNacgtRYKM".to_vec()), 0..60),
                0..40),
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            let enc = encode_hashed(&lens, &seq, 11, 13, 6).expect("encode");
            let (out_lens, out_seq) = decode(&enc).expect("decode");
            proptest::prop_assert_eq!(out_lens, lens);
            proptest::prop_assert_eq!(out_seq, seq);
        }

        // Arbitrary bytes must never panic or abort the decoder — only Ok/Err.
        #[test]
        fn decode_never_aborts_on_garbage(bytes in proptest::collection::vec(0u8..=255, 0..256)) {
            let _ = decode(&bytes);
        }
    }

    fn roundtrip_hashed(lens: &[u32], seq: &[u8], order: usize, h: usize, hb: u32) {
        let enc = encode_hashed(lens, seq, order, h, hb).expect("encode");
        let (out_lens, out_seq) = decode(&enc).expect("decode");
        assert_eq!(out_lens, lens, "lengths mismatch (h={h} hb={hb})");
        assert_eq!(
            out_seq, seq,
            "sequence mismatch (order {order} h={h} hb={hb})"
        );
    }

    #[test]
    fn hashed_tier_roundtrips() {
        let seq: Vec<u8> = (0..5000u32)
            .map(|i| SYM2BASE[((i * 7 + i / 3) % 4) as usize])
            .collect();
        // Various table sizes, including tiny ones that force heavy collisions.
        for hb in [4u32, 8, 12, 16] {
            roundtrip_hashed(&[100; 50], &seq, 11, 13, hb);
        }
        // With N/exception bytes mixed in.
        roundtrip_hashed(&[20], b"ACGTNNRYACGTacgtACGT", 11, 13, 8);
    }

    #[test]
    fn hashed_disabled_matches_plain() {
        // hash_order <= order, or hash_bits == 0, must be byte-identical to plain
        // encode (the "levels 1-7 unchanged" guarantee).
        let seq: Vec<u8> = (0..3000u32)
            .map(|i| SYM2BASE[((i * 5 + 1) % 4) as usize])
            .collect();
        let lens = [150; 20];
        let plain = encode(&lens, &seq, 11).expect("plain");
        assert_eq!(encode_hashed(&lens, &seq, 11, 0, 0).expect("h0"), plain);
        assert_eq!(encode_hashed(&lens, &seq, 11, 13, 0).expect("hb0"), plain);
        assert_eq!(encode_hashed(&lens, &seq, 11, 11, 20).expect("h<=k"), plain);
    }

    #[test]
    fn hashed_encode_is_deterministic() {
        let seq: Vec<u8> = (0..8000u32)
            .map(|i| SYM2BASE[((i * 11 + i / 7) % 4) as usize])
            .collect();
        let a = encode_hashed(&[200; 40], &seq, 11, 13, 14).expect("a");
        let b = encode_hashed(&[200; 40], &seq, 11, 13, 14).expect("b");
        assert_eq!(a, b, "hashed encode must be reproducible");
    }

    #[test]
    fn hashed_beats_order11_on_order13_dependency() {
        // The next base is determined by the bases 12-13 back (the shared "AA"/"CC"
        // prefix) but NOT by the last 11 (a fixed 11-mer), so order-11 must pay ~1
        // bit at each such base while the order-13 hashed tier resolves it. This is
        // exactly the win the context-masking bug silently lost, so the hashed
        // stream must be clearly smaller than plain order-11.
        let p11: &[u8] = b"ACGTAACCGGT"; // 11 bases
        let mut seq = Vec::new();
        let mut lens = Vec::new();
        for i in 0..3000u32 {
            let bit = (i.wrapping_mul(2_654_435_761) >> 20) & 1 == 1;
            let (pre, nxt) = if bit { (b'A', b'G') } else { (b'C', b'T') };
            let mut read = vec![pre, pre];
            read.extend_from_slice(p11);
            read.push(nxt);
            lens.push(read.len() as u32);
            seq.extend_from_slice(&read);
        }
        let plain = encode(&lens, &seq, 11).expect("plain").len();
        let hashed = encode_hashed(&lens, &seq, 11, 13, 16)
            .expect("hashed")
            .len();
        assert!(
            hashed < plain,
            "order-13 hashed tier should beat order-11 here: {hashed} vs {plain}"
        );
        // Still byte-exact.
        let enc = encode_hashed(&lens, &seq, 11, 13, 16).unwrap();
        assert_eq!(decode(&enc).unwrap().1, seq);
    }

    fn push_varint(out: &mut Vec<u8>, mut v: u64) {
        loop {
            let b = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(b);
                break;
            }
            out.push(b | 0x80);
        }
    }

    // A hostile length header must produce a clean `Err`, never abort the
    // process on an impossible allocation.
    #[test]
    fn rejects_huge_length_count() {
        let mut buf = vec![FORMAT_VERSION, 1, 0, 0]; // version, order k=1, no hash tier
        push_varint(&mut buf, u64::MAX >> 8); // n: absurd length count
        buf.push(1); // fixed = true -> resize(n, f) path
        push_varint(&mut buf, 100); // f
        assert!(matches!(decode(&buf), Err(Error::Malformed(_))));
    }

    #[test]
    fn rejects_huge_total_length() {
        let mut buf = vec![FORMAT_VERSION, 1, 0, 0];
        push_varint(&mut buf, 1000); // n reads
        buf.push(1); // fixed = true
        push_varint(&mut buf, u32::MAX as u64); // each read u32::MAX long
        assert!(matches!(decode(&buf), Err(Error::Malformed(_))));
    }
}
