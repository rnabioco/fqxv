//! AVX2 order-0 rANS decoder and encoder.
//!
//! Produces byte-identical output to [`crate::scalar`]; only throughput
//! differs. The 32 interleaved states are processed 8 at a time (one AVX2
//! 256-bit register). Decode walks rounds in ascending state order (the order
//! the scalar reference consumes renorm words); encode walks them in descending
//! order (the order it emits, reversed by `finish`), so both share one stream
//! format. The encoder is described at [`encode_order0`].
//!
//! Two design choices keep this competitive with scalar on micro-architectures
//! with usable gather (Intel Broadwell/Skylake/Cascade Lake):
//!
//! - **One gather per group.** The frequency, cumulative-offset, and symbol
//!   lookups are packed into a single `u32` slot table (`(freq-1)<<20 |
//!   offset<<8 | sym`), so a group resolves with one `vpgatherdd` over a 16 KiB
//!   table that fits L1, instead of three gathers over 48 KiB.
//! - **Branchless renorm.** The Nx16 invariant guarantees at most one 16-bit
//!   word is consumed per state per step, so the data-dependent renorm becomes
//!   a compare + a permute that compacts the freshly read words into the lanes
//!   that need them (indexed by a 256-entry shuffle table), replacing the
//!   per-lane scalar loop.

use std::arch::x86_64::*;

use crate::model::{Model, N_STATES, RANS_L, SCALE_BITS, TOTFREQ};
use crate::scalar::Reader;
use crate::{Error, Result};

/// One `u32` per slot packing everything a decode step needs:
/// `((freq - 1) << 20) | (offset << 8) | sym`, where `offset = slot - cum[sym]`
/// is the position within the symbol's frequency range.
///
/// `freq` reaches `TOTFREQ` (4096, 13 bits) when one symbol owns the whole
/// alphabet, so we store `freq - 1` (12 bits) and add one back after unpacking.
/// `offset` and `sym` are 12 and 8 bits, so all three share a single word and
/// resolve with one gather over a 16 KiB (L1-resident) table.
pub(crate) struct SlotTable {
    pub(crate) packed: Vec<i32>,
}

impl SlotTable {
    pub(crate) fn from_model(model: &Model) -> Self {
        let n = TOTFREQ as usize;
        let mut packed = vec![0i32; n];
        for (slot, p) in packed.iter_mut().enumerate() {
            let s = model.slot2sym[slot];
            let f = u32::from(model.freq[s as usize]);
            let c = u32::from(model.cum[s as usize]);
            let offset = slot as u32 - c;
            *p = (((f - 1) << 20) | (offset << 8) | u32::from(s)) as i32;
        }
        SlotTable { packed }
    }
}

/// Decode an order-0 stream using AVX2. Falls back to per-symbol scalar work
/// for the tail that doesn't fill a 32-state round.
pub(crate) fn decode_order0(src: &[u8]) -> Result<Vec<u8>> {
    // Parse header (matches scalar's order-0 layout).
    let mut r = Reader::new(src);
    let order = r.u8()?;
    let n = r.u64()? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    debug_assert_eq!(order, 0, "avx2::decode_order0 called on non-order-0 stream");

    let mut freq = [0u16; 256];
    for f in &mut freq {
        *f = r.u16()?;
    }
    let model = Model::from_freqs(freq)?;
    let tables = SlotTable::from_model(&model);

    let mut states = [0u32; N_STATES];
    for x in &mut states {
        *x = r.u32()?;
    }

    // The renorm region is whatever follows the header.
    let mut sp = r.pos();
    let mut out = Model::alloc_output(n)?;

    let full_rounds = n / N_STATES;
    // SAFETY: guarded by the AVX2 feature detection in `crate::decode`.
    unsafe {
        decode_rounds(&tables, &mut states, &mut out, src, &mut sp, full_rounds)?;
    }

    // Scalar tail for the final partial round.
    for i in (full_rounds * N_STATES)..n {
        let x = &mut states[i % N_STATES];
        let slot = (*x & (TOTFREQ - 1)) as usize;
        let s = model.slot2sym[slot];
        out[i] = s;
        let f = u32::from(model.freq[s as usize]);
        let c = u32::from(model.cum[s as usize]);
        let mut v = f * (*x >> SCALE_BITS) + slot as u32 - c;
        while v < RANS_L {
            v = (v << 16) | u32::from(read_word(src, &mut sp)?);
        }
        *x = v;
    }
    Ok(out)
}

/// For each 8-bit "needs-renorm" mask, the `permutevar8x32` index vector that
/// gathers the consumed words into their lanes: lane `i` (if it needs a word)
/// reads word `popcount(mask & ((1<<i) - 1))`; other lanes read word 0 and are
/// discarded by the blend. 256 entries × 8 lanes = 8 KiB, built once per decode.
fn renorm_perm_table() -> [[i32; 8]; 256] {
    let mut table = [[0i32; 8]; 256];
    let mut mask = 0usize;
    while mask < 256 {
        let mut rank = 0i32;
        let mut i = 0;
        while i < 8 {
            if (mask >> i) & 1 == 1 {
                table[mask][i] = rank;
                rank += 1;
            }
            i += 1;
        }
        mask += 1;
    }
    table
}

#[target_feature(enable = "avx2")]
unsafe fn decode_rounds(
    tables: &SlotTable,
    states: &mut [u32; N_STATES],
    out: &mut [u8],
    src: &[u8],
    sp: &mut usize,
    full_rounds: usize,
) -> Result<()> {
    let perm = renorm_perm_table();
    let mask12 = _mm256_set1_epi32((TOTFREQ - 1) as i32);
    let mask_off = _mm256_set1_epi32(0xfff);
    let mask_sym = _mm256_set1_epi32(0xff);
    let one = _mm256_set1_epi32(1);
    let packed_ptr = tables.packed.as_ptr();

    for round in 0..full_rounds {
        let base = round * N_STATES;
        // Four 8-lane groups → 32 states, in ascending state order.
        for g in 0..(N_STATES / 8) {
            let off = g * 8;
            // SAFETY (memory intrinsics below): reached only under
            // `target_feature = "avx2"`. Loads/stores are the unaligned variants
            // over valid 8-lane slices, and the gather index is an in-bounds slot
            // (< TOTFREQ). Pure register intrinsics need no `unsafe` here.
            let v = unsafe { _mm256_loadu_si256(states.as_ptr().add(off).cast()) };
            let slot = _mm256_and_si256(v, mask12);

            // One gather resolves freq, offset, and symbol.
            let p = unsafe { _mm256_i32gather_epi32::<4>(packed_ptr, slot) };
            let f = _mm256_add_epi32(_mm256_srli_epi32::<20>(p), one);
            let offset = _mm256_and_si256(_mm256_srli_epi32::<8>(p), mask_off);
            let sym = _mm256_and_si256(p, mask_sym);

            // v' = freq * (v >> SCALE_BITS) + offset  (mod 2^32, as scalar does).
            let vsh = _mm256_srli_epi32::<{ SCALE_BITS as i32 }>(v);
            let v2 = _mm256_add_epi32(_mm256_mullo_epi32(f, vsh), offset);

            // Renormalize: lanes with v' < RANS_L (top 16 bits zero) pull one
            // fresh 16-bit word. Compact the consumed words into those lanes.
            let need = _mm256_cmpeq_epi32(_mm256_srli_epi32::<16>(v2), _mm256_setzero_si256());
            let m = _mm256_movemask_ps(_mm256_castsi256_ps(need)) as usize;
            let k = m.count_ones() as usize;

            let renormed = if k == 0 {
                v2
            } else {
                // SAFETY: reached only under `target_feature = "avx2"`.
                let words = unsafe { load_words(src, sp, k)? };
                // SAFETY: `perm[m]` is a valid 8-lane `[i32; 8]` slice.
                let idx = unsafe { _mm256_loadu_si256(perm[m].as_ptr().cast()) };
                let placed = _mm256_permutevar8x32_epi32(words, idx);
                // v' = (v' << 16) | word, only where needed.
                let shifted = _mm256_or_si256(_mm256_slli_epi32::<16>(v2), placed);
                _mm256_blendv_epi8(v2, shifted, need)
            };

            // SAFETY: `states[off..off+8]` and `symtmp` are valid 8-lane targets.
            unsafe { _mm256_storeu_si256(states.as_mut_ptr().add(off).cast(), renormed) };

            // Pack the 8 symbol lanes down to bytes and store them.
            let mut symtmp = [0i32; 8];
            unsafe { _mm256_storeu_si256(symtmp.as_mut_ptr().cast(), sym) };
            for l in 0..8 {
                out[base + off + l] = symtmp[l] as u8;
            }
        }
    }
    Ok(())
}

/// Load `k` (1..=8) consecutive little-endian 16-bit words from `src[*sp..]`,
/// zero-extended into the low `k` lanes of a `__m256i`, and advance `*sp` by
/// `2 * k`. Reads up to 16 bytes; callers stay `2*8` bytes inside `src`.
#[target_feature(enable = "avx2")]
unsafe fn load_words(src: &[u8], sp: &mut usize, k: usize) -> Result<__m256i> {
    let end = *sp + 2 * k;
    if end > src.len() {
        return Err(Error::Malformed("truncated renorm stream"));
    }
    // Fast path: a full 16-byte load when the eight words are all in bounds.
    let words = if *sp + 16 <= src.len() {
        // SAFETY: `*sp + 16 <= src.len()`, so 16 bytes from `*sp` are in bounds.
        let lo = unsafe { _mm_loadu_si128(src.as_ptr().add(*sp).cast()) };
        _mm256_cvtepu16_epi32(lo)
    } else {
        let mut buf = [0u16; 8];
        for (j, w) in buf.iter_mut().enumerate().take(k) {
            let o = *sp + 2 * j;
            *w = u16::from_le_bytes([src[o], src[o + 1]]);
        }
        // SAFETY: `buf` is a valid 8-lane (16-byte) `[u16; 8]`.
        _mm256_cvtepu16_epi32(unsafe { _mm_loadu_si128(buf.as_ptr().cast()) })
    };
    *sp = end;
    Ok(words)
}

#[inline]
fn read_word(src: &[u8], sp: &mut usize) -> Result<u16> {
    let b = src
        .get(*sp..*sp + 2)
        .ok_or(Error::Malformed("truncated renorm stream"))?;
    *sp += 2;
    Ok(u16::from_le_bytes([b[0], b[1]]))
}

/// Per-symbol encode constants, indexed by symbol byte (256-entry L1 tables):
/// `magic` is a division-by-`freq` reciprocal, `fc` packs `(freq << 12) | cum`.
///
/// The rANS encode map is always `x' = (v / freq) << SCALE_BITS + v % freq + cum`
/// regardless of how the quotient is obtained, so instead of the scalar coder's
/// reciprocal-or-division split we compute `q = v / freq` exactly for every lane
/// with one estimate-and-correct step: `magic = floor((2^32 - 1) / freq)` gives
/// `q0 = mulhi(v, magic) ∈ {q - 1, q}` (never an overestimate), and a single
/// conditional `+1` corrects it. Branchless, and byte-identical to scalar.
pub(crate) struct EncTables {
    pub(crate) magic: Vec<i32>,
    pub(crate) fc: Vec<i32>,
}

impl EncTables {
    pub(crate) fn from_model(model: &Model) -> Self {
        let mut magic = vec![0i32; 256];
        let mut fc = vec![0i32; 256];
        for s in 0..256 {
            let f = u32::from(model.freq[s]);
            if f == 0 {
                continue; // absent symbol: never encoded, placeholders are fine.
            }
            magic[s] = (((1u64 << 32) - 1) / u64::from(f)) as i32;
            fc[s] = ((f << 12) | u32::from(model.cum[s])) as i32;
        }
        EncTables { magic, fc }
    }
}

/// Encode an order-0 stream using AVX2. Byte-identical to
/// [`crate::scalar::encode_order0`]; the partial tail round is coded with the
/// scalar coder so the interleaved states stay in step.
pub(crate) fn encode_order0(src: &[u8]) -> Vec<u8> {
    let mut counts = [0u32; 256];
    for &b in src {
        counts[b as usize] += 1;
    }
    let model = Model::from_counts(&counts);
    let enc = model.enc_table();
    let tables = EncTables::from_model(&model);

    let n = src.len();
    let mut states = [RANS_L; N_STATES];
    let mut renorm: Vec<u16> = Vec::with_capacity(n);

    // The scalar coder walks symbols in reverse; the SIMD path must push renorm
    // words in that same descending-position order so the final reverse in
    // `finish` yields the order the forward decoder consumes. So: the partial
    // tail (highest positions) first, then full rounds descending.
    let full_rounds = n / N_STATES;
    for i in (full_rounds * N_STATES..n).rev() {
        let s = src[i] as usize;
        let x = &mut states[i % N_STATES];
        crate::scalar::encode_symbol(x, &enc[s], &mut renorm);
    }
    // SAFETY: guarded by the AVX2 feature detection in `crate::encode`.
    unsafe {
        encode_rounds(&tables, &mut states, src, &mut renorm, full_rounds);
    }

    let mut out = Vec::new();
    out.push(crate::scalar::ORDER0);
    out.extend_from_slice(&(n as u64).to_le_bytes());
    if n != 0 {
        for s in 0..256 {
            out.extend_from_slice(&model.freq[s].to_le_bytes());
        }
        crate::scalar::finish(&mut out, &states, &mut renorm);
    }
    out
}

/// High 32 bits of each lane's unsigned 32×32→64 product (`mulhi_epu32`), built
/// from the two `vpmuludq` halves (even lanes, then odd lanes shifted into place).
#[target_feature(enable = "avx2")]
#[inline]
fn mulhi_epu32(a: __m256i, b: __m256i) -> __m256i {
    let evn = _mm256_mul_epu32(a, b);
    let a_odd = _mm256_srli_epi64::<32>(a);
    let b_odd = _mm256_srli_epi64::<32>(b);
    let odd = _mm256_mul_epu32(a_odd, b_odd);
    let evn_hi = _mm256_srli_epi64::<32>(evn);
    let odd_hi = _mm256_slli_epi64::<32>(_mm256_srli_epi64::<32>(odd));
    _mm256_blend_epi32::<0b1010_1010>(evn_hi, odd_hi)
}

#[target_feature(enable = "avx2")]
unsafe fn encode_rounds(
    tables: &EncTables,
    states: &mut [u32; N_STATES],
    src: &[u8],
    renorm: &mut Vec<u16>,
    full_rounds: usize,
) {
    let magic_ptr = tables.magic.as_ptr();
    let fc_ptr = tables.fc.as_ptr();
    let mask12 = _mm256_set1_epi32(0xfff);
    let mask16 = _mm256_set1_epi32(0xffff);
    let one = _mm256_set1_epi32(1);
    let sign = _mm256_set1_epi32(i32::MIN); // 0x8000_0000 for unsigned compare

    let mut wordtmp = [0i32; 8];

    // Descending rounds, descending groups, descending lanes — reverse of the
    // decoder's consume order (see `encode_order0`).
    for round in (0..full_rounds).rev() {
        let base = round * N_STATES;
        for g in (0..N_STATES / 8).rev() {
            let off = g * 8;

            // The 8 symbol bytes at these positions become gather indices.
            // SAFETY: `base + off + 8 <= n`; 8 bytes are in bounds.
            let sym_bytes = unsafe { _mm_loadl_epi64(src.as_ptr().add(base + off).cast()) };
            let syms = _mm256_cvtepu8_epi32(sym_bytes);

            // SAFETY: symbol indices are < 256, in-bounds for the 256-entry tables.
            let m = unsafe { _mm256_i32gather_epi32::<4>(magic_ptr, syms) };
            let fc = unsafe { _mm256_i32gather_epi32::<4>(fc_ptr, syms) };
            let freq = _mm256_srli_epi32::<12>(fc);
            let cum = _mm256_and_si256(fc, mask12);

            // SAFETY: `states[off..off+8]` is a valid 8-lane slice.
            let v0 = unsafe { _mm256_loadu_si256(states.as_ptr().add(off).cast()) };

            // Renorm: emit one 16-bit word from each lane whose state exceeds
            // x_max - 1 = (freq << 20) - 1, then shift that lane down by 16.
            let xmax_m1 = _mm256_sub_epi32(_mm256_slli_epi32::<20>(freq), one);
            let need =
                _mm256_cmpgt_epi32(_mm256_xor_si256(v0, sign), _mm256_xor_si256(xmax_m1, sign));
            let words = _mm256_and_si256(v0, mask16);
            let v = _mm256_blendv_epi8(v0, _mm256_srli_epi32::<16>(v0), need);

            // Emit the surviving words in descending lane order (matching the
            // scalar reverse walk). At most one word per lane (Nx16 invariant);
            // iterate only the set lanes (high bit first) rather than all eight.
            let mut mask = _mm256_movemask_ps(_mm256_castsi256_ps(need)) as u32;
            if mask != 0 {
                // SAFETY: `wordtmp` is a valid 8-lane target.
                unsafe { _mm256_storeu_si256(wordtmp.as_mut_ptr().cast(), words) };
                while mask != 0 {
                    let l = 31 - mask.leading_zeros() as usize; // highest set lane
                    renorm.push(wordtmp[l] as u16);
                    mask &= !(1 << l);
                }
            }

            // Apply the rANS map: q = v / freq (estimate, then one correction),
            // r = v - q*freq, x' = (q << SCALE_BITS) + r + cum.
            let q0 = mulhi_epu32(v, m);
            let r0 = _mm256_sub_epi32(v, _mm256_mullo_epi32(q0, freq));
            let ge = _mm256_cmpgt_epi32(r0, _mm256_sub_epi32(freq, one));
            let q = _mm256_sub_epi32(q0, ge); // ge == -1 → +1
            let r = _mm256_sub_epi32(r0, _mm256_and_si256(ge, freq));
            let x = _mm256_add_epi32(
                _mm256_add_epi32(_mm256_slli_epi32::<{ SCALE_BITS as i32 }>(q), r),
                cum,
            );

            // SAFETY: `states[off..off+8]` is a valid 8-lane target.
            unsafe { _mm256_storeu_si256(states.as_mut_ptr().add(off).cast(), x) };
        }
    }
}
