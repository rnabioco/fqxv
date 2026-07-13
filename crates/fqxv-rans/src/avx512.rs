//! AVX-512 order-0 rANS decoder and encoder.
//!
//! Byte-identical to [`crate::scalar`]; only throughput differs. The 32
//! interleaved states map onto two 512-bit registers per round (16 lanes each),
//! so this shares the AVX2 stream format and reuses its packed slot / magic
//! tables ([`crate::avx2::SlotTable`], [`crate::avx2::EncTables`]).
//!
//! AVX-512 buys two things over the AVX2 backend:
//!
//! - **Native `k`-mask compares.** `vpcmp*_epu32` yields an unsigned-compare
//!   mask directly, dropping the AVX2 XOR-with-sign-bit trick.
//! - **`vpexpandd` renorm.** Decode distributes the freshly read words into the
//!   lanes that need them with one `vpexpandd` (mask-driven expand), replacing
//!   the AVX2 precomputed permute table. Encode still emits with a scalar
//!   bit-scan because 16-bit `vpcompressw` needs AVX-512-VBMI2 (Ice Lake+),
//!   which the Skylake/Cascade Lake baseline lacks.

use std::arch::x86_64::*;

use crate::avx2::{EncTables, SlotTable};
use crate::model::{Model, N_STATES, RANS_L, SCALE_BITS, TOTFREQ};
use crate::scalar::Reader;
use crate::{Error, Result};

const LANES: usize = 16;

/// Decode an order-0 stream using AVX-512. The partial tail round is coded
/// scalar so the interleaved states stay in step with the encoder.
pub(crate) fn decode_order0(src: &[u8]) -> Result<Vec<u8>> {
    let mut r = Reader::new(src);
    let order = r.u8()?;
    let n = r.u64()? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    debug_assert_eq!(
        order, 0,
        "avx512::decode_order0 called on non-order-0 stream"
    );

    let mut freq = [0u16; 256];
    for f in &mut freq {
        *f = r.u16()?;
    }
    let model = Model::from_freqs(freq);
    let tables = SlotTable::from_model(&model);

    let mut states = [0u32; N_STATES];
    for x in &mut states {
        *x = r.u32()?;
    }

    let mut sp = r.pos();
    let mut out = vec![0u8; n];
    let full_rounds = n / N_STATES;
    // SAFETY: guarded by the AVX-512 feature detection in `crate::decode`.
    unsafe {
        decode_rounds(&tables, &mut states, &mut out, src, &mut sp, full_rounds)?;
    }

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

#[target_feature(enable = "avx512f")]
unsafe fn decode_rounds(
    tables: &SlotTable,
    states: &mut [u32; N_STATES],
    out: &mut [u8],
    src: &[u8],
    sp: &mut usize,
    full_rounds: usize,
) -> Result<()> {
    let mask12 = _mm512_set1_epi32((TOTFREQ - 1) as i32);
    let mask_off = _mm512_set1_epi32(0xfff);
    let mask_sym = _mm512_set1_epi32(0xff);
    let one = _mm512_set1_epi32(1);
    let rans_l = _mm512_set1_epi32(RANS_L as i32);
    let base_ptr = tables.packed.as_ptr();

    for round in 0..full_rounds {
        let base = round * N_STATES;
        for g in 0..(N_STATES / LANES) {
            let off = g * LANES;
            // SAFETY: `states[off..off+16]` is a valid 16-lane slice.
            let v = unsafe { _mm512_loadu_si512(states.as_ptr().add(off).cast()) };
            let slot = _mm512_and_si512(v, mask12);

            // One gather resolves freq, offset, and symbol.
            // SAFETY: slot indices are < TOTFREQ, in bounds for the packed table.
            let p = unsafe { _mm512_i32gather_epi32::<4>(slot, base_ptr) };
            let f = _mm512_add_epi32(_mm512_srli_epi32::<20>(p), one);
            let offset = _mm512_and_si512(_mm512_srli_epi32::<8>(p), mask_off);
            let sym = _mm512_and_si512(p, mask_sym);

            let vsh = _mm512_srli_epi32::<{ SCALE_BITS }>(v);
            let v2 = _mm512_add_epi32(_mm512_mullo_epi32(f, vsh), offset);

            // Lanes with v' < RANS_L pull one fresh 16-bit word; `vpexpandd`
            // scatters the k consumed words into exactly those lanes (ascending).
            let need = _mm512_cmplt_epu32_mask(v2, rans_l);
            let k = need.count_ones() as usize;
            let renormed = if k == 0 {
                v2
            } else {
                let words = unsafe { load_words(src, sp, k)? };
                let expanded = _mm512_maskz_expand_epi32(need, words);
                let shifted = _mm512_or_si512(_mm512_slli_epi32::<16>(v2), expanded);
                _mm512_mask_blend_epi32(need, v2, shifted)
            };

            // SAFETY: `states[off..off+16]` is a valid 16-lane target.
            unsafe { _mm512_storeu_si512(states.as_mut_ptr().add(off).cast(), renormed) };

            // Narrow the 16 symbol lanes to 16 bytes and store them.
            let sym_bytes = _mm512_cvtepi32_epi8(sym);
            // SAFETY: `base + off + 16 <= full_rounds*32 <= n`.
            unsafe { _mm_storeu_si128(out.as_mut_ptr().add(base + off).cast(), sym_bytes) };
        }
    }
    Ok(())
}

/// Load `k` (1..=16) consecutive little-endian 16-bit words from `src[*sp..]`,
/// zero-extended into the low `k` lanes of a `__m512i`, and advance `*sp` by
/// `2 * k`. Reads up to 32 bytes; falls back to a scalar fill near the end.
#[target_feature(enable = "avx512f")]
unsafe fn load_words(src: &[u8], sp: &mut usize, k: usize) -> Result<__m512i> {
    let end = *sp + 2 * k;
    if end > src.len() {
        return Err(Error::Malformed("truncated renorm stream"));
    }
    let words = if *sp + 32 <= src.len() {
        // SAFETY: 32 bytes from `*sp` are in bounds.
        let lo = unsafe { _mm256_loadu_si256(src.as_ptr().add(*sp).cast()) };
        _mm512_cvtepu16_epi32(lo)
    } else {
        let mut buf = [0u16; 16];
        for (j, w) in buf.iter_mut().enumerate().take(k) {
            let o = *sp + 2 * j;
            *w = u16::from_le_bytes([src[o], src[o + 1]]);
        }
        // SAFETY: `buf` is a valid 16-lane (32-byte) `[u16; 16]`.
        _mm512_cvtepu16_epi32(unsafe { _mm256_loadu_si256(buf.as_ptr().cast()) })
    };
    *sp = end;
    Ok(words)
}

/// Encode an order-0 stream using AVX-512. Byte-identical to
/// [`crate::scalar::encode_order0`]; the partial tail round is coded scalar.
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

    let full_rounds = n / N_STATES;
    for i in (full_rounds * N_STATES..n).rev() {
        let s = src[i] as usize;
        let x = &mut states[i % N_STATES];
        crate::scalar::encode_symbol(x, &enc[s], &mut renorm);
    }
    // SAFETY: guarded by the AVX-512 feature detection in `crate::encode`.
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

/// High 32 bits of each lane's unsigned 32×32→64 product, from the two
/// `vpmuludq` halves (even lanes, then odd lanes shifted into place).
#[target_feature(enable = "avx512f")]
#[inline]
fn mulhi_epu32(a: __m512i, b: __m512i) -> __m512i {
    let evn = _mm512_mul_epu32(a, b);
    let odd = _mm512_mul_epu32(_mm512_srli_epi64::<32>(a), _mm512_srli_epi64::<32>(b));
    let evn_hi = _mm512_srli_epi64::<32>(evn);
    let odd_hi = _mm512_slli_epi64::<32>(_mm512_srli_epi64::<32>(odd));
    _mm512_mask_blend_epi32(0xAAAA, evn_hi, odd_hi)
}

#[target_feature(enable = "avx512f")]
unsafe fn encode_rounds(
    tables: &EncTables,
    states: &mut [u32; N_STATES],
    src: &[u8],
    renorm: &mut Vec<u16>,
    full_rounds: usize,
) {
    let magic_ptr = tables.magic.as_ptr();
    let fc_ptr = tables.fc.as_ptr();
    let mask12 = _mm512_set1_epi32(0xfff);
    let mask16 = _mm512_set1_epi32(0xffff);
    let one = _mm512_set1_epi32(1);
    let mut wordtmp = [0i32; 16];

    for round in (0..full_rounds).rev() {
        let base = round * N_STATES;
        for g in (0..N_STATES / LANES).rev() {
            let off = g * LANES;

            // The 16 symbol bytes at these positions become gather indices.
            // SAFETY: `base + off + 16 <= n`.
            let sym_bytes = unsafe { _mm_loadu_si128(src.as_ptr().add(base + off).cast()) };
            let syms = _mm512_cvtepu8_epi32(sym_bytes);

            // SAFETY: symbol indices are < 256, in bounds for the 256-entry tables.
            let m = unsafe { _mm512_i32gather_epi32::<4>(syms, magic_ptr) };
            let fc = unsafe { _mm512_i32gather_epi32::<4>(syms, fc_ptr) };
            let freq = _mm512_srli_epi32::<12>(fc);
            let cum = _mm512_and_si512(fc, mask12);

            // SAFETY: `states[off..off+16]` is a valid 16-lane slice.
            let v0 = unsafe { _mm512_loadu_si512(states.as_ptr().add(off).cast()) };

            // Emit one 16-bit word from each lane whose state exceeds
            // x_max - 1 = (freq << 20) - 1, then shift that lane down by 16.
            let xmax_m1 = _mm512_sub_epi32(_mm512_slli_epi32::<20>(freq), one);
            let need = _mm512_cmpgt_epu32_mask(v0, xmax_m1);
            let words = _mm512_and_si512(v0, mask16);
            let v = _mm512_mask_srli_epi32(v0, need, v0, 16);

            // Emit surviving words in descending lane order (scalar bit-scan;
            // 16-bit vpcompressw would need VBMI2). At most one word per lane.
            let mut mask = need as u32;
            if mask != 0 {
                // SAFETY: `wordtmp` is a valid 16-lane target.
                unsafe { _mm512_storeu_si512(wordtmp.as_mut_ptr().cast(), words) };
                while mask != 0 {
                    let l = 31 - mask.leading_zeros() as usize;
                    renorm.push(wordtmp[l] as u16);
                    mask &= !(1 << l);
                }
            }

            // q = v / freq (estimate, then one correction), r = v - q*freq,
            // x' = (q << SCALE_BITS) + r + cum.
            let q0 = mulhi_epu32(v, m);
            let r0 = _mm512_sub_epi32(v, _mm512_mullo_epi32(q0, freq));
            let ge = _mm512_cmpge_epu32_mask(r0, freq);
            let q = _mm512_mask_add_epi32(q0, ge, q0, one);
            let r = _mm512_mask_sub_epi32(r0, ge, r0, freq);
            let x = _mm512_add_epi32(
                _mm512_add_epi32(_mm512_slli_epi32::<{ SCALE_BITS }>(q), r),
                cum,
            );

            // SAFETY: `states[off..off+16]` is a valid 16-lane target.
            unsafe { _mm512_storeu_si512(states.as_mut_ptr().add(off).cast(), x) };
        }
    }
}

#[inline]
fn read_word(src: &[u8], sp: &mut usize) -> Result<u16> {
    let b = src
        .get(*sp..*sp + 2)
        .ok_or(Error::Malformed("truncated renorm stream"))?;
    *sp += 2;
    Ok(u16::from_le_bytes([b[0], b[1]]))
}
