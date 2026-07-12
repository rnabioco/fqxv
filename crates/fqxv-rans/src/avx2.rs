//! AVX2 order-0 rANS decoder.
//!
//! Produces byte-identical output to [`crate::scalar`]; only throughput
//! differs. The 32 interleaved states are processed 8 at a time (one AVX2
//! 256-bit register). Each round decodes 32 symbols in ascending state order,
//! which is exactly the order the scalar reference consumes renorm words, so
//! the two share one stream format.
//!
//! The gather-heavy work — slot masking, and the frequency / cumulative /
//! symbol table lookups — is vectorized; the data-dependent 16-bit renorm read
//! is done with a short scalar loop over the 8 lanes (ascending), which keeps
//! the stream cursor exactly in step with the encoder. A fully vectorized
//! renorm compaction is a future optimization.

use std::arch::x86_64::*;

use crate::model::{Model, N_STATES, RANS_L, SCALE_BITS, TOTFREQ};
use crate::scalar::Reader;
use crate::{Error, Result};

/// Slot-indexed lookup tables, so a single gather per table resolves a state.
struct SlotTables {
    freq: Vec<i32>, // freq[slot]
    cum: Vec<i32>,  // cum[slot]
    sym: Vec<i32>,  // symbol byte, zero-extended
}

impl SlotTables {
    fn from_model(model: &Model) -> Self {
        let n = TOTFREQ as usize;
        let mut freq = vec![0i32; n];
        let mut cum = vec![0i32; n];
        let mut sym = vec![0i32; n];
        for slot in 0..n {
            let s = model.slot2sym[slot];
            sym[slot] = i32::from(s);
            freq[slot] = i32::from(model.freq[s as usize]);
            cum[slot] = i32::from(model.cum[s as usize]);
        }
        SlotTables { freq, cum, sym }
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
    let model = Model::from_freqs(freq);
    let tables = SlotTables::from_model(&model);

    let mut states = [0u32; N_STATES];
    for x in &mut states {
        *x = r.u32()?;
    }

    // The renorm region is whatever follows the header.
    let mut sp = r.pos();
    let mut out = vec![0u8; n];

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

#[target_feature(enable = "avx2")]
unsafe fn decode_rounds(
    tables: &SlotTables,
    states: &mut [u32; N_STATES],
    out: &mut [u8],
    src: &[u8],
    sp: &mut usize,
    full_rounds: usize,
) -> Result<()> {
    let mask = _mm256_set1_epi32((TOTFREQ - 1) as i32);
    let mut tmp = [0u32; 8];
    let mut symtmp = [0i32; 8];

    for round in 0..full_rounds {
        let base = round * N_STATES;
        // Four 8-lane groups → 32 states, in ascending state order.
        for g in 0..(N_STATES / 8) {
            let off = g * 8;
            // SAFETY: reached only under `target_feature = "avx2"`. Loads/stores
            // are the unaligned variants over valid 8-lane slices, and gather
            // indices are in-bounds slots (< TOTFREQ).
            unsafe {
                let v = _mm256_loadu_si256(states.as_ptr().add(off).cast());
                let slot = _mm256_and_si256(v, mask);

                // One gather per slot-indexed table.
                let f = _mm256_i32gather_epi32::<4>(tables.freq.as_ptr(), slot);
                let c = _mm256_i32gather_epi32::<4>(tables.cum.as_ptr(), slot);
                let sym = _mm256_i32gather_epi32::<4>(tables.sym.as_ptr(), slot);

                // v' = f * (v >> SCALE_BITS) + slot - c
                let vsh = _mm256_srli_epi32::<{ SCALE_BITS as i32 }>(v);
                let fm = _mm256_mullo_epi32(f, vsh);
                let v2 = _mm256_add_epi32(_mm256_sub_epi32(fm, c), slot);

                _mm256_storeu_si256(tmp.as_mut_ptr().cast(), v2);
                _mm256_storeu_si256(symtmp.as_mut_ptr().cast(), sym);
            }

            // Scalar renorm per lane (ascending), then write states + output.
            for l in 0..8 {
                let mut x = tmp[l];
                while x < RANS_L {
                    x = (x << 16) | u32::from(read_word(src, sp)?);
                }
                states[off + l] = x;
                out[base + off + l] = symtmp[l] as u8;
            }
        }
    }
    Ok(())
}

#[inline]
fn read_word(src: &[u8], sp: &mut usize) -> Result<u16> {
    let b = src
        .get(*sp..*sp + 2)
        .ok_or(Error::Malformed("truncated renorm stream"))?;
    *sp += 2;
    Ok(u16::from_le_bytes([b[0], b[1]]))
}
