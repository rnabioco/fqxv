//! Binary-decomposition logistic-mixing quality coder for long reads.
//!
//! Codes each long-read quality symbol as `d = ceil(log2(k))` binary decisions down
//! a bit-tree instead of one k-way distribution over the alphabet. Each bit is
//! predicted by several context tiers of increasing richness whose
//! bit-probabilities are mixed in the logit domain (stretch/squash) with weights
//! adapted per bit-position — the lpaq/zpaq design. This is ~an order of magnitude
//! less work per symbol than the softmax path while keeping the same sub-CoLoRd
//! ratio (simulated: 3.15 vs CoLoRd 3.19 on full-alphabet HiFi).
//!
//! Fixed-point/integer throughout (12-bit probabilities, integer stretch/squash
//! tables), so encode and decode are bit-identical across platforms.

use fqxv_range::{Decoder, Encoder};

use crate::{Error, Result};

/// Probability precision: 12-bit, `PONE` = 1.0. Also the range-coder total (< 2^16).
const PBITS: u32 = 12;
const PONE: u32 = 1 << PBITS; // 4096
/// Number of context tiers mixed (coarse/mid/rich).
const NMODELS: usize = 3;
/// Rich context hashed into this many slots to bound the table.
const RICH_BITS: u32 = 20;
/// Mixer weight scale (Q16) and learning-rate shift. Weight update:
/// `w += (err·stretch) >> LR_SHIFT`. 11 is the joint HiFi/ONT optimum from a
/// real-byte sweep.
const W_ONE: i32 = 1 << 16;
const LR_SHIFT: u32 = 11;
/// Per-tier bit-probability adaptation rate: `p += ((bit<<12) - p) >> PRATE_SHIFT`
/// (1/32, the sim optimum).
const PRATE_SHIFT: u32 = 5;

/// Integer stretch/squash tables (lpaq): `squash(x)=4096/(1+e^(-x/256))` for
/// `x ∈ [-2047, 2047]`, and `stretch` its inverse mapping a 12-bit prob to `x`.
struct Logistic {
    squash: Vec<u16>,  // index x+2048 -> p in [1, 4095]
    stretch: Vec<i16>, // index p (0..4096) -> x in [-2047, 2047]
}

impl Logistic {
    fn build() -> Self {
        let mut squash = vec![0u16; 4096];
        for (i, s) in squash.iter_mut().enumerate() {
            let x = i as f64 - 2048.0;
            let p = (f64::from(PONE) / (1.0 + (-x / 256.0).exp())).round();
            *s = p.clamp(1.0, (PONE - 1) as f64) as u16;
        }
        // stretch = inverse of squash: for each prob, the smallest x with squash(x) >= p.
        let mut stretch = vec![0i16; (PONE + 1) as usize];
        let mut pi = 0usize;
        for x in 0..4096usize {
            let p = squash[x] as usize;
            while pi <= p {
                stretch[pi] = (x as i32 - 2048) as i16;
                pi += 1;
            }
        }
        while pi <= PONE as usize {
            stretch[pi] = 2047;
            pi += 1;
        }
        Logistic { squash, stretch }
    }

    #[inline]
    fn stretch(&self, p: u16) -> i32 {
        // SAFETY: `p` is a 12-bit probability, so `p <= PONE < stretch.len() = PONE+1`.
        i32::from(unsafe { *self.stretch.get_unchecked(p as usize) })
    }
    #[inline]
    fn squash(&self, x: i32) -> u16 {
        let xi = (x.clamp(-2047, 2047) + 2048) as usize;
        // SAFETY: `xi ∈ [1, 4095] < squash.len() = 4096`.
        unsafe { *self.squash.get_unchecked(xi) }
    }
}

/// A tier's per-context bit-tree of 12-bit probabilities: `probs[ctx*nnodes + node]`.
struct Tier {
    probs: Vec<u16>,
    nnodes: usize,
    mask: u32,
    hashed: bool,
}

impl Tier {
    fn new(bits: u32, nnodes: usize, hashed: bool) -> Self {
        let n_ctx = 1usize << bits;
        Tier {
            probs: vec![(PONE / 2) as u16; n_ctx * nnodes],
            nnodes,
            mask: (n_ctx as u32) - 1,
            hashed,
        }
    }
    #[inline]
    fn base(&self, key: u32) -> usize {
        let idx = if self.hashed {
            (key.wrapping_mul(2654435761) >> (32 - self.mask.count_ones())) & self.mask
        } else {
            key & self.mask
        };
        idx as usize * self.nnodes
    }
}

/// Binary-mixing coder state, shared by encode and decode.
struct BinMixer {
    tiers: [Tier; NMODELS],
    logistic: Logistic,
    weights: Vec<[i32; NMODELS]>, // per bit-position gate
    d: u32,
}

impl BinMixer {
    fn new(k: usize) -> Self {
        let d = if k <= 1 {
            1
        } else {
            u32::BITS - (k as u32 - 1).leading_zeros()
        };
        let nnodes = 1usize << d;
        let tiers = [
            Tier::new(16, nnodes, false),
            Tier::new(18, nnodes, false),
            Tier::new(RICH_BITS, nnodes, true),
        ];
        BinMixer {
            tiers,
            logistic: Logistic::build(),
            weights: vec![[W_ONE / NMODELS as i32; NMODELS]; d as usize],
            d,
        }
    }

    /// Predict P(bit=1) as a 12-bit prob for tree `node` under `keys`, and the
    /// per-tier stretched inputs (needed for the weight update). Returns
    /// `(p, [stretch_t], gate)`.
    /// The three tier base offsets for a context (constant across a symbol's bits;
    /// only the tree `node` changes), computed once per symbol.
    #[inline]
    fn bases(&self, keys: &[u32; NMODELS]) -> [usize; NMODELS] {
        core::array::from_fn(|t| self.tiers[t].base(keys[t]))
    }

    #[inline]
    fn predict(
        &self,
        bases: &[usize; NMODELS],
        node: usize,
        bpos: u32,
    ) -> (u16, [i32; NMODELS], usize) {
        let gate = bpos as usize;
        let mut st = [0i32; NMODELS];
        let mut x = 0i64;
        // SAFETY: `gate = bpos < d = weights.len()`; and `bases[t] + node <
        // tiers[t].probs.len()` — `bases[t] ≤ (n_ctx-1)*nnodes` and the tree node
        // stays in `[1, nnodes-1]` (`node = 2*node + bit` from 1, ≤ 2*(nnodes/2-1)+1).
        let w = unsafe { self.weights.get_unchecked(gate) };
        for t in 0..NMODELS {
            let p = unsafe { *self.tiers[t].probs.get_unchecked(bases[t] + node) };
            let s = self.logistic.stretch(p);
            st[t] = s;
            x += i64::from(w[t]) * i64::from(s);
        }
        let p = self.logistic.squash((x >> 16) as i32);
        (p.clamp(1, (PONE - 1) as u16), st, gate)
    }

    /// Adapt the mixer weights and each tier's node probability after coding `bit`.
    #[inline]
    fn update(
        &mut self,
        bases: &[usize; NMODELS],
        node: usize,
        bit: u32,
        p: u16,
        st: &[i32; NMODELS],
        gate: usize,
    ) {
        let err = (bit << PBITS) as i32 - i32::from(p); // [-4095, 4095]
        let target = (bit << PBITS) as i32;
        // SAFETY: same index invariants as `predict` (gate < weights.len(),
        // bases[t] + node < probs.len()).
        let w = unsafe { self.weights.get_unchecked_mut(gate) };
        for t in 0..NMODELS {
            w[t] += (err * st[t]) >> LR_SHIFT;
        }
        for t in 0..NMODELS {
            let idx = bases[t] + node;
            let slot = unsafe { self.tiers[t].probs.get_unchecked_mut(idx) };
            let cur = i32::from(*slot);
            *slot = (cur + ((target - cur) >> PRATE_SHIFT)) as u16;
        }
    }

    /// Code one symbol (`dv`, the dense index) as `d` bits, MSB first.
    #[inline]
    fn encode_sym(&mut self, enc: &mut Encoder, keys: &[u32; NMODELS], dv: usize) {
        let bases = self.bases(keys);
        let mut node = 1usize;
        for bpos in 0..self.d {
            let bit = (dv as u32 >> (self.d - 1 - bpos)) & 1;
            let (p, st, gate) = self.predict(&bases, node, bpos);
            // bit==1 occupies [0, p); bit==0 occupies [p, PONE).
            if bit == 1 {
                enc.encode(0, u32::from(p), PONE);
            } else {
                enc.encode(u32::from(p), PONE - u32::from(p), PONE);
            }
            self.update(&bases, node, bit, p, &st, gate);
            node = 2 * node + bit as usize;
        }
    }

    /// Decode one symbol, returning its dense index.
    #[inline]
    fn decode_sym(&mut self, dec: &mut Decoder<'_>, keys: &[u32; NMODELS]) -> usize {
        let bases = self.bases(keys);
        let mut node = 1usize;
        let mut dv = 0u32;
        for bpos in 0..self.d {
            let (p, st, gate) = self.predict(&bases, node, bpos);
            let target = dec.freq(PONE);
            let bit = if target < u32::from(p) { 1u32 } else { 0 };
            if bit == 1 {
                dec.decode(0, u32::from(p));
            } else {
                dec.decode(u32::from(p), PONE - u32::from(p));
            }
            self.update(&bases, node, bit, p, &st, gate);
            node = 2 * node + bit as usize;
            dv = (dv << 1) | bit;
        }
        dv as usize
    }
}

/// 2-bit base code, matching [`crate::base_code`].
#[inline]
fn bcode(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Per-block quantizer for the three recent-quality context fields, indexed by the
/// quality value on the encoder's `cv = byte - qmin` scale.
///
/// The three fields feed the coarse/mid/rich tiers at their fixed bit widths
/// (`q1` 6 bits, `q2` 3 bits, `q3` 2 bits). [`QCtx::flat`] reproduces the original
/// fixed `>>1`/`>>3`/`>>4` shifts exactly (byte-identical `MODE_SEQ_BINMIX`); a
/// per-block quantizer table (built by the caller from the block's histogram)
/// instead spends full resolution where quality actually varies — the fqzcomp
/// `qtab` idea (`MODE_SEQ_BINMIX_Q`). Because encode trials both and keeps the
/// smaller, the table can only match or shrink a block.
pub(crate) struct QCtx {
    g1: [u8; 256],
    g2: [u8; 256],
    g3: [u8; 256],
}

impl QCtx {
    /// The original fixed shifts, as tables: `g1[cv]=cv>>1`, `g2[cv]=cv>>3`,
    /// `g3[cv]=cv>>4`. [`keys`] masks each field, so this is byte-identical to the
    /// pre-quantizer coder.
    pub(crate) fn flat() -> Self {
        let mut g1 = [0u8; 256];
        let mut g2 = [0u8; 256];
        let mut g3 = [0u8; 256];
        for cv in 0..256usize {
            g1[cv] = (cv >> 1) as u8;
            g2[cv] = (cv >> 3) as u8;
            g3[cv] = (cv >> 4) as u8;
        }
        QCtx { g1, g2, g3 }
    }

    /// Build a quantizer from explicit per-cv bucket tables (see
    /// `crate::build_quant_ctx`).
    pub(crate) fn from_tables(g1: [u8; 256], g2: [u8; 256], g3: [u8; 256]) -> Self {
        QCtx { g1, g2, g3 }
    }
}

/// Context keys (coarse 16b, mid 18b, rich 22b) from the base window and recent
/// qualities — base identity and homopolymer runs, where long-read quality lives.
/// The three recent-quality fields arrive already quantized (see [`QCtx`]); [`keys`]
/// only masks them to their tier bit widths.
#[inline]
#[allow(clippy::too_many_arguments)] // a context-key builder over one feature per arg
fn keys(
    qc1: u32,
    qc2: u32,
    qc3: u32,
    base: usize,
    next: usize,
    next2: usize,
    prevbase: usize,
    hp: usize,
) -> [u32; 3] {
    let f1 = qc1 & 0x3F;
    let q2c = qc2 & 0x7;
    let q3c = qc3 & 0x3;
    let (b, nx, n2, pb) = (
        base as u32 & 3,
        next as u32 & 3,
        next2 as u32 & 3,
        prevbase as u32 & 3,
    );
    let h = hp.min(7) as u32;
    let coarse = f1 | (q2c << 6) | (b << 9) | (nx << 11) | (h << 13);
    let mid = coarse | (pb << 16);
    let rich = mid | (n2 << 18) | (q3c << 20);
    [coarse, mid, rich]
}

/// Encode per-read qualities with the binary-mixing coder under quantizer `qc`.
pub(crate) fn encode(
    lens: &[u32],
    binned: &[u8],
    seq: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    qc: &QCtx,
) -> Vec<u8> {
    let mut mx = BinMixer::new(k);
    let mut enc = Encoder::new();
    let mut rest = binned;
    let mut srest = seq;
    for &l in lens {
        let (read, tail) = rest.split_at(l as usize);
        rest = tail;
        let (sread, stail) = srest.split_at(l as usize);
        srest = stail;
        let (mut q1, mut q2, mut q3) = (0u8, 0u8, 0u8);
        let mut prev_base = u8::MAX;
        let mut run = 0usize;
        for (pos, &b) in read.iter().enumerate() {
            let dv = dense[b as usize] as usize;
            let base = sread[pos];
            let next = sread.get(pos + 1).copied().unwrap_or(u8::MAX);
            let next2 = sread.get(pos + 2).copied().unwrap_or(u8::MAX);
            run = if base == prev_base { run + 1 } else { 1 };
            let kk = keys(
                qc.g1[q1 as usize] as u32,
                qc.g2[q2 as usize] as u32,
                qc.g3[q3 as usize] as u32,
                bcode(base),
                bcode(next),
                bcode(next2),
                bcode(prev_base),
                run,
            );
            prev_base = base;
            mx.encode_sym(&mut enc, &kk, dv);
            let cv = b - qmin;
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
    enc.finish()
}

/// Decode a binary-mixing payload into `quals` under quantizer `qc`.
#[allow(clippy::too_many_arguments)] // parallel of `encode`, one arg per coder input
pub(crate) fn decode(
    lens: &[u32],
    payload: &[u8],
    seq: &[u8],
    syms: &[u8],
    qmin: u8,
    k: usize,
    qc: &QCtx,
    quals: &mut Vec<u8>,
) -> Result<()> {
    let mut mx = BinMixer::new(k);
    let mut dec = Decoder::new(payload);
    let mut srest = seq;
    for &l in lens {
        if srest.len() < l as usize {
            return Err(Error::Malformed("sequence shorter than quality lengths"));
        }
        let (sread, stail) = srest.split_at(l as usize);
        srest = stail;
        let (mut q1, mut q2, mut q3) = (0u8, 0u8, 0u8);
        let mut prev_base = u8::MAX;
        let mut run = 0usize;
        for pos in 0..l as usize {
            let base = sread[pos];
            let next = sread.get(pos + 1).copied().unwrap_or(u8::MAX);
            let next2 = sread.get(pos + 2).copied().unwrap_or(u8::MAX);
            run = if base == prev_base { run + 1 } else { 1 };
            let kk = keys(
                qc.g1[q1 as usize] as u32,
                qc.g2[q2 as usize] as u32,
                qc.g3[q3 as usize] as u32,
                bcode(base),
                bcode(next),
                bcode(next2),
                bcode(prev_base),
                run,
            );
            prev_base = base;
            let dv = mx.decode_sym(&mut dec, &kk);
            let b = *syms
                .get(dv)
                .ok_or(Error::Malformed("decoded symbol outside alphabet"))?;
            quals.push(b);
            let cv = b - qmin;
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense_of(quals: &[u8]) -> (Vec<u8>, [u8; 256], u8, usize) {
        let mut present = [false; 256];
        for &b in quals {
            present[b as usize] = true;
        }
        let syms: Vec<u8> = (0..=255u8).filter(|&b| present[b as usize]).collect();
        let mut map = [0u8; 256];
        for (i, &b) in syms.iter().enumerate() {
            map[b as usize] = i as u8;
        }
        (syms.clone(), map, syms[0], syms.len())
    }

    #[test]
    fn roundtrip_binmix() {
        let bases = b"ACGTACGTAAAACCCCGGGGTTTTACGTTTTTAAAACGCG";
        let (mut seq, mut quals, mut lens) = (Vec::new(), Vec::new(), Vec::new());
        for r in 0..6u32 {
            let l = 9 + r as usize;
            for i in 0..l {
                let bb = bases[(r as usize * 5 + i) % bases.len()];
                seq.push(bb);
                quals.push(33 + ((bcode(bb) * 9 + i * 2 + r as usize) % 42) as u8);
            }
            lens.push(l as u32);
        }
        let (syms, dense, qmin, k) = dense_of(&quals);
        let qc = QCtx::flat();
        let payload = encode(&lens, &quals, &seq, &dense, qmin, k, &qc);
        let mut out = Vec::new();
        decode(&lens, &payload, &seq, &syms, qmin, k, &qc, &mut out).unwrap();
        assert_eq!(out, quals, "binmix round-trip must be lossless");
    }

    #[test]
    fn roundtrip_binmix_tabled_quantizer() {
        // A non-flat quantizer must round-trip identically to the flat one — the
        // g-tables only relabel the context, never the coded symbol.
        let bases = b"ACGTACGTAAAACCCCGGGGTTTTACGTTTTTAAAACGCG";
        let (mut seq, mut quals, mut lens) = (Vec::new(), Vec::new(), Vec::new());
        for r in 0..6u32 {
            let l = 9 + r as usize;
            for i in 0..l {
                let bb = bases[(r as usize * 5 + i) % bases.len()];
                seq.push(bb);
                quals.push(33 + ((bcode(bb) * 9 + i * 2 + r as usize) % 42) as u8);
            }
            lens.push(l as u32);
        }
        let (syms, dense, qmin, k) = dense_of(&quals);
        // A hand-built tabled quantizer (dense rank for q1, coarse for q2/q3).
        let mut g1 = [0u8; 256];
        let mut g2 = [0u8; 256];
        let mut g3 = [0u8; 256];
        for (i, &b) in syms.iter().enumerate() {
            let cv = (b - qmin) as usize;
            g1[cv] = i as u8;
            g2[cv] = (i as u8) >> 2;
            g3[cv] = (i as u8) >> 3;
        }
        let qc = QCtx::from_tables(g1, g2, g3);
        let payload = encode(&lens, &quals, &seq, &dense, qmin, k, &qc);
        let mut out = Vec::new();
        decode(&lens, &payload, &seq, &syms, qmin, k, &qc, &mut out).unwrap();
        assert_eq!(out, quals, "tabled-quantizer round-trip must be lossless");
    }

    #[test]
    fn roundtrip_binmix_single_symbol() {
        let seq = vec![b'A'; 40];
        let quals = vec![b'~'; 40];
        let lens = vec![40u32];
        let (syms, dense, qmin, k) = dense_of(&quals);
        let qc = QCtx::flat();
        let payload = encode(&lens, &quals, &seq, &dense, qmin, k, &qc);
        let mut out = Vec::new();
        decode(&lens, &payload, &seq, &syms, qmin, k, &qc, &mut out).unwrap();
        assert_eq!(out, quals);
    }
}
