//! Fixed-point logistic (geometric) context mixing for long-read quality.
//!
//! The single-context adaptive model ([`crate::context_lr`]) is at its ceiling on
//! HiFi/ONT quality: adding context bits dilutes the per-block model faster than it
//! adds signal. This module instead mixes several context models of increasing
//! richness. For each quality symbol it forms, from every model `i`, a log-domain
//! feature `f_i[s] = ln P_i(s)`, blends them with adaptive weights
//! `logit(s) = Σ_i w[g][i]·f_i[s]`, takes the softmax as the coding distribution,
//! range-codes the symbol, then nudges the weights down the coding-loss gradient.
//! The weight set is selected by a gate `g` = how well-trained the richest model is
//! (secondary estimation): mixing trusts the rich model only where it has evidence.
//!
//! An entropy/coder simulation on full-alphabet HiFi put this below CoLoRd
//! (~3.16 vs ~3.19 bits/qual) where the single context plateaus at ~3.29.
//!
//! **Determinism.** Everything here is integer/fixed-point with lookup tables — no
//! float at runtime — so encode and decode reproduce bit-identical distributions on
//! any platform, and a future SIMD backend can match the scalar path exactly (the
//! crate's `SIMD ≡ scalar` invariant). Floating point is used only to *build* the
//! constant tables once.

use fqxv_range::{Decoder, Encoder};

use crate::{Error, Result};

/// Per-context model total-frequency cap (matches [`fqxv_range::SimpleModel`]).
const MODEL_MAX_TOT: u32 = 1 << 13;
/// Frequency increment on the observed symbol.
const STEP: u16 = 16;
/// Range-coder total for the mixed distribution. Must be `< 2^16` (the coder's
/// `BOT`); 2^15 gives ~15-bit probability precision, far finer than the mixer's
/// edge over the single context (~0.12 bits/qual).
const TOTAL: u32 = 1 << 15;

/// Fixed-point scale for natural-log feature values (`ln` in units of `1/LN_ONE`).
const LN_ONE: i32 = 1 << 16;
/// Fixed-point scale for mixing weights.
const W_ONE: i64 = 1 << 16;
/// Weight-update learning-rate numerator/shift: `w += (grad·LR_NUM) >> LR_SHIFT`.
/// `LR_NUM/2^LR_SHIFT ≈ 0.001`, the simulated optimum (lower lr won the sweep).
const LR_NUM: i64 = 1049;
const LR_SHIFT: u32 = 20;
/// Fixed-point shift for the softmax normalization reciprocal (`inv = spread·2^S/z`).
/// 32 bits keeps `e·inv` within `u64` and the quantization below one coder unit.
const RECIP_SHIFT: u32 = 32;
/// Weight clamp so a runaway gradient can't destabilize the mix.
const W_MIN: i64 = 0;
const W_MAX: i64 = 8 * W_ONE;

/// Number of gate buckets (weight sets), indexed by the richest model's training.
const GATES: usize = 8;

/// Rich-tier hash table size (log2). The 22-bit rich context is hashed into this
/// many slots to bound memory/init cost; collisions add mild noise the coarse
/// fallback absorbs.
const RICH_BITS: u32 = 20;

/// Exp table: `EXP[q] = round(exp(-(q<<EXP_SHIFT)/LN_ONE) · EXP_ONE)`. Covers the
/// range where the softmax weight is non-negligible; past it the symbol gets the
/// floor frequency 1.
const EXP_ONE: u32 = 1 << 15;
const EXP_SHIFT: u32 = 8;
const EXP_SIZE: usize = 2817; // ceil(11·LN_ONE / 2^EXP_SHIFT): exp(-11) < 1/EXP_ONE

/// Precomputed constant tables (built once from float; used only via integer ops).
struct Tables {
    /// `LN[v] = round(ln(v)·LN_ONE)` for `v` in `0..=MODEL_MAX_TOT` (`LN[0]` unused).
    ln: Vec<i32>,
    /// `EXP[q]` softmax weights (see [`EXP_SIZE`]).
    exp: Vec<u32>,
}

impl Tables {
    fn build() -> Self {
        let mut ln = vec![0i32; MODEL_MAX_TOT as usize + STEP as usize + 1];
        for (v, slot) in ln.iter_mut().enumerate().skip(1) {
            *slot = ((v as f64).ln() * f64::from(LN_ONE)).round() as i32;
        }
        let mut exp = vec![0u32; EXP_SIZE];
        for (q, slot) in exp.iter_mut().enumerate() {
            let d = (q << EXP_SHIFT) as f64 / f64::from(LN_ONE);
            *slot = ((-d).exp() * f64::from(EXP_ONE)).round() as u32;
        }
        Tables { ln, exp }
    }
}

/// One tier of context models: a flat frequency table over `n_ctx` contexts, each a
/// `k`-symbol adaptive histogram, plus a per-context total. `mask`/`hash` map a raw
/// context key into `0..n_ctx`.
struct Tier {
    freq: Vec<u16>, // n_ctx * k
    tot: Vec<u32>,  // n_ctx
    k: usize,
    mask: u32,
    hashed: bool,
}

impl Tier {
    fn new(bits: u32, k: usize, hashed: bool) -> Self {
        let n_ctx = 1usize << bits;
        Tier {
            freq: vec![1u16; n_ctx * k],
            tot: vec![k as u32; n_ctx],
            k,
            mask: (n_ctx as u32) - 1,
            hashed,
        }
    }

    #[inline]
    fn slot(&self, key: u32) -> usize {
        let idx = if self.hashed {
            // multiplicative (Fibonacci) hash, then keep the high bits.
            (key.wrapping_mul(2654435761) >> (32 - self.mask.count_ones())) & self.mask
        } else {
            key & self.mask
        };
        idx as usize * self.k
    }

    #[inline]
    fn ctx(&self, key: u32) -> usize {
        // context index (not the freq base) for the tot array
        if self.hashed {
            ((key.wrapping_mul(2654435761) >> (32 - self.mask.count_ones())) & self.mask) as usize
        } else {
            (key & self.mask) as usize
        }
    }

    #[inline]
    fn update(&mut self, key: u32, sym: usize) {
        let base = self.slot(key);
        let ci = self.ctx(key);
        self.freq[base + sym] += STEP;
        self.tot[ci] += u32::from(STEP);
        if self.tot[ci] > MODEL_MAX_TOT {
            let row = &mut self.freq[base..base + self.k];
            let mut t = 0u32;
            for f in row.iter_mut() {
                *f = (*f + 1) >> 1;
                t += u32::from(*f);
            }
            self.tot[ci] = t;
        }
    }
}

/// Long-read quality mixer: the model tiers, adaptive gated weights, and scratch
/// buffers reused across symbols. Encode and decode share this exact state machine.
struct Mixer {
    tiers: Vec<Tier>,
    /// weights[gate][tier], fixed-point `W_ONE`.
    weights: Vec<[i64; NMODELS]>,
    tables: Tables,
    k: usize,
    // scratch, reused each symbol to avoid per-symbol allocation. `feat` is
    // tier-major (`feat[t*k + s]`) so each tier's k features are contiguous —
    // the layout a SIMD backend loads 8/16 lanes at a time.
    feat: Vec<i32>, // NMODELS * k
    logit: Vec<i32>,
    e: Vec<u32>,
    freq_rc: Vec<u32>,
}

/// Number of mixed model tiers. Simulation: gated 3-tier {coarse16b, mid18b,
/// rich22b} at lr≈0.003 already beats CoLoRd; a 4th tier adds little.
const NMODELS: usize = 3;

impl Mixer {
    fn new(k: usize) -> Self {
        let tiers = vec![
            Tier::new(16, k, false), // coarse: q1>>1|q2>>3|base|next|hp
            Tier::new(18, k, false), // mid: coarse + prevbase
            Tier::new(RICH_BITS, k, true), // rich: mid + next2 + q3 (hashed)
        ];
        let init = W_ONE / NMODELS as i64;
        Mixer {
            tiers,
            weights: vec![[init; NMODELS]; GATES],
            tables: Tables::build(),
            k,
            feat: vec![0i32; NMODELS * k],
            logit: vec![0i32; k],
            e: vec![0u32; k],
            freq_rc: vec![0u32; k],
        }
    }

    /// Build the mixed coding distribution into `self.freq_rc` (summing to `TOTAL`),
    /// returning the gate index used (needed for the post-coding weight update).
    /// `keys[t]` is tier `t`'s raw context key for this symbol.
    #[inline]
    fn predict(&mut self, keys: &[u32; NMODELS]) -> usize {
        let k = self.k;
        // Per-tier log-prob features (tier-major) and the gate from the richest
        // tier's training. The k features of a tier are contiguous — a SIMD-shaped
        // layout, and the `tier.freq[base..base+k]` slice is one sequential scan.
        let mut gate = 0usize;
        for (t, tier) in self.tiers.iter().enumerate() {
            let base = tier.slot(keys[t]);
            let ci = tier.ctx(keys[t]);
            let ln_tot = self.tables.ln[tier.tot[ci] as usize];
            let dst = &mut self.feat[t * k..t * k + k];
            let row = &tier.freq[base..base + k];
            for (d, &f) in dst.iter_mut().zip(row) {
                *d = self.tables.ln[f as usize] - ln_tot;
            }
            if t == NMODELS - 1 {
                // log2(tot) bucketed: tot in [k..MODEL_MAX_TOT] -> 0..GATES-1.
                let lg = 31 - tier.tot[ci].leading_zeros();
                gate = (lg.saturating_sub(6) as usize).min(GATES - 1);
            }
        }
        let w = &self.weights[gate];
        // logit[s] = (Σ_t w[t]·feat[t*k+s]) >> 16 (back to Q16 nat-log units).
        let mut max_logit = i32::MIN;
        for s in 0..k {
            let mut l = 0i64;
            for t in 0..NMODELS {
                l += w[t] * i64::from(self.feat[t * k + s]);
            }
            let lg = (l >> 16) as i32;
            self.logit[s] = lg;
            if lg > max_logit {
                max_logit = lg;
            }
        }
        // softmax weights via the exp table
        let mut z = 0u64;
        for s in 0..k {
            let d = (max_logit - self.logit[s]) as u32; // >= 0, Q16
            let q = (d >> EXP_SHIFT) as usize;
            let e = if q < EXP_SIZE { self.tables.exp[q] } else { 0 };
            self.e[s] = e;
            z += u64::from(e);
        }
        // Normalize to TOTAL by reciprocal-multiply — no per-symbol divide (~93
        // integer divides/symbol were the dominant cost), and the exact integer
        // form a SIMD backend reproduces bit-for-bit. Floor keeps Σ ≤ TOTAL, so
        // the residue on the peak symbol is non-negative.
        let spread = u64::from(TOTAL) - k as u64;
        let z = z.max(1);
        let inv = (spread << RECIP_SHIFT) / z;
        let mut sum = 0u32;
        let mut argmax = 0usize;
        let mut max_e = 0u32;
        for s in 0..k {
            let f = ((u64::from(self.e[s]) * inv) >> RECIP_SHIFT) as u32 + 1;
            self.freq_rc[s] = f;
            sum += f;
            if self.e[s] > max_e {
                max_e = self.e[s];
                argmax = s;
            }
        }
        self.freq_rc[argmax] += TOTAL - sum;
        gate
    }

    /// After coding `sym`, adapt the weights (gradient of `-ln P(sym)`) and the
    /// per-tier models. Must be called identically on encode and decode.
    #[inline]
    fn update(&mut self, keys: &[u32; NMODELS], sym: usize, gate: usize) {
        let k = self.k;
        // E_P[f_i] = Σ_s P(s)·f_i[s] = (Σ_s e[s]·f_i[s]) / Z
        let z: u64 = self.e.iter().map(|&e| u64::from(e)).sum::<u64>().max(1);
        let w = &mut self.weights[gate];
        for t in 0..NMODELS {
            let mut acc = 0i64;
            let ft = &self.feat[t * k..t * k + k];
            for s in 0..k {
                acc += i64::from(self.e[s]) * i64::from(ft[s]);
            }
            let ep = acc / z as i64; // E_P[f_t], Q16
            let grad = i64::from(ft[sym]) - ep; // Q16
            let nw = w[t] + ((grad * LR_NUM) >> LR_SHIFT);
            w[t] = nw.clamp(W_MIN, W_MAX);
        }
        for (t, tier) in self.tiers.iter_mut().enumerate() {
            tier.update(keys[t], sym);
        }
    }
}

/// Compute the three context keys for a quality position from the sequence window
/// and previous-quality state. Mirrors [`crate::context_lr`]'s feature set, split
/// across the coarse/mid/rich tiers.
#[inline]
fn keys(q1: u8, q2: u8, q3: u8, base: usize, next: usize, next2: usize, prevbase: usize, hp: usize) -> [u32; 3] {
    let f1 = (q1 as u32 >> 1) & 0x3F;
    let q2c = (q2 as u32 >> 3) & 0x7;
    let q3c = (q3 as u32 >> 4) & 0x3;
    let b = base as u32 & 0x3;
    let nx = next as u32 & 0x3;
    let n2 = next2 as u32 & 0x3;
    let pb = prevbase as u32 & 0x3;
    let h = hp.min(7) as u32;
    let coarse = f1 | (q2c << 6) | (b << 9) | (nx << 11) | (h << 13); // 16 bits
    let mid = coarse | (pb << 16); // 18 bits
    let rich = mid | (n2 << 18) | (q3c << 20); // 22 bits (hashed by the tier)
    [coarse, mid, rich]
}

/// 2-bit base code (A/C/G/T -> 0/1/2/3; else 0), matching [`crate::base_code`].
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

/// Encode per-read qualities under the mixer, conditioning on the reads' bases.
/// `binned` is the (possibly binned) quality bytes; `seq` the concatenated bases in
/// lockstep; `dense` maps a quality byte to its dense symbol index; `qmin` the
/// context origin. Returns the range-coded payload.
pub(crate) fn encode(lens: &[u32], binned: &[u8], seq: &[u8], dense: &[u8; 256], qmin: u8, k: usize) -> Vec<u8> {
    let mut mx = Mixer::new(k);
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
            let kk = keys(q1, q2, q3, bcode(base), bcode(next), bcode(next2), bcode(prev_base), run);
            prev_base = base;
            let gate = mx.predict(&kk);
            let mut cum = 0u32;
            for s in 0..dv {
                cum += mx.freq_rc[s];
            }
            enc.encode(cum, mx.freq_rc[dv], TOTAL);
            mx.update(&kk, dv, gate);
            let cv = b - qmin;
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
    enc.finish()
}

/// Decode the mixer payload into `quals`, reconstructing the same features from the
/// caller-supplied decoded `seq`. `syms` maps a dense index back to its quality byte.
pub(crate) fn decode(
    lens: &[u32],
    payload: &[u8],
    seq: &[u8],
    syms: &[u8],
    qmin: u8,
    k: usize,
    quals: &mut Vec<u8>,
) -> Result<()> {
    let mut mx = Mixer::new(k);
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
            let kk = keys(q1, q2, q3, bcode(base), bcode(next), bcode(next2), bcode(prev_base), run);
            prev_base = base;
            let gate = mx.predict(&kk);
            let target = dec.freq(TOTAL);
            let mut cum = 0u32;
            let mut dv = 0usize;
            while dv + 1 < k && cum + mx.freq_rc[dv] <= target {
                cum += mx.freq_rc[dv];
                dv += 1;
            }
            dec.decode(cum, mx.freq_rc[dv]);
            mx.update(&kk, dv, gate);
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

    // Build a dense alphabet like the parent codec.
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
        let qmin = syms[0];
        let k = syms.len();
        (syms, map, qmin, k)
    }

    #[test]
    fn roundtrip_mixed() {
        // A few reads with base-correlated quality so the mixer has real signal.
        let bases = b"ACGTACGTAAAACCCCGGGGTTTTACGTACGTTTTTAAAA";
        let mut seq = Vec::new();
        let mut quals = Vec::new();
        let mut lens = Vec::new();
        for r in 0..5u32 {
            let l = 8 + r as usize;
            for i in 0..l {
                let bb = bases[(r as usize * 3 + i) % bases.len()];
                seq.push(bb);
                // quality correlates with base + a little variation
                let q = 33 + ((bcode(bb) * 7 + i + r as usize) % 40) as u8;
                quals.push(q);
            }
            lens.push(l as u32);
        }
        let (syms, dense, qmin, k) = dense_of(&quals);
        let payload = encode(&lens, &quals, &seq, &dense, qmin, k);
        let mut out = Vec::new();
        decode(&lens, &payload, &seq, &syms, qmin, k, &mut out).unwrap();
        assert_eq!(out, quals, "mixer round-trip must be lossless");
    }

    #[test]
    fn roundtrip_single_symbol() {
        // Degenerate alphabet (k=1): every freq_rc collapses to TOTAL on one symbol.
        let seq = vec![b'A'; 50];
        let quals = vec![b'I'; 50];
        let lens = vec![50u32];
        let (syms, dense, qmin, k) = dense_of(&quals);
        let payload = encode(&lens, &quals, &seq, &dense, qmin, k);
        let mut out = Vec::new();
        decode(&lens, &payload, &seq, &syms, qmin, k, &mut out).unwrap();
        assert_eq!(out, quals);
    }
}
