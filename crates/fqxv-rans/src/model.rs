//! Order-0 frequency model.
//!
//! Byte counts are normalized so the frequencies sum to exactly [`TOTFREQ`]
//! (a 12-bit total), which is the invariant the rANS coder relies on. We also
//! build the cumulative table and the `slot -> symbol` lookup the decoder uses.

/// Number of bits in the frequency total (the rANS "scale").
pub(crate) const SCALE_BITS: u32 = 12;
/// The frequency total; every model sums to this.
pub(crate) const TOTFREQ: u32 = 1 << SCALE_BITS;
/// Lower bound of the normalized rANS state interval `[RANS_L, RANS_L << 16)`.
///
/// With 16-bit renormalization this bound guarantees at most one 16-bit word is
/// emitted/consumed per state per step — the property the SIMD backends rely on.
pub(crate) const RANS_L: u32 = 1 << 16;
/// Number of interleaved rANS states (maps to SIMD lanes later).
pub(crate) const N_STATES: usize = 32;

/// A normalized order-0 model over the 256-byte alphabet.
#[derive(Debug)]
pub(crate) struct Model {
    /// Per-symbol frequency, summing to [`TOTFREQ`] (0 for absent symbols).
    pub(crate) freq: [u16; 256],
    /// Cumulative frequency; `cum[s]` is the start slot of symbol `s`,
    /// `cum[256] == TOTFREQ`.
    pub(crate) cum: [u16; 257],
    /// `slot -> symbol` reverse lookup, length [`TOTFREQ`].
    pub(crate) slot2sym: Vec<u8>,
}

/// Precomputed per-symbol constants for division-free rANS encoding.
///
/// The rANS encode map is `x' = (x / f) << SCALE_BITS + (x % f) + c`. The naive
/// form runs an integer `DIV`+`MOD` per symbol (tens of cycles on the critical
/// path). We replace `x / f` with an Alverson reciprocal multiply — precompute a
/// fixed-point reciprocal per symbol and turn the division into a multiply-high
/// plus shift. This is the ryg_rans `RansEncSymbolInit` construction and yields
/// output byte-identical to the division form (see the equivalence test below).
#[derive(Clone, Copy, Debug)]
pub(crate) struct EncSym {
    /// Renormalization bound: emit 16-bit words while state `>= x_max`.
    pub(crate) x_max: u64,
    rcp_freq: u32,
    rcp_shift: u32,
    bias: u32,
    cmpl_freq: u32,
}

impl EncSym {
    /// Build the reciprocal constants for a symbol with frequency `freq` and
    /// cumulative start `cum`. Absent symbols (`freq == 0`) are never encoded;
    /// their constants are harmless placeholders.
    pub(crate) fn new(freq: u32, cum: u32) -> Self {
        // Same renorm bound the division form uses: at most one 16-bit word per
        // step given `RANS_L`'s 16-bit renormalization.
        let x_max = (u64::from(RANS_L >> SCALE_BITS) << 16) * u64::from(freq);
        if freq < 2 {
            // freq == 1 (or 0): the reciprocal degenerates; ryg's bias trick
            // recovers the exact map, `q = mulhi(x, ~0)` giving `x - 1`.
            EncSym {
                x_max,
                rcp_freq: !0,
                rcp_shift: 0,
                bias: cum + (1 << SCALE_BITS) - 1,
                cmpl_freq: (1u32 << SCALE_BITS).wrapping_sub(freq),
            }
        } else {
            // Smallest `shift` with `2^shift >= freq`.
            let mut shift = 0u32;
            while freq > (1u32 << shift) {
                shift += 1;
            }
            let rcp_freq = (1u64 << (shift + 31)).div_ceil(u64::from(freq)) as u32;
            EncSym {
                x_max,
                rcp_freq,
                rcp_shift: shift - 1,
                bias: cum,
                cmpl_freq: (1u32 << SCALE_BITS) - freq,
            }
        }
    }

    /// Apply the rANS symbol map to a renormalized state `v` (`v < 2^32`),
    /// division-free. Equivalent to `(v / f) << SCALE_BITS + (v % f) + c`.
    #[inline]
    pub(crate) fn apply(&self, v: u32) -> u32 {
        let q = (((u64::from(v) * u64::from(self.rcp_freq)) >> 32) as u32) >> self.rcp_shift;
        v.wrapping_add(self.bias)
            .wrapping_add(q.wrapping_mul(self.cmpl_freq))
    }
}

impl Model {
    /// Build a model from raw byte counts.
    pub(crate) fn from_counts(counts: &[u32; 256]) -> Self {
        Self::from_freqs(normalize(counts))
    }

    /// Per-symbol division-free encode constants for this model.
    pub(crate) fn enc_table(&self) -> [EncSym; 256] {
        let mut table = [EncSym::new(0, 0); 256];
        for s in 0..256 {
            table[s] = EncSym::new(u32::from(self.freq[s]), u32::from(self.cum[s]));
        }
        table
    }

    /// Build a model from an already-normalized frequency table (decoder path).
    pub(crate) fn from_freqs(freq: [u16; 256]) -> Self {
        let mut cum = [0u16; 257];
        for s in 0..256 {
            cum[s + 1] = cum[s] + freq[s];
        }
        let mut slot2sym = vec![0u8; TOTFREQ as usize];
        for s in 0..256 {
            slot2sym[cum[s] as usize..cum[s + 1] as usize].fill(s as u8);
        }
        Model {
            freq,
            cum,
            slot2sym,
        }
    }
}

/// Normalize byte counts to frequencies summing to exactly [`TOTFREQ`].
///
/// Every symbol that occurs gets a frequency of at least 1; symbols that never
/// occur stay at 0. Rounding error is absorbed by the largest bucket(s), never
/// driving an occurring symbol below 1.
fn normalize(counts: &[u32; 256]) -> [u16; 256] {
    let total: u64 = counts.iter().map(|&c| u64::from(c)).sum();
    let mut freq = [0u16; 256];
    if total == 0 {
        return freq;
    }

    let mut sum: i64 = 0;
    let (mut max_s, mut max_f) = (0usize, 0u32);
    for s in 0..256 {
        if counts[s] == 0 {
            continue;
        }
        let mut f = ((u64::from(counts[s]) * u64::from(TOTFREQ)) / total) as u32;
        if f == 0 {
            f = 1;
        }
        freq[s] = f as u16;
        sum += i64::from(f);
        if f > max_f {
            max_f = f;
            max_s = s;
        }
    }

    // Absorb the rounding residual so the table sums to exactly TOTFREQ.
    let mut diff = i64::from(TOTFREQ) - sum;
    while diff != 0 {
        if diff > 0 {
            // Growing the largest bucket is always safe.
            freq[max_s] += diff as u16;
            diff = 0;
        } else {
            // Shrink the largest bucket, but never below 1; spill to the next.
            let take = (i64::from(freq[max_s]) - 1).clamp(0, -diff);
            freq[max_s] -= take as u16;
            diff += take;
            if diff != 0 {
                let (mut nm, mut nf) = (usize::MAX, 1i64);
                for s in 0..256 {
                    if i64::from(freq[s]) > nf {
                        nf = i64::from(freq[s]);
                        nm = s;
                    }
                }
                debug_assert!(nm != usize::MAX, "no bucket > 1 to absorb residual");
                max_s = nm;
            }
        }
    }
    freq
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_sum(counts: &[u32; 256]) {
        let f = normalize(counts);
        let sum: u32 = f.iter().map(|&x| u32::from(x)).sum();
        assert_eq!(sum, TOTFREQ, "freqs must sum to TOTFREQ");
        for s in 0..256 {
            if counts[s] > 0 {
                assert!(f[s] >= 1, "occurring symbol {s} dropped to 0");
            } else {
                assert_eq!(f[s], 0, "absent symbol {s} got frequency");
            }
        }
    }

    #[test]
    fn single_symbol() {
        let mut c = [0u32; 256];
        c[65] = 100;
        check_sum(&c);
        assert_eq!(normalize(&c)[65], TOTFREQ as u16);
    }

    #[test]
    fn uniform_alphabet() {
        let c = [7u32; 256];
        check_sum(&c);
    }

    #[test]
    fn highly_skewed() {
        let mut c = [1u32; 256];
        c[0] = 10_000_000;
        check_sum(&c);
    }

    #[test]
    fn two_symbols() {
        let mut c = [0u32; 256];
        c[0] = 3;
        c[255] = 1;
        check_sum(&c);
    }

    #[test]
    fn reciprocal_matches_division() {
        // The reciprocal encode map must equal the DIV/MOD form for every valid
        // frequency and every renormalized state `v < x_max` — that equivalence
        // is what keeps output byte-identical across the optimization.
        for freq in [
            1u32, 2, 3, 4, 7, 15, 16, 17, 255, 256, 1000, 2048, 4095, 4096,
        ] {
            let cum = 123u32.min(TOTFREQ - freq); // any legal start
            let sym = EncSym::new(freq, cum);
            let x_max = sym.x_max;
            // Sample states across the renorm interval [x_max>>16, x_max).
            let lo = (x_max >> 16).max(1);
            let step = ((x_max - lo) / 997).max(1);
            let mut v = lo;
            while v < x_max {
                let vv = v as u32;
                let expect = ((vv / freq) << SCALE_BITS) + (vv % freq) + cum;
                assert_eq!(sym.apply(vv), expect, "freq={freq} v={vv}");
                v += step;
            }
        }
    }
}
