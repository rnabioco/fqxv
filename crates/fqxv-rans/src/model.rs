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

impl Model {
    /// Build a model from raw byte counts.
    pub(crate) fn from_counts(counts: &[u32; 256]) -> Self {
        Self::from_freqs(normalize(counts))
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
}
