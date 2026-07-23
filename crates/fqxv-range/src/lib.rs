//! Binary range coder and adaptive frequency models.
//!
//! This is the serial entropy backend for the fqzcomp quality model
//! ([`fqxv-fqzcomp`](https://docs.rs/fqxv-fqzcomp)). The range coder is the
//! classic Subbotin carryless design (public domain) used by the CRAM codecs;
//! [`SimpleModel`] is an adaptive frequency table over a small alphabet.
//!
//! ```
//! use fqxv_range::{Encoder, Decoder, SimpleModel};
//! let data = [3usize, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
//! let mut enc = Encoder::new();
//! let mut m = SimpleModel::<10>::new();
//! for &s in &data { m.encode(&mut enc, s); }
//! let bytes = enc.finish();
//!
//! let mut dec = Decoder::new(&bytes);
//! let mut m = SimpleModel::<10>::new();
//! let out: Vec<usize> = (0..data.len()).map(|_| m.decode(&mut dec)).collect();
//! assert_eq!(out, data);
//! ```

use thiserror::Error;

/// Errors returned by the range coder.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed range-coded stream: {0}")]
    Malformed(&'static str),
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

// Subbotin carryless range coder constants.
const TOP: u32 = 1 << 24;
const BOT: u32 = 1 << 16;

/// Range encoder. Feed cumulative-frequency intervals via [`Encoder::encode`],
/// then [`Encoder::finish`] to flush the tail.
#[derive(Debug)]
pub struct Encoder {
    low: u32,
    range: u32,
    out: Vec<u8>,
}

impl Default for Encoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Encoder {
    /// Create an empty encoder.
    #[must_use]
    pub fn new() -> Self {
        Encoder {
            low: 0,
            range: u32::MAX,
            out: Vec::new(),
        }
    }

    /// Encode a symbol occupying `[cum, cum + freq)` out of `tot`.
    ///
    /// Requires `freq > 0`, `cum + freq <= tot`, and `tot < BOT` (2^16).
    #[inline]
    pub fn encode(&mut self, cum: u32, freq: u32, tot: u32) {
        self.range /= tot;
        self.low = self.low.wrapping_add(cum * self.range);
        self.range *= freq;
        self.renorm();
    }

    #[inline]
    fn renorm(&mut self) {
        while (self.low ^ self.low.wrapping_add(self.range)) < TOP
            || (self.range < BOT && {
                self.range = self.low.wrapping_neg() & (BOT - 1);
                true
            })
        {
            self.out.push((self.low >> 24) as u8);
            self.low <<= 8;
            self.range <<= 8;
        }
    }

    /// Flush the coder and return the compressed bytes.
    #[must_use]
    pub fn finish(mut self) -> Vec<u8> {
        for _ in 0..4 {
            self.out.push((self.low >> 24) as u8);
            self.low <<= 8;
        }
        self.out
    }
}

/// Range decoder, the inverse of [`Encoder`]. Use [`Decoder::freq`] to find the
/// target within `tot`, look up the symbol, then [`Decoder::decode`] to consume
/// its interval.
#[derive(Debug)]
pub struct Decoder<'a> {
    low: u32,
    range: u32,
    code: u32,
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Decoder<'a> {
    /// Create a decoder over `buf` (the output of [`Encoder::finish`]).
    #[must_use]
    pub fn new(buf: &'a [u8]) -> Self {
        let mut d = Decoder {
            low: 0,
            range: u32::MAX,
            code: 0,
            buf,
            pos: 0,
        };
        for _ in 0..4 {
            let b = d.next_byte();
            d.code = (d.code << 8) | u32::from(b);
        }
        d
    }

    #[inline]
    fn next_byte(&mut self) -> u8 {
        let b = self.buf.get(self.pos).copied().unwrap_or(0);
        self.pos += 1;
        b
    }

    /// Return the target frequency in `[0, tot)` for the next symbol. The caller
    /// maps it to a symbol interval `[cum, cum + freq)` and calls [`decode`].
    ///
    /// [`decode`]: Decoder::decode
    #[inline]
    pub fn freq(&mut self, tot: u32) -> u32 {
        self.range /= tot;
        self.code.wrapping_sub(self.low) / self.range
    }

    /// Consume the interval `[cum, cum + freq)` identified from [`freq`].
    ///
    /// [`freq`]: Decoder::freq
    #[inline]
    pub fn decode(&mut self, cum: u32, freq: u32) {
        self.low = self.low.wrapping_add(cum * self.range);
        self.range *= freq;
        self.renorm();
    }

    #[inline]
    fn renorm(&mut self) {
        while (self.low ^ self.low.wrapping_add(self.range)) < TOP
            || (self.range < BOT && {
                self.range = self.low.wrapping_neg() & (BOT - 1);
                true
            })
        {
            let b = self.next_byte();
            self.code = (self.code << 8) | u32::from(b);
            self.low <<= 8;
            self.range <<= 8;
        }
    }
}

/// An adaptive order-0 frequency model over `N` symbols, driving the range
/// coder. Frequencies start uniform and grow toward observed symbols, with
/// periodic halving to bound the total and keep the model adaptive.
#[derive(Debug, Clone)]
pub struct SimpleModel<const N: usize> {
    freq: [u16; N],
    tot: u32,
}

impl<const N: usize> Default for SimpleModel<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> SimpleModel<N> {
    /// Total frequency cap; kept well below `BOT` so `range / tot` stays valid.
    const MAX_TOT: u32 = 1 << 13;
    /// Frequency increment applied to the observed symbol.
    const STEP: u16 = 16;

    /// Create a model with uniform frequencies over all `N` symbols.
    #[must_use]
    pub fn new() -> Self {
        SimpleModel {
            freq: [1; N],
            tot: N as u32,
        }
    }

    /// Create a model over only the first `active` symbols (`1..=N`); the rest
    /// carry zero frequency and are never coded. Use this when the true alphabet
    /// is smaller than the compile-time capacity `N`: the unused ("phantom")
    /// slots otherwise sit at frequency 1 forever, permanently taxing every coded
    /// symbol — on a skewed 4-symbol quality stream coded with `N = 94`, those 90
    /// phantoms cap the dominant symbol's probability and cost several percent.
    ///
    /// # Panics
    /// If `active` is 0 or greater than `N`.
    #[must_use]
    pub fn with_active(active: usize) -> Self {
        assert!(active >= 1 && active <= N, "active {active} out of 1..={N}");
        let mut freq = [0u16; N];
        for f in &mut freq[..active] {
            *f = 1;
        }
        SimpleModel {
            freq,
            tot: active as u32,
        }
    }

    /// Encode symbol `sym` (`0..N`) and adapt.
    pub fn encode(&mut self, enc: &mut Encoder, sym: usize) {
        let mut cum = 0u32;
        for f in &self.freq[..sym] {
            cum += u32::from(*f);
        }
        enc.encode(cum, u32::from(self.freq[sym]), self.tot);
        self.update(sym);
    }

    /// Decode and adapt, returning the symbol (`0..N`).
    pub fn decode(&mut self, dec: &mut Decoder<'_>) -> usize {
        let target = dec.freq(self.tot);
        let mut cum = 0u32;
        let mut sym = 0;
        while sym + 1 < N && cum + u32::from(self.freq[sym]) <= target {
            cum += u32::from(self.freq[sym]);
            sym += 1;
        }
        dec.decode(cum, u32::from(self.freq[sym]));
        self.update(sym);
        sym
    }

    #[inline]
    fn update(&mut self, sym: usize) {
        self.freq[sym] += Self::STEP;
        self.tot += u32::from(Self::STEP);
        if self.tot > Self::MAX_TOT {
            let mut t = 0u32;
            for f in &mut self.freq {
                *f = (*f + 1) >> 1; // halve, keep >= 1
                t += u32::from(*f);
            }
            self.tot = t;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip<const N: usize>(data: &[usize]) {
        let mut enc = Encoder::new();
        let mut m = SimpleModel::<N>::new();
        for &s in data {
            assert!(s < N);
            m.encode(&mut enc, s);
        }
        let bytes = enc.finish();

        let mut dec = Decoder::new(&bytes);
        let mut m = SimpleModel::<N>::new();
        let out: Vec<usize> = (0..data.len()).map(|_| m.decode(&mut dec)).collect();
        assert_eq!(out, data, "range-coder round-trip mismatch");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip::<4>(&[]);
    }

    #[test]
    fn roundtrip_single() {
        roundtrip::<8>(&[5]);
    }

    #[test]
    fn roundtrip_binary_skewed() {
        // Mostly zeros: should compress and round-trip.
        let mut data = vec![0usize; 10_000];
        for i in (0..data.len()).step_by(97) {
            data[i] = 1;
        }
        roundtrip::<2>(&data);
    }

    #[test]
    fn skewed_compresses() {
        let data = vec![0usize; 100_000];
        let mut enc = Encoder::new();
        let mut m = SimpleModel::<4>::new();
        for &s in &data {
            m.encode(&mut enc, s);
        }
        let bytes = enc.finish();
        // 100k near-deterministic symbols should be tiny.
        assert!(
            bytes.len() < 200,
            "expected strong compression, got {}",
            bytes.len()
        );
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(data in proptest::collection::vec(0usize..64, 0..5000)) {
            roundtrip::<64>(&data);
        }

        #[test]
        fn roundtrip_full_alphabet(data in proptest::collection::vec(0usize..256, 0..3000)) {
            roundtrip::<256>(&data);
        }

        /// Arbitrary bytes — never a valid encoding — must never panic or abort
        /// the decoder; a `Decoder` built over garbage and driven through a
        /// bounded number of `SimpleModel` decodes only ever yields (meaningless)
        /// in-range symbols. `next_byte` reads past the buffer as zeros, so there
        /// is no length header and no allocation to bound; the guard is that the
        /// coder never overruns, divides by zero, or indexes out of range.
        #[test]
        fn decode_never_aborts_on_garbage(bytes in proptest::collection::vec(0u8..=255, 0..256)) {
            let mut dec = Decoder::new(&bytes);
            let mut m = SimpleModel::<64>::new();
            for _ in 0..1000 {
                let sym = m.decode(&mut dec);
                proptest::prop_assert!(sym < 64);
            }
        }
    }
}
