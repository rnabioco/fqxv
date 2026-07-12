//! rANS Nx16 entropy coder — a clean-room implementation of the CRAM 3.1 codec.
//!
//! The coder maintains 32 interleaved rANS states so that renormalization and
//! symbol lookup vectorize cleanly. Three backends live behind one API:
//!
//! - [`Backend::Scalar`] — always available, the correctness reference.
//! - [`Backend::Sse42`] — x86-64 SSE4.2.
//! - [`Backend::Avx2`] — x86-64 AVX2.
//!
//! The best available backend is chosen at runtime via
//! `is_x86_feature_detected!`; encoded output is identical across backends
//! (only throughput differs).
//!
//! Implemented from the CRAM codecs specification
//! (<https://samtools.github.io/hts-specs/CRAMcodecs.pdf>); see
//! `THIRD-PARTY-NOTICES.md`.
//!
//! Status: **scaffold.** The scalar reference is implemented next (M1).

#![doc(html_root_url = "https://docs.rs/fqxv-rans")]

mod model;
mod scalar;

use thiserror::Error;

/// Errors returned by the rANS codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The input was longer than the coder supports in a single block.
    #[error("input too large for a single rANS block: {0} bytes")]
    InputTooLarge(usize),
    /// The compressed stream was malformed or truncated.
    #[error("malformed rANS stream: {0}")]
    Malformed(&'static str),
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// The context order of the frequency model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    /// Order-0: symbol frequencies are independent of context.
    Zero,
    /// Order-1: frequencies are conditioned on the previous byte.
    One,
}

/// Which vector backend performs the coding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    /// Portable scalar reference (always available).
    Scalar,
    /// x86-64 SSE4.2.
    Sse42,
    /// x86-64 AVX2.
    Avx2,
}

impl Backend {
    /// Select the fastest backend supported by the running CPU.
    #[must_use]
    pub fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx2") {
                return Backend::Avx2;
            }
            if std::is_x86_feature_detected!("sse4.2") {
                return Backend::Sse42;
            }
        }
        Backend::Scalar
    }
}

/// Encode `src` with the given model order.
///
/// The output is byte-identical regardless of backend; today all orders run on
/// the scalar coder while the SIMD backends are brought up (M1).
pub fn encode(src: &[u8], order: Order) -> Result<Vec<u8>> {
    match order {
        Order::Zero => Ok(scalar::encode_order0(src)),
        Order::One => Err(Error::NotImplemented("rANS order-1 encode (M1 follow-up)")),
    }
}

/// Decode a stream produced by [`encode`].
pub fn decode(src: &[u8]) -> Result<Vec<u8>> {
    scalar::decode(src)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backend_detect_is_available() {
        // Whatever we detect must be a real variant; smoke test only.
        let b = Backend::detect();
        assert!(matches!(b, Backend::Scalar | Backend::Sse42 | Backend::Avx2));
    }

    fn roundtrip(src: &[u8]) {
        let enc = encode(src, Order::Zero).expect("encode");
        let dec = decode(&enc).expect("decode");
        assert_eq!(dec, src, "round-trip mismatch (len {})", src.len());
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(b"");
    }

    #[test]
    fn roundtrip_single_byte() {
        roundtrip(b"Q");
    }

    #[test]
    fn roundtrip_all_same() {
        roundtrip(&[b'A'; 1000]);
    }

    #[test]
    fn roundtrip_short_and_odd_lengths() {
        // Lengths straddling the 32-state interleave boundary.
        for n in [1usize, 2, 31, 32, 33, 63, 64, 65, 100] {
            let data: Vec<u8> = (0..n).map(|i| (i * 7 + 3) as u8).collect();
            roundtrip(&data);
        }
    }

    #[test]
    fn roundtrip_full_alphabet() {
        let data: Vec<u8> = (0..=255u8).cycle().take(10_000).collect();
        roundtrip(&data);
    }

    #[test]
    fn compresses_skewed_data() {
        // Quality-like data with a small skewed alphabet should shrink well.
        let mut data = Vec::new();
        for i in 0..50_000u32 {
            data.push(b'I' - (i % 8) as u8 * (i % 3 == 0) as u8);
        }
        let enc = encode(&data, Order::Zero).expect("encode");
        assert!(
            enc.len() < data.len() / 2,
            "expected >2x on skewed data, got {} -> {}",
            data.len(),
            enc.len()
        );
        assert_eq!(decode(&enc).expect("decode"), data);
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(data: Vec<u8>) {
            let enc = encode(&data, Order::Zero).expect("encode");
            let dec = decode(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, data);
        }

        #[test]
        fn roundtrip_small_alphabet(
            data in proptest::collection::vec(0u8..6, 0..2000)
        ) {
            let enc = encode(&data, Order::Zero).expect("encode");
            let dec = decode(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, data);
        }
    }
}
