//! rANS Nx16 entropy coder — a clean-room implementation of the CRAM 3.1 codec.
//!
//! The coder maintains 32 interleaved rANS states (16-bit renormalization) so
//! that renormalization and symbol lookup vectorize cleanly. Order-0 and
//! order-1 models are supported.
//!
//! Decode backends live behind one API and are chosen at runtime via
//! `is_x86_feature_detected!`; the decoded output is identical whichever runs:
//!
//! - **scalar** — always available, the correctness reference (all orders).
//! - **AVX2** — order-0 decode on x86-64.
//!
//! An SSE4.2 backend and vectorized order-1/encode are future work; [`Backend`]
//! already reports the detected tier.
//!
//! Implemented from the CRAM codecs specification
//! (<https://samtools.github.io/hts-specs/CRAMcodecs.pdf>); see
//! `THIRD-PARTY-NOTICES.md`.

#![doc(html_root_url = "https://docs.rs/fqxv-rans")]

#[cfg(target_arch = "x86_64")]
mod avx2;
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
        Order::One => Ok(scalar::encode_order1(src)),
    }
}

/// Internal entry points exposed only for benchmarking (the `bench` feature).
/// Not part of the public API and may change at any time.
#[cfg(feature = "bench")]
#[doc(hidden)]
pub mod bench_api {
    use crate::Result;

    /// Decode an order-0/order-1 stream with the portable scalar backend.
    pub fn decode_scalar(src: &[u8]) -> Result<Vec<u8>> {
        crate::scalar::decode(src)
    }

    /// Decode an order-0 stream with the AVX2 backend.
    #[cfg(target_arch = "x86_64")]
    pub fn decode_avx2(src: &[u8]) -> Result<Vec<u8>> {
        crate::avx2::decode_order0(src)
    }
}

/// Decode a stream produced by [`encode`].
///
/// Defaults to the scalar backend. On AMD Zen 3 (this cluster's EPYC 7713P) the
/// AVX2 gather-based order-0 decoder measured *slower* than scalar (≈208 vs
/// ≈250 MiB/s) because `vpgatherdd` is microcoded — the same reason htscodecs
/// falls back from gather-heavy SIMD. We therefore dispatch on measured, not
/// assumed, speed: scalar is the default, and the AVX2 path stays available via
/// [`bench_api`] / [`Backend`] for micro-architectures with fast gather and for
/// a future gather-free SIMD design. Order-0 output is byte-identical either way.
pub fn decode(src: &[u8]) -> Result<Vec<u8>> {
    scalar::decode(src)
}

/// Decode using an explicitly chosen [`Backend`].
///
/// [`Backend::Avx2`] uses the AVX2 order-0 decoder when the stream is order-0
/// and the CPU supports AVX2; it falls back to scalar otherwise (including for
/// order-1). Provided for benchmarking and for callers who know their
/// micro-architecture has fast gather — see [`decode`] for why scalar is the
/// default. Output is identical regardless of backend.
pub fn decode_with(src: &[u8], backend: Backend) -> Result<Vec<u8>> {
    #[cfg(target_arch = "x86_64")]
    {
        if backend == Backend::Avx2
            && src.first() == Some(&0)
            && std::is_x86_feature_detected!("avx2")
        {
            return avx2::decode_order0(src);
        }
    }
    let _ = backend;
    scalar::decode(src)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backend_detect_is_available() {
        // Whatever we detect must be a real variant; smoke test only.
        let b = Backend::detect();
        assert!(matches!(
            b,
            Backend::Scalar | Backend::Sse42 | Backend::Avx2
        ));
    }

    fn roundtrip(src: &[u8]) {
        for order in [Order::Zero, Order::One] {
            let enc = encode(src, order).expect("encode");
            let dec = decode(&enc).expect("decode");
            assert_eq!(
                dec,
                src,
                "round-trip mismatch (order {order:?}, len {})",
                src.len()
            );
        }
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
            for order in [Order::Zero, Order::One] {
                let enc = encode(&data, order).expect("encode");
                let dec = decode(&enc).expect("decode");
                proptest::prop_assert_eq!(&dec, &data);
            }
        }

        #[test]
        fn roundtrip_arbitrary_order1(data: Vec<u8>) {
            let enc = encode(&data, Order::One).expect("encode");
            let dec = decode(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, data);
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_matches_scalar_order0() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        // Sizes straddling the 32-state round boundary, plus large multi-round.
        for n in [0usize, 1, 5, 31, 32, 33, 64, 65, 1000, 100_003] {
            let data: Vec<u8> = (0..n).map(|i| (i * 131 % 251) as u8).collect();
            let enc = encode(&data, Order::Zero).expect("encode");
            let scalar_dec = crate::scalar::decode(&enc).expect("scalar decode");
            let avx2_dec = crate::avx2::decode_order0(&enc).expect("avx2 decode");
            assert_eq!(scalar_dec, data, "scalar mismatch at n={n}");
            assert_eq!(avx2_dec, data, "avx2 mismatch at n={n}");
        }
    }

    #[test]
    fn order1_beats_order0_on_correlated_data() {
        // Runs of correlated bytes: order-1 should model the transitions better.
        let mut data = Vec::new();
        let mut x = 0u8;
        for i in 0..100_000u32 {
            x = x.wrapping_add((i % 4) as u8);
            data.push(b'A' + (x % 4));
        }
        let o0 = encode(&data, Order::Zero).expect("o0").len();
        let o1 = encode(&data, Order::One).expect("o1").len();
        assert!(o1 < o0, "order-1 ({o1}) should beat order-0 ({o0}) here");
    }
}
