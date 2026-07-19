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
//! - **AVX2** — order-0 decode and encode (8-wide). One L1-resident gather per
//!   group plus a branchless renorm; encode is gather-bound (see [`encode`]).
//! - **AVX-512** — order-0 decode and encode (16-wide, two groups per round).
//!   Adds native `k`-mask compares and a `vpexpandd` renorm; the widest
//!   detected path wins, so it supersedes AVX2 on Intel (decode ≈5×, encode
//!   ≈1.4× scalar on Cascade Lake).
//!
//! There is deliberately no SSE4.2 vector path: gather (`vpgatherdd`) requires
//! AVX2, and without it a 4-lane path degenerates to per-lane scalar loads with
//! no advantage over the scalar coder. Order-1 stays scalar (its context is the
//! previous decoded symbol — a cross-lane serial dependency). [`Backend`]
//! reports the detected tier; [`Backend::Sse42`] runs the scalar path.
//!
//! Implemented from the CRAM codecs specification
//! (<https://samtools.github.io/hts-specs/CRAMcodecs.pdf>); see
//! `THIRD-PARTY-NOTICES.md`.

#![doc(html_root_url = "https://docs.rs/fqxv-rans")]

#[cfg(target_arch = "x86_64")]
mod avx2;
#[cfg(target_arch = "x86_64")]
mod avx512;
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
    /// x86-64 AVX-512 (F).
    Avx512,
}

impl Backend {
    /// Select the fastest backend supported by the running CPU.
    #[must_use]
    pub fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512f") {
                return Backend::Avx512;
            }
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
/// The output is byte-identical regardless of backend — a stream written on an
/// AVX2 host decodes bit-for-bit the same as one written on a scalar host, and
/// two hosts encode the same input to the same bytes.
///
/// Order-1 is always scalar. Order-0 picks the widest vector encoder that wins
/// on the host: AVX-512 where present (≈233 vs ≈163 MiB/s scalar on Cascade
/// Lake — the 16-wide gather turns the AVX2 wash into a win), otherwise the AVX2
/// encoder on Broadwell-class cores (≈206 vs ≈137). Both are gather-bound, so an
/// AVX-512-class core that also has AVX2 skips AVX2 (which only ties there) in
/// favor of AVX-512. Output is byte-identical to scalar whichever path runs.
pub fn encode(src: &[u8], order: Order) -> Result<Vec<u8>> {
    match order {
        Order::Zero => {
            #[cfg(target_arch = "x86_64")]
            {
                if std::is_x86_feature_detected!("avx512f") {
                    return Ok(avx512::encode_order0(src));
                }
                if std::is_x86_feature_detected!("avx2") {
                    return Ok(avx2::encode_order0(src));
                }
            }
            Ok(scalar::encode_order0(src))
        }
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

    /// Encode an order-0 stream with the portable scalar backend.
    pub fn encode_order0_scalar(src: &[u8]) -> Vec<u8> {
        crate::scalar::encode_order0(src)
    }

    /// Encode an order-0 stream with the AVX2 backend.
    #[cfg(target_arch = "x86_64")]
    pub fn encode_order0_avx2(src: &[u8]) -> Vec<u8> {
        crate::avx2::encode_order0(src)
    }

    /// Decode an order-0 stream with the AVX-512 backend.
    #[cfg(target_arch = "x86_64")]
    pub fn decode_avx512(src: &[u8]) -> Result<Vec<u8>> {
        crate::avx512::decode_order0(src)
    }

    /// Encode an order-0 stream with the AVX-512 backend.
    #[cfg(target_arch = "x86_64")]
    pub fn encode_order0_avx512(src: &[u8]) -> Vec<u8> {
        crate::avx512::encode_order0(src)
    }
}

/// Decode a stream produced by [`encode`].
///
/// Order-0 streams take the widest vector decoder the CPU supports — AVX-512,
/// else AVX2; order-1 (and non-x86) fall back to scalar. Both vector decoders
/// resolve each group with one L1-resident gather and a branchless renorm.
/// Measured on Intel Cascade Lake (8 MiB order-0): scalar ≈163, AVX2 ≈315,
/// AVX-512 ≈817 MiB/s (5×) — AVX-512 adds native `k`-mask compares and a
/// single `vpexpandd` for the renorm distribution. Output is byte-identical to
/// scalar whichever path runs.
///
/// Historical note: an earlier three-gather AVX2 decoder ran *slower* than
/// scalar (microcoded `vpgatherdd`, worst on AMD Zen 3), which is why scalar was
/// once the default. The single-gather rewrite reverses that on Intel; on Zen it
/// is strictly less gather-bound than the version that regressed, but re-measure
/// there before assuming a win. [`decode_with`] forces a specific [`Backend`].
pub fn decode(src: &[u8]) -> Result<Vec<u8>> {
    #[cfg(target_arch = "x86_64")]
    {
        // Order tag 0 == order-0; the only vectorized path today.
        if src.first() == Some(&0) {
            if std::is_x86_feature_detected!("avx512f") {
                return avx512::decode_order0(src);
            }
            if std::is_x86_feature_detected!("avx2") {
                return avx2::decode_order0(src);
            }
        }
    }
    scalar::decode(src)
}

/// The declared output length in a stream header, without decoding.
///
/// The header is `[u8 order][u64 n]` little-endian (see [`scalar`]); `n` is the
/// byte count the stream reconstructs to. Returns `None` if `src` is too short
/// to hold the header.
fn peek_output_len(src: &[u8]) -> Option<u64> {
    let bytes = src.get(1..9)?;
    Some(u64::from_le_bytes(bytes.try_into().unwrap()))
}

/// Decode `src`, rejecting a declared output length above `max_len`.
///
/// The output length is a `u64` in the stream header. A corrupt or hostile value
/// can be far larger than the stream could ever legitimately reconstruct: a
/// single-symbol model emits its entire output from the header alone, consuming
/// no renorm input, so there is no input-relative bound that separates a valid
/// high-compression stream from a decompression bomb. `try_reserve` does not stop
/// the allocation under memory overcommit — a multi-gigabyte `n` reserves fine
/// and then stalls or is OOM-killed while the buffer is zeroed.
///
/// Callers decoding untrusted bytes therefore pass the largest length they could
/// legitimately expect and this rejects anything beyond it *before* allocating.
/// [`decode`] imposes no bound and is only safe on data whose integrity is
/// already established (e.g. a CRC-checked container block).
pub fn decode_bounded(src: &[u8], max_len: usize) -> Result<Vec<u8>> {
    if peek_output_len(src).is_some_and(|n| n > max_len as u64) {
        return Err(Error::Malformed("rANS output length exceeds caller bound"));
    }
    decode(src)
}

/// Decode using an explicitly chosen [`Backend`].
///
/// [`Backend::Avx512`] / [`Backend::Avx2`] use the corresponding order-0 vector
/// decoder when the stream is order-0 and the CPU supports that feature; any
/// other case (order-1, unsupported CPU, [`Backend::Sse42`], [`Backend::Scalar`])
/// falls back to scalar. Provided for benchmarking and for pinning a backend;
/// [`decode`] already selects the widest supported path. Output is identical
/// regardless of backend.
pub fn decode_with(src: &[u8], backend: Backend) -> Result<Vec<u8>> {
    #[cfg(target_arch = "x86_64")]
    {
        if src.first() == Some(&0) {
            if backend == Backend::Avx512 && std::is_x86_feature_detected!("avx512f") {
                return avx512::decode_order0(src);
            }
            if backend == Backend::Avx2 && std::is_x86_feature_detected!("avx2") {
                return avx2::decode_order0(src);
            }
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
            Backend::Scalar | Backend::Sse42 | Backend::Avx2 | Backend::Avx512
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
    fn decode_bounded_rejects_inflated_length() {
        // A single-symbol stream emits its whole output from the `n` header with
        // no input to consume, so an inflated `n` is a decompression bomb. The
        // dangerous band is a length that `try_reserve` accepts under overcommit
        // (so `decode` would stall/OOM) but that exceeds the caller's bound.
        let enc = encode(&[b'A'; 64], Order::Zero).unwrap();
        let mut bomb = enc.clone();
        bomb[1..9].copy_from_slice(&(1u64 << 31).to_le_bytes()); // 2 GiB
        assert!(matches!(
            decode_bounded(&bomb, 1 << 20),
            Err(Error::Malformed(_))
        ));
        // The unmodified stream decodes fine when its true length is within bound.
        assert_eq!(decode_bounded(&enc, 1 << 20).unwrap(), vec![b'A'; 64]);
        // A header too short to hold the length field still errors (in decode),
        // rather than being mistaken for an unbounded stream.
        assert!(decode_bounded(&enc[..1], 1 << 20).is_err());
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

        // The AVX2 order-0 decoder must be byte-identical to scalar for every
        // input. A small skewed alphabet over multi-round lengths stresses the
        // vectorized renorm compaction (varied per-group renorm masks) and the
        // 32-state round boundary far harder than the fixed-size smoke test.
        #[cfg(target_arch = "x86_64")]
        #[test]
        fn avx2_matches_scalar_multiround(
            data in proptest::collection::vec(0u8..40, 0..5000)
        ) {
            if !std::is_x86_feature_detected!("avx2") {
                return Ok(());
            }
            let enc = encode(&data, Order::Zero).expect("encode");
            let scalar_dec = crate::scalar::decode(&enc).expect("scalar decode");
            let avx2_dec = crate::avx2::decode_order0(&enc).expect("avx2 decode");
            proptest::prop_assert_eq!(&avx2_dec, &scalar_dec);
            proptest::prop_assert_eq!(&avx2_dec, &data);
        }

        // The AVX2 order-0 *encoder* must emit the exact same bytes as scalar —
        // otherwise the same input compresses differently on different hosts.
        #[cfg(target_arch = "x86_64")]
        #[test]
        fn avx2_encode_matches_scalar(
            data in proptest::collection::vec(0u8..40, 0..5000)
        ) {
            if !std::is_x86_feature_detected!("avx2") {
                return Ok(());
            }
            let avx2_enc = crate::avx2::encode_order0(&data);
            let scalar_enc = crate::scalar::encode_order0(&data);
            proptest::prop_assert_eq!(&avx2_enc, &scalar_enc);
            proptest::prop_assert_eq!(crate::scalar::decode(&avx2_enc).expect("decode"), data);
        }

        // The AVX-512 order-0 decode and encode must both be byte-identical to
        // scalar: decode reproduces the bytes, encode emits the same stream.
        #[cfg(target_arch = "x86_64")]
        #[test]
        fn avx512_matches_scalar(
            data in proptest::collection::vec(0u8..40, 0..5000)
        ) {
            if !std::is_x86_feature_detected!("avx512f") {
                return Ok(());
            }
            let scalar_enc = crate::scalar::encode_order0(&data);
            let avx512_enc = crate::avx512::encode_order0(&data);
            proptest::prop_assert_eq!(&avx512_enc, &scalar_enc);
            let avx512_dec = crate::avx512::decode_order0(&scalar_enc).expect("avx512 decode");
            proptest::prop_assert_eq!(&avx512_dec, &data);
        }
    }

    // Frequency regimes the branchless SIMD reciprocal must get exactly right:
    // a dominant symbol (freq > 2048 → the scalar coder's exact-division path)
    // and a single-symbol alphabet (freq == TOTFREQ, the never-renorm case).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_encode_matches_scalar_skewed() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let cases: [Vec<u8>; 3] = [
            std::iter::repeat_n(b'A', 60_000)
                .chain((0..4000).map(|i| b'B' + (i % 20) as u8))
                .collect(),
            vec![b'Q'; 40_000],
            (0..80_000u32)
                .map(|i| b'!' + (i % 3 == 0) as u8 * 7)
                .collect(),
        ];
        for data in cases {
            assert_eq!(
                crate::avx2::encode_order0(&data),
                crate::scalar::encode_order0(&data),
                "avx2/scalar encode diverged (len {})",
                data.len()
            );
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

    #[test]
    fn decode_rejects_malformed_freq_table_without_aborting() {
        // order-0 stream, n=1, with a frequency table that doesn't sum to TOTFREQ.
        // Pre-hardening `Model::from_freqs` sliced `slot2sym` out of bounds and
        // panicked (a process abort under panic=abort); now `decode` returns Err.
        let mut src = vec![0u8]; // order 0
        src.extend_from_slice(&1u64.to_le_bytes()); // n = 1
        let mut freq = [0u16; 256];
        freq[0] = crate::model::TOTFREQ as u16 + 100; // exceeds TOTFREQ
        for f in freq {
            src.extend_from_slice(&f.to_le_bytes());
        }
        assert!(matches!(decode(&src), Err(Error::Malformed(_))));
    }

    #[test]
    fn decode_rejects_huge_output_length_without_aborting() {
        // Valid single-symbol freq table but a corrupt, enormous output length.
        // The output reservation must fail cleanly instead of aborting on a
        // multi-exabyte infallible allocation.
        let mut src = vec![0u8]; // order 0
        src.extend_from_slice(&u64::MAX.to_le_bytes()); // n = huge
        let mut freq = [0u16; 256];
        freq[0] = crate::model::TOTFREQ as u16; // sums to exactly TOTFREQ
        for f in freq {
            src.extend_from_slice(&f.to_le_bytes());
        }
        // Encoder states (values irrelevant — the allocation is attempted first).
        for _ in 0..crate::model::N_STATES {
            src.extend_from_slice(&0u32.to_le_bytes());
        }
        assert!(matches!(decode(&src), Err(Error::Malformed(_))));
    }

    /// Apply byte substitutions, an optional 0xFF window (to drive a length/count
    /// field high), and an optional truncation to `data`.
    fn corrupt(
        mut data: Vec<u8>,
        subs: &[(usize, u8)],
        wipe: Option<(usize, usize)>,
        trunc: Option<usize>,
    ) -> Vec<u8> {
        for &(pos, val) in subs {
            if !data.is_empty() {
                let i = pos % data.len();
                data[i] = val;
            }
        }
        if let (Some((pos, w)), false) = (wipe, data.is_empty()) {
            let i = pos % data.len();
            for b in data.iter_mut().skip(i).take(w) {
                *b = 0xFF;
            }
        }
        if let Some(t) = trunc
            && !data.is_empty()
        {
            data.truncate(t % data.len());
        }
        data
    }

    proptest::proptest! {
        /// Encode valid data, corrupt the stream, then decode: a mutated rANS
        /// stream must never panic or abort — only return Ok/Err. Covers both
        /// orders and whichever backend `decode` dispatches to on this CPU.
        ///
        /// Uses `decode_bounded`, because that is the entry point this scenario
        /// describes. `corrupt` rewrites arbitrary offsets, including the 8-byte
        /// output length at `src[1..9]`, and plain `decode` is documented to
        /// allocate whatever that says — so calling it here tested the one API
        /// that is explicitly not for untrusted bytes, and duly allocated on
        /// them: 15.6 GB resident and a 120s CI timeout, from a seed rare enough
        /// that 20000 cases under a 4 GB cap never hit it. The bound is the only
        /// thing separating "handles corruption" from "obeys corruption".
        #[test]
        fn decode_survives_mutation(
            data in proptest::collection::vec(0u8..40, 0..2000),
            zero_order in proptest::bool::ANY,
            subs in proptest::collection::vec((proptest::num::usize::ANY, proptest::num::u8::ANY), 0..12),
            wipe in proptest::option::of((proptest::num::usize::ANY, 1usize..9)),
            trunc in proptest::option::of(proptest::num::usize::ANY),
        ) {
            let order = if zero_order { Order::Zero } else { Order::One };
            let enc = encode(&data, order).expect("encode");
            // What a caller legitimately expects back: the data we encoded.
            let _ = decode_bounded(&corrupt(enc, &subs, wipe, trunc), data.len());
        }
    }
}
