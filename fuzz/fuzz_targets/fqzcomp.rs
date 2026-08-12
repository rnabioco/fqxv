#![no_main]
//! Fuzz the quality-score (fqzcomp) decoder.

use libfuzzer_sys::fuzz_target;

/// A constant sequence long enough to satisfy any plausible fuzzed length
/// array, so the sequence-conditioned modes are reachable (below).
static SEQ: [u8; 1 << 16] = [b'A'; 1 << 16];

fuzz_target!(|data: &[u8]| {
    let _ = fqxv_fqzcomp::decode(data);
    // The sequence-conditioned modes (3/4 and the chunked 5/6) reject any
    // decode whose sequence length differs from the declared total, so the
    // bare `decode` above never reaches their payload paths — including the
    // mode-5/6 chunk-table parse and segment fan. Offer a few fixed sequence
    // lengths; the fuzzer learns length arrays summing to one of them.
    for n in [64usize, 4096, SEQ.len()] {
        let _ = fqxv_fqzcomp::decode_seq(data, &SEQ[..n]);
    }
});
