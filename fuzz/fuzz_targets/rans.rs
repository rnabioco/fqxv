#![no_main]
//! Fuzz the rANS entropy decoder. It reads an output-length header and a
//! frequency table straight from the stream, so it must reject a malformed
//! table or an absurd length rather than abort.

use libfuzzer_sys::fuzz_target;

/// Cap the fuzzed output size. The stream's second field is a `u64` output length
/// (`[u8 order][u64 n]…`) that the standalone decoder trusts — a tiny stream can
/// declare a multi-GB output (a decompression bomb the *container* bounds
/// structurally, not this raw API). That is resource exhaustion, not the
/// panic/abort robustness this target checks, so skip it: logic bugs still fire at
/// small `n`.
const MAX_OUTPUT: u64 = 64 << 20; // 64 MiB

fuzz_target!(|data: &[u8]| {
    if let Some(len) = data.get(1..9) {
        if u64::from_le_bytes(len.try_into().unwrap()) > MAX_OUTPUT {
            return;
        }
    }
    let _ = fqxv_rans::decode(data);
});
