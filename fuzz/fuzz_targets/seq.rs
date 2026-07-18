#![no_main]
//! Fuzz the order-k sequence decoder (context model + exception list).

use libfuzzer_sys::fuzz_target;

// The order-k context table (`4^k` models) and hashed table (`2^hb` slots) are
// now bounded on decode: `decode` refits both against the declared output length
// (mirroring the encoder's `fit_order`/`fit_bits`), so an oversized `k`/`hb` in a
// tiny stream is rejected as malformed rather than allocated. That closes the
// former parameter-driven OOM (#142), so the target no longer skips large
// `k`/`hb` — libFuzzer is free to explore them, with a tight `-rss_limit_mb` to
// catch any regression of the bound.
fuzz_target!(|data: &[u8]| {
    let _ = fqxv_seq::decode(data);
});
