#![no_main]
//! Fuzz the order-k sequence decoder (context model + exception list).

use libfuzzer_sys::fuzz_target;

// The order-k context table (`4^k` models) and hashed table (`2^hb` slots) scale
// with the header's `k` (byte 1) and `hb` (byte 3) — both valid parameters, so a
// large one is memory pressure, not a logic bug. Skip them so the target explores
// decode logic without parameter-driven OOM.
fuzz_target!(|data: &[u8]| {
    if data.get(1).is_some_and(|&k| k > 8) || data.get(3).is_some_and(|&hb| hb > 20) {
        return;
    }
    let _ = fqxv_seq::decode(data);
});
