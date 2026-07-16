#![no_main]
//! Fuzz the quality-score (fqzcomp) decoder.

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = fqxv_fqzcomp::decode(data);
});
