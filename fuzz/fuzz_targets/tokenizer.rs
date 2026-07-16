#![no_main]
//! Fuzz the read-name tokenizer decoder — the op-RLE, delta-plane, and
//! name-template paths all size allocations from stream-supplied counts.

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = fqxv_tokenizer::decode(data);
});
