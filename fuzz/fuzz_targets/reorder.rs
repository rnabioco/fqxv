#![no_main]
//! Fuzz the reorder (SPRING-style) decoders: the clustered read decoders (via
//! the version-dispatching `decode_clustered_auto`) and every global-reference
//! method. Each reads per-read lengths, block counts, and reference offsets from
//! the stream, so all must reject a corrupt value rather than abort.

use fqxv_reorder::GlobalReference;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = fqxv_reorder::decode_clustered_auto(data);
    let _ = GlobalReference::decode(data);
    let _ = GlobalReference::decode_blocked(data);
    let _ = GlobalReference::decode_lzma(data);
    let _ = GlobalReference::decode_packed(data);
});
