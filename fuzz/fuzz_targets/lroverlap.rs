#![no_main]
//! Fuzz the long-read overlap decoders. Each must return `Ok`/`Err` on any
//! input — never panic, abort, or OOM.
//!
//! `tile_decode` is the multi-reference tiler (`SEQ_METHOD_TILE`): its block
//! self-describes read lengths, per-tile neighbour deltas/offsets, and framed
//! entropy streams, all of which a mutated input can make inconsistent. It must
//! fail closed (a strictly-earlier neighbour id, an in-range offset, a tile whose
//! ops produce exactly its declared length, matching framed symbol counts) rather
//! than index out of bounds or allocate on an attacker-controlled length. `decode`
//! is the consensus codec's self-contained block; both share the same framing
//! primitives, so fuzzing them together exercises that shared decode surface.

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = fqxv_lroverlap::tile_decode(data);
    // NOTE: `fqxv_lroverlap::decode` (the consensus codec) is intentionally *not*
    // fuzzed here yet. Adding this first lroverlap target surfaced that its
    // `reconstruct`/`read_header` still size allocations (`vec![0u8; total_bases]`,
    // `Vec::with_capacity(n_contigs)`, per-insertion `Vec::with_capacity(k)`) from
    // unvalidated header varints — the #142 class of decode alloc bomb, but on a
    // path that predates this codec and was never fuzzed. `tile_decode` bounds all
    // of those; hardening the consensus decoder the same way is a tracked follow-up.
});
