#![no_main]
//! Fuzz the container decode entries. Every one must return `Ok`/`Err` on any
//! input, never panic or abort.
//!
//! `verify` is intentionally omitted: its reorder path is the same decode
//! `decompress` already drives, its plain path is a whole-file CRC (not a
//! decode), and it would build an all-cores rayon pool per case. The entries
//! here all stay on a single-thread pool. Seed the corpus with real `.fqxv`
//! files (e.g. the ones in the repo root) so mutation gets past the magic + CRC.

use libfuzzer_sys::fuzz_target;
use std::io::{sink, Cursor};

fuzz_target!(|data: &[u8]| {
    let _ = fqxv::decompress(Cursor::new(data), sink(), 1);
    let _ = fqxv::decompress_recover(Cursor::new(data), sink(), 1);
    let _ = fqxv::inspect(Cursor::new(data));
    let _ = fqxv::expected_reads(Cursor::new(data));
});
