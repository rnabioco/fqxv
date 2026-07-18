//! Corruption-robustness harness for the container decode surface.
//!
//! Every public decode entry point must return `Ok` or a clean `Err` on any
//! input — a corrupt or truncated archive must never panic or abort the process
//! (the release profile is `panic = "abort"`, so a panic is a crash). These tests
//! feed mutated and arbitrary bytes to `decompress`, `decompress_recover`,
//! `verify`, `inspect`, and `expected_reads` and assert only that none of them
//! panic; the returned `Result` is intentionally ignored.
//!
//! Note the container verifies a block/frame CRC before handing a payload to a
//! sub-codec, so random mutation here rarely reaches a codec's internals — that
//! surface is fuzzed directly in each codec crate's own tests. This harness
//! guards the container's own framing/index/header parsing.

use std::io::{sink, Cursor};
use std::sync::LazyLock;

use fqxv::{compress, decompress, decompress_recover, expected_reads, inspect, verify, Params};
use proptest::prelude::*;

/// Synthetic interleaved-ish FASTQ with a little duplication so the reorder
/// layout has something to cluster.
fn synth_fastq(n: usize) -> Vec<u8> {
    let seqs = [
        b"ACGTACGTACGTACGT",
        b"TTGGCCAATTGGCCAA",
        b"ACGTACGTACGTACGT",
    ];
    let mut out = Vec::new();
    for i in 0..n {
        let s = seqs[i % seqs.len()];
        out.extend_from_slice(format!("@read.{i} lane:1:{i}\n").as_bytes());
        out.extend_from_slice(s);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&vec![b'I'; s.len()]);
        out.push(b'\n');
    }
    out
}

/// Compress `fastq` with `params` into an in-memory archive.
fn archive(fastq: &[u8], params: Params) -> Vec<u8> {
    let mut out = Vec::new();
    compress(fastq, &mut out, params).expect("compress a valid archive");
    out
}

/// The cheap, pool-free structural readers — safe to call in a tight exhaustive
/// loop. Each must not panic on `bytes`.
fn probe_structural(bytes: &[u8]) {
    let _ = verify(Cursor::new(bytes), 1);
    let _ = inspect(Cursor::new(bytes));
    let _ = expected_reads(Cursor::new(bytes));
}

/// The full decode entries (these spin a rayon pool), for the bounded proptest.
fn probe_decode(bytes: &[u8]) {
    let _ = decompress(Cursor::new(bytes), sink(), 1);
    let _ = decompress_recover(Cursor::new(bytes), sink(), 1);
}

fn plain_params() -> Params {
    Params {
        threads: 1,
        ..Params::default()
    }
}

fn reorder_params() -> Params {
    Params {
        reorder: true,
        threads: 1,
        ..Params::default()
    }
}

#[test]
fn exhaustive_single_byte_and_truncation_never_panics() {
    // Plain layout only: its structural readers are pool-free and, crucially,
    // `verify` is a cheap whole-file CRC pass (the reorder layout's `verify` runs
    // a full decode, far too slow to call once per byte). Reorder-layout
    // robustness is covered by `mutated_archive_never_panics` below. Flip every
    // byte, then truncate at every length, and confirm nothing panics.
    let good = archive(&synth_fastq(6), plain_params());
    for i in 0..good.len() {
        let mut m = good.clone();
        m[i] ^= 0xFF;
        probe_structural(&m);
    }
    for len in 0..=good.len() {
        probe_structural(&good[..len]);
    }
}

// Base archives compressed once (the reorder path is slow), then cloned and
// mutated per proptest case so each case is just mutate + decode.
static PLAIN_ARCHIVE: LazyLock<Vec<u8>> =
    LazyLock::new(|| archive(&synth_fastq(10), plain_params()));
static REORDER_ARCHIVE: LazyLock<Vec<u8>> =
    LazyLock::new(|| archive(&synth_fastq(10), reorder_params()));

proptest! {
    // Each case spins a decode + recovery pool; a modest count is plenty for a
    // per-push regression net (proptest explores more across CI runs), and the
    // exhaustive test above already sweeps every byte of a plain archive.
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Mutate a valid archive (byte substitutions, including runs that can push a
    /// length/count field to an absurd value, plus a truncation) and confirm no
    /// decode entry panics.
    #[test]
    fn mutated_archive_never_panics(
        reorder in any::<bool>(),
        subs in prop::collection::vec((any::<usize>(), any::<u8>()), 0..12),
        wipe in proptest::option::of((any::<usize>(), 1usize..9)),
        truncate in proptest::option::of(any::<usize>()),
    ) {
        let mut bytes = if reorder {
            REORDER_ARCHIVE.clone()
        } else {
            PLAIN_ARCHIVE.clone()
        };
        for (pos, val) in subs {
            if !bytes.is_empty() {
                let i = pos % bytes.len();
                bytes[i] = val;
            }
        }
        // Set a small window to 0xFF — a cheap way to drive a length field to a
        // huge value, the exact shape that used to trigger an abort.
        if let (Some((pos, w)), false) = (wipe, bytes.is_empty()) {
            let i = pos % bytes.len();
            for b in bytes.iter_mut().skip(i).take(w) {
                *b = 0xFF;
            }
        }
        if let Some(t) = truncate {
            if !bytes.is_empty() {
                bytes.truncate(t % bytes.len());
            }
        }
        probe_structural(&bytes);
        probe_decode(&bytes);
    }

    /// Arbitrary bytes straight into every entry point (shallow but free).
    #[test]
    fn arbitrary_bytes_never_panic(bytes in prop::collection::vec(any::<u8>(), 0..512)) {
        probe_structural(&bytes);
        probe_decode(&bytes);
    }
}
