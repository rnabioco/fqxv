//! Throughput benchmark for the read-reordering planner.
//!
//! `plan` computes a canonical minimizer per read and sorts to cluster
//! reverse-complement-aware — the stage whose sort keys hold the oriented
//! sequence. Reads are sampled from a shared reference (so minimizers genuinely
//! cluster) with ~half reverse-complemented (so both the borrow and the
//! `revcomp`-allocate key paths are exercised).
//!
//! Run on a compute node: `cargo bench -p fqxv-reorder`.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use fqxv_reorder::{plan, revcomp, DEFAULT_K};

/// Sample `n_reads` reads of `read_len` from a shared pseudo-random reference,
/// reverse-complementing about half. Returns `(lens, concatenated seq)`.
fn make_reads(n_reads: usize, read_len: usize) -> (Vec<u32>, Vec<u8>) {
    const LUT: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let ref_len = 1usize << 20; // 1 MiB reference
    let mut state = 0x1234_5678u32;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        state
    };
    let mut reference = Vec::with_capacity(ref_len);
    for _ in 0..ref_len {
        reference.push(LUT[(next() & 3) as usize]);
    }

    let mut seq = Vec::with_capacity(n_reads * read_len);
    let mut lens = Vec::with_capacity(n_reads);
    for _ in 0..n_reads {
        let pos = (next() as usize) % (ref_len - read_len);
        let window = &reference[pos..pos + read_len];
        if next() & 1 == 0 {
            seq.extend_from_slice(&revcomp(window));
        } else {
            seq.extend_from_slice(window);
        }
        lens.push(read_len as u32);
    }
    (lens, seq)
}

fn bench_reorder(c: &mut Criterion) {
    let read_len = 150usize;
    let n_reads = 50_000;
    let (lens, seq) = make_reads(n_reads, read_len);

    let mut g = c.benchmark_group("reorder_plan_50k");
    g.throughput(Throughput::Bytes(seq.len() as u64));
    g.bench_function("plan", |b| {
        b.iter(|| plan(&lens, std::hint::black_box(&seq), DEFAULT_K))
    });
    g.finish();
}

criterion_group!(benches, bench_reorder);
criterion_main!(benches);
