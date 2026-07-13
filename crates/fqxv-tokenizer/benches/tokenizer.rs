//! Throughput benchmarks for the read-name tokenizer.
//!
//! Realistic Illumina-style names with incrementing tile/x/y coordinates — the
//! case the per-column delta model and `natural_digits` are tuned for (each name
//! carries several numeric tokens, so the digit-counting path is hot).
//!
//! Run on a compute node: `cargo bench -p fqxv-tokenizer`.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use fqxv_tokenizer::{decode, encode};

/// Deterministic Illumina-like names: fixed instrument/run/flowcell prefix, then
/// lane/tile plus incrementing x/y coordinates (the numeric-token-heavy shape).
fn illumina_names(n: usize) -> Vec<Vec<u8>> {
    let mut v = Vec::with_capacity(n);
    let mut state = 0x9e37_79b9u32;
    let mut x = 1000u32;
    let mut y = 1000u32;
    for i in 0..n {
        // xorshift for cheap deterministic jitter in the coordinate steps.
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        x = x.wrapping_add(1 + (state & 7));
        y = y.wrapping_add(1 + ((state >> 3) & 7));
        let lane = 1 + (i as u32 % 4);
        let tile = 1101 + (i as u32 / 5000) % 100;
        let name = format!("SIM01:47:000000000-A1B2C:{lane}:{tile}:{x}:{y}");
        v.push(name.into_bytes());
    }
    v
}

fn bench_tokenizer(c: &mut Criterion) {
    let n = 50_000;
    let names_owned = illumina_names(n);
    let names: Vec<&[u8]> = names_owned.iter().map(Vec::as_slice).collect();
    let total: usize = names_owned.iter().map(Vec::len).sum();
    let enc = encode(&names).expect("encode");
    assert_eq!(decode(&enc).unwrap(), names_owned, "roundtrip");

    let mut g = c.benchmark_group("tokenizer_names");
    g.throughput(Throughput::Bytes(total as u64));
    g.bench_function("encode", |b| {
        b.iter(|| encode(std::hint::black_box(&names)).unwrap())
    });
    g.bench_function("decode", |b| {
        b.iter(|| decode(std::hint::black_box(&enc)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, bench_tokenizer);
criterion_main!(benches);
