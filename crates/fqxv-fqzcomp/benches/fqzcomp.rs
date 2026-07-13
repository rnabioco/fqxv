//! Throughput benchmarks for the quality-score context model.
//!
//! `encode_lossless` is the default path (borrows the input, no binning copy);
//! `encode_bin8` exercises the lossy binning map; `decode` is the range-decode
//! side. Data is a skewed ~40-symbol alphabet with runs, so the model is
//! realistic (not uniform) without any RNG dependency.
//!
//! Run on a compute node: `cargo bench -p fqxv-fqzcomp`.

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use fqxv_fqzcomp::{decode, encode, QualityBinning};

/// Deterministic quality-like bytes: a skewed alphabet with short runs.
fn quality_like(n: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(n);
    let mut state = 0x1234_5678u32;
    let mut q = 30i32;
    for _ in 0..n {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        let step = (state % 7) as i32 - 3;
        q = (q + step).clamp(2, 41);
        v.push(b'!' + q as u8);
    }
    v
}

fn bench_fqzcomp(c: &mut Criterion) {
    let read_len = 150usize;
    let n_reads = (2 << 20) / read_len; // ~2 MiB of qualities
    let total = n_reads * read_len;
    let quals = quality_like(total);
    let lens = vec![read_len as u32; n_reads];
    let enc = encode(&lens, &quals, QualityBinning::Lossless).expect("encode");
    assert_eq!(decode(&enc).unwrap().1, quals, "roundtrip");

    let mut g = c.benchmark_group("fqzcomp_quality_2MiB");
    g.throughput(Throughput::Bytes(total as u64));
    g.bench_function("encode_lossless", |b| {
        b.iter(|| {
            encode(
                &lens,
                std::hint::black_box(&quals),
                QualityBinning::Lossless,
            )
            .unwrap()
        })
    });
    g.bench_function("encode_bin8", |b| {
        b.iter(|| encode(&lens, std::hint::black_box(&quals), QualityBinning::Bin8).unwrap())
    });
    g.bench_function("decode", |b| {
        b.iter(|| decode(std::hint::black_box(&enc)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, bench_fqzcomp);
criterion_main!(benches);
