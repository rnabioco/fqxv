//! Throughput benchmarks: scalar vs AVX2 order-0 decode, plus encode.
//!
//! Run on a compute node: `cargo bench -p fqxv-rans --features bench`.
//! (The `bench` feature exposes the internal scalar/AVX2 decoders.)

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use fqxv_rans::{Order, bench_api, encode};

/// Deterministic quality-like bytes: a skewed ~40-symbol alphabet with runs,
/// so the model is realistic (not uniform) without any RNG dependency.
fn quality_like(n: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(n);
    let mut state = 0x1234_5678u32;
    let mut q = 30i32;
    for _ in 0..n {
        // xorshift for cheap deterministic noise
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        let step = (state % 7) as i32 - 3;
        q = (q + step).clamp(2, 41);
        v.push(b'!' + q as u8);
    }
    v
}

fn bench_decode(c: &mut Criterion) {
    let n = 8 << 20; // 8 MiB
    let data = quality_like(n);
    let enc = encode(&data, Order::Zero).expect("encode");
    assert_eq!(bench_api::decode_scalar(&enc).unwrap(), data);

    let mut g = c.benchmark_group("order0_decode_8MiB");
    g.throughput(Throughput::Bytes(n as u64));
    g.bench_function(BenchmarkId::new("scalar", n), |b| {
        b.iter(|| bench_api::decode_scalar(std::hint::black_box(&enc)).unwrap())
    });
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        assert_eq!(bench_api::decode_avx2(&enc).unwrap(), data);
        g.bench_function(BenchmarkId::new("avx2", n), |b| {
            b.iter(|| bench_api::decode_avx2(std::hint::black_box(&enc)).unwrap())
        });
    }
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx512f") {
        assert_eq!(bench_api::decode_avx512(&enc).unwrap(), data);
        g.bench_function(BenchmarkId::new("avx512", n), |b| {
            b.iter(|| bench_api::decode_avx512(std::hint::black_box(&enc)).unwrap())
        });
    }
    g.finish();
}

fn bench_encode(c: &mut Criterion) {
    let n = 8 << 20;
    let data = quality_like(n);
    // Scalar and AVX2 order-0 encoders must agree byte-for-byte.
    assert_eq!(
        bench_api::encode_order0_scalar(&data),
        encode(&data, Order::Zero).unwrap()
    );
    let mut g = c.benchmark_group("encode_8MiB");
    g.throughput(Throughput::Bytes(n as u64));
    g.bench_function("order0_scalar", |b| {
        b.iter(|| bench_api::encode_order0_scalar(std::hint::black_box(&data)))
    });
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        g.bench_function("order0_avx2", |b| {
            b.iter(|| bench_api::encode_order0_avx2(std::hint::black_box(&data)))
        });
    }
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx512f") {
        assert_eq!(
            bench_api::encode_order0_avx512(&data),
            encode(&data, Order::Zero).unwrap()
        );
        g.bench_function("order0_avx512", |b| {
            b.iter(|| bench_api::encode_order0_avx512(std::hint::black_box(&data)))
        });
    }
    g.bench_function("order1", |b| {
        b.iter(|| encode(std::hint::black_box(&data), Order::One).unwrap())
    });
    g.finish();
}

criterion_group!(benches, bench_decode, bench_encode);
criterion_main!(benches);
