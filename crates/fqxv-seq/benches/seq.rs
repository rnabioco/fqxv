//! Throughput benchmarks for the order-k ACGT sequence context model.
//!
//! Encode/decode at the maximum context order (the high-compression path). The
//! bases are deterministic pseudo-random ACGT with a sprinkling of `N`, so the
//! exception path is exercised without dominating.
//!
//! Run on a compute node: `cargo bench -p fqxv-seq`.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use fqxv_seq::{decode, encode};

/// Deterministic ACGT bytes with rare `N` (~1 in 4096), no RNG dependency.
fn acgt_like(n: usize) -> Vec<u8> {
    const LUT: [u8; 4] = *b"ACGT";
    let mut v = Vec::with_capacity(n);
    let mut state = 0x2545_f491u32;
    for _ in 0..n {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        if state & 0xfff == 0 {
            v.push(b'N');
        } else {
            v.push(LUT[(state & 3) as usize]);
        }
    }
    v
}

fn bench_seq(c: &mut Criterion) {
    let read_len = 150usize;
    let n_reads = (2 << 20) / read_len; // ~2 MiB of bases
    let total = n_reads * read_len;
    let seq = acgt_like(total);
    let lens = vec![read_len as u32; n_reads];
    // Order 8: exercises the per-base context loop without the ~42 MB model
    // allocation of MAX_ORDER dominating each iteration.
    let order = 8;
    let enc = encode(&lens, &seq, order).expect("encode");
    assert_eq!(decode(&enc).unwrap().1, seq, "roundtrip");

    let mut g = c.benchmark_group("seq_acgt_2MiB");
    g.throughput(Throughput::Bytes(total as u64));
    g.bench_function("encode_order8", |b| {
        b.iter(|| encode(&lens, std::hint::black_box(&seq), order).unwrap())
    });
    g.bench_function("decode", |b| {
        b.iter(|| decode(std::hint::black_box(&enc)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, bench_seq);
criterion_main!(benches);
