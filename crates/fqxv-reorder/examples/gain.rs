//! Measure the sequence-compression gain from reordering.
//!
//! Usage: `cargo run --release -p fqxv-reorder --example gain -- <seq-file> <read-len> [k] [order]`
//! Compresses the sequence stream with fqxv-seq before and after minimizer
//! clustering, and reports bits/base plus the permutation cost.

use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: gain <file> <read-len> [k] [order]");
    let read_len: usize = args.next().expect("read-len").parse().unwrap();
    let k: usize = args.next().unwrap_or_else(|| "15".into()).parse().unwrap();
    let order: usize = args.next().unwrap_or_else(|| "11".into()).parse().unwrap();

    let seq = std::fs::read(&path).expect("read");
    let n = seq.len() / read_len;
    let seq = &seq[..n * read_len];
    let lens = vec![read_len as u32; n];

    let baseline = fqxv_seq::encode(&lens, seq, order).unwrap().len();

    let t0 = Instant::now();
    let p = fqxv_reorder::plan(&lens, seq, k);
    let plan_s = t0.elapsed().as_secs_f64();

    // Build the reordered, canonically-oriented sequence stream.
    let mut reordered = Vec::with_capacity(seq.len());
    for &oi in &p.order {
        let read = &seq[oi as usize * read_len..(oi as usize + 1) * read_len];
        if p.flip[oi as usize] {
            reordered.extend_from_slice(&fqxv_reorder::revcomp(read));
        } else {
            reordered.extend_from_slice(read);
        }
    }
    let reordered_c = fqxv_seq::encode(&lens, &reordered, order).unwrap().len();

    // Order-preserving cost: store the permutation (~log2(n!) bits; here as raw
    // 4-byte indices compressed by the same order-0 path, a loose upper bound).
    let perm_bytes: Vec<u8> = p.order.iter().flat_map(|x| x.to_le_bytes()).collect();
    let perm_c = fqxv_seq::encode(&[perm_bytes.len() as u32], &perm_bytes, 4)
        .map(|v| v.len())
        .unwrap_or(perm_bytes.len());

    let bpb = |bytes: usize| (bytes * 8) as f64 / (n * read_len) as f64;
    println!("reads {n}  read_len {read_len}  k {k}  order {order}");
    println!(
        "  baseline (no reorder)  {baseline:>10} bytes  {:.3} bits/base",
        bpb(baseline)
    );
    println!(
        "  reordered              {reordered_c:>10} bytes  {:.3} bits/base",
        bpb(reordered_c)
    );
    println!(
        "  reordered + perm       {:>10} bytes  {:.3} bits/base  (order-preserving)",
        reordered_c + perm_c,
        bpb(reordered_c + perm_c)
    );
    println!(
        "  gain: {:.1}% (reorder-free) / {:.1}% (order-preserving);  plan {plan_s:.1}s",
        100.0 * (baseline as f64 - reordered_c as f64) / baseline as f64,
        100.0 * (baseline as f64 - (reordered_c + perm_c) as f64) / baseline as f64,
    );
}
