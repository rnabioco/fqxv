//! fqxv-seq round-trip + ratio on a raw sequence stream.
//!
//! Usage: `cargo run --release -p fqxv-seq --example seq -- <seq-file> <read-len> [order]`

use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: seq <file> <read-len> [order]");
    let read_len: u32 = args.next().expect("read-len").parse().expect("int");
    let order: usize = args
        .next()
        .unwrap_or_else(|| "8".into())
        .parse()
        .expect("int");

    let seq = std::fs::read(&path).expect("read input");
    let n = seq.len() / read_len as usize;
    let seq = &seq[..n * read_len as usize];
    let lens = vec![read_len; n];

    let t0 = Instant::now();
    let enc = fqxv_seq::encode(&lens, seq, order).expect("encode");
    let c = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let (_l, dec) = fqxv_seq::decode(&enc).expect("decode");
    let d = t1.elapsed().as_secs_f64();
    assert_eq!(dec, seq, "round-trip mismatch!");

    let bpb = (enc.len() * 8) as f64 / seq.len() as f64;
    let mbps = |s: f64| (seq.len() as f64 / 1e6) / s;
    println!(
        "order-{order:<2}  {} -> {} bytes  {:.3} bits/base  enc {:.0} MB/s  dec {:.0} MB/s  OK",
        seq.len(),
        enc.len(),
        bpb,
        mbps(c),
        mbps(d)
    );
}
