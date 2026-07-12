//! Tokenizer round-trip + ratio on a file of newline-separated read names.
//!
//! Usage: `cargo run --release -p fqxv-tokenizer --example names -- <names-file>`

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: names <file>");
    let data = std::fs::read(&path).expect("read input");
    let names: Vec<&[u8]> = data
        .split(|&b| b == b'\n')
        .filter(|l| !l.is_empty())
        .collect();
    let raw: usize = names.iter().map(|n| n.len() + 1).sum();

    let t0 = Instant::now();
    let enc = fqxv_tokenizer::encode(&names).expect("encode");
    let c = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let out = fqxv_tokenizer::decode(&enc).expect("decode");
    let d = t1.elapsed().as_secs_f64();
    let expect: Vec<Vec<u8>> = names.iter().map(|n| n.to_vec()).collect();
    assert_eq!(out, expect, "round-trip mismatch!");

    let mbps = |s: f64| (raw as f64 / 1e6) / s;
    println!(
        "{} names  {} -> {} bytes  {:.3}x  ({:.3} bytes/name)  enc {:.0} MB/s  dec {:.0} MB/s  OK",
        names.len(),
        raw,
        enc.len(),
        raw as f64 / enc.len() as f64,
        enc.len() as f64 / names.len() as f64,
        mbps(c),
        mbps(d),
    );
}
