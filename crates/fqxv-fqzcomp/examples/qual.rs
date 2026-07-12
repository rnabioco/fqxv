//! fqzcomp quality model round-trip + ratio on a raw quality stream.
//!
//! Usage: `cargo run --release -p fqxv-fqzcomp --example qual -- <qual-file> <read-len>`
//! The file is the concatenated quality bytes (no newlines); read-len is the
//! fixed per-read length.

use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: qual <file> <read-len>");
    let read_len: u32 = args
        .next()
        .expect("read-len")
        .parse()
        .expect("read-len int");

    let quals = std::fs::read(&path).expect("read input");
    let n_reads = quals.len() / read_len as usize;
    let quals = &quals[..n_reads * read_len as usize];
    let lens = vec![read_len; n_reads];

    let t0 = Instant::now();
    let enc =
        fqxv_fqzcomp::encode(&lens, quals, fqxv_fqzcomp::QualityBinning::Lossless).expect("encode");
    let c_secs = t0.elapsed().as_secs_f64();

    let t1 = Instant::now();
    let (_lens, dec) = fqxv_fqzcomp::decode(&enc).expect("decode");
    let d_secs = t1.elapsed().as_secs_f64();
    assert_eq!(dec, quals, "round-trip mismatch!");

    let bpb = (enc.len() * 8) as f64 / quals.len() as f64;
    let ratio = quals.len() as f64 / enc.len() as f64;
    let mbps = |s: f64| (quals.len() as f64 / 1e6) / s;
    println!("file        {path}  ({} reads x {read_len}bp)", n_reads);
    println!("input       {} bytes", quals.len());
    println!("fqzcomp     {} bytes", enc.len());
    println!("ratio       {ratio:.3}x  ({bpb:.3} bits/byte)");
    println!("encode      {c_secs:.2}s  ({:.0} MB/s)", mbps(c_secs));
    println!("decode      {d_secs:.2}s  ({:.0} MB/s)", mbps(d_secs));
    println!("round-trip  OK");
}
