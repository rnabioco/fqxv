//! Order-0 rANS round-trip + ratio on a raw byte file.
//!
//! Usage: `cargo run --release -p fqxv-rans --example entropy -- <file>`
//! Reports input/output sizes, ratio, bits/byte, and verifies the round-trip.
//! Handy for pointing at an extracted FASTQ quality or sequence stream.

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: entropy <file>");
    let data = std::fs::read(&path).expect("read input");

    let t0 = Instant::now();
    let enc = fqxv_rans::encode(&data, fqxv_rans::Order::Zero).expect("encode");
    let c_secs = t0.elapsed().as_secs_f64();

    let t1 = Instant::now();
    let dec = fqxv_rans::decode(&enc).expect("decode");
    let d_secs = t1.elapsed().as_secs_f64();

    assert_eq!(dec, data, "round-trip mismatch!");

    let ratio = data.len() as f64 / enc.len() as f64;
    let bits_per_byte = (enc.len() * 8) as f64 / data.len() as f64;
    let mbps = |secs: f64| (data.len() as f64 / 1e6) / secs;
    println!("file           {path}");
    println!("input          {} bytes", data.len());
    println!("rANS order-0   {} bytes", enc.len());
    println!("ratio          {ratio:.3}x  ({bits_per_byte:.3} bits/byte)");
    println!("encode         {c_secs:.2}s  ({:.0} MB/s)", mbps(c_secs));
    println!("decode         {d_secs:.2}s  ({:.0} MB/s)", mbps(d_secs));
    println!("round-trip     OK");
}
