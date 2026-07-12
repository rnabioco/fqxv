//! Order-0 rANS round-trip + ratio on a raw byte file.
//!
//! Usage: `cargo run --release -p fqxv-rans --example entropy -- <file>`
//! Reports input/output sizes, ratio, bits/byte, and verifies the round-trip.
//! Handy for pointing at an extracted FASTQ quality or sequence stream.

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: entropy <file>");
    let data = std::fs::read(&path).expect("read input");
    println!("file  {path}  ({} bytes)", data.len());
    println!(
        "{:<8} {:>12} {:>8} {:>10} {:>9} {:>9}",
        "order", "bytes", "ratio", "bits/byte", "enc MB/s", "dec MB/s"
    );
    for (label, order) in [
        ("order-0", fqxv_rans::Order::Zero),
        ("order-1", fqxv_rans::Order::One),
    ] {
        let t0 = Instant::now();
        let enc = fqxv_rans::encode(&data, order).expect("encode");
        let c_secs = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let dec = fqxv_rans::decode(&enc).expect("decode");
        let d_secs = t1.elapsed().as_secs_f64();
        assert_eq!(dec, data, "round-trip mismatch ({label})!");

        let ratio = data.len() as f64 / enc.len() as f64;
        let bpb = (enc.len() * 8) as f64 / data.len() as f64;
        let mbps = |s: f64| (data.len() as f64 / 1e6) / s;
        println!(
            "{label:<8} {:>12} {ratio:>7.3}x {bpb:>10.3} {:>9.0} {:>9.0}",
            enc.len(),
            mbps(c_secs),
            mbps(d_secs)
        );
    }
}
