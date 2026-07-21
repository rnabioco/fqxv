//! Throwaway harness for the reflzma window experiment (#197 follow-up): parse a
//! FASTQ, extract real per-read lengths + concatenated bases, run only the raw
//! LZMA sequence encode, and report bits/base and encode time. Isolates the LZMA
//! codec from the container so match-finder parameters can be A/B'd fast.
//!
//! `cargo run --release --example lzma_probe -- <fastq> [max_reads]`
//! `FQXV_LZMA_MAX_CHAIN=<n>` overrides the match-finder chain depth for sweeps.

use std::fs::File;
use std::io::BufReader;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: lzma_probe <fastq> [max_reads]");
    let limit: usize = args.next().map_or(usize::MAX, |s| s.parse().unwrap());

    let mut reader = noodles_fastq::io::Reader::new(BufReader::new(File::open(&path).unwrap()));
    let mut rec = noodles_fastq::Record::default();
    let mut lens: Vec<u32> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    while lens.len() < limit {
        if reader.read_record(&mut rec).unwrap() == 0 {
            break;
        }
        let s = rec.sequence();
        lens.push(s.len() as u32);
        seq.extend_from_slice(s);
    }
    let bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();

    let t0 = std::time::Instant::now();
    let coded = fqxv_seq::lzma::encode(&lens, &seq).unwrap();
    let dt = t0.elapsed();
    eprintln!(
        "reads={} bases={} coded={} b/base={:.3} encode={:.1}s",
        lens.len(),
        bases,
        coded.len(),
        coded.len() as f64 * 8.0 / bases as f64,
        dt.as_secs_f64(),
    );

    let (dl, ds) = fqxv_seq::lzma::decode(&coded).unwrap();
    assert_eq!(dl, lens, "length round-trip");
    assert_eq!(ds, seq, "sequence round-trip");
    eprintln!("round-trip OK");
}
