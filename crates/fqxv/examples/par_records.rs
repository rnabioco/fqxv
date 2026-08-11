//! Record-consumer scaling probe: the serial [`fqxv::RecordReader`] funnel vs
//! the parallel [`fqxv::decompress_records_par`] visitor on the same archive.
//!
//! Runs a cheap record consumer — count records, sum bases, fold a 64-bit
//! per-record checksum — the shape of a k-mer counter or filter's hot loop, and
//! reports wall time per repetition. The checksum makes the two paths
//! cross-check each other: it is order-independent (a wrapping sum of per-record
//! xxh3 digests), so serial and parallel runs must agree exactly.
//!
//! Usage: cargo run --release -p fqxv --example par_records -- <archive.fqxv> <serial|par> [threads] [reps]

use std::time::Instant;
use xxhash_rust::xxh3::Xxh3;

/// The per-record fold both modes run: records, bases, order-free checksum.
#[derive(Default, Clone, Copy)]
struct Acc {
    records: u64,
    bases: u64,
    checksum: u64,
}

impl Acc {
    fn push(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) {
        let mut h = Xxh3::new();
        h.update(name);
        h.update(seq);
        h.update(qual);
        self.records += 1;
        self.bases += seq.len() as u64;
        self.checksum = self.checksum.wrapping_add(h.digest());
    }
    fn merge(self, other: Acc) -> Acc {
        Acc {
            records: self.records + other.records,
            bases: self.bases + other.bases,
            checksum: self.checksum.wrapping_add(other.checksum),
        }
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let usage = "usage: par_records <archive.fqxv> <serial|par> [threads] [reps]";
    let path = args.next().expect(usage);
    let mode = args.next().expect(usage);
    let threads: usize = args.next().map_or(0, |s| s.parse().expect("threads"));
    let reps: usize = args.next().map_or(3, |s| s.parse().expect("reps"));

    for rep in 0..reps {
        let file = std::fs::File::open(&path).expect("open archive");
        let t0 = Instant::now();
        let acc = match mode.as_str() {
            "serial" => {
                // The pull iterator: block-parallel decode behind a bounded
                // channel, records consumed by this one thread.
                let mut acc = Acc::default();
                for rec in fqxv::RecordReader::new(file, threads) {
                    let rec = rec.expect("record");
                    acc.push(&rec.name, &rec.seq, &rec.qual);
                }
                acc
            }
            "par" => {
                let (acc, _stats) = fqxv::decompress_records_par(
                    file,
                    threads,
                    Acc::default,
                    |acc, _index, rec| {
                        acc.push(rec.name, rec.seq, rec.qual);
                        Ok(())
                    },
                    Acc::merge,
                )
                .expect("decompress_records_par");
                acc
            }
            other => panic!("unknown mode {other:?}; {usage}"),
        };
        let dt = t0.elapsed();
        println!(
            "{mode}\tthreads={threads}\trep={rep}\t{:.3}s\trecords={}\tbases={}\tchecksum={:#018x}",
            dt.as_secs_f64(),
            acc.records,
            acc.bases,
            acc.checksum
        );
    }
}
