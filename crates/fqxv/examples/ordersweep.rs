//! Sweep the sequence-model order over one FASTQ file, holding the block size
//! fixed (a single block) so the only variable is `seq_order`. Prints the
//! compressed size, ratio, and throughput per order — the data behind the
//! default-order choice. Run: `cargo run --profile profiling --example ordersweep -- reads.fastq`.

use std::fs::File;
use std::io;
use std::time::Instant;

use fqxv::Params;

fn main() {
    let path = std::env::args().nth(1).expect("usage: ordersweep <fastq>");
    let in_size = std::fs::metadata(&path).expect("stat input").len();
    println!("order    out_bytes   ratio    secs   in-MB/s");
    for order in [7u8, 8, 9, 10, 11] {
        let params = Params {
            seq_order: order,
            block_reads: 8_000_000, // one block for a 2M-read input
            threads: 8,
            ..Default::default()
        };
        let f = File::open(&path).expect("open input");
        let t = Instant::now();
        let stats = fqxv::compress(f, io::sink(), params).expect("compress");
        let secs = t.elapsed().as_secs_f64();
        println!(
            "{:>5}  {:>11}  {:>6.3}  {:>6.1}  {:>7.1}",
            order,
            stats.out_bytes,
            in_size as f64 / stats.out_bytes as f64,
            secs,
            in_size as f64 / 1e6 / secs,
        );
    }
}
