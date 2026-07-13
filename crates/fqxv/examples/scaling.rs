//! Block-level rayon scaling probe: compress the same FASTQ at increasing
//! thread counts and report speedup vs one thread.
//!
//! Blocks are fqxv's unit of parallelism, so `block_reads` must be small enough
//! that block_count >> max_threads or the curve is starved, not stalled. The
//! input is read into memory once; each trial compresses it from a cursor into
//! a sink, so only the codec work is timed (no disk/gzip in the loop).
//!
//! Usage:
//!   cargo run --release -p fqxv --example scaling -- <fastq> [block_reads] [maxthreads]

use std::io::Cursor;
use std::time::Instant;

use fqxv::{compress, Params};

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: scaling <fastq> [block_reads] [maxthreads]");
    let block_reads: usize = args
        .next()
        .map_or(32_000, |s| s.parse().expect("block_reads"));
    let max_threads: usize = args.next().map_or_else(
        || std::thread::available_parallelism().map_or(64, |n| n.get()),
        |s| s.parse().expect("maxthreads"),
    );

    let bytes = std::fs::read(&path).expect("read fastq");
    eprintln!(
        "{path}: {:.2} GiB in memory, block_reads={block_reads}, up to {max_threads} threads",
        bytes.len() as f64 / (1u64 << 30) as f64
    );

    let mk_params = |threads: usize| Params {
        block_reads,
        threads,
        ..Params::default()
    };

    // Warm up + report block count / ratio once.
    let warm = compress(Cursor::new(&bytes), std::io::sink(), mk_params(0)).expect("warm compress");
    eprintln!(
        "{} reads, {} blocks, ratio {:.2}x\n",
        warm.reads,
        warm.blocks,
        bytes.len() as f64 / warm.out_bytes as f64
    );

    let mut threads = vec![1usize];
    let mut t = 2;
    while t <= max_threads {
        threads.push(t);
        t *= 2;
    }
    if *threads.last().unwrap() != max_threads {
        threads.push(max_threads);
    }

    let mb = bytes.len() as f64 / 1e6;
    let mut base = 0.0f64;
    println!(
        "{:>7} {:>10} {:>10} {:>9} {:>7} {:>6}",
        "threads", "wall_s", "in MB/s", "speedup", "ideal", "eff%"
    );
    for &n in &threads {
        // Best of 2 runs to damp scheduler noise on a shared node.
        let mut best = f64::INFINITY;
        for _ in 0..2 {
            let t = Instant::now();
            let s = compress(Cursor::new(&bytes), std::io::sink(), mk_params(n)).expect("compress");
            std::hint::black_box(s);
            best = best.min(t.elapsed().as_secs_f64());
        }
        if n == 1 {
            base = best;
        }
        println!(
            "{:>7} {:>10.3} {:>10.0} {:>8.2}x {:>6}x {:>5.0}%",
            n,
            best,
            mb / best,
            base / best,
            n,
            100.0 * (base / best) / n as f64,
        );
    }
}
