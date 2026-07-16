//! Find read-to-read overlaps in a FASTQ and report what was found.
//!
//! Validates the detector against real data (the unit tests are synthetic) and
//! is the harness for comparing recall against minimap2's all-vs-all.
//!
//! ```text
//! cargo run --release --example overlaps -- reads.fastq [ont|hifi]
//! ```

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use rayon::prelude::*;

use fqxv_lroverlap::{find_overlaps, ChainOpts, Index, Repeat, Sketch};

/// Minimal FASTQ reader: returns concatenated sequence plus per-read lengths.
/// Non-ACGT is passed through — `minimizers` breaks k-mer runs on it.
fn read_fastq(path: &str) -> (Vec<u32>, Vec<u8>) {
    let f = BufReader::with_capacity(1 << 20, File::open(path).expect("open fastq"));
    let (mut lens, mut seq) = (Vec::new(), Vec::new());
    for (i, line) in f.lines().enumerate() {
        // Record is 4 lines; the sequence is line 1 of each.
        if i % 4 != 1 {
            continue;
        }
        let l = line.expect("read line");
        lens.push(l.len() as u32);
        seq.extend_from_slice(l.as_bytes());
    }
    (lens, seq)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: overlaps <reads.fastq> [ont|hifi]");
        std::process::exit(2);
    }
    let sketch = match args.get(2).map(String::as_str) {
        Some("hifi") => Sketch::hifi(),
        _ => Sketch::ont(),
    };

    let t0 = Instant::now();
    let (lens, seq) = read_fastq(&args[1]);
    let total: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    println!(
        "reads {} · bases {} · mean len {} · parsed in {:.1}s",
        lens.len(),
        total,
        total / lens.len().max(1) as u64,
        t0.elapsed().as_secs_f64()
    );
    println!("sketch: w={} k={}", sketch.w, sketch.k);

    let t1 = Instant::now();
    let idx = Index::build(&lens, &seq, sketch, Repeat::default()).expect("build index");
    println!(
        "index: {} occurrences · {} distinct minimizers dropped as repetitive · \
         busiest kept {}x · {:.1}s",
        idx.len(),
        idx.dropped(),
        idx.max_kept(),
        t1.elapsed().as_secs_f64()
    );

    // Per-read byte offsets.
    let mut offs = Vec::with_capacity(lens.len() + 1);
    let mut acc = 0usize;
    for &l in &lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    let t2 = Instant::now();
    // (chain records, distinct partner reads). These differ: a pair whose
    // alignment fragments into several chains yields several records but is
    // still one usable coding partner, and the codec only ever wants the best
    // chain per partner — so `partners` is the number that actually matters.
    let counts: Vec<(usize, usize)> = (0..lens.len())
        .into_par_iter()
        .map(|r| {
            let s = &seq[offs[r]..offs[r + 1]];
            let ov = find_overlaps(&idx, r as u32, s, ChainOpts::default());
            let mut t: Vec<u32> = ov.iter().map(|o| o.target).collect();
            t.sort_unstable();
            t.dedup();
            (ov.len(), t.len())
        })
        .collect();
    let elapsed = t2.elapsed().as_secs_f64();

    let total_ov: usize = counts.iter().map(|c| c.0).sum();
    let total_partners: usize = counts.iter().map(|c| c.1).sum();
    let orphans = counts.iter().filter(|c| c.1 == 0).count();
    let mut sorted: Vec<usize> = counts.iter().map(|c| c.1).collect();
    sorted.sort_unstable();
    let median = sorted.get(sorted.len() / 2).copied().unwrap_or(0);
    let n = lens.len().max(1) as f64;

    println!("--- overlaps ---");
    println!(
        "  chain records    : {total_ov} (mean {:.1}/read)",
        total_ov as f64 / n
    );
    println!(
        "  distinct partners: {total_partners} (mean {:.1}/read)  <- what the codec uses",
        total_partners as f64 / n
    );
    println!(
        "  records/partner  : {:.2}  (>1 = alignments fragmenting into several chains)",
        total_ov as f64 / total_partners.max(1) as f64
    );
    println!("  median partners  : {median}");
    println!("  max partners     : {}", sorted.last().copied().unwrap_or(0));
    println!(
        "  reads with none  : {orphans} ({:.2}%)  <- these fall back to order-k coding",
        100.0 * orphans as f64 / n
    );
    println!(
        "  search time      : {elapsed:.1}s ({:.1} Mbase/s)",
        (total as f64 / 1e6) / elapsed
    );
}
