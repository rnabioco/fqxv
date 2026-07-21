//! Sweep the multi-reference tiler's two ratio levers — the alignment band and
//! the best-of-N reference fan-out (`tile_band` × `tile_max_refs`) — over one
//! long-read block, reporting bits/base and encode wall-time for each point.
//!
//! This is the fast inner loop for the CoLoRd-parity follow-up: it calls
//! [`fqxv_lroverlap::tile_encode`] directly on a real ONT block held in memory, so
//! a whole grid runs in one process without the container overhead. The container
//! only adopts the tiling candidate when it wins `keep_smaller`, so the number to
//! watch is bits/base of the coded block versus encode time.
//!
//! ```text
//! cargo run --release -p fqxv-lroverlap --example tile_sweep -- reads.fastq [max_mb]
//! # grid overrides (space/comma-separated):
//! FQXV_TILE_BANDS="256 512 768 1024" FQXV_TILE_REFS="1 2 3 4 6" \
//!   cargo run --release -p fqxv-lroverlap --example tile_sweep -- reads.fastq 220
//! ```
//!
//! `max_mb` caps the block at roughly that many megabases (whole reads only),
//! standing in for one container block; omit it to tile the whole file.

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use fqxv_lroverlap::{EncodeOpts, Sketch, tile_encode};

/// Read a FASTQ's sequence lines into `(lens, seq)`, stopping once the running
/// base count would exceed `max_bases` (0 = no cap). Whole reads only, so the
/// truncated block still starts on a read boundary.
fn read_fastq(path: &str, max_bases: u64) -> (Vec<u32>, Vec<u8>) {
    let f = BufReader::with_capacity(1 << 20, File::open(path).expect("open fastq"));
    let (mut lens, mut seq) = (Vec::new(), Vec::new());
    let mut total: u64 = 0;
    for (i, line) in f.lines().enumerate() {
        if i % 4 != 1 {
            continue;
        }
        let l = line.expect("line");
        if max_bases != 0 && total + l.len() as u64 > max_bases && !lens.is_empty() {
            break;
        }
        total += l.len() as u64;
        lens.push(l.len() as u32);
        seq.extend_from_slice(l.as_bytes());
    }
    (lens, seq)
}

/// Parse a whitespace/comma-separated env list of `usize`, or a default.
fn grid(var: &str, default: &[usize]) -> Vec<usize> {
    match env::var(var) {
        Ok(v) => v
            .split(|c: char| c.is_whitespace() || c == ',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().expect("grid value"))
            .collect(),
        Err(_) => default.to_vec(),
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: tile_sweep <reads.fastq> [max_mb]");
        std::process::exit(2);
    }
    let max_bases = args
        .get(2)
        .and_then(|s| s.parse::<u64>().ok())
        .map_or(0, |mb| mb * 1_000_000);

    let (lens, seq) = read_fastq(&args[1], max_bases);
    let total_bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    let sketch = Sketch::ont();
    println!(
        "block: {} reads · {:.1} Mbase · sketch w={} k={}",
        lens.len(),
        total_bases as f64 / 1e6,
        sketch.w,
        sketch.k,
    );

    let bands = grid("FQXV_TILE_BANDS", &[256, 512, 768]);
    let refs = grid("FQXV_TILE_REFS", &[1, 2, 3, 4]);

    println!(
        "{:>6} {:>5} {:>14} {:>10} {:>9} {:>8}",
        "band", "refs", "coded_B", "b/base", "vs 1ref", "enc_s",
    );
    // Baseline (band 256, 1 ref) if it is in the grid, else the first point coded.
    let mut baseline: Option<f64> = None;
    for &band in &bands {
        for &nref in &refs {
            let opts = EncodeOpts {
                sketch,
                tile_band: band,
                tile_max_refs: nref,
                ..EncodeOpts::default()
            };
            let t = Instant::now();
            let coded = tile_encode(&lens, &seq, &opts).expect("tile_encode");
            let secs = t.elapsed().as_secs_f64();
            let bpb = (coded.len() as f64 * 8.0) / total_bases as f64;
            let rel = match baseline {
                Some(b) => format!("{:+.2}%", (bpb / b - 1.0) * 100.0),
                None => {
                    baseline = Some(bpb);
                    "  base".to_string()
                }
            };
            println!(
                "{band:>6} {nref:>5} {:>14} {bpb:>10.4} {rel:>9} {secs:>8.1}",
                coded.len(),
            );
        }
    }
}
