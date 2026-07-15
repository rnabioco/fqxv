//! Measure the fqxv quality codec directly on a FASTQ's quality stream.
//!
//! Usage: `qual_bench <reads.fastq>`
//!
//! Isolates the quality codec from all container/CLI overhead: reads the
//! quality lines, encodes them losslessly, and reports the compressed size
//! and bits/symbol, then verifies an exact round-trip.

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::ExitCode;

use fqxv_fqzcomp::{decode, encode, QualityBinning};

fn main() -> ExitCode {
    let path = match env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("usage: qual_bench <reads.fastq>");
            return ExitCode::FAILURE;
        }
    };

    let file = match File::open(&path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error opening {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let reader = BufReader::new(file);

    let mut quals: Vec<u8> = Vec::new();
    let mut lens: Vec<u32> = Vec::new();

    for (i, line) in reader.split(b'\n').enumerate() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("error reading {path}: {e}");
                return ExitCode::FAILURE;
            }
        };
        // Quality is line 3 (0-indexed) of each 4-line record.
        if i % 4 == 3 {
            lens.push(line.len() as u32);
            quals.extend_from_slice(&line);
        }
    }

    let n_reads = lens.len();
    let total = quals.len();

    let encoded = match encode(&lens, &quals, QualityBinning::Lossless) {
        Ok(e) => e,
        Err(e) => {
            eprintln!("encode failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let compressed = encoded.len();
    let bits_per_symbol = if total > 0 {
        (compressed as f64 * 8.0) / total as f64
    } else {
        0.0
    };

    println!("reads:            {n_reads}");
    println!("total qual bytes: {total}");
    println!("compressed bytes: {compressed}");
    println!("bits/symbol:      {bits_per_symbol:.4}");

    let (dec_lens, dec_quals) = match decode(&encoded) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("decode failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    assert_eq!(dec_lens, lens, "lengths did not round-trip");
    assert_eq!(dec_quals, quals, "qualities did not round-trip");
    println!("roundtrip OK");

    ExitCode::SUCCESS
}
