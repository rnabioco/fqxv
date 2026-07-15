//! Assembler A/B harness: compares fqxv's rightward-only assembler (+ overlap
//! merge) against the SPRING-style prototypes — [`assemble_global_bidir`]
//! (leftward extension) and [`assemble_global_refdriven`] (reference-driven
//! greedy) — by reference contig count and base count. This is how the assembler
//! lever (the real gap to SPRING's ~55 MB / 12.5 MB reference) was measured; see
//! the design issue for the findings.
//!
//! Usage: cargo run --release -p fqxv --example containment_scan -- <fastq>

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

const REORDER_K: usize = 15;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: containment_scan <fastq>");

    let mut reader =
        noodles_fastq::io::Reader::new(BufReader::new(File::open(&path).expect("open fastq")));
    let mut rec = noodles_fastq::Record::default();
    let mut lens: Vec<u32> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    while reader.read_record(&mut rec).expect("read") != 0 {
        let s: &[u8] = rec.sequence();
        lens.push(s.len() as u32);
        seq.extend_from_slice(s);
    }
    eprintln!("{path}: {} reads", lens.len());

    let plan = fqxv_reorder::plan(&lens, &seq, REORDER_K);
    let n = lens.len();
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in &lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);
    let cl_reads: Vec<Vec<u8>> = plan
        .order
        .iter()
        .map(|&oi| {
            let oi = oi as usize;
            let s = &seq[offs[oi]..offs[oi + 1]];
            if plan.flip[oi] {
                fqxv_reorder::revcomp(s)
            } else {
                s.to_vec()
            }
        })
        .collect();
    let cl_anchors: Vec<u32> = plan
        .order
        .iter()
        .map(|&oi| plan.anchor[oi as usize])
        .collect();
    let refs_all: Vec<&[u8]> = cl_reads.iter().map(Vec::as_slice).collect();

    let mb = |b: usize| b as f64 / 1e6;

    // Baseline: rightward-only assembler (+ overlap-merge).
    let t = Instant::now();
    let (reference, places) = fqxv_reorder::assemble_global(&refs_all, &cl_anchors);
    let base_asm = (reference.n_contigs(), reference.total_bases());
    let (merged, _pl) = fqxv_reorder::merge_reference(&refs_all, &reference, &places);
    let base_merged = (merged.n_contigs(), merged.total_bases());
    eprintln!("baseline assemble+merge: {:.1}s", t.elapsed().as_secs_f64());

    // Prototype: SPRING-style bidirectional assembler (leftward extension).
    let t = Instant::now();
    let bref = fqxv_reorder::assemble_global_bidir(&refs_all, &cl_anchors);
    let bidir_asm = (bref.n_contigs(), bref.total_bases());
    eprintln!("bidir assemble: {:.1}s", t.elapsed().as_secs_f64());

    // Prototype: SPRING-style reference-driven greedy assembler (rightward-only).
    let t = Instant::now();
    let rdref = fqxv_reorder::assemble_global_refdriven(&refs_all, &cl_anchors);
    let rd_asm = (rdref.n_contigs(), rdref.total_bases());
    eprintln!("refdriven assemble: {:.1}s", t.elapsed().as_secs_f64());

    println!("\n== assembler A/B (rightward-only vs bidirectional) ==");
    println!(
        "  rightward assemble : {:>8} contigs  {:>7.1} MB bases",
        base_asm.0,
        mb(base_asm.1)
    );
    println!(
        "  rightward +merge   : {:>8} contigs  {:>7.1} MB bases",
        base_merged.0,
        mb(base_merged.1)
    );
    println!(
        "  BIDIR assemble     : {:>8} contigs  {:>7.1} MB bases   ({:+.1}% bases vs rightward assemble)",
        bidir_asm.0,
        mb(bidir_asm.1),
        100.0 * (bidir_asm.1 as f64 - base_asm.1 as f64) / base_asm.1 as f64,
    );
    println!(
        "  REFDRIVEN assemble : {:>8} contigs  {:>7.1} MB bases   ({:+.1}% bases vs rightward assemble)",
        rd_asm.0,
        mb(rd_asm.1),
        100.0 * (rd_asm.1 as f64 - base_asm.1 as f64) / base_asm.1 as f64,
    );
    println!("  (SPRING reference ~55 MB bases / 12.5 MB coded; merge further reduces both.)");
}
