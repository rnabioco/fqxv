//! Op-mix of the clustered contig-assembly sequence codec, on real data.
//!
//! This answers a single question about the reorder path: when reads are
//! globally minimizer-clustered and coded against a growing consensus contig,
//! what fraction land in each op —
//!
//!   * `MATCH`   — byte-identical to the previous read (free),
//!   * `CONTIG`  — placed on the consensus, coded as offset + a few mismatches
//!                 (cheap; only the novel tail is context-coded), or
//!   * `LITERAL` — seeds a new contig, coded from scratch by the `fqxv_seq`
//!                 order-k context model (expensive).
//!
//! It runs the same global clustering (`fqxv_reorder::plan`) and per-block
//! classification (`fqxv_reorder::op_stats`) the container uses on the
//! `--order any` path, so the numbers match what the encoder actually does.
//! A high LITERAL share — in reads or, more tellingly, in *bases* — is the
//! signal that single-minimizer clustering is fragmenting contigs and leaving
//! cross-read redundancy that SPRING's assembly would capture.
//!
//! Usage: cargo run --release -p fqxv --example op_mix -- <fastq> [k]

use std::fs::File;
use std::io::BufReader;

/// Minimizer k for clustering — matches the container's `REORDER_K`.
const REORDER_K: usize = 15;
/// Reads per block — matches the container's `REORDER_BLOCK_READS` (the contig
/// resets at each block boundary, so op-mix must block the same way).
const REORDER_BLOCK_READS: usize = 1 << 18;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: op_mix <fastq> [k]");
    let k: usize = args.next().map_or(REORDER_K, |s| s.parse().expect("k"));

    // Sequence stream only — op-mix is a property of the sequence codec.
    let mut reader =
        noodles_fastq::io::Reader::new(BufReader::new(File::open(&path).expect("open fastq")));
    let mut rec = noodles_fastq::Record::default();
    let mut lens: Vec<u32> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    while reader.read_record(&mut rec).expect("read record") != 0 {
        let s: &[u8] = rec.sequence();
        lens.push(s.len() as u32);
        seq.extend_from_slice(s);
    }
    let n = lens.len();
    eprintln!("{path}: {n} reads, {} seq bases, clustering k={k}", seq.len());

    // Global minimizer clustering — the whole-file plan the container computes.
    let plan = fqxv_reorder::plan(&lens, &seq, k);

    // Cumulative offsets into the concatenated seq.
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in &lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // Clustered, oriented reads + their minimizer anchors, exactly as the
    // container builds them before calling `encode_clustered`.
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
    let cl_anchors: Vec<u32> = plan.order.iter().map(|&oi| plan.anchor[oi as usize]).collect();

    // Replay classification per block (contig resets each block, as in encode).
    let bsz = REORDER_BLOCK_READS.max(1);
    let mut agg = fqxv_reorder::OpStats::default();
    let mut s = 0usize;
    while s < n {
        let e = (s + bsz).min(n);
        let refs: Vec<&[u8]> = cl_reads[s..e].iter().map(Vec::as_slice).collect();
        agg.merge(&fqxv_reorder::op_stats(&refs, &cl_anchors[s..e]));
        s = e;
    }

    let reads = agg.reads.max(1) as f64;
    let bases = agg.total_bases.max(1) as f64;
    // Bases the context coder actually sees from scratch: whole literals plus the
    // novel tails of contig reads. This — not the read-count split — is what
    // drives the sequence-stream size.
    let scratch_bases = agg.literal_bases + agg.novel_tail_bases;

    println!("\nop-mix (share of reads):");
    println!("  {:<8} {:>12} {:>8}", "op", "reads", "%");
    let row = |name: &str, cnt: usize| {
        println!(
            "  {:<8} {:>12} {:>7.2}%",
            name,
            cnt,
            100.0 * cnt as f64 / reads
        );
    };
    row("MATCH", agg.matches);
    row("CONTIG", agg.contigs);
    row("LITERAL", agg.literals);

    println!("\nbase accounting (share of {} total bases):", agg.total_bases);
    let brow = |name: &str, b: u64| {
        println!(
            "  {:<22} {:>14} {:>7.2}%",
            name,
            b,
            100.0 * b as f64 / bases
        );
    };
    brow("MATCH (free)", agg.match_bases);
    brow("CONTIG overlap (diff)", agg.contig_overlap_bases);
    brow("CONTIG novel tail", agg.novel_tail_bases);
    brow("LITERAL (from scratch)", agg.literal_bases);
    brow("-> context-coded total", scratch_bases);

    let contigs = agg.contigs.max(1) as f64;
    println!(
        "\nCONTIG reads: {:.3} mismatches/read, {:.1} overlap bases/read, {:.1} novel tail/read",
        agg.contig_mismatches as f64 / contigs,
        agg.contig_overlap_bases as f64 / contigs,
        agg.novel_tail_bases as f64 / contigs,
    );
    println!(
        "Headline: {:.1}% of bases are context-coded from scratch \
         (LITERAL {:.1}% + novel tails {:.1}%).",
        100.0 * scratch_bases as f64 / bases,
        100.0 * agg.literal_bases as f64 / bases,
        100.0 * agg.novel_tail_bases as f64 / bases,
    );
}
