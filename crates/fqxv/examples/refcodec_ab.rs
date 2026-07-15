//! A/B of the global-reference coder: the order-k model (`encode_blocked`,
//! REF_METHOD_SEQPAR) vs the SPRING-style 2-bit-pack + LZMA coder (`encode_packed`,
//! REF_METHOD_PACK). The reorder container codes the reference once per file and
//! every read is a position on it, so the reference-frame bytes are the *only*
//! thing this codec changes — the seq-stream and whole-archive delta equal the
//! reference-bytes delta measured here. Directly quantifies the #52 lever: how
//! much repeat structure the order-k model was leaving on the table.
//!
//! Runs the container's real pipeline: global minimizer clustering, orientation,
//! `assemble_global`, then the default `merge_reference` pass. Both codings are
//! round-tripped on the merged reference to prove byte-exactness.
//!
//! Usage: cargo run --release -p fqxv --example refcodec_ab -- <fastq> [order]

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use fqxv_reorder::GlobalReference;

const REORDER_K: usize = 15;
const SEQ_ORDER: usize = 11;
const REF_SEQ_BLOCKS: usize = 64; // matches the container's REF_SEQ_BLOCKS

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: refcodec_ab <fastq> [order]");
    let order: usize = args.next().map_or(SEQ_ORDER, |s| s.parse().expect("order"));

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
    eprintln!("{path}: {n} reads, {} seq bases", seq.len());

    let t = Instant::now();
    let plan = fqxv_reorder::plan(&lens, &seq, REORDER_K);
    eprintln!("plan: {:.1}s", t.elapsed().as_secs_f64());

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

    let t = Instant::now();
    let (reference, places) = fqxv_reorder::assemble_global(&refs_all, &cl_anchors);
    eprintln!(
        "assemble_global: {:.1}s ({} contigs)",
        t.elapsed().as_secs_f64(),
        reference.n_contigs()
    );
    let t = Instant::now();
    let (merged, _mplaces) = fqxv_reorder::merge_reference(&refs_all, &reference, &places);
    eprintln!(
        "merge_reference: {:.1}s ({} contigs, {} bases)",
        t.elapsed().as_secs_f64(),
        merged.n_contigs(),
        merged.total_bases(),
    );

    let mb = |b: usize| b as f64 / 1e6;

    // A: shipped block-parallel order-k model.
    let t = Instant::now();
    let seqpar = merged
        .encode_blocked(order, REF_SEQ_BLOCKS)
        .expect("encode_blocked");
    let seqpar_s = t.elapsed().as_secs_f64();
    let rt = GlobalReference::decode_blocked(&seqpar).expect("decode_blocked");
    assert_eq!(rt.raw_bases(), merged.raw_bases(), "seqpar round-trip");
    assert_eq!(rt.contig_lens(), merged.contig_lens(), "seqpar lens");

    let raw = merged.total_bases();
    let bpb = |bytes: usize| 8.0 * bytes as f64 / raw as f64;
    let d = |c: usize| 100.0 * (c as f64 - seqpar.len() as f64) / seqpar.len() as f64;

    println!("\n== global reference coder A/B (order {order}) ==");
    println!(
        "  merged reference : {} contigs, {} bases",
        merged.n_contigs(),
        raw
    );
    println!(
        "  SEQPAR (order-k) : {:>7.2} MB  ({:.3} b/base)  encode {:.1}s   (baseline)",
        mb(seqpar.len()),
        bpb(seqpar.len()),
        seqpar_s
    );

    // SPRING-faithful 2-bit-pack + LZMA (the measured winner).
    let t = Instant::now();
    let packed = merged.encode_packed().expect("encode_packed");
    let packed_s = t.elapsed().as_secs_f64();
    let rt = GlobalReference::decode_packed(&packed).expect("decode_packed");
    assert_eq!(rt.raw_bases(), merged.raw_bases(), "packed round-trip");
    assert_eq!(rt.contig_lens(), merged.contig_lens(), "packed lens");
    println!(
        "  PACK+LZMA (SPRING): {:>7.2} MB  ({:.3} b/base)  encode {:.1}s   {:+.1}%",
        mb(packed.len()),
        bpb(packed.len()),
        packed_s,
        d(packed.len()),
    );
    println!("  (deltas == whole-archive seq delta; all round-trip OK. 2bit+xz(SPRING proxy) ~20.42 MB = -7.2%; raw xz ceiling ~19.87 = -9.7%.)");

    // Optional: dump the raw merged-reference bases so external reference
    // compressors (zstd --long, xz -9e) can re-measure the true LZ headroom vs
    // today's SEQPAR baseline — the memory's -16.6% projection predates the
    // current fqxv_seq coder, so it needs re-checking before investing in a
    // stronger clean-room LZ.
    if let Some(dump) = args.next() {
        std::fs::write(&dump, merged.raw_bases()).expect("dump raw bases");
        eprintln!(
            "wrote {} raw reference bases to {dump}",
            merged.total_bases()
        );
    }
}
