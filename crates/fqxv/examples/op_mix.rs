//! A/B of the reorder sequence codec: single-contig (v2) vs literal-rescue (v3).
//!
//! Runs the container's global minimizer clustering (`fqxv_reorder::plan`),
//! orients reads, then per 256Ki-read block (the container's
//! `REORDER_BLOCK_READS`) measures, for BOTH codecs:
//!
//!   * the op-mix — MATCH (exact dup, free), CONTIG (offset + a few mismatches),
//!     LITERAL (context-coded from scratch), plus the share of bases coded from
//!     scratch (literals + novel tails), and
//!   * the actual compressed sequence-stream bytes (the real ratio signal).
//!
//! The rescue codec is also round-tripped on the real data to confirm it is
//! byte-exact. A high LITERAL / from-scratch share under v2 that drops under v3
//! — and a smaller v3 byte total — is the payoff of re-attaching would-be
//! literals to earlier contigs (the SPRING-assembly lever).
//!
//! Usage: cargo run --release -p fqxv --example op_mix -- <fastq> [k]

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use fqxv_reorder::OpStats;

/// Minimizer k for clustering — matches the container's `REORDER_K`.
const REORDER_K: usize = 15;
/// Reads per block — matches the container's `REORDER_BLOCK_READS`.
const REORDER_BLOCK_READS: usize = 1 << 18;
/// Sequence context order the container passes to the codec by default.
const SEQ_ORDER: usize = 11;

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
    let total_bases = seq.len();
    eprintln!("{path}: {n} reads, {total_bases} seq bases, clustering k={k}");

    // Global minimizer clustering — the whole-file plan the container computes.
    let t = Instant::now();
    let plan = fqxv_reorder::plan(&lens, &seq, k);
    eprintln!("plan: {:.1}s", t.elapsed().as_secs_f64());

    // Cumulative offsets into the concatenated seq.
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in &lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // Clustered, oriented reads + their minimizer anchors, exactly as the
    // container builds them before calling the codec.
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

    // Per-block: op-mix + compressed bytes for both codecs, round-trip v3.
    let bsz = REORDER_BLOCK_READS.max(1);
    let mut v2 = OpStats::default();
    let mut v3 = OpStats::default();
    let mut v2_bytes = 0usize;
    let mut v3_bytes = 0usize;
    let mut enc_v2_s = 0.0f64;
    let mut enc_v3_s = 0.0f64;
    let mut s = 0usize;
    while s < n {
        let e = (s + bsz).min(n);
        let refs: Vec<&[u8]> = cl_reads[s..e].iter().map(Vec::as_slice).collect();
        let anch = &cl_anchors[s..e];

        v2.merge(&fqxv_reorder::op_stats(&refs, anch));
        v3.merge(&fqxv_reorder::op_stats_rescue(&refs, anch));

        let ta = Instant::now();
        let e2 = fqxv_reorder::encode_clustered(&refs, anch, SEQ_ORDER).expect("v2 encode");
        enc_v2_s += ta.elapsed().as_secs_f64();
        let tb = Instant::now();
        let e3 = fqxv_reorder::encode_clustered_rescue(&refs, anch, SEQ_ORDER).expect("v3 encode");
        enc_v3_s += tb.elapsed().as_secs_f64();
        v2_bytes += e2.len();
        v3_bytes += e3.len();

        // Byte-exactness on real data.
        let dec = fqxv_reorder::decode_clustered_rescue(&e3).expect("v3 decode");
        assert_eq!(dec.len(), refs.len(), "v3 read count");
        for (a, b) in dec.iter().zip(refs.iter()) {
            assert_eq!(a.as_slice(), *b, "v3 round-trip mismatch");
        }
        s = e;
    }

    let pf = |x: u64, tot: u64| 100.0 * x as f64 / tot.max(1) as f64;
    let report = |tag: &str, st: &OpStats, bytes: usize, enc_s: f64| {
        let scratch = st.literal_bases + st.novel_tail_bases;
        let reads = st.reads.max(1) as u64;
        let cts = st.contigs.max(1) as f64;
        println!(
            "\n== {tag} ==\n  \
             reads : MATCH {:.1}%  CONTIG {:.1}%  LITERAL {:.1}%\n  \
             bases : MATCH {:.1}%  overlap {:.1}%  novel-tail {:.1}%  LITERAL {:.1}%\n  \
             from-scratch bases (LITERAL + novel tail): {:.1}%\n  \
             CONTIG: {:.3} mismatch/read, {:.1} overlap/read\n  \
             seq bytes: {} ({:.1} MB)  |  encode {:.1}s",
            pf(st.matches as u64, reads),
            pf(st.contigs as u64, reads),
            pf(st.literals as u64, reads),
            pf(st.match_bases, st.total_bases),
            pf(st.contig_overlap_bases, st.total_bases),
            pf(st.novel_tail_bases, st.total_bases),
            pf(st.literal_bases, st.total_bases),
            pf(scratch, st.total_bases),
            st.contig_mismatches as f64 / cts,
            st.contig_overlap_bases as f64 / cts,
            bytes,
            bytes as f64 / 1e6,
            enc_s,
        );
    };

    report("v2  single-contig", &v2, v2_bytes, enc_v2_s);
    report("v3  literal-rescue", &v3, v3_bytes, enc_v3_s);

    let lit2 = pf(v2.literals as u64, v2.reads.max(1) as u64);
    let lit3 = pf(v3.literals as u64, v3.reads.max(1) as u64);
    let sc2 = pf(v2.literal_bases + v2.novel_tail_bases, v2.total_bases);
    let sc3 = pf(v3.literal_bases + v3.novel_tail_bases, v3.total_bases);
    let byte_delta = 100.0 * (v3_bytes as f64 - v2_bytes as f64) / v2_bytes.max(1) as f64;
    println!(
        "\n== delta (v2 -> v3) ==\n  \
         LITERAL reads      : {:.1}% -> {:.1}%\n  \
         from-scratch bases : {:.1}% -> {:.1}%\n  \
         seq bytes          : {:.1} MB -> {:.1} MB  ({:+.1}%)\n  \
         round-trip (v3, real data): OK  ({} reads)",
        lit2,
        lit3,
        sc2,
        sc3,
        v2_bytes as f64 / 1e6,
        v3_bytes as f64 / 1e6,
        byte_delta,
        n,
    );
}
