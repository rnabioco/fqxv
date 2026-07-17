//! Issue #102: does `fqxv-reorder` strand Illumina reads as literals *because it
//! has no indel op*?
//!
//! `fqxv-reorder` codes a read against a contig with `OP_CONTIG (offdelta, slen,
//! nmis, pos, subs)` — substitutions only. There is no insert or delete op
//! anywhere in the crate, and `try_place` is an ungapped compare. So a single
//! indel shifts every base after it, the read mismatches its way past the
//! `overlap / 4` budget, and an otherwise-perfect neighbour is stranded as a
//! LITERAL and coded from scratch. The cost of hitting one is total.
//!
//! #102 asks for the cheapest experiment that settles whether that matters on
//! real Illumina data, rather than porting the long-read overlap stack on the
//! strength of it working for long reads. This is that experiment: run the
//! container's real clustering, run the real assembler, and for every read it
//! strands, re-test the SAME candidates it considered with a small band. A read
//! counts only if a gapped compare lands inside the same budget the ungapped one
//! blew, using at least one indel — i.e. the missing indel op is the reason.
//!
//! If `indel_rescuable` is ~0, close #102: there is nothing to win and the
//! machinery is not worth the complexity. If it is material, the fix is narrow —
//! give the short-read codec an indel op, reusing `align_banded`.
//!
//! Usage: cargo run --release -p fqxv --example indel_probe -- <fastq> [band]

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use fqxv_reorder::IndelProbe;

/// Minimizer k for clustering — matches the container's `REORDER_K`.
const REORDER_K: usize = 15;
/// Reads per block — matches the container's `REORDER_BLOCK_READS`.
const REORDER_BLOCK_READS: usize = 1 << 18;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: indel_probe <fastq> [band]");
    let band: usize = args.next().map_or(4, |s| s.parse().expect("band"));

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
    eprintln!(
        "{path}: {n} reads, {} seq bases, clustering k={REORDER_K}, band={band}",
        seq.len()
    );

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

    // Clustered, oriented reads + anchors — exactly what the container hands the
    // codec. The probe has to see the codec's own input or it measures something
    // else.
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

    let t = Instant::now();
    let mut tot = IndelProbe::default();
    let bsz = REORDER_BLOCK_READS.max(1);
    let mut s = 0usize;
    while s < n {
        let e = (s + bsz).min(n);
        let refs: Vec<&[u8]> = cl_reads[s..e].iter().map(Vec::as_slice).collect();
        let st = fqxv_reorder::indel_probe(&refs, &cl_anchors[s..e], band);
        tot.reads += st.reads;
        tot.literals += st.literals;
        tot.indel_rescuable += st.indel_rescuable;
        tot.truly_novel += st.truly_novel;
        tot.rescuable_bases += st.rescuable_bases;
        tot.literal_bases += st.literal_bases;
        s = e;
    }
    eprintln!("probe: {:.1}s", t.elapsed().as_secs_f64());

    let pct = |a: u64, b: u64| {
        if b == 0 {
            0.0
        } else {
            a as f64 * 100.0 / b as f64
        }
    };
    println!("--- issue #102: indel-driven literals ---");
    println!("  reads placed or stranded   {}", tot.reads);
    println!(
        "  LITERAL                    {} ({:.2}% of reads, {} bases)",
        tot.literals,
        pct(tot.literals, tot.reads),
        tot.literal_bases
    );
    println!(
        "  ...rescuable by an indel   {} ({:.2}% of literals, {:.4}% of reads)",
        tot.indel_rescuable,
        pct(tot.indel_rescuable, tot.literals),
        pct(tot.indel_rescuable, tot.reads)
    );
    println!(
        "  ...genuinely novel         {} ({:.2}% of literals)",
        tot.truly_novel,
        pct(tot.truly_novel, tot.literals)
    );
    println!(
        "  bases an indel op reclaims {} ({:.4}% of all literal bases)",
        tot.rescuable_bases,
        pct(tot.rescuable_bases, tot.literal_bases)
    );
    println!(
        "\n  Verdict: an indel op would move {:.4}% of reads off the literal path.",
        pct(tot.indel_rescuable, tot.reads)
    );
}
