//! TEMPORARY A/B harness (remove with `assemble_global_rightward`): compares the
//! rightward-only assembler against the bidirectional one on a FASTQ file, both
//! before and after the overlap-merge, reporting contigs / reference bases /
//! encoded reference-frame bytes. Answers: does leftward extension still reduce
//! the reference AFTER `merge_reference`, or does the merge subsume it?
//!
//! Usage: `cargo run --release -p fqxv-reorder --example bidir_ab -- <file.fastq>`

use std::io::{BufRead, BufReader};

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: bidir_ab <file.fastq>");
    let f = std::fs::File::open(&path).expect("open");
    let mut rdr = BufReader::new(f);

    // Minimal 4-line FASTQ reader → concatenated seq + lengths.
    let (mut seq, mut lens) = (Vec::<u8>::new(), Vec::<u32>::new());
    let mut line = String::new();
    let mut li = 0u64;
    loop {
        line.clear();
        if rdr.read_line(&mut line).expect("read") == 0 {
            break;
        }
        if li % 4 == 1 {
            let s = line.trim_end();
            seq.extend_from_slice(s.as_bytes());
            lens.push(s.len() as u32);
        }
        li += 1;
    }
    eprintln!("reads: {}", lens.len());

    let p = fqxv_reorder::plan(&lens, &seq, fqxv_reorder::DEFAULT_K);
    let mut offs = vec![0usize];
    for &l in &lens {
        offs.push(offs.last().unwrap() + l as usize);
    }
    let cl: Vec<Vec<u8>> = p
        .order
        .iter()
        .map(|&oi| {
            let oi = oi as usize;
            let s = &seq[offs[oi]..offs[oi + 1]];
            if p.flip[oi] {
                fqxv_reorder::revcomp(s)
            } else {
                s.to_vec()
            }
        })
        .collect();
    let refs: Vec<&[u8]> = cl.iter().map(Vec::as_slice).collect();
    let anchors: Vec<u32> = p.order.iter().map(|&oi| p.anchor[oi as usize]).collect();

    for (tag, assemble) in [
        (
            "rightward",
            fqxv_reorder::assemble_global_rightward
                as fn(
                    &[&[u8]],
                    &[u32],
                )
                    -> (fqxv_reorder::GlobalReference, Vec<fqxv_reorder::Place4>),
        ),
        ("bidir", fqxv_reorder::assemble_global),
    ] {
        let (reference, places) = assemble(&refs, &anchors);
        let pre_c = reference.n_contigs();
        let pre_b = reference.total_bases();
        let pre_enc = reference.encode(11, 0, 0).expect("enc").len();
        let (merged, mp) = fqxv_reorder::merge_reference(&refs, &reference, &places);
        let post_c = merged.n_contigs();
        let post_b = merged.total_bases();
        let post_enc = merged.encode(11, 0, 0).expect("enc").len();
        let post_pack = merged.encode_packed().expect("pack").len();
        // Per-read placement cost against the merged reference (one block).
        let block = fqxv_reorder::encode_global_block(&refs, &mp, &merged)
            .expect("block")
            .len();
        println!(
            "{tag:>10}  pre: {pre_c:>8} contigs {pre_b:>10} bases {pre_enc:>10} enc  |  \
             post-merge: {post_c:>8} contigs {post_b:>10} bases  ref[order11 {post_enc} / pack {post_pack}]  block {block}  TOTAL {}",
            post_pack + block
        );
    }
}
