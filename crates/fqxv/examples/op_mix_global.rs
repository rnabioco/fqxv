//! A/B of the reorder sequence codec: block-local literal-rescue (v3) vs the
//! two-pass global-reference codec (v4). Confirms — before any container work —
//! whether the SPRING-style global reference actually shrinks the NovaSeq
//! sequence stream (issue #52), and where v4's bytes go (the `contig_id` stream
//! is the crux that killed naive global assembly).
//!
//! Runs the container's global minimizer clustering (`fqxv_reorder::plan`),
//! orients reads, then:
//!   * v3: per 256Ki-read block, `encode_clustered_rescue` — the block-local
//!     assembler the container ships today (issue baseline 37.21 MB).
//!   * v4: `assemble_global` once over ALL clustered reads → one frozen
//!     reference (coded once via `fqxv_seq`), then `encode_global_block` per
//!     256Ki-read block against that shared reference. The v4 total is the
//!     reference frame PLUS every per-block payload.
//!
//! Both codecs are round-tripped on the real data to prove byte-exactness.
//!
//! Usage: cargo run --release -p fqxv --example op_mix_global -- <fastq> [k]

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use fqxv_reorder::{GlobalReference, Place4};
use rayon::prelude::*;

/// Minimizer k for clustering — matches the container's `REORDER_K`.
const REORDER_K: usize = 15;
/// Reads per block — matches the container's `REORDER_BLOCK_READS`.
const REORDER_BLOCK_READS: usize = 1 << 18;
/// Sequence context order the container passes to the codec by default.
const SEQ_ORDER: usize = 11;

/// Sum the length-prefixed sub-stream sizes of a `[version][n][len][bytes]...`
/// reorder block payload, so the per-stream breakdown can be reported.
fn stream_sizes(payload: &[u8]) -> Vec<usize> {
    let mut pos = 1usize; // skip version byte
    // read varint n (reads)
    read_varint(payload, &mut pos);
    let mut sizes = Vec::new();
    while pos < payload.len() {
        let len = read_varint(payload, &mut pos) as usize;
        sizes.push(len);
        pos += len;
    }
    sizes
}

/// Minimal varint reader (LEB128) matching `fqxv_bytes::write_varint`.
fn read_varint(src: &[u8], pos: &mut usize) -> u64 {
    let mut v = 0u64;
    let mut shift = 0u32;
    loop {
        let b = src[*pos];
        *pos += 1;
        v |= u64::from(b & 0x7f) << shift;
        if b & 0x80 == 0 {
            break;
        }
        shift += 7;
    }
    v
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: op_mix_global <fastq> [k]");
    let k: usize = args.next().map_or(REORDER_K, |s| s.parse().expect("k"));

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

    let t = Instant::now();
    let plan = fqxv_reorder::plan(&lens, &seq, k);
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

    let bsz = REORDER_BLOCK_READS.max(1);
    let ranges: Vec<(usize, usize)> = (0..n).step_by(bsz).map(|s| (s, (s + bsz).min(n))).collect();

    // ---- v3: block-local literal-rescue (the current container codec) --------
    let t = Instant::now();
    let mut v3_bytes = 0usize;
    for &(s, e) in &ranges {
        let e3 =
            fqxv_reorder::encode_clustered_rescue(&refs_all[s..e], &cl_anchors[s..e], SEQ_ORDER)
                .expect("v3 encode");
        // Byte-exactness.
        let dec = fqxv_reorder::decode_clustered_rescue(&e3).expect("v3 decode");
        for (a, b) in dec.iter().zip(refs_all[s..e].iter()) {
            assert_eq!(a.as_slice(), *b, "v3 round-trip mismatch");
        }
        v3_bytes += e3.len();
    }
    let v3_s = t.elapsed().as_secs_f64();

    // ---- v4: two-pass global reference ---------------------------------------
    let t = Instant::now();
    let (reference, places) = fqxv_reorder::assemble_global(&refs_all, &cl_anchors);
    let asm_s = t.elapsed().as_secs_f64();

    // The reference is coded ONCE for the whole file (74% of v4's bytes on
    // NovaSeq), so sweep the hashed high-order tier and keep the smallest — a
    // per-block model can't afford this, but a one-shot whole-file blob can.
    let configs: &[(usize, usize, u32)] = &[
        (SEQ_ORDER, 0, 0),   // plain dense order-11 (baseline)
        (SEQ_ORDER, 13, 25), // the --max non-reorder tier
        (SEQ_ORDER, 16, 26),
        (SEQ_ORDER, 20, 26),
        (SEQ_ORDER, 24, 26),
    ];
    println!("\n== reference coding sweep (order, hash_order, hash_bits) ==");
    let mut best: Option<(Vec<u8>, (usize, usize, u32))> = None;
    let t = Instant::now();
    for &(o, ho, hb) in configs {
        let tc = Instant::now();
        let payload = reference.encode(o, ho, hb).expect("reference encode");
        println!(
            "  o{o} h{ho} b{hb}: {:.2} MB   ({:.1}s)",
            payload.len() as f64 / 1e6,
            tc.elapsed().as_secs_f64()
        );
        if best.as_ref().is_none_or(|(b, _)| payload.len() < b.len()) {
            best = Some((payload, (o, ho, hb)));
        }
    }
    let (ref_payload, ref_cfg) = best.expect("at least one config");
    let ref_s = t.elapsed().as_secs_f64();
    // The chosen reference must round-trip.
    let ref_reload = GlobalReference::decode(&ref_payload).expect("reference decode");
    println!("  -> chosen o{} h{} b{}", ref_cfg.0, ref_cfg.1, ref_cfg.2);

    let t = Instant::now();
    let mut v4_block_bytes = 0usize;
    // Accumulate the per-stream breakdown across blocks (ops,cid,off,slen,nmis,pos,subs,tail).
    let mut agg = vec![0usize; 8];
    for &(s, e) in &ranges {
        let payload =
            fqxv_reorder::encode_global_block(&refs_all[s..e], &places[s..e], &ref_reload)
                .expect("v4 encode");
        let dec = fqxv_reorder::decode_global_block(&payload, &ref_reload).expect("v4 decode");
        for (a, b) in dec.iter().zip(refs_all[s..e].iter()) {
            assert_eq!(a.as_slice(), *b, "v4 round-trip mismatch");
        }
        for (slot, sz) in agg.iter_mut().zip(stream_sizes(&payload)) {
            *slot += sz;
        }
        v4_block_bytes += payload.len();
    }
    let v4_s = t.elapsed().as_secs_f64();
    let v4_bytes = v4_block_bytes + ref_payload.len();

    let mb = |b: usize| b as f64 / 1e6;
    let labels = ["ops", "cid", "off", "slen", "nmis", "pos", "subs", "tail"];
    println!("\n== v4 per-block stream breakdown (summed over blocks) ==");
    for (name, sz) in labels.iter().zip(&agg) {
        println!("  {name:>5}: {:.2} MB", mb(*sz));
    }
    println!(
        "\n== reference frame ==\n  contigs: {}\n  reference bases: {} ({:.1} MB raw)\n  reference payload (fqxv_seq): {:.2} MB",
        reference.n_contigs(),
        reference.total_bases(),
        mb(reference.total_bases()),
        mb(ref_payload.len()),
    );

    let delta = 100.0 * (v4_bytes as f64 - v3_bytes as f64) / v3_bytes.max(1) as f64;
    println!(
        "\n== v3 (block-local rescue) vs v4 (global reference) ==\n  \
         v3 seq bytes : {:.2} MB   (encode {:.1}s)\n  \
         v4 seq bytes : {:.2} MB   = reference {:.2} MB + blocks {:.2} MB\n  \
                        (assemble {:.1}s, ref-code {:.1}s, blocks {:.1}s)\n  \
         delta        : {:+.1}%\n  \
         round-trip (v3 + v4, real data): OK  ({n} reads)",
        mb(v3_bytes),
        v3_s,
        mb(v4_bytes),
        mb(ref_payload.len()),
        mb(v4_block_bytes),
        asm_s,
        ref_s,
        v4_s,
        delta,
    );

    // ---- overlap-merge refinement (A): threshold sweep ----------------------
    //
    // The greedy reference fragments into many short contigs (~1.4 reads each);
    // `merge_reference` chains contigs whose suffix overlaps another's prefix into
    // longer super-contigs, storing shared sequence once and remapping the reads.
    // A single pass captures ~all the gain (iterating to convergence adds <0.3%),
    // but it merges only ~18% of contigs — the other 82% never find a qualifying
    // successor. This sweeps the overlap-search thresholds to tell WHY: if looser
    // criteria (shorter min overlap, higher mismatch budget, wider prefix window)
    // merge substantially more with the ratio still improving, the default
    // thresholds are the limiter; if not, the remaining contigs are genuinely
    // distinct sequence and need a different mechanism (containment absorption).
    // Each config is ONE pass, format-transparent, round-tripped on real data.
    let global_ref_plain = reference
        .encode(SEQ_ORDER, 0, 0)
        .expect("plain ref enc")
        .len();
    let global_total_plain = global_ref_plain + v4_block_bytes;
    let default_cfg = fqxv_reorder::MergeConfig::default();
    let sweep: [(&str, fqxv_reorder::MergeConfig); 6] = [
        ("default    ", default_cfg),
        (
            "ovl16      ",
            fqxv_reorder::MergeConfig {
                min_ovl: 16,
                ..default_cfg
            },
        ),
        (
            "mism/5     ",
            fqxv_reorder::MergeConfig {
                mism_div: 5,
                ..default_cfg
            },
        ),
        (
            "prefix128  ",
            fqxv_reorder::MergeConfig {
                prefix: 128,
                ..default_cfg
            },
        ),
        (
            "fanout32   ",
            fqxv_reorder::MergeConfig {
                fanout: 32,
                ..default_cfg
            },
        ),
        (
            "loose-all  ",
            fqxv_reorder::MergeConfig {
                min_ovl: 16,
                mism_div: 5,
                prefix: 128,
                fanout: 32,
                ..default_cfg
            },
        ),
    ];
    println!("\n== overlap-merge refinement (A): single-pass threshold sweep ==");
    println!("  config        contigs     %merged    refMB    blkMB   v4totMB   vs-v3    time(s)");
    println!(
        "  {:12}  {:>7}   {:>7}   {:>8.2} {:>7.2} {:>8.2}   {:+5.1}%   {:>6}",
        "none(v4)",
        reference.n_contigs(),
        "-",
        mb(global_ref_plain),
        mb(v4_block_bytes),
        mb(global_total_plain),
        100.0 * (global_total_plain as f64 - v3_bytes as f64) / v3_bytes.max(1) as f64,
        "-",
    );
    let base_contigs = reference.n_contigs();
    for (name, cfg) in sweep {
        let t = Instant::now();
        let (merged, mplaces) =
            fqxv_reorder::merge_reference_with(&refs_all, &reference, &places, cfg);
        let merge_s = t.elapsed().as_secs_f64();
        let mref_payload = merged.encode(SEQ_ORDER, 0, 0).expect("merged ref enc");
        let mref = GlobalReference::decode(&mref_payload).expect("merged ref dec");
        let t = Instant::now();
        let mut merged_block_bytes = 0usize;
        for &(s, e) in &ranges {
            let payload = fqxv_reorder::encode_global_block(&refs_all[s..e], &mplaces[s..e], &mref)
                .expect("merged v4 enc");
            let dec = fqxv_reorder::decode_global_block(&payload, &mref).expect("merged v4 dec");
            for (a, b) in dec.iter().zip(refs_all[s..e].iter()) {
                assert_eq!(a.as_slice(), *b, "merged round-trip mismatch");
            }
            merged_block_bytes += payload.len();
        }
        let code_s = t.elapsed().as_secs_f64();
        let merged_total = merged_block_bytes + mref_payload.len();
        let contigs = merged.n_contigs();
        let pct_merged = 100.0 * (base_contigs - contigs) as f64 / base_contigs.max(1) as f64;
        println!(
            "  {name}  {:>7}   {:>6.1}%   {:>8.2} {:>7.2} {:>8.2}   {:+5.1}%   {:>6.1}",
            contigs,
            pct_merged,
            mb(mref_payload.len()),
            mb(merged_block_bytes),
            mb(merged_total),
            100.0 * (merged_total as f64 - v3_bytes as f64) / v3_bytes.max(1) as f64,
            merge_s + code_s,
        );
    }

    // ---- windowed parallel assembly sweep -----------------------------------
    //
    // Pass-1 global assembly is an inherently-sequential fold (each read placed
    // against the growing reference), so it is the compress-time floor. Splitting
    // the clustered reads into W non-overlapping windows lets the windows assemble
    // CONCURRENTLY across cores — and each window's smaller working set stays in
    // cache — at the cost of losing cross-WINDOW overlaps (a read near a window
    // edge can't attach to a contig in another window). W=1 is the current global
    // codec. This measures the ratio/speed trade so the container can pick W.
    let window_counts = [1usize, 2, 4, 8, 16, 32];
    println!("\n== windowed parallel assembly (W windows, assembled concurrently) ==");
    println!("   W   windows   contigs      refMB    blkMB    totMB     Δ-ratio    wall(s)");
    let mut global_tot = 0f64;
    for (wi, &w) in window_counts.iter().enumerate() {
        let wsz = n.div_ceil(w);
        let wins: Vec<(usize, usize)> = (0..w)
            .map(|i| (i * wsz, ((i + 1) * wsz).min(n)))
            .filter(|(s, e)| s < e)
            .collect();
        let t = Instant::now();
        // Each window is independent: assemble it, code its reference once, then
        // code its reads in 256Ki blocks against that reference (round-tripping).
        let per_window: Vec<(usize, usize, usize)> = wins
            .par_iter()
            .map(|&(ws, we)| {
                let wrefs = &refs_all[ws..we];
                let wanch = &cl_anchors[ws..we];
                let (reference, places) = fqxv_reorder::assemble_global(wrefs, wanch);
                let ref_payload = reference.encode(SEQ_ORDER, 0, 0).expect("ref enc");
                let refl = GlobalReference::decode(&ref_payload).expect("ref dec");
                let m = we - ws;
                let mut blk_bytes = 0usize;
                let mut s = 0usize;
                while s < m {
                    let e = (s + REORDER_BLOCK_READS).min(m);
                    let payload =
                        fqxv_reorder::encode_global_block(&wrefs[s..e], &places[s..e], &refl)
                            .expect("v4 enc");
                    let dec = fqxv_reorder::decode_global_block(&payload, &refl).expect("v4 dec");
                    for (a, b) in dec.iter().zip(wrefs[s..e].iter()) {
                        assert_eq!(a.as_slice(), *b, "windowed v4 round-trip mismatch");
                    }
                    blk_bytes += payload.len();
                    s = e;
                }
                (reference.n_contigs(), ref_payload.len(), blk_bytes)
            })
            .collect();
        let wall = t.elapsed().as_secs_f64();
        let contigs: usize = per_window.iter().map(|x| x.0).sum();
        let ref_b: usize = per_window.iter().map(|x| x.1).sum();
        let blk_b: usize = per_window.iter().map(|x| x.2).sum();
        let tot = mb(ref_b + blk_b);
        if wi == 0 {
            global_tot = tot;
        }
        let dratio = 100.0 * (tot - global_tot) / global_tot.max(1e-9);
        println!(
            "  {w:>3}   {:>7}   {contigs:>7}   {:>8.2} {:>8.2} {:>8.2}   {:+7.2}%   {wall:>7.1}",
            wins.len(),
            mb(ref_b),
            mb(blk_b),
            tot,
            dratio,
        );
    }
    // Keep `Place4` import meaningful even if the compiler prunes it otherwise.
    let _ = std::mem::size_of::<Place4>();
}
