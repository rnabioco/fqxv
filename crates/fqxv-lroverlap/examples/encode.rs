//! End-to-end: reads -> overlaps -> layout -> consensus -> edit scripts ->
//! entropy-coded streams, reporting the achieved DNA-stream bits/base.
//!
//! This closes the loop. M1b measured 0.0401 bits/base on HiFi using the TRUE
//! reference and minimap2's placement — an oracle. This runs the same
//! measurement on the crate's own assembly and its own placement, so the number
//! it prints is one the codec can actually deliver.
//!
//! Targets: HiFi <= 0.068 beats CoLoRd; ~0.040 matches the oracle bound.
//!
//! ```text
//! cargo run --release --example encode -- reads.fastq [ont|hifi]
//! ```

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use rayon::prelude::*;

use fqxv_lroverlap::{
    align_banded, consensus, find_overlaps, layout, place_against, Anchored, ChainOpts,
    ConsensusOpts, Index, LayoutOpts, Op, Repeat, Sketch,
};

fn read_fastq(path: &str) -> (Vec<u32>, Vec<u8>) {
    let f = BufReader::with_capacity(1 << 20, File::open(path).expect("open fastq"));
    let (mut lens, mut seq) = (Vec::new(), Vec::new());
    for (i, line) in f.lines().enumerate() {
        if i % 4 != 1 {
            continue;
        }
        let l = line.expect("line");
        lens.push(l.len() as u32);
        seq.extend_from_slice(l.as_bytes());
    }
    (lens, seq)
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

fn push_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let b = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(b);
            break;
        }
        out.push(b | 0x80);
    }
}

fn base_code(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// The streams a real codec would emit, one per kind so each gets its own model.
#[derive(Default)]
struct Streams {
    ops: Vec<u8>,
    runs: Vec<u8>,
    subs: Vec<u8>,
    ins_bases: Vec<u8>,
    indel_lens: Vec<u8>,
    placements: Vec<u8>,
    /// One byte per placed read: does it align reverse-complemented to the
    /// reference. Its own stream because it is one bit of real information that
    /// the decoder cannot derive, and rANS takes it back down to ~that. The
    /// layout-framed encode could skip it (reads were pre-oriented by the
    /// layout); placing raw reads against a reference cannot.
    flips: Vec<u8>,
    /// Bases of reads that could not be placed — coded standalone.
    literals: Vec<u8>,
}

impl Streams {
    fn merge(&mut self, o: Streams) {
        self.ops.extend(o.ops);
        self.runs.extend(o.runs);
        self.subs.extend(o.subs);
        self.ins_bases.extend(o.ins_bases);
        self.indel_lens.extend(o.indel_lens);
        self.placements.extend(o.placements);
        self.flips.extend(o.flips);
        self.literals.extend(o.literals);
    }
}

/// Entropy-code one stream, returning its compressed size.
fn encoded_len(buf: &[u8]) -> usize {
    if buf.is_empty() {
        return 0;
    }
    fqxv_rans::encode(buf, fqxv_rans::Order::One)
        .expect("rans")
        .len()
}

fn report(label: &str, raw: usize, enc: usize, total_bases: u64) {
    if raw == 0 {
        println!("    {label:<12}      empty");
        return;
    }
    println!(
        "    {label:<12} {raw:>12} B raw -> {enc:>11} B  ({:.3} b/sym, {:.4} bits/base)",
        (enc as f64 * 8.0) / raw as f64,
        (enc as f64 * 8.0) / total_bases as f64
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: encode <reads.fastq> [ont|hifi]");
        std::process::exit(2);
    }
    let sketch = match args.get(2).map(String::as_str) {
        Some("hifi") => Sketch::hifi(),
        _ => Sketch::ont(),
    };

    let (lens, seq) = read_fastq(&args[1]);
    let total_bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    println!(
        "reads {} · bases {} · sketch w={} k={}",
        lens.len(),
        total_bases,
        sketch.w,
        sketch.k
    );

    let mut offs = Vec::with_capacity(lens.len() + 1);
    let mut acc = 0usize;
    for &l in &lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);
    let read_at = |r: usize| &seq[offs[r]..offs[r + 1]];

    // ---- SUBSAMPLE FOR LAYOUT ------------------------------------------
    // The layout is fed a SUBSAMPLE, not every read. Two problems have the same
    // answer. (1) `overlaps` is quadratic in coverage: 636.8s of a 770s run at
    // 300x, already saturating 64 cores. (2) The layout starves at depth — at
    // 300x a read's ~800 partners are mostly COINCIDENT reads, which carry no
    // layout information yet crowd out the few that extend it, shattering the
    // layout into 43287 contigs where 40x gives ONE (see `layout.rs`).
    //
    // So do not fix the layout for 300x — feed it 40x, where it already works,
    // and place the other 260x against the consensus it produces.
    // `place_against` is one hop from a fixed frame and never reads a layout
    // offset, so the deep vote and the encode are untouched. The reference is
    // then amortised over every read rather than over the subsample, which is
    // where most of the predicted win comes from.
    //
    // A stride, not a coverage target: the genome size is not ours to assume.
    // FQXV_LAYOUT_STRIDE=8 takes ecoli_hifi's 300x down to ~40x. Default 1 is
    // the old behaviour exactly.
    let stride: usize = env::var("FQXV_LAYOUT_STRIDE")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&s| s >= 1)
        .unwrap_or(1);
    let sub: Vec<u32> = (0..lens.len() as u32).step_by(stride).collect();
    let sub_lens: Vec<u32> = sub.iter().map(|&r| lens[r as usize]).collect();
    let mut sub_seq: Vec<u8> = Vec::new();
    for &r in &sub {
        sub_seq.extend_from_slice(read_at(r as usize));
    }
    let mut sub_offs = Vec::with_capacity(sub.len() + 1);
    let mut acc = 0usize;
    for &l in &sub_lens {
        sub_offs.push(acc);
        acc += l as usize;
    }
    sub_offs.push(acc);
    let sub_read_at = |i: usize| &sub_seq[sub_offs[i]..sub_offs[i + 1]];
    if stride > 1 {
        let sub_bases: u64 = sub_lens.iter().map(|&l| u64::from(l)).sum();
        println!(
            "layout subsample: {} of {} reads (stride {}) · {:.0} of {:.0} Mbase",
            sub.len(),
            lens.len(),
            stride,
            sub_bases as f64 / 1e6,
            total_bases as f64 / 1e6
        );
    }

    let t = Instant::now();
    let idx = Index::build(&sub_lens, &sub_seq, sketch, Repeat::default()).expect("index");
    println!(
        "index: {} occs · {:.1}s",
        idx.len(),
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let ovs: Vec<Vec<_>> = (0..sub.len())
        .into_par_iter()
        .map(|r| find_overlaps(&idx, r as u32, sub_read_at(r), ChainOpts::default()))
        .collect();
    println!("overlaps: {:.1}s", t.elapsed().as_secs_f64());

    let t = Instant::now();
    let contigs = layout(&sub_lens, &ovs, LayoutOpts::default());
    let placed: usize = contigs
        .iter()
        .filter(|c| c.reads.len() > 1)
        .map(|c| c.reads.len())
        .sum();
    let singletons = contigs.iter().filter(|c| c.reads.len() == 1).count();
    // SPAN, not just count. A contig count alone is not a correctness check: a
    // layout that collapses distinct loci onto each other reports FEWER contigs
    // while getting worse. The number that catches both failure modes is total
    // span against the genome -- fragmentation inflates it (43.9 Mb on a 5.16 Mb
    // genome), collapse deflates it (1.5 Mb). Only ~1x is right.
    let mut spans: Vec<u64> = contigs
        .iter()
        .filter(|c| c.reads.len() > 1)
        .map(|c| u64::from(c.len))
        .collect();
    spans.sort_unstable_by(|a, b| b.cmp(a));
    let total_span: u64 = spans.iter().sum();
    let half = total_span / 2;
    let mut acc = 0u64;
    let n50 = spans
        .iter()
        .find(|&&s| {
            acc += s;
            acc >= half
        })
        .copied()
        .unwrap_or(0);
    println!(
        "layout: {} contigs · {} reads on multi-read contigs · {} singletons · {:.1}s",
        contigs.len(),
        placed,
        singletons,
        t.elapsed().as_secs_f64()
    );
    println!(
        "layout span: {:.2} Mb total · largest {:.2} Mb · N50 {:.2} Mb  (1x genome is the target)",
        total_span as f64 / 1e6,
        spans.first().copied().unwrap_or(0) as f64 / 1e6,
        n50 as f64 / 1e6,
    );

    // ---- CONSENSUS over the subsample ----------------------------------
    // Only multi-read contigs get a reference; a singleton is one read and has
    // nothing to vote with.
    let t = Instant::now();
    let consensi: Vec<Vec<u8>> = contigs
        .par_iter()
        .filter(|c| c.reads.len() > 1)
        .filter_map(|c| {
            // Orient reads as the layout placed them, straight into the table
            // `consensus` wants — which it indexes by read id, so it is
            // full-length over the SUBSAMPLE's local ids, the ids the layout
            // speaks. Built in one pass rather than oriented into a Vec and then
            // cloned into this one: that held two full copies of every read on
            // the contig, times however many contigs rayon ran at once.
            let mut by_id: Vec<Vec<u8>> = vec![Vec::new(); sub.len()];
            for p in &c.reads {
                let r = sub_read_at(p.read as usize);
                by_id[p.read as usize] = if p.flip { revcomp(r) } else { r.to_vec() };
            }
            // The sketch must match the platform: the draft is built and voted
            // by chaining, so a sketch too sparse for the error rate loses reads
            // from both.
            let cons = consensus(
                c,
                &by_id,
                ConsensusOpts {
                    sketch,
                    ..ConsensusOpts::default()
                },
            );
            // Dump the consensus so its quality can be measured directly against
            // the true genome, rather than inferred from the bits/base it
            // produces. FQXV_DUMP_CONS=<dir> writes one FASTA per contig.
            if let Ok(dir) = std::env::var("FQXV_DUMP_CONS") {
                use std::io::Write;
                if let Ok(mut f) = File::create(format!("{dir}/contig{}.fa", c.id)) {
                    let _ = writeln!(f, ">contig{} reads={}", c.id, c.reads.len());
                    for chunk in cons.seq.chunks(80) {
                        let _ = f.write_all(chunk);
                        let _ = f.write_all(b"\n");
                    }
                }
            }
            if cons.is_empty() {
                None
            } else {
                Some(cons.seq)
            }
        })
        .collect();
    let ref_bases: usize = consensi.iter().map(Vec::len).sum();
    println!(
        "consensus: {} references · {:.2} Mb · {:.1}s",
        consensi.len(),
        ref_bases as f64 / 1e6,
        t.elapsed().as_secs_f64()
    );

    // ---- PLACE EVERY READ ----------------------------------------------
    // Including the ones the layout never saw. This is the whole point of the
    // subsample: the reference is built from 1/stride of the data and paid for
    // by all of it.
    let t = Instant::now();
    let all_refs: Vec<&[u8]> = (0..lens.len()).map(read_at).collect();
    let placed: Vec<Vec<Option<Anchored>>> = consensi
        .iter()
        .map(|cs| place_against(cs, &all_refs, sketch, ChainOpts::default()))
        .collect();
    // Each read goes to its best-scoring reference. Score DESC, then reference
    // index, so the assignment cannot depend on iteration order.
    let mut best: Vec<Option<(usize, Anchored)>> = vec![None; lens.len()];
    for (ci, per_read) in placed.iter().enumerate() {
        for (r, a) in per_read.iter().enumerate() {
            let Some(a) = a else { continue };
            let take = match best[r] {
                None => true,
                Some((_, b)) => a.score > b.score,
            };
            if take {
                best[r] = Some((ci, *a));
            }
        }
    }
    let n_placed = best.iter().filter(|b| b.is_some()).count();
    println!(
        "placement: {} of {} reads on a reference · {:.1}s",
        n_placed,
        lens.len(),
        t.elapsed().as_secs_f64()
    );

    // ---- CODE EACH READ AGAINST ITS REFERENCE --------------------------
    let t = Instant::now();
    // Group by reference, ordered by offset then read: `placements` is delta-
    // coded on the start, so offset order is what makes it small, and the
    // tie-break keeps the order total.
    let mut by_ref: Vec<Vec<(u32, u32)>> = vec![Vec::new(); consensi.len()];
    for (r, b) in best.iter().enumerate() {
        if let Some((ci, a)) = b {
            by_ref[*ci].push((a.offset, r as u32));
        }
    }
    for v in &mut by_ref {
        v.sort_unstable();
    }

    // The band the read-vs-reference alignment is allowed.
    //
    // This was 96, and 96 was throttling the codec. Ground truth (minimap2) puts
    // a read at 0.00374 edits/base against these consensi — 0.00235 of that is
    // the read's own error against the true genome, the rest the consensus's
    // 0.00155 — while the codec was emitting 0.0078. A banded DP cannot find a
    // path that leaves the band, so it was paying for edits that are not there.
    // Every one of them also splits a match run, which is why `match-runs` is
    // half the archive. Measured, ecoli_hifi 300x stride 8:
    //
    //     band    ops raw    bits/base   scripts
    //       96      22.8M       0.0835       94s
    //      384      17.2M       0.0704     ~350s
    //      768      17.0M       0.0699      717s
    //
    // It plateaus at ~384: the drift a read needs is real and BOUNDED, not
    // runaway placement error (which would keep improving). 768 buys 0.7% for
    // double the work; 384 buys 16% for four times it, which an archiver should
    // take — compression is paid once.
    //
    // 384 is still a workaround. `align_banded` asks callers to "size `band` from
    // the chain's observed diagonal drift", and this hardcodes a constant instead
    // — it cannot do better, because `Anchored` carries an offset but not the
    // drift the chain saw. Widening it globally pays that drift's WORST case on
    // every read. Threading the per-read drift out of `place_against` should get
    // 0.0699 nearer to band-96 cost.
    // Per-read: the chain's own drift plus a margin for the unanchored ends,
    // capped so one pathological chain cannot allocate a quadratic table. A
    // fixed band charges every read the tail's worst case; this charges each
    // read what its own chain says it needs. FQXV_BAND overrides with a constant
    // for A/B.
    let band_margin: usize = env::var("FQXV_BAND_MARGIN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(96);
    let band_fixed: Option<usize> = env::var("FQXV_BAND").ok().and_then(|v| v.parse().ok());
    let band_cap: usize = env::var("FQXV_BAND_CAP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2048);
    let mut per_contig: Vec<Streams> = consensi
        .par_iter()
        .enumerate()
        .map(|(ci, cs)| {
            let mut s = Streams::default();
            // Align in PARALLEL, then emit serially. The alignment is ~10 ms per
            // read (a 12.9 kb read at band 96 is ~2.5M cells) and there are
            // thousands per reference; the emission is microseconds. Doing both
            // in one loop made the whole thing serial, because delta-coded
            // placements depend on the previous read's start — so the cheap
            // order-dependent step was forcing the expensive order-independent
            // one to run on a single core.
            let aligned: Vec<Option<(usize, Vec<Op>)>> = by_ref[ci]
                .par_iter()
                .map(|&(_, r)| {
                    let a = best[r as usize].expect("assigned").1;
                    let raw = read_at(r as usize);
                    let read = if a.flip { revcomp(raw) } else { raw.to_vec() };
                    let start = a.offset as usize;
                    let end = (start + read.len() + 64).min(cs.len());
                    if start >= end {
                        return None;
                    }
                    let band = band_fixed
                        .unwrap_or_else(|| a.drift as usize + band_margin)
                        .min(band_cap);
                    Some((start, align_banded(&cs[start..end], &read, band).ops))
                })
                .collect();

            let mut prev_start: i64 = 0;
            for (idx, slot) in aligned.into_iter().enumerate() {
                let r = by_ref[ci][idx].1 as usize;
                let Some((start, ops)) = slot else {
                    // Unplaceable -> code standalone.
                    s.literals.extend(read_at(r).iter().map(|&b| base_code(b)));
                    continue;
                };
                s.flips.push(u8::from(best[r].expect("assigned").1.flip));
                // Placement: delta-coded start. Order-dependent, hence serial.
                let d = start as i64 - prev_start;
                prev_start = start as i64;
                push_varint(&mut s.placements, ((d << 1) ^ (d >> 63)) as u64);

                for op in &ops {
                    match op {
                        Op::Match(n) => {
                            s.ops.push(0);
                            push_varint(&mut s.runs, u64::from(*n));
                        }
                        Op::Sub(b) => {
                            s.ops.push(1);
                            s.subs.push(base_code(*b));
                        }
                        Op::Ins(bs) => {
                            s.ops.push(2);
                            push_varint(&mut s.indel_lens, bs.len() as u64);
                            s.ins_bases.extend(bs.iter().map(|&b| base_code(b)));
                        }
                        Op::Del(n) => {
                            s.ops.push(3);
                            push_varint(&mut s.indel_lens, u64::from(*n));
                        }
                    }
                }
            }
            s
        })
        .collect();

    // Reads that landed on no reference code standalone — the graceful
    // degradation path, and the number to watch: at 300x the unsubsampled layout
    // stranded 40146 of them as 496 Mbase of literals.
    let mut orphans = Streams::default();
    for (r, b) in best.iter().enumerate() {
        if b.is_none() {
            orphans
                .literals
                .extend(read_at(r).iter().map(|&x| base_code(x)));
        }
    }
    per_contig.push(orphans);
    println!("scripts: {:.1}s", t.elapsed().as_secs_f64());

    let mut all = Streams::default();
    for s in per_contig {
        all.merge(s);
    }

    println!("--- streams (order-1 rANS) ---");
    // Encode the seven streams CONCURRENTLY — they are independent, and doing
    // them one after another left ~10s of a 67s run on a single core.
    let named: Vec<(&str, &[u8])> = vec![
        ("ops", &all.ops),
        ("match-runs", &all.runs),
        ("sub-bases", &all.subs),
        ("ins-bases", &all.ins_bases),
        ("indel-lens", &all.indel_lens),
        ("placements", &all.placements),
        ("flips", &all.flips),
        ("literals", &all.literals),
    ];
    let sizes: Vec<usize> = named.par_iter().map(|(_, b)| encoded_len(b)).collect();
    let mut total = 0usize;
    for ((label, buf), enc) in named.iter().zip(&sizes) {
        report(label, buf.len(), *enc, total_bases);
        total += *enc;
    }

    // Reference: 2-bit packed (refpack's 2bit+LZMA would do better).
    let ref_bytes = ref_bases.div_ceil(4);
    println!(
        "    reference    {ref_bases} bases -> {ref_bytes} B (2-bit, conservative) \
         ({:.4} bits/base)",
        (ref_bytes as f64 * 8.0) / total_bases as f64
    );

    let grand = total + ref_bytes;
    println!("  ---");
    println!(
        "  TOTAL: {grand} B = {:.4} bits/base   (fqxv-seq today 0.653 · CoLoRd 0.068 · \
         M1b oracle 0.040)",
        (grand as f64 * 8.0) / total_bases as f64
    );
}
