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
    align_banded, consensus, find_overlaps, layout, place_against, ChainOpts, ConsensusOpts, Index,
    LayoutOpts, Op, Repeat, Sketch,
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

    let t = Instant::now();
    let idx = Index::build(&lens, &seq, sketch, Repeat::default()).expect("index");
    println!(
        "index: {} occs · {:.1}s",
        idx.len(),
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let ovs: Vec<Vec<_>> = (0..lens.len())
        .into_par_iter()
        .map(|r| find_overlaps(&idx, r as u32, read_at(r), ChainOpts::default()))
        .collect();
    println!("overlaps: {:.1}s", t.elapsed().as_secs_f64());

    let t = Instant::now();
    let contigs = layout(&lens, &ovs, LayoutOpts::default());
    let placed: usize = contigs
        .iter()
        .filter(|c| c.reads.len() > 1)
        .map(|c| c.reads.len())
        .sum();
    let singletons = contigs.iter().filter(|c| c.reads.len() == 1).count();
    println!(
        "layout: {} contigs · {} reads on multi-read contigs · {} singletons · {:.1}s",
        contigs.len(),
        placed,
        singletons,
        t.elapsed().as_secs_f64()
    );

    // Per contig: consensus, then code every read against it.
    let t = Instant::now();
    let per_contig: Vec<(Streams, usize)> = contigs
        .par_iter()
        .map(|c| {
            let mut s = Streams::default();
            // Orient reads as the layout placed them.
            let oriented: Vec<Vec<u8>> = c
                .reads
                .iter()
                .map(|p| {
                    let r = read_at(p.read as usize);
                    if p.flip {
                        revcomp(r)
                    } else {
                        r.to_vec()
                    }
                })
                .collect();
            // A singleton has no reference: its bases are literals.
            if c.reads.len() == 1 {
                s.literals.extend(oriented[0].iter().map(|&b| base_code(b)));
                return (s, 0);
            }
            // `consensus` indexes by read id, so build a full-length table.
            let mut by_id: Vec<Vec<u8>> = vec![Vec::new(); lens.len()];
            for (i, p) in c.reads.iter().enumerate() {
                by_id[p.read as usize] = oriented[i].clone();
            }
            let cons = consensus(c, &by_id, ConsensusOpts::default());
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
                for o in &oriented {
                    s.literals.extend(o.iter().map(|&b| base_code(b)));
                }
                return (s, 0);
            }

            // M3 REFINEMENT. Do NOT reuse the layout's offsets here. They are
            // composed hop-by-hop from chain point estimates, so their indel
            // error accumulates across the contig with nothing re-anchoring it —
            // measured at 0.036 subs/base emitted where 0.0025 exist. Re-place
            // every read against the finished consensus instead: one hop from a
            // fixed frame, so error cannot compound.
            let refs: Vec<&[u8]> = oriented.iter().map(|v| v.as_slice()).collect();
            let anchored = place_against(&cons.seq, &refs, sketch, ChainOpts::default());

            // Align in PARALLEL, then emit serially. The alignment is ~10 ms per
            // read (a 12.9 kb read at band 96 is ~2.5M cells) and there are
            // thousands per contig; the emission is microseconds. Doing both in
            // one loop made the whole thing serial, because delta-coded
            // placements depend on the previous read's start — so the cheap
            // order-dependent step was forcing the expensive order-independent
            // one to run on a single core.
            let aligned: Vec<Option<(usize, Vec<Op>)>> = (0..c.reads.len())
                .into_par_iter()
                .map(|i| {
                    let read = &oriented[i];
                    let a = anchored[i]?;
                    // `oriented` is already in the layout's frame, so a further
                    // flip here means the layout and the consensus disagree about
                    // orientation. Treat as unplaceable rather than guess.
                    if a.flip {
                        return None;
                    }
                    let start = a.offset as usize;
                    let end = (start + read.len() + 64).min(cons.seq.len());
                    if start >= end {
                        return None;
                    }
                    Some((start, align_banded(&cons.seq[start..end], read, 96).ops))
                })
                .collect();

            let mut prev_start: i64 = 0;
            for (i, slot) in aligned.into_iter().enumerate() {
                let Some((start, ops)) = slot else {
                    // Unplaceable -> code standalone.
                    s.literals.extend(oriented[i].iter().map(|&b| base_code(b)));
                    continue;
                };
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
            (s, cons.seq.len())
        })
        .collect();
    println!("consensus + scripts: {:.1}s", t.elapsed().as_secs_f64());

    let mut all = Streams::default();
    let mut ref_bases = 0usize;
    for (s, rb) in per_contig {
        all.merge(s);
        ref_bases += rb;
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
