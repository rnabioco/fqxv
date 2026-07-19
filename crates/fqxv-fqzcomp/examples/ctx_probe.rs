//! Conditional-entropy probe for long-read quality context designs.
//!
//! Prints `H(quality | context)` in bits/qual for several candidate context
//! packings over the SAME positions of a FASTQ, so a lower number means the
//! design's features carry more signal. Use it to *rank* candidate features
//! before spending real coder bits on them.
//!
//! ```text
//! cargo run --release -p fqxv-fqzcomp --example ctx_probe -- reads.fastq [Nreads]
//! ```
//!
//! IMPORTANT — empirical entropy is a ranking tool, not a verdict. It assumes
//! each context's symbol distribution is known perfectly. The real codec uses a
//! per-context *adaptive* model ([`fqxv_range::SimpleModel`]) that must learn each
//! distribution online, paying a cold-start cost per context. Adding context bits
//! multiplies that cost, and past a point it swamps the entropy gain — so a design
//! that wins here can still lose on real bytes. A HiFi sweep found exactly this:
//! the richer contexts below have lower entropy yet produced *larger* archives
//! than the shipped 16-bit context. Always confirm a winner by measuring the real
//! `quality` stream (`fqxv info --json`) at full block size before shipping it.

use std::fs::File;
use std::io::{BufRead, BufReader};

/// Phred range 0..=93.
const QN: usize = 94;

/// 2-bit base code (A/C/G/T -> 0/1/2/3; anything else -> 0), matching the codec.
fn bc(b: u8) -> u32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// One candidate context design: a name, the number of packed bits (so the table
/// is sized correctly), and a running per-context symbol histogram.
struct Design {
    name: &'static str,
    counts: Vec<u32>,
}

impl Design {
    fn new(name: &'static str, bits: u32) -> Self {
        Design {
            name,
            counts: vec![0u32; (1usize << bits) * QN],
        }
    }

    fn bump(&mut self, ctx: u32, q: u32) {
        self.counts[(ctx as usize) * QN + q as usize] += 1;
    }

    /// `(H in bits/qual, populated context count)`.
    fn entropy(&self) -> (f64, usize) {
        let mut acc = 0.0f64;
        let mut tot = 0u64;
        let mut nctx = 0usize;
        for row in self.counts.chunks_exact(QN) {
            let n: u64 = row.iter().map(|&c| u64::from(c)).sum();
            if n == 0 {
                continue;
            }
            nctx += 1;
            tot += n;
            let nf = n as f64;
            for &c in row {
                if c != 0 {
                    let cf = f64::from(c);
                    acc += -cf * (cf / nf).log2();
                }
            }
        }
        (acc / tot as f64, nctx)
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: ctx_probe <fastq> [Nreads]");
    let nmax: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);

    // Candidates. `shipped` is the current long-read context; the rest add features
    // ranked by an earlier sweep (q3 = third-previous quality, next2 = next-next
    // base, prevbase = base at i-1). Bit widths must bound the packed context.
    let mut designs = vec![
        Design::new("shipped  q1>>1,q2>>3,base,next,hp           (16b)", 16),
        Design::new("+prevbase                                   (18b)", 18),
        Design::new("+prevbase +next2                            (20b)", 20),
        Design::new("+next2 +q3>>4                               (20b)", 20),
        Design::new("+prevbase +next2 +q3>>4                     (22b)", 22),
    ];

    let mut qseen = [false; 256];
    let mut totbases = 0u64;
    let mut nreads = 0usize;
    let mut lines = BufReader::new(File::open(path).expect("open fastq")).lines();
    while nreads < nmax {
        let Some(_h) = lines.next() else { break };
        let seq = lines.next().expect("seq line").unwrap();
        let _plus = lines.next().expect("+ line").unwrap();
        let qual = lines.next().expect("qual line").unwrap();
        nreads += 1;
        let sb = seq.as_bytes();
        let qb = qual.as_bytes();
        totbases += qb.len() as u64;
        for &b in qb {
            qseen[b as usize] = true;
        }
        let (mut q1, mut q2, mut q3) = (0u32, 0u32, 0u32);
        let mut run = 0u32;
        let mut prev = 255u8;
        let mut prevbase = 0u32;
        for i in 0..qb.len() {
            let q = u32::from(qb[i].wrapping_sub(33)).min(93);
            let base = if i < sb.len() { bc(sb[i]) } else { 0 };
            let next = if i + 1 < sb.len() { bc(sb[i + 1]) } else { 0 };
            let next2 = if i + 2 < sb.len() { bc(sb[i + 2]) } else { 0 };
            let cur = if i < sb.len() { sb[i] } else { 0 };
            run = if i > 0 && cur == prev { run + 1 } else { 1 };
            let hp = run.min(7);
            let f1 = q1 >> 1; // fine q1, 6 bits
            let q2c = (q2 >> 3) & 0x7; // coarse q2, 3 bits
            let q3c = (q3 >> 4) & 0x3; // coarse q3, 2 bits

            let core = f1 | (q2c << 6) | (base << 9) | (next << 11) | (hp << 13);
            designs[0].bump(core, q);
            designs[1].bump(core | (prevbase << 16), q);
            designs[2].bump(core | (prevbase << 16) | (next2 << 18), q);
            designs[3].bump(core | (next2 << 16) | (q3c << 18), q);
            designs[4].bump(core | (prevbase << 16) | (next2 << 18) | (q3c << 20), q);

            q3 = q2;
            q2 = q1;
            q1 = q;
            prev = cur;
            prevbase = base;
        }
    }

    let k = qseen.iter().filter(|&&x| x).count();
    let mean_len = if nreads > 0 { totbases / nreads as u64 } else { 0 };
    println!("reads={nreads}  distinct_qual_k={k}  mean_read_len={mean_len}");
    let base_h = designs[0].entropy().0;
    for d in &designs {
        let (h, nc) = d.entropy();
        println!("H = {h:.4} bits/qual  d-vs-shipped {:+.4}  [ctxs={nc:>9}]  {}", h - base_h, d.name);
    }
}
