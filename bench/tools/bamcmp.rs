//! bamcmp — order-independent alignment digests + quality distortion, for the
//! BAM round-trip proof (bam_identity.sh). Replaces a slow `samtools sort -n |
//! cut | sort | md5sum` per variant and a pure-Python per-base quality loop.
//!
//! Two subcommands, both single O(n) streaming passes with bounded memory:
//!
//!   digest    read headerless SAM (from `samtools view`) on stdin, emit three
//!             multiset digests. Each hashes a field range of every record and
//!             *sums* the hashes — addition is commutative, so the total is
//!             independent of record order and survives read reordering AND
//!             renaming, while preserving multiplicity (a true multiset digest):
//!               content = whole record        (QNAME + FLAG..QUAL + tags)
//!               body    = record minus QNAME  (FLAG..QUAL + tags)
//!               place   = FLAG..SEQ only       (no QNAME, no QUAL, no tags)
//!             `place` is invariant under reorder, rename, AND lossy quality
//!             binning (reads don't move); `body` adds QUAL so it changes only
//!             when quality changes; `content` adds QNAME+order sensitivity for
//!             the strict lossless check. Same digest ⇔ same record multiset.
//!
//!   qualdelta A.fastq B.fastq   per-base Phred distortion between two FASTQs in
//!             the SAME order (original vs a lossy round-trip). Prints:
//!               <nbases> <changed> <mean_abs> <rmse> <max>
//!
//! Not fqxv's own hash — a standalone harness tool. Build (same as fqdigest):
//!   rustc -O --edition 2021 bench/bamcmp.rs -o <path>/bamcmp

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

// wyhash-style multiply-fold primitives, identical to fqdigest.rs: fast, good
// avalanche, far cheaper than std SipHash for many small writes.
const P0: u64 = 0xa076_1d64_78bd_642f;
const P1: u64 = 0xe703_7ed1_a0b4_28db;
const P2: u64 = 0x8ebc_6af0_9c88_c6e3;

#[inline(always)]
fn mum(a: u64, b: u64) -> u64 {
    let r = (a as u128).wrapping_mul(b as u128);
    (r as u64) ^ ((r >> 64) as u64)
}

#[inline(always)]
fn fold(mut h: u64, data: &[u8]) -> u64 {
    let mut chunks = data.chunks_exact(8);
    for c in &mut chunks {
        let v = u64::from_le_bytes(c.try_into().unwrap());
        h = mum(h ^ P1, v ^ P0);
    }
    let rem = chunks.remainder();
    if !rem.is_empty() {
        let mut last = [0u8; 8];
        last[..rem.len()].copy_from_slice(rem);
        h = mum(h ^ P1, u64::from_le_bytes(last) ^ P0);
    }
    mum(h, (data.len() as u64) ^ P2)
}

/// One 64-bit lane of a multiset sum: fold `data` into a domain-seeded state.
#[inline(always)]
fn lane(domain: u64, data: &[u8]) -> u64 {
    mum(fold(domain ^ P0, data) ^ P0, P1)
}

/// A 128-bit commutative multiset accumulator (two independent lanes).
#[derive(Default, Clone, Copy)]
struct Acc {
    lo: u64,
    hi: u64,
    n: u64,
}
impl Acc {
    #[inline(always)]
    fn add(&mut self, data: &[u8]) {
        self.lo = self.lo.wrapping_add(lane(1, data));
        self.hi = self.hi.wrapping_add(lane(2, data));
        self.n += 1;
    }
}

fn trim_len(l: &[u8]) -> usize {
    let mut end = l.len();
    if end > 0 && l[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && l[end - 1] == b'\r' {
        end -= 1;
    }
    end
}

/// Byte offset of the nth (0-based) tab in `line`, or None if there are fewer.
#[inline]
fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
    line.iter()
        .enumerate()
        .filter(|&(_, &b)| b == b'\t')
        .nth(n)
        .map(|(i, _)| i)
}

fn digest<R: Read>(r: R) -> io::Result<(Acc, Acc, Acc)> {
    let mut br = BufReader::with_capacity(1 << 20, r);
    let (mut content, mut body, mut place) = (Acc::default(), Acc::default(), Acc::default());
    let mut line = Vec::new();
    loop {
        line.clear();
        if br.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let rec = &line[..trim_len(&line)];
        if rec.is_empty() {
            continue;
        }
        content.add(rec);
        // body = everything after the first tab (drop QNAME).
        if let Some(t0) = nth_tab(rec, 0) {
            body.add(&rec[t0 + 1..]);
            // place = FLAG..SEQ = [after tab0 .. tab9] (tab before QUAL, field 10).
            if let Some(t9) = nth_tab(rec, 9) {
                place.add(&rec[t0 + 1..t9]);
            } else {
                place.add(&rec[t0 + 1..]); // <11 fields: fall back to the body span
            }
        } else {
            body.add(rec);
            place.add(rec);
        }
    }
    Ok((content, body, place))
}

fn qualdelta(a: &str, b: &str) -> io::Result<(u64, u64, f64, f64, u64)> {
    let mut ra = BufReader::with_capacity(1 << 20, File::open(a)?);
    let mut rb = BufReader::with_capacity(1 << 20, File::open(b)?);
    let (mut la, mut lb) = (Vec::new(), Vec::new());
    let (mut n, mut changed, mut sabs, mut ssq, mut mx) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let mut i = 0usize;
    loop {
        la.clear();
        lb.clear();
        let ga = ra.read_until(b'\n', &mut la)?;
        let gb = rb.read_until(b'\n', &mut lb)?;
        if ga == 0 || gb == 0 {
            break;
        }
        if i % 4 == 3 {
            let qa = &la[..trim_len(&la)];
            let qb = &lb[..trim_len(&lb)];
            for (&ca, &cb) in qa.iter().zip(qb.iter()) {
                let d = (ca as i32 - cb as i32).unsigned_abs() as u64;
                n += 1;
                if d != 0 {
                    changed += 1;
                    sabs += d;
                    ssq += d * d;
                    if d > mx {
                        mx = d;
                    }
                }
            }
        }
        i += 1;
    }
    let (mean_abs, rmse) = if n == 0 {
        (0.0, 0.0)
    } else {
        (sabs as f64 / n as f64, (ssq as f64 / n as f64).sqrt())
    };
    Ok((n, changed, mean_abs, rmse, mx))
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out = io::stdout().lock();
    match args.first().map(String::as_str) {
        Some("digest") => {
            let (c, b, p) = digest(io::stdin().lock())?;
            writeln!(out, "content {:016x}{:016x} {}", c.lo, c.hi, c.n)?;
            writeln!(out, "body {:016x}{:016x} {}", b.lo, b.hi, b.n)?;
            writeln!(out, "place {:016x}{:016x} {}", p.lo, p.hi, p.n)?;
        }
        Some("qualdelta") if args.len() == 3 => {
            let (n, changed, mean_abs, rmse, mx) = qualdelta(&args[1], &args[2])?;
            writeln!(out, "{} {} {:.4} {:.4} {}", n, changed, mean_abs, rmse, mx)?;
        }
        _ => {
            eprintln!("usage:\n  samtools view IN.bam | bamcmp digest\n  bamcmp qualdelta ORIG.fastq BINNED.fastq");
            std::process::exit(2);
        }
    }
    Ok(())
}
