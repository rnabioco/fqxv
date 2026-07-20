//! Diagnostic: compare raw-value vs splitmix64-hashed minimizer selection for
//! the short-read clustering key in `plan::min_canonical`.
//!
//! The library selects the read's representative k-mer by comparing raw 2-bit
//! packed values, so the all-A k-mer (which packs to 0) wins any read that
//! contains one and selection is biased toward A-rich k-mers generally.
//! `fqxv-lroverlap` hashes for exactly this reason. This measures whether the
//! bias actually degrades clustering on real short-read data.
//!
//! Usage: `keystats <in.fastq> [max_reads]`

use std::collections::HashMap;
use std::env;
use std::fs;

const K: usize = 15;

#[inline]
fn code_fold(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 255,
    }
}

/// splitmix64, as used by `fqxv-lroverlap::minimizer::hash64`.
#[inline]
fn hash64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Minimum canonical k-mer of `read`. When `hash` is set, selection compares
/// splitmix64 of the canonical value; the returned key is the canonical k-mer
/// either way so cluster identity is comparable between the two schemes.
/// Mirrors `plan::min_canonical`.
fn min_canonical(read: &[u8], k: usize, hash: bool) -> u64 {
    if read.len() < k || k == 0 || k > 32 {
        return u64::MAX;
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k as u64 - 1);
    let (mut fwd, mut rc, mut valid) = (0u64, 0u64, 0usize);
    let (mut best_score, mut best_key) = (u64::MAX, u64::MAX);
    for &b in read {
        let c = code_fold(b);
        if c == 255 {
            fwd = 0;
            rc = 0;
            valid = 0;
            continue;
        }
        let c = u64::from(c);
        fwd = ((fwd << 2) | c) & mask;
        rc = ((rc >> 2) | ((3 - c) << shift)) & mask;
        valid += 1;
        if valid >= k {
            // Strand choice stays a raw comparison (revcomp-symmetric); only the
            // *winner* selection changes.
            let canon = if fwd <= rc { fwd } else { rc };
            let score = if hash { hash64(canon) } else { canon };
            if score < best_score {
                best_score = score;
                best_key = canon;
            }
        }
    }
    best_key
}

/// Base composition of a packed k-mer: returns the count of the most frequent
/// base, used as a crude low-complexity signal.
fn max_base_run(mut key: u64, k: usize) -> u32 {
    let mut counts = [0u32; 4];
    for _ in 0..k {
        counts[(key & 3) as usize] += 1;
        key >>= 2;
    }
    *counts.iter().max().unwrap()
}

struct Stats {
    label: &'static str,
    reads: usize,
    distinct: usize,
    reads_in_singletons: usize,
    largest: usize,
    top: Vec<usize>,
    reads_in_top_1k: usize,
}

fn summarize(label: &'static str, keys: &[u64]) -> Stats {
    let mut counts: HashMap<u64, usize> = HashMap::with_capacity(keys.len());
    for &key in keys {
        *counts.entry(key).or_insert(0) += 1;
    }
    let mut sizes: Vec<usize> = counts.values().copied().collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    let reads_in_singletons = sizes.iter().filter(|&&s| s == 1).count();
    let reads_in_top_1k: usize = sizes.iter().take(1000).sum();
    Stats {
        label,
        reads: keys.len(),
        distinct: sizes.len(),
        reads_in_singletons,
        largest: sizes.first().copied().unwrap_or(0),
        top: sizes.iter().take(10).copied().collect(),
        reads_in_top_1k,
    }
}

fn report(s: &Stats) {
    let pct = |n: usize| 100.0 * n as f64 / s.reads as f64;
    println!("  [{}]", s.label);
    println!("    distinct keys        {:>12}", s.distinct);
    println!(
        "    reads in singletons  {:>12}  ({:.2}%)   <- cannot cluster at all",
        s.reads_in_singletons,
        pct(s.reads_in_singletons)
    );
    println!(
        "    reads in clusters>=2 {:>12}  ({:.2}%)",
        s.reads - s.reads_in_singletons,
        pct(s.reads - s.reads_in_singletons)
    );
    println!(
        "    largest cluster      {:>12}  ({:.2}%)",
        s.largest,
        pct(s.largest)
    );
    println!(
        "    reads in top-1k clus {:>12}  ({:.2}%)   <- concentration",
        s.reads_in_top_1k,
        pct(s.reads_in_top_1k)
    );
    println!("    top-10 cluster sizes {:?}", s.top);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let path = args.get(1).expect("usage: keystats <in.fastq> [max_reads]");
    let max: usize = args
        .get(2)
        .map_or(usize::MAX, |s| s.parse().expect("max_reads"));

    eprintln!("reading {path} ...");
    let data = fs::read(path).expect("read fastq");

    // Every 4th line starting at line 1 is a sequence.
    let mut raw = Vec::new();
    let mut hashed = Vec::new();
    let mut nolow = 0usize; // winning k-mer is low-complexity (>=80% one base)
    let mut polya = 0usize; // winning k-mer is a pure homopolymer
    let mut nolow_h = 0usize;
    let mut polya_h = 0usize;
    let mut invalid = 0usize;

    for (i, line) in data.split(|&b| b == b'\n').enumerate() {
        if i % 4 != 1 {
            continue;
        }
        if raw.len() >= max {
            break;
        }
        if line.is_empty() {
            continue;
        }
        let kr = min_canonical(line, K, false);
        let kh = min_canonical(line, K, true);
        if kr == u64::MAX {
            invalid += 1;
            continue;
        }
        let thresh = (K as u32 * 4).div_ceil(5); // 80%
        if max_base_run(kr, K) >= thresh {
            nolow += 1;
        }
        if max_base_run(kr, K) == K as u32 {
            polya += 1;
        }
        if max_base_run(kh, K) >= thresh {
            nolow_h += 1;
        }
        if max_base_run(kh, K) == K as u32 {
            polya_h += 1;
        }
        raw.push(kr);
        hashed.push(kh);
    }

    let n = raw.len();
    println!("\n=== {path} ===");
    println!("reads analyzed {n}  (k={K}, no-valid-kmer skipped: {invalid})");
    let pct = |x: usize| 100.0 * x as f64 / n as f64;
    println!("\n  selected-k-mer composition (the bias itself):");
    println!(
        "    raw:    low-complexity (>=80% one base) {:>10} ({:.2}%)   pure homopolymer {:>8} ({:.3}%)",
        nolow,
        pct(nolow),
        polya,
        pct(polya)
    );
    println!(
        "    hashed: low-complexity (>=80% one base) {:>10} ({:.2}%)   pure homopolymer {:>8} ({:.3}%)",
        nolow_h,
        pct(nolow_h),
        polya_h,
        pct(polya_h)
    );
    println!("\n  clustering quality:");
    report(&summarize("raw (current)", &raw));
    report(&summarize("hashed", &hashed));
}
