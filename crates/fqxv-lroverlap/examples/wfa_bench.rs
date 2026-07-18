//! WFA vs banded DP: is score-proportional alignment faster for this codec?
//!
//! The codec aligns each long read against a consensus. The regimes are HiFi
//! (~13 kb at ~0.5% error, substitution-dominant) and ONT (~13 kb at ~10% error,
//! indel-heavy). `align_banded` is `O(n · band)` regardless of similarity;
//! `wfa_align` is score-proportional. This measures both at matched settings
//! across a divergence sweep, reports ns/alignment and WFA's `O(s²)` traceback
//! footprint, and prints the speed crossover.
//!
//! Both aligners are given the SAME band: the true diagonal drift of the applied
//! edits plus a margin, which is the minimal band that still lets the DP find the
//! optimum — so the two find the same edit distance and the comparison is
//! apples-to-apples (same answer, different work).
//!
//! ```text
//! cargo run --release --example wfa_bench
//! ```

use std::time::{Duration, Instant};

use fqxv_lroverlap::{align_banded, apply, wfa_align, wfa_cells};

/// Deterministic PRNG (SplitMix64) — no `rand`, fully seeded and reproducible.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

const BASES: &[u8; 4] = b"ACGT";

/// A synthetic reference and a mutated copy, plus the true diagonal-drift band
/// the optimal path needs (max |insertions − deletions| seen, + margin).
struct Pair {
    refr: Vec<u8>,
    query: Vec<u8>,
    band: usize,
}

/// Generate `refr` of length `len`, then walk it emitting `query`: with
/// probability `p` apply an edit, split into substitution / insertion / deletion
/// by the given weights (indel-heavy for ONT, substitution-heavy for HiFi).
fn make_pair(seed: u64, len: usize, p: f64, w_sub: u32, w_ins: u32, w_del: u32) -> Pair {
    let mut rng = Rng(seed);
    let refr: Vec<u8> = (0..len).map(|_| BASES[rng.below(4) as usize]).collect();

    let scale = 1_000_000u64;
    let thresh = (p * scale as f64) as u64;
    let w_tot = w_sub + w_ins + w_del;

    let mut query = Vec::with_capacity(len);
    let mut bal: i64 = 0; // insertions − deletions so far
    let mut max_drift: i64 = 0;
    let mut i = 0usize;
    while i < refr.len() {
        if rng.below(scale) >= thresh {
            query.push(refr[i]);
            i += 1;
            continue;
        }
        let r = rng.below(u64::from(w_tot)) as u32;
        if r < w_sub {
            // Substitution: a different base.
            let mut b = BASES[rng.below(4) as usize];
            while b == refr[i] {
                b = BASES[rng.below(4) as usize];
            }
            query.push(b);
            i += 1;
        } else if r < w_sub + w_ins {
            // Insertion: emit an extra base, consume none of `refr`.
            query.push(BASES[rng.below(4) as usize]);
            bal += 1;
        } else {
            // Deletion: skip a reference base.
            i += 1;
            bal -= 1;
        }
        max_drift = max_drift.max(bal.abs());
    }
    let band = max_drift as usize + 8;
    Pair { refr, query, band }
}

/// Run `f` in a loop until at least `budget` has elapsed, returning ns/call.
fn bench<F: FnMut()>(budget: Duration, mut f: F) -> f64 {
    // Warm up and estimate.
    let t0 = Instant::now();
    f();
    let est = t0.elapsed().max(Duration::from_nanos(1));
    let reps = ((budget.as_nanos() / est.as_nanos()).max(1)) as u64;
    let t = Instant::now();
    for _ in 0..reps {
        f();
    }
    t.elapsed().as_nanos() as f64 / reps as f64
}

fn main() {
    // (label, error%, length, sub/ins/del weights)
    let regimes: &[(&str, f64, usize, u32, u32, u32)] = &[
        ("HiFi   0.5%", 0.005, 13_000, 90, 5, 5),
        ("mid    2%  ", 0.02, 13_000, 70, 15, 15),
        ("mid    5%  ", 0.05, 13_000, 60, 20, 20),
        ("ONT    10% ", 0.10, 13_000, 40, 30, 30),
        ("ONT    14% ", 0.14, 13_000, 40, 30, 30),
    ];
    // Distinct seeded pairs per regime, averaged, so a single lucky layout does
    // not decide the number.
    let pairs_per_regime = 4;
    let budget = Duration::from_millis(150);

    println!(
        "{:<12} {:>7} {:>6} {:>7} {:>7} {:>13} {:>13} {:>8} {:>10}",
        "regime", "dist", "band", "wfaMB", "dpMB", "wfa ns", "dp ns", "speedup", "dp/wfa mem"
    );
    println!("{}", "-".repeat(96));

    for &(label, p, len, ws, wi, wd) in regimes {
        let pairs: Vec<Pair> = (0..pairs_per_regime)
            .map(|s| {
                make_pair(
                    0xC0FFEE ^ (s as u64) << 20 ^ (len as u64),
                    len,
                    p,
                    ws,
                    wi,
                    wd,
                )
            })
            .collect();

        // Sanity: both aligners round-trip and agree on the optimal distance.
        let mut dist_sum = 0u64;
        let mut band_sum = 0usize;
        let mut wfa_cells_sum = 0usize;
        let mut dp_cells_sum = 0usize;
        for pr in &pairs {
            let big = (pr.refr.len() + pr.query.len()) as u32;
            let w = wfa_align(&pr.refr, &pr.query, big);
            let d = align_banded(&pr.refr, &pr.query, pr.band);
            assert_eq!(apply(&pr.refr, &w.ops), pr.query, "WFA must round-trip");
            assert_eq!(apply(&pr.refr, &d.ops), pr.query, "DP must round-trip");
            assert_eq!(
                w.dist, d.dist,
                "WFA and exact-band DP must find the same edit distance"
            );
            dist_sum += u64::from(w.dist);
            band_sum += pr.band;
            wfa_cells_sum += wfa_cells(&pr.refr, &pr.query, big);
            // DP footprint: (n+1)·stride cells; traceback packs 2 bits/cell, so
            // ~stride/4 bytes/row, plus two i32 score rows. The traceback matrix
            // dominates, so report it as the comparable stored bytes.
            let stride = (2 * pr.band + 1).min(pr.query.len() + 1);
            dp_cells_sum += (pr.refr.len() + 1) * stride;
        }
        let n = pairs.len();
        let dist = dist_sum as f64 / n as f64;
        let band = band_sum / n;
        // WFA stores i32 offsets (4 B); DP traceback packs 2 bits/cell (0.25 B).
        let wfa_mb = (wfa_cells_sum as f64 * 4.0) / (n as f64 * 1e6);
        let dp_mb = (dp_cells_sum as f64 * 0.25) / (n as f64 * 1e6);

        let wfa_ns = bench(budget, || {
            for pr in &pairs {
                let big = (pr.refr.len() + pr.query.len()) as u32;
                std::hint::black_box(wfa_align(&pr.refr, &pr.query, big));
            }
        }) / n as f64;
        let dp_ns = bench(budget, || {
            for pr in &pairs {
                std::hint::black_box(align_banded(&pr.refr, &pr.query, pr.band));
            }
        }) / n as f64;

        println!(
            "{:<12} {:>7.0} {:>6} {:>7.2} {:>7.3} {:>13.0} {:>13.0} {:>7.2}x {:>9.1}x",
            label,
            dist,
            band,
            wfa_mb,
            dp_mb,
            wfa_ns,
            dp_ns,
            dp_ns / wfa_ns,
            dp_mb / wfa_mb,
        );
    }

    println!();
    println!("speedup > 1 means WFA is faster; dp/wfa mem > 1 means WFA uses more memory.");
}
