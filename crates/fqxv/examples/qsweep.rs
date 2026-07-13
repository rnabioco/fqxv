//! Quality-context sweep: measure the compressed size of the quality stream
//! under several context models, on real data, to guide the fqzcomp context.
//!
//! The shipping context (fqxv-fqzcomp) is `q1(6) | q2>>2(4) | delta(2) |
//! q3>>4(2) | pos-bucket(4)` — it treats the previous quality as an offset from
//! qmin (up to 6 bits) even when the alphabet is tiny, and coarsens q2/q3. For
//! binned data (few distinct values) a DENSE remap lets several full-resolution
//! previous qualities fit in the same budget, freeing bits for position — which
//! fqzcomp's presets weight heavily. This harness ranks such variants by bytes.
//!
//! Each config codes the same dense symbols, so sizes are directly comparable;
//! all use `SimpleModel::<64>` like the shipping codec, so per-context warmup is
//! identical and only the context structure differs.
//!
//! Usage: cargo run --release -p fqxv --example qsweep -- <fastq>

use std::fs::File;
use std::io::BufReader;

use fqxv_range::{Encoder, SimpleModel};

const QMAX: usize = 64;
/// Max previous qualities any config uses.
const KMAX: usize = 6;

/// A context configuration over dense (remapped) qualities.
struct Cfg {
    name: &'static str,
    /// Number of previous qualities in the context.
    k_prev: usize,
    /// Position cap (position feature = min(pos, pos_cap); pos_cap+1 buckets).
    pos_cap: usize,
    /// Delta (transition-count) cap (delta feature = min(delta, delta_cap)).
    delta_cap: usize,
}

impl Cfg {
    fn table_size(&self, nsym: usize) -> usize {
        nsym.pow(self.k_prev as u32) * (self.pos_cap + 1) * (self.delta_cap + 1)
    }
    #[inline]
    fn ctx(&self, prev: &[u8; KMAX], delta: u8, pos: usize, nsym: usize) -> usize {
        let mut idx = 0usize;
        for j in 0..self.k_prev {
            idx = idx * nsym + prev[j] as usize;
        }
        idx = idx * (self.pos_cap + 1) + pos.min(self.pos_cap);
        idx * (self.delta_cap + 1) + (delta as usize).min(self.delta_cap)
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: qsweep <fastq>");
    let mut reader =
        noodles_fastq::io::Reader::new(BufReader::new(File::open(&path).expect("open fastq")));
    let mut rec = noodles_fastq::Record::default();
    let mut lens: Vec<u32> = Vec::new();
    let mut quals: Vec<u8> = Vec::new();
    while reader.read_record(&mut rec).expect("read record") != 0 {
        let q = rec.quality_scores();
        lens.push(q.len() as u32);
        quals.extend_from_slice(q);
    }

    // Dense remap of the distinct quality values.
    let mut present = [false; 256];
    for &b in &quals {
        present[b as usize] = true;
    }
    let mut qmap = [0u8; 256];
    let mut nsym = 0usize;
    for (v, &p) in present.iter().enumerate() {
        if p {
            qmap[v] = nsym as u8;
            nsym += 1;
        }
    }
    let qmin = present.iter().position(|&p| p).unwrap_or(0) as u8;
    let qsize = present
        .iter()
        .rposition(|&p| p)
        .map_or(1, |m| m + 1 - qmin as usize);
    eprintln!(
        "{path}: {} reads, {} quals, {nsym} distinct values (qmin={qmin}, qsize={qsize})",
        lens.len(),
        quals.len()
    );

    // Baseline: replicate the shipping fqxv-fqzcomp context exactly.
    let baseline = encode_fqxv(&lens, &quals, qmin, qsize);
    report("fqxv-current (18b)", baseline, quals.len());

    let configs = [
        Cfg {
            name: "q1 dense + pos150 + d4",
            k_prev: 1,
            pos_cap: 149,
            delta_cap: 3,
        },
        Cfg {
            name: "q1q2 + pos150 + d4",
            k_prev: 2,
            pos_cap: 149,
            delta_cap: 3,
        },
        Cfg {
            name: "q1q2q3 + pos150 + d4",
            k_prev: 3,
            pos_cap: 149,
            delta_cap: 3,
        },
        Cfg {
            name: "q1q2q3q4 + pos150 + d4",
            k_prev: 4,
            pos_cap: 149,
            delta_cap: 3,
        },
        Cfg {
            name: "q1q2q3q4 + pos150 + d0",
            k_prev: 4,
            pos_cap: 149,
            delta_cap: 0,
        },
        Cfg {
            name: "q1q2q3q4q5 + pos150 + d4",
            k_prev: 5,
            pos_cap: 149,
            delta_cap: 3,
        },
        Cfg {
            name: "q1q2q3 + pos256 + d8",
            k_prev: 3,
            pos_cap: 255,
            delta_cap: 7,
        },
        Cfg {
            name: "q1q2q3q4 + pos48 + d4",
            k_prev: 4,
            pos_cap: 47,
            delta_cap: 3,
        },
        // Larger contexts — more training data on the full file may favor these.
        Cfg {
            name: "q1q2q3q4 + pos256 + d8",
            k_prev: 4,
            pos_cap: 255,
            delta_cap: 7,
        },
        Cfg {
            name: "q1q2q3 + pos256 + d16",
            k_prev: 3,
            pos_cap: 255,
            delta_cap: 15,
        },
        Cfg {
            name: "q1q2q3q4q5 + pos256 + d4",
            k_prev: 5,
            pos_cap: 255,
            delta_cap: 3,
        },
        Cfg {
            name: "q1q2q3q4q5q6 + pos150 + d4",
            k_prev: 6,
            pos_cap: 149,
            delta_cap: 3,
        },
    ];
    for c in &configs {
        let ts = c.table_size(nsym);
        if ts > (1 << 24) {
            println!("  {:<28} skipped (table {} too large)", c.name, ts);
            continue;
        }
        let bytes = encode_cfg(c, &lens, &quals, &qmap, nsym, ts);
        report(&format!("{} [{}ctx]", c.name, ts), bytes, quals.len());
    }
}

fn report(name: &str, bytes: usize, n_quals: usize) {
    println!(
        "  {:<40} {:>12} B  {:>6.1} MB  {:>6.3} bits/q",
        name,
        bytes,
        bytes as f64 / 1e6,
        8.0 * bytes as f64 / n_quals as f64,
    );
}

/// Exact replica of the shipping fqxv-fqzcomp context, for a fair baseline.
fn encode_fqxv(lens: &[u32], quals: &[u8], qmin: u8, _qsize: usize) -> usize {
    let pos_bucket = |pos: usize| -> usize {
        if pos < 16 {
            pos >> 1
        } else {
            (8 + (pos >> 5)).min(15)
        }
    };
    let ctx = |q1: u8, q2: u8, q3: u8, delta: u8, pos: usize| -> usize {
        (q1 as usize)
            | ((q2 as usize >> 2) << 6)
            | ((delta as usize) << 10)
            | ((q3 as usize >> 4) << 12)
            | (pos_bucket(pos) << 14)
    };
    let mut models = vec![SimpleModel::<QMAX>::new(); 1 << 18];
    let mut enc = Encoder::new();
    let mut rest = quals;
    for &l in lens {
        let (read, tail) = rest.split_at(l as usize);
        rest = tail;
        let (mut q1, mut q2, mut q3, mut delta) = (0u8, 0u8, 0u8, 0u8);
        for (pos, &b) in read.iter().enumerate() {
            let sym = b - qmin;
            let c = ctx(q1, q2, q3, delta, pos);
            models[c].encode(&mut enc, sym as usize);
            if pos > 0 && sym != q1 {
                delta = (delta + 1).min(3);
            }
            q3 = q2;
            q2 = q1;
            q1 = sym;
        }
    }
    enc.finish().len()
}

fn encode_cfg(
    c: &Cfg,
    lens: &[u32],
    quals: &[u8],
    qmap: &[u8; 256],
    nsym: usize,
    ts: usize,
) -> usize {
    let mut models = vec![SimpleModel::<QMAX>::new(); ts];
    let mut enc = Encoder::new();
    let mut rest = quals;
    for &l in lens {
        let (read, tail) = rest.split_at(l as usize);
        rest = tail;
        let mut prev = [0u8; KMAX];
        let mut delta = 0u8;
        for (pos, &b) in read.iter().enumerate() {
            let d = qmap[b as usize];
            let cx = c.ctx(&prev, delta, pos, nsym);
            models[cx].encode(&mut enc, d as usize);
            if pos > 0 && d != prev[0] {
                delta = delta.saturating_add(1);
            }
            for j in (1..KMAX).rev() {
                prev[j] = prev[j - 1];
            }
            prev[0] = d;
        }
    }
    enc.finish().len()
}
