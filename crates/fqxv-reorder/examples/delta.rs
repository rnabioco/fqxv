//! Measure a reorder-paired differential sequence coder against fqxv-seq.
//!
//! Usage: `cargo run --release -p fqxv-reorder --example delta -- <seq-file> <read-len> [k]`
//!
//! After minimizer clustering, each read is coded relative to the previous read
//! in the reordered stream: MATCH (identical), DELTA (same length, few
//! mismatches), or LITERAL. The op / mismatch / literal streams are then
//! entropy-coded with fqxv-rans. This is the SPRING-style mechanism — cheap
//! per-read coding that leans on clustering to expose duplicate redundancy.

use fqxv_rans::{encode as rans, Order};
use std::time::Instant;

const MATCH: u8 = 0;
const DELTA: u8 = 1;
const LITERAL: u8 = 2;

fn base2(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0, // N and friends fold to A here (measurement only)
    }
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            break;
        }
        out.push(byte | 0x80);
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: delta <file> <read-len> [k]");
    let read_len: usize = args.next().expect("read-len").parse().unwrap();
    let k: usize = args.next().unwrap_or_else(|| "15".into()).parse().unwrap();

    let raw = std::fs::read(&path).expect("read");
    let n = raw.len() / read_len;
    let seq = &raw[..n * read_len];
    let lens = vec![read_len as u32; n];

    let baseline = fqxv_seq::encode(&lens, seq, 11).unwrap().len();

    let t0 = Instant::now();
    let plan = fqxv_reorder::plan(&lens, seq, k);
    // Materialize reordered, canonically-oriented reads.
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n);
    for &oi in &plan.order {
        let r = &seq[oi as usize * read_len..(oi as usize + 1) * read_len];
        reads.push(if plan.flip[oi as usize] {
            fqxv_reorder::revcomp(r)
        } else {
            r.to_vec()
        });
    }

    // Differential coding vs the previous read. Literal reads are held aside and
    // coded with the fqxv-seq context model (not raw 2-bit), so dedup and the
    // context model stack.
    let (mut ops, mut nmis, mut pos, mut subs) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut lit_seq, mut lit_lens): (Vec<u8>, Vec<u32>) = (Vec::new(), Vec::new());
    let thresh = read_len / 4;
    let (mut n_match, mut n_delta, mut n_lit) = (0u64, 0u64, 0u64);
    for i in 0..reads.len() {
        let cur = &reads[i];
        let prev: &[u8] = if i > 0 { &reads[i - 1] } else { &[] };
        if prev == cur.as_slice() {
            ops.push(MATCH);
            n_match += 1;
            continue;
        }
        if prev.len() == cur.len() {
            let mism: Vec<usize> = (0..cur.len()).filter(|&j| prev[j] != cur[j]).collect();
            if mism.len() <= thresh {
                ops.push(DELTA);
                put_varint(&mut nmis, mism.len() as u64);
                let mut last = 0usize;
                for &m in &mism {
                    put_varint(&mut pos, (m - last) as u64);
                    last = m;
                    subs.push(base2(cur[m]));
                }
                n_delta += 1;
                continue;
            }
        }
        ops.push(LITERAL);
        lit_seq.extend_from_slice(cur);
        lit_lens.push(cur.len() as u32);
        n_lit += 1;
    }
    let plan_s = t0.elapsed().as_secs_f64();

    let c = |v: &[u8], o| rans(v, o).map(|x| x.len()).unwrap_or(v.len());
    let lits = fqxv_seq::encode(&lit_lens, &lit_seq, 11)
        .map(|v| v.len())
        .unwrap_or(lit_seq.len());
    let total = c(&ops, Order::One)
        + c(&nmis, Order::Zero)
        + c(&pos, Order::Zero)
        + c(&subs, Order::One)
        + lits;

    let bpb = |bytes: usize| (bytes * 8) as f64 / (n * read_len) as f64;
    println!("reads {n}  read_len {read_len}  k {k}  (plan {plan_s:.1}s)");
    println!("  reads: {n_match} match / {n_delta} delta / {n_lit} literal");
    println!(
        "  fqxv-seq order-11 (baseline)      {baseline:>10} bytes  {:.3} bits/base",
        bpb(baseline)
    );
    println!(
        "  reorder + delta + ctx-literals    {total:>10} bytes  {:.3} bits/base",
        bpb(total)
    );
    println!(
        "    ops {} nmis {} pos {} subs {} lits(ctx) {}",
        c(&ops, Order::One),
        c(&nmis, Order::Zero),
        c(&pos, Order::Zero),
        c(&subs, Order::One),
        lits,
    );
    println!(
        "  gain vs baseline: {:.1}% (order not preserved)",
        100.0 * (baseline as f64 - total as f64) / baseline as f64
    );
}
