//! A/B of the **lroverlap** shared-reference frame coder.
//!
//! `Reference::encode` currently stores the consensus bases with order-0/order-1
//! rANS over the *raw* 4-symbol bases. `fqxv-reorder` already established that
//! this is the wrong representation for a consensus reference: 2-bit-packing
//! first and running a byte-domain LZ over the packed bytes beats it, because
//! the packing is a hard 2 bits/base floor and the LZ then reaches the
//! long-range near-duplicate-contig repeats a context model cannot see (see
//! `fqxv-reorder::refpack`, REF_METHOD_PACK, −5.9% there).
//!
//! This measures the same swap for lroverlap (issue #197).
//!
//! **Result on PacBio Revio WGS (SRR36938642, 40k reads): the lever is not
//! there.** 2-bit+LZMA saves only **3.07%** of the frame (5,670,352 →
//! 5,496,143 B), against reorder's −5.9%, and raw LZMA is 11.5% *worse* than the
//! rANS coder shipping today:
//!
//! ```text
//! rANS raw (now)        5,670,352 B   1.921 b/base    0.00%
//! LZMA raw              6,323,130 B   2.143 b/base  +11.51%
//! 2bit+LZMA (refpack)   5,496,143 B   1.862 b/base   -3.07%
//! ```
//!
//! The reason is visible in the b/base column: the consensus is already at
//! **1.92 bits/base**, essentially the 2-bit floor for unique sequence. An
//! lroverlap consensus over a real genome is mostly *unique* sequence, so there
//! is no long-range repeat structure for an LZ to find — unlike reorder's
//! reference, which is built by clustering reads and therefore carries
//! near-duplicate contigs. Same codec, different input statistics.
//!
//! This also bounds the whole idea: the shared-reference and per-block layouts
//! are within **4 bytes** of each other on this data, so saving X bytes on the
//! frame saves X−4 bytes on the archive. −3.07% of the frame is 174 KB, i.e.
//! **0.46%** of the sequence stream. The "15%" ceiling only obtains if the frame
//! were free.
//!
//! The packed side is measured with `fqxv-reorder`'s public `GlobalReference`
//! (`from_lens_seq` + `encode_packed`), which is the exact 2-bit+LZMA codec in
//! question. That is a measurement convenience only — lroverlap must not depend
//! on reorder, so shipping this would mean hoisting the codec into a crate both
//! can use. At 0.46% that refactor does not pay for itself; kept here so the
//! negative result is reproducible on other data (amplicon or higher-coverage
//! references may still have the repeat structure this one lacks).
//!
//! Usage: cargo run --release -p fqxv --example lroverlap_refcodec_ab -- <fastq> [ont|hifi]

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use fqxv_lroverlap::{EncodeOpts, Sketch, build_reference};
use fqxv_reorder::GlobalReference;

fn mib(n: usize) -> f64 {
    n as f64 / (1024.0 * 1024.0)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: lroverlap_refcodec_ab <fastq> [ont|hifi]");
    let sketch = match args.next().as_deref() {
        Some("hifi") => Sketch::hifi(),
        _ => Sketch::ont(),
    };

    let mut reader =
        noodles_fastq::io::Reader::new(BufReader::new(File::open(&path).expect("open fastq")));
    let mut rec = noodles_fastq::Record::default();
    let (mut lens, mut seq) = (Vec::<u32>::new(), Vec::<u8>::new());
    while reader.read_record(&mut rec).expect("read record") != 0 {
        let s: &[u8] = rec.sequence();
        lens.push(s.len() as u32);
        seq.extend_from_slice(s);
    }
    eprintln!("{path}: {} reads, {} bases", lens.len(), seq.len());

    let t = Instant::now();
    let opts = EncodeOpts {
        sketch,
        ..EncodeOpts::default()
    };
    let reference = build_reference(&lens, &seq, &opts).expect("build reference");
    eprintln!(
        "reference: {} contigs, {} bases, built in {:.1}s",
        reference.len(),
        reference.total_bases(),
        t.elapsed().as_secs_f64()
    );
    if reference.is_empty() {
        eprintln!("no contigs — nothing to compare");
        return;
    }

    // A: what ships today — rANS order-0/1 over the raw consensus bases.
    let cur = reference.encode().expect("encode current");

    // B: 2-bit-pack + LZMA over the same contigs, via reorder's codec.
    let ref_lens: Vec<u32> = reference
        .consensi()
        .iter()
        .map(|c| c.len() as u32)
        .collect();
    let ref_bases: Vec<u8> = reference.consensi().concat();
    let gref = GlobalReference::from_lens_seq(&ref_lens, ref_bases).expect("global reference");
    let packed = gref.encode_packed().expect("encode packed");
    let lzma = gref.encode_lzma().expect("encode lzma");

    // Round-trip both alternatives so a size win is not a correctness loss.
    let rt_packed = GlobalReference::decode_packed(&packed).expect("decode packed");
    assert_eq!(
        rt_packed.raw_bases(),
        gref.raw_bases(),
        "packed round-trip mismatch"
    );
    assert_eq!(
        rt_packed.contig_lens(),
        gref.contig_lens(),
        "packed lens mismatch"
    );

    let bases = reference.total_bases().max(1);
    let bpb = |n: usize| (n as f64 * 8.0) / bases as f64;
    let delta = |n: usize| (n as f64 - cur.len() as f64) / cur.len() as f64 * 100.0;

    println!("\nreference frame, {bases} consensus bases:");
    println!(
        "  {:<22} {:>12} {:>10} {:>9}",
        "coder", "bytes", "b/base", "vs now"
    );
    println!("  {:-<56}", "");
    println!(
        "  {:<22} {:>12} {:>10.3} {:>8.2}%",
        "rANS raw (now)",
        cur.len(),
        bpb(cur.len()),
        0.0
    );
    println!(
        "  {:<22} {:>12} {:>10.3} {:>8.2}%",
        "LZMA raw",
        lzma.len(),
        bpb(lzma.len()),
        delta(lzma.len())
    );
    println!(
        "  {:<22} {:>12} {:>10.3} {:>8.2}%",
        "2bit+LZMA (refpack)",
        packed.len(),
        bpb(packed.len()),
        delta(packed.len())
    );
    println!(
        "\n  frame saving vs now: {:.2} MiB",
        mib(cur.len()) - mib(packed.len())
    );
}
