//! Per-stream cost breakdown of the three FASTQ codecs, on real data.
//!
//! Reads a plain FASTQ, splits it into the three streams the container codes
//! independently — names (fqxv-tokenizer → rANS), sequence (fqxv-seq → range
//! coder), quality (fqxv-fqzcomp → range coder) — and times each codec's
//! encode + decode single-threaded. The point is the *relative* split: how much
//! of the work lands on the rANS path (GPU-parallelizable) versus the serial
//! adaptive range coder (not). Parallelism scales all three together, so the
//! single-thread ratio is the Amdahl ceiling for GPU-accelerating rANS.
//!
//! Usage: `cargo run --release -p fqxv --example stream_split -- <fastq> [order]`

use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use fqxv_fqzcomp::QualityBinning;

struct Row {
    stream: &'static str,
    codec: &'static str,
    family: &'static str,
    raw: usize,
    comp: usize,
    enc_s: f64,
    dec_s: f64,
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: stream_split <fastq> [seq_order]");
    let order: usize = args.next().map_or(11, |s| s.parse().expect("order"));

    // Collect the three streams exactly as the codecs consume them.
    let mut reader =
        noodles_fastq::io::Reader::new(BufReader::new(File::open(&path).expect("open fastq")));
    let mut rec = noodles_fastq::Record::default();
    let mut headers: Vec<Vec<u8>> = Vec::new();
    let mut lens: Vec<u32> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut quals: Vec<u8> = Vec::new();
    while reader.read_record(&mut rec).expect("read record") != 0 {
        let name: &[u8] = rec.name().as_ref();
        let desc: &[u8] = rec.description().as_ref();
        let mut h = name.to_vec();
        if !desc.is_empty() {
            h.push(b' ');
            h.extend_from_slice(desc);
        }
        headers.push(h);
        let s: &[u8] = rec.sequence();
        let q: &[u8] = rec.quality_scores();
        lens.push(s.len() as u32);
        seq.extend_from_slice(s);
        quals.extend_from_slice(q);
    }
    let name_refs: Vec<&[u8]> = headers.iter().map(|h| h.as_slice()).collect();
    let names_raw: usize = headers.iter().map(|h| h.len()).sum();
    let n_reads = lens.len();
    eprintln!(
        "{path}: {n_reads} reads, seq order {order}\n\
         names {names_raw} B, seq {} B, qual {} B",
        seq.len(),
        quals.len()
    );

    let time = |f: &dyn Fn() -> Vec<u8>| {
        let t = Instant::now();
        let out = f();
        (out, t.elapsed().as_secs_f64())
    };

    // NAMES — tokenizer, rANS entropy backend.
    let (names_c, names_enc) = time(&|| fqxv_tokenizer::encode(&name_refs).expect("names enc"));
    let t = Instant::now();
    let names_d = fqxv_tokenizer::decode(&names_c).expect("names dec");
    let names_dec = t.elapsed().as_secs_f64();
    assert_eq!(names_d, headers, "names round-trip");

    // SEQUENCE — order-k adaptive context model, range coded.
    let (seq_c, seq_enc) = time(&|| fqxv_seq::encode(&lens, &seq, order).expect("seq enc"));
    let t = Instant::now();
    let (_, seq_d) = fqxv_seq::decode(&seq_c).expect("seq dec");
    let seq_dec = t.elapsed().as_secs_f64();
    assert_eq!(seq_d, seq, "seq round-trip");

    // QUALITY — fqzcomp context model, range coded (lossless).
    let (qual_c, qual_enc) =
        time(&|| fqxv_fqzcomp::encode(&lens, &quals, QualityBinning::Lossless).expect("qual enc"));
    let t = Instant::now();
    let (_, qual_d) = fqxv_fqzcomp::decode(&qual_c).expect("qual dec");
    let qual_dec = t.elapsed().as_secs_f64();
    assert_eq!(qual_d, quals, "qual round-trip");

    let rows = [
        Row {
            stream: "names",
            codec: "tokenizer",
            family: "rANS",
            raw: names_raw,
            comp: names_c.len(),
            enc_s: names_enc,
            dec_s: names_dec,
        },
        Row {
            stream: "sequence",
            codec: "seq",
            family: "range",
            raw: seq.len(),
            comp: seq_c.len(),
            enc_s: seq_enc,
            dec_s: seq_dec,
        },
        Row {
            stream: "quality",
            codec: "fqzcomp",
            family: "range",
            raw: quals.len(),
            comp: qual_c.len(),
            enc_s: qual_enc,
            dec_s: qual_dec,
        },
    ];

    let tot_enc: f64 = rows.iter().map(|r| r.enc_s).sum();
    let tot_dec: f64 = rows.iter().map(|r| r.dec_s).sum();
    let mbps = |bytes: usize, s: f64| (bytes as f64 / 1e6) / s;

    println!(
        "\n{:<9} {:<10} {:<6} {:>10} {:>10} {:>6} {:>9} {:>9} {:>7} {:>7}",
        "stream",
        "codec",
        "family",
        "raw MB",
        "comp MB",
        "ratio",
        "enc MB/s",
        "dec MB/s",
        "enc %",
        "dec %"
    );
    for r in &rows {
        println!(
            "{:<9} {:<10} {:<6} {:>10.1} {:>10.1} {:>5.2}x {:>9.0} {:>9.0} {:>6.1}% {:>6.1}%",
            r.stream,
            r.codec,
            r.family,
            r.raw as f64 / 1e6,
            r.comp as f64 / 1e6,
            r.raw as f64 / r.comp as f64,
            mbps(r.raw, r.enc_s),
            mbps(r.raw, r.dec_s),
            100.0 * r.enc_s / tot_enc,
            100.0 * r.dec_s / tot_dec,
        );
    }

    let rans_enc: f64 = rows
        .iter()
        .filter(|r| r.family == "rANS")
        .map(|r| r.enc_s)
        .sum();
    let rans_dec: f64 = rows
        .iter()
        .filter(|r| r.family == "rANS")
        .map(|r| r.dec_s)
        .sum();
    println!(
        "\nAmdahl split (share of codec time):\n  \
         rANS  (GPU-parallelizable): enc {:.1}%  dec {:.1}%\n  \
         range (serial arithmetic) : enc {:.1}%  dec {:.1}%",
        100.0 * rans_enc / tot_enc,
        100.0 * rans_dec / tot_dec,
        100.0 * (tot_enc - rans_enc) / tot_enc,
        100.0 * (tot_dec - rans_dec) / tot_dec,
    );
    println!(
        "\nCeiling: even an infinitely fast GPU rANS caps end-to-end codec speedup\n\
         at enc {:.2}x / dec {:.2}x (1 / range-fraction). Production seq uses\n\
         fqxv-reorder (clustering + range) — even more serial than plain seq here.",
        tot_enc / (tot_enc - rans_enc),
        tot_dec / (tot_dec - rans_dec),
    );
}
