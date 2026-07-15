//! fqdigest — order-independent content digest of FASTQ, for the corpus
//! round-trip check (corpus.sh). Replaces a slow `awk | sort | md5sum`.
//!
//! The check needs to answer "do two FASTQs hold the same multiset of
//! (name, sequence, quality) records, ignoring order?" The awk version sorted
//! every record (O(n log n)) and, for binned modes, walked each quality byte in
//! interpreted awk. Here we hash each record and *sum* the hashes: addition is
//! commutative, so the total is independent of read order, and it preserves
//! multiplicity (a duplicated record adds twice), giving a true multiset digest
//! in one O(n) streaming pass with no sort and no whole-file buffering.
//!
//! The `+` separator line (record line 3) is excluded, matching fqxv's
//! documented `+`-normalization. With `--bin`, each quality byte is passed
//! through fqxv's bin table first — the expected content of a correct lossy
//! round-trip (mirrors QualityBinning::apply in fqxv-fqzcomp).
//!
//! Usage: fqdigest [--bin bin8|bin4|bin2] [--no-qual] [--no-names] FILE...  (or `-`/none = stdin)
//! `--no-qual` hashes (name, sequence) only, for tools that drop quality.
//! `--no-names` hashes (sequence, quality) only, for tools that renumber reads.
//! Output: a 32-hex-digit digest. Same digest ⇔ same record multiset.
//!
//! Not fqxv's own hash — a standalone harness tool. Build:
//!   rustc -O bench/fqdigest.rs -o <path>/fqdigest

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

/// Illumina quality bin tables. Must mirror QualityBinning::apply.
fn bin_q(scheme: Bin, q: u8) -> u8 {
    match scheme {
        Bin::None => q,
        Bin::Bin8 => match q {
            0..=1 => q,
            2..=9 => 6,
            10..=19 => 15,
            20..=24 => 22,
            25..=29 => 27,
            30..=34 => 33,
            35..=39 => 37,
            _ => 40,
        },
        Bin::Bin4 => match q {
            0..=2 => 2,
            3..=17 => 12,
            18..=29 => 24,
            _ => 40,
        },
        Bin::Bin2 => {
            if q <= 24 {
                15
            } else {
                37
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Bin {
    None,
    Bin8,
    Bin4,
    Bin2,
}

// wyhash-style multiply-fold primitives. Fast (multi-GB/s) with good avalanche;
// std's SipHash was ~10x slower here and per-`write` call overhead dominated.
const P0: u64 = 0xa076_1d64_78bd_642f;
const P1: u64 = 0xe703_7ed1_a0b4_28db;
const P2: u64 = 0x8ebc_6af0_9c88_c6e3;

#[inline(always)]
fn mum(a: u64, b: u64) -> u64 {
    let r = (a as u128).wrapping_mul(b as u128);
    (r as u64) ^ ((r >> 64) as u64)
}

/// Fold one byte slice into a running hash state (8 bytes at a time + tail).
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

/// One 64-bit lane of the multiset sum: chain-fold the three record fields into
/// a seed derived from `domain`. Chaining keeps field boundaries significant
/// without materializing a concatenated buffer.
#[inline(always)]
fn lane(domain: u64, name: &[u8], seq: &[u8], qual: &[u8]) -> u64 {
    let mut h = domain ^ P0;
    h = fold(h, name);
    h = fold(h ^ P1, seq);
    h = fold(h ^ P2, qual);
    mum(h ^ P0, P1)
}

fn digest_reader<R: Read>(
    r: R,
    scheme: Bin,
    no_qual: bool,
    no_names: bool,
    lo: &mut u64,
    hi: &mut u64,
) -> io::Result<()> {
    // Stream line-by-line with reused buffers: memory is bounded to four records'
    // worth regardless of file size (a multi-GB run must not be slurped whole).
    let mut br = BufReader::with_capacity(1 << 20, r);
    let (mut name, mut seq, mut plus, mut qual) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    loop {
        name.clear();
        if br.read_until(b'\n', &mut name)? == 0 {
            break; // clean EOF at a record boundary
        }
        seq.clear();
        plus.clear();
        qual.clear();
        // A short final record (truncated file) reads 0 here and is dropped.
        if br.read_until(b'\n', &mut seq)? == 0 {
            break;
        }
        br.read_until(b'\n', &mut plus)?; // line 3: excluded
        if br.read_until(b'\n', &mut qual)? == 0 {
            break;
        }

        let s = trim(&seq);
        if scheme != Bin::None && !no_qual {
            let end = trim_len(&qual);
            for b in &mut qual[..end] {
                *b = bin_q(scheme, b.wrapping_sub(33)) + 33;
            }
        }
        // `--no-qual`: hash (name, sequence) only, for tools that drop quality.
        // `--no-names`: hash (sequence, quality) only, for tools that renumber
        // reads (e.g. `fqxv --order shuffle`, SPRING-style) — the round-trip is
        // compared on the retained content.
        let n: &[u8] = if no_names { &[] } else { trim(&name) };
        let q: &[u8] = if no_qual { &[] } else { trim(&qual) };
        *lo = lo.wrapping_add(lane(1, n, s, q));
        *hi = hi.wrapping_add(lane(2, n, s, q));
    }
    Ok(())
}

/// Length of a line with any trailing `\n` and `\r` removed.
#[inline]
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

#[inline]
fn trim(l: &[u8]) -> &[u8] {
    &l[..trim_len(l)]
}

fn main() -> io::Result<()> {
    let mut scheme = Bin::None;
    let mut no_qual = false;
    let mut no_names = false;
    let mut files: Vec<String> = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--bin" => {
                scheme = match args.next().as_deref() {
                    Some("bin8") => Bin::Bin8,
                    Some("bin4") => Bin::Bin4,
                    Some("bin2") => Bin::Bin2,
                    Some("none") | Some("lossless") => Bin::None,
                    other => {
                        eprintln!("fqdigest: unknown --bin {:?}", other);
                        std::process::exit(2);
                    }
                };
            }
            "--no-qual" => no_qual = true,
            "--no-names" => no_names = true,
            "-h" | "--help" => {
                eprintln!(
                    "usage: fqdigest [--bin bin8|bin4|bin2] [--no-qual] [--no-names] FILE...  (- or none = stdin)"
                );
                return Ok(());
            }
            _ => files.push(a),
        }
    }

    let (mut lo, mut hi) = (0u64, 0u64);
    if files.is_empty() || files == ["-"] {
        digest_reader(io::stdin().lock(), scheme, no_qual, no_names, &mut lo, &mut hi)?;
    } else {
        for f in &files {
            let file = File::open(f)?;
            digest_reader(file, scheme, no_qual, no_names, &mut lo, &mut hi)?;
        }
    }

    let mut out = io::stdout().lock();
    writeln!(out, "{:016x}{:016x}", lo, hi)
}
