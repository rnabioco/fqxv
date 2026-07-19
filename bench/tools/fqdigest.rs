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
//! Usage: fqdigest [--bin bin8|bin4|bin2|ont|hifi] [--no-qual] [--no-names] FILE...  (or `-`/none = stdin)
//! `--no-qual` hashes (name, sequence) only, for tools that drop quality.
//! `--no-names` hashes (sequence, quality) only, for tools that renumber reads.
//! Output: a 32-hex-digit digest. Same digest ⇔ same record multiset.
//!
//! `--distort ORIG RT` is a separate mode: per-base quality distortion of a lossy
//! round-trip `RT` vs the original `ORIG`. Prints `mae rmse pct` (mean abs error,
//! RMSE, % of bases whose quality changed), or `-1 -1 -1` if no records matched.
//! Replaces the old interpreted per-byte awk loop. It diffs in lockstep (O(1)
//! memory) when reads kept their order — the common case, since binning doesn't
//! reorder — and only falls back to a name->qual map when a reordering tool
//! (reorder-bin*, SPRING) shuffled the reads.
//!
//! Not fqxv's own hash — a standalone harness tool. Build:
//!   rustc -O bench/fqdigest.rs -o <path>/fqdigest

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

/// Quality bin tables. Must mirror fqxv-fqzcomp `QualityBinning::apply`, which
/// covers Illumina (bin8/bin4/bin2) and the long-read schemes (ONT, HiFi).
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
        // CoLoRd ONT 4-level table (representatives 3/10/18/35).
        Bin::Ont => match q {
            0..=6 => 3,
            7..=13 => 10,
            14..=25 => 18,
            _ => 35,
        },
        // CoLoRd HiFi 5-level table: as ONT, but the top Q93 symbol is kept exactly.
        Bin::Hifi => match q {
            0..=6 => 3,
            7..=13 => 10,
            14..=25 => 18,
            26..=92 => 35,
            _ => 93,
        },
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Bin {
    None,
    Bin8,
    Bin4,
    Bin2,
    Ont,
    Hifi,
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

/// Read one FASTQ record into the reused buffers. Returns false at a clean EOF
/// (or a truncated final record), matching `digest_reader`'s framing: line 3
/// (`+`) is consumed but ignored.
fn read_record<R: BufRead>(
    br: &mut R,
    name: &mut Vec<u8>,
    seq: &mut Vec<u8>,
    plus: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) -> io::Result<bool> {
    name.clear();
    if br.read_until(b'\n', name)? == 0 {
        return Ok(false);
    }
    seq.clear();
    plus.clear();
    qual.clear();
    if br.read_until(b'\n', seq)? == 0 {
        return Ok(false);
    }
    br.read_until(b'\n', plus)?;
    if br.read_until(b'\n', qual)? == 0 {
        return Ok(false);
    }
    Ok(true)
}

/// Running distortion accumulators: sum |Δ|, sum Δ², count changed, count total.
type Stats = (u64, u64, u64, u64);

/// Fold one original/round-trip quality pair (already trimmed) into the stats.
/// Phred+33 cancels in the difference, so raw bytes are compared directly.
#[inline]
fn accum(oq: &[u8], rq: &[u8], s: &mut Stats) {
    let m = oq.len().min(rq.len());
    for i in 0..m {
        let d = (rq[i] as i32 - oq[i] as i32).unsigned_abs() as u64;
        s.0 += d;
        s.1 += d * d;
        if d > 0 {
            s.2 += 1;
        }
        s.3 += 1;
    }
}

/// Fast path: when orig and rt hold the same records in the same order (the common
/// case — quality binning preserves read order), diff them in lockstep with no map
/// and O(1) memory. Returns `None` the moment a positional name mismatch or a
/// record-count mismatch shows the round-trip reordered reads, so the caller can
/// fall back to the name-matched map path (which then reproduces the awk exactly).
fn distort_lockstep(orig_path: &str, rt_path: &str) -> io::Result<Option<Stats>> {
    let mut bo = BufReader::with_capacity(1 << 20, File::open(orig_path)?);
    let mut br = BufReader::with_capacity(1 << 20, File::open(rt_path)?);
    let (mut on, mut os, mut op, mut oq) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut rn, mut rs, mut rp, mut rq) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut s: Stats = (0, 0, 0, 0);
    loop {
        let a = read_record(&mut bo, &mut on, &mut os, &mut op, &mut oq)?;
        let b = read_record(&mut br, &mut rn, &mut rs, &mut rp, &mut rq)?;
        match (a, b) {
            (false, false) => return Ok(Some(s)), // both hit EOF together, all names matched
            (true, true) => {
                if trim(&on) != trim(&rn) {
                    return Ok(None); // reordered → fall back to the map path
                }
                accum(trim(&oq), trim(&rq), &mut s);
            }
            _ => return Ok(None), // differing record counts → let the map path decide
        }
    }
}

/// Fallback: records matched by name via a name->qual map of the original — the
/// order-independent path for read-reordering lossy tools (reorder-bin*, SPRING).
/// Same memory profile as the awk it replaces, but one compiled O(bases) pass.
fn distort_mapped(orig_path: &str, rt_path: &str) -> io::Result<Stats> {
    let mut map: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    let (mut name, mut seq, mut plus, mut qual) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());

    let mut br = BufReader::with_capacity(1 << 20, File::open(orig_path)?);
    while read_record(&mut br, &mut name, &mut seq, &mut plus, &mut qual)? {
        map.insert(trim(&name).to_vec(), trim(&qual).to_vec());
    }

    let mut s: Stats = (0, 0, 0, 0);
    let mut br = BufReader::with_capacity(1 << 20, File::open(rt_path)?);
    while read_record(&mut br, &mut name, &mut seq, &mut plus, &mut qual)? {
        if let Some(oq) = map.get(trim(&name)) {
            accum(oq, trim(&qual), &mut s);
        }
    }
    Ok(s)
}

/// Per-base quality distortion of round-trip `rt_path` vs original `orig_path`.
/// Prints "mae rmse pct" (or "-1 -1 -1" if nothing matched). Tries the O(1)-memory
/// lockstep pass first (correct whenever reads kept their order) and falls back to
/// the name-matched map only when reordering is detected.
fn distort(orig_path: &str, rt_path: &str) -> io::Result<()> {
    let (sum_abs, sum_sq, changed, n) = match distort_lockstep(orig_path, rt_path)? {
        Some(s) => s,
        None => distort_mapped(orig_path, rt_path)?,
    };

    let mut out = io::stdout().lock();
    if n == 0 {
        writeln!(out, "-1 -1 -1")
    } else {
        let (nf, mae) = (n as f64, sum_abs as f64 / n as f64);
        let rmse = (sum_sq as f64 / nf).sqrt();
        let pct = 100.0 * changed as f64 / nf;
        writeln!(out, "{:.4} {:.4} {:.4}", mae, rmse, pct)
    }
}

fn main() -> io::Result<()> {
    let mut scheme = Bin::None;
    let mut no_qual = false;
    let mut no_names = false;
    let mut distort_orig: Option<String> = None;
    let mut files: Vec<String> = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--distort" => {
                distort_orig = Some(args.next().unwrap_or_else(|| {
                    eprintln!("fqdigest: --distort needs an ORIG file");
                    std::process::exit(2);
                }));
            }
            "--bin" => {
                scheme = match args.next().as_deref() {
                    Some("bin8") => Bin::Bin8,
                    Some("bin4") => Bin::Bin4,
                    Some("bin2") => Bin::Bin2,
                    Some("ont") => Bin::Ont,
                    Some("hifi") => Bin::Hifi,
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
                    "usage: fqdigest [--bin bin8|bin4|bin2|ont|hifi] [--no-qual] [--no-names] FILE...  (- or none = stdin)\n       fqdigest --distort ORIG RT   (per-base quality distortion: mae rmse pct)"
                );
                return Ok(());
            }
            _ => files.push(a),
        }
    }

    // Distortion mode: `--distort ORIG` consumes ORIG; the lone positional is RT.
    if let Some(orig) = distort_orig {
        match files.as_slice() {
            [rt] => return distort(&orig, rt),
            _ => {
                eprintln!("fqdigest: --distort ORIG takes exactly one RT file");
                std::process::exit(2);
            }
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
