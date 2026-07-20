//! fqsim — synthetic FASTQ generator with per-platform length, error and
//! quality models, for test data that exercises the real codec paths.
//!
//! The point is *structured* redundancy. Reads drawn independently at random
//! share nothing, so they compress to near the entropy floor and tell you
//! nothing about the codecs that exist to exploit cross-read similarity —
//! fqxv-reorder's minimizer clustering and fqxv-lroverlap's overlap coding both
//! come out looking pointless. So fqsim builds a random reference genome once,
//! then samples reads *from* it at a chosen coverage. Overlapping reads then
//! carry genuine shared sequence, the reorder and overlap paths have something
//! to find, and ratios move the way they do on real data.
//!
//! Quality is correlated with error, not independent of it: a base that was
//! mutated gets a low score drawn from the bottom of the platform's range.
//! That is what makes the quality stream compressible in the way fqzcomp
//! assumes, and an uncorrelated model understates it.
//!
//! Platforms set read length, error mix, quality model and name format:
//!
//! * `novaseq` — 150 bp, low substitution rate, 4-level binned quality
//!   ({2,12,23,37}, the RTA3 bins). The binning is the point: it is why
//!   NovaSeq quality compresses so much better than HiSeq.
//! * `hiseq`   — 150 bp, low substitution rate, full 2..41 continuous quality.
//! * `ont`     — lognormal ~8 kb, ~5% error dominated by indels and
//!   homopolymer miscalls, quality centred near Q12.
//! * `hifi`    — lognormal ~15 kb, ~0.1% error, quality near Q32.
//!
//! Long-read error is homopolymer-aware: inside a run of identical bases the
//! indel rate is raised sharply, which is the dominant real ONT/PacBio failure
//! mode and the thing homopolymer-context quality models key on.
//!
//! Output is plain FASTQ. Pipe to `bgzip`/`gzip` for compressed input; keeping
//! the tool dependency-free is worth more than built-in gzip.
//!
//! Deterministic: the same `--seed` gives byte-identical output, so a
//! regression run is reproducible without staging data files.
//!
//! Usage:
//!   fqsim --platform novaseq --reads 1000000 -o reads.fastq
//!   fqsim --platform novaseq --reads 1000000 --paired sample   # sample_1/_2.fastq
//!   fqsim --platform ont --reads 20000 --genome 5000000 -o ont.fastq
//!   fqsim --platform hifi --reads 50000 --coverage 30 -o hifi.fastq
//!
//! Options:
//!   --platform P     novaseq | hiseq | ont | hifi        (default novaseq)
//!   --reads N        number of reads (spots, when --paired)
//!   --genome BP      reference length (default 1000000; ignored with --coverage)
//!   --coverage X     size the genome so reads give ~X-fold coverage instead
//!   --len L          override mean read length
//!   --seed S         PRNG seed (default 1)
//!   --sub-rate R     override substitution rate (per base)
//!   --ins-rate R     override insertion rate
//!   --del-rate R     override deletion rate
//!   --frag-len L     mean paired fragment length (default 350)
//!   -o FILE          output ("-" or omitted = stdout)
//!   --paired PREFIX  write PREFIX_1.fastq and PREFIX_2.fastq
//!
//! Build: rustc -O -o fqsim bench/tools/fqsim.rs

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::process::ExitCode;

// ---------------------------------------------------------------------------
// PRNG: splitmix64 to seed, xoshiro256++ for the stream. Deliberately
// self-contained — a data generator that pulls in a crate graph is a worse
// tool, and neither algorithm is subtle enough to be worth a dependency.
// ---------------------------------------------------------------------------

struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            x ^ (x >> 31)
        };
        Self {
            s: [next(), next(), next(), next()],
        }
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let r = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }

    /// Uniform in [0, 1).
    #[inline]
    fn f64(&mut self) -> f64 {
        // 53 significant bits, the most an f64 mantissa holds exactly.
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in [0, n).
    #[inline]
    fn below(&mut self, n: u64) -> u64 {
        // Lemire's multiply-shift. Biased by at most 2^-64 for our n, which is
        // far below anything a compression benchmark could resolve.
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }

    /// Standard normal, via Box-Muller (one of the pair; the other is dropped —
    /// caching it would save a call we are not bottlenecked on).
    fn normal(&mut self) -> f64 {
        let u1 = self.f64().max(f64::MIN_POSITIVE);
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Lognormal length with the given mean, clamped to a sane floor. `sigma`
    /// is on the log scale: long-read length distributions are heavy-tailed,
    /// and a normal would put implausible mass near zero.
    fn lognormal_len(&mut self, mean: f64, sigma: f64) -> usize {
        let mu = mean.ln() - 0.5 * sigma * sigma;
        let v = (mu + sigma * self.normal()).exp();
        v.max(64.0) as usize
    }
}

// ---------------------------------------------------------------------------
// Platform profiles
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum Platform {
    NovaSeq,
    HiSeq,
    Ont,
    HiFi,
}

#[derive(Clone, Copy, PartialEq)]
enum QualModel {
    /// NovaSeq RTA3: every score is one of four bins.
    Binned4,
    /// Continuous range, degrading toward the 3' end.
    Continuous,
    /// Long-read: centred on a mean with a wide spread, no positional ramp.
    LongRead,
}

struct Profile {
    name: &'static str,
    mean_len: f64,
    /// Log-scale spread; 0 means fixed-length (short-read).
    len_sigma: f64,
    sub_rate: f64,
    ins_rate: f64,
    del_rate: f64,
    /// Multiplier on indel rates inside a homopolymer run — the dominant
    /// long-read error mode, and ~1.0 (off) for short reads.
    homopolymer_boost: f64,
    qual: QualModel,
    /// Centre of the quality distribution; meaning depends on `qual`.
    qual_mean: f64,
    qual_sd: f64,
}

impl Platform {
    fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "novaseq" => Some(Self::NovaSeq),
            "hiseq" => Some(Self::HiSeq),
            "ont" | "nanopore" => Some(Self::Ont),
            "hifi" | "pacbio" => Some(Self::HiFi),
            _ => None,
        }
    }

    fn profile(self) -> Profile {
        match self {
            Self::NovaSeq => Profile {
                name: "novaseq",
                mean_len: 150.0,
                len_sigma: 0.0,
                sub_rate: 0.001,
                ins_rate: 0.000_02,
                del_rate: 0.000_02,
                homopolymer_boost: 1.0,
                qual: QualModel::Binned4,
                qual_mean: 0.0,
                qual_sd: 0.0,
            },
            Self::HiSeq => Profile {
                name: "hiseq",
                mean_len: 150.0,
                len_sigma: 0.0,
                sub_rate: 0.002,
                ins_rate: 0.000_05,
                del_rate: 0.000_05,
                homopolymer_boost: 1.0,
                qual: QualModel::Continuous,
                qual_mean: 34.0,
                qual_sd: 4.0,
            },
            Self::Ont => Profile {
                name: "ont",
                mean_len: 8000.0,
                len_sigma: 0.9,
                sub_rate: 0.020,
                ins_rate: 0.015,
                del_rate: 0.015,
                homopolymer_boost: 6.0,
                qual: QualModel::LongRead,
                qual_mean: 12.0,
                qual_sd: 4.0,
            },
            Self::HiFi => Profile {
                name: "hifi",
                mean_len: 15000.0,
                len_sigma: 0.35,
                sub_rate: 0.000_4,
                ins_rate: 0.000_3,
                del_rate: 0.000_3,
                homopolymer_boost: 4.0,
                qual: QualModel::LongRead,
                qual_mean: 32.0,
                qual_sd: 5.0,
            },
        }
    }
}

/// NovaSeq RTA3 quality bins.
const BINS4: [u8; 4] = [2, 12, 23, 37];

const BASES: [u8; 4] = *b"ACGT";

// ---------------------------------------------------------------------------
// Generation
// ---------------------------------------------------------------------------

/// A random reference the reads are sampled from. Uniform base composition —
/// real genomes have structure (repeats, GC skew) that would make the sequence
/// stream *more* compressible, so this is the conservative choice: a codec that
/// wins here wins by exploiting cross-read overlap, not reference structure.
fn make_genome(rng: &mut Rng, len: usize) -> Vec<u8> {
    let mut g = vec![0u8; len];
    for b in g.iter_mut() {
        *b = BASES[rng.below(4) as usize];
    }
    g
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Copy `genome[start..start+len]` into `out`, reverse-complemented when `rc`.
/// Wraps at the end of the genome so a long read near the boundary still gets a
/// full-length draw rather than a short one that would skew the length model.
fn draw_template(genome: &[u8], start: usize, len: usize, rc: bool, out: &mut Vec<u8>) {
    out.clear();
    let n = genome.len();
    if rc {
        for i in (0..len).rev() {
            out.push(complement(genome[(start + i) % n]));
        }
    } else {
        for i in 0..len {
            out.push(genome[(start + i) % n]);
        }
    }
}

/// Apply the platform error model to `template`, writing bases to `seq` and the
/// matching quality scores to `qual`.
///
/// Errors and quality are produced together on purpose. A substituted or
/// inserted base gets a score from the bottom of the range while a correct one
/// gets a draw from the platform's normal distribution, so the quality stream
/// carries real information about the sequence — the correlation fqzcomp's
/// context model (and the homopolymer-aware long-read context) is built to use.
#[allow(clippy::too_many_arguments)]
fn apply_errors(
    rng: &mut Rng,
    p: &Profile,
    template: &[u8],
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) {
    seq.clear();
    qual.clear();
    let mut run_base = 0u8;
    let mut run_len = 0usize;

    for (pos, &base) in template.iter().enumerate() {
        // Track the homopolymer run we are inside, to raise the indel rate the
        // way real long-read basecallers miscall run lengths.
        if base == run_base {
            run_len += 1;
        } else {
            run_base = base;
            run_len = 1;
        }
        let boost = if run_len >= 3 { p.homopolymer_boost } else { 1.0 };

        // Deletion: emit nothing for this template base.
        if rng.f64() < p.del_rate * boost {
            continue;
        }

        // Insertion: emit a random base before the real one, scored low.
        if rng.f64() < p.ins_rate * boost {
            seq.push(BASES[rng.below(4) as usize]);
            qual.push(error_qual(rng, p));
        }

        if rng.f64() < p.sub_rate {
            // Substitution: any of the three other bases, equally likely.
            let orig = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            let alt = (orig + 1 + rng.below(3) as usize) % 4;
            seq.push(BASES[alt]);
            qual.push(error_qual(rng, p));
        } else {
            seq.push(base);
            qual.push(good_qual(rng, p, pos, template.len()));
        }
    }

    // A read that lost everything to deletions is not a useful record; a single
    // base keeps the FASTQ well-formed.
    if seq.is_empty() {
        seq.push(b'A');
        qual.push(error_qual(rng, p));
    }
}

/// Quality for a base the model got wrong: the bottom bin, or a low draw.
fn error_qual(rng: &mut Rng, p: &Profile) -> u8 {
    let q = match p.qual {
        QualModel::Binned4 => BINS4[rng.below(2) as usize] as f64,
        QualModel::Continuous => 2.0 + rng.f64() * 10.0,
        QualModel::LongRead => (p.qual_mean * 0.35).max(2.0) + rng.f64() * 4.0,
    };
    phred_to_char(q)
}

/// Quality for a correctly-called base. Short-read models ramp down toward the
/// 3' end, which is the dominant positional effect on Illumina and the reason
/// position is part of fqzcomp's context.
fn good_qual(rng: &mut Rng, p: &Profile, pos: usize, len: usize) -> u8 {
    let frac = if len > 1 {
        pos as f64 / (len - 1) as f64
    } else {
        0.0
    };
    let q = match p.qual {
        QualModel::Binned4 => {
            // Later cycles slide probability mass toward the lower bins.
            let r = rng.f64() + frac * 0.35;
            if r < 0.72 {
                BINS4[3] as f64
            } else if r < 0.92 {
                BINS4[2] as f64
            } else if r < 0.99 {
                BINS4[1] as f64
            } else {
                BINS4[0] as f64
            }
        }
        QualModel::Continuous => p.qual_mean - frac * 8.0 + rng.normal() * p.qual_sd,
        QualModel::LongRead => p.qual_mean + rng.normal() * p.qual_sd,
    };
    phred_to_char(q)
}

#[inline]
fn phred_to_char(q: f64) -> u8 {
    // Clamp to the Sanger range fqxv's quality codec expects (0..=93).
    let q = q.round().clamp(0.0, 93.0) as u8;
    b'!' + q
}

/// Platform-shaped read name. Name structure is not cosmetic here: the
/// tokenizer codes names by column, so a realistic field layout is what makes
/// the names stream behave like it does on real data.
fn write_name(out: &mut Vec<u8>, p: &Profile, plat: Platform, idx: u64, mate: usize, rng: &mut Rng) {
    out.clear();
    match plat {
        Platform::NovaSeq | Platform::HiSeq => {
            let lane = 1 + (idx % 4);
            let tile = 1101 + (idx / 4) % 78;
            let x = 1000 + rng.below(30000);
            let y = 1000 + rng.below(30000);
            let inst = if plat == Platform::NovaSeq {
                "A00123:45:HXXXXDSXX"
            } else {
                "D00360:95:C7T3KANXX"
            };
            let _ = write!(
                out,
                "{inst}:{lane}:{tile}:{x}:{y} {mate}:N:0:ATCACG",
                mate = mate.max(1)
            );
        }
        Platform::Ont => {
            // A UUID-shaped read id plus the run metadata MinKNOW emits.
            let a = rng.next_u64();
            let b = rng.next_u64();
            let _ = write!(
                out,
                "{:08x}-{:04x}-{:04x}-{:04x}-{:012x} runid=fqsim{} read={} ch={} start_time=2026-07-19T00:00:00Z",
                (a >> 32) as u32,
                (a >> 16) as u16,
                a as u16,
                (b >> 48) as u16,
                b & 0xffff_ffff_ffff,
                p.name.len(),
                idx,
                1 + (idx % 512),
            );
        }
        Platform::HiFi => {
            let _ = write!(out, "m64011_190830_220126/{}/ccs", 1000 + idx);
        }
    }
}

fn write_record(w: &mut impl Write, name: &[u8], seq: &[u8], qual: &[u8]) -> io::Result<()> {
    w.write_all(b"@")?;
    w.write_all(name)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

struct Args {
    platform: Platform,
    reads: u64,
    genome: usize,
    coverage: Option<f64>,
    len: Option<f64>,
    seed: u64,
    sub: Option<f64>,
    ins: Option<f64>,
    del: Option<f64>,
    frag: f64,
    out: Option<String>,
    paired: Option<String>,
}

fn usage() -> &'static str {
    "usage: fqsim [--platform novaseq|hiseq|ont|hifi] --reads N [--genome BP] [--coverage X]\n\
     \x20            [--len L] [--seed S] [--sub-rate R] [--ins-rate R] [--del-rate R]\n\
     \x20            [--frag-len L] [-o FILE] [--paired PREFIX]"
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        platform: Platform::NovaSeq,
        reads: 0,
        genome: 1_000_000,
        coverage: None,
        len: None,
        seed: 1,
        sub: None,
        ins: None,
        del: None,
        frag: 350.0,
        out: None,
        paired: None,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let next = |i: &mut usize, what: &str| -> Result<String, String> {
            *i += 1;
            argv.get(*i)
                .cloned()
                .ok_or_else(|| format!("{what} needs a value"))
        };
        match argv[i].as_str() {
            "--platform" => {
                let v = next(&mut i, "--platform")?;
                a.platform =
                    Platform::parse(&v).ok_or_else(|| format!("unknown platform '{v}'"))?;
            }
            "--reads" => a.reads = next(&mut i, "--reads")?.parse().map_err(|_| "bad --reads")?,
            "--genome" => {
                a.genome = next(&mut i, "--genome")?
                    .parse()
                    .map_err(|_| "bad --genome")?
            }
            "--coverage" => {
                a.coverage = Some(
                    next(&mut i, "--coverage")?
                        .parse()
                        .map_err(|_| "bad --coverage")?,
                )
            }
            "--len" => a.len = Some(next(&mut i, "--len")?.parse().map_err(|_| "bad --len")?),
            "--seed" => a.seed = next(&mut i, "--seed")?.parse().map_err(|_| "bad --seed")?,
            "--sub-rate" => {
                a.sub = Some(
                    next(&mut i, "--sub-rate")?
                        .parse()
                        .map_err(|_| "bad --sub-rate")?,
                )
            }
            "--ins-rate" => {
                a.ins = Some(
                    next(&mut i, "--ins-rate")?
                        .parse()
                        .map_err(|_| "bad --ins-rate")?,
                )
            }
            "--del-rate" => {
                a.del = Some(
                    next(&mut i, "--del-rate")?
                        .parse()
                        .map_err(|_| "bad --del-rate")?,
                )
            }
            "--frag-len" => {
                a.frag = next(&mut i, "--frag-len")?
                    .parse()
                    .map_err(|_| "bad --frag-len")?
            }
            "-o" | "--out" => a.out = Some(next(&mut i, "-o")?),
            "--paired" => a.paired = Some(next(&mut i, "--paired")?),
            "-h" | "--help" => return Err(usage().to_string()),
            other => return Err(format!("unknown argument '{other}'\n{}", usage())),
        }
        i += 1;
    }
    if a.reads == 0 {
        return Err(format!("--reads is required\n{}", usage()));
    }
    Ok(a)
}

fn open_out(path: Option<&str>) -> io::Result<Box<dyn Write>> {
    match path {
        None | Some("-") => Ok(Box::new(BufWriter::with_capacity(
            1 << 20,
            io::stdout().lock(),
        ))),
        Some(p) => Ok(Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?))),
    }
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("fqsim: {e}");
            return ExitCode::from(2);
        }
    };
    match run(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("fqsim: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Args) -> io::Result<()> {
    let mut p = args.platform.profile();
    if let Some(l) = args.len {
        p.mean_len = l;
    }
    if let Some(r) = args.sub {
        p.sub_rate = r;
    }
    if let Some(r) = args.ins {
        p.ins_rate = r;
    }
    if let Some(r) = args.del {
        p.del_rate = r;
    }

    let paired = args.paired.is_some();
    // Paired mode draws a fragment per spot and reads both ends, so a spot
    // consumes two reads' worth of bases.
    let bases_per_spot = if paired {
        p.mean_len * 2.0
    } else {
        p.mean_len
    };
    let genome_len = match args.coverage {
        Some(x) if x > 0.0 => ((args.reads as f64 * bases_per_spot) / x).max(1000.0) as usize,
        _ => args.genome,
    };

    let mut rng = Rng::new(args.seed);
    let genome = make_genome(&mut rng, genome_len);

    let (mut w1, mut w2): (Box<dyn Write>, Option<Box<dyn Write>>) = match &args.paired {
        Some(prefix) => (
            open_out(Some(&format!("{prefix}_1.fastq")))?,
            Some(open_out(Some(&format!("{prefix}_2.fastq")))?),
        ),
        None => (open_out(args.out.as_deref())?, None),
    };

    let mut template = Vec::new();
    let mut seq = Vec::new();
    let mut qual = Vec::new();
    let mut name = Vec::new();

    for idx in 0..args.reads {
        if let Some(w2) = w2.as_mut() {
            // One fragment, both ends: R1 forward from the start, R2 reverse-
            // complemented from the end — the orientation real paired-end data
            // has, and what fqxv's interleaving and the reorder codec's
            // revcomp-aware clustering expect to see.
            let frag = p.mean_len.max(args.frag) as usize;
            let start = rng.below(genome.len() as u64) as usize;
            let rlen = p.mean_len as usize;

            draw_template(&genome, start, rlen, false, &mut template);
            apply_errors(&mut rng, &p, &template, &mut seq, &mut qual);
            write_name(&mut name, &p, args.platform, idx, 1, &mut rng);
            write_record(&mut w1, &name, &seq, &qual)?;

            let r2_start = start + frag.saturating_sub(rlen);
            draw_template(&genome, r2_start, rlen, true, &mut template);
            apply_errors(&mut rng, &p, &template, &mut seq, &mut qual);
            write_name(&mut name, &p, args.platform, idx, 2, &mut rng);
            write_record(w2, &name, &seq, &qual)?;
        } else {
            let rlen = if p.len_sigma > 0.0 {
                rng.lognormal_len(p.mean_len, p.len_sigma)
            } else {
                p.mean_len as usize
            };
            let start = rng.below(genome.len() as u64) as usize;
            let rc = rng.next_u64() & 1 == 0;

            draw_template(&genome, start, rlen, rc, &mut template);
            apply_errors(&mut rng, &p, &template, &mut seq, &mut qual);
            write_name(&mut name, &p, args.platform, idx, 1, &mut rng);
            write_record(&mut w1, &name, &seq, &qual)?;
        }
    }

    w1.flush()?;
    if let Some(w2) = w2.as_mut() {
        w2.flush()?;
    }

    let cov = (args.reads as f64 * bases_per_spot) / genome_len as f64;
    eprintln!(
        "fqsim: {} reads, platform {}, genome {} bp (~{:.1}x coverage), seed {}",
        args.reads, p.name, genome_len, cov, args.seed
    );
    Ok(())
}
