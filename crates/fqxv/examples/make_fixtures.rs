//! Regenerate the golden archive fixtures in `crates/fqxv/tests/fixtures/`.
//!
//! ```text
//! cargo run --release -p fqxv --example make_fixtures -- crates/fqxv/tests/fixtures
//! ```
//!
//! The fixtures exist so a later release can prove it still decodes archives an
//! earlier one wrote — the one compatibility property the rest of the suite cannot
//! check, because every other test compares a build against itself.
//!
//! **Regenerating them is almost always the wrong fix for a failing golden test.**
//! A mismatch means this build no longer reproduces what a previous one wrote,
//! which is exactly the regression the fixtures exist to catch. Regenerate only
//! when deliberately adding a fixture, or when a format change has been made
//! knowingly and gated (see the evolution policy in `docs/design/container.md`) —
//! and in that case keep the old fixtures too, since they must still decode.
//!
//! Inputs are generated here rather than checked in: only the archives and a
//! manifest of the expected decode (length + CRC-32C) are stored, which keeps the
//! directory small enough that each input can have the coverage its codec needs.
//! A fixture whose input is too thin silently falls back to the plain order-k path
//! and pins nothing the `plain_se` fixture does not already cover.

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use fqxv::{Params, compress, compress_multi};

/// Deterministic, dependency-free PRNG (xorshift64*). Reproducibility across
/// platforms matters more here than statistical quality — the fixtures must come
/// out byte-identical wherever they are regenerated.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    /// Roughly normal via the sum of three uniforms; enough structure that the
    /// quality model has something real to model.
    fn quality(&mut self, lo: u8, hi: u8) -> u8 {
        let span = u64::from(hi - lo);
        let s = (self.next() % span + self.next() % span + self.next() % span) / 3;
        lo + s as u8
    }
}

const ACGT: [u8; 4] = *b"ACGT";

fn genome(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| ACGT[rng.below(4)]).collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn mutate(rng: &mut Rng, s: &[u8], per_mille: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len());
    for &b in s {
        match rng.next() % 1000 {
            r if r < per_mille * 6 / 10 => out.push(ACGT[rng.below(4)]),
            r if r < per_mille * 8 / 10 => {}
            r if r < per_mille => {
                out.push(b);
                out.push(ACGT[rng.below(4)]);
            }
            _ => out.push(b),
        }
    }
    out
}

fn illumina(rng: &mut Rng, n: usize, rlen: usize, reference: &[u8], mate: Option<u8>) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..n {
        let p = rng.below(reference.len() - rlen);
        let mut s = mutate(rng, &reference[p..p + rlen], 2);
        s.resize(rlen, b'A');
        let q: Vec<u8> = (0..s.len()).map(|_| 33 + rng.quality(20, 41)).collect();
        let mut name = format!(
            "@INST7:41:HFXXX:2:{}:{}:{}",
            1101 + i % 8,
            1000 + i,
            2000 + i * 3
        );
        if let Some(m) = mate {
            let _ = write!(name, " {m}:N:0:ATCACG");
        }
        out.extend_from_slice(name.as_bytes());
        out.push(b'\n');
        out.extend_from_slice(&s);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&q);
        out.push(b'\n');
    }
    out
}

fn longread(rng: &mut Rng, n: usize, lo: usize, hi: usize, reference: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..n {
        let len = lo + rng.below(hi - lo);
        let p = rng.below(reference.len().saturating_sub(len).max(1));
        let mut s = reference[p..(p + len).min(reference.len())].to_vec();
        if rng.next().is_multiple_of(2) {
            s = revcomp(&s);
        }
        let s = mutate(rng, &s, 60);
        let q: Vec<u8> = (0..s.len()).map(|_| 33 + rng.quality(8, 30)).collect();
        out.extend_from_slice(format!("@read_{i:05} runid=fixture\n").as_bytes());
        out.extend_from_slice(&s);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&q);
        out.push(b'\n');
    }
    out
}

/// CRC-32C (Castagnoli), the same polynomial the container uses. Duplicated here
/// rather than exposed from the crate: the manifest is test scaffolding, not API.
fn crc32c(data: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &b in data {
        crc ^= u32::from(b);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0x82F6_3B78 & (!(crc & 1)).wrapping_add(1));
        }
    }
    !crc
}

fn interleave(a: &[u8], b: &[u8]) -> Vec<u8> {
    let recs = |s: &[u8]| -> Vec<Vec<u8>> {
        let lines: Vec<&[u8]> = s.split(|&c| c == b'\n').collect();
        lines
            .chunks(4)
            .filter(|c| c.len() == 4)
            .map(|c| {
                let mut r = c.join(&b'\n');
                r.push(b'\n');
                r
            })
            .collect()
    };
    let (ra, rb) = (recs(a), recs(b));
    let mut out = Vec::new();
    for (x, y) in ra.iter().zip(&rb) {
        out.extend_from_slice(x);
        out.extend_from_slice(y);
    }
    out
}

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("usage: make_fixtures <output-dir>");
        std::process::exit(2)
    });
    let dir = Path::new(&dir);
    fs::create_dir_all(dir).expect("create fixture dir");

    let mut rng = Rng(0x0F9C_0DEC_5EED_1234);
    let ref_short = genome(&mut rng, 20_000);
    let ref_long = genome(&mut rng, 6_000);
    let ref_dense = genome(&mut rng, 1_200);
    let mut manifest = String::from("# name\texpected_len\texpected_crc32c\tnote\n");

    // Record what the archive actually decodes to, not what went in, and check the
    // two agree wherever the layout promises to preserve order. The reorder-lossy
    // layout deliberately renumbers, so there the decoded stream *is* the contract
    // and the input is not — recording the input would pin a value no build ever
    // produced.
    let mut emit =
        |name: &str, archive: Vec<u8>, input: &[u8], order_preserved: bool, note: &str| {
            let mut decoded = Vec::new();
            fqxv::decompress(&archive[..], &mut decoded, 1)
                .unwrap_or_else(|e| panic!("fixture `{name}` does not round-trip: {e}"));
            if order_preserved {
                assert!(
                    decoded == input,
                    "fixture `{name}` claims to preserve order but did not round-trip"
                );
            }
            fs::write(dir.join(format!("{name}.fqxv")), &archive).expect("write fixture");
            let _ = writeln!(
                manifest,
                "{name}\t{}\t{}\t{note}",
                decoded.len(),
                crc32c(&decoded)
            );
            println!("  {name:<16} {:>8} B archive   {note}", archive.len());
        };

    // 1. Plain layout, order-k sequence codec — the default path.
    let se = illumina(&mut rng, 400, 100, &ref_short, None);
    let mut out = Vec::new();
    compress(&se[..], &mut out, Params::default()).expect("plain_se");
    emit("plain_se", out, &se, true, "plain layout, order-k sequence");

    // 2. G=2 interleaved, and the only user of the header extension region: the
    //    member-label record, non-critical tag 0x01.
    let r1 = illumina(&mut rng, 200, 100, &ref_short, Some(1));
    let r2 = illumina(&mut rng, 200, 100, &ref_short, Some(2));
    let mut out = Vec::new();
    compress_multi(
        vec![Box::new(&r1[..]), Box::new(&r2[..])],
        &mut out,
        Params {
            member_labels: vec!["R1".into(), "R2".into()],
            ..Params::default()
        },
    )
    .expect("paired_labels");
    emit(
        "paired_labels",
        out,
        &interleave(&r1, &r2),
        true,
        "G=2 interleaved, header extension tag 0x01 (member labels)",
    );

    // 3. The whole-file reorder layout — a different header and block structure
    //    entirely. Order is discarded here because that is the only configuration
    //    in which the layout wins its never-worse gate at a size worth checking in;
    //    with the permutation stored, plain is genuinely smaller and the gate
    //    rightly says so.
    let ro = illumina(&mut rng, 4_000, 100, &ref_dense, None);
    let mut out = Vec::new();
    compress(
        &ro[..],
        &mut out,
        Params {
            reorder: true,
            keep_order: false,
            // What `--order shuffle` sets. Regenerating names is most of why this
            // layout wins here: it drops the name stream entirely rather than
            // coding it, which is what tips the never-worse gate at fixture scale.
            regenerate_names: true,
            rescue: true,
            ..Params::default()
        },
    )
    .expect("reorder");
    emit(
        "reorder",
        out,
        &ro,
        false,
        "whole-file global-cluster reorder layout (order discarded)",
    );

    // 4. Long reads at real coverage, so the cross-read overlap codec wins its
    //    gate and the fixture pins a sequence method byte other than 0.
    let lr = longread(&mut rng, 120, 1_200, 2_500, &ref_long);
    let mut out = Vec::new();
    compress(&lr[..], &mut out, Params::default()).expect("longread");
    emit(
        "longread",
        out,
        &lr,
        true,
        "long-read cross-read overlap sequence codec",
    );

    fs::write(dir.join("manifest.tsv"), manifest).expect("write manifest");
    println!("wrote manifest.tsv");
}
