//! End-to-end CLI test for opt-in lossy quality binning: `compress
//! --quality-bin bin4` must change qualities to the RTA4 4-bin table, keep
//! names/bases exact, and be reported by `info`.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// The fqxv binary under test (set by Cargo for integration tests).
const FQXV: &str = env!("CARGO_BIN_EXE_fqxv");

/// A private temp dir Cargo manages for this test binary.
fn tmp(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    p.push(name);
    p
}

/// NovaSeq X / RTA4 v1.2 4-bin table, replicated here so the test pins the
/// expected mapping independently of the library.
fn bin4(byte: u8) -> u8 {
    let q = byte.saturating_sub(33);
    let b = match q {
        0..=2 => 2,
        3..=17 => 12,
        18..=29 => 24,
        _ => 40,
    };
    33 + b
}

fn run(args: &[&str]) -> std::process::Output {
    let out = Command::new(FQXV).args(args).output().expect("spawn fqxv");
    assert!(
        out.status.success(),
        "fqxv {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

/// Every 4th FASTQ line starting at `which` (1 = sequence, 3 = quality),
/// concatenated in order.
fn record_line(fastq: &[u8], which: usize) -> Vec<u8> {
    fastq
        .split(|&b| b == b'\n')
        .enumerate()
        .filter(|(i, l)| i % 4 == which && !l.is_empty())
        .flat_map(|(_, l)| l.iter().copied())
        .collect()
}

#[test]
fn compress_quality_bin4_roundtrip() {
    // Qualities span every RTA4 band so the mapping is exercised end to end.
    let input = b"\
@read.1 lane\n\
ACGTAC\n\
+\n\
!5?JI#\n\
@read.2 lane\n\
NNGGCC\n\
+\n\
0:F,+.\n";

    let in_path = tmp("qb_in.fastq");
    let arc_path = tmp("qb.fqxv");
    let rt_path = tmp("qb_rt.fastq");
    fs::write(&in_path, input).unwrap();

    // Compress with 4-bin lossy quality.
    run(&[
        "compress",
        in_path.to_str().unwrap(),
        "-o",
        arc_path.to_str().unwrap(),
        "--force",
        "--quality-bin",
        "bin4",
        "--threads",
        "1",
    ]);

    // `info --tsv` reports the binning tag (column 6 = quality_binning = 2).
    let info = run(&["info", arc_path.to_str().unwrap(), "--tsv"]);
    let stdout = String::from_utf8(info.stdout).unwrap();
    let mut lines = stdout.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let data: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = header.iter().position(|&h| h == "quality_binning").unwrap();
    assert_eq!(data[col], "2", "info should report the 4-bin tag");

    // Decompress and check the lossy contract: bases exact, qualities binned.
    run(&[
        "decompress",
        arc_path.to_str().unwrap(),
        "-o",
        rt_path.to_str().unwrap(),
        "--force",
        "--threads",
        "1",
    ]);
    let rt = fs::read(&rt_path).unwrap();

    assert_eq!(
        record_line(&rt, 1),
        record_line(input, 1),
        "bases must survive lossy quality binning exactly"
    );
    let want: Vec<u8> = record_line(input, 3).iter().map(|&b| bin4(b)).collect();
    assert_eq!(
        record_line(&rt, 3),
        want,
        "decompressed qualities must equal input through the RTA4 4-bin table"
    );
}
