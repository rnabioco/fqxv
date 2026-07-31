//! Golden-archive compatibility: every checked-in fixture must still decode to
//! the bytes it decoded to when it was written.
//!
//! This is the only test in the workspace that compares this build against a
//! *previous* one. Every other test — the round-trips, the proptests, the fuzz
//! targets, the thread-determinism byte-comparisons, the accession corpus —
//! compares a build against itself, and is therefore invariant to any change that
//! moves the encoder and decoder together. That is precisely the change that
//! breaks archives already on disk, so without these fixtures the project's
//! central promise is asserted and never checked.
//!
//! A failure here is not a flaky test. It means this build cannot correctly read
//! an archive an earlier release wrote. Regenerating the fixtures makes the red go
//! away and ships the regression; see `examples/make_fixtures.rs`.

use std::fs;
use std::path::PathBuf;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// CRC-32C (Castagnoli) — mirrors `examples/make_fixtures.rs`, which produced the
/// manifest values.
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

struct Fixture {
    name: String,
    expected_len: usize,
    expected_crc: u32,
    note: String,
}

fn manifest() -> Vec<Fixture> {
    let path = fixture_dir().join("manifest.tsv");
    let raw =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    raw.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 4, "malformed manifest row: {l:?}");
            Fixture {
                name: f[0].to_string(),
                expected_len: f[1].parse().expect("expected_len"),
                expected_crc: f[2].parse().expect("expected_crc32c"),
                note: f[3].to_string(),
            }
        })
        .collect()
}

#[test]
fn golden_archives_still_decode() {
    let fixtures = manifest();
    assert!(
        !fixtures.is_empty(),
        "no golden fixtures found — the compatibility guarantee would be untested"
    );

    for f in &fixtures {
        let path = fixture_dir().join(format!("{}.fqxv", f.name));
        let archive = fs::read(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));

        let mut got = Vec::new();
        fqxv::decompress(&archive[..], &mut got, 1).unwrap_or_else(|e| {
            panic!(
                "BACKWARD COMPATIBILITY BROKEN: fixture `{}` ({}) no longer decodes: {e}\n\
                 This archive was written by an earlier release and must stay readable.",
                f.name, f.note
            )
        });

        // Length first: it localizes a truncation or an extra/missing record far
        // better than a checksum mismatch can.
        assert_eq!(
            got.len(),
            f.expected_len,
            "BACKWARD COMPATIBILITY BROKEN: fixture `{}` ({}) decoded to {} bytes, expected {}",
            f.name,
            f.note,
            got.len(),
            f.expected_len
        );
        assert_eq!(
            crc32c(&got),
            f.expected_crc,
            "BACKWARD COMPATIBILITY BROKEN: fixture `{}` ({}) decoded to the right length but \
             different bytes — a codec now decodes an old archive incorrectly",
            f.name,
            f.note
        );
    }
}

/// The fixtures are only worth their bytes if they cover distinct paths. This
/// pins the set so a future edit cannot quietly drop one — losing, say, the only
/// reorder-layout fixture would leave that layout with no compatibility coverage
/// while this file still looked green.
#[test]
fn golden_fixtures_cover_every_layout() {
    let names: Vec<String> = manifest().into_iter().map(|f| f.name).collect();
    for required in [
        "plain_se",      // plain layout, order-k sequence
        "paired_labels", // grouped, and the header extension region
        "reorder",       // the whole-file reorder layout
        "longread",      // the cross-read overlap sequence codec
    ] {
        assert!(
            names.iter().any(|n| n == required),
            "golden fixture `{required}` is missing; every on-disk layout needs one"
        );
    }
}

/// `inspect` reads the header and footer without decoding any block, so it has its
/// own parsing path — one that random access and the remote/range clients share.
/// A fixture that decodes but no longer inspects would be a real regression for
/// those callers.
#[test]
fn golden_archives_still_inspect() {
    for f in manifest() {
        let path = fixture_dir().join(format!("{}.fqxv", f.name));
        let archive = fs::read(&path).expect("read fixture");
        let info = fqxv::inspect(std::io::Cursor::new(&archive)).unwrap_or_else(|e| {
            panic!("fixture `{}` ({}) no longer inspects: {e}", f.name, f.note)
        });
        assert_eq!(
            info.format_version >> 8,
            u16::from(fqxv::FORMAT_MAJOR),
            "fixture `{}` reports a foreign format major",
            f.name
        );
        assert!(
            info.reads > 0,
            "fixture `{}` inspected as an empty archive",
            f.name
        );
    }
}
