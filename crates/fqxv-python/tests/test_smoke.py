"""End-to-end smoke tests for the fqxv Python bindings.

Run after `maturin develop`:  pytest crates/fqxv-python/tests/

The tests build a tiny `.fqxv` fixture with the `fqxv` CLI, so they finish in
well under a second and never latch onto whatever large archive happens to sit
at the repo root. Overrides:
  FQXV_FIXTURE  path to an existing archive to test instead of the built one
  FQXV_BIN      path to the `fqxv` CLI (else PATH, else target/{release,debug})
"""
import os
import pathlib
import shutil
import subprocess

import pytest

import fqxv

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]


def _fqxv_cli():
    """Locate the fqxv CLI binary: FQXV_BIN, then PATH, then a built target/."""
    env = os.environ.get("FQXV_BIN")
    if env:
        return env
    found = shutil.which("fqxv")
    if found:
        return found
    for profile in ("release", "debug"):
        cand = REPO_ROOT / "target" / profile / "fqxv"
        if cand.exists():
            return str(cand)
    return None


def _write_fastq(dest_dir):
    """Write a tiny synthetic FASTQ and return its path."""
    fastq = dest_dir / "tiny.fastq"
    lines = []
    for i in range(8):
        lines += [f"@read{i} desc", "ACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIII"]
    fastq.write_text("\n".join(lines) + "\n")
    return fastq


def _build_fixture(dest_dir):
    """Compress a tiny synthetic FASTQ into a `.fqxv` the tests can read."""
    fastq = _write_fastq(dest_dir)
    cli = _fqxv_cli()
    if cli is None:
        pytest.skip("fqxv CLI not found (set FQXV_BIN or build the workspace)")
    archive = dest_dir / "tiny.fqxv"
    subprocess.run(
        [cli, "--quiet", "compress", str(fastq), "-o", str(archive)],
        check=True,
        capture_output=True,
    )
    return archive


@pytest.fixture(scope="session")
def archive(tmp_path_factory):
    """A tiny `.fqxv` fixture, built once per test session."""
    env = os.environ.get("FQXV_FIXTURE")
    if env and pathlib.Path(env).exists():
        return env
    return str(_build_fixture(tmp_path_factory.mktemp("fqxv")))


@pytest.fixture(scope="session")
def fastq(tmp_path_factory):
    """A tiny raw FASTQ fixture. `estimate` reads FASTQ input, not an archive."""
    return str(_write_fastq(tmp_path_factory.mktemp("fqxv_fastq")))


def test_iterate_records(archive):
    recs = list(fqxv.open(archive))
    assert recs, "expected at least one record"
    r0 = recs[0]
    assert isinstance(r0.name, bytes)
    assert isinstance(r0.sequence, bytes)
    assert isinstance(r0.quality, bytes)
    assert len(r0.sequence) == len(r0.quality) == len(r0)


def test_bytes_input_matches_path(archive):
    from_path = list(fqxv.open(archive))
    data = pathlib.Path(archive).read_bytes()
    from_bytes = list(fqxv.open(data))
    assert len(from_path) == len(from_bytes)
    assert from_path[0].sequence == from_bytes[0].sequence


def test_break_does_not_hang(archive):
    it = fqxv.open(archive)
    first = next(it)
    assert first is not None
    del it  # early drop must not deadlock the decode thread


def test_inspect_matches_iteration(archive):
    info = fqxv.inspect(archive)
    n = sum(1 for _ in fqxv.open(archive))
    assert info.reads == n
    # `format_version` is packed `(major << 8) | minor`; the major must be v1.
    assert info.format_version >> 8 == 1


def test_decompress_to_bytes_roundtrips(archive):
    raw = fqxv.decompress_to_bytes(archive)
    # Four lines per record.
    assert raw.count(b"\n") == 4 * sum(1 for _ in fqxv.open(archive))


def test_projection_matches_iteration(archive):
    info = fqxv.inspect(archive)
    if info.reordered:
        pytest.skip("projection unavailable for reordered archives")
    seqs = fqxv.read_sequences(archive)
    recs = list(fqxv.open(archive))
    assert len(seqs) == len(recs)
    assert seqs[0] == recs[0].sequence
    idx = fqxv.open_index(archive)
    assert idx.total_reads == len(recs)
    block0 = fqxv.read_block(archive, 0)
    assert block0[0].sequence == recs[0].sequence


def test_estimate_projects_size(fastq):
    est = fqxv.estimate(fastq)
    # The whole 8-read fixture fits in one sample, so the numbers are exact.
    assert est.sample_reads == 8
    assert est.exhausted
    assert est.archive_bytes > 0
    assert est.ratio > 0.0
    # archive_bytes is the three streams plus per-block frame overhead.
    assert est.names_bytes + est.sequence_bytes + est.quality_bytes <= est.archive_bytes


def test_estimate_accepts_bytes_and_is_deterministic(fastq):
    from_path = fqxv.estimate(fastq)
    from_bytes = fqxv.estimate(pathlib.Path(fastq).read_bytes())
    assert from_bytes.sample_reads == from_path.sample_reads
    assert from_bytes.archive_bytes == from_path.archive_bytes


def test_estimate_rejects_unknown_binning(fastq):
    with pytest.raises(ValueError):
        fqxv.estimate(fastq, quality_binning="nonsense")


def test_estimate_accepts_paired_sources(fastq):
    """A list/tuple of mates estimates one archive: the per-stream sample sizes
    sum, so a pair reports twice a single mate's reads and bytes."""
    single = fqxv.estimate(fastq)
    pair = fqxv.estimate([fastq, fastq])
    assert pair.sample_reads == 2 * single.sample_reads
    assert pair.sample_bases == 2 * single.sample_bases
    assert pair.archive_bytes == 2 * single.archive_bytes
    assert pair.exhausted
    # A tuple is treated identically to a list.
    assert fqxv.estimate((fastq, fastq)).archive_bytes == pair.archive_bytes
    # bytes is a single source, not an iterable of sources.
    data = pathlib.Path(fastq).read_bytes()
    assert fqxv.estimate(data).sample_reads == single.sample_reads


def test_estimate_rejects_empty_source_list():
    with pytest.raises(ValueError):
        fqxv.estimate([])


@pytest.fixture(scope="session")
def illumina_fastq(tmp_path_factory):
    """Illumina-style names (colon-delimited instrument coordinates)."""
    d = tmp_path_factory.mktemp("il")
    p = d / "il.fastq"
    p.write_text(
        "".join(
            f"@INST:1:FC:1:{i}:100:200 1:N:0:1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
            for i in range(400)
        )
    )
    return str(p)


@pytest.fixture(scope="session")
def nanopore_fastq(tmp_path_factory):
    """Nanopore-style names (UUID read id, runid= description tag)."""
    d = tmp_path_factory.mktemp("nano")
    p = d / "nano.fastq"

    def rec(i):
        seq = ("ACGTTGCAAGTC" * 30)[: 300 + i % 40]
        qual = "".join(chr(33 + ((j * 7 + i) % 40)) for j in range(len(seq)))
        uuid = f"{i:08x}-1234-5678-9abc-def012345678"
        return f"@{uuid} runid=abc ch={i}\n{seq}\n+\n{qual}\n"

    p.write_text("".join(rec(i) for i in range(400)))
    return str(p)


def test_estimate_reports_platform(illumina_fastq, nanopore_fastq):
    assert fqxv.estimate(illumina_fastq).platform == "illumina"
    assert fqxv.estimate(nanopore_fastq).platform == "nanopore"
    # A matched pair keeps the shared platform.
    assert fqxv.estimate([illumina_fastq, illumina_fastq]).platform == "illumina"


def test_estimate_rejects_cross_platform_group(illumina_fastq, nanopore_fastq):
    """An accidental Illumina + Nanopore group is an error, not a silent sum."""
    with pytest.raises(fqxv.FqxvError, match="multiple platforms"):
        fqxv.estimate([illumina_fastq, nanopore_fastq])


def test_verify_passes_on_good_archive(archive):
    # Returns None on success; raises on failure.
    assert fqxv.verify(archive) is None


def test_verify_rejects_non_archive():
    with pytest.raises(Exception):
        fqxv.verify(b"definitely not an fqxv archive")


def test_version_is_exposed():
    # The package re-exports its installed version so callers can log/branch on it.
    assert isinstance(fqxv.__version__, str)
    assert fqxv.__version__  # non-empty


# --------------------------------------------------------------------------- #
# Long-read (Nanopore/PacBio) fixtures. Reads averaging >500 bp make the quality
# coder condition each score on its base (fqzcomp MODE_SEQ), so the quality stream
# can't be projected alone — the decoder must be handed that group's decoded
# sequence. Short-read fixtures never hit this path, so it needs its own coverage.
# --------------------------------------------------------------------------- #
# Quality that tracks the base makes the sequence context win the size trial, so
# MODE_SEQ (quality-needs-sequence) is reliably selected.
_QMAP = {"A": "I", "C": "F", "G": ";", "T": "5"}


def _write_longread_fastq(dest_dir):
    fastq = dest_dir / "longread.fastq"
    import random

    rng = random.Random(7)
    with fastq.open("w") as fh:
        for i in range(300):
            length = 750 + (i % 120)  # mean ~810 bp, comfortably over 500
            seq = "".join(rng.choices("ACGT", k=length))
            qual = "".join(_QMAP[b] for b in seq)
            fh.write(
                f"@ln{i:08x}-1111-2222-3333-444455556666 runid=x ch={i}\n"
                f"{seq}\n+\n{qual}\n"
            )
    return fastq


@pytest.fixture(scope="session")
def longread_fastq(tmp_path_factory):
    return str(_write_longread_fastq(tmp_path_factory.mktemp("fqxv_lr_fastq")))


@pytest.fixture(scope="session")
def longread_archive(tmp_path_factory):
    """A multi-group Nanopore archive whose quality is coded against the sequence.
    Explicit ``--platform nanopore`` keeps each block's sequence self-contained (no
    shared reference), so a projected quality read can refetch and decode it."""
    d = tmp_path_factory.mktemp("fqxv_lr")
    fastq = _write_longread_fastq(d)
    cli = _fqxv_cli()
    if cli is None:
        pytest.skip("fqxv CLI not found (set FQXV_BIN or build the workspace)")
    archive = d / "longread.fqxv"
    # --block-reads 100 → 3 row groups over 300 reads.
    subprocess.run(
        [
            cli, "--quiet", "compress", "--platform", "nanopore",
            "--block-reads", "100", str(fastq), "-o", str(archive),
        ],
        check=True,
        capture_output=True,
    )
    return str(archive)


def test_longread_quality_is_coded_against_sequence(longread_archive):
    info = fqxv.inspect(longread_archive)
    assert not info.reordered, "long-read fixture must be projectable (plain layout)"

    idx = fqxv.open_index(longread_archive)
    assert idx.num_groups >= 2, "fixture should span multiple row groups"

    # The quality column of group 0 must be sequence-conditioned — otherwise this
    # test is silently exercising the short-read path, not the one it exists for.
    raw = pathlib.Path(longread_archive).read_bytes()
    start, end = idx.stream_range(0, "quality")
    assert fqxv.quality_needs_sequence_bytes(raw[start:end])

    # The projection paths must decode it correctly: read_qualities internally
    # refetches each group's sequence to decode its quality.
    local = list(fqxv.open(longread_archive))
    assert fqxv.read_qualities(longread_archive) == [r.quality for r in local]
    assert fqxv.read_sequences(longread_archive) == [r.sequence for r in local]
    block0 = fqxv.read_block(longread_archive, 0)
    assert [r.quality for r in block0] == [r.quality for r in local[: len(block0)]]


@pytest.mark.parametrize("fixture_name", ["fastq", "longread_fastq"])
@pytest.mark.parametrize("level", [1, 5, 7, 8, 9])
def test_estimate_archive_bytes_matches_cli(request, fixture_name, level):
    """`estimate` re-implements the CLI's ``level_to_*`` effort mapping; a copy that
    drifts would silently mis-size. Lock them together: the projected archive size
    must equal ``fqxv compress --level N --estimate tsv`` on the same input. The
    fixtures are small enough that both fully consume them, so the size is exact.
    The long-read fixture additionally exercises the tiler knobs (`level_to_tile`),
    which short reads never touch."""
    cli = _fqxv_cli()
    if cli is None:
        pytest.skip("fqxv CLI not found (set FQXV_BIN or build the workspace)")
    fastq = request.getfixturevalue(fixture_name)
    out = subprocess.run(
        [cli, "compress", "--level", str(level), "--estimate", "tsv", fastq],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip().splitlines()
    # TSV: "file\tinput_bytes\test_fqxv_bytes\tratio" then one data row.
    cli_bytes = int(out[1].split("\t")[2])
    py_bytes = fqxv.estimate(fastq, level=level).archive_bytes
    assert py_bytes == cli_bytes, (
        f"level {level} on {fixture_name}: python={py_bytes} cli={cli_bytes} — "
        "the level_to_* mapping in fqxv-python drifted from fqxv-cli"
    )
