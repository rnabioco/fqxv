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


def _build_fixture(dest_dir):
    """Compress a tiny synthetic FASTQ into a `.fqxv` the tests can read."""
    fastq = dest_dir / "tiny.fastq"
    lines = []
    for i in range(8):
        lines += [f"@read{i} desc", "ACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIII"]
    fastq.write_text("\n".join(lines) + "\n")
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
