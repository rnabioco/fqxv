"""End-to-end smoke tests for the fqxv Python bindings.

Run after `maturin develop`:  pytest crates/fqxv-python/tests/
Requires a `.fqxv` fixture; pass one via FQXV_FIXTURE or the tests build one
from the repo's checked-in sample if the `fqxv` CLI is on PATH.
"""
import os
import subprocess
import shutil
import pathlib

import pytest

import fqxv


def _fixture(tmp_path):
    env = os.environ.get("FQXV_FIXTURE")
    if env and pathlib.Path(env).exists():
        return env
    # Fall back to a checked-in archive at the repo root.
    root = pathlib.Path(__file__).resolve().parents[3]
    for name in ("SRR2584863.fqxv", "SRR28588231.fqxv"):
        cand = root / name
        if cand.exists():
            return str(cand)
    pytest.skip("no .fqxv fixture available (set FQXV_FIXTURE)")


def test_iterate_records(tmp_path):
    path = _fixture(tmp_path)
    recs = list(fqxv.open(path))
    assert recs, "expected at least one record"
    r0 = recs[0]
    assert isinstance(r0.name, bytes)
    assert isinstance(r0.sequence, bytes)
    assert isinstance(r0.quality, bytes)
    assert len(r0.sequence) == len(r0.quality) == len(r0)


def test_bytes_input_matches_path(tmp_path):
    path = _fixture(tmp_path)
    from_path = list(fqxv.open(path))
    data = pathlib.Path(path).read_bytes()
    from_bytes = list(fqxv.open(data))
    assert len(from_path) == len(from_bytes)
    assert from_path[0].sequence == from_bytes[0].sequence


def test_break_does_not_hang(tmp_path):
    path = _fixture(tmp_path)
    it = fqxv.open(path)
    first = next(it)
    assert first is not None
    del it  # early drop must not deadlock the decode thread


def test_inspect_matches_iteration(tmp_path):
    path = _fixture(tmp_path)
    info = fqxv.inspect(path)
    n = sum(1 for _ in fqxv.open(path))
    assert info.reads == n
    assert info.format_version >= 3


def test_decompress_to_bytes_roundtrips(tmp_path):
    path = _fixture(tmp_path)
    raw = fqxv.decompress_to_bytes(path)
    # Four lines per record.
    assert raw.count(b"\n") == 4 * sum(1 for _ in fqxv.open(path))


def test_projection_matches_iteration(tmp_path):
    path = _fixture(tmp_path)
    info = fqxv.inspect(path)
    if info.reordered:
        pytest.skip("projection unavailable for reordered archives")
    seqs = fqxv.read_sequences(path)
    recs = list(fqxv.open(path))
    assert len(seqs) == len(recs)
    assert seqs[0] == recs[0].sequence
    idx = fqxv.open_index(path)
    assert idx.total_reads == len(recs)
    block0 = fqxv.read_block(path, 0)
    assert block0[0].sequence == recs[0].sequence
