"""Stream-selective decode (`open(streams=...)`) and FASTA output (`fasta=True`).

The selection contract under test: selected fields decode byte-identically to a
full decode, deselected fields come back as empty ``bytes`` (their coded streams
are skipped, never entropy-decoded), and ``fasta=True`` matches the CLI's
``decompress --fasta`` byte-for-byte. Long-read archives exercise the one
cross-stream dependency (quality is coded against the bases), and a globally
reordered archive exercises the whole-file layout's selection paths.

Run after ``maturin develop``:  pytest crates/fqxv-python/tests/test_select.py
"""
import itertools
import os
import pathlib
import random
import shutil
import subprocess

import pytest

import fqxv

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]


def _fqxv_cli():
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


def _compress(fastq, archive, *args):
    cli = _fqxv_cli()
    if cli is None:
        pytest.skip("fqxv CLI not found (set FQXV_BIN or build the workspace)")
    subprocess.run(
        [cli, "--quiet", "compress", *args, str(fastq), "-o", str(archive)],
        check=True,
        capture_output=True,
    )
    return str(archive)


@pytest.fixture(scope="session")
def archive(tmp_path_factory):
    """A small plain-layout archive spanning several blocks, with per-read
    distinct names, varied lengths, and non-constant quality — so a stream
    swapped or dropped by selection cannot silently equal another read's."""
    d = tmp_path_factory.mktemp("fqxv_select")
    fastq = d / "reads.fastq"
    rng = random.Random(42)
    with fastq.open("w") as fh:
        for i in range(64):
            seq = "".join(rng.choices("ACGT", k=20 + i % 21))
            qual = "".join(chr(33 + ((j * 3 + i) % 40)) for j in range(len(seq)))
            fh.write(f"@sel{i} x={i}\n{seq}\n+\n{qual}\n")
    return _compress(fastq, d / "reads.fqxv", "--block-reads", "16")


@pytest.fixture(scope="session")
def full(archive):
    """The full-decode reference all selections are compared against."""
    return list(fqxv.open(archive))


# Every subset of the three streams, including empty (decode nothing).
SELECTIONS = [
    sel
    for n in range(4)
    for sel in itertools.combinations(("names", "sequence", "quality"), n)
]


@pytest.mark.parametrize("selected", SELECTIONS, ids=lambda s: "+".join(s) or "none")
def test_selection_matches_full_decode(archive, full, selected):
    recs = list(fqxv.open(archive, streams=selected))
    assert len(recs) == len(full)
    for got, want in zip(recs, full):
        assert got.name == (want.name if "names" in selected else b"")
        assert got.sequence == (want.sequence if "sequence" in selected else b"")
        assert got.quality == (want.quality if "quality" in selected else b"")


def test_single_name_and_aliases(archive, full):
    seqs = [r.sequence for r in full]
    # A bare string is one stream name, not an iterable of characters.
    assert [r.sequence for r in fqxv.open(archive, streams="sequence")] == seqs
    assert [r.sequence for r in fqxv.open(archive, streams="seq")] == seqs
    # The short aliases mirror Index.stream_range's accepted names.
    recs = list(fqxv.open(archive, streams=("name", "qual")))
    assert [r.name for r in recs] == [r.name for r in full]
    assert [r.quality for r in recs] == [r.quality for r in full]
    assert all(r.sequence == b"" for r in recs)


def test_empty_selection_counts_reads(archive, full):
    """streams=() decodes nothing: every field empty, one record per read."""
    recs = list(fqxv.open(archive, streams=()))
    assert len(recs) == len(full)
    assert all(r.name == r.sequence == r.quality == b"" for r in recs)
    assert all(len(r) == 0 for r in recs)  # len() is the sequence length


def test_selection_works_on_bytes_source(archive, full):
    data = pathlib.Path(archive).read_bytes()
    got = [r.sequence for r in fqxv.open(data, streams=("seq",))]
    assert got == [r.sequence for r in full]


def test_bad_selection_raises(archive):
    with pytest.raises(ValueError, match="unknown stream"):
        fqxv.open(archive, streams="nonsense")
    with pytest.raises(ValueError, match="unknown stream"):
        fqxv.open(archive, streams=("names", "nonsense"))
    with pytest.raises(TypeError):
        fqxv.open(archive, streams=123)
    with pytest.raises(TypeError):
        fqxv.open(archive, streams=(1, 2))


# --------------------------------------------------------------------------- #
# FASTA output (fasta=True), mirroring the CLI's `decompress --fasta`.
# --------------------------------------------------------------------------- #
def _fasta_of(records):
    return b"".join(b">" + r.name + b"\n" + r.sequence + b"\n" for r in records)


def test_fasta_bytes_matches_full_decode(archive, full):
    assert fqxv.decompress_to_bytes(archive, fasta=True) == _fasta_of(full)


def test_fasta_path_matches_cli(archive, full, tmp_path):
    cli = _fqxv_cli()
    if cli is None:
        pytest.skip("fqxv CLI not found (set FQXV_BIN or build the workspace)")
    py_out = tmp_path / "py.fasta"
    n = fqxv.decompress_to_path(archive, py_out, fasta=True)
    assert n == len(full)
    cli_out = tmp_path / "cli.fasta"
    subprocess.run(
        [cli, "--quiet", "decompress", "--fasta", archive, "-o", str(cli_out)],
        check=True,
        capture_output=True,
    )
    assert py_out.read_bytes() == cli_out.read_bytes()


# --------------------------------------------------------------------------- #
# Long-read archives: the one cross-stream dependency. Quality is coded against
# the bases (fqzcomp MODE_SEQ), so selecting quality decodes the sequence stream
# internally — but still returns it empty unless selected.
# --------------------------------------------------------------------------- #
_QMAP = {"A": "I", "C": "F", "G": ";", "T": "5"}


@pytest.fixture(scope="session")
def longread_archive(tmp_path_factory):
    d = tmp_path_factory.mktemp("fqxv_select_lr")
    fastq = d / "lr.fastq"
    rng = random.Random(7)
    with fastq.open("w") as fh:
        for i in range(300):
            seq = "".join(rng.choices("ACGT", k=750 + (i % 120)))
            qual = "".join(_QMAP[b] for b in seq)
            fh.write(f"@ln{i:08x}-1111-2222 ch={i}\n{seq}\n+\n{qual}\n")
    return _compress(
        fastq, d / "lr.fqxv", "--platform", "nanopore", "--block-reads", "100"
    )


def test_longread_quality_selection_decodes_sequence_internally(longread_archive):
    # Sanity: the fixture must really condition quality on sequence, or this
    # test silently exercises the short-read path.
    idx = fqxv.open_index(longread_archive)
    raw = pathlib.Path(longread_archive).read_bytes()
    start, end = idx.stream_range(0, "quality")
    assert fqxv.quality_needs_sequence_bytes(raw[start:end])

    full = list(fqxv.open(longread_archive))
    recs = list(fqxv.open(longread_archive, streams=("quality",)))
    assert [r.quality for r in recs] == [r.quality for r in full]
    assert all(r.name == b"" and r.sequence == b"" for r in recs)

    # Sequence-only is the fast path the selection exists for.
    srecs = list(fqxv.open(longread_archive, streams="seq"))
    assert [r.sequence for r in srecs] == [r.sequence for r in full]
    assert all(r.name == b"" and r.quality == b"" for r in srecs)

    assert fqxv.decompress_to_bytes(longread_archive, fasta=True) == _fasta_of(full)


# --------------------------------------------------------------------------- #
# Globally reordered archive (--order shuffle): the whole-file layout's
# selection paths (independent per-stream frames + flip bitmap).
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="session")
def shuffle_archive(tmp_path_factory):
    """Overlapping reads sampled from one reference, half reverse-complemented,
    so the reorder layout genuinely wins the size contest and is kept."""
    d = tmp_path_factory.mktemp("fqxv_select_shuffle")
    fastq = d / "overlap.fastq"
    rng = random.Random(20260810)
    ref = "".join(rng.choice("ACGT") for _ in range(50_000))
    comp = str.maketrans("ACGT", "TGCA")
    with fastq.open("w") as fh:
        for i in range(4000):
            start = rng.randrange(len(ref) - 150)
            seq = ref[start : start + 150]
            if i % 2:
                seq = seq.translate(comp)[::-1]
            qual = "".join(chr(35 + ((j + i) % 30)) for j in range(len(seq)))
            fh.write(f"@ov{i}\n{seq}\n+\n{qual}\n")
    return _compress(fastq, d / "overlap.fqxv", "--order", "shuffle")


def test_shuffle_selection_matches_full_decode(shuffle_archive):
    # The fixture is built to make the reorder layout win; if a codec change
    # ever re-routes it to plain, this assert says so rather than silently
    # testing the wrong layout.
    assert fqxv.inspect(shuffle_archive).reordered
    full = list(fqxv.open(shuffle_archive))
    assert len(full) == 4000
    for selected in [("names",), ("sequence",), ("quality",), ("names", "sequence")]:
        recs = list(fqxv.open(shuffle_archive, streams=selected))
        assert len(recs) == len(full)
        for got, want in zip(recs, full):
            assert got.name == (want.name if "names" in selected else b"")
            assert got.sequence == (want.sequence if "sequence" in selected else b"")
            assert got.quality == (want.quality if "quality" in selected else b"")
    assert fqxv.decompress_to_bytes(shuffle_archive, fasta=True) == _fasta_of(full)
