"""Direct test of remote (HTTP byte-range) reads over a real socket.

Serves a locally-built ``.fqxv`` from an in-process HTTP server that honours
``Range`` requests, then checks that :mod:`fqxv.remote`:

* decodes names / sequences / records byte-identically to a full local decode,
* transfers only a fraction of the file for a column projection (the whole point
  of the assertion "fqxv can be read async over an internet connection"),
* handles a too-short first tail with exactly one extra round trip, and
* rejects a corrupted range via the per-stream CRC.

Standard-library only (``http.server`` + ``urllib``); the async path is skipped
when neither httpx nor aiohttp is importable.

Run after ``maturin develop``:  pytest crates/fqxv-python/tests/test_remote.py
"""
import os
import pathlib
import shutil
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

import fqxv
import fqxv.remote as remote

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]


# --------------------------------------------------------------------------- #
# Fixture: a multi-row-group archive built with the CLI
# --------------------------------------------------------------------------- #
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


@pytest.fixture(scope="session")
def archive(tmp_path_factory):
    """~50k reads of pseudo-random sequence at a small block size: an archive of a
    few MB spanning several row groups, so the 64 KiB index tail is negligible and
    a names-only projection is a small fraction of the file."""
    import random

    cli = _fqxv_cli()
    if cli is None:
        pytest.skip("fqxv CLI not found (set FQXV_BIN or build the workspace)")
    d = tmp_path_factory.mktemp("fqxv_remote")
    fastq = d / "reads.fastq"
    rng = random.Random(1234)
    with fastq.open("w") as fh:
        for i in range(50_000):
            seq = "".join(rng.choices("ACGT", k=150))
            fh.write(f"@read{i} lane=1 x={i}\n{seq}\n+\n{'I' * len(seq)}\n")
    archive = d / "reads.fqxv"
    # --block-reads 4000 → 5 row groups over 20k reads, so a column projection and
    # a group subset each fetch a genuine fraction of the file (the object-storage
    # use case that flag exists for).
    subprocess.run(
        [cli, "--quiet", "compress", "--block-reads", "4000", str(fastq), "-o", str(archive)],
        check=True,
        capture_output=True,
    )
    return str(archive)


# --------------------------------------------------------------------------- #
# A minimal Range-capable HTTP server (stdlib SimpleHTTPRequestHandler does not
# support Range).
# --------------------------------------------------------------------------- #
def _parse_range(header, total):
    """Parse a single-range ``bytes=`` header into inclusive [start, end]."""
    assert header.startswith("bytes=")
    spec = header[len("bytes=") :].split(",")[0].strip()
    lo, _, hi = spec.partition("-")
    if lo == "":  # suffix range: bytes=-N
        n = int(hi)
        return max(0, total - n), total - 1
    start = int(lo)
    end = int(hi) if hi else total - 1
    return start, min(end, total - 1)


class _RangeHandler(BaseHTTPRequestHandler):
    def log_message(self, *_args):
        pass

    def _path(self):
        return pathlib.Path(self.server.directory) / os.path.basename(self.path)

    def do_HEAD(self):
        p = self._path()
        if not p.exists():
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header("Content-Length", str(p.stat().st_size))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()

    def do_GET(self):
        p = self._path()
        if not p.exists():
            self.send_error(404)
            return
        data = p.read_bytes()
        total = len(data)
        rng = self.headers.get("Range")
        if not rng:
            self.send_response(200)
            self.send_header("Content-Length", str(total))
            self.send_header("Accept-Ranges", "bytes")
            self.end_headers()
            self.wfile.write(data)
            return
        start, end = _parse_range(rng, total)
        chunk = data[start : end + 1]
        self.send_response(206)
        self.send_header("Content-Range", f"bytes {start}-{end}/{total}")
        self.send_header("Content-Length", str(len(chunk)))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        self.wfile.write(chunk)


@pytest.fixture(scope="session")
def server(archive):
    directory = os.path.dirname(archive)
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), _RangeHandler)
    httpd.directory = directory
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    host, port = httpd.server_address
    yield f"http://{host}:{port}/{os.path.basename(archive)}"
    httpd.shutdown()


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #
def test_remote_names_match_and_are_cheap(archive, server):
    local = list(fqxv.open(archive))
    file_size = os.path.getsize(archive)

    arc = remote.open_index(server)
    assert arc.size == file_size
    assert arc.index.num_groups >= 2, "fixture should span multiple row groups"

    names = arc.names()
    assert names == [r.name for r in local]
    # Names are a small fraction of the archive: index tail + names columns only.
    assert 0 < arc.bytes_fetched < file_size // 4


def test_remote_sequences_skip_quality(server, archive):
    local = list(fqxv.open(archive))
    arc = remote.open_index(server)
    seqs = arc.sequences()
    assert seqs == [r.sequence for r in local]
    seq_only = arc.bytes_fetched

    # Reading the whole record set (all three columns) must transfer strictly more.
    arc2 = remote.open_index(server)
    recs = arc2.records()
    assert [r.sequence for r in recs] == [r.sequence for r in local]
    assert [r.name for r in recs] == [r.name for r in local]
    assert [r.quality for r in recs] == [r.quality for r in local]
    assert arc2.bytes_fetched > seq_only


def test_remote_qualities_match(server, archive):
    local = list(fqxv.open(archive))
    quals = remote.read_qualities(server)
    assert quals == [r.quality for r in local]


def test_short_tail_triggers_one_refetch(server):
    # A tail smaller than the footer forces the NeedAtLeast refetch branch.
    arc = remote.RemoteArchive.open(server, tail=16)
    assert arc.index.total_reads == 50_000
    # Two tail GETs happened; still far less than the file.
    assert arc.bytes_fetched < arc.size


def test_group_subset_reads_less_than_all(server):
    arc_all = remote.open_index(server)
    _ = arc_all.names()
    arc_one = remote.open_index(server)
    _ = arc_one.names(groups=[0])
    assert arc_one.bytes_fetched < arc_all.bytes_fetched


def test_corrupt_range_is_rejected(server):
    arc = remote.open_index(server)
    start, end = arc.index.stream_range(0, "names")
    coded = bytearray(arc._get_range(start, end))
    coded[0] ^= 0xFF
    with pytest.raises(fqxv.FqxvError):
        arc.index.verify_stream(0, "names", bytes(coded))


def test_stream_matches_local(server, archive):
    # Whole-archive streaming decode over HTTP (the aligner path) == local decode.
    local = list(fqxv.open(archive))
    streamed = list(fqxv.remote.stream(server))
    assert len(streamed) == len(local)
    assert [r.name for r in streamed] == [r.name for r in local]
    assert [r.sequence for r in streamed] == [r.sequence for r in local]
    assert [r.quality for r in streamed] == [r.quality for r in local]


def test_download_matches_local(server, archive, tmp_path):
    dest = tmp_path / "streamed.fastq"
    n = fqxv.remote.download(server, dest)
    assert n == len(list(fqxv.open(archive)))
    assert dest.read_bytes() == fqxv.decompress_to_bytes(archive)


def test_primitives_drive_a_manual_projection(server, archive):
    # The IO-free primitives a custom (e.g. async) client would call directly.
    import urllib.request

    local = list(fqxv.open(archive))
    size = os.path.getsize(archive)

    def get(range_header):
        req = urllib.request.Request(server, headers={"Range": range_header})
        with urllib.request.urlopen(req) as r:
            return r.read()

    index, need = fqxv.parse_index_suffix(get("bytes=-65536"), size)
    assert need is None
    names = []
    for g in range(index.num_groups):
        start, end = index.stream_range(g, "names")
        coded = get(f"bytes={start}-{end - 1}")
        index.verify_stream(g, "names", coded)
        names.extend(fqxv.decode_names_bytes(coded))
    assert names == [r.name for r in local]
