"""Read an ``.fqxv`` archive over HTTP without downloading the whole file.

Two capabilities, both driven by the container's footer index (per-stream
``(offset, length, crc32c)`` reachable from a fixed 12-byte EOF trailer):

* **Streaming** the whole archive — :func:`stream` / :func:`download` feed an HTTP
  response straight into the decoder, so reads flow out as blocks arrive. This is
  the whole-file path; for a private object presign the URL or pass an
  ``Authorization`` header. (Equivalent to the CLI's
  ``aws s3 cp s3://… - | fqxv decompress -``.)
* **Column projection** — :class:`RemoteArchive` fetches just the footer index from
  the archive tail, then only the column(s) you ask for with HTTP ``Range``
  requests (read names are typically <1% of the file), CRC-verified per stream. A
  copy tool can't do this — it downloads the whole object — which is why this lives
  here rather than delegating to ``aws s3 cp``.

The HTTP client is the standard library (``urllib``), so this module has no
dependencies. Everything it does rests on IO-free primitives in the compiled core
(``fqxv.parse_index_suffix``, ``Index.stream_range`` / ``verify_stream``,
``fqxv.decode_{names,sequences,qualities}_bytes``); to drive it with a different
client — an async ``httpx``/``aiohttp`` session for concurrent range fetches, or a
``boto3`` client — call those primitives directly (see :class:`RemoteArchive` for
the shape). A short async example:

    idx_bytes = (await client.get(url, headers={"Range": "bytes=-65536"})).content
    index, need = fqxv.parse_index_suffix(idx_bytes, total_size)
    start, end = index.stream_range(0, "names")
    coded = (await client.get(url, headers={"Range": f"bytes={start}-{end-1}"})).content
    index.verify_stream(0, "names", coded)
    names = fqxv.decode_names_bytes(coded)

The globally-reordered layout (``--order shuffle``) has no footer index and cannot
be projected; :func:`fqxv.parse_index_suffix` raises for it. Use :func:`stream`.
"""

from __future__ import annotations

import urllib.request
from collections import namedtuple
from typing import Iterable, List, Optional, Tuple

from . import _fqxv

__all__ = [
    "RemoteRecord",
    "RemoteArchive",
    "open_index",
    "read_names",
    "read_sequences",
    "read_qualities",
    "stream",
    "download",
]

# A projected record; field names mirror the compiled `Record`.
RemoteRecord = namedtuple("RemoteRecord", ["name", "sequence", "quality"])

# First tail fetch: big enough that the footer almost always fits in one round trip.
DEFAULT_TAIL = 64 << 10


def _total_size(content_range: Optional[str], body_len: int) -> int:
    """Total archive size from a 206 ``Content-Range`` (``bytes 100-199/12345``),
    or the body length when the server ignored the range and sent the whole file."""
    if content_range:
        tail = content_range.rsplit("/", 1)[-1].strip()
        if tail.isdigit():
            return int(tail)
    return body_len


def _decode_column(stream_name: str, coded: bytes, seq: Optional[bytes] = None) -> List[bytes]:
    if stream_name == "names":
        return _fqxv.decode_names_bytes(coded)
    if stream_name == "sequence":
        return _fqxv.decode_sequences_bytes(coded)
    if stream_name == "quality":
        return _fqxv.decode_qualities_bytes(coded, seq)
    raise ValueError(f"unknown stream {stream_name!r}; expected names, sequence, or quality")


# --------------------------------------------------------------------------- #
# Column projection over HTTP Range (synchronous, standard library only)
# --------------------------------------------------------------------------- #
class RemoteArchive:
    """A remote reader over one archive URL, holding the parsed footer index so each
    column read costs one range request (no re-fetch of the tail). ``bytes_fetched``
    tracks the total transferred, which stays a small fraction of ``size`` for a
    projection."""

    def __init__(self, url: str, index, size: int, *, headers=None, bytes_fetched: int = 0):
        self.url = url
        self.index = index
        self.size = size
        self.headers = dict(headers or {})
        self.bytes_fetched = bytes_fetched

    def _get(self, range_header: str) -> Tuple[bytes, Optional[str]]:
        headers = dict(self.headers)
        headers["Range"] = range_header
        req = urllib.request.Request(self.url, headers=headers)
        with urllib.request.urlopen(req) as resp:
            return resp.read(), resp.headers.get("Content-Range")

    def _get_range(self, start: int, end: int) -> bytes:
        # HTTP byte ranges are inclusive on both ends; ours are half-open.
        body, _ = self._get(f"bytes={start}-{end - 1}")
        self.bytes_fetched += len(body)
        return body

    @classmethod
    def open(cls, url: str, *, tail: int = DEFAULT_TAIL, headers=None) -> "RemoteArchive":
        """Open ``url``, fetching only the archive tail to parse the footer index."""
        arc = cls(url, index=None, size=0, headers=headers)
        body, content_range = arc._get(f"bytes=-{tail}")
        arc.bytes_fetched += len(body)
        size = _total_size(content_range, len(body))
        index, need = _fqxv.parse_index_suffix(body, size)
        if need is not None:
            body, content_range = arc._get(f"bytes=-{need}")
            arc.bytes_fetched += len(body)
            size = _total_size(content_range, len(body))
            index, need = _fqxv.parse_index_suffix(body, size)
            if need is not None:
                raise _fqxv.FqxvError("footer did not fit in the refetched tail")
        arc.index, arc.size = index, size
        return arc

    def _column(self, stream_name: str, group: int) -> bytes:
        start, end = self.index.stream_range(group, stream_name)
        coded = self._get_range(start, end)
        self.index.verify_stream(group, stream_name, coded)
        return coded

    def _groups(self, groups: Optional[Iterable[int]]) -> List[int]:
        return list(range(self.index.num_groups)) if groups is None else list(groups)

    def names(self, groups: Optional[Iterable[int]] = None) -> List[bytes]:
        """Fetch and decode the read names (one range GET per row group)."""
        out: List[bytes] = []
        for g in self._groups(groups):
            out.extend(_decode_column("names", self._column("names", g)))
        return out

    def sequences(self, groups: Optional[Iterable[int]] = None) -> List[bytes]:
        """Fetch and decode the sequences, skipping the quality stream entirely."""
        out: List[bytes] = []
        for g in self._groups(groups):
            out.extend(_decode_column("sequence", self._column("sequence", g)))
        return out

    def qualities(self, groups: Optional[Iterable[int]] = None) -> List[bytes]:
        """Fetch and decode the qualities. Long-read quality is coded against the
        sequence, so those groups also fetch their sequence column."""
        out: List[bytes] = []
        for g in self._groups(groups):
            qcoded = self._column("quality", g)
            seq = None
            if _fqxv.quality_needs_sequence_bytes(qcoded):
                seq = b"".join(_decode_column("sequence", self._column("sequence", g)))
            out.extend(_decode_column("quality", qcoded, seq))
        return out

    def records(self, groups: Optional[Iterable[int]] = None) -> List[RemoteRecord]:
        """Fetch all three columns for the selected groups and zip them into
        :class:`RemoteRecord` tuples (this reads the whole archive body)."""
        out: List[RemoteRecord] = []
        for g in self._groups(groups):
            names = _decode_column("names", self._column("names", g))
            seqs = _decode_column("sequence", self._column("sequence", g))
            qcoded = self._column("quality", g)
            seq = b"".join(seqs) if _fqxv.quality_needs_sequence_bytes(qcoded) else None
            quals = _decode_column("quality", qcoded, seq)
            out.extend(RemoteRecord(n, s, q) for n, s, q in zip(names, seqs, quals))
        return out


def open_index(url: str, *, tail: int = DEFAULT_TAIL, headers=None) -> RemoteArchive:
    """Open a remote archive and parse its footer index (one tail GET)."""
    return RemoteArchive.open(url, tail=tail, headers=headers)


def read_names(url: str, groups: Optional[Iterable[int]] = None, **kw) -> List[bytes]:
    return RemoteArchive.open(url, **kw).names(groups)


def read_sequences(url: str, groups: Optional[Iterable[int]] = None, **kw) -> List[bytes]:
    return RemoteArchive.open(url, **kw).sequences(groups)


def read_qualities(url: str, groups: Optional[Iterable[int]] = None, **kw) -> List[bytes]:
    return RemoteArchive.open(url, **kw).qualities(groups)


# --------------------------------------------------------------------------- #
# Streaming whole-archive decode (mirrors `fqxv decompress -`)
# --------------------------------------------------------------------------- #
def _urlopen(url: str, headers=None):
    return urllib.request.urlopen(urllib.request.Request(url, headers=dict(headers or {})))


def stream(url: str, *, headers=None, threads: int = 0):
    """Iterate records from a remote archive without staging it to disk. The HTTP
    response feeds straight into the streaming decoder, so records flow out as
    blocks arrive. Yields fqxv ``Record`` objects in original order.

    For a private object presign the URL or pass an ``Authorization`` header via
    ``headers``. For ``s3://`` with the AWS credential chain, hand a boto3 body to
    :func:`fqxv.open` directly (``fqxv.open(s3.get_object(...)["Body"])``) or use the
    CLI (``aws s3 cp s3://… - | fqxv decompress -``)."""
    return _fqxv.open(_urlopen(url, headers), threads=threads)


def download(url: str, dest, *, headers=None, threads: int = 0) -> int:
    """Stream a remote archive and decode it to interleaved FASTQ at ``dest``
    without buffering the archive in memory. Returns the read count."""
    return _fqxv.decompress_to_path(_urlopen(url, headers), str(dest), threads=threads)
