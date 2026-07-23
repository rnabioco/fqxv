"""Read-only Python bindings for the ``fqxv`` FASTQ archiver.

The compiled core lives in :mod:`fqxv._fqxv`; this package re-exports it and adds
:mod:`fqxv.remote`, which reads an archive over HTTP using byte-range requests —
fetching just the footer index and then only the column(s) you ask for, rather
than downloading the whole file.
"""

from ._fqxv import (
    Estimate,
    FqxvError,
    GroupLoc,
    Index,
    Info,
    Reader,
    Record,
    decode_names_bytes,
    decode_qualities_bytes,
    decode_sequences_bytes,
    decompress_to_bytes,
    decompress_to_path,
    estimate,
    inspect,
    open,
    open_index,
    parse_index_suffix,
    quality_needs_sequence_bytes,
    read_block,
    read_names,
    read_qualities,
    read_sequences,
    verify,
)
from . import remote

__all__ = [
    # Classes
    "Estimate",
    "FqxvError",
    "GroupLoc",
    "Index",
    "Info",
    "Reader",
    "Record",
    # Whole-archive / streaming / projection (local path or bytes)
    "open",
    "decompress_to_path",
    "decompress_to_bytes",
    "inspect",
    "open_index",
    "read_names",
    "read_sequences",
    "read_qualities",
    "read_block",
    "estimate",
    "verify",
    # Random-access primitives (used by fqxv.remote; rarely called directly)
    "parse_index_suffix",
    "decode_names_bytes",
    "decode_sequences_bytes",
    "decode_qualities_bytes",
    "quality_needs_sequence_bytes",
    # Remote (HTTP Range) reader
    "remote",
]
