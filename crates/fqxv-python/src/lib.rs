//! Read-only Python bindings for the `fqxv` FASTQ archiver.
//!
//! Exposes three things: a streaming record iterator ([`open`] → [`Reader`]),
//! whole-archive convenience ([`decompress_to_path`], [`decompress_to_bytes`],
//! [`inspect`]), and column projection / random access over the v3 footer index
//! ([`open_index`], [`read_names`], [`read_sequences`], [`read_qualities`],
//! [`read_block`]). Every entry point accepts a filesystem path (`str` /
//! `os.PathLike`) or in-memory `bytes`.
#![allow(missing_docs)] // PyO3 glue; user-facing docs live in the .pyi / Python layer.

use std::fs::File;
use std::io::{self, Cursor, Read, Seek, SeekFrom};
use std::path::PathBuf;

use pyo3::create_exception;
use pyo3::exceptions::{PyException, PyIOError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use fqxv_core::{
    decode_block_contents, decode_names, decode_quality, decode_sequence, Index, Stream,
};

create_exception!(
    fqxv,
    FqxvError,
    PyException,
    "Base error for fqxv archive failures."
);

/// Map an [`fqxv_core::Error`] onto a Python exception. I/O failures become
/// `OSError`; everything else becomes `fqxv.FqxvError`. The enum is
/// `#[non_exhaustive]`, hence the catch-all.
fn map_err(e: fqxv_core::Error) -> PyErr {
    match e {
        fqxv_core::Error::Io(io) => PyIOError::new_err(io.to_string()),
        other => FqxvError::new_err(other.to_string()),
    }
}

/// A concrete seekable reader over either an in-memory buffer or a file, so the
/// generic decode helpers monomorphize once instead of over a trait object.
enum SeekReader {
    Bytes(Cursor<Vec<u8>>),
    File(File),
}

impl Read for SeekReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            SeekReader::Bytes(c) => c.read(buf),
            SeekReader::File(f) => f.read(buf),
        }
    }
}

impl Seek for SeekReader {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match self {
            SeekReader::Bytes(c) => c.seek(pos),
            SeekReader::File(f) => f.seek(pos),
        }
    }
}

/// Open a Python `source` (`bytes`, `str`, or `os.PathLike`) as a seekable reader.
fn open_source(obj: &Bound<'_, PyAny>) -> PyResult<SeekReader> {
    if let Ok(b) = obj.downcast::<PyBytes>() {
        return Ok(SeekReader::Bytes(Cursor::new(b.as_bytes().to_vec())));
    }
    let path: PathBuf = obj
        .extract()
        .map_err(|_| PyTypeError::new_err("source must be bytes, a str path, or an os.PathLike"))?;
    let file =
        File::open(&path).map_err(|e| PyIOError::new_err(format!("{}: {e}", path.display())))?;
    Ok(SeekReader::File(file))
}

/// One decoded FASTQ record. `name` excludes the leading `@`; `sequence` and
/// `quality` are raw bytes with no line endings.
#[pyclass(name = "Record", frozen)]
struct PyRecord {
    inner: fqxv_core::Record,
}

#[pymethods]
impl PyRecord {
    #[getter]
    fn name<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.inner.name)
    }
    #[getter]
    fn sequence<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.inner.seq)
    }
    #[getter]
    fn quality<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.inner.qual)
    }
    fn __len__(&self) -> usize {
        self.inner.seq.len()
    }
    fn __repr__(&self) -> String {
        format!(
            "Record(name={:?}, len={})",
            String::from_utf8_lossy(&self.inner.name),
            self.inner.seq.len()
        )
    }
}

/// A streaming iterator over an archive's records, in original file order.
/// Handles every layout, including globally-reordered archives.
///
/// `unsendable`: the reader owns a channel `Receiver` (`!Sync`), so the object is
/// pinned to the thread that created it — the natural usage for an iterator.
#[pyclass(name = "Reader", unsendable)]
struct PyReader {
    inner: fqxv_core::RecordReader,
}

#[pymethods]
impl PyReader {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyRecord>> {
        let py = slf.py();
        // Release the GIL while blocking on the decode thread's channel. The
        // outstanding PyRefMut guard still blocks any concurrent borrow of this
        // Reader, so the &mut stays exclusive.
        let inner = &mut slf.inner;
        match py.allow_threads(move || inner.next()) {
            None => Ok(None),
            Some(Ok(rec)) => Ok(Some(PyRecord { inner: rec })),
            Some(Err(e)) => Err(map_err(e)),
        }
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    #[pyo3(signature = (*_args))]
    fn __exit__(&self, _args: &Bound<'_, PyAny>) -> bool {
        false
    }
}

/// Metadata about an archive, from its header and footer index (no payload decode).
#[pyclass(name = "Info")]
struct PyInfo {
    inner: fqxv_core::Info,
}

#[pymethods]
impl PyInfo {
    #[getter]
    fn reads(&self) -> u64 {
        self.inner.reads
    }
    #[getter]
    fn blocks(&self) -> u64 {
        self.inner.blocks
    }
    #[getter]
    fn group_size(&self) -> u8 {
        self.inner.group_size
    }
    #[getter]
    fn reordered(&self) -> bool {
        self.inner.reordered
    }
    #[getter]
    fn keep_order(&self) -> bool {
        self.inner.keep_order
    }
    #[getter]
    fn regenerated_names(&self) -> bool {
        self.inner.regenerated_names
    }
    #[getter]
    fn plus_normalized(&self) -> bool {
        self.inner.plus_normalized
    }
    #[getter]
    fn format_version(&self) -> u16 {
        self.inner.format_version
    }
    #[getter]
    fn seq_order(&self) -> u8 {
        self.inner.seq_order
    }
    #[getter]
    fn quality_binning(&self) -> u8 {
        self.inner.quality_binning
    }
    #[getter]
    fn names_bytes(&self) -> u64 {
        self.inner.names_bytes
    }
    #[getter]
    fn sequence_bytes(&self) -> u64 {
        self.inner.seq_bytes
    }
    #[getter]
    fn quality_bytes(&self) -> u64 {
        self.inner.qual_bytes
    }
    #[getter]
    fn platform(&self) -> &'static str {
        self.inner.platform.label()
    }
    #[getter]
    fn whole_file_crc(&self) -> Option<u32> {
        self.inner.whole_file_crc
    }
    fn __repr__(&self) -> String {
        format!(
            "Info(reads={}, blocks={}, format_version={}, platform={:?})",
            self.inner.reads,
            self.inner.blocks,
            self.inner.format_version,
            self.inner.platform.label(),
        )
    }
}

/// On-disk location of one row group.
#[pyclass(name = "GroupLoc", frozen)]
struct PyGroupLoc {
    #[pyo3(get)]
    block_offset: u64,
    #[pyo3(get)]
    read_count: u32,
}

/// A parsed footer row-group index (plain layout only).
#[pyclass(name = "Index")]
struct PyIndex {
    total_reads: u64,
    whole_file_crc: u32,
    groups: Vec<PyGroupLoc>,
}

#[pymethods]
impl PyIndex {
    #[getter]
    fn total_reads(&self) -> u64 {
        self.total_reads
    }
    #[getter]
    fn whole_file_crc(&self) -> u32 {
        self.whole_file_crc
    }
    #[getter]
    fn num_groups(&self) -> usize {
        self.groups.len()
    }
    fn groups(&self, py: Python<'_>) -> Vec<Py<PyGroupLoc>> {
        self.groups
            .iter()
            .map(|g| {
                Py::new(
                    py,
                    PyGroupLoc {
                        block_offset: g.block_offset,
                        read_count: g.read_count,
                    },
                )
            })
            .collect::<PyResult<Vec<_>>>()
            .expect("GroupLoc alloc")
    }
    fn __repr__(&self) -> String {
        format!(
            "Index(total_reads={}, num_groups={})",
            self.total_reads,
            self.groups.len()
        )
    }
}

/// Slice a concatenated stream into per-read byte vectors using `lengths`.
fn slice_reads(lengths: &[u32], data: &[u8], out: &mut Vec<Vec<u8>>) {
    let mut off = 0usize;
    for &l in lengths {
        let l = l as usize;
        out.push(data[off..off + l].to_vec());
        off += l;
    }
}

/// Fetch and decode one projected stream across the requested groups (or all).
fn project(
    mut r: SeekReader,
    stream: Stream,
    groups: Option<Vec<usize>>,
) -> Result<Vec<Vec<u8>>, fqxv_core::Error> {
    let index = Index::read(&mut r)?;
    let selected: Vec<usize> = groups.unwrap_or_else(|| (0..index.groups().len()).collect());
    let ranges = index.byte_ranges(&selected, stream)?;
    let mut out = Vec::new();
    for (&g, range) in selected.iter().zip(ranges) {
        r.seek(SeekFrom::Start(range.start))?;
        let mut buf = vec![0u8; (range.end - range.start) as usize];
        r.read_exact(&mut buf)?;
        index.verify_stream(g, stream, &buf)?;
        match stream {
            Stream::Names => out.extend(decode_names(&buf)?),
            Stream::Sequence => {
                let (lens, bases) = decode_sequence(&buf)?;
                slice_reads(&lens, &bases, &mut out);
            }
            Stream::Quality => {
                let (lens, qual) = decode_quality(&buf)?;
                slice_reads(&lens, &qual, &mut out);
            }
        }
    }
    Ok(out)
}

/// Fetch and decode one whole block into records.
fn read_one_block(
    mut r: SeekReader,
    group: usize,
) -> Result<Vec<fqxv_core::Record>, fqxv_core::Error> {
    let index = Index::read(&mut r)?;
    let loc = *index
        .groups()
        .get(group)
        .ok_or(fqxv_core::Error::Malformed("row-group index out of range"))?;
    r.seek(SeekFrom::Start(loc.block_offset))?;
    // Block frame head: [4 BLOCK_MAGIC][8 payload_len][4 crc]. decode_block_contents
    // runs its own content-digest check, so we only need the payload length here.
    let mut head = [0u8; 16];
    r.read_exact(&mut head)?;
    let len = u64::from_le_bytes(head[4..12].try_into().unwrap()) as usize;
    let mut payload = vec![0u8; len];
    r.read_exact(&mut payload)?;
    let bc = decode_block_contents(&payload)?;
    let mut recs = Vec::with_capacity(bc.names.len());
    let mut off = 0usize;
    for (i, name) in bc.names.into_iter().enumerate() {
        let l = bc.lengths[i] as usize;
        recs.push(fqxv_core::Record {
            name,
            seq: bc.sequence[off..off + l].to_vec(),
            qual: bc.quality[off..off + l].to_vec(),
        });
        off += l;
    }
    Ok(recs)
}

/// Convert a list of raw byte vectors into a Python `list[bytes]`.
fn to_bytes_list<'py>(py: Python<'py>, items: Vec<Vec<u8>>) -> Vec<Bound<'py, PyBytes>> {
    items.iter().map(|v| PyBytes::new(py, v)).collect()
}

/// Open an archive for streaming record iteration.
#[pyfunction]
#[pyo3(signature = (source, *, threads = 0))]
fn open(source: &Bound<'_, PyAny>, threads: usize) -> PyResult<PyReader> {
    let reader = open_source(source)?;
    Ok(PyReader {
        inner: fqxv_core::RecordReader::new(reader, threads),
    })
}

/// Decompress an archive to interleaved FASTQ at `dest`. Returns the read count.
#[pyfunction]
#[pyo3(signature = (source, dest, *, threads = 0))]
fn decompress_to_path(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    dest: PathBuf,
    threads: usize,
) -> PyResult<u64> {
    let reader = open_source(source)?;
    let file =
        File::create(&dest).map_err(|e| PyIOError::new_err(format!("{}: {e}", dest.display())))?;
    let stats = py
        .allow_threads(move || fqxv_core::decompress(reader, io::BufWriter::new(file), threads))
        .map_err(map_err)?;
    Ok(stats.reads)
}

/// Decompress an archive and return interleaved FASTQ as `bytes`.
#[pyfunction]
#[pyo3(signature = (source, *, threads = 0))]
fn decompress_to_bytes<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    threads: usize,
) -> PyResult<Bound<'py, PyBytes>> {
    let reader = open_source(source)?;
    let out = py
        .allow_threads(move || {
            let mut buf = Vec::new();
            fqxv_core::decompress(reader, &mut buf, threads).map(|_| buf)
        })
        .map_err(map_err)?;
    Ok(PyBytes::new(py, &out))
}

/// Read an archive's metadata (header + footer index; no payload decode).
#[pyfunction]
fn inspect(py: Python<'_>, source: &Bound<'_, PyAny>) -> PyResult<PyInfo> {
    let reader = open_source(source)?;
    let info = py
        .allow_threads(move || fqxv_core::inspect(reader))
        .map_err(map_err)?;
    Ok(PyInfo { inner: info })
}

/// Parse the footer row-group index (plain layout only; rejects reorder archives).
#[pyfunction]
fn open_index(py: Python<'_>, source: &Bound<'_, PyAny>) -> PyResult<PyIndex> {
    let mut reader = open_source(source)?;
    let index = py
        .allow_threads(move || Index::read(&mut reader))
        .map_err(map_err)?;
    let groups = index
        .groups()
        .iter()
        .map(|g| PyGroupLoc {
            block_offset: g.block_offset,
            read_count: g.read_count,
        })
        .collect();
    Ok(PyIndex {
        total_reads: index.total_reads(),
        whole_file_crc: index.whole_file_crc(),
        groups,
    })
}

/// Decode just the read names for the given row groups (or all). `list[bytes]`.
#[pyfunction]
#[pyo3(signature = (source, groups = None))]
fn read_names<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    groups: Option<Vec<usize>>,
) -> PyResult<Vec<Bound<'py, PyBytes>>> {
    let reader = open_source(source)?;
    let data = py
        .allow_threads(move || project(reader, Stream::Names, groups))
        .map_err(map_err)?;
    Ok(to_bytes_list(py, data))
}

/// Decode just the sequences for the given row groups (or all). `list[bytes]`.
#[pyfunction]
#[pyo3(signature = (source, groups = None))]
fn read_sequences<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    groups: Option<Vec<usize>>,
) -> PyResult<Vec<Bound<'py, PyBytes>>> {
    let reader = open_source(source)?;
    let data = py
        .allow_threads(move || project(reader, Stream::Sequence, groups))
        .map_err(map_err)?;
    Ok(to_bytes_list(py, data))
}

/// Decode just the qualities for the given row groups (or all). `list[bytes]`.
#[pyfunction]
#[pyo3(signature = (source, groups = None))]
fn read_qualities<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    groups: Option<Vec<usize>>,
) -> PyResult<Vec<Bound<'py, PyBytes>>> {
    let reader = open_source(source)?;
    let data = py
        .allow_threads(move || project(reader, Stream::Quality, groups))
        .map_err(map_err)?;
    Ok(to_bytes_list(py, data))
}

/// Random access: decode one whole row group into a `list[Record]`.
#[pyfunction]
fn read_block(py: Python<'_>, source: &Bound<'_, PyAny>, group: usize) -> PyResult<Vec<PyRecord>> {
    let reader = open_source(source)?;
    let recs = py
        .allow_threads(move || read_one_block(reader, group))
        .map_err(map_err)?;
    Ok(recs.into_iter().map(|inner| PyRecord { inner }).collect())
}

#[pymodule]
fn fqxv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRecord>()?;
    m.add_class::<PyReader>()?;
    m.add_class::<PyInfo>()?;
    m.add_class::<PyIndex>()?;
    m.add_class::<PyGroupLoc>()?;
    m.add("FqxvError", m.py().get_type::<FqxvError>())?;
    m.add_function(wrap_pyfunction!(open, m)?)?;
    m.add_function(wrap_pyfunction!(decompress_to_path, m)?)?;
    m.add_function(wrap_pyfunction!(decompress_to_bytes, m)?)?;
    m.add_function(wrap_pyfunction!(inspect, m)?)?;
    m.add_function(wrap_pyfunction!(open_index, m)?)?;
    m.add_function(wrap_pyfunction!(read_names, m)?)?;
    m.add_function(wrap_pyfunction!(read_sequences, m)?)?;
    m.add_function(wrap_pyfunction!(read_qualities, m)?)?;
    m.add_function(wrap_pyfunction!(read_block, m)?)?;
    Ok(())
}
