//! Record-level decoding: pull one `(name, sequence, quality)` at a time instead
//! of a serialized FASTQ byte stream.
//!
//! [`decompress`] writes interleaved FASTQ *text* to a [`Write`] sink; nothing in
//! that path hands back structured records. This module adds that surface without
//! any codec-specific logic: it drives the ordinary [`decompress`] into a
//! [`Write`] sink that reassembles the emitted four-line records (the same trick
//! [`content_stats`](super::content_stats) uses via `StatsSink`). Because it rides
//! the canonical decode path it is **layout-complete for free** — plain, grouped,
//! and both globally-clustered reorder layouts all funnel through the same sink —
//! whereas decoding blocks directly (via `decode_block_contents`) would not handle
//! the footer-less reorder layout.
//!
//! Two entry points share one sink:
//! - [`decompress_records`] — a push primitive; invokes a closure per record.
//! - [`RecordReader`] — a pull [`Iterator`]; runs the decode on a background
//!   thread feeding a bounded channel, so records stream out with bounded memory.

use super::*;
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread::{self, JoinHandle};

/// One decoded FASTQ record.
///
/// `name` is the header with the leading `@` stripped (description included);
/// `seq` and `qual` are the raw sequence and quality bytes with no line endings.
/// These are exactly the bytes [`decompress`] would emit for the record, split
/// back apart — the `+` separator line is dropped (as it is on the byte path).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Read name and description (no leading `@`).
    pub name: Vec<u8>,
    /// Sequence bases.
    pub seq: Vec<u8>,
    /// Quality scores (same length as `seq`).
    pub qual: Vec<u8>,
}

/// Bounded channel depth between the decode thread and a [`RecordReader`]. Large
/// enough to keep the decoder busy across a rayon batch, small enough to bound
/// memory; `send` blocks when full, so this is the backpressure knob.
const CHANNEL_CAP: usize = 1024;

/// A [`Write`] sink that reassembles decoded interleaved FASTQ into [`Record`]s
/// and hands each to `emit`. Lines cycle name → sequence → `+` → quality; a record
/// is emitted once its quality line completes. Line reassembly buffers across
/// arbitrary chunk boundaries exactly like `StatsSink`.
struct RecordSink<F> {
    emit: F,
    /// Bytes of the current line seen so far (newline excluded).
    line: Vec<u8>,
    /// Which line of the current record: 0 name, 1 seq, 2 `+`, 3 qual.
    line_no: u8,
    /// Name/sequence stashed until the record's quality line completes it.
    name: Vec<u8>,
    seq: Vec<u8>,
}

impl<F: FnMut(Record) -> io::Result<()>> RecordSink<F> {
    fn new(emit: F) -> Self {
        RecordSink {
            emit,
            line: Vec::new(),
            line_no: 0,
            name: Vec::new(),
            seq: Vec::new(),
        }
    }

    /// Commit the current line at its record position, emitting a record when the
    /// quality line closes it. Buffers are moved (not copied) into the record.
    fn commit_line(&mut self) -> io::Result<()> {
        // Defensive, mirroring `StatsSink`: `write_record` emits `\n`-terminated
        // lines with no `\r`, so this is a no-op in practice.
        if self.line.last() == Some(&b'\r') {
            self.line.pop();
        }
        match self.line_no {
            0 => {
                let mut name = std::mem::take(&mut self.line);
                if name.first() == Some(&b'@') {
                    name.remove(0);
                }
                self.name = name;
            }
            1 => self.seq = std::mem::take(&mut self.line),
            2 => self.line.clear(), // `+` separator line
            _ => {
                let qual = std::mem::take(&mut self.line);
                let record = Record {
                    name: std::mem::take(&mut self.name),
                    seq: std::mem::take(&mut self.seq),
                    qual,
                };
                (self.emit)(record)?;
            }
        }
        self.line_no = (self.line_no + 1) % 4;
        Ok(())
    }
}

impl<F: FnMut(Record) -> io::Result<()>> Write for RecordSink<F> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut rest = buf;
        while let Some(nl) = rest.iter().position(|&b| b == b'\n') {
            self.line.extend_from_slice(&rest[..nl]);
            self.commit_line()?;
            rest = &rest[nl + 1..];
        }
        self.line.extend_from_slice(rest);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Decode an archive and invoke `on_record` for every record, in original file
/// order. Layout-agnostic (drives [`decompress`]). Returns the decode [`Stats`].
///
/// ```no_run
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let file = std::fs::File::open("reads.fqxv")?;
/// let mut n = 0u64;
/// fqxv::decompress_records(file, 0, |_rec| n += 1)?;
/// println!("{n} reads");
/// # Ok(()) }
/// ```
pub fn decompress_records<R: Read>(
    reader: R,
    threads: usize,
    mut on_record: impl FnMut(Record),
) -> Result<Stats> {
    let sink = RecordSink::new(move |rec| {
        on_record(rec);
        Ok(())
    });
    decompress(reader, sink, threads)
}

/// A pull [`Iterator`] over an archive's records.
///
/// The decode runs on a background thread that feeds a bounded channel; [`next`]
/// pulls from it, so records stream with bounded memory regardless of archive
/// size. Handles every layout (it drives [`decompress`]).
///
/// A decode error surfaces as the iterator's final `Some(Err(_))` — because the
/// terminal result (including a trailing-block corruption caught only at
/// end-of-stream) arrives when the decode thread finishes, not through the record
/// channel. Dropping the reader early (e.g. breaking out of a loop) cleanly stops
/// and joins the decode thread.
///
/// [`next`]: Iterator::next
#[derive(Debug)]
pub struct RecordReader {
    /// `None` once the channel has disconnected and the thread been joined.
    rx: Option<Receiver<Record>>,
    handle: Option<JoinHandle<Result<Stats>>>,
    /// Terminal `Ok` stats, captured on join for [`RecordReader::finish`].
    stats: Option<Stats>,
}

impl RecordReader {
    /// Start decoding `reader` on a background thread. `threads` matches
    /// [`decompress`] (0 = a default pool). The reader is moved into the thread,
    /// hence the `Send + 'static` bound; a `Box<dyn Read + Send>`, [`std::fs::File`],
    /// or `Cursor<Vec<u8>>` all satisfy it.
    pub fn new<R: Read + Send + 'static>(reader: R, threads: usize) -> Self {
        let (tx, rx) = sync_channel::<Record>(CHANNEL_CAP);
        let handle = thread::spawn(move || {
            let sink = RecordSink::new(|rec| {
                tx.send(rec).map_err(|_| {
                    io::Error::new(io::ErrorKind::BrokenPipe, "record receiver dropped")
                })
            });
            decompress(reader, sink, threads)
        });
        RecordReader {
            rx: Some(rx),
            handle: Some(handle),
            stats: None,
        }
    }

    /// Join the decode thread, returning its terminal [`Stats`] (or the decode
    /// error). Drains and discards any records not yet pulled. Intended as an
    /// alternative to iterating; after a decode error has already surfaced through
    /// [`Iterator::next`], this reports the last known [`Stats`] (default if none).
    pub fn finish(mut self) -> Result<Stats> {
        if let Some(rx) = self.rx.take() {
            while rx.recv().is_ok() {}
        }
        self.join()
            .unwrap_or_else(|| Ok(self.stats.unwrap_or_default()))
    }

    /// Join the decode thread if it is still owned, mapping a thread panic to an
    /// error. Returns `None` if already joined.
    fn join(&mut self) -> Option<Result<Stats>> {
        self.handle.take().map(|h| {
            h.join()
                .unwrap_or(Err(Error::Malformed("record decode thread panicked")))
        })
    }
}

impl Iterator for RecordReader {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        let rx = self.rx.as_ref()?;
        match rx.recv() {
            Ok(record) => Some(Ok(record)),
            // Channel closed: the decode thread finished (cleanly or with an
            // error that only the JoinHandle carries). Join and, on error, yield
            // it once — subsequent calls see `rx == None` and return `None`.
            Err(_) => {
                self.rx = None;
                match self.join() {
                    Some(Ok(stats)) => {
                        self.stats = Some(stats);
                        None
                    }
                    Some(Err(e)) => Some(Err(e)),
                    None => None,
                }
            }
        }
    }
}

impl Drop for RecordReader {
    fn drop(&mut self) {
        // Close the channel BEFORE joining: a decode thread blocked in `send`
        // unblocks (BrokenPipe → `decompress` unwinds), so an early break never
        // deadlocks the join.
        drop(self.rx.take());
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Build a `.fqxv` archive from raw FASTQ bytes under the given params.
    fn compress_bytes(input: &[u8], params: Params) -> Vec<u8> {
        let mut out = Vec::new();
        compress(input, &mut out, params).expect("compress");
        out
    }

    /// Reference records: what `decompress` emits, parsed back into `(name, seq, qual)`.
    fn expected_records(archive: &[u8]) -> Vec<Record> {
        let mut fastq = Vec::new();
        decompress(archive, &mut fastq, 1).expect("decompress");
        let mut recs = Vec::new();
        let mut lines = fastq.split(|&b| b == b'\n');
        while let Some(name) = lines.next() {
            if name.is_empty() {
                break; // trailing newline
            }
            let seq = lines.next().unwrap();
            let _plus = lines.next().unwrap();
            let qual = lines.next().unwrap();
            recs.push(Record {
                name: name[1..].to_vec(), // strip '@'
                seq: seq.to_vec(),
                qual: qual.to_vec(),
            });
        }
        recs
    }

    const SAMPLE: &[u8] = b"\
@r1 desc one\nACGTACGT\n+\nIIIIFFF#\n\
@r2 desc two\nNNGGCCTA\n+\n###IIIFF\n\
@r3\nTTAA\n+\nIIII\n";

    fn params_plain() -> Params {
        Params::default()
    }
    fn params_reorder_keep() -> Params {
        Params {
            reorder: true,
            keep_order: true,
            ..Params::default()
        }
    }
    fn params_reorder_discard() -> Params {
        Params {
            reorder: true,
            keep_order: false,
            ..Params::default()
        }
    }

    /// Collect via the pull iterator and the push callback; assert they agree with
    /// each other and with the reference `decompress` output.
    fn assert_agrees(archive: &[u8], check_names: bool) {
        let expected = expected_records(archive);

        let iter_recs: Vec<Record> = RecordReader::new(io::Cursor::new(archive.to_vec()), 1)
            .map(|r| r.expect("record"))
            .collect();
        let mut push_recs = Vec::new();
        decompress_records(archive, 1, |r| push_recs.push(r)).expect("decompress_records");

        assert_eq!(iter_recs, push_recs, "iterator and callback disagree");
        assert_eq!(iter_recs.len(), expected.len(), "record count");
        for (got, want) in iter_recs.iter().zip(&expected) {
            assert_eq!(got.seq, want.seq, "sequence");
            assert_eq!(got.qual, want.qual, "quality");
            if check_names {
                assert_eq!(got.name, want.name, "name");
            }
        }
    }

    #[test]
    fn roundtrip_plain() {
        assert_agrees(&compress_bytes(SAMPLE, params_plain()), true);
    }

    #[test]
    fn roundtrip_reorder_keep_order() {
        assert_agrees(&compress_bytes(SAMPLE, params_reorder_keep()), true);
    }

    #[test]
    fn roundtrip_reorder_discard_order() {
        // Discard-order may permute reads and regenerate names, so only the record
        // count and the (seq, qual) multiset are guaranteed; skip name equality.
        assert_agrees(&compress_bytes(SAMPLE, params_reorder_discard()), false);
    }

    #[test]
    fn records_reserialize_to_decompress_output() {
        let archive = compress_bytes(SAMPLE, params_plain());
        let mut reserialized = Vec::new();
        decompress_records(&archive[..], 1, |r| {
            write_record(&mut reserialized, &r.name, &r.seq, &r.qual)
        })
        .unwrap();
        let mut direct = Vec::new();
        decompress(&archive[..], &mut direct, 1).unwrap();
        assert_eq!(
            reserialized, direct,
            "records must re-serialize byte-identically"
        );
    }

    #[test]
    fn early_drop_does_not_hang() {
        // Many small blocks so the decode thread is mid-stream (would block on a
        // full channel) when we drop after one record.
        let mut input = Vec::new();
        for i in 0..5000u32 {
            input.extend_from_slice(format!("@r{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").as_bytes());
        }
        let archive = compress_bytes(
            &input,
            Params {
                block_reads: 16,
                ..Params::default()
            },
        );
        let mut reader = RecordReader::new(io::Cursor::new(archive), 1);
        let first = reader.next().expect("at least one record").expect("ok");
        assert_eq!(first.seq, b"ACGTACGTACGT");
        drop(reader); // Drop must not deadlock on the still-running decode thread.
    }

    #[test]
    fn corrupted_archive_yields_one_error_then_none() {
        let mut archive = compress_bytes(SAMPLE, params_plain());
        // Corrupt deep in the payload (past the 10-byte header) to trip a CRC.
        let mid = archive.len() / 2;
        archive[mid] ^= 0xFF;
        let mut reader = RecordReader::new(io::Cursor::new(archive), 1);
        let mut errors = 0;
        let mut sawn_none_after_error = false;
        while let Some(item) = reader.next() {
            if item.is_err() {
                errors += 1;
                // The very next pull must terminate the iterator.
                sawn_none_after_error = reader.next().is_none();
                break;
            }
        }
        assert!(errors >= 1, "corruption must surface as an error");
        assert!(sawn_none_after_error, "iterator must end after the error");
    }

    proptest! {
        #[test]
        fn prop_roundtrip_plain(
            reads in proptest::collection::vec(
                (
                    "[A-Za-z0-9:_.-]{1,20}",
                    proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')], 1..40),
                ),
                1..30,
            )
        ) {
            let mut input = Vec::new();
            for (i, (name, seq)) in reads.iter().enumerate() {
                input.extend_from_slice(b"@");
                input.extend_from_slice(name.as_bytes());
                // Force unique names so the tokenizer round-trips them verbatim.
                input.extend_from_slice(format!("_{i}\n").as_bytes());
                input.extend_from_slice(seq);
                input.push(b'\n');
                input.extend_from_slice(b"+\n");
                // Deterministic quality of the matching length.
                input.extend(seq.iter().map(|_| b'I'));
                input.push(b'\n');
            }
            let archive = compress_bytes(&input, params_plain());
            let recs: Vec<Record> = RecordReader::new(io::Cursor::new(archive.clone()), 1)
                .map(|r| r.expect("record"))
                .collect();
            let expected = expected_records(&archive);
            prop_assert_eq!(recs, expected);
        }
    }
}
