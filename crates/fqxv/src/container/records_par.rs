//! Parallel record visiting: a map-reduce hook over an archive's decoded
//! records that scales with the container's block parallelism.
//!
//! The serial record surface ([`decompress_records`](super::decompress_records),
//! [`RecordReader`](super::RecordReader)) funnels every record through one
//! consumer — the decode fans out across blocks on rayon, then re-serializes
//! into a single callback or channel, so a record-level consumer (k-mer
//! counting, filtering, per-read stats) is capped at one core no matter how
//! wide the decode ran. This module adds [`decompress_records_par`], a
//! per-worker-state map-reduce visitor in the shape BINSEQ showed scales
//! near-linearly: blocks decode in parallel and each block's records are
//! visited *inside* its decode task, so the visitor work scales with the block
//! parallelism instead of collapsing onto a channel drain.
//!
//! Records are delivered as borrowed [`RecordRef`] views into the decoded
//! block's buffers — no per-record allocation on the delivery path — each
//! tagged with its record index, so an order-sensitive consumer can restore
//! the serial order afterwards while an order-free reduction never pays for it.
//!
//! [`decompress_records_par_select`] composes the visitor with a
//! [`StreamSelection`] (issue #272): deselected streams are genuinely skipped
//! *inside* the per-block decode tasks — never entropy-decoded — so the
//! per-stream savings measured for the serial selective decode (~12x
//! sequence-only on ONT, where the sequence-conditioned quality model dominates
//! decode; ~1.3-1.6x single-threaded on short-read Illumina) stack with the
//! visitor's block-parallel scaling.

use super::*;
use rayon::prelude::*;
use std::sync::{Mutex, PoisonError};

/// A borrowed view of one decoded FASTQ record, valid only for the duration of
/// a visitor call.
///
/// Field semantics match [`Record`]: `name` is the header with the leading `@`
/// stripped (description included); `seq` and `qual` are the raw sequence and
/// quality bytes (equal length, no line endings). The `+` separator line is
/// normalized away, exactly as on every other decode path. Borrowing is what
/// keeps the parallel visitor allocation-free per record; call
/// [`RecordRef::to_record`] to keep a copy past the visit.
#[derive(Debug, Clone, Copy)]
pub struct RecordRef<'a> {
    /// Read name and description (no leading `@`).
    pub name: &'a [u8],
    /// Sequence bases.
    pub seq: &'a [u8],
    /// Quality scores (same length as `seq`).
    pub qual: &'a [u8],
}

impl RecordRef<'_> {
    /// Copy this view into an owned [`Record`].
    pub fn to_record(&self) -> Record {
        Record {
            name: self.name.to_vec(),
            seq: self.seq.to_vec(),
            qual: self.qual.to_vec(),
        }
    }
}

impl From<RecordRef<'_>> for Record {
    fn from(r: RecordRef<'_>) -> Record {
        r.to_record()
    }
}

/// Minimum records per parallel visit task on the whole-file reorder path.
/// Small enough that even a modest archive fans out, large enough that the
/// per-task state checkout is noise against the visit work.
const MIN_VISIT_CHUNK: usize = 256;

/// A lock-guarded pool of visitor states, checked out per decode task and
/// merged once at the end. States are created lazily, so at most one exists
/// per concurrently running task (bounded by the thread count) — the
/// per-worker-state model — rather than one per block.
struct StatePool<'a, S, I> {
    states: Mutex<Vec<S>>,
    init: &'a I,
}

impl<'a, S, I: Fn() -> S> StatePool<'a, S, I> {
    fn new(init: &'a I) -> Self {
        StatePool {
            states: Mutex::new(Vec::new()),
            init,
        }
    }

    /// Check a state out, creating one if none is idle. The lock is held only
    /// for the pop, never across visitor calls.
    fn take(&self) -> S {
        self.states
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .pop()
            .unwrap_or_else(self.init)
    }

    /// Return a state after a task finishes with it.
    fn put(&self, s: S) {
        self.states
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .push(s);
    }

    /// Fold every state into one with `merge` (serially, on the caller's
    /// thread). An archive that delivered no records yields a fresh state.
    fn into_merged<M: FnMut(S, S) -> S>(self, merge: M) -> S {
        self.states
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner)
            .into_iter()
            .reduce(merge)
            .unwrap_or_else(self.init)
    }
}

/// Decode an archive and fan its records out to a map-reduce visitor in
/// parallel, returning the merged state and the decode [`Stats`].
///
/// Equivalent to [`decompress_records_par_select`] with
/// [`StreamSelection::ALL`] (use that entry point to skip whole streams — e.g.
/// a parallel sequence-only pass); the one difference is bookkeeping:
/// [`Stats::out_bytes`] here is the serialized-FASTQ byte count, as documented
/// below.
///
/// Three hooks define the reduction, mirroring rayon's `fold`/`reduce` shape
/// (chosen over a visitor trait because every other hook in this API surface is
/// a closure, and a closure triple needs no impl boilerplate for the common
/// counting/filtering consumers):
///
/// - `init` builds one worker state. It is called lazily, at most once per
///   concurrently running task (bounded by the thread count), and each state is
///   reused across many blocks — so a large state (a k-mer table, a histogram)
///   is paid per *worker*, not per block or per record.
/// - `visit` is called once per record with the worker state, the record's
///   index, and a borrowed [`RecordRef`] valid only for that call. Returning an
///   error aborts the decode and propagates; construct one via
///   `fqxv::Error::from(std::io::Error::other(..))` for consumer-side failures.
/// - `merge` folds two states into one; it runs serially at the end, on the
///   caller's thread.
///
/// `threads` sizes the rayon pool as in [`decompress`] (0 = all cores).
///
/// # Index and ordering contract
///
/// The index passed to `visit` is the record's position in the **serial decode
/// order** — exactly the position the same record occupies in a
/// [`RecordReader`](super::RecordReader) iteration of the same archive. For
/// plain, grouped, and keep-order reordered archives that is the original input
/// order; for discard-order (`--order shuffle`) archives it is the archive's
/// clustered output order, the only order such an archive can reproduce. For a
/// grouped archive (group size `G` > 1, see [`Stats::group_size`]) records are
/// spot-interleaved, so member identity is `index % G`.
///
/// Records arrive as **disjoint contiguous runs of ascending indices** (a whole
/// block on the plain layout; a fixed index range on the reorder layouts). Runs
/// are visited concurrently and complete in arbitrary order; each state sees
/// its runs strictly sequentially, so worker state needs no internal
/// synchronization. Every index in `0..stats.reads` is visited exactly once on
/// a successful return.
///
/// The partitioning of records across states is **not deterministic** run to
/// run (it follows rayon's scheduling), so the merged result is deterministic
/// only when the reduction is insensitive to partitioning — sums, min/max,
/// mergeable tables, and other commutative folds all qualify. An
/// order-sensitive consumer should instead collect `(index, record)` pairs and
/// sort, or use the serial [`RecordReader`](super::RecordReader).
///
/// # Layout behavior
///
/// Every layout is delivered in parallel; none falls back to a serial funnel:
///
/// - **Plain / grouped / long-read blocks** (footer-carrying layout): blocks
///   decode in parallel and each block's records are visited inside its decode
///   task. Per-block stream digests are verified before that block's records
///   are visited.
/// - **Whole-file reorder layouts** (both keep-order and discard-order): the
///   entropy decode is block-parallel and the per-record normalization
///   (un-permuting, un-flipping, name regeneration) fans out too; records are
///   then visited in parallel over index ranges. The archive's whole-output
///   digest is verified before the first visit — a bounded, hash-only serial
///   pass, the layout's one serial section.
///
/// # Errors and panics
///
/// A decode or visitor error aborts the run and propagates. Because delivery is
/// concurrent, the visitor may already have observed an arbitrary subset of
/// records when the error surfaces (the same partial-delivery caveat the push
/// [`decompress_records`](super::decompress_records) has); no record is visited
/// twice regardless. A panicking visitor propagates the panic to the caller.
///
/// The returned [`Stats`] mirror a full decode: `reads`/`blocks` count the
/// archive, `out_bytes` is the interleaved-FASTQ byte count an equivalent
/// [`decompress`] would have written, and `group_size` is the archive's
/// recorded interleaving (1 for single-end).
///
/// ```no_run
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let file = std::fs::File::open("reads.fqxv")?;
/// let (bases, stats) = fqxv::decompress_records_par(
///     file,
///     0,
///     || 0u64,
///     |acc, _index, rec| {
///         *acc += rec.seq.len() as u64;
///         Ok(())
///     },
///     |a, b| a + b,
/// )?;
/// println!("{bases} bases in {} reads", stats.reads);
/// # Ok(()) }
/// ```
pub fn decompress_records_par<R, S, I, V, M>(
    reader: R,
    threads: usize,
    init: I,
    visit: V,
    merge: M,
) -> Result<(S, Stats)>
where
    R: Read,
    S: Send,
    I: Fn() -> S + Sync,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()> + Sync,
    M: FnMut(S, S) -> S,
{
    let (state, mut stats) =
        decompress_records_par_select(reader, threads, StreamSelection::ALL, init, visit, merge)?;
    // The selective entry point counts emitted field bytes (its `out_bytes`
    // contract, matching `decompress_records_select`); this API documents the
    // FASTQ byte count instead. With everything selected the two differ by
    // exactly the per-record framing — `@`, the `+` line, and three newlines: 6
    // bytes — since `qual` and `seq` are equal length, so
    // `6 + name + 2*seq == 6 + name + seq + qual` per record.
    stats.out_bytes += 6 * stats.reads;
    Ok((state, stats))
}

/// [`decompress_records_par`] with a [`StreamSelection`]: the parallel
/// map-reduce visitor over only the streams the caller asks for. Deselected
/// [`RecordRef`] fields come back **empty** and their coded streams are
/// genuinely skipped *inside* the per-block decode tasks — seeked past, never
/// entropy-decoded — the same contract as the serial
/// [`decompress_records_select`](super::decompress_records_select). The
/// per-stream savings stack with the visitor's block-parallel scaling: a
/// sequence-only pass on a long-read (ONT) archive measured ~12x at any thread
/// count (quality dominates long-read decode), and short-read Illumina saves a
/// further ~1.3-1.6x of decode CPU on top of the visitor's scaling (see the
/// [`stream_select`](super::stream_select) module docs for the measured
/// baselines).
///
/// The hooks, index/ordering contract, per-layout delivery, and error/panic
/// behavior are exactly [`decompress_records_par`]'s; with
/// [`StreamSelection::ALL`] the visited records are identical to that
/// function's. Two things change under a *partial* selection:
///
/// - **Integrity**: every stream that *is* decoded keeps its content-digest
///   check (per-block stream digests on the plain layout) and all frames read
///   keep their CRCs, but a skipped stream's digest cannot be checked, and the
///   whole-file reorder layout's joint output digest — which spans all three
///   streams — is skipped for any partial selection, exactly as on the serial
///   selective path.
/// - **[`Stats::out_bytes`]** counts the bytes of the emitted (selected)
///   fields, matching [`decompress_records_select`](super::decompress_records_select)
///   — not the serialized-FASTQ byte count [`decompress_records_par`] reports.
///
/// The long-read cross-stream dependency carries over: quality coded against
/// the bases decodes the sequence stream internally when only quality is
/// selected (still delivered empty unless selected), and a keep-order reorder
/// archive still decodes the (tiny) permutation and flip bitmap to restore its
/// *sequences* to original order.
///
/// ```no_run
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let file = std::fs::File::open("reads.fqxv")?;
/// // Parallel k-mer-counter shape: sequence only, quality and names skipped.
/// let (bases, stats) = fqxv::decompress_records_par_select(
///     file,
///     0,
///     fqxv::StreamSelection::SEQUENCE_ONLY,
///     || 0u64,
///     |acc, _index, rec| {
///         *acc += rec.seq.len() as u64; // rec.name and rec.qual are empty
///         Ok(())
///     },
///     |a, b| a + b,
/// )?;
/// println!("{bases} bases in {} reads", stats.reads);
/// # Ok(()) }
/// ```
pub fn decompress_records_par_select<R, S, I, V, M>(
    reader: R,
    threads: usize,
    selection: StreamSelection,
    init: I,
    visit: V,
    merge: M,
) -> Result<(S, Stats)>
where
    R: Read,
    S: Send,
    I: Fn() -> S + Sync,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()> + Sync,
    M: FnMut(S, S) -> S,
{
    let pool = build_pool(threads)?;
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    let states = StatePool::new(&init);
    let stats = if header.flags & FLAG_GLOBAL_REORDER != 0 {
        par_visit_reordered(r, &pool, &header, selection, &states, &visit)?
    } else {
        par_visit_plain(r, &pool, &header, selection, &states, &visit)?
    };
    Ok((states.into_merged(merge), stats))
}

/// The `n_reads` field of a raw block payload (`[3 x 8 digests][4 n_reads]…`),
/// read without decoding: the per-block record index bases are needed *before*
/// the batch fans out. [`decode_block_parts`] reads the same field and
/// cross-checks it against every decoded stream's length, so a value that lies
/// cannot survive to a successful return.
fn peek_block_reads(payload: &[u8]) -> Result<u32> {
    let bytes = payload
        .get(STREAM_DIGESTS_LEN..STREAM_DIGESTS_LEN + 4)
        .ok_or(Error::Truncated)?;
    Ok(u32::from_le_bytes(bytes.try_into().unwrap()))
}

/// Plain-layout parallel visit: read block batches serially (as [`decompress`]
/// does), decode the batch's blocks in parallel under `sel`
/// ([`decode_block_parts_select`], which skips deselected streams), and visit
/// each block's records inside its decode task. `r` is positioned just past the
/// header.
fn par_visit_plain<R, S, I, V>(
    mut r: BufReader<R>,
    pool: &rayon::ThreadPool,
    header: &Header,
    sel: StreamSelection,
    states: &StatePool<'_, S, I>,
    visit: &V,
) -> Result<Stats>
where
    R: Read,
    S: Send,
    I: Fn() -> S + Sync,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()> + Sync,
{
    // The whole-file shared reference frame feeds only the sequence decoder
    // (and, transitively, sequence-conditioned long-read quality), so a
    // names-only or empty selection skips its decode too.
    let reference = if header.flags & FLAG_GLOBAL_REFERENCE != 0 && !(sel.sequence || sel.quality) {
        skip_framed(&mut r)?;
        None
    } else {
        read_reference_frame(&mut r, header.flags)?
    };
    let reference = reference.as_ref();
    let batch = pool.current_num_threads().max(1);
    let mut stats = Stats {
        group_size: header.group_size.max(1),
        ..Stats::default()
    };
    let mut base = 0u64;
    for_each_block_batch(&mut r, batch, |raw_blocks| {
        // Record index base of each block in the batch, from the payloads'
        // n_reads fields (verified against the decode inside the task).
        let mut bases = Vec::with_capacity(raw_blocks.len());
        for b in raw_blocks {
            bases.push(base);
            base += u64::from(peek_block_reads(b)?);
        }
        let decoded: Vec<Result<(u64, u64)>> = pool.install(|| {
            (0..raw_blocks.len())
                .into_par_iter()
                .map(|k| visit_block(&raw_blocks[k], bases[k], reference, sel, states, visit))
                .collect()
        });
        for d in decoded {
            let (reads, bytes) = d?;
            stats.reads += reads;
            stats.blocks += 1;
            stats.out_bytes += bytes;
        }
        Ok(())
    })?;
    Ok(stats)
}

/// Decode one plain-layout block under `sel` and visit its records in index
/// order with a checked-out state; deselected [`RecordRef`] fields are empty.
/// Returns `(reads, emitted field bytes)`.
fn visit_block<S, I, V>(
    payload: &[u8],
    base: u64,
    reference: Option<&fqxv_lroverlap::Reference>,
    sel: StreamSelection,
    states: &StatePool<'_, S, I>,
    visit: &V,
) -> Result<(u64, u64)>
where
    I: Fn() -> S,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()>,
{
    let (n_reads, names, lens, seq, qual) = decode_block_parts_select(payload, reference, sel)?;
    let mut st = states.take();
    let mut off = 0usize;
    let mut bytes = 0u64;
    for i in 0..n_reads {
        // Per-read lengths exist whenever sequence or quality was decoded (they
        // share them); with neither selected every slice is empty.
        let l = if sel.sequence || sel.quality {
            lens[i] as usize
        } else {
            0
        };
        let end = off
            .checked_add(l)
            .ok_or(Error::Malformed("read length overflow"))?;
        let s: &[u8] = if sel.sequence {
            seq.get(off..end).ok_or(Error::Malformed(
                "sequence shorter than declared read lengths",
            ))?
        } else {
            &[]
        };
        let q: &[u8] = if sel.quality {
            qual.get(off..end).ok_or(Error::Malformed(
                "quality shorter than declared read lengths",
            ))?
        } else {
            &[]
        };
        let name: &[u8] = if sel.names { &names[i] } else { &[] };
        bytes += (name.len() + s.len() + q.len()) as u64;
        visit(
            &mut st,
            base + i as u64,
            RecordRef {
                name,
                seq: s,
                qual: q,
            },
        )?;
        off = end;
    }
    states.put(st);
    Ok((n_reads as u64, bytes))
}

/// Whole-file reorder parallel visit: decode the archive into output-order
/// records ([`decode_reordered_records_select`] — block-parallel; under a full
/// selection the joint output digest is verified first, under a partial one it
/// is unverifiable and skipped), then visit contiguous index ranges in
/// parallel. `r` is positioned just past the header.
fn par_visit_reordered<R, S, I, V>(
    r: BufReader<R>,
    pool: &rayon::ThreadPool,
    header: &Header,
    sel: StreamSelection,
    states: &StatePool<'_, S, I>,
    visit: &V,
) -> Result<Stats>
where
    R: Read,
    S: Send,
    I: Fn() -> S + Sync,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()> + Sync,
{
    let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
    let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
    let rr = decode_reordered_records_select(r, pool, keep_order, has_reference, sel)?;
    let n = rr.n;
    let chunk = n
        .div_ceil(pool.current_num_threads().max(1) * 8)
        .max(MIN_VISIT_CHUNK);
    let ranges: Vec<(usize, usize)> = (0..n)
        .step_by(chunk)
        .map(|s| (s, (s + chunk).min(n)))
        .collect();
    let out_bytes = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(a, b)| -> Result<u64> {
                let mut st = states.take();
                let mut bytes = 0u64;
                for i in a..b {
                    // A deselected stream's vectors are empty (not
                    // per-record-empty) — deliver empty slices for it.
                    let name: &[u8] = if rr.names.is_empty() {
                        &[]
                    } else {
                        &rr.names[i]
                    };
                    let seq: &[u8] = if rr.seqs.is_empty() { &[] } else { &rr.seqs[i] };
                    let qual: &[u8] = if rr.qoffs.is_empty() {
                        &[]
                    } else {
                        &rr.quals[rr.qoffs[i]..rr.qoffs[i + 1]]
                    };
                    bytes += (name.len() + seq.len() + qual.len()) as u64;
                    visit(&mut st, i as u64, RecordRef { name, seq, qual })?;
                }
                states.put(st);
                Ok(bytes)
            })
            .try_reduce(|| 0u64, |a, b| Ok(a + b))
    })?;
    Ok(Stats {
        reads: n as u64,
        blocks: rr.n_blocks as u64,
        out_bytes,
        group_size: header.group_size.max(1),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::records::{RecordReader, decompress_records};

    /// Thread counts every equality test runs at, to shake out scheduling- and
    /// partitioning-dependent assumptions.
    const THREADS: [usize; 4] = [1, 2, 4, 8];

    /// Reference sequence: the records the serial reader yields, in order.
    fn serial_records(archive: &[u8]) -> Vec<Record> {
        RecordReader::new(io::Cursor::new(archive.to_vec()), 1)
            .map(|r| r.expect("serial record"))
            .collect()
    }

    /// Collect `(index, record)` pairs from the parallel visitor.
    fn par_records(archive: &[u8], threads: usize) -> (Vec<(u64, Record)>, Stats) {
        let (mut pairs, stats) = decompress_records_par(
            archive,
            threads,
            Vec::new,
            |v: &mut Vec<(u64, Record)>, i, rec| {
                v.push((i, rec.to_record()));
                Ok(())
            },
            |mut a, mut b| {
                a.append(&mut b);
                a
            },
        )
        .expect("decompress_records_par");
        pairs.sort_by_key(|&(i, _)| i);
        (pairs, stats)
    }

    /// Assert the parallel visitor delivers exactly the serial sequence — every
    /// index once, each record equal to the serial reader's record at that
    /// index — at several thread counts.
    fn assert_par_matches_serial(archive: &[u8]) {
        let expected = serial_records(archive);
        for threads in THREADS {
            let (pairs, stats) = par_records(archive, threads);
            assert_eq!(
                pairs.len(),
                expected.len(),
                "record count at {threads} threads"
            );
            assert_eq!(stats.reads, expected.len() as u64, "stats.reads");
            for (k, (i, rec)) in pairs.iter().enumerate() {
                assert_eq!(*i, k as u64, "indices must be dense at {threads} threads");
                assert_eq!(rec, &expected[k], "record {k} at {threads} threads");
            }
        }
    }

    /// A plain single-end fixture large enough for several blocks.
    fn plain_input(n: u32) -> Vec<u8> {
        let mut input = Vec::new();
        for i in 0..n {
            let seq = match i % 4 {
                0 => "ACGTACGTACGTACGT",
                1 => "TTGGCCAATTGGCCAA",
                2 => "GATTACAGATTACANN",
                _ => "CCCCGGGGAAAATTTT",
            };
            input.extend_from_slice(
                format!("@read.{i} desc{}\n{seq}\n+\nIIIIFFFF####IIII\n", i % 7).as_bytes(),
            );
        }
        input
    }

    fn archive_plain(n: u32, block_reads: usize) -> Vec<u8> {
        let mut archive = Vec::new();
        compress(
            &plain_input(n)[..],
            &mut archive,
            Params {
                block_reads,
                seq_order: 4,
                ..Params::default()
            },
        )
        .expect("compress");
        archive
    }

    #[test]
    fn plain_matches_serial_across_blocks_and_threads() {
        assert_par_matches_serial(&archive_plain(3000, 64));
    }

    #[test]
    fn plain_single_block_matches_serial() {
        assert_par_matches_serial(&archive_plain(50, 1 << 16));
    }

    #[test]
    fn empty_archive_yields_init_state() {
        let mut archive = Vec::new();
        compress(&b""[..], &mut archive, Params::default()).expect("compress empty");
        let (state, stats) = decompress_records_par(
            &archive[..],
            2,
            || 41u64,
            |_s, _i, _r| panic!("no records to visit"),
            |a, _b| a,
        )
        .expect("empty archive");
        assert_eq!(state, 41, "state must be a fresh init()");
        assert_eq!(stats.reads, 0);
    }

    /// Paired-spot grouped archive (G=2); small blocks so several blocks each
    /// hold whole spots.
    fn archive_grouped(n_spots: u32) -> Vec<u8> {
        let mut input = Vec::new();
        for i in 0..n_spots {
            for m in 1..=2u32 {
                input.extend_from_slice(format!("@sp.{i}/{m}\nACGTTGCA\n+\nIIIIFFFF\n").as_bytes());
            }
        }
        let mut archive = Vec::new();
        compress_interleaved(
            &input[..],
            &mut archive,
            Params {
                block_reads: 64,
                seq_order: 4,
                ..Params::default()
            },
            2,
        )
        .expect("compress_interleaved");
        archive
    }

    #[test]
    fn grouped_interleaved_matches_serial_and_member_index() {
        let archive = archive_grouped(800);
        assert_par_matches_serial(&archive);
        // Member identity is index % G: the fixture encodes it in the name.
        let (pairs, stats) = par_records(&archive, 4);
        assert_eq!(stats.group_size, 2);
        for (i, rec) in &pairs {
            let member = i % 2 + 1;
            assert!(
                rec.name.ends_with(format!("/{member}").as_bytes()),
                "record {i} must be member {member}, name {:?}",
                String::from_utf8_lossy(&rec.name)
            );
        }
    }

    /// Short redundant reads that cluster well, with counter-style names.
    fn reorder_input(n: u32) -> Vec<u8> {
        let motifs = [
            "ACGTACGTACGTACGTACGTACGT",
            "TTTTGGGGCCCCAAAATTTTGGGG",
            "GATTACAGATTACAGATTACAGAT",
        ];
        let mut input = Vec::new();
        for i in 0..n {
            let seq = motifs[(i % 3) as usize];
            input.extend_from_slice(
                format!("@r.{i}\n{seq}\n+\nIIIIIIIIIIIIFFFFFFFF####\n").as_bytes(),
            );
        }
        input
    }

    fn archive_clustered(n: u32, keep_order: bool, regenerate_names: bool) -> Vec<u8> {
        let params = Params {
            reorder: true,
            keep_order,
            regenerate_names,
            seq_order: 4,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress_clustered_forced(&reorder_input(n)[..], &mut archive, params, 1)
            .expect("forced clustered compress");
        assert_ne!(
            archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
            0,
            "fixture must be the clustered layout"
        );
        archive
    }

    #[test]
    fn reorder_keep_order_matches_serial() {
        assert_par_matches_serial(&archive_clustered(2000, true, false));
    }

    #[test]
    fn reorder_discard_order_matches_serial() {
        // Discard-order emits clustered order; the serial reader defines that
        // order and the parallel indices must reproduce it exactly.
        assert_par_matches_serial(&archive_clustered(2000, false, false));
    }

    #[test]
    fn reorder_regenerated_names_matches_serial() {
        assert_par_matches_serial(&archive_clustered(2000, false, true));
    }

    /// 600 bp overlapping reads clear the long-read gate, so blocks run the
    /// long-read candidate contest (order-k / overlap / LZMA — whichever wins,
    /// the visitor rides the same decode dispatch).
    fn archive_long_read() -> Vec<u8> {
        let genome: Vec<u8> = (0..4000u32)
            .map(|i| b"ACGT"[((i.wrapping_mul(2_654_435_761) >> 13) & 3) as usize])
            .collect();
        let mut input = Vec::new();
        for i in 0..60u32 {
            let start = (i as usize * 53) % (genome.len() - 600);
            let mut s = genome[start..start + 600].to_vec();
            s[77] = b'N';
            let qual = vec![b'I'; s.len()];
            input.extend_from_slice(format!("@lr.{i}\n").as_bytes());
            input.extend_from_slice(&s);
            input.extend_from_slice(b"\n+\n");
            input.extend_from_slice(&qual);
            input.push(b'\n');
        }
        let mut archive = Vec::new();
        compress(
            &input[..],
            &mut archive,
            Params {
                block_reads: 16,
                ..Params::default()
            },
        )
        .expect("compress long reads");
        archive
    }

    #[test]
    fn long_read_blocks_match_serial() {
        assert_par_matches_serial(&archive_long_read());
    }

    /// Selections the parallel selective visitor is exercised at, per layout:
    /// the analysis-shaped ones, plus `QUALITY_ONLY` (which on long-read
    /// archives exercises the sequence-conditioned quality dependency).
    const PAR_SELECTIONS: [StreamSelection; 4] = [
        StreamSelection::SEQUENCE_ONLY,
        StreamSelection::NAMES_AND_SEQUENCE,
        StreamSelection::NAMES_ONLY,
        StreamSelection::QUALITY_ONLY,
    ];

    /// Reference records for a selection: the serial selective decode.
    fn serial_select_records(archive: &[u8], sel: StreamSelection) -> Vec<Record> {
        let mut recs = Vec::new();
        decompress_records_select(archive, 1, sel, |r| recs.push(r)).expect("serial select");
        recs
    }

    /// Collect `(index, record)` pairs from the parallel selective visitor,
    /// sorted by index.
    fn par_select_records(
        archive: &[u8],
        threads: usize,
        sel: StreamSelection,
    ) -> (Vec<(u64, Record)>, Stats) {
        let (mut pairs, stats) = decompress_records_par_select(
            archive,
            threads,
            sel,
            Vec::new,
            |v: &mut Vec<(u64, Record)>, i, rec| {
                v.push((i, rec.to_record()));
                Ok(())
            },
            |mut a, mut b| {
                a.append(&mut b);
                a
            },
        )
        .expect("decompress_records_par_select");
        pairs.sort_by_key(|&(i, _)| i);
        (pairs, stats)
    }

    /// Assert the parallel selective visitor delivers exactly the serial
    /// selective decode's records — every index once, equal records (deselected
    /// fields empty on both sides), `out_bytes` counting the emitted fields —
    /// for every selection and thread count.
    fn assert_par_select_matches_serial(archive: &[u8], what: &str) {
        for sel in PAR_SELECTIONS {
            let expected = serial_select_records(archive, sel);
            let field_bytes: u64 = expected
                .iter()
                .map(|r| (r.name.len() + r.seq.len() + r.qual.len()) as u64)
                .sum();
            for threads in THREADS {
                let (pairs, stats) = par_select_records(archive, threads, sel);
                assert_eq!(
                    pairs.len(),
                    expected.len(),
                    "{what} {sel:?}: record count at {threads} threads"
                );
                assert_eq!(stats.reads, expected.len() as u64, "{what} {sel:?}: reads");
                assert_eq!(
                    stats.out_bytes, field_bytes,
                    "{what} {sel:?}: out_bytes must count the emitted fields"
                );
                for (k, (i, rec)) in pairs.iter().enumerate() {
                    assert_eq!(
                        *i, k as u64,
                        "{what} {sel:?}: indices must be dense at {threads} threads"
                    );
                    assert_eq!(
                        rec, &expected[k],
                        "{what} {sel:?}: record {k} at {threads} threads"
                    );
                }
            }
        }
    }

    #[test]
    fn select_plain_matches_serial() {
        assert_par_select_matches_serial(&archive_plain(3000, 64), "plain multi-block");
    }

    #[test]
    fn select_grouped_matches_serial() {
        assert_par_select_matches_serial(&archive_grouped(800), "grouped G=2");
    }

    #[test]
    fn select_reorder_keep_order_matches_serial() {
        assert_par_select_matches_serial(
            &archive_clustered(2000, true, false),
            "reorder keep-order",
        );
    }

    #[test]
    fn select_reorder_discard_order_matches_serial() {
        assert_par_select_matches_serial(
            &archive_clustered(2000, false, false),
            "reorder discard-order",
        );
    }

    #[test]
    fn select_reorder_regenerated_names_matches_serial() {
        assert_par_select_matches_serial(
            &archive_clustered(2000, false, true),
            "reorder regen-names",
        );
    }

    #[test]
    fn select_long_read_blocks_match_serial() {
        assert_par_select_matches_serial(&archive_long_read(), "long-read plain");
    }

    #[test]
    fn full_selection_is_identical_to_decompress_records_par() {
        // `StreamSelection::ALL` must ride the exact full-decode path
        // (`decompress_records_par` is now a wrapper over it): same records,
        // same indices, same stats — with `out_bytes` differing only by the
        // documented per-record FASTQ framing (6 bytes per record).
        for archive in [
            archive_plain(1000, 64),
            archive_grouped(400),
            archive_clustered(1000, true, false),
            archive_clustered(1000, false, true),
        ] {
            let (all_pairs, all_stats) = par_select_records(&archive, 4, StreamSelection::ALL);
            let (par_pairs, par_stats) = par_records(&archive, 4);
            assert_eq!(all_pairs, par_pairs, "ALL must deliver identical records");
            assert_eq!(all_stats.reads, par_stats.reads);
            assert_eq!(all_stats.blocks, par_stats.blocks);
            assert_eq!(all_stats.group_size, par_stats.group_size);
            assert_eq!(
                par_stats.out_bytes,
                all_stats.out_bytes + 6 * all_stats.reads,
                "FASTQ byte count = field bytes + 6 per record"
            );
        }
    }

    #[test]
    fn none_selection_visits_every_index_with_empty_fields() {
        for (archive, want_reads) in [
            (archive_plain(500, 64), 500u64),
            (archive_clustered(1000, false, false), 1000),
        ] {
            let (pairs, stats) = par_select_records(&archive, 4, StreamSelection::NONE);
            assert_eq!(stats.reads, want_reads);
            assert_eq!(pairs.len() as u64, want_reads);
            assert_eq!(stats.out_bytes, 0);
            for (k, (i, rec)) in pairs.iter().enumerate() {
                assert_eq!(*i, k as u64, "indices must be dense");
                assert!(
                    rec.name.is_empty() && rec.seq.is_empty() && rec.qual.is_empty(),
                    "NONE must deliver empty fields"
                );
            }
        }
    }

    /// Order-independent reduction: (records, bases, wrapping checksum) must
    /// match the same reduction computed serially, at every thread count.
    #[test]
    fn order_free_reduction_matches_serial() {
        let archive = archive_plain(3000, 64);
        let fold = |acc: &mut (u64, u64, u64), name: &[u8], seq: &[u8], qual: &[u8]| {
            let mut h = Xxh3::new();
            h.update(name);
            h.update(seq);
            h.update(qual);
            acc.0 += 1;
            acc.1 += seq.len() as u64;
            acc.2 = acc.2.wrapping_add(h.digest());
        };
        let mut want = (0u64, 0u64, 0u64);
        decompress_records(&archive[..], 1, |r| {
            fold(&mut want, &r.name, &r.seq, &r.qual)
        })
        .expect("serial");
        for threads in THREADS {
            let (got, stats) = decompress_records_par(
                &archive[..],
                threads,
                || (0u64, 0u64, 0u64),
                |acc, _i, rec| {
                    fold(acc, rec.name, rec.seq, rec.qual);
                    Ok(())
                },
                |a, b| (a.0 + b.0, a.1 + b.1, a.2.wrapping_add(b.2)),
            )
            .expect("parallel");
            assert_eq!(got, want, "reduction at {threads} threads");
            assert_eq!(stats.reads, want.0);
        }
    }

    #[test]
    fn stats_match_decompress_output_bytes() {
        for archive in [
            archive_plain(500, 64),
            archive_clustered(1000, true, false),
            archive_clustered(1000, false, false),
        ] {
            let mut fastq = Vec::new();
            decompress(&archive[..], &mut fastq, 1).expect("decompress");
            let (_, stats) =
                decompress_records_par(&archive[..], 4, || (), |_, _, _| Ok(()), |a, _| a)
                    .expect("par");
            assert_eq!(
                stats.out_bytes,
                fastq.len() as u64,
                "out_bytes must equal the serialized FASTQ size"
            );
        }
    }

    #[test]
    fn visitor_error_propagates() {
        let archive = archive_plain(3000, 64);
        for threads in [1, 4] {
            let err = decompress_records_par(
                &archive[..],
                threads,
                || (),
                |_s, i, _rec| {
                    if i == 1234 {
                        Err(Error::from(io::Error::other("visitor abort")))
                    } else {
                        Ok(())
                    }
                },
                |a, _| a,
            );
            assert!(
                matches!(err, Err(Error::Io(_))),
                "visitor error must propagate at {threads} threads, got {err:?}"
            );
        }
    }

    #[test]
    fn visitor_panic_propagates() {
        let archive = archive_plain(1000, 64);
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = decompress_records_par(
                &archive[..],
                4,
                || (),
                |_s, i, _rec| {
                    if i == 500 {
                        panic!("visitor panic");
                    }
                    Ok(())
                },
                |a, _| a,
            );
        }));
        assert!(panicked.is_err(), "visitor panic must propagate");
    }

    #[test]
    fn corrupted_archive_errors() {
        let mut archive = archive_plain(500, 64);
        let mid = archive.len() / 2;
        archive[mid] ^= 0xFF;
        let res = decompress_records_par(&archive[..], 4, || (), |_, _, _| Ok(()), |a, _| a);
        assert!(res.is_err(), "corruption must surface as an error");
    }
}
