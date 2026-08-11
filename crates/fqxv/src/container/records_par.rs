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
    let pool = build_pool(threads)?;
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    let states = StatePool::new(&init);
    let stats = if header.flags & FLAG_GLOBAL_REORDER != 0 {
        let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
        let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
        par_visit_reordered(
            r,
            &pool,
            keep_order,
            header.group_size,
            has_reference,
            &states,
            &visit,
        )?
    } else {
        par_visit_plain(r, &pool, &header, &states, &visit)?
    };
    Ok((states.into_merged(merge), stats))
}

/// The serialized-FASTQ byte count of one record: `@` + name + newline, the
/// sequence line, the normalized `+` line, and the quality line.
fn fastq_bytes(name_len: usize, seq_len: usize) -> u64 {
    (6 + name_len + 2 * seq_len) as u64
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
/// does), decode the batch's blocks in parallel, and visit each block's records
/// inside its decode task. `r` is positioned just past the header.
fn par_visit_plain<R, S, I, V>(
    mut r: BufReader<R>,
    pool: &rayon::ThreadPool,
    header: &Header,
    states: &StatePool<'_, S, I>,
    visit: &V,
) -> Result<Stats>
where
    R: Read,
    S: Send,
    I: Fn() -> S + Sync,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()> + Sync,
{
    let reference = read_reference_frame(&mut r, header.flags)?;
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
                .map(|k| visit_block(&raw_blocks[k], bases[k], reference, states, visit))
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

/// Decode one plain-layout block and visit its records in index order with a
/// checked-out state. Returns `(reads, serialized FASTQ bytes)`.
fn visit_block<S, I, V>(
    payload: &[u8],
    base: u64,
    reference: Option<&fqxv_lroverlap::Reference>,
    states: &StatePool<'_, S, I>,
    visit: &V,
) -> Result<(u64, u64)>
where
    I: Fn() -> S,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()>,
{
    let (n_reads, names, lens, seq, qual) = decode_block_parts(payload, reference)?;
    let mut st = states.take();
    let mut off = 0usize;
    let mut bytes = 0u64;
    for i in 0..n_reads {
        let l = lens[i] as usize;
        let (s, q) = read_slices(&seq, &qual, off, l)?;
        bytes += fastq_bytes(names[i].len(), l);
        visit(
            &mut st,
            base + i as u64,
            RecordRef {
                name: &names[i],
                seq: s,
                qual: q,
            },
        )?;
        off += l;
    }
    states.put(st);
    Ok((n_reads as u64, bytes))
}

/// Whole-file reorder parallel visit: decode the archive into output-order
/// records ([`decode_reordered_records`], block-parallel and digest-verified),
/// then visit contiguous index ranges in parallel. `r` is positioned just past
/// the header.
fn par_visit_reordered<R, S, I, V>(
    r: BufReader<R>,
    pool: &rayon::ThreadPool,
    keep_order: bool,
    group_size: u8,
    has_reference: bool,
    states: &StatePool<'_, S, I>,
    visit: &V,
) -> Result<Stats>
where
    R: Read,
    S: Send,
    I: Fn() -> S + Sync,
    V: Fn(&mut S, u64, RecordRef<'_>) -> Result<()> + Sync,
{
    let rr = decode_reordered_records(r, pool, keep_order, has_reference)?;
    let n = rr.names.len();
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
                    let qual = &rr.quals[rr.qoffs[i]..rr.qoffs[i + 1]];
                    bytes += fastq_bytes(rr.names[i].len(), rr.seqs[i].len());
                    visit(
                        &mut st,
                        i as u64,
                        RecordRef {
                            name: &rr.names[i],
                            seq: &rr.seqs[i],
                            qual,
                        },
                    )?;
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
        group_size: group_size.max(1),
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

    #[test]
    fn grouped_interleaved_matches_serial_and_member_index() {
        // Paired spots; small blocks so several blocks each hold whole spots.
        let mut input = Vec::new();
        for i in 0..800u32 {
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

    #[test]
    fn long_read_blocks_match_serial() {
        // 600 bp overlapping reads clear the long-read gate, so blocks run the
        // long-read candidate contest (order-k / overlap / LZMA — whichever
        // wins, the visitor rides the same decode dispatch).
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
        assert_par_matches_serial(&archive);
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
