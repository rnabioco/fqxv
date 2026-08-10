//! Stream-selective decode: reconstruct only the per-read streams a caller asks
//! for, *skipping* — never entropy-decoding — the rest.
//!
//! Every layout stores names, sequence, and quality as independently framed
//! coded streams (the plain layout length-prefixes them inside each block
//! payload; the whole-file reorder layout frames each one separately), so a
//! deselected stream can be seeked past without running its codec. The range
//! coders are ~all of decode compute, so skipping a stream saves its full
//! decode cost. How much that buys depends on where the compute sits: on
//! long-read archives the sequence-conditioned quality model dominates and a
//! sequence-only decode measured **~12x** faster (ONT, 1 and 16 threads); on
//! short-read Illumina the sequence context model shares the cost, for a
//! ~1.3-1.6x single-threaded win (NovaSeq/MiSeq) — either way, fast enough to
//! use an archive directly as k-mer-counting or alignment input.
//!
//! What can be skipped, per layout:
//!
//! - **Plain (and grouped `G > 1`)**: every stream decodes independently — with
//!   one exception: long-read quality is coded *against the bases*
//!   ([`quality_needs_sequence`](super::random_access::quality_needs_sequence)),
//!   so selecting quality there also decodes the sequence stream internally
//!   (still returned empty unless selected). The whole-file reference frame is
//!   only decoded when sequence or quality is selected (a names-only pass skips
//!   it).
//! - **Reorder, order preserved (`--order any`)**: names and quality are coded
//!   in original order and decode independently; the sequence partition needs
//!   the permutation and flip bitmap (both tiny) to restore original order.
//! - **Reorder, order discarded (`--order shuffle`)**: reads emerge in
//!   clustered order; each stream decodes independently (quality additionally
//!   reads the flip bitmap to un-reverse flipped reads; regenerated names come
//!   from the stored template, near-free).
//!
//! Integrity under selection: each *decoded* stream is still verified against
//! its stored per-stream content digest (plain layout) and the frame CRCs of
//! everything read are still checked. A skipped stream's digest cannot be
//! verified (its content is never reconstructed); in the plain layout its
//! stored bytes are still covered by the block frame CRC, while in the reorder
//! layout a skipped frame's CRC is not checked at all. The reorder layout's
//! whole-output digest spans all three streams jointly, so any partial
//! selection skips that check too — a full decode ([`StreamSelection::ALL`])
//! keeps every check, exactly as [`decompress`](super::decompress).

use super::records::RecordSink;
use super::*;
use rayon::prelude::*;

/// Which of the three per-read streams to decode. Deselected streams are
/// skipped (their coded bytes are never entropy-decoded) and the matching
/// [`Record`] fields come back empty.
///
/// See the [module docs](self) for what skipping means per archive layout and
/// which integrity checks a partial selection retains.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StreamSelection {
    /// Decode read names and descriptions.
    pub names: bool,
    /// Decode sequence bases.
    pub sequence: bool,
    /// Decode quality scores. On long-read archives the quality stream is coded
    /// against the bases, so selecting quality also costs the sequence decode
    /// (the sequence still comes back empty unless selected).
    pub quality: bool,
}

impl StreamSelection {
    /// Decode everything — equivalent to a full decode, with every integrity
    /// check intact.
    pub const ALL: StreamSelection = StreamSelection {
        names: true,
        sequence: true,
        quality: true,
    };
    /// Decode nothing: every record comes back with all fields empty. Useful to
    /// count reads at framing speed (no stream is entropy-decoded).
    pub const NONE: StreamSelection = StreamSelection {
        names: false,
        sequence: false,
        quality: false,
    };
    /// Sequence bases only — k-mer counting / classification input.
    pub const SEQUENCE_ONLY: StreamSelection = StreamSelection {
        names: false,
        sequence: true,
        quality: false,
    };
    /// Read names only.
    pub const NAMES_ONLY: StreamSelection = StreamSelection {
        names: true,
        sequence: false,
        quality: false,
    };
    /// Quality scores only.
    pub const QUALITY_ONLY: StreamSelection = StreamSelection {
        names: false,
        sequence: false,
        quality: true,
    };
    /// Names + sequence — what FASTA output needs ([`decompress_fasta`]).
    pub const NAMES_AND_SEQUENCE: StreamSelection = StreamSelection {
        names: true,
        sequence: true,
        quality: false,
    };
}

impl Default for StreamSelection {
    /// Defaults to [`StreamSelection::ALL`].
    fn default() -> Self {
        StreamSelection::ALL
    }
}

/// Decode an archive invoking `on_record` per record — [`decompress_records`]
/// with a [`StreamSelection`]: deselected [`Record`] fields come back empty and
/// their coded streams are genuinely skipped, not decoded-and-discarded, so a
/// sequence-only pass saves the skipped streams' full decode cost (~12x
/// measured on ONT, where quality dominates; ~1.3-1.6x single-threaded on
/// short-read Illumina — see the [module docs](self)).
///
/// Records arrive in the same order [`decompress`] would emit them (original
/// file order; clustered order for an `--order shuffle` archive), interleaved
/// for grouped archives. With [`StreamSelection::ALL`] this is exactly
/// [`decompress_records`]. The returned [`Stats::out_bytes`] counts the bytes
/// of the emitted (selected) fields.
///
/// ```no_run
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let file = std::fs::File::open("reads.fqxv")?;
/// let mut bases = 0u64;
/// fqxv::decompress_records_select(file, 0, fqxv::StreamSelection::SEQUENCE_ONLY, |rec| {
///     bases += rec.seq.len() as u64; // rec.name and rec.qual are empty
/// })?;
/// println!("{bases} bases");
/// # Ok(()) }
/// ```
pub fn decompress_records_select<R: Read>(
    reader: R,
    threads: usize,
    selection: StreamSelection,
    mut on_record: impl FnMut(Record),
) -> Result<Stats> {
    decompress_select(reader, threads, selection, move |rec| {
        on_record(rec);
        Ok(())
    })
}

/// Decompress an archive to single-line FASTA (`>name` + sequence) on `writer`,
/// skipping the quality stream entirely — it is never entropy-decoded, which
/// makes this faster than [`decompress`] by the quality stream's full decode
/// cost: ~12x measured on long-read (ONT) archives, where the
/// sequence-conditioned quality model dominates decode, and ~1.3-1.6x
/// single-threaded on short-read Illumina (see the [module docs](self)).
/// Grouped archives emit interleaved records, like [`decompress`].
///
/// One caveat: long-read archives whose quality would need the sequence are
/// unaffected (quality is simply skipped), but the sequence stream itself is
/// always decoded — this is a *quality*-skipping fast path, not a projection of
/// an unrelated stream. [`Stats::out_bytes`] counts the FASTA bytes written.
///
/// ```no_run
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let file = std::fs::File::open("reads.fqxv")?;
/// let out = std::fs::File::create("reads.fasta")?;
/// fqxv::decompress_fasta(file, out, 0)?;
/// # Ok(()) }
/// ```
pub fn decompress_fasta<R: Read, W: Write>(reader: R, writer: W, threads: usize) -> Result<Stats> {
    let mut w = BufWriter::new(writer);
    let mut fasta_bytes = 0u64;
    let mut stats = decompress_select(
        reader,
        threads,
        StreamSelection::NAMES_AND_SEQUENCE,
        |rec| {
            fasta_bytes += (rec.name.len() + rec.seq.len() + 3) as u64;
            w.write_all(b">")?;
            w.write_all(&rec.name)?;
            w.write_all(b"\n")?;
            w.write_all(&rec.seq)?;
            w.write_all(b"\n")
        },
    )?;
    w.flush()?;
    stats.out_bytes = fasta_bytes;
    Ok(stats)
}

/// Shared driver behind every stream-selective entry point: dispatch on the
/// archive layout and push [`Record`]s (deselected fields empty) into `emit`.
/// A full selection delegates to the canonical [`decompress`] path so every
/// integrity check (including the reorder layout's whole-output digest) stays
/// exactly as strong there. `emit` errors abort the decode as [`Error::Io`].
pub(crate) fn decompress_select<R: Read>(
    reader: R,
    threads: usize,
    sel: StreamSelection,
    mut emit: impl FnMut(Record) -> io::Result<()>,
) -> Result<Stats> {
    if sel == StreamSelection::ALL {
        // Full selection: ride the canonical full-decode path (all layouts, all
        // digest checks), reassembling records from its text output.
        let mut field_bytes = 0u64;
        let sink = RecordSink::new(|rec: Record| {
            field_bytes += (rec.name.len() + rec.seq.len() + rec.qual.len()) as u64;
            emit(rec)
        });
        let mut stats = decompress(reader, sink, threads)?;
        stats.out_bytes = field_bytes;
        return Ok(stats);
    }

    let pool = build_pool(threads)?;
    let batch = pool.current_num_threads().max(1);
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;

    let mut out_bytes = 0u64;
    let mut stats;
    {
        let mut emit2 = |rec: Record| -> Result<()> {
            out_bytes += (rec.name.len() + rec.seq.len() + rec.qual.len()) as u64;
            emit(rec).map_err(Error::Io)
        };
        if header.flags & FLAG_GLOBAL_REORDER != 0 {
            let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
            let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
            stats = decode_reordered_select(
                &mut r,
                &pool,
                keep_order,
                has_reference,
                header.group_size,
                sel,
                &mut emit2,
            )?;
        } else {
            // Whole-file shared reference frame (plain layout): only the
            // sequence decoder (and, transitively, sequence-conditioned
            // quality) indexes it, so a names-only pass skips its decode too.
            let reference =
                if header.flags & FLAG_GLOBAL_REFERENCE != 0 && !(sel.sequence || sel.quality) {
                    skip_framed(&mut r)?;
                    None
                } else {
                    read_reference_frame(&mut r, header.flags)?
                };
            let reference = reference.as_ref();
            stats = Stats::default();
            for_each_block_batch(&mut r, batch, |raw_blocks| {
                let decoded: Vec<Result<(u64, Vec<Record>)>> = pool.install(|| {
                    raw_blocks
                        .par_iter()
                        .map(|b| decode_block_records_select(b, reference, sel))
                        .collect()
                });
                for d in decoded {
                    let (reads, recs) = d?;
                    for rec in recs {
                        emit2(rec)?;
                    }
                    stats.reads += reads;
                    stats.blocks += 1;
                }
                Ok(())
            })?;
        }
    }
    stats.out_bytes = out_bytes;
    Ok(stats)
}

/// Decode one plain-layout block into per-read [`Record`]s under `sel`
/// (deselected fields empty). Wraps [`decode_block_parts_select`] with
/// bounds-checked per-read slicing.
fn decode_block_records_select(
    buf: &[u8],
    shared_ref: Option<&fqxv_lroverlap::Reference>,
    sel: StreamSelection,
) -> Result<(u64, Vec<Record>)> {
    let (n_reads, mut names, lens, seq, qual) = decode_block_parts_select(buf, shared_ref, sel)?;
    let mut recs = Vec::with_capacity(n_reads);
    let mut off = 0usize;
    for i in 0..n_reads {
        // Per-read length: present whenever sequence or quality was decoded
        // (validated to have exactly n_reads entries); irrelevant otherwise.
        let l = if sel.sequence || sel.quality {
            lens[i] as usize
        } else {
            0
        };
        let end = off
            .checked_add(l)
            .ok_or(Error::Malformed("read length overflow"))?;
        let s = if sel.sequence {
            seq.get(off..end)
                .ok_or(Error::Malformed(
                    "sequence shorter than declared read lengths",
                ))?
                .to_vec()
        } else {
            Vec::new()
        };
        let q = if sel.quality {
            qual.get(off..end)
                .ok_or(Error::Malformed(
                    "quality shorter than declared read lengths",
                ))?
                .to_vec()
        } else {
            Vec::new()
        };
        let name = if sel.names {
            std::mem::take(&mut names[i])
        } else {
            Vec::new()
        };
        recs.push(Record {
            name,
            seq: s,
            qual: q,
        });
        off = end;
    }
    Ok((n_reads as u64, recs))
}

/// Stream-selective decode of the whole-file globally-clustered reorder layout.
/// `r` is positioned just past the header. Frames a deselected stream would
/// feed are skipped via [`skip_framed`]; in particular the permutation is only
/// decoded when a `keep_order` archive needs its *sequences* un-permuted
/// (names/quality are stored in original order there), and the shared global
/// reference only when sequences are selected. The whole-output digest spans
/// all three streams jointly, so it cannot be verified under a partial
/// selection (full selections never reach here — see [`decompress_select`]).
fn decode_reordered_select<R: Read>(
    mut r: R,
    pool: &rayon::ThreadPool,
    keep_order: bool,
    has_reference: bool,
    group_size: u8,
    sel: StreamSelection,
    emit: &mut dyn FnMut(Record) -> Result<()>,
) -> Result<Stats> {
    let mut n_buf = [0u8; 8];
    r.read_exact(&mut n_buf)?;
    let n = u64::from_le_bytes(n_buf) as usize;

    // Flip bitmap (tiny): sequences un-flip with it in both modes; clustered
    // quality was stored reversed for flipped reads and un-reverses with it too.
    let need_flip = sel.sequence || (!keep_order && sel.quality);
    let flip = if need_flip {
        read_framed(&mut r, "reorder flip bitmap")?
    } else {
        skip_framed(&mut r)?;
        Vec::new()
    };
    // Permutation: restores original positions of the clustered sequences in a
    // keep-order archive. Names/quality are coded in original order there and
    // never need it; discard-order archives store none (empty frame).
    let need_perm = keep_order && sel.sequence;
    let perm_c = if need_perm {
        read_framed(&mut r, "reorder permutation")?
    } else {
        skip_framed(&mut r)?;
        Vec::new()
    };
    let tmpl_bytes = read_framed(&mut r, "reorder name template")?;
    let template = if tmpl_bytes.is_empty() {
        None
    } else {
        Some(fqxv_tokenizer::NameTemplate::from_bytes(&tmpl_bytes)?)
    };
    let regen = template.is_some();
    // Shared global reference: only the version-4 sequence blocks index it.
    let reference = if has_reference {
        if sel.sequence {
            let ref_bytes = read_framed(&mut r, "reorder global reference")?;
            Some(decode_reference_frame(&ref_bytes)?)
        } else {
            skip_framed(&mut r)?;
            None
        }
    } else {
        None
    };
    let mut nb = [0u8; 4];
    r.read_exact(&mut nb)?;
    let n_blocks = u32::from_le_bytes(nb) as usize;

    let mut seq_payloads: Vec<Vec<u8>> = Vec::with_capacity(n_blocks.min(1 << 20));
    for i in 0..n_blocks {
        if sel.sequence {
            seq_payloads.push(read_framed(&mut r, &format!("reorder sequence block {i}"))?);
        } else {
            skip_framed(&mut r)?;
        }
    }
    // Regenerated names come from the template, not the (empty) name frames.
    let need_names = sel.names && !regen;
    let mut name_payloads: Vec<Vec<u8>> = Vec::with_capacity(n_blocks.min(1 << 20));
    let mut qual_payloads: Vec<Vec<u8>> = Vec::with_capacity(n_blocks.min(1 << 20));
    for i in 0..n_blocks {
        if need_names {
            name_payloads.push(read_framed(&mut r, &format!("reorder name block {i}"))?);
        } else {
            skip_framed(&mut r)?;
        }
        if sel.quality {
            qual_payloads.push(read_framed(&mut r, &format!("reorder quality block {i}"))?);
        } else {
            skip_framed(&mut r)?;
        }
    }
    // The whole-output digest folds name+sequence+quality jointly, so a partial
    // selection cannot recompute it; consume the frame to keep the stream in
    // sync (and its own CRC checked).
    let _digest = read_framed(&mut r, "reorder output digest")?;

    // Entropy-decode the selected partitions in parallel.
    let reference_ref = reference.as_ref();
    let mut cl_reads: Vec<Vec<u8>> = if sel.sequence {
        let blocks: Vec<Vec<Vec<u8>>> = pool.install(|| {
            seq_payloads
                .par_iter()
                .map(|p| Ok(fqxv_reorder::decode_clustered_any(p, reference_ref)?))
                .collect::<Result<_>>()
        })?;
        let mut v = reserve_checked(n)?;
        for b in blocks {
            v.extend(b);
        }
        if v.len() != n {
            return Err(Error::Malformed("reordered stream length disagreement"));
        }
        v
    } else {
        Vec::new()
    };
    let mut names: Vec<Vec<u8>> = if need_names {
        let blocks: Vec<Vec<Vec<u8>>> = pool.install(|| {
            name_payloads
                .par_iter()
                .map(|p| Ok(fqxv_tokenizer::decode(p)?))
                .collect::<Result<_>>()
        })?;
        let mut v = reserve_checked(n)?;
        for b in blocks {
            v.extend(b);
        }
        if v.len() != n {
            return Err(Error::Malformed("reordered stream length disagreement"));
        }
        v
    } else {
        Vec::new()
    };
    let (lens, quals): (Vec<u32>, Vec<u8>) = if sel.quality {
        let blocks: Vec<(Vec<u32>, Vec<u8>)> = pool.install(|| {
            qual_payloads
                .par_iter()
                .map(|p| Ok(fqxv_fqzcomp::decode(p)?))
                .collect::<Result<_>>()
        })?;
        let mut lens: Vec<u32> = reserve_checked(n)?;
        let mut quals: Vec<u8> = Vec::new();
        for (ls, qs) in blocks {
            lens.extend(ls);
            quals.extend(qs);
        }
        if lens.len() != n {
            return Err(Error::Malformed("reordered stream length disagreement"));
        }
        (lens, quals)
    } else {
        (Vec::new(), Vec::new())
    };

    let mut qoff = 0usize;
    if keep_order {
        // Sequences un-permute into original order; names/quality are already
        // stored in original order. (A keep-order archive never regenerates
        // names, so `regen` is always false here; guarded anyway.)
        let mut seq_orig = if sel.sequence {
            unpermute_reads(cl_reads, &flip, &perm_c, n)?
        } else {
            Vec::new()
        };
        for i in 0..n {
            let seq = if sel.sequence {
                std::mem::take(&mut seq_orig[i])
            } else {
                Vec::new()
            };
            let qual = if sel.quality {
                let l = lens[i] as usize;
                let q = quals
                    .get(qoff..qoff + l)
                    .ok_or(Error::Malformed("quality underrun"))?
                    .to_vec();
                qoff += l;
                q
            } else {
                Vec::new()
            };
            if sel.sequence && sel.quality && seq.len() != qual.len() {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            let name = if need_names {
                std::mem::take(&mut names[i])
            } else {
                Vec::new()
            };
            emit(Record { name, seq, qual })?;
        }
    } else {
        // Reads emerge in clustered order; un-flip each selected stream back to
        // the record's original content.
        for j in 0..n {
            let seq = if sel.sequence {
                let s = std::mem::take(&mut cl_reads[j]);
                if flip_bit(&flip, j) {
                    fqxv_reorder::revcomp(&s)
                } else {
                    s
                }
            } else {
                Vec::new()
            };
            let qual = if sel.quality {
                let l = lens[j] as usize;
                let mut q = quals
                    .get(qoff..qoff + l)
                    .ok_or(Error::Malformed("quality underrun"))?
                    .to_vec();
                qoff += l;
                if flip_bit(&flip, j) {
                    q.reverse();
                }
                q
            } else {
                Vec::new()
            };
            if sel.sequence && sel.quality && seq.len() != qual.len() {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            let name = if sel.names {
                if let Some(t) = &template {
                    t.regenerate(j)
                } else {
                    std::mem::take(&mut names[j])
                }
            } else {
                Vec::new()
            };
            emit(Record { name, seq, qual })?;
        }
    }
    Ok(Stats {
        reads: n as u64,
        blocks: n_blocks as u64,
        out_bytes: 0,
        group_size: group_size.max(1),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn compress_bytes(input: &[u8], params: Params) -> Vec<u8> {
        let mut out = Vec::new();
        compress(input, &mut out, params).expect("compress");
        out
    }

    /// Reference records from the canonical full decode.
    fn full_records(archive: &[u8]) -> Vec<Record> {
        let mut recs = Vec::new();
        decompress_records(archive, 1, |r| recs.push(r)).expect("full decode");
        recs
    }

    /// Records from the stream-selective decode.
    fn select_records(archive: &[u8], sel: StreamSelection) -> Vec<Record> {
        let mut recs = Vec::new();
        decompress_records_select(archive, 1, sel, |r| recs.push(r)).expect("select decode");
        recs
    }

    /// The selections worth exercising: every single stream, the useful pairs,
    /// nothing, and everything (the full-path delegation).
    const SELECTIONS: [StreamSelection; 8] = [
        StreamSelection::SEQUENCE_ONLY,
        StreamSelection::NAMES_ONLY,
        StreamSelection::QUALITY_ONLY,
        StreamSelection::NAMES_AND_SEQUENCE,
        StreamSelection {
            names: false,
            sequence: true,
            quality: true,
        },
        StreamSelection {
            names: true,
            sequence: false,
            quality: true,
        },
        StreamSelection::NONE,
        StreamSelection::ALL,
    ];

    /// Assert every selection of `archive` agrees with the full decode: selected
    /// fields byte-identical, deselected fields empty, record counts equal.
    fn assert_selections_agree(archive: &[u8], what: &str) {
        let expected = full_records(archive);
        for sel in SELECTIONS {
            let got = select_records(archive, sel);
            assert_eq!(got.len(), expected.len(), "{what} {sel:?}: record count");
            for (i, (g, e)) in got.iter().zip(&expected).enumerate() {
                if sel.names {
                    assert_eq!(g.name, e.name, "{what} {sel:?} record {i}: name");
                } else {
                    assert!(
                        g.name.is_empty(),
                        "{what} {sel:?} record {i}: name not empty"
                    );
                }
                if sel.sequence {
                    assert_eq!(g.seq, e.seq, "{what} {sel:?} record {i}: sequence");
                } else {
                    assert!(g.seq.is_empty(), "{what} {sel:?} record {i}: seq not empty");
                }
                if sel.quality {
                    assert_eq!(g.qual, e.qual, "{what} {sel:?} record {i}: quality");
                } else {
                    assert!(
                        g.qual.is_empty(),
                        "{what} {sel:?} record {i}: qual not empty"
                    );
                }
            }
        }
    }

    const SAMPLE: &[u8] = b"\
@r1 desc one\nACGTACGT\n+\nIIIIFFF#\n\
@r2 desc two\nNNGGCCTA\n+\n###IIIFF\n\
@r3\nTTAA\n+\nIIII\n\
@r4\nACACACAC\n+\nHHHHGGGG\n";

    #[test]
    fn plain_layout_selections_agree() {
        assert_selections_agree(&compress_bytes(SAMPLE, Params::default()), "plain");
    }

    #[test]
    fn plain_layout_multiple_blocks_selections_agree() {
        // Several small blocks so the batched block loop is exercised.
        let mut input = Vec::new();
        for i in 0..300u32 {
            input.extend_from_slice(
                format!("@read.{i} extra\nACGTACGTNACG\n+\nIIIIFFFF##II\n").as_bytes(),
            );
        }
        let archive = compress_bytes(
            &input,
            Params {
                block_reads: 32,
                ..Params::default()
            },
        );
        assert_selections_agree(&archive, "plain multi-block");
    }

    #[test]
    fn grouped_layout_selections_agree() {
        let r1: &[u8] = b"@spot1/1\nACGTACGT\n+\nIIIIFFFF\n@spot2/1\nTTGGCCAA\n+\nHHHHGGGG\n";
        let r2: &[u8] = b"@spot1/2\nCCCCGGGG\n+\nJJJJKKKK\n@spot2/2\nAATTCCGG\n+\nFFFFIIII\n";
        let readers: Vec<Box<dyn io::Read + Send>> = vec![Box::new(r1), Box::new(r2)];
        let mut archive = Vec::new();
        compress_multi(readers, &mut archive, Params::default()).unwrap();
        assert_selections_agree(&archive, "grouped");
    }

    #[test]
    fn reorder_keep_order_selections_agree() {
        // Force the clustered layout (the never-worse gate would pick plain on a
        // fixture this small) with the permutation kept.
        let mut archive = Vec::new();
        compress_clustered_forced(
            SAMPLE,
            &mut archive,
            Params {
                reorder: true,
                keep_order: true,
                ..Params::default()
            },
            1,
        )
        .unwrap();
        assert_selections_agree(&archive, "reorder keep-order");
    }

    #[test]
    fn reorder_discard_order_selections_agree() {
        // Clustered output order: the full decode emits clustered order too, so
        // positional comparison against it is exact.
        let mut archive = Vec::new();
        compress_clustered_forced(
            SAMPLE,
            &mut archive,
            Params {
                reorder: true,
                keep_order: false,
                ..Params::default()
            },
            1,
        )
        .unwrap();
        assert_selections_agree(&archive, "reorder discard-order");
    }

    #[test]
    fn reorder_regenerated_names_selections_agree() {
        // Discard-order with template-regenerated names: the name frames are
        // empty and names come from the stored template at each output position.
        let mut input = Vec::new();
        for i in 1..=20u32 {
            input
                .extend_from_slice(format!("@sra.{i} {i}\nACGTACGTAA\n+\nIIIIFFFFHH\n").as_bytes());
        }
        let mut archive = Vec::new();
        compress_clustered_forced(
            &input[..],
            &mut archive,
            Params {
                reorder: true,
                keep_order: false,
                regenerate_names: true,
                ..Params::default()
            },
            1,
        )
        .unwrap();
        assert_selections_agree(&archive, "reorder regen-names");
    }

    #[test]
    fn reorder_grouped_selections_agree() {
        let r1: &[u8] = b"@spot1/1\nACGTACGT\n+\nIIIIFFFF\n@spot2/1\nTTGGCCAA\n+\nHHHHGGGG\n";
        let r2: &[u8] = b"@spot1/2\nCCCCGGGG\n+\nJJJJKKKK\n@spot2/2\nAATTCCGG\n+\nFFFFIIII\n";
        let mut interleaved = Vec::new();
        // Interleave the two mates spot by spot, as compress_multi would.
        let (m1, m2) = (full_records_of(r1), full_records_of(r2));
        for (a, b) in m1.iter().zip(&m2) {
            write_record(&mut interleaved, &a.name, &a.seq, &a.qual);
            write_record(&mut interleaved, &b.name, &b.seq, &b.qual);
        }
        let mut archive = Vec::new();
        compress_clustered_forced(
            &interleaved[..],
            &mut archive,
            Params {
                reorder: true,
                ..Params::default()
            },
            2,
        )
        .unwrap();
        assert_selections_agree(&archive, "reorder grouped");
    }

    /// Parse raw FASTQ text into records (test-local helper).
    fn full_records_of(fastq: &[u8]) -> Vec<Record> {
        let mut recs = Vec::new();
        let mut lines = fastq.split(|&b| b == b'\n');
        while let Some(name) = lines.next() {
            if name.is_empty() {
                break;
            }
            let seq = lines.next().unwrap();
            let _plus = lines.next().unwrap();
            let qual = lines.next().unwrap();
            recs.push(Record {
                name: name[1..].to_vec(),
                seq: seq.to_vec(),
                qual: qual.to_vec(),
            });
        }
        recs
    }

    #[test]
    fn lossy_binning_selections_agree() {
        // A lossy-quality archive: the select decode must reproduce the stored
        // (post-binning) quality bytes, exactly as the full decode does.
        let archive = compress_bytes(
            SAMPLE,
            Params {
                quality_binning: QualityBinning::Bin8,
                ..Params::default()
            },
        );
        assert_selections_agree(&archive, "lossy bin8");
    }

    #[test]
    fn fasta_output_matches_full_decode() {
        let archive = compress_bytes(SAMPLE, Params::default());
        let mut fasta = Vec::new();
        let stats = decompress_fasta(&archive[..], &mut fasta, 1).unwrap();
        let mut expected = Vec::new();
        for rec in full_records(&archive) {
            expected.push(b'>');
            expected.extend_from_slice(&rec.name);
            expected.push(b'\n');
            expected.extend_from_slice(&rec.seq);
            expected.push(b'\n');
        }
        assert_eq!(fasta, expected, "FASTA must be >name + sequence per record");
        assert_eq!(stats.reads, 4);
        assert_eq!(
            stats.out_bytes,
            fasta.len() as u64,
            "out_bytes = FASTA bytes"
        );
    }

    #[test]
    fn record_reader_with_selection_agrees() {
        let archive = compress_bytes(SAMPLE, Params::default());
        for sel in SELECTIONS {
            let via_iter: Vec<Record> =
                RecordReader::with_selection(io::Cursor::new(archive.clone()), 1, sel)
                    .map(|r| r.expect("record"))
                    .collect();
            let via_push = select_records(&archive, sel);
            assert_eq!(via_iter, via_push, "{sel:?}: iterator vs callback");
        }
    }

    #[test]
    fn none_selection_counts_reads_without_decoding() {
        let archive = compress_bytes(SAMPLE, Params::default());
        let stats = decompress_records_select(&archive[..], 1, StreamSelection::NONE, |rec| {
            assert!(rec.name.is_empty() && rec.seq.is_empty() && rec.qual.is_empty());
        })
        .unwrap();
        assert_eq!(stats.reads, 4);
        assert_eq!(stats.out_bytes, 0);
    }

    /// Byte offset of the first *quality* frame's payload in a whole-file
    /// reorder archive, by walking its self-describing layout. Test-only: used
    /// to corrupt a deselected stream and prove it is genuinely skipped.
    fn reorder_first_qual_payload_offset(archive: &[u8]) -> usize {
        let mut c = io::Cursor::new(archive);
        let header = read_header(&mut c).expect("header");
        assert_ne!(
            header.flags & FLAG_GLOBAL_REORDER,
            0,
            "not a reorder archive"
        );
        let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
        let mut pos = header.header_len as usize + 8; // past [8] n
        let frame_len =
            |pos: usize| u32::from_le_bytes(archive[pos..pos + 4].try_into().unwrap()) as usize;
        let skip = |pos: &mut usize| *pos += 8 + frame_len(*pos);
        skip(&mut pos); // flip
        skip(&mut pos); // perm
        skip(&mut pos); // template
        if has_reference {
            skip(&mut pos); // shared reference
        }
        let n_blocks = u32::from_le_bytes(archive[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        for _ in 0..n_blocks {
            skip(&mut pos); // sequence block
        }
        skip(&mut pos); // first name block
        pos + 8 // past the qual frame's [len][crc] head
    }

    #[test]
    fn skipped_streams_are_genuinely_skipped_not_decoded() {
        // Corrupt a byte inside the first QUALITY frame of a reorder archive.
        // Its frame CRC now fails, so any decode that reads that stream errors —
        // while a selection that skips quality must still succeed byte-exact,
        // proving the frame was never read into the decoder at all.
        let mut archive = Vec::new();
        compress_clustered_forced(
            SAMPLE,
            &mut archive,
            Params {
                reorder: true,
                keep_order: true,
                ..Params::default()
            },
            1,
        )
        .unwrap();
        let expected = full_records(&archive);
        let off = reorder_first_qual_payload_offset(&archive);
        archive[off] ^= 0xFF;

        // Full decode reads the quality frame: must fail its CRC.
        let mut sink = Vec::new();
        assert!(
            decompress(&archive[..], &mut sink, 1).is_err(),
            "full decode must fail on the corrupted quality frame"
        );
        // Quality-selecting decode also fails.
        assert!(
            decompress_records_select(&archive[..], 1, StreamSelection::QUALITY_ONLY, |_| {})
                .is_err(),
            "quality-selecting decode must fail on the corrupted quality frame"
        );
        // A quality-skipping selection never touches the frame: succeeds, exact.
        let got = select_records(&archive, StreamSelection::NAMES_AND_SEQUENCE);
        assert_eq!(got.len(), expected.len());
        for (g, e) in got.iter().zip(&expected) {
            assert_eq!(g.name, e.name);
            assert_eq!(g.seq, e.seq);
            assert!(g.qual.is_empty());
        }
    }

    proptest! {
        // Mirrors records.rs's framing property at the same cost budget: each
        // case runs one compress and two decodes over a random corpus and a
        // random selection.
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(48))]

        #[test]
        fn prop_selections_agree_with_full_decode(
            reads in proptest::collection::vec(
                (
                    "[A-Za-z0-9:_.-]{1,20}",
                    proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')], 1..40),
                ),
                1..30,
            ),
            names in any::<bool>(),
            sequence in any::<bool>(),
            quality in any::<bool>(),
        ) {
            let mut input = Vec::new();
            for (i, (name, seq)) in reads.iter().enumerate() {
                input.extend_from_slice(b"@");
                input.extend_from_slice(name.as_bytes());
                input.extend_from_slice(format!("_{i}\n").as_bytes());
                input.extend_from_slice(seq);
                input.push(b'\n');
                input.extend_from_slice(b"+\n");
                input.extend(seq.iter().map(|_| b'I'));
                input.push(b'\n');
            }
            // Low order for the same reason as records.rs's proptest: this is a
            // framing/selection property, not a sequence-coder one.
            let archive = compress_bytes(&input, Params { seq_order: 2, ..Params::default() });
            let sel = StreamSelection { names, sequence, quality };
            let expected = full_records(&archive);
            let got = select_records(&archive, sel);
            prop_assert_eq!(got.len(), expected.len());
            let empty: Vec<u8> = Vec::new();
            for (g, e) in got.iter().zip(&expected) {
                prop_assert_eq!(&g.name, if sel.names { &e.name } else { &empty });
                prop_assert_eq!(&g.seq, if sel.sequence { &e.seq } else { &empty });
                prop_assert_eq!(&g.qual, if sel.quality { &e.qual } else { &empty });
            }
        }
    }
}
