//! Reference-free read reordering by canonical-minimizer clustering.
//!
//! Each read is reduced to its minimum canonical k-mer (the smaller of a k-mer
//! and its reverse complement, over the whole read). Sorting reads by that key
//! (then by oriented sequence) brings exact and reverse-complement duplicates —
//! and near-duplicates sharing a minimizer — next to each other, so the order-k
//! sequence model in `fqxv-seq` sees long runs of near-identical reads. This is
//! the cross-read redundancy lever (PgRC/SPRING-class), the piece that plain
//! per-read context modeling can't reach.
//!
//! [`plan`] returns the emission [`Plan::order`] and per-read [`Plan::flip`]
//! (whether a read is stored reverse-complemented so its minimizer is in
//! canonical orientation). The caller reorders names/sequence/quality by `order`,
//! reverse-complements + reverses the quality of flipped reads, and stores the
//! permutation to restore the original order.
//!
//! ```
//! use fqxv_reorder::{plan, revcomp};
//! // read 1 is the reverse complement of read 0 — they should cluster.
//! let a = b"ACGTTTGACCGATT";
//! let ra = revcomp(a);
//! let mut seq = a.to_vec();
//! seq.extend_from_slice(&ra);
//! let lens = [a.len() as u32, ra.len() as u32];
//! let p = plan(&lens, &seq, 7);
//! assert_eq!(p.order.len(), 2);
//! // exactly one of the mates is flipped so both store identically.
//! assert_ne!(p.flip[0], p.flip[1]);
//! ```

use std::borrow::Cow;
use std::collections::HashMap;

use rayon::prelude::*;
use thiserror::Error;

/// Errors returned by the reordering engine.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed clustered stream: {0}")]
    Malformed(&'static str),
    /// Underlying rANS coder failure.
    #[error(transparent)]
    Rans(#[from] fqxv_rans::Error),
    /// Underlying sequence codec failure.
    #[error(transparent)]
    Seq(#[from] fqxv_seq::Error),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// A reordering: the order to emit reads, and which reads are stored
/// reverse-complemented.
#[derive(Debug, Clone, Default)]
pub struct Plan {
    /// Original read indices in the order they should be emitted.
    pub order: Vec<u32>,
    /// Per original-read flag: store the read reverse-complemented.
    pub flip: Vec<bool>,
    /// Per original-read: start position of the clustering minimizer in the
    /// oriented read. Adjacent clustered reads share the minimizer, so the
    /// difference of anchors is their alignment shift (for overlap coding).
    pub anchor: Vec<u32>,
}

/// Default minimizer k-mer length.
pub const DEFAULT_K: usize = 15;

/// byte -> 2-bit code, 255 for non-ACGT.
#[inline]
fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 255,
    }
}

/// Reverse complement (complementing A/C/G/T in either case; other bytes,
/// including `N`, are passed through), reversed.
#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// Minimum canonical k-mer of `read` and whether that minimizer sits on the
/// reverse strand (i.e. the read should be flipped to canonicalize it).
/// Returns `(min_canonical_kmer, flip, anchor)` where `anchor` is the start
/// position of that minimizer k-mer in the *oriented* read (the read as stored:
/// reverse-complemented iff `flip`). The anchor lets clustered reads be aligned
/// by their shared minimizer for shifted-overlap coding.
fn min_canonical(read: &[u8], k: usize) -> (u64, bool, u32) {
    if read.len() < k || k == 0 || k > 32 {
        return (u64::MAX, false, 0);
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k as u64 - 1);
    let (mut fwd, mut rc, mut valid) = (0u64, 0u64, 0usize);
    let (mut best, mut best_flip, mut best_end) = (u64::MAX, false, 0usize);
    for (idx, &b) in read.iter().enumerate() {
        let c = code(b);
        if c == 255 {
            fwd = 0;
            rc = 0;
            valid = 0;
            continue;
        }
        let c = u64::from(c);
        fwd = ((fwd << 2) | c) & mask;
        rc = ((rc >> 2) | ((3 - c) << shift)) & mask;
        valid += 1;
        if valid >= k {
            let (canon, is_rc) = if fwd <= rc { (fwd, false) } else { (rc, true) };
            if canon < best {
                best = canon;
                best_flip = is_rc;
                best_end = idx; // last base of the minimizing k-mer
            }
        }
    }
    // No valid (N-free) k-mer found — no minimizer, no anchor.
    if best == u64::MAX {
        return (u64::MAX, false, 0);
    }
    // Anchor = start of the minimizer k-mer in the oriented read.
    let len = read.len();
    let anchor = if best_flip {
        (len - 1 - best_end) as u32
    } else {
        (best_end + 1 - k) as u32
    };
    (best, best_flip, anchor)
}

/// Build a clustering [`Plan`] for the reads in `seq` (lengths in `lens`).
///
/// `k` is the minimizer length (see [`DEFAULT_K`]); it is clamped to `1..=32`.
#[must_use]
pub fn plan(lens: &[u32], seq: &[u8], k: usize) -> Plan {
    let k = k.clamp(1, 32);
    let n = lens.len();

    // Byte offset of each read (so the key build can run in parallel).
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // (canonical minimizer, oriented sequence, original index, flip, anchor).
    // Building each key is independent, so it runs across cores. The oriented
    // sequence is only a sort tiebreak, so borrow the input for the common
    // non-flipped case and allocate (via `revcomp`) only when a read flips.
    type Key<'a> = (u64, Cow<'a, [u8]>, u32, bool, u32);
    let mut keys: Vec<Key<'_>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let read = &seq[offs[i]..offs[i + 1]];
            let (canon, flip, anchor) = min_canonical(read, k);
            let oriented: Cow<'_, [u8]> = if flip {
                Cow::Owned(revcomp(read))
            } else {
                Cow::Borrowed(read)
            };
            (canon, oriented, i as u32, flip, anchor)
        })
        .collect();

    // Parallel sort: cluster by minimizer, then within a cluster order by
    // anchor DESCENDING. Higher anchor = the shared minimizer sits later in the
    // read = the read starts earlier on the shared coordinate, so reads emerge
    // left-to-right and the contig assembler grows a reference rightward without
    // ever extending left. The oriented-sequence + original-index tie-breaks
    // make the comparator a TOTAL order, so the output is byte-identical
    // regardless of thread count (the determinism invariant).
    keys.par_sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| b.4.cmp(&a.4))
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
    });

    let order = keys.iter().map(|key| key.2).collect();
    let mut flip = vec![false; n];
    let mut anchor = vec![0u32; n];
    for key in &keys {
        flip[key.2 as usize] = key.3;
        anchor[key.2 as usize] = key.4;
    }
    Plan {
        order,
        flip,
        anchor,
    }
}

// --- clustered contig-assembly codec ----------------------------------------

/// Exact duplicate of the previous emitted read.
const OP_MATCH: u8 = 0;
/// A member of the current contig: placed at an offset on the growing consensus
/// reference, storing only the overlap's mismatches and the novel tail bases.
const OP_CONTIG: u8 = 1;
/// Seeds a new contig — coded as a literal via the `fqxv_seq` context model.
const OP_LITERAL: u8 = 2;
/// Minimum overlap with the contig to bother placing a read on it.
const MIN_CONTIG_OVERLAP: usize = 16;
/// Half-width of the offset search around the anchor-implied placement. The
/// shared minimizer fixes the offset exactly for substitution errors (an error
/// inside the minimizer k-mer would move the read to a different cluster), so
/// this window only has to absorb the small shifts that indels introduce. The
/// offset is stored explicitly, so widening the search is purely an
/// encoder-side choice the decoder never sees.
const OFF_SEARCH: i64 = 8;

fn zigzag(d: i64) -> u64 {
    ((d << 1) ^ (d >> 63)) as u64
}
fn unzigzag(z: u64) -> i64 {
    ((z >> 1) as i64) ^ -((z & 1) as i64)
}

/// A contig column: per-base A/C/G/T vote counts plus the current consensus
/// byte (the plurality base, or a first-seen non-ACGT byte until an ACGT wins).
#[derive(Clone)]
struct Column {
    votes: [u32; 4],
    base: u8,
}

/// Fold base `b` into a contig column, updating the consensus to the plurality
/// (ties resolve to the lowest A<C<G<T so encode and decode always agree).
/// Non-ACGT bytes don't vote, so a column keeps its first-seen byte until an
/// ACGT base wins — the same rule on both sides keeps the reference in sync.
#[inline]
fn cast_vote(col: &mut Column, b: u8) {
    let c = code(b);
    if c < 4 {
        col.votes[c as usize] += 1;
        // Highest count wins; on a tie the lowest base index wins.
        let best = (0..4)
            .max_by_key(|&i| (col.votes[i], std::cmp::Reverse(i)))
            .unwrap();
        col.base = b"ACGT"[best];
    }
}

/// Seed a fresh contig column from the first read to cover a position.
#[inline]
fn seed_column(b: u8) -> Column {
    let mut col = Column {
        votes: [0; 4],
        base: b, // first-seen, until an ACGT vote takes over
    };
    cast_vote(&mut col, b);
    col
}

/// Decide whether read `cur` (with minimizer `anchor`) can be placed on the
/// current `contig` (whose seed read sat at `ref_anchor`). Returns
/// `Some((offset, overlap, mismatch_positions))` when the placement is cheaper
/// than a literal, else `None`. Pure — no mutation — so both [`encode_clustered`]
/// and [`op_stats`] share one source of truth for the classification.
fn place_on_contig(
    contig: &[Column],
    cur: &[u8],
    anchor: u32,
    ref_anchor: u32,
) -> Option<(usize, usize, Vec<usize>)> {
    if contig.is_empty() || cur.is_empty() {
        return None;
    }
    // The shared-minimizer anchor gives the structurally-correct offset, which is
    // exact for substitution errors (an error inside the minimizer k-mer would
    // move the read to another cluster). Try that offset first — the common path.
    // Only if it fails acceptance do we search a small window around it to rescue
    // reads an indel has shifted off the anchor. The chosen offset is stored
    // explicitly, so the search is invisible to the decoder.
    let center = ref_anchor as i64 - anchor as i64;
    let try_off = |off: usize| -> Option<(usize, usize, Vec<usize>)> {
        let overlap = cur.len().min(contig.len() - off);
        if overlap == 0 || overlap < MIN_CONTIG_OVERLAP.min(cur.len()) {
            return None;
        }
        // Mismatches vs the CONSENSUS of all reads placed so far.
        let mism: Vec<usize> = (0..overlap)
            .filter(|&j| cur[j] != contig[off + j].base)
            .collect();
        let novel_n = cur.len() - overlap;
        // Cheap enough to be a real overlap, and smaller than a literal.
        (mism.len() <= overlap / 4 && novel_n + mism.len() * 2 < cur.len())
            .then_some((off, overlap, mism))
    };
    let anchor_ok = (center >= 0 && center as usize <= contig.len())
        .then(|| try_off(center as usize))
        .flatten();
    anchor_ok.or_else(|| {
        // Anchor offset was rejected: scan the window for the placement with the
        // fewest mismatches (ties nearest the anchor).
        let lo = (center - OFF_SEARCH).max(0);
        let hi = (center + OFF_SEARCH).min(contig.len() as i64);
        let mut best: Option<(usize, usize, Vec<usize>)> = None;
        let mut best_key = (usize::MAX, i64::MAX);
        for off in lo..=hi {
            if off == center {
                continue; // already tried
            }
            if let Some((o, ov, mism)) = try_off(off as usize) {
                let key = (mism.len(), (off - center).abs());
                if key < best_key {
                    best_key = key;
                    best = Some((o, ov, mism));
                }
            }
        }
        best
    })
}

/// Assemble reads that are already in clustered, left-to-right order (via
/// [`plan`]) into contigs and code them against a growing consensus reference.
///
/// A read is one of: `MATCH` (identical to the previous read); `CONTIG` (it
/// overlaps the current contig at the shift implied by the shared minimizer
/// [`Plan::anchor`] — store the overlap's mismatches and the novel tail, which
/// extends the reference); or `LITERAL` (seeds a new contig, coded with the
/// [`fqxv_seq`] model). This captures the shifted overlaps of deep coverage —
/// a read matches the consensus of all reads before it on the contig, not just
/// its immediate predecessor. Byte-exact.
pub fn encode_clustered(reads: &[&[u8]], anchors: &[u32], seq_order: usize) -> Result<Vec<u8>> {
    let mut ops = Vec::with_capacity(reads.len());
    let (mut offdelta, mut slen) = (Vec::new(), Vec::new());
    let (mut nmis, mut pos, mut subs) = (Vec::new(), Vec::new(), Vec::new());
    let (mut novel, mut lit_seq, mut lit_lens): (Vec<u8>, Vec<u8>, Vec<u32>) =
        (Vec::new(), Vec::new(), Vec::new());

    // The current contig: a growing plurality-consensus reference (one voting
    // `Column` per position). `ref_anchor` is the shared minimizer's position in
    // it (the seed read's anchor); `prev_off` delta-codes offsets.
    let mut contig: Vec<Column> = Vec::new();
    let mut ref_anchor: u32 = 0;
    let mut prev_off: usize = 0;

    for (i, &cur) in reads.iter().enumerate() {
        if i > 0 && cur == reads[i - 1] {
            ops.push(OP_MATCH);
            continue;
        }
        // Place `cur` on the contig (shared-minimizer anchor, small indel-rescue
        // window). See [`place_on_contig`].
        let placed = place_on_contig(&contig, cur, anchors[i], ref_anchor);
        match placed {
            Some((off, overlap, mism)) => {
                ops.push(OP_CONTIG);
                write_varint(&mut offdelta, zigzag(off as i64 - prev_off as i64));
                write_varint(&mut slen, cur.len() as u64);
                write_varint(&mut nmis, mism.len() as u64);
                let mut last = 0usize;
                for &m in &mism {
                    write_varint(&mut pos, (m - last) as u64);
                    last = m;
                    subs.push(cur[m]);
                }
                novel.extend_from_slice(&cur[overlap..]);
                // Fold this read into the consensus for the reads that follow.
                for (j, &b) in cur.iter().enumerate().take(overlap) {
                    cast_vote(&mut contig[off + j], b);
                }
                for &b in &cur[overlap..] {
                    contig.push(seed_column(b));
                }
                prev_off = off;
            }
            None => {
                ops.push(OP_LITERAL);
                lit_seq.extend_from_slice(cur);
                lit_lens.push(cur.len() as u32);
                contig = cur.iter().map(|&b| seed_column(b)).collect();
                ref_anchor = anchors[i];
                prev_off = 0;
            }
        }
    }

    let ops_c = fqxv_rans::encode(&ops, fqxv_rans::Order::One)?;
    let offdelta_c = fqxv_rans::encode(&offdelta, fqxv_rans::Order::Zero)?;
    let slen_c = fqxv_rans::encode(&slen, fqxv_rans::Order::Zero)?;
    let nmis_c = fqxv_rans::encode(&nmis, fqxv_rans::Order::Zero)?;
    let pos_c = fqxv_rans::encode(&pos, fqxv_rans::Order::Zero)?;
    let subs_c = fqxv_rans::encode(&subs, fqxv_rans::Order::One)?;
    let novel_c = fqxv_seq::encode(&[novel.len() as u32], &novel, seq_order)?;
    let lit_c = fqxv_seq::encode(&lit_lens, &lit_seq, seq_order)?;

    let mut out = Vec::new();
    out.push(2u8); // version 2: contig-assembly layout
    write_varint(&mut out, reads.len() as u64);
    for s in [
        &ops_c,
        &offdelta_c,
        &slen_c,
        &nmis_c,
        &pos_c,
        &subs_c,
        &novel_c,
        &lit_c,
    ] {
        write_varint(&mut out, s.len() as u64);
        out.extend_from_slice(s);
    }
    Ok(out)
}

/// Op-mix tally from the clustered contig-assembly codec — a diagnostic that
/// replays [`encode_clustered`]'s classification and consensus updates exactly
/// (via the shared [`place_on_contig`]) but skips entropy coding, so the counts
/// reflect what the real encoder does. A high `literals` / `literal_bases` share
/// is the signal that clustering is leaving cross-read redundancy uncaptured.
///
/// Call it per block on the same clustered, oriented slices the container feeds
/// [`encode_clustered`]; the contig resets at each call, matching the per-block
/// encoding. Fields are additive across blocks (see [`OpStats::merge`]).
#[derive(Debug, Default, Clone)]
pub struct OpStats {
    /// Reads seen.
    pub reads: usize,
    /// Reads coded as `MATCH` (byte-identical to the previous read).
    pub matches: usize,
    /// Reads placed on a contig (`CONTIG`).
    pub contigs: usize,
    /// Reads that seeded a fresh contig (`LITERAL`) — context-coded from scratch.
    pub literals: usize,
    /// Total substitution mismatches across all `CONTIG` reads.
    pub contig_mismatches: u64,
    /// Total overlap bases coded differentially (as offset + mismatches).
    pub contig_overlap_bases: u64,
    /// Total novel-tail bases (the `CONTIG` overhang past the contig) — these go
    /// to the `fqxv_seq` context model, so they cost like literal bases.
    pub novel_tail_bases: u64,
    /// Total bases in `LITERAL` reads — context-coded from scratch.
    pub literal_bases: u64,
    /// Total bases in `MATCH` reads — coded for free (one op symbol).
    pub match_bases: u64,
    /// All bases seen (overlap + novel tail + literal + match).
    pub total_bases: u64,
}

impl OpStats {
    /// Add another block's tally into this one.
    pub fn merge(&mut self, o: &OpStats) {
        self.reads += o.reads;
        self.matches += o.matches;
        self.contigs += o.contigs;
        self.literals += o.literals;
        self.contig_mismatches += o.contig_mismatches;
        self.contig_overlap_bases += o.contig_overlap_bases;
        self.novel_tail_bases += o.novel_tail_bases;
        self.literal_bases += o.literal_bases;
        self.match_bases += o.match_bases;
        self.total_bases += o.total_bases;
    }
}

/// Classify a clustered, oriented block of reads exactly as [`encode_clustered`]
/// would and return the [`OpStats`] tally — no entropy coding, no output. `reads`
/// and `anchors` are the same slices the container passes to `encode_clustered`.
pub fn op_stats(reads: &[&[u8]], anchors: &[u32]) -> OpStats {
    let mut st = OpStats::default();
    let mut contig: Vec<Column> = Vec::new();
    let mut ref_anchor: u32 = 0;
    for (i, &cur) in reads.iter().enumerate() {
        st.reads += 1;
        st.total_bases += cur.len() as u64;
        if i > 0 && cur == reads[i - 1] {
            st.matches += 1;
            st.match_bases += cur.len() as u64;
            continue;
        }
        match place_on_contig(&contig, cur, anchors[i], ref_anchor) {
            Some((off, overlap, mism)) => {
                st.contigs += 1;
                st.contig_mismatches += mism.len() as u64;
                st.contig_overlap_bases += overlap as u64;
                st.novel_tail_bases += (cur.len() - overlap) as u64;
                for (j, &b) in cur.iter().enumerate().take(overlap) {
                    cast_vote(&mut contig[off + j], b);
                }
                for &b in &cur[overlap..] {
                    contig.push(seed_column(b));
                }
            }
            None => {
                st.literals += 1;
                st.literal_bases += cur.len() as u64;
                contig = cur.iter().map(|&b| seed_column(b)).collect();
                ref_anchor = anchors[i];
            }
        }
    }
    st
}

/// Decode a stream produced by [`encode_clustered`], returning the reads in the
/// same (clustered) order.
pub fn decode_clustered(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != 2 {
        return Err(Error::Malformed("unsupported version"));
    }
    let n = r.varint()? as usize;
    let ops = fqxv_rans::decode(r.take_stream()?)?;
    let offdelta = fqxv_rans::decode(r.take_stream()?)?;
    let slen = fqxv_rans::decode(r.take_stream()?)?;
    let nmis = fqxv_rans::decode(r.take_stream()?)?;
    let pos = fqxv_rans::decode(r.take_stream()?)?;
    let subs = fqxv_rans::decode(r.take_stream()?)?;
    let (_, novel) = fqxv_seq::decode(r.take_stream()?)?;
    let (lit_lens, lit_seq) = fqxv_seq::decode(r.take_stream()?)?;

    let mut c_offdelta = Cursor::new(&offdelta);
    let mut c_slen = Cursor::new(&slen);
    let mut c_nmis = Cursor::new(&nmis);
    let mut c_pos = Cursor::new(&pos);
    let (mut subs_pos, mut lit_pos, mut lit_idx, mut novel_pos) = (0usize, 0usize, 0usize, 0usize);
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n.min(1 << 22));

    // The current contig, voted identically to the encoder.
    let mut contig: Vec<Column> = Vec::new();
    let mut prev_off: usize = 0;

    for i in 0..n {
        let op = *ops.get(i).ok_or(Error::Malformed("op underrun"))?;
        match op {
            OP_MATCH => {
                let read = reads
                    .last()
                    .ok_or(Error::Malformed("MATCH with no previous"))?
                    .clone();
                reads.push(read);
            }
            OP_LITERAL => {
                let l = *lit_lens
                    .get(lit_idx)
                    .ok_or(Error::Malformed("lit len underrun"))? as usize;
                lit_idx += 1;
                let bytes = lit_seq
                    .get(lit_pos..lit_pos + l)
                    .ok_or(Error::Malformed("lit data underrun"))?
                    .to_vec();
                lit_pos += l;
                contig = bytes.iter().map(|&b| seed_column(b)).collect();
                prev_off = 0;
                reads.push(bytes);
            }
            OP_CONTIG => {
                let off = usize::try_from(prev_off as i64 + unzigzag(c_offdelta.varint()?))
                    .map_err(|_| Error::Malformed("bad contig offset"))?;
                if off > contig.len() {
                    return Err(Error::Malformed("contig offset past reference"));
                }
                let cur_len = c_slen.varint()? as usize;
                let overlap = cur_len.min(contig.len() - off);
                let mut read = vec![0u8; cur_len];
                for (j, slot) in read.iter_mut().enumerate().take(overlap) {
                    *slot = contig[off + j].base; // consensus of prior reads
                }
                let m = c_nmis.varint()? as usize;
                let mut p = 0usize;
                for _ in 0..m {
                    p += c_pos.varint()? as usize;
                    let b = *subs
                        .get(subs_pos)
                        .ok_or(Error::Malformed("subs underrun"))?;
                    subs_pos += 1;
                    *read
                        .get_mut(p)
                        .ok_or(Error::Malformed("mismatch position out of range"))? = b;
                }
                for slot in read.iter_mut().skip(overlap) {
                    *slot = *novel
                        .get(novel_pos)
                        .ok_or(Error::Malformed("novel underrun"))?;
                    novel_pos += 1;
                }
                // Fold this read into the consensus, exactly as the encoder did.
                for (j, &b) in read.iter().enumerate().take(overlap) {
                    cast_vote(&mut contig[off + j], b);
                }
                for &b in &read[overlap..] {
                    contig.push(seed_column(b));
                }
                prev_off = off;
                reads.push(read);
            }
            _ => return Err(Error::Malformed("unknown op")),
        }
    }
    Ok(reads)
}

// --- literal-rescue contig-assembly codec (prototype, version 3) -------------
//
// The version-2 codec above keeps a SINGLE active contig: a read that fails to
// place on it seeds a fresh contig and the old one is discarded. On deep data
// that strands ~15% of reads as LITERALs (context-coded from scratch) even
// though they overlap reads on an *earlier* contig — the redundancy SPRING's
// assembly captures. This codec keeps every contig alive and, before a read
// becomes a literal, looks it up against a k-mer index of all contigs so it can
// attach to whichever one it overlaps. The index is ENCODER-ONLY: each CONTIG
// read stores the contig it landed on (a small back-reference) plus its offset,
// so the decoder never searches — it just replays votes into the same contigs.

/// k-mer length for the rescue index (matches the clustering minimizer k).
const RESCUE_K: usize = DEFAULT_K;

/// Forward 2-bit k-mer packed from `seq[start..start+k]`, or `None` if the
/// window runs off the end or contains a non-ACGT byte. `k <= 32`.
#[inline]
fn kmer_at(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }
    let mut v = 0u64;
    for &b in &seq[start..start + k] {
        let c = code(b);
        if c >= 4 {
            return None;
        }
        v = (v << 2) | u64::from(c);
    }
    Some(v)
}

/// A chosen placement of a read on a contig.
struct Placement {
    ci: usize,
    off: usize,
    overlap: usize,
    mism: Vec<usize>,
}

/// Multi-contig assembler with an encoder-side k-mer index. Shared by the
/// rescue encoder and its op-mix diagnostic so both make identical decisions.
#[derive(Default)]
struct Assembler {
    contigs: Vec<Vec<Column>>,
    ref_anchors: Vec<u32>,
    /// k-mer -> (contig index, column position); most-recent occurrence wins.
    index: HashMap<u64, (u32, u32)>,
}

/// Acceptance test: can `cur` sit on `contig` at `off`? Returns
/// `(overlap, mismatch_positions)` when it is cheaper than a literal.
fn try_place(contig: &[Column], cur: &[u8], off: usize) -> Option<(usize, Vec<usize>)> {
    if off > contig.len() {
        return None;
    }
    let overlap = cur.len().min(contig.len() - off);
    if overlap == 0 || overlap < MIN_CONTIG_OVERLAP.min(cur.len()) {
        return None;
    }
    let mism: Vec<usize> = (0..overlap).filter(|&j| cur[j] != contig[off + j].base).collect();
    let novel_n = cur.len() - overlap;
    (mism.len() <= overlap / 4 && novel_n + mism.len() * 2 < cur.len()).then_some((overlap, mism))
}

impl Assembler {
    /// Index every k-mer starting in `[from, to)` of contig `ci`'s consensus.
    /// Only called on freshly-appended columns, so cost is linear in new bases;
    /// overlap columns whose consensus later shifts are left stale on purpose —
    /// the index only proposes candidates, [`try_place`] validates against the
    /// live consensus, so staleness costs recall, never correctness.
    fn index_range(&mut self, ci: usize, from: usize, to: usize) {
        let n = self.contigs[ci].len();
        let hi = to.min(n.saturating_sub(RESCUE_K - 1));
        for start in from..hi {
            let mut v = 0u64;
            let mut ok = true;
            for j in 0..RESCUE_K {
                let c = code(self.contigs[ci][start + j].base);
                if c >= 4 {
                    ok = false;
                    break;
                }
                v = (v << 2) | u64::from(c);
            }
            if ok {
                self.index.insert(v, (ci as u32, start as u32));
            }
        }
    }

    /// Best placement of `cur` (minimizer at `anchor`) across all contigs, or
    /// `None` if it should seed a new one. Candidates come from the most-recent
    /// contig at the anchor-implied offset (the v2 fast path) plus every contig a
    /// sampled read k-mer points at. Deterministic: candidates are deduped and
    /// scored by (mismatches, recency, offset), independent of hash iteration.
    fn place(&self, cur: &[u8], anchor: u32) -> Option<Placement> {
        if cur.is_empty() || self.contigs.is_empty() {
            return None;
        }
        let mut cands: Vec<(usize, usize)> = Vec::new();
        let last = self.contigs.len() - 1;
        let center = self.ref_anchors[last] as i64 - anchor as i64;
        if center >= 0 && center as usize <= self.contigs[last].len() {
            cands.push((last, center as usize));
        }
        // Non-overlapping k-mers cover every base, so a read with a few errors
        // still has clean k-mers to match on.
        let mut start = 0;
        while start + RESCUE_K <= cur.len() {
            if let Some(code) = kmer_at(cur, start, RESCUE_K) {
                if let Some(&(ci, cpos)) = self.index.get(&code) {
                    let off = cpos as i64 - start as i64;
                    if off >= 0 && off as usize <= self.contigs[ci as usize].len() {
                        cands.push((ci as usize, off as usize));
                    }
                }
            }
            start += RESCUE_K;
        }
        cands.sort_unstable();
        cands.dedup();

        let mut best: Option<Placement> = None;
        let mut best_key = (usize::MAX, usize::MAX, usize::MAX);
        for (ci, off) in cands {
            if let Some((overlap, mism)) = try_place(&self.contigs[ci], cur, off) {
                let key = (mism.len(), self.contigs.len() - 1 - ci, off);
                if key < best_key {
                    best_key = key;
                    best = Some(Placement { ci, off, overlap, mism });
                }
            }
        }
        best
    }

    /// Fold a placed read into contig `ci`'s consensus, extending it and
    /// indexing the newly-appended columns.
    fn commit(&mut self, ci: usize, cur: &[u8], off: usize, overlap: usize) {
        let old_len = self.contigs[ci].len();
        for (j, &b) in cur.iter().enumerate().take(overlap) {
            cast_vote(&mut self.contigs[ci][off + j], b);
        }
        for &b in &cur[overlap..] {
            self.contigs[ci].push(seed_column(b));
        }
        let new_len = self.contigs[ci].len();
        if new_len > old_len {
            let from = old_len.saturating_sub(RESCUE_K - 1);
            self.index_range(ci, from, new_len);
        }
    }

    /// Seed a fresh contig from a literal read and index all its k-mers.
    fn seed(&mut self, cur: &[u8], anchor: u32) {
        let ci = self.contigs.len();
        self.contigs.push(cur.iter().map(|&b| seed_column(b)).collect());
        self.ref_anchors.push(anchor);
        self.index_range(ci, 0, cur.len());
    }
}

/// Literal-rescue variant of [`encode_clustered`]: keeps every contig alive and
/// attaches would-be literals to any contig they overlap (see the module note on
/// the version-3 codec). Byte-exactly reversible by [`decode_clustered_rescue`].
pub fn encode_clustered_rescue(
    reads: &[&[u8]],
    anchors: &[u32],
    seq_order: usize,
) -> Result<Vec<u8>> {
    let mut ops = Vec::with_capacity(reads.len());
    let (mut cref, mut offdelta, mut slen) = (Vec::new(), Vec::new(), Vec::new());
    let (mut nmis, mut pos, mut subs) = (Vec::new(), Vec::new(), Vec::new());
    let (mut novel, mut lit_seq, mut lit_lens): (Vec<u8>, Vec<u8>, Vec<u32>) =
        (Vec::new(), Vec::new(), Vec::new());

    let mut asm = Assembler::default();
    // Per-contig previous offset, for delta-coding offsets within a contig.
    let mut last_off: Vec<usize> = Vec::new();

    for (i, &cur) in reads.iter().enumerate() {
        if i > 0 && cur == reads[i - 1] {
            ops.push(OP_MATCH);
            continue;
        }
        match asm.place(cur, anchors[i]) {
            Some(p) => {
                ops.push(OP_CONTIG);
                // Back-reference: contigs ago (0 = most recent). Small under
                // clustered order, so it entropy-codes cheaply.
                write_varint(&mut cref, (asm.contigs.len() - 1 - p.ci) as u64);
                write_varint(&mut offdelta, zigzag(p.off as i64 - last_off[p.ci] as i64));
                write_varint(&mut slen, cur.len() as u64);
                write_varint(&mut nmis, p.mism.len() as u64);
                let mut last = 0usize;
                for &m in &p.mism {
                    write_varint(&mut pos, (m - last) as u64);
                    last = m;
                    subs.push(cur[m]);
                }
                novel.extend_from_slice(&cur[p.overlap..]);
                last_off[p.ci] = p.off;
                asm.commit(p.ci, cur, p.off, p.overlap);
            }
            None => {
                ops.push(OP_LITERAL);
                lit_seq.extend_from_slice(cur);
                lit_lens.push(cur.len() as u32);
                asm.seed(cur, anchors[i]);
                last_off.push(0);
            }
        }
    }

    let ops_c = fqxv_rans::encode(&ops, fqxv_rans::Order::One)?;
    let cref_c = fqxv_rans::encode(&cref, fqxv_rans::Order::Zero)?;
    let offdelta_c = fqxv_rans::encode(&offdelta, fqxv_rans::Order::Zero)?;
    let slen_c = fqxv_rans::encode(&slen, fqxv_rans::Order::Zero)?;
    let nmis_c = fqxv_rans::encode(&nmis, fqxv_rans::Order::Zero)?;
    let pos_c = fqxv_rans::encode(&pos, fqxv_rans::Order::Zero)?;
    let subs_c = fqxv_rans::encode(&subs, fqxv_rans::Order::One)?;
    let novel_c = fqxv_seq::encode(&[novel.len() as u32], &novel, seq_order)?;
    let lit_c = fqxv_seq::encode(&lit_lens, &lit_seq, seq_order)?;

    let mut out = Vec::new();
    out.push(3u8); // version 3: literal-rescue contig-assembly layout
    write_varint(&mut out, reads.len() as u64);
    for s in [
        &ops_c, &cref_c, &offdelta_c, &slen_c, &nmis_c, &pos_c, &subs_c, &novel_c, &lit_c,
    ] {
        write_varint(&mut out, s.len() as u64);
        out.extend_from_slice(s);
    }
    Ok(out)
}

/// Decode a stream from [`encode_clustered_rescue`], returning the reads in
/// clustered order. Maintains the same set of contigs the encoder built (no
/// k-mer index needed — each read carries its contig back-reference and offset).
pub fn decode_clustered_rescue(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != 3 {
        return Err(Error::Malformed("unsupported version"));
    }
    let n = r.varint()? as usize;
    let ops = fqxv_rans::decode(r.take_stream()?)?;
    let cref = fqxv_rans::decode(r.take_stream()?)?;
    let offdelta = fqxv_rans::decode(r.take_stream()?)?;
    let slen = fqxv_rans::decode(r.take_stream()?)?;
    let nmis = fqxv_rans::decode(r.take_stream()?)?;
    let pos = fqxv_rans::decode(r.take_stream()?)?;
    let subs = fqxv_rans::decode(r.take_stream()?)?;
    let (_, novel) = fqxv_seq::decode(r.take_stream()?)?;
    let (lit_lens, lit_seq) = fqxv_seq::decode(r.take_stream()?)?;

    let mut c_cref = Cursor::new(&cref);
    let mut c_offdelta = Cursor::new(&offdelta);
    let mut c_slen = Cursor::new(&slen);
    let mut c_nmis = Cursor::new(&nmis);
    let mut c_pos = Cursor::new(&pos);
    let (mut subs_pos, mut lit_pos, mut lit_idx, mut novel_pos) = (0usize, 0usize, 0usize, 0usize);
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n.min(1 << 22));

    let mut contigs: Vec<Vec<Column>> = Vec::new();
    let mut last_off: Vec<usize> = Vec::new();

    for i in 0..n {
        let op = *ops.get(i).ok_or(Error::Malformed("op underrun"))?;
        match op {
            OP_MATCH => {
                let read = reads
                    .last()
                    .ok_or(Error::Malformed("MATCH with no previous"))?
                    .clone();
                reads.push(read);
            }
            OP_LITERAL => {
                let l = *lit_lens
                    .get(lit_idx)
                    .ok_or(Error::Malformed("lit len underrun"))? as usize;
                lit_idx += 1;
                let bytes = lit_seq
                    .get(lit_pos..lit_pos + l)
                    .ok_or(Error::Malformed("lit data underrun"))?
                    .to_vec();
                lit_pos += l;
                contigs.push(bytes.iter().map(|&b| seed_column(b)).collect());
                last_off.push(0);
                reads.push(bytes);
            }
            OP_CONTIG => {
                let back = c_cref.varint()? as usize;
                let ci = contigs
                    .len()
                    .checked_sub(1 + back)
                    .ok_or(Error::Malformed("contig back-reference out of range"))?;
                let off = usize::try_from(last_off[ci] as i64 + unzigzag(c_offdelta.varint()?))
                    .map_err(|_| Error::Malformed("bad contig offset"))?;
                if off > contigs[ci].len() {
                    return Err(Error::Malformed("contig offset past reference"));
                }
                let cur_len = c_slen.varint()? as usize;
                let overlap = cur_len.min(contigs[ci].len() - off);
                let mut read = vec![0u8; cur_len];
                for (j, slot) in read.iter_mut().enumerate().take(overlap) {
                    *slot = contigs[ci][off + j].base;
                }
                let m = c_nmis.varint()? as usize;
                let mut p = 0usize;
                for _ in 0..m {
                    p += c_pos.varint()? as usize;
                    let b = *subs.get(subs_pos).ok_or(Error::Malformed("subs underrun"))?;
                    subs_pos += 1;
                    *read
                        .get_mut(p)
                        .ok_or(Error::Malformed("mismatch position out of range"))? = b;
                }
                for slot in read.iter_mut().skip(overlap) {
                    *slot = *novel.get(novel_pos).ok_or(Error::Malformed("novel underrun"))?;
                    novel_pos += 1;
                }
                for (j, &b) in read.iter().enumerate().take(overlap) {
                    cast_vote(&mut contigs[ci][off + j], b);
                }
                for &b in &read[overlap..] {
                    contigs[ci].push(seed_column(b));
                }
                last_off[ci] = off;
                reads.push(read);
            }
            _ => return Err(Error::Malformed("unknown op")),
        }
    }
    Ok(reads)
}

/// Op-mix tally for the literal-rescue codec — the [`op_stats`] analogue for
/// [`encode_clustered_rescue`], driving the same [`Assembler`] so the counts
/// match the encoder. Lets the diagnostic measure how many literals the rescue
/// pass recovers.
pub fn op_stats_rescue(reads: &[&[u8]], anchors: &[u32]) -> OpStats {
    let mut st = OpStats::default();
    let mut asm = Assembler::default();
    for (i, &cur) in reads.iter().enumerate() {
        st.reads += 1;
        st.total_bases += cur.len() as u64;
        if i > 0 && cur == reads[i - 1] {
            st.matches += 1;
            st.match_bases += cur.len() as u64;
            continue;
        }
        match asm.place(cur, anchors[i]) {
            Some(p) => {
                st.contigs += 1;
                st.contig_mismatches += p.mism.len() as u64;
                st.contig_overlap_bases += p.overlap as u64;
                st.novel_tail_bases += (cur.len() - p.overlap) as u64;
                asm.commit(p.ci, cur, p.off, p.overlap);
            }
            None => {
                st.literals += 1;
                st.literal_bases += cur.len() as u64;
                asm.seed(cur, anchors[i]);
            }
        }
    }
    st
}

fn write_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            break;
        }
        out.push(byte | 0x80);
    }
}

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Cursor { buf, pos: 0 }
    }
    fn u8(&mut self) -> Result<u8> {
        let b = *self
            .buf
            .get(self.pos)
            .ok_or(Error::Malformed("truncated"))?;
        self.pos += 1;
        Ok(b)
    }
    fn varint(&mut self) -> Result<u64> {
        let (mut v, mut shift) = (0u64, 0u32);
        loop {
            let b = self.u8()?;
            v |= u64::from(b & 0x7f) << shift;
            if b & 0x80 == 0 {
                return Ok(v);
            }
            shift += 7;
            if shift >= 64 {
                return Err(Error::Malformed("varint too long"));
            }
        }
    }
    fn take_stream(&mut self) -> Result<&'a [u8]> {
        let n = self.varint()? as usize;
        let end = self.pos + n;
        let s = self
            .buf
            .get(self.pos..end)
            .ok_or(Error::Malformed("stream truncated"))?;
        self.pos = end;
        Ok(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_basic() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
        assert_eq!(revcomp(b"AACG"), b"CGTT");
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
    }

    #[test]
    fn order_is_a_permutation() {
        let reads: Vec<&[u8]> = vec![b"ACGTACGT", b"TTTTGGGG", b"ACGTACGT", b"CCCCAAAA"];
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.concat();
        let p = plan(&lens, &seq, 5);
        let mut sorted = p.order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3]);
        assert_eq!(p.flip.len(), 4);
    }

    #[test]
    fn duplicates_cluster_adjacently() {
        // Two copies of read A and two of read B, shuffled; A's and B's should
        // each become adjacent after planning.
        let a: &[u8] = b"ACGTTTGACCGATTGCA";
        let b: &[u8] = b"GGGGCCCCAAAATTTTG";
        let reads = [a, b, a, b, a];
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.concat();
        let p = plan(&lens, &seq, 9);
        // Map emitted order back to which template (a=0/b=1) each read is.
        let tmpl: Vec<u8> = p
            .order
            .iter()
            .map(|&i| u8::from(reads[i as usize] == b))
            .collect();
        // Count adjacency changes; clustered => few transitions (<= 1).
        let transitions = tmpl.windows(2).filter(|w| w[0] != w[1]).count();
        assert!(transitions <= 1, "reads not clustered: {tmpl:?}");
    }

    #[test]
    fn revcomp_duplicate_flips_to_match() {
        let a = b"ACGTTTGACCGATTGCA";
        let ra = revcomp(a);
        let seq: Vec<u8> = a.iter().chain(ra.iter()).copied().collect();
        let lens = [a.len() as u32, ra.len() as u32];
        let p = plan(&lens, &seq, 9);
        // The two mates share a canonical minimizer; exactly one is flipped so
        // both are stored in the same orientation.
        assert_ne!(p.flip[0], p.flip[1]);
        // After orienting, the stored sequences are identical.
        let s0 = if p.flip[0] { revcomp(a) } else { a.to_vec() };
        let s1 = if p.flip[1] { revcomp(&ra) } else { ra.clone() };
        assert_eq!(s0, s1);
    }

    #[test]
    fn handles_short_and_n_reads() {
        let reads: Vec<&[u8]> = vec![b"AC", b"NNNNNNN", b"", b"ACGTACGTAC"];
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.concat();
        let p = plan(&lens, &seq, 5);
        assert_eq!(p.order.len(), 4);
    }

    #[test]
    fn clustered_roundtrip() {
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGTACGT", // literal (first)
            b"ACGTACGTACGT", // match
            b"ACGTAAGTACGT", // delta (1 mismatch)
            b"ACGTACGTNCGT", // delta including an N
            b"TTTTGGGGCCCC", // literal
            b"",             // empty read
            b"",             // match (empty == empty)
        ];
        let enc = encode_clustered(&reads, &vec![0u32; reads.len()], 4).expect("encode");
        let dec = decode_clustered(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn clustered_empty() {
        let enc = encode_clustered(&[], &[], 4).unwrap();
        assert!(decode_clustered(&enc).unwrap().is_empty());
    }

    #[test]
    fn clustered_deduplicates() {
        // Many copies of one read collapse to MATCH ops; the encoded size is
        // dominated by fixed per-stream overhead, not the read count.
        let read = b"ACGTTTGACCGATTGCAACGTTTGACCGATTGCA";
        let reads: Vec<&[u8]> = vec![&read[..]; 10_000];
        let raw = read.len() * reads.len();
        let enc = encode_clustered(&reads, &vec![0u32; reads.len()], 6).unwrap();
        assert!(
            enc.len() < raw / 50,
            "expected heavy dedup, got {} for {raw} raw",
            enc.len()
        );
        assert_eq!(decode_clustered(&enc).unwrap().len(), 10_000);
    }

    #[test]
    fn clustered_shift_overlap_roundtrips() {
        // Overlapping windows of a reference share sequence at a shift, which
        // should trigger SHIFT ops (and round-trip exactly).
        let reference = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGCATCGATCGTAGCTAGCAT";
        let win = 30usize;
        let (mut lens, mut seq) = (Vec::new(), Vec::new());
        for start in 0..=(reference.len() - win) {
            seq.extend_from_slice(&reference[start..start + win]);
            lens.push(win as u32);
        }
        let p = plan(&lens, &seq, DEFAULT_K);
        let mut offs = vec![0usize];
        for &l in &lens {
            offs.push(offs.last().unwrap() + l as usize);
        }
        let cl: Vec<Vec<u8>> = p
            .order
            .iter()
            .map(|&oi| {
                let oi = oi as usize;
                let s = &seq[offs[oi]..offs[oi + 1]];
                if p.flip[oi] {
                    revcomp(s)
                } else {
                    s.to_vec()
                }
            })
            .collect();
        let refs: Vec<&[u8]> = cl.iter().map(Vec::as_slice).collect();
        let anchors: Vec<u32> = p.order.iter().map(|&oi| p.anchor[oi as usize]).collect();
        let enc = encode_clustered(&refs, &anchors, 8).unwrap();
        let expect: Vec<Vec<u8>> = refs.iter().map(|r| r.to_vec()).collect();
        assert_eq!(decode_clustered(&enc).unwrap(), expect);
    }

    #[test]
    fn rescue_roundtrip() {
        // The same mix the v2 codec round-trips must also round-trip under rescue.
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGTACGT",
            b"ACGTACGTACGT", // match
            b"ACGTAAGTACGT", // 1 mismatch
            b"ACGTACGTNCGT", // mismatch incl. N
            b"TTTTGGGGCCCC", // literal
            b"",             // empty
            b"",             // match (empty == empty)
        ];
        let enc = encode_clustered_rescue(&reads, &vec![0u32; reads.len()], 4).expect("encode");
        let dec = decode_clustered_rescue(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn rescue_attaches_to_earlier_contig() {
        // Two unrelated references interleaved: A, B, A, B, A. After B seeds the
        // "current" contig, the v2 codec would strand the next A as a LITERAL;
        // the rescue index lets it re-attach to the earlier A contig. This asserts
        // the round-trip; op_stats_rescue reports the recovered literal.
        let a = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGC";
        let b = b"TTAGGCCATTACAGGTACCATGACATTGGACATTACAGGTTCAAGT";
        let reads: Vec<&[u8]> = vec![&a[..], &b[..], &a[..], &b[..], &a[..]];
        let anchors = vec![0u32; reads.len()];
        let enc = encode_clustered_rescue(&reads, &anchors, 6).expect("encode");
        let dec = decode_clustered_rescue(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
        // The third A (index 2) should be rescued onto A's contig, not stranded:
        // only the first A and the first B are literals.
        let st = op_stats_rescue(&reads, &anchors);
        assert_eq!(st.literals, 2, "expected only the two seeds to be literals");
    }

    proptest::proptest! {
        #[test]
        fn rescue_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..30),
                0..50)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let enc = encode_clustered_rescue(&refs, &vec![0u32; refs.len()], 4).expect("encode");
            let dec = decode_clustered_rescue(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, reads);
        }

        #[test]
        fn clustered_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..30),
                0..50)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let enc = encode_clustered(&refs, &vec![0u32; refs.len()], 4).expect("encode");
            let dec = decode_clustered(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, reads);
        }

        #[test]
        fn plan_is_valid_permutation(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..40),
                0..50),
            k in 1usize..=20,
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            let p = plan(&lens, &seq, k);
            let mut sorted = p.order.clone();
            sorted.sort_unstable();
            let expect: Vec<u32> = (0..reads.len() as u32).collect();
            proptest::prop_assert_eq!(sorted, expect);
            proptest::prop_assert_eq!(p.flip.len(), reads.len());
        }
    }
}
