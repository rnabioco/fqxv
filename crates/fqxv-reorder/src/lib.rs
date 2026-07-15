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
use std::hash::{BuildHasherDefault, Hasher};

use fqxv_bytes::{unzigzag, write_varint, zigzag};
use rayon::prelude::*;
use thiserror::Error;

/// A minimal integer hasher for the assembly maps. Their keys are already
/// well-mixed — 2-bit-packed k-mers and dense contig ids — so a single
/// multiplicative (Fibonacci) mix beats SipHash on the ~10^8 inserts/probes the
/// global assembler drives, and it is the throughput bottleneck of the v4 encode
/// path. Byte-output-preserving: these maps are only ever probed by key (never
/// iterated), and callers sort any candidate set deterministically, so the hash
/// choice cannot change the encoded stream.
#[derive(Default)]
struct IntHasher(u64);

impl Hasher for IntHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write_u64(&mut self, v: u64) {
        self.0 = v.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
    #[inline]
    fn write_u32(&mut self, v: u32) {
        self.0 = u64::from(v).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        // Only u32/u64 keys are hashed in this crate; keep a correct fallback so
        // the impl is total regardless of future key types.
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(b);
        }
        self.0 = self.0.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

/// A `HashMap` over integer keys using [`IntHasher`].
type IntMap<K, V> = HashMap<K, V, BuildHasherDefault<IntHasher>>;

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
    index: IntMap<u64, (u32, u32)>,
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
    let mism: Vec<usize> = (0..overlap)
        .filter(|&j| cur[j] != contig[off + j].base)
        .collect();
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
                    best = Some(Placement {
                        ci,
                        off,
                        overlap,
                        mism,
                    });
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
        self.contigs
            .push(cur.iter().map(|&b| seed_column(b)).collect());
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
        &ops_c,
        &cref_c,
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

/// Decode a clustered sequence block written by [`encode_clustered`] (version 2),
/// [`encode_clustered_rescue`] (version 3), or [`encode_global_block`] (version
/// 4), dispatching on the leading version byte. A version-4 block references the
/// shared frozen [`GlobalReference`], so `reference` must be `Some` for it;
/// versions 2/3 are self-contained and ignore it. Blocks may mix versions freely
/// within one archive.
pub fn decode_clustered_any(
    src: &[u8],
    reference: Option<&GlobalReference>,
) -> Result<Vec<Vec<u8>>> {
    match src.first() {
        Some(2) => decode_clustered(src),
        Some(3) => decode_clustered_rescue(src),
        Some(4) => {
            let r = reference.ok_or(Error::Malformed("version-4 block without reference"))?;
            decode_global_block(src, r)
        }
        _ => Err(Error::Malformed("unsupported version")),
    }
}

/// Back-compat shim: dispatch a version-2/3 block with no shared reference.
/// Equivalent to [`decode_clustered_any`] with `None`; version-4 blocks error.
pub fn decode_clustered_auto(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    decode_clustered_any(src, None)
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

// --- global-reference contig-assembly codec (prototype, version 4) -----------
//
// The v3 codec keeps assembly BLOCK-LOCAL: its multi-contig `Assembler` resets
// at every 256Ki-read block, so cross-block overlaps are lost, and enlarging the
// block only trades that gain against an exploding per-read `cref` recency
// back-reference over the growing contig set (see issue #52). v4 inverts the
// structure, SPRING-style: assemble ONE global reference over all clustered
// reads, freeze its final consensus, store it once (context-coded via
// `fqxv_seq`, deduplicated by construction), and code every read as a *position*
// on that shared reference — `(contig_id, offset, few mismatches)` — with the
// contig id DELTA-coded in clustered order rather than a global recency
// back-reference. Because clustering keeps same-contig reads adjacent, the id
// delta is mostly zero with rare jumps for reads a k-mer rescued onto a far
// contig; that is the lever that kills the `cref` blowup.
//
// Unlike v2/v3 a v4 block is NOT self-contained: it references the frozen
// global reference, which lives once at the whole-file level. Encoding is a
// two-pass whole-file mode: [`assemble_global`] builds+freezes the reference and
// the per-read placements, then [`encode_global_block`] codes each (parallel)
// block against the frozen reference. [`decode_global_block`] replays reads
// against the same reference — no vote/consensus reconstruction needed, so
// decode is a straight slice-and-patch.

/// The frozen global reference produced by [`assemble_global`]: the final
/// plurality-consensus bytes of every contig, concatenated, with per-contig
/// offsets. Reads are coded as positions on it and decoded by slicing it.
#[derive(Debug, Default, Clone)]
pub struct GlobalReference {
    /// Concatenated final consensus bytes of all contigs.
    seq: Vec<u8>,
    /// Byte offset of each contig in `seq`; `offs.len() == n_contigs + 1`.
    offs: Vec<usize>,
}

impl GlobalReference {
    /// Number of contigs in the reference.
    #[must_use]
    pub fn n_contigs(&self) -> usize {
        self.offs.len().saturating_sub(1)
    }

    /// Total reference bytes (the from-scratch content, stored once).
    #[must_use]
    pub fn total_bases(&self) -> usize {
        self.seq.len()
    }

    /// Consensus bytes of contig `ci`.
    fn contig(&self, ci: usize) -> &[u8] {
        &self.seq[self.offs[ci]..self.offs[ci + 1]]
    }

    /// The concatenated consensus of every contig — the exact bytes [`encode`]
    /// codes. Exposed for analysis (e.g. comparing the order-k coder against a
    /// long-range compressor on the raw reference).
    ///
    /// [`encode`]: GlobalReference::encode
    #[must_use]
    pub fn raw_bases(&self) -> &[u8] {
        &self.seq
    }

    /// Serialize the reference: contig count, then the concatenated consensus
    /// context-coded by [`fqxv_seq`] with per-contig lengths (so contigs are
    /// deduplicated and modeled as sequence, not stored raw). The reference is
    /// coded once for the whole file, so it is worth an aggressive hashed
    /// high-order tier (`hash_order`/`hash_bits`, as in [`fqxv_seq::encode_hashed`]);
    /// pass `hash_order == 0` for the plain dense order-`seq_order` model.
    pub fn encode(&self, seq_order: usize, hash_order: usize, hash_bits: u32) -> Result<Vec<u8>> {
        let lens: Vec<u32> = (0..self.n_contigs())
            .map(|c| (self.offs[c + 1] - self.offs[c]) as u32)
            .collect();
        let coded = fqxv_seq::encode_hashed(&lens, &self.seq, seq_order, hash_order, hash_bits)?;
        let mut out = Vec::new();
        write_varint(&mut out, self.n_contigs() as u64);
        write_varint(&mut out, coded.len() as u64);
        out.extend_from_slice(&coded);
        Ok(out)
    }

    /// Per-contig consensus lengths, in contig order. Together with
    /// [`raw_bases`](GlobalReference::raw_bases) these fully describe the
    /// reference, so an external coder can compress the bases and round-trip via
    /// [`from_lens_seq`](GlobalReference::from_lens_seq).
    #[must_use]
    pub fn contig_lens(&self) -> Vec<u32> {
        (0..self.n_contigs())
            .map(|c| (self.offs[c + 1] - self.offs[c]) as u32)
            .collect()
    }

    /// Rebuild a reference from per-contig lengths and the concatenated
    /// consensus (the inverse of [`contig_lens`](GlobalReference::contig_lens) +
    /// [`raw_bases`](GlobalReference::raw_bases)). Errors if the lengths do not
    /// sum to `seq.len()`.
    pub fn from_lens_seq(lens: &[u32], seq: Vec<u8>) -> Result<GlobalReference> {
        let mut offs = Vec::with_capacity(lens.len() + 1);
        let mut acc = 0usize;
        offs.push(0);
        for &l in lens {
            acc += l as usize;
            offs.push(acc);
        }
        if acc != seq.len() {
            return Err(Error::Malformed("reference length disagreement"));
        }
        Ok(GlobalReference { seq, offs })
    }

    /// Reverse of [`GlobalReference::encode`].
    pub fn decode(src: &[u8]) -> Result<GlobalReference> {
        let mut r = Cursor::new(src);
        let n = r.varint()? as usize;
        let coded = r.take_stream()?;
        let (lens, seq) = fqxv_seq::decode(coded)?;
        if lens.len() != n {
            return Err(Error::Malformed("reference contig count mismatch"));
        }
        let mut offs = Vec::with_capacity(n + 1);
        let mut acc = 0usize;
        offs.push(0);
        for l in &lens {
            acc += *l as usize;
            offs.push(acc);
        }
        if acc != seq.len() {
            return Err(Error::Malformed("reference length disagreement"));
        }
        Ok(GlobalReference { seq, offs })
    }

    /// Block-parallel variant of [`encode`](GlobalReference::encode): split the
    /// contigs into up to `n_blocks` contiguous groups (by contig index — fixed,
    /// so the output is byte-identical regardless of thread count) and code each
    /// group's `(lens, bases)` with a plain order-`seq_order` [`fqxv_seq`] model
    /// **in parallel**. Far faster than one whole-reference pass (and than xz) at
    /// a small ratio cost from the per-block context resets. Frame:
    /// `[varint n_blocks]` then, per block, `[varint n_contigs][varint len][coded]`.
    pub fn encode_blocked(&self, seq_order: usize, n_blocks: usize) -> Result<Vec<u8>> {
        let nc = self.n_contigs();
        let per = nc.div_ceil(n_blocks.clamp(1, nc.max(1)));
        let bounds: Vec<(usize, usize)> = (0..nc)
            .step_by(per.max(1))
            .map(|s| (s, (s + per).min(nc)))
            .collect();
        let coded: Vec<Vec<u8>> = bounds
            .par_iter()
            .map(|&(s, e)| -> Result<Vec<u8>> {
                let lens: Vec<u32> = (s..e)
                    .map(|c| (self.offs[c + 1] - self.offs[c]) as u32)
                    .collect();
                Ok(fqxv_seq::encode(
                    &lens,
                    &self.seq[self.offs[s]..self.offs[e]],
                    seq_order,
                )?)
            })
            .collect::<Result<_>>()?;
        let mut out = Vec::new();
        write_varint(&mut out, bounds.len() as u64);
        for (&(s, e), c) in bounds.iter().zip(&coded) {
            write_varint(&mut out, (e - s) as u64);
            write_varint(&mut out, c.len() as u64);
            out.extend_from_slice(c);
        }
        Ok(out)
    }

    /// Reverse of [`encode_blocked`](GlobalReference::encode_blocked).
    pub fn decode_blocked(src: &[u8]) -> Result<GlobalReference> {
        let mut r = Cursor::new(src);
        let nb = r.varint()? as usize;
        let mut blocks: Vec<(usize, &[u8])> = Vec::with_capacity(nb);
        for _ in 0..nb {
            let ncb = r.varint()? as usize;
            blocks.push((ncb, r.take_stream()?));
        }
        let decoded: Vec<(Vec<u32>, Vec<u8>)> = blocks
            .par_iter()
            .map(|&(ncb, coded)| -> Result<(Vec<u32>, Vec<u8>)> {
                let (lens, seq) = fqxv_seq::decode(coded)?;
                if lens.len() != ncb {
                    return Err(Error::Malformed("blocked reference contig count mismatch"));
                }
                Ok((lens, seq))
            })
            .collect::<Result<_>>()?;
        let mut lens = Vec::new();
        let mut seq = Vec::new();
        for (bl, bs) in decoded {
            lens.extend_from_slice(&bl);
            seq.extend_from_slice(&bs);
        }
        Self::from_lens_seq(&lens, seq)
    }
}

/// Where one read sits on the frozen reference: contig `ci`, starting at column
/// `off`. The read length (hence overlap) comes from the read itself, so this is
/// all the placement state a read needs.
#[derive(Debug, Clone, Copy, Default)]
pub struct Place4 {
    /// Contig index in the [`GlobalReference`].
    pub ci: u32,
    /// Start column of the read on that contig.
    pub off: u32,
}

/// Pass 1 of the v4 codec: assemble ALL clustered reads into one global set of
/// contigs (the multi-contig [`Assembler`], never reset), freeze the final
/// consensus into a [`GlobalReference`], and record each read's placement.
///
/// Exact duplicates of the previous read are not re-folded (they inherit the
/// previous read's placement), matching v3's `MATCH` short-circuit so the
/// reference structure is the global analogue of v3's per-block contigs. Every
/// read gets a valid `(ci, off)` so a read that lands at a parallel-block
/// boundary in pass 2 still has a reference position even when it can't be a
/// block-local `MATCH`. Deterministic: a sequential fold over the deterministic
/// clustered order.
#[must_use]
pub fn assemble_global(reads: &[&[u8]], anchors: &[u32]) -> (GlobalReference, Vec<Place4>) {
    assemble_window(reads, anchors)
}

/// The serial greedy fold over one window of reads: place each read on the
/// growing multi-contig assembly (or seed a new contig), then freeze the
/// consensus. Contig ids in the returned placements are local to this window.
fn assemble_window(reads: &[&[u8]], anchors: &[u32]) -> (GlobalReference, Vec<Place4>) {
    let mut asm = Assembler::default();
    let mut places: Vec<Place4> = Vec::with_capacity(reads.len());
    for (i, &cur) in reads.iter().enumerate() {
        if i > 0 && cur == reads[i - 1] {
            places.push(places[i - 1]);
            continue;
        }
        match asm.place(cur, anchors[i]) {
            Some(p) => {
                places.push(Place4 {
                    ci: p.ci as u32,
                    off: p.off as u32,
                });
                asm.commit(p.ci, cur, p.off, p.overlap);
            }
            None => {
                let ci = asm.contigs.len();
                asm.seed(cur, anchors[i]);
                places.push(Place4 {
                    ci: ci as u32,
                    off: 0,
                });
            }
        }
    }
    // Freeze: concatenate every contig's final consensus byte.
    let total: usize = asm.contigs.iter().map(Vec::len).sum();
    let mut seq = Vec::with_capacity(total);
    let mut offs = Vec::with_capacity(asm.contigs.len() + 1);
    offs.push(0);
    for c in &asm.contigs {
        for col in c {
            seq.push(col.base);
        }
        offs.push(seq.len());
    }
    (GlobalReference { seq, offs }, places)
}

/// Parallel windowed assembly: split the clustered reads into `n_windows`
/// contiguous windows (by read index — fixed, so the result is byte-identical
/// regardless of thread count), assemble each **in parallel** with the serial
/// [`assemble_window`], then concatenate their frozen references (remapping each
/// window's local contig ids by a running offset). Windowing costs cross-window
/// deduplication, but a following [`merge_reference`] recovers most of it by
/// chaining duplicate contigs — so this is a near-ratio-neutral speedup of the
/// otherwise-serial [`assemble_global`] fold. `n_windows == 1` reproduces
/// [`assemble_global`] exactly.
#[must_use]
pub fn assemble_global_windowed(
    reads: &[&[u8]],
    anchors: &[u32],
    n_windows: usize,
) -> (GlobalReference, Vec<Place4>) {
    let n = reads.len();
    if n == 0 {
        return (
            GlobalReference {
                seq: Vec::new(),
                offs: vec![0],
            },
            Vec::new(),
        );
    }
    let per = n.div_ceil(n_windows.clamp(1, n));
    let ranges: Vec<(usize, usize)> = (0..n)
        .step_by(per.max(1))
        .map(|s| (s, (s + per).min(n)))
        .collect();
    let windows: Vec<(GlobalReference, Vec<Place4>)> = ranges
        .par_iter()
        .map(|&(s, e)| assemble_window(&reads[s..e], &anchors[s..e]))
        .collect();

    let mut seq = Vec::new();
    let mut offs = vec![0usize];
    let mut places = Vec::with_capacity(n);
    let mut contig_off = 0u32;
    for (gref, wplaces) in windows {
        seq.extend_from_slice(&gref.seq);
        for w in 1..gref.offs.len() {
            offs.push(offs[offs.len() - 1] + (gref.offs[w] - gref.offs[w - 1]));
        }
        for p in wplaces {
            places.push(Place4 {
                ci: p.ci + contig_off,
                off: p.off,
            });
        }
        contig_off += gref.n_contigs() as u32;
    }
    (GlobalReference { seq, offs }, places)
}

// --- PROTOTYPE: SPRING-style bidirectional assembly (measurement only) --------
//
// fqxv's [`assemble_global`] grows each contig RIGHTWARD only: `place` requires a
// non-negative offset, so a read that overlaps a contig's *left* end (it would
// start before the contig) can't attach and seeds a near-duplicate new contig.
// SPRING (`reorder.h`) instead extends a contig rightward until dry, then
// LEFTWARD (its `left_search`), and only then starts a new contig — which is why
// its reference has fewer/longer contigs (~55M bases vs fqxv's 88.9M). This
// prototype adds that leftward extension to measure the contig/base-count
// reduction. It tracks NO per-read placements (not container-wired) — it only
// freezes the assembled contigs so the reference base count can be compared.

/// A contig that can grow at BOTH ends. Columns live in a `VecDeque` (cheap front
/// growth); `left` counts columns prepended so the k-mer index, keyed by
/// ORIGIN-relative position (`deque_index - left`), survives prepends.
type Contig = std::collections::VecDeque<Column>;

#[derive(Default)]
struct BidirAssembler {
    contigs: Vec<Contig>,
    left: Vec<usize>,
    ref_anchors: Vec<u32>,
    index: IntMap<u64, (u32, i32)>,
}

/// k-mer at deque position `start` of a contig, or `None` (off-end / non-ACGT).
#[inline]
fn contig_kmer(dq: &Contig, start: usize) -> Option<u64> {
    if start + RESCUE_K > dq.len() {
        return None;
    }
    let mut v = 0u64;
    for j in 0..RESCUE_K {
        let c = code(dq[start + j].base);
        if c >= 4 {
            return None;
        }
        v = (v << 2) | u64::from(c);
    }
    Some(v)
}

/// Acceptance test for placing `cur` on `dq` at deque shift `s` (the deque index
/// of `cur[0]`; may be negative — leftward). Returns `(overlap, mismatches)` when
/// cheaper than a literal, mirroring [`try_place`]'s thresholds.
fn try_place_shift(dq: &Contig, cur: &[u8], s: i64) -> Option<(usize, usize)> {
    let len = dq.len() as i64;
    let ov_start = s.max(0);
    let ov_end = (s + cur.len() as i64).min(len);
    if ov_end <= ov_start {
        return None;
    }
    let overlap = (ov_end - ov_start) as usize;
    if overlap < MIN_CONTIG_OVERLAP.min(cur.len()) {
        return None;
    }
    let mut mism = 0usize;
    for c in ov_start..ov_end {
        let r = (c - s) as usize;
        if cur[r] != dq[c as usize].base {
            mism += 1;
        }
    }
    let novel = cur.len() - overlap;
    (mism <= overlap / 4 && novel + mism * 2 < cur.len()).then_some((overlap, mism))
}

impl BidirAssembler {
    /// Index k-mers for deque positions `[from, to)` of contig `ci`, keyed by
    /// origin-relative position (survives future prepends).
    fn index_range(&mut self, ci: usize, from: usize, to: usize) {
        let Self {
            contigs,
            left,
            index,
            ..
        } = self;
        let dq = &contigs[ci];
        let l = left[ci] as i32;
        let hi = to.min(dq.len().saturating_sub(RESCUE_K - 1));
        for start in from..hi {
            if let Some(v) = contig_kmer(dq, start) {
                index.insert(v, (ci as u32, start as i32 - l));
            }
        }
    }

    /// Best placement of `cur` across contigs as `(ci, shift)` (shift may be
    /// negative), or `None` to seed. Same candidate sources as [`Assembler::place`]
    /// (anchor fast-path on the last contig + sampled read k-mers) plus negative
    /// shifts. Deterministic: candidates deduped and scored by (mism, recency).
    fn place(&self, cur: &[u8], anchor: u32) -> Option<(usize, i64)> {
        if cur.is_empty() || self.contigs.is_empty() {
            return None;
        }
        let mut cands: Vec<(usize, i64)> = Vec::new();
        let last = self.contigs.len() - 1;
        // Anchor fast-path on the last contig (in deque frame: + left).
        let center = self.ref_anchors[last] as i64 - anchor as i64 + self.left[last] as i64;
        cands.push((last, center));
        let mut start = 0usize;
        while start + RESCUE_K <= cur.len() {
            if let Some(k) = kmer_at(cur, start, RESCUE_K) {
                if let Some(&(ci, orp)) = self.index.get(&k) {
                    let vpos = orp as i64 + self.left[ci as usize] as i64;
                    cands.push((ci as usize, vpos - start as i64));
                }
            }
            start += RESCUE_K;
        }
        cands.sort_unstable();
        cands.dedup();

        let mut best: Option<(usize, i64)> = None;
        let mut best_key = (usize::MAX, usize::MAX);
        for (ci, s) in cands {
            if let Some((_ov, mism)) = try_place_shift(&self.contigs[ci], cur, s) {
                let key = (mism, self.contigs.len() - 1 - ci);
                if key < best_key {
                    best_key = key;
                    best = Some((ci, s));
                }
            }
        }
        // SPRING-style EXHAUSTIVE shift search on the current (last) contig, as a
        // fallback when the k-mer index found no placement: try every shift within
        // ±readlen/2 of the anchor-implied centre and accept the best within the
        // mismatch budget. This catches reads that DO align to the current contig
        // but whose sampled k-mers all carried an error (so the index missed
        // them) — SPRING searches all shifts (`0..maxshift`) rather than only
        // k-mer-implied offsets, which is a chunk of why it places more reads.
        if best.is_none() {
            let last = self.contigs.len() - 1;
            let center = self.ref_anchors[last] as i64 - anchor as i64 + self.left[last] as i64;
            let maxshift = (cur.len() / 2) as i64;
            for s in (center - maxshift)..=(center + maxshift) {
                if let Some((_ov, mism)) = try_place_shift(&self.contigs[last], cur, s) {
                    if (mism, 0usize) < best_key {
                        best_key = (mism, 0);
                        best = Some((last, s));
                    }
                }
            }
        }
        best
    }

    /// Fold `cur` into contig `ci` at shift `s`, extending right (append) and/or
    /// left (prepend) as needed, then re-index the changed ends.
    fn commit(&mut self, ci: usize, cur: &[u8], s: i64) {
        if s >= 0 {
            let off = s as usize;
            let old_len = self.contigs[ci].len();
            let ov = cur.len().min(old_len - off);
            for k in 0..ov {
                cast_vote(&mut self.contigs[ci][off + k], cur[k]);
            }
            for &b in &cur[ov..] {
                self.contigs[ci].push_back(seed_column(b));
            }
            if self.contigs[ci].len() > old_len {
                self.index_range(
                    ci,
                    old_len.saturating_sub(RESCUE_K - 1),
                    self.contigs[ci].len(),
                );
            }
        } else {
            let p = (-s) as usize;
            // Prepend cur[0..p] (front-most is cur[0]).
            for &b in cur[..p].iter().rev() {
                self.contigs[ci].push_front(seed_column(b));
            }
            self.left[ci] += p;
            // cur[p..] now aligns to deque[p..].
            let ov = (cur.len() - p).min(self.contigs[ci].len() - p);
            for k in 0..ov {
                cast_vote(&mut self.contigs[ci][p + k], cur[p + k]);
            }
            for &b in &cur[p + ov..] {
                self.contigs[ci].push_back(seed_column(b));
            }
            // Re-index the whole (usually short) contig — positions in deque frame
            // shifted, though origin-relative keys for old columns are unchanged;
            // re-indexing the front join + tail is enough, done generously here.
            let len = self.contigs[ci].len();
            self.index_range(ci, 0, (p + RESCUE_K).min(len));
            if ov < cur.len() - p {
                self.index_range(ci, len.saturating_sub(RESCUE_K - 1), len);
            }
        }
    }

    fn seed(&mut self, cur: &[u8], anchor: u32) {
        let ci = self.contigs.len();
        self.contigs
            .push(cur.iter().map(|&b| seed_column(b)).collect());
        self.left.push(0);
        self.ref_anchors.push(anchor);
        self.index_range(ci, 0, cur.len());
    }
}

/// PROTOTYPE (measurement only): SPRING-style BIDIRECTIONAL global assembly.
/// Same as [`assemble_global`] but reads that overlap a contig's LEFT end prepend
/// to it instead of seeding a near-duplicate contig. Returns only the frozen
/// reference (no placements) — for measuring the contig/base-count reduction.
#[must_use]
pub fn assemble_global_bidir(reads: &[&[u8]], anchors: &[u32]) -> GlobalReference {
    let mut asm = BidirAssembler::default();
    for (i, &cur) in reads.iter().enumerate() {
        if cur.is_empty() {
            continue;
        }
        if i > 0 && cur == reads[i - 1] {
            continue; // exact duplicate of previous — inherits its placement
        }
        match asm.place(cur, anchors[i]) {
            Some((ci, s)) => asm.commit(ci, cur, s),
            None => asm.seed(cur, anchors[i]),
        }
    }
    let mut seq = Vec::new();
    let mut offs = vec![0usize];
    for dq in &asm.contigs {
        for col in dq {
            seq.push(col.base);
        }
        offs.push(seq.len());
    }
    GlobalReference { seq, offs }
}

// --- PROTOTYPE: fully-integrated SPRING-style reference-driven assembly ---------
//
// fqxv's assembler is READ-DRIVEN: it walks reads in clustered order and pushes
// each onto its best existing contig, so deep-coverage reads that don't attach
// seed near-duplicate contigs (88.9M reference bases vs SPRING's ~55M). SPRING is
// REFERENCE-DRIVEN (`reorder.h`): it seeds one contig, then greedily PULLS IN
// every remaining read that overlaps the current reference — extending RIGHTWARD
// until dry, then LEFTWARD (its `left_search`) until dry — before starting the
// next contig, so a whole region collapses into ONE contig. This is the full
// integrated version of that: bidirectional, max-overlap ("all-shift") read
// selection so the end sweeps slowly and skips nothing, pulling every matching
// read. Reads are already minimizer-oriented, so reverse-complement search is
// omitted (a read that is RC relative to a cross-cluster contig is not pulled — a
// known undercount vs SPRING). Returns only the frozen reference (no placements);
// this measures the base-count reduction from the assembly STRATEGY.

/// Read-index fanout, sampling, and per-end scan window for [`assemble_global_refdriven`].
const RD_FANOUT: usize = 24;
const RD_SAMPLE: usize = 4;
const RD_WINDOW: usize = 200;

/// PROTOTYPE (measurement): the reference-driven bidirectional assembler above.
///
/// **Does NOT converge yet** (measured on 4M NovaSeq — see the design issue): all
/// tried extension-selection rules produce MORE reference bases than fqxv's
/// baseline, not fewer, because the pieces are not tuned together:
/// - max-overlap (min advance): contigs *creep* one base per read and stall at
///   ~200 bp (~+84% bases, 887k contigs);
/// - max-reach (max advance, current): chains through spurious short overlaps into
///   *chimeric* contigs that duplicate regions (~+328% bases, 1.2M contigs).
/// Reaching SPRING's ~55 MB needs coupled tuning (consensus quality, overlap +
/// identity thresholds, the sliding-window/shift bookkeeping, chimera guards) with
/// instrumentation — a real assembler-engineering project, not this prototype. The
/// separately-measured leftward lever ([`assemble_global_bidir`], −8.8% bases) is
/// the tractable win.
#[must_use]
pub fn assemble_global_refdriven(reads: &[&[u8]], _anchors: &[u32]) -> GlobalReference {
    const ORIENT_BIT: u32 = 1 << 31;
    let n = reads.len();
    // Reverse complements, so a read that is RC relative to a cross-cluster contig
    // can still be pulled in (SPRING searches both orientations at each step).
    let rc_reads: Vec<Vec<u8>> = reads.iter().map(|r| revcomp(r)).collect();
    // The oriented sequence for a candidate `(rid, orient)`.
    let oriented = |rid: usize, orient: u32| -> &[u8] {
        if orient == 0 {
            reads[rid]
        } else {
            &rc_reads[rid]
        }
    };
    // Read index (built once): k-mer -> [(read id, read pos | orient-bit)]. Both
    // the forward read and its reverse complement are indexed.
    let mut index: IntMap<u64, Vec<(u32, u32)>> = IntMap::default();
    for rid in 0..n {
        for (orient, seq) in [(0u32, reads[rid]), (ORIENT_BIT, rc_reads[rid].as_slice())] {
            let mut p = 0usize;
            while p + RESCUE_K <= seq.len() {
                if let Some(k) = kmer_at(seq, p, RESCUE_K) {
                    let e = index.entry(k).or_default();
                    if e.len() < RD_FANOUT {
                        e.push((rid as u32, p as u32 | orient));
                    }
                }
                p += RD_SAMPLE;
            }
        }
    }

    let mut used = vec![false; n];
    let mut contigs: Vec<Contig> = Vec::new();
    // Per-step candidate map (rid -> (best read_start, orient)), reused.
    let mut cand: IntMap<u32, (i64, u32)> = IntMap::default();

    for seed in 0..n {
        if used[seed] || reads[seed].is_empty() {
            continue;
        }
        used[seed] = true;
        let mut c: Contig = reads[seed].iter().map(|&b| seed_column(b)).collect();

        // RIGHTWARD sweep: pull the max-overlap read that extends the right end.
        loop {
            let l = c.len();
            cand.clear();
            let lo = l.saturating_sub(RD_WINDOW);
            let hi = l.saturating_sub(RESCUE_K);
            for p in lo..=hi {
                let Some(k) = contig_kmer(&c, p) else {
                    continue;
                };
                let Some(cs) = index.get(&k) else { continue };
                for &(rid, packed) in cs {
                    if used[rid as usize] {
                        continue;
                    }
                    let orient = packed & ORIENT_BIT;
                    let rstart = p as i64 - (packed & !ORIENT_BIT) as i64;
                    // Keep the smallest read_start (largest overlap) per read.
                    cand.entry(rid)
                        .and_modify(|e| {
                            if rstart < e.0 {
                                *e = (rstart, orient);
                            }
                        })
                        .or_insert((rstart, orient));
                }
            }
            // Among valid extenders, take the one reaching FURTHEST right (max
            // `rs + read_len`) so the contig grows by the largest reliable step and
            // a whole region collapses into one contig, rather than creeping one
            // base at a time. `reach` is the new contig length if this read is
            // taken; the hamming gate keeps the overlap reliable.
            let mut best: Option<(u32, usize, usize, u32)> = None; // (rid, rs, reach, orient)
            for (&rid, &(rstart, orient)) in &cand {
                if rstart < 0 {
                    continue;
                }
                let rs = rstart as usize;
                let r = oriented(rid as usize, orient);
                let reach = rs + r.len();
                if reach <= l {
                    continue; // no rightward extension
                }
                let overlap = l - rs;
                // read[j] vs contig column rs+j.
                if overlap < MIN_CONTIG_OVERLAP || !hamming_ok(r, &c, 0, rs, overlap) {
                    continue;
                }
                if best.is_none_or(|(brid, _, br, _)| reach > br || (reach == br && rid < brid)) {
                    best = Some((rid, rs, reach, orient));
                }
            }
            let Some((rid, rs, _, orient)) = best else {
                break;
            };
            used[rid as usize] = true;
            let r = oriented(rid as usize, orient);
            for cc in rs..l {
                cast_vote(&mut c[cc], r[cc - rs]);
            }
            for &b in &r[l - rs..] {
                c.push_back(seed_column(b));
            }
        }

        // LEFTWARD sweep: pull the max-overlap read that extends the left end.
        loop {
            let l = c.len();
            cand.clear();
            let hi = RD_WINDOW.min(l.saturating_sub(RESCUE_K));
            for p in 0..=hi {
                let Some(k) = contig_kmer(&c, p) else {
                    continue;
                };
                let Some(cs) = index.get(&k) else { continue };
                for &(rid, packed) in cs {
                    if used[rid as usize] {
                        continue;
                    }
                    let orient = packed & ORIENT_BIT;
                    let rstart = p as i64 - (packed & !ORIENT_BIT) as i64;
                    // Keep the smallest (most-negative) read_start = reaches furthest left.
                    cand.entry(rid)
                        .and_modify(|e| {
                            if rstart < e.0 {
                                *e = (rstart, orient);
                            }
                        })
                        .or_insert((rstart, orient));
                }
            }
            // Take the read reaching FURTHEST left (largest prepend), the mirror of
            // the rightward furthest-reach rule.
            let mut best: Option<(u32, usize, usize, u32)> = None; // (rid, prepend, pfx, orient)
            for (&rid, &(rstart, orient)) in &cand {
                if rstart >= 0 {
                    continue; // must extend left
                }
                let r = oriented(rid as usize, orient);
                let oend = (rstart + r.len() as i64).min(l as i64);
                if oend <= 0 {
                    continue;
                }
                let overlap = oend as usize; // overlap region is columns [0, overlap)
                let pfx = (-rstart) as usize;
                // read[pfx+j] vs contig column j.
                if overlap < MIN_CONTIG_OVERLAP || !hamming_ok(r, &c, pfx, 0, overlap) {
                    continue;
                }
                if best.is_none_or(|(brid, _, bp, _)| pfx > bp || (pfx == bp && rid < brid)) {
                    best = Some((rid, pfx, pfx, orient));
                }
            }
            let Some((rid, pfx, _, orient)) = best else {
                break;
            };
            used[rid as usize] = true;
            let r = oriented(rid as usize, orient);
            for &b in r[..pfx].iter().rev() {
                c.push_front(seed_column(b));
            }
            // The read now aligns at deque index 0; vote read[pfx..] into c[pfx..].
            let l2 = c.len();
            let ov = (r.len() - pfx).min(l2 - pfx);
            for k in 0..ov {
                cast_vote(&mut c[pfx + k], r[pfx + k]);
            }
            for &b in &r[pfx + ov..] {
                c.push_back(seed_column(b));
            }
        }

        contigs.push(c);
    }

    let mut seq = Vec::new();
    let mut offs = vec![0usize];
    for c in &contigs {
        for col in c {
            seq.push(col.base);
        }
        offs.push(seq.len());
    }
    GlobalReference { seq, offs }
}

/// Hamming check over `overlap` columns: `read[read_off + j]` vs contig column
/// `col_off + j`, accepting up to `overlap / 4` mismatches (mirrors [`try_place`]).
#[inline]
fn hamming_ok(read: &[u8], c: &Contig, read_off: usize, col_off: usize, overlap: usize) -> bool {
    let budget = overlap / 4;
    let mut mism = 0usize;
    for j in 0..overlap {
        if read[read_off + j] != c[col_off + j].base {
            mism += 1;
            if mism > budget {
                return false;
            }
        }
    }
    true
}

/// Pass 2 of the v4 codec: code one block of clustered reads as positions on the
/// frozen `reference`, using the placements from [`assemble_global`]. Each read
/// is `MATCH` (byte-identical to the block-previous read) or `CONTIG` — a
/// delta-coded contig id, a per-contig delta-coded offset, and the substitutions
/// versus the frozen consensus. `places` is the slice for this block (same range
/// as `reads`). Byte-exactly reversible by [`decode_global_block`] given the
/// same reference.
pub fn encode_global_block(
    reads: &[&[u8]],
    places: &[Place4],
    reference: &GlobalReference,
) -> Result<Vec<u8>> {
    let mut ops = Vec::with_capacity(reads.len());
    let (mut cid, mut offdelta, mut slen) = (Vec::new(), Vec::new(), Vec::new());
    let (mut nmis, mut pos, mut subs, mut tail) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());

    let mut prev_cid: i64 = 0;
    // Per-contig previous offset (delta-coded within a contig). Bounded by the
    // distinct contigs a single block references, so a map stays small.
    let mut last_off: IntMap<u32, i64> = IntMap::default();

    for (i, &cur) in reads.iter().enumerate() {
        if i > 0 && cur == reads[i - 1] {
            ops.push(OP_MATCH);
            continue;
        }
        ops.push(OP_CONTIG);
        let p = places[i];
        let ci = p.ci as i64;
        write_varint(&mut cid, zigzag(ci - prev_cid));
        prev_cid = ci;
        let off = p.off as usize;
        let lo = last_off.entry(p.ci).or_insert(0);
        write_varint(&mut offdelta, zigzag(off as i64 - *lo));
        *lo = off as i64;
        write_varint(&mut slen, cur.len() as u64);

        let contig = reference.contig(p.ci as usize);
        let overlap = cur.len().min(contig.len().saturating_sub(off));
        let mism: Vec<usize> = (0..overlap)
            .filter(|&j| cur[j] != contig[off + j])
            .collect();
        write_varint(&mut nmis, mism.len() as u64);
        let mut last = 0usize;
        for &m in &mism {
            write_varint(&mut pos, (m - last) as u64);
            last = m;
            subs.push(cur[m]);
        }
        // On real short-read data every placed read fits within its frozen
        // contig, so `tail` stays empty; keep it as a safety valve for edge
        // cases (short reference slices) so the codec never loses bytes.
        tail.extend_from_slice(&cur[overlap..]);
    }

    let ops_c = fqxv_rans::encode(&ops, fqxv_rans::Order::One)?;
    let cid_c = fqxv_rans::encode(&cid, fqxv_rans::Order::Zero)?;
    let offdelta_c = fqxv_rans::encode(&offdelta, fqxv_rans::Order::Zero)?;
    let slen_c = fqxv_rans::encode(&slen, fqxv_rans::Order::Zero)?;
    let nmis_c = fqxv_rans::encode(&nmis, fqxv_rans::Order::Zero)?;
    let pos_c = fqxv_rans::encode(&pos, fqxv_rans::Order::Zero)?;
    let subs_c = fqxv_rans::encode(&subs, fqxv_rans::Order::One)?;
    let tail_c = fqxv_rans::encode(&tail, fqxv_rans::Order::One)?;

    let mut out = Vec::new();
    out.push(4u8); // version 4: global-reference layout
    write_varint(&mut out, reads.len() as u64);
    for s in [
        &ops_c,
        &cid_c,
        &offdelta_c,
        &slen_c,
        &nmis_c,
        &pos_c,
        &subs_c,
        &tail_c,
    ] {
        write_varint(&mut out, s.len() as u64);
        out.extend_from_slice(s);
    }
    Ok(out)
}

/// Decode a block written by [`encode_global_block`] against the same frozen
/// `reference`, returning the reads in clustered order. No consensus is rebuilt:
/// each read is a slice of the reference with its substitutions patched in.
pub fn decode_global_block(src: &[u8], reference: &GlobalReference) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != 4 {
        return Err(Error::Malformed("unsupported version"));
    }
    let n = r.varint()? as usize;
    let ops = fqxv_rans::decode(r.take_stream()?)?;
    let cid = fqxv_rans::decode(r.take_stream()?)?;
    let offdelta = fqxv_rans::decode(r.take_stream()?)?;
    let slen = fqxv_rans::decode(r.take_stream()?)?;
    let nmis = fqxv_rans::decode(r.take_stream()?)?;
    let pos = fqxv_rans::decode(r.take_stream()?)?;
    let subs = fqxv_rans::decode(r.take_stream()?)?;
    let tail = fqxv_rans::decode(r.take_stream()?)?;

    let mut c_cid = Cursor::new(&cid);
    let mut c_offdelta = Cursor::new(&offdelta);
    let mut c_slen = Cursor::new(&slen);
    let mut c_nmis = Cursor::new(&nmis);
    let mut c_pos = Cursor::new(&pos);
    let (mut subs_pos, mut tail_pos) = (0usize, 0usize);
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n.min(1 << 22));

    let mut prev_cid: i64 = 0;
    let mut last_off: IntMap<u32, i64> = IntMap::default();

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
            OP_CONTIG => {
                let ci_i = prev_cid + unzigzag(c_cid.varint()?);
                prev_cid = ci_i;
                let ci = u32::try_from(ci_i).map_err(|_| Error::Malformed("bad contig id"))?;
                if ci as usize >= reference.n_contigs() {
                    return Err(Error::Malformed("contig id out of range"));
                }
                let lo = last_off.entry(ci).or_insert(0);
                let off = usize::try_from(*lo + unzigzag(c_offdelta.varint()?))
                    .map_err(|_| Error::Malformed("bad contig offset"))?;
                *lo = off as i64;
                let cur_len = c_slen.varint()? as usize;
                let contig = reference.contig(ci as usize);
                if off > contig.len() {
                    return Err(Error::Malformed("contig offset past reference"));
                }
                let overlap = cur_len.min(contig.len() - off);
                let mut read = vec![0u8; cur_len];
                read[..overlap].copy_from_slice(&contig[off..off + overlap]);
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
                    *slot = *tail
                        .get(tail_pos)
                        .ok_or(Error::Malformed("tail underrun"))?;
                    tail_pos += 1;
                }
                reads.push(read);
            }
            _ => return Err(Error::Malformed("unknown op")),
        }
    }
    Ok(reads)
}

// --- overlap-merge assembler refinement (prototype) --------------------------
//
// The greedy [`assemble_global`] pass never compares contigs to EACH OTHER, so
// deep short-read data fragments into many contigs barely longer than one read
// (on 4M NovaSeq reads: 492K contigs averaging ~204 bp, ~1.4 reads each). The
// reference — which stores that content once — is then most of the v4 seq bytes.
// [`merge_reference`] is an overlap-layout refinement (OLC-lite): chain contigs
// whose suffix overlaps another contig's PREFIX into longer super-contigs, store
// the shared overlap once, and remap every read's placement onto the merged
// reference. It is format-transparent — [`encode_global_block`] recomputes each
// read's mismatches against whatever reference it is handed — so it is a pure
// encoder-side swap the decoder never sees. Deterministic.

/// k-mer length for detecting contig overlaps (matches the assembly minimizer).
const MERGE_K: usize = RESCUE_K;
/// Shortest contig-contig overlap worth merging.
const MIN_MERGE_OVL: usize = 24;
/// Index each contig's first `MERGE_PREFIX` bases; a successor's start must land
/// within here for the overlap to be found. Bounds the index to ~one short-read
/// worth of prefix per contig.
const MERGE_PREFIX: usize = 64;
/// Probe each contig's last `MERGE_SUFFIX` bases for overlaps into a successor.
const MERGE_SUFFIX: usize = 220;
/// Cap the candidates kept per k-mer so a repetitive k-mer can't blow up cost.
const MERGE_FANOUT: usize = 16;

/// Union-find root with path halving.
fn uf_find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize];
        x = parent[x as usize];
    }
    x
}

/// Tunable thresholds for [`merge_reference_with`]. [`MergeConfig::default`]
/// reproduces [`merge_reference`] byte-for-byte, so sweeping these is a pure
/// encoder-side experiment — the decoder never sees the reference shape.
#[derive(Debug, Clone, Copy)]
pub struct MergeConfig {
    /// Shortest contig-contig overlap worth merging.
    pub min_ovl: usize,
    /// Index each contig's first `prefix` bases as successor entry points.
    pub prefix: usize,
    /// Probe each contig's last `suffix` bases for overlaps into a successor.
    pub suffix: usize,
    /// Cap candidates kept per k-mer so a repetitive k-mer can't blow up cost.
    pub fanout: usize,
    /// Mismatch budget for an overlap is `overlap / mism_div` (larger = stricter).
    pub mism_div: usize,
}

impl Default for MergeConfig {
    fn default() -> Self {
        Self {
            min_ovl: MIN_MERGE_OVL,
            prefix: MERGE_PREFIX,
            suffix: MERGE_SUFFIX,
            fanout: MERGE_FANOUT,
            mism_div: 8,
        }
    }
}

/// Overlap-merge with default thresholds ([`MergeConfig::default`]). See
/// [`merge_reference_with`] for the full semantics.
#[must_use]
pub fn merge_reference(
    reads: &[&[u8]],
    reference: &GlobalReference,
    places: &[Place4],
) -> (GlobalReference, Vec<Place4>) {
    merge_reference_with(reads, reference, places, MergeConfig::default())
}

/// Number of contig chunks the merge k-mer index is built over in parallel.
/// The combined index is invariant to this (chunks are combined in contig order),
/// so it affects only parallelism, not the output.
const MERGE_INDEX_CHUNKS: usize = 64;

/// Build the prefix-k-mer index for [`merge_reference_with`], in parallel. Each
/// contig chunk builds a partial (fan-out-capped) map; the partials are combined
/// in contig order and re-capped, giving the exact same first-N entries per
/// k-mer as a serial build — so the result is deterministic regardless of thread
/// count. This is the merge's hottest step, so parallelizing it matters most.
fn build_merge_index(
    contigs: &[&[u8]],
    prefix: usize,
    fanout: usize,
) -> IntMap<u64, Vec<(u32, u32)>> {
    let nc = contigs.len();
    let chunk = nc.div_ceil(MERGE_INDEX_CHUNKS.clamp(1, nc.max(1))).max(1);
    let partials: Vec<IntMap<u64, Vec<(u32, u32)>>> = (0..nc)
        .step_by(chunk)
        .collect::<Vec<_>>()
        .par_iter()
        .map(|&start| {
            let end = (start + chunk).min(nc);
            let mut m: IntMap<u64, Vec<(u32, u32)>> = IntMap::default();
            for ci in start..end {
                let c = contigs[ci];
                let hi = c.len().min(prefix);
                let mut s = 0usize;
                while s + MERGE_K <= hi {
                    if let Some(code) = kmer_at(c, s, MERGE_K) {
                        let e = m.entry(code).or_default();
                        if e.len() < fanout {
                            e.push((ci as u32, s as u32));
                        }
                    }
                    s += 1;
                }
            }
            m
        })
        .collect();
    let mut index: IntMap<u64, Vec<(u32, u32)>> = IntMap::default();
    for part in partials {
        for (code, list) in part {
            let e = index.entry(code).or_default();
            for item in list {
                if e.len() < fanout {
                    e.push(item);
                } else {
                    break;
                }
            }
        }
    }
    index
}

/// Overlap-merge a greedy reference (see the module note): returns a new
/// `(reference, placements)` with fewer, longer contigs, usable by
/// [`encode_global_block`] unchanged. After chaining, the merged consensus is
/// RE-VOTED from every read at its remapped position, so overlap columns reflect
/// all contributing reads (not just the earliest contig's bytes) — that keeps the
/// per-read mismatch cost down. `reads` are the clustered, oriented reads that
/// produced `places` (read `i` has placement `places[i]`). Purely additive
/// refinement — never splits a contig, so every read keeps a valid placement.
/// `cfg` tunes the overlap-search thresholds ([`MergeConfig`]).
#[must_use]
pub fn merge_reference_with(
    reads: &[&[u8]],
    reference: &GlobalReference,
    places: &[Place4],
    cfg: MergeConfig,
) -> (GlobalReference, Vec<Place4>) {
    let nc = reference.n_contigs();
    if nc < 2 {
        return (reference.clone(), places.to_vec());
    }
    let contigs: Vec<&[u8]> = (0..nc).map(|c| reference.contig(c)).collect();

    // 1. Index each contig's PREFIX k-mers -> [(contig, pos)] (capped fan-out).
    // Built over contig CHUNKS in parallel and combined in contig order, so the
    // fan-out cap keeps the same first-N entries as a serial build — the combined
    // index is independent of the chunk count, hence of the thread count.
    let index = build_merge_index(&contigs, cfg.prefix, cfg.fanout);

    // 2. For each contig A, find its best successor B: A's suffix overlaps B's
    //    prefix (B starts at offset `s` inside A, overlap = A.len - s reaches A's
    //    end and matches B[0..overlap] within a small mismatch budget). Prefer the
    //    longest overlap, then fewest mismatches, then smallest ids (determinism).
    // Each contig's best successor depends only on the immutable `contigs` and
    // `index`, so compute them in parallel — this is the merge's hottest loop
    // (per-contig suffix probing + mismatch scans). `best_key` is a total order
    // ((MAX-ovl, mism, bi, s) minimised), so the winner — and the whole result —
    // is independent of thread count. `succ[ai] = (contig B, shift s)`.
    let succ: Vec<Option<(u32, u32)>> = (0..nc)
        .into_par_iter()
        .map(|ai| {
            let a = contigs[ai];
            if a.len() < MERGE_K {
                return None;
            }
            let lo = a.len().saturating_sub(cfg.suffix);
            let mut best_key = (usize::MAX, usize::MAX, usize::MAX, usize::MAX);
            let mut best: Option<(u32, u32)> = None;
            let mut pos_a = lo;
            while pos_a + MERGE_K <= a.len() {
                if let Some(code) = kmer_at(a, pos_a, MERGE_K) {
                    if let Some(list) = index.get(&code) {
                        for &(bi_u, pos_b_u) in list {
                            let bi = bi_u as usize;
                            if bi == ai {
                                continue;
                            }
                            let pos_b = pos_b_u as usize;
                            if pos_a < pos_b {
                                continue;
                            }
                            let s = pos_a - pos_b;
                            if s == 0 || s >= a.len() {
                                continue;
                            }
                            let ovl = a.len() - s;
                            let b = contigs[bi];
                            if ovl < cfg.min_ovl || ovl > b.len() {
                                continue;
                            }
                            let budget = ovl / cfg.mism_div;
                            let mut mism = 0usize;
                            for t in 0..ovl {
                                if a[s + t] != b[t] {
                                    mism += 1;
                                    if mism > budget {
                                        break;
                                    }
                                }
                            }
                            if mism > budget {
                                continue;
                            }
                            let key = (usize::MAX - ovl, mism, bi, s);
                            if key < best_key {
                                best_key = key;
                                best = Some((bi as u32, s as u32));
                            }
                        }
                    }
                }
                pos_a += 1;
            }
            best
        })
        .collect();

    // 3. Resolve successor edges into simple chains: each contig gets at most one
    //    successor and one predecessor, no cycles (union-find). Deterministic:
    //    accept edges in contig order.
    let mut parent: Vec<u32> = (0..nc as u32).collect();
    let mut pred_taken = vec![false; nc];
    let mut chosen: Vec<Option<(u32, u32)>> = vec![None; nc];
    for ai in 0..nc {
        if let Some((bi, s)) = succ[ai] {
            let b = bi as usize;
            if pred_taken[b] {
                continue;
            }
            if uf_find(&mut parent, ai as u32) == uf_find(&mut parent, bi) {
                continue; // would close a cycle
            }
            chosen[ai] = Some((bi, s));
            pred_taken[b] = true;
            let ra = uf_find(&mut parent, ai as u32);
            let rb = uf_find(&mut parent, bi);
            parent[ra as usize] = rb;
        }
    }

    // 4. Walk each chain head (no predecessor) into a super-contig, recording each
    //    original contig's (super id, offset). Overlap bytes come from the earlier
    //    contig; only a successor's non-overlapping tail is appended.
    let mut super_id = vec![u32::MAX; nc];
    let mut super_off = vec![0u32; nc];
    let mut new_seq: Vec<u8> = Vec::with_capacity(reference.total_bases());
    let mut new_offs: Vec<usize> = vec![0];
    let mut sid = 0u32;
    for head in 0..nc {
        if pred_taken[head] {
            continue;
        }
        let super_start = new_seq.len();
        new_seq.extend_from_slice(contigs[head]);
        super_id[head] = sid;
        super_off[head] = 0;
        let mut cur = head;
        let mut base = 0usize; // super-offset of `cur` relative to super_start
        while let Some((bi, s)) = chosen[cur] {
            let bi = bi as usize;
            let bbase = base + s as usize;
            super_id[bi] = sid;
            super_off[bi] = bbase as u32;
            let cur_super_len = new_seq.len() - super_start;
            let b = contigs[bi];
            if bbase + b.len() > cur_super_len {
                let new_from = cur_super_len - bbase; // first novel base of B
                new_seq.extend_from_slice(&b[new_from..]);
            }
            base = bbase;
            cur = bi;
        }
        new_offs.push(new_seq.len());
        sid += 1;
    }

    // Remap each read onto its merged super-contig.
    let new_places: Vec<Place4> = places
        .iter()
        .map(|p| {
            let oc = p.ci as usize;
            Place4 {
                ci: super_id[oc],
                off: super_off[oc] + p.off,
            }
        })
        .collect();

    // Re-vote the merged consensus: fold every read into its remapped position and
    // take the per-column plurality (ties to the lowest base, matching the greedy
    // assembler). Overlap columns now reflect all reads, so the reads mismatch the
    // reference less — recovering most of the block-byte cost the layout-only merge
    // would otherwise add. Columns with no ACGT vote keep their laid-down byte
    // (preserving non-ACGT reference content).
    let mut votes = vec![[0u32; 4]; new_seq.len()];
    for (r, pl) in reads.iter().zip(&new_places) {
        let start = new_offs[pl.ci as usize] + pl.off as usize;
        for (j, &byte) in r.iter().enumerate() {
            let c = code(byte);
            if c < 4 {
                votes[start + j][c as usize] += 1;
            }
        }
    }
    // Per-column plurality is independent per position, so resolve in parallel.
    // Deterministic: each output byte is a pure function of that column's votes
    // (ties to the lowest base via `Reverse(k)`).
    new_seq
        .par_iter_mut()
        .zip(votes.par_iter())
        .for_each(|(byte, v)| {
            if v.iter().any(|&x| x > 0) {
                let best = (0..4)
                    .max_by_key(|&k| (v[k], std::cmp::Reverse(k)))
                    .unwrap();
                *byte = b"ACGT"[best];
            }
        });

    let merged = GlobalReference {
        seq: new_seq,
        offs: new_offs,
    };
    (merged, new_places)
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
        fqxv_bytes::read_varint(self.buf, &mut self.pos).ok_or(Error::Malformed("varint too long"))
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

    /// Encode a whole read set with the two-pass v4 codec as a single block and
    /// decode it against the frozen reference; must be byte-exact.
    fn v4_roundtrip_one_block(reads: &[&[u8]], anchors: &[u32]) -> Vec<Vec<u8>> {
        let (reference, places) = assemble_global(reads, anchors);
        // Reference must serialize/deserialize losslessly too.
        let ref_bytes = reference.encode(6, 0, 0).expect("ref encode");
        let reference = GlobalReference::decode(&ref_bytes).expect("ref decode");
        let enc = encode_global_block(reads, &places, &reference).expect("v4 encode");
        // Determinism: re-encoding the same input is byte-identical.
        assert_eq!(
            encode_global_block(reads, &places, &reference).expect("v4 again"),
            enc,
            "v4 encode not deterministic"
        );
        decode_global_block(&enc, &reference).expect("v4 decode")
    }

    /// Assemble, overlap-merge the reference, then encode/decode every read
    /// against the MERGED reference — the placement remap must stay byte-exact.
    fn v4_merged_roundtrip(reads: &[&[u8]], anchors: &[u32]) -> (usize, usize, Vec<Vec<u8>>) {
        let (reference, places) = assemble_global(reads, anchors);
        let before = reference.n_contigs();
        let (merged, mplaces) = merge_reference(reads, &reference, &places);
        let after = merged.n_contigs();
        // Merged reference must serialize/reload and every read must round-trip.
        let ref_bytes = merged.encode(6, 0, 0).expect("ref encode");
        let merged = GlobalReference::decode(&ref_bytes).expect("ref decode");
        let enc = encode_global_block(reads, &mplaces, &merged).expect("v4 encode");
        let dec = decode_global_block(&enc, &merged).expect("v4 decode");
        (before, after, dec)
    }

    #[test]
    fn merge_reference_roundtrips_and_shrinks() {
        // Overlapping windows of one reference fragment into several contigs under
        // the greedy pass; overlap-merge should chain them into fewer contigs and
        // still round-trip byte-exactly.
        let reference_seq = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGCATCGATCGTAGCTAGCATTTACAGGTACCATGACATTGG";
        let win = 40usize;
        let (mut lens, mut seq) = (Vec::new(), Vec::new());
        for start in (0..=(reference_seq.len() - win)).step_by(3) {
            seq.extend_from_slice(&reference_seq[start..start + win]);
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
        let (before, after, dec) = v4_merged_roundtrip(&refs, &anchors);
        let expect: Vec<Vec<u8>> = refs.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect, "merged reference must round-trip");
        assert!(after <= before, "merge must not increase contig count");
    }

    #[test]
    fn global_roundtrip() {
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGTACGT",
            b"ACGTACGTACGT", // match
            b"ACGTAAGTACGT", // 1 mismatch
            b"ACGTACGTNCGT", // mismatch incl. N
            b"TTTTGGGGCCCC", // seeds a second contig
            b"",             // empty
            b"",             // match (empty == empty)
        ];
        let dec = v4_roundtrip_one_block(&reads, &vec![0u32; reads.len()]);
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn global_attaches_to_earlier_contig() {
        // The A,B,A,B,A interleave that strands the third A as a LITERAL under
        // v2: v4's global reference has one A contig every A read maps onto, so
        // there are exactly two contigs (one A, one B) — the reference dedups.
        let a = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGC";
        let b = b"TTAGGCCATTACAGGTACCATGACATTGGACATTACAGGTTCAAGT";
        let reads: Vec<&[u8]> = vec![&a[..], &b[..], &a[..], &b[..], &a[..]];
        let anchors = vec![0u32; reads.len()];
        let (reference, _places) = assemble_global(&reads, &anchors);
        assert_eq!(
            reference.n_contigs(),
            2,
            "reference should hold two contigs"
        );
        let dec = v4_roundtrip_one_block(&reads, &anchors);
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn global_multi_block_shares_reference() {
        // Assemble globally, then code in several small blocks against the one
        // frozen reference (the container's parallel-block shape). Reads at block
        // boundaries can't be block-local MATCH, so this exercises the `(ci, off)`
        // fallback placement for every read.
        let reference_seq = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGCATCGATCGTAGCTAGCAT";
        let win = 30usize;
        let (mut lens, mut seq) = (Vec::new(), Vec::new());
        for start in 0..=(reference_seq.len() - win) {
            seq.extend_from_slice(&reference_seq[start..start + win]);
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

        let (reference, places) = assemble_global(&refs, &anchors);
        let ref_bytes = reference.encode(8, 0, 0).expect("ref encode");
        let reference = GlobalReference::decode(&ref_bytes).expect("ref decode");
        // Cut into 4-read blocks and code/decode each against the shared ref.
        let mut got: Vec<Vec<u8>> = Vec::new();
        let mut s = 0usize;
        while s < refs.len() {
            let e = (s + 4).min(refs.len());
            let enc = encode_global_block(&refs[s..e], &places[s..e], &reference).expect("enc");
            got.extend(decode_global_block(&enc, &reference).expect("dec"));
            s = e;
        }
        let expect: Vec<Vec<u8>> = refs.iter().map(|r| r.to_vec()).collect();
        assert_eq!(got, expect);
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

        // Harden the rescue codec against the way the container actually drives it:
        // longer reads, a realistic `seq_order`, and NON-uniform anchors (the
        // anchors above are all zero, so the anchor-implied offset path — the v2
        // fast candidate inside `Assembler::place` — went largely unexercised).
        // Also pins encode determinism, which the container's determinism invariant
        // relies on now that rescue is on by default.
        #[test]
        fn rescue_roundtrip_varied_anchors(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..80)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let anchors: Vec<u32> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| ((r.len() * 7 + i * 3) % 41) as u32)
                .collect();
            for order in [1usize, 8] {
                let v3 = encode_clustered_rescue(&refs, &anchors, order).expect("v3 encode");
                proptest::prop_assert_eq!(
                    decode_clustered_rescue(&v3).expect("v3 decode"),
                    reads.clone()
                );
                // Deterministic: re-encoding the same input yields identical bytes.
                proptest::prop_assert_eq!(
                    encode_clustered_rescue(&refs, &anchors, order).expect("v3 again"),
                    v3
                );
                // v2 must also round-trip under the same non-uniform anchors.
                let v2 = encode_clustered(&refs, &anchors, order).expect("v2 encode");
                proptest::prop_assert_eq!(
                    decode_clustered(&v2).expect("v2 decode"),
                    reads.clone()
                );
            }
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

        // v4 two-pass: assemble globally, serialize+reload the reference, then
        // code the whole set as one block and as several small blocks (the
        // parallel-decode shape) — both against the shared frozen reference. Uses
        // non-uniform anchors so the anchor-implied candidate path is exercised.
        #[test]
        fn global_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..80)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let anchors: Vec<u32> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| ((r.len() * 7 + i * 3) % 41) as u32)
                .collect();
            let (reference, places) = assemble_global(&refs, &anchors);
            let ref_bytes = reference.encode(4, 0, 0).expect("ref encode");
            proptest::prop_assert_eq!(
                reference.encode(4, 0, 0).expect("ref again"), ref_bytes.clone(),
                "reference encode not deterministic"
            );
            let reference = GlobalReference::decode(&ref_bytes).expect("ref decode");
            // Single block.
            let enc = encode_global_block(&refs, &places, &reference).expect("v4 encode");
            proptest::prop_assert_eq!(
                decode_global_block(&enc, &reference).expect("v4 decode"),
                reads.clone()
            );
            // Several small blocks sharing the reference.
            let mut got: Vec<Vec<u8>> = Vec::new();
            let mut s = 0usize;
            while s < refs.len() {
                let e = (s + 7).min(refs.len());
                let b = encode_global_block(&refs[s..e], &places[s..e], &reference).expect("enc");
                got.extend(decode_global_block(&b, &reference).expect("dec"));
                s = e;
            }
            proptest::prop_assert_eq!(got, reads);
        }

        // The overlap-merge refinement must preserve the v4 round-trip: assemble,
        // merge the reference (remapping placements), then encode/decode every read
        // against the merged reference. Also pins merge determinism and that it
        // never grows the contig set.
        #[test]
        fn merge_reference_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..80)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let anchors: Vec<u32> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| ((r.len() * 5 + i * 7) % 37) as u32)
                .collect();
            let (reference, places) = assemble_global(&refs, &anchors);
            let (merged, mplaces) = merge_reference(&refs, &reference, &places);
            proptest::prop_assert!(merged.n_contigs() <= reference.n_contigs());
            // Merge is deterministic.
            let (merged2, _) = merge_reference(&refs, &reference, &places);
            proptest::prop_assert_eq!(merged.total_bases(), merged2.total_bases());
            let ref_bytes = merged.encode(4, 0, 0).expect("ref encode");
            let merged = GlobalReference::decode(&ref_bytes).expect("ref decode");
            let enc = encode_global_block(&refs, &mplaces, &merged).expect("v4 encode");
            proptest::prop_assert_eq!(
                decode_global_block(&enc, &merged).expect("v4 decode"),
                reads.clone()
            );
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
