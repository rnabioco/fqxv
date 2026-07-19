//! Rescue assembler: gapped placement of reads the clustered codec stranded.

use super::*;

/// k-mer length for the rescue index (matches the clustering minimizer k).
pub(crate) const RESCUE_K: usize = DEFAULT_K;

/// Forward 2-bit k-mer packed from `seq[start..start+k]`, or `None` if the
/// window runs off the end or contains a non-ACGT byte. `k <= 32`.
#[inline]
pub(crate) fn kmer_at(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }
    let mut v = 0u64;
    for &b in &seq[start..start + k] {
        let c = code_fold(b);
        if c >= 4 {
            return None;
        }
        v = (v << 2) | u64::from(c);
    }
    Some(v)
}

/// A chosen placement of a read on a contig.
pub(crate) struct Placement {
    pub(crate) ci: usize,
    pub(crate) off: usize,
    pub(crate) overlap: usize,
    pub(crate) mism: Vec<usize>,
}

/// Multi-contig assembler with an encoder-side k-mer index. Shared by the
/// rescue encoder and its op-mix diagnostic so both make identical decisions.
#[derive(Default)]
pub(crate) struct Assembler {
    pub(crate) contigs: Vec<Vec<Column>>,
    pub(crate) ref_anchors: Vec<u32>,
    /// k-mer -> (contig index, column position); most-recent occurrence wins.
    pub(crate) index: IntMap<u64, (u32, u32)>,
}

/// Acceptance test: can `cur` sit on `contig` at `off`? Returns
/// `(overlap, mismatch_positions)` when it is cheaper than a literal.
pub(crate) fn try_place(contig: &[Column], cur: &[u8], off: usize) -> Option<(usize, Vec<usize>)> {
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

/// Cheapest gapped answer to "would one small indel have saved this read?".
///
/// A banded Levenshtein between `cur` and the contig's consensus at `off`,
/// returning `(edits, indels)` for the best alignment, or `None` if the band
/// cannot reach the end. `band` is tiny (a few bases): the question is whether a
/// *small* indel explains an ungapped rejection, not whether some arbitrary
/// gapped alignment exists.
///
/// Deliberately its own routine rather than `fqxv-lroverlap`'s `align_banded`:
/// that crate is a sibling in the DAG, and a diagnostic is not a reason to add an
/// edge. If the answer here says indels are worth coding, the real fix should
/// reuse that aligner properly — see issue #102.
pub(crate) fn gapped_edits(
    contig: &[Column],
    cur: &[u8],
    off: usize,
    band: usize,
) -> Option<(usize, usize)> {
    if off > contig.len() {
        return None;
    }
    let refr: Vec<u8> = contig[off..(off + cur.len() + band).min(contig.len())]
        .iter()
        .map(|c| c.base)
        .collect();
    let (n, m) = (refr.len(), cur.len());
    if n == 0 || m == 0 {
        return None;
    }
    // (edits, indels) per cell, minimising edits then indels. Full table: reads
    // are ~150 bp here, so this is a few thousand cells and clarity wins.
    //
    // The reference's 3' end is FREE: the read may stop anywhere in it, so the
    // answer is the best "read fully consumed" cell across ALL rows, not the last
    // one. A global alignment would force the read to consume the whole
    // reference — `band` extra deletions it never made — which both inflates the
    // edit count and makes `indels > 0` vacuously true, so every literal would
    // look indel-rescuable. The read's 5' end stays anchored at `off`, with
    // ordinary deletion cost for shifting, because that is what `try_place`
    // itself assumes.
    let inf = (usize::MAX / 4, usize::MAX / 4);
    let mut prev: Vec<(usize, usize)> = (0..=m).map(|j| (j, j)).collect();
    let mut cur_row = vec![inf; m + 1];
    let mut best = prev[m];
    for i in 1..=n {
        cur_row[0] = (i, i);
        for j in 1..=m {
            let cost = usize::from(refr[i - 1] != cur[j - 1]);
            let diag = (prev[j - 1].0 + cost, prev[j - 1].1);
            let up = (prev[j].0 + 1, prev[j].1 + 1); // consume ref = deletion
            let left = (cur_row[j - 1].0 + 1, cur_row[j - 1].1 + 1); // insertion
            cur_row[j] = diag.min(up).min(left);
        }
        best = best.min(cur_row[m]);
        std::mem::swap(&mut prev, &mut cur_row);
        cur_row.fill(inf);
    }
    (best.0 < usize::MAX / 4).then_some(best)
}

/// Issue #102: does `fqxv-reorder` strand Illumina reads as literals *because it
/// has no indel op*?
///
/// `try_place` is an ungapped compare, so one indel shifts every base after it
/// and the read mismatches its way past the `overlap / 4` budget — stranding a
/// read that is otherwise a perfect neighbour. This counts how often that
/// actually happens, by re-testing every literal against the SAME candidates the
/// assembler considered, with a small band.
///
/// A read is counted `indel_rescuable` only if a gapped compare lands inside the
/// same budget the ungapped one blew, using at least one indel. That is the
/// narrowest reading of the question: not "could some aligner place this" but
/// "was the *lack of an indel op* the reason it was refused".
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct IndelProbe {
    /// Reads examined (excludes exact-duplicate MATCH reads).
    pub reads: u64,
    /// Reads the assembler stranded as literals.
    pub literals: u64,
    /// Literals a gapped compare would have placed within the same mismatch
    /// budget, using at least one indel.
    pub indel_rescuable: u64,
    /// Literals that stay unplaceable even gapped — genuinely novel sequence.
    pub truly_novel: u64,
    /// Bases held in `indel_rescuable` literals: what an indel op could reclaim.
    pub rescuable_bases: u64,
    /// Total literal bases, for scale.
    pub literal_bases: u64,
}

/// Run the [`IndelProbe`] over one block of clustered reads.
#[must_use]
pub fn indel_probe(reads: &[&[u8]], anchors: &[u32], band: usize) -> IndelProbe {
    let mut st = IndelProbe::default();
    let mut asm = Assembler::default();
    for (i, &cur) in reads.iter().enumerate() {
        if i > 0 && cur == reads[i - 1] {
            continue; // MATCH: never reaches the placer
        }
        st.reads += 1;
        match asm.place(cur, anchors[i]) {
            Some(p) => asm.commit(p.ci, cur, p.off, p.overlap),
            None => {
                st.literals += 1;
                st.literal_bases += cur.len() as u64;
                // Re-test the assembler's own candidates, gapped.
                let budget = cur.len() / 4;
                let rescued = asm
                    .candidates(cur, anchors[i])
                    .into_iter()
                    .filter_map(|(ci, off)| gapped_edits(&asm.contigs[ci], cur, off, band))
                    .any(|(edits, indels)| edits <= budget && indels > 0);
                if rescued {
                    st.indel_rescuable += 1;
                    st.rescuable_bases += cur.len() as u64;
                } else {
                    st.truly_novel += 1;
                }
                asm.seed(cur, anchors[i]);
            }
        }
    }
    st
}

impl Assembler {
    /// Index every k-mer starting in `[from, to)` of contig `ci`'s consensus.
    /// Only called on freshly-appended columns, so cost is linear in new bases;
    /// overlap columns whose consensus later shifts are left stale on purpose —
    /// the index only proposes candidates, [`try_place`] validates against the
    /// live consensus, so staleness costs recall, never correctness.
    fn index_range(&mut self, ci: usize, from: usize, to: usize) {
        // Low 2*RESCUE_K bits: the rolling window of the last RESCUE_K bases.
        const MASK: u64 = (1u64 << (2 * RESCUE_K)) - 1;
        let Self { contigs, index, .. } = self;
        let contig = &contigs[ci];
        let n = contig.len();
        let hi = to.min(n.saturating_sub(RESCUE_K - 1)); // k-mer starts are `< hi`
        if from >= hi {
            return;
        }
        // Roll a 2-bit-packed k-mer across the columns instead of recomputing all
        // RESCUE_K bases at every start (O(bases) not O(bases * k)). A non-ACGT
        // base resets the run, so only all-ACGT windows are indexed — identical
        // keys and insert order to the per-start recompute, so output is unchanged.
        let last = hi + RESCUE_K - 1; // exclusive column bound, `<= n`
        let mut v = 0u64;
        let mut run = 0usize;
        for p in from..last {
            let c = code_fold(contig[p].base);
            if c >= 4 {
                run = 0;
                continue;
            }
            v = ((v << 2) | u64::from(c)) & MASK;
            run += 1;
            if run >= RESCUE_K {
                index.insert(v, (ci as u32, (p + 1 - RESCUE_K) as u32));
            }
        }
    }

    /// Best placement of `cur` (minimizer at `anchor`) across all contigs, or
    /// `None` if it should seed a new one. Candidates come from the most-recent
    /// contig at the anchor-implied offset (the v2 fast path) plus every contig a
    /// sampled read k-mer points at. Deterministic: candidates are deduped and
    /// scored by (mismatches, recency, offset), independent of hash iteration.
    /// The (contig, offset) pairs worth validating for `cur`, from its anchor
    /// against the newest contig plus every k-mer index hit. Split out of
    /// [`Assembler::place`] so a diagnostic can re-test the SAME candidates the
    /// codec saw — a probe that generated its own candidates would be measuring
    /// its own recall, not the codec's.
    fn candidates(&self, cur: &[u8], anchor: u32) -> Vec<(usize, usize)> {
        let mut cands: Vec<(usize, usize)> = Vec::new();
        if cur.is_empty() || self.contigs.is_empty() {
            return cands;
        }
        let last = self.contigs.len() - 1;
        let center = self.ref_anchors[last] as i64 - anchor as i64;
        if center >= 0 && center as usize <= self.contigs[last].len() {
            cands.push((last, center as usize));
        }
        // Non-overlapping k-mers cover every base, so a read with a few errors
        // still has clean k-mers to match on.
        let mut start = 0;
        while start + RESCUE_K <= cur.len() {
            if let Some(code) = kmer_at(cur, start, RESCUE_K)
                && let Some(&(ci, cpos)) = self.index.get(&code)
            {
                let off = cpos as i64 - start as i64;
                if off >= 0 && off as usize <= self.contigs[ci as usize].len() {
                    cands.push((ci as usize, off as usize));
                }
            }
            start += RESCUE_K;
        }
        cands.sort_unstable();
        cands.dedup();
        cands
    }

    pub(crate) fn place(&self, cur: &[u8], anchor: u32) -> Option<Placement> {
        let mut best: Option<Placement> = None;
        let mut best_key = (usize::MAX, usize::MAX, usize::MAX);
        for (ci, off) in self.candidates(cur, anchor) {
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
    pub(crate) fn commit(&mut self, ci: usize, cur: &[u8], off: usize, overlap: usize) {
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
    pub(crate) fn seed(&mut self, cur: &[u8], anchor: u32) {
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
    let n = check_n(r.varint()? as usize)?;
    let s_ops = r.take_stream()?;
    let s_cref = r.take_stream()?;
    let s_offdelta = r.take_stream()?;
    let s_slen = r.take_stream()?;
    let s_nmis = r.take_stream()?;
    let s_pos = r.take_stream()?;
    let s_subs = r.take_stream()?;
    let s_novel = r.take_stream()?;
    let s_lit = r.take_stream()?;

    let ops = fqxv_rans::decode_bounded(s_ops, n)?;
    let cref = fqxv_rans::decode_bounded(s_cref, per_read_varints(n))?;
    let offdelta = fqxv_rans::decode_bounded(s_offdelta, per_read_varints(n))?;
    let slen = fqxv_rans::decode_bounded(s_slen, per_read_varints(n))?;
    let nmis = fqxv_rans::decode_bounded(s_nmis, per_read_varints(n))?;
    let subs = fqxv_rans::decode_bounded(s_subs, MAX_DECODED_BASES)?;
    let pos = fqxv_rans::decode_bounded(s_pos, subs.len().saturating_mul(10))?;
    let (_, novel) = fqxv_seq::decode(s_novel)?;
    let (lit_lens, lit_seq) = fqxv_seq::decode(s_lit)?;

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
                let mut read = alloc_read(cur_len)?;
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
/// [`encode_clustered_rescue`], driving the same `Assembler` so the counts
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
