//! Clustered contig-assembly sequence codec (the `encode_clustered` path).

use super::*;

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
        // window). See `place_on_contig`.
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
/// (via the shared `place_on_contig`) but skips entropy coding, so the counts
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

/// Allocate a `len`-byte read buffer fallibly. `len` is a per-read length decoded
/// from an untrusted stream, so a corrupt value must error rather than abort the
/// process on a huge infallible allocation.
pub(crate) fn alloc_read(len: usize) -> Result<Vec<u8>> {
    // Reorder is a short-read layout (mean length <= `REORDER_MAX_MEAN_LEN`), so a
    // single read this large is a corrupt length, not real data. Reject it before
    // allocating — a bomb declaring a multi-GB read would otherwise reserve (and
    // fill) that much before any downstream check.
    const MAX_READ_LEN: usize = 1 << 24; // 16 MiB — far above any real read
    if len > MAX_READ_LEN {
        return Err(Error::Malformed("read length implausibly large"));
    }
    let mut read = Vec::new();
    read.try_reserve_exact(len)
        .map_err(|_| Error::Malformed("read length too large to allocate"))?;
    read.resize(len, 0);
    Ok(read)
}

/// Decode a stream produced by [`encode_clustered`], returning the reads in the
/// same (clustered) order.
pub fn decode_clustered(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != 2 {
        return Err(Error::Malformed("unsupported version"));
    }
    let n = check_n(r.varint()? as usize)?;
    // Borrow every stream in the order the encoder wrote them, then decode in an
    // order that lets each bound the next. `take_stream` borrows from `src`, not
    // from the cursor, so the two orders need not agree.
    let s_ops = r.take_stream()?;
    let s_offdelta = r.take_stream()?;
    let s_slen = r.take_stream()?;
    let s_nmis = r.take_stream()?;
    let s_pos = r.take_stream()?;
    let s_subs = r.take_stream()?;
    let s_novel = r.take_stream()?;
    let s_lit = r.take_stream()?;

    let ops = fqxv_rans::decode_bounded(s_ops, n)?; // exactly one op per read
    let offdelta = fqxv_rans::decode_bounded(s_offdelta, per_read_varints(n))?;
    let slen = fqxv_rans::decode_bounded(s_slen, per_read_varints(n))?;
    let nmis = fqxv_rans::decode_bounded(s_nmis, per_read_varints(n))?;
    // `subs` is one byte per mismatch, so decoding it first SIZES `pos`, which is
    // one varint per mismatch. That is a real bound rather than a guess.
    let subs = fqxv_rans::decode_bounded(s_subs, MAX_DECODED_BASES)?;
    let pos = fqxv_rans::decode_bounded(s_pos, subs.len().saturating_mul(10))?;
    let (_, novel) = fqxv_seq::decode(s_novel)?;
    let (lit_lens, lit_seq) = fqxv_seq::decode(s_lit)?;

    let mut c_offdelta = Cursor::new(&offdelta);
    let mut c_slen = Cursor::new(&slen);
    let mut c_nmis = Cursor::new(&nmis);
    let mut c_pos = Cursor::new(&pos);
    let (mut subs_pos, mut lit_pos, mut lit_idx, mut novel_pos) = (0usize, 0usize, 0usize, 0usize);
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n.min(1 << 22));

    // The current contig, voted identically to the encoder.
    let mut contig: Vec<Column> = Vec::new();
    let mut prev_off: usize = 0;
    // Running total of reconstructed bases, bounded below (see the accumulation).
    let mut out_bases = 0usize;

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
                let mut read = alloc_read(cur_len)?;
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
        // MATCH clones and CONTIG copies can expand a few KB of coded streams into
        // unbounded output. Each per-read length is already capped (`alloc_read`,
        // the literal streams), but the aggregate is not — bound it against the
        // same per-block ceiling the rANS streams use. Reorder is a short-read
        // layout written at <= `REORDER_BLOCK_READS` reads/block (~128 MiB honest
        // max, half this ceiling), so only an amplification bomb trips it.
        out_bases = out_bases
            .checked_add(reads.last().map_or(0, Vec::len))
            .filter(|&t| t <= MAX_DECODED_BASES)
            .ok_or(Error::Malformed(
                "reconstructed output exceeds decode limit",
            ))?;
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
