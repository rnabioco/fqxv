//! Global-reference assembly: one frozen consensus reference per block.

use super::*;
use rayon::prelude::*;

/// The frozen global reference produced by [`assemble_global`]: the final
/// plurality-consensus bytes of every contig, concatenated, with per-contig
/// offsets. Reads are coded as positions on it and decoded by slicing it.
#[derive(Debug, Default, Clone)]
pub struct GlobalReference {
    /// Concatenated final consensus bytes of all contigs.
    pub(crate) seq: Vec<u8>,
    /// Byte offset of each contig in `seq`; `offs.len() == n_contigs + 1`.
    pub(crate) offs: Vec<usize>,
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
    pub(crate) fn contig(&self, ci: usize) -> &[u8] {
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
        // `nb` is an untrusted block count. Reserving it FALLIBLY is not enough:
        // `try_reserve` still asks the allocator for the whole thing, and under
        // memory overcommit that request is granted and the process is OOM-killed
        // later while the buffer is touched (the same reason `decode_bounded`
        // exists — see its docs). The `reorder` fuzz target found this: a 5-byte
        // input whose leading varint reads 32905357309 produced a single
        // malloc(789728575416), 24 bytes per `(usize, &[u8])`.
        //
        // Unlike an rANS stream, this one has an exact structural bound. Every
        // block costs at least two varints here (`n_contigs` and `len`), so a
        // stream of `src.len()` bytes cannot describe more than half that many
        // blocks. Check it before the reserve, and the count can never multiply
        // into an allocation the input did not pay for.
        if nb > src.len() / 2 {
            return Err(Error::Malformed(
                "reference block count exceeds stream size",
            ));
        }
        let mut blocks: Vec<(usize, &[u8])> = Vec::new();
        blocks
            .try_reserve_exact(nb)
            .map_err(|_| Error::Malformed("reference block count too large to allocate"))?;
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

    /// LZMA-class variant: explicit LZ77 matching (hash-chain finder) with LZMA
    /// entropy coding — context literals with matched-byte prediction, a length
    /// coder, position-slot + aligned distance coding, and rep0–3 short codes —
    /// over the [`fqxv_range`](fqxv_range) range coder. This *copies* the
    /// long-range near-duplicate contigs the order-k model and BWT can only model,
    /// the redundancy that separates xz (~1.79 b/base here) from the context model
    /// (~1.98). Coded over one whole-reference window (a serial decode, fine since
    /// the reference is coded once per file); deterministic, so thread-count
    /// independent. See the `fqxv_seq::lzma` module. Gated never-worse by the caller.
    pub fn encode_lzma(&self) -> Result<Vec<u8>> {
        Ok(fqxv_seq::lzma::encode(&self.contig_lens(), &self.seq)?)
    }

    /// Reverse of [`encode_lzma`](GlobalReference::encode_lzma).
    pub fn decode_lzma(src: &[u8]) -> Result<GlobalReference> {
        let (lens, seq) = fqxv_seq::lzma::decode(src)?;
        Self::from_lens_seq(&lens, seq)
    }

    /// SPRING-faithful reference coder: **2-bit-pack the ACGT consensus (4
    /// bases/byte), then LZMA the packed bytes** — exactly SPRING's
    /// `pack_compress_seq` (2-bit pack + BSC). The packing is a hard 2 bits/base
    /// floor and a byte-domain LZ then captures the long-range near-duplicate
    /// repeats; on real references this beats the order-k model on raw bases by
    /// ~7%. Non-ACGT bytes (rare in a plurality consensus) are exception-coded so
    /// it stays byte-exact. See the `refpack` module. Gated never-worse.
    pub fn encode_packed(&self) -> Result<Vec<u8>> {
        refpack::encode(&self.contig_lens(), &self.seq)
    }

    /// Reverse of [`encode_packed`](GlobalReference::encode_packed).
    pub fn decode_packed(src: &[u8]) -> Result<GlobalReference> {
        let (lens, seq) = refpack::decode(src)?;
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
/// contigs (the multi-contig `Assembler`, never reset), freeze the final
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
pub(crate) fn assemble_window(reads: &[&[u8]], anchors: &[u32]) -> (GlobalReference, Vec<Place4>) {
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
/// `assemble_window`, then concatenate their frozen references (remapping each
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
    let n = check_n(r.varint()? as usize)?;
    let s_ops = r.take_stream()?;
    let s_cid = r.take_stream()?;
    let s_offdelta = r.take_stream()?;
    let s_slen = r.take_stream()?;
    let s_nmis = r.take_stream()?;
    let s_pos = r.take_stream()?;
    let s_subs = r.take_stream()?;
    let s_tail = r.take_stream()?;

    let ops = fqxv_rans::decode_bounded(s_ops, n)?;
    let cid = fqxv_rans::decode_bounded(s_cid, per_read_varints(n))?;
    let offdelta = fqxv_rans::decode_bounded(s_offdelta, per_read_varints(n))?;
    let slen = fqxv_rans::decode_bounded(s_slen, per_read_varints(n))?;
    let nmis = fqxv_rans::decode_bounded(s_nmis, per_read_varints(n))?;
    let subs = fqxv_rans::decode_bounded(s_subs, MAX_DECODED_BASES)?;
    let pos = fqxv_rans::decode_bounded(s_pos, subs.len().saturating_mul(10))?;
    // Novel tail bases: one byte each, so the block's sequence budget bounds it.
    let tail = fqxv_rans::decode_bounded(s_tail, MAX_DECODED_BASES)?;

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
                let mut read = alloc_read(cur_len)?;
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
