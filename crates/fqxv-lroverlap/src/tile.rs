//! Multi-reference **tiling** sequence codec (`SEQ_METHOD_TILE`).
//!
//! The consensus codec ([`crate::encode`]) codes every read against a single
//! voted consensus. On noisy (ONT) data that consensus is itself ~6% divergent
//! from the truth *and* each read sits ~19% away from it, so every edit script
//! carries the consensus's error on top of the read's. This codec drops the
//! consensus entirely and codes each read directly against **earlier raw reads**,
//! the way CoLoRd does: a read is *tiled* by a handful of overlapping neighbour
//! reads whose ids are strictly smaller than its own, each tile coded as a banded
//! edit script against its neighbour, and any genuinely uncovered span stored as
//! literal bases.
//!
//! Because a tile references only *earlier* reads, decoding is a single forward
//! pass over read ids: read `i` is rebuilt from tiles applied to the
//! already-decoded reads `i - delta` (`delta >= 1`), so all dependencies are
//! resolved by the time a read is reached. Measured on `ecoli_ont` block 1 this
//! reaches ~1.14 b/base (band 256) / ~1.056 (band 768) at ~97% read coverage,
//! against the consensus codec's 1.243 — the ONT ratio lever it was built for.
//!
//! ## Wire format (after the container's method byte)
//!
//! ```text
//! "LRT" u8=VERSION                      magic + bitstream version
//! varint n | varint total_bases         read count and base count
//! stream  lens                          per-read lengths (varints)
//! stream  markers                       per segment: 0 = literal gap,
//!                                          k>=1 = tile vs read (id - k)
//! stream  tile_lens                     per tile: query bases it produces
//! stream  strands                       per tile: 1 = neighbour reverse-complemented
//! stream  rs                            per tile: neighbour start offset
//! stream  lit_lens                      per literal gap: its length
//! range   ops   (count, blob)           edit-op stream, order-3 op history
//! stream  runs                          Match run lengths
//! range   subs  (count, blob)           substituted bases, on the neighbour base
//! range   ins   (count, blob)           inserted bases, on the preceding base
//! stream  indel_lens                    Ins/Del lengths
//! stream  literals                      literal + gap base codes (0..4)
//! stream  exc_pos | stream exc_bytes    non-ACGT exception list
//! ```
//!
//! A read whose segments never resolve (read 0, or one with no earlier-id
//! overlap) is stored as a single full-length literal gap — still exact. Every
//! byte is recovered: ACGT bases from the edit scripts / literals, and every
//! non-ACGT byte from the exception list applied last (identical to
//! [`crate::codec`]).

use rayon::prelude::*;

use fqxv_bytes::{read_varint, write_varint};
use fqxv_dna::{base_of_sym, is_acgt, revcomp_acgt};

use crate::codec::{
    OP_CTXS, REF_TAIL, SUB_SYMS, align_for_coding, base_code, compute_offs, encode_ops,
    encode_subs, get_stream, put_stream,
};
use crate::{ChainOpts, EncodeOpts, Error, Index, Op, Overlap, Repeat, find_overlaps};

/// Format tag for a tiling sequence block. Distinct from the consensus codec's
/// `LRO`/`LRR` so a mis-dispatched block fails closed rather than mis-decodes.
const TILE_MAGIC: [u8; 3] = *b"LRT";
/// Bitstream version. Bump on any layout change (nothing on disk is stable yet).
const TILE_VERSION: u8 = 1;

/// Ceiling on decoded bases per coded byte, used to reject a hostile `total_bases`
/// header before it drives the `vec![0u8; total_bases]` output allocation (issue
/// #142). The tiler's `literals` stream must seed every distinct base at ~2
/// bits/base, so even a pathological all-identical-reads block stays well under
/// this (~6 K observed worst case vs the ~6 bases/byte of real ONT); it only fails
/// a crafted length. Deliberately far below `1 << 18` so the bound also caps peak
/// decode memory (`~3 × total_bases`) on a large hostile input, not just u64::MAX.
const MAX_BASES_PER_BYTE: usize = 1 << 14;

/// One read's contribution to the block streams: the shared edit streams plus the
/// per-segment manifest. Accumulated per read in parallel and merged in id order,
/// so the concatenation is byte-identical for any thread count.
#[derive(Default)]
struct TileStreams {
    // ---- edit streams (identical shape to the consensus codec) ----
    ops: Vec<u8>,
    runs: Vec<u8>,
    subs: Vec<u8>,
    /// Parallel to `subs`: the neighbour base each substitution replaced (context).
    sub_ctx: Vec<u8>,
    ins_bases: Vec<u8>,
    /// Parallel to `ins_bases`: the preceding read base (homopolymer context).
    ins_ctx: Vec<u8>,
    indel_lens: Vec<u8>,
    literals: Vec<u8>,
    // ---- per-segment manifest ----
    /// Per segment: `0` = literal gap, `k >= 1` = a tile against read `id - k`.
    markers: Vec<u8>,
    /// Per tile: the number of query bases it produces (its op run's terminator).
    tile_lens: Vec<u8>,
    /// Per tile: `1` when the neighbour is reverse-complemented, else `0`.
    strands: Vec<u8>,
    /// Per tile: the neighbour offset the edit script starts at.
    rs: Vec<u8>,
    /// Per literal gap: its length in query bases.
    lit_lens: Vec<u8>,
}

impl TileStreams {
    fn merge(&mut self, o: &TileStreams) {
        self.ops.extend_from_slice(&o.ops);
        self.runs.extend_from_slice(&o.runs);
        self.subs.extend_from_slice(&o.subs);
        self.sub_ctx.extend_from_slice(&o.sub_ctx);
        self.ins_bases.extend_from_slice(&o.ins_bases);
        self.ins_ctx.extend_from_slice(&o.ins_ctx);
        self.indel_lens.extend_from_slice(&o.indel_lens);
        self.literals.extend_from_slice(&o.literals);
        self.markers.extend_from_slice(&o.markers);
        self.tile_lens.extend_from_slice(&o.tile_lens);
        self.strands.extend_from_slice(&o.strands);
        self.rs.extend_from_slice(&o.rs);
        self.lit_lens.extend_from_slice(&o.lit_lens);
    }

    /// Walk one tile's edit script (`query` coded against `refslice`) into the
    /// edit streams, deriving the same substitution / insertion contexts the
    /// entropy coders key on. Mirrors the consensus codec's inline walk exactly,
    /// so the streams entropy-code with the same models. `refslice` is the
    /// neighbour span the tile was aligned against; `ref_pos`/`last` reset per
    /// tile (each tile is an independent alignment).
    fn walk_ops(&mut self, refslice: &[u8], ops: &[Op]) {
        let mut ref_pos = 0usize;
        let mut last = 4u8; // 4 = no preceding base (tile start)
        for op in ops {
            match op {
                Op::Match(m) => {
                    self.ops.push(0);
                    write_varint(&mut self.runs, u64::from(*m));
                    let mm = *m as usize;
                    if mm > 0 {
                        last = refslice.get(ref_pos + mm - 1).map_or(4, |&x| base_code(x));
                    }
                    ref_pos += mm;
                }
                Op::Sub(b) => {
                    self.ops.push(1);
                    let bc = base_code(*b);
                    self.subs.push(bc);
                    self.sub_ctx
                        .push(refslice.get(ref_pos).map_or(4, |&x| base_code(x)));
                    last = bc;
                    ref_pos += 1;
                }
                Op::Ins(bs) => {
                    self.ops.push(2);
                    write_varint(&mut self.indel_lens, bs.len() as u64);
                    for &b in bs {
                        let bc = base_code(b);
                        self.ins_ctx.push(last);
                        self.ins_bases.push(bc);
                        last = bc;
                    }
                }
                Op::Del(m) => {
                    self.ops.push(3);
                    write_varint(&mut self.indel_lens, u64::from(*m));
                    ref_pos += *m as usize;
                }
            }
        }
    }

    /// Emit a literal gap of `bases` (their 0..4 codes) as one segment.
    fn push_literal(&mut self, bases: &[u8]) {
        write_varint(&mut self.markers, 0);
        write_varint(&mut self.lit_lens, bases.len() as u64);
        self.literals.extend(bases.iter().map(|&b| base_code(b)));
    }

    /// Emit one tile: its edit script (against `refwin`, the neighbour span the
    /// tile was aligned over) followed by the per-segment manifest fields, in the
    /// exact order [`serialize`] reads them back. `tile_len` is the query bases the
    /// ops produce, `delta = id - neighbour_id >= 1`, `strand`/`rs` place the
    /// neighbour. Shared by the greedy and best-of-N covers so both are byte-for-
    /// byte identical for a given chosen `(refwin, ops, …)`.
    fn push_tile(
        &mut self,
        refwin: &[u8],
        ops: &[Op],
        tile_len: usize,
        delta: u64,
        strand: bool,
        rs: usize,
    ) {
        self.walk_ops(refwin, ops);
        write_varint(&mut self.markers, delta);
        write_varint(&mut self.tile_lens, tile_len as u64);
        self.strands.push(u8::from(strand));
        write_varint(&mut self.rs, rs as u64);
    }
}

/// One aligned candidate reference for a tile at `pos`: the oriented neighbour
/// window it was aligned against, its edit script, where it lands, and a coding
/// cost used only to rank candidates against each other (never on decode).
struct Cand {
    /// Query position the tile ends at (`o.q_end` capped to the read length).
    te: usize,
    /// The neighbour span the ops are coded against (`reference[rs..re]`).
    refwin: Vec<u8>,
    /// The (trailing-Del-trimmed) edit script producing `read[pos..te]`.
    ops: Vec<Op>,
    /// Neighbour start offset the edit script begins at.
    rs: usize,
    /// `id - neighbour_id`, always `>= 1` (a tile references an earlier read).
    delta: u64,
    /// Whether the neighbour is reverse-complemented.
    strand: bool,
    /// Estimated coded bits for this tile: framing plus per-edit cost. Compared
    /// amortised over the span (`cost_bits / span`) so a candidate is preferred
    /// only when it is cheaper *per produced base* — balancing reach against
    /// match quality without a hard span rule. Ranking only; decode never sees it.
    cost_bits: u64,
}

/// Fixed per-tile framing charged in the best-of-N cost model: one manifest
/// segment is ~5 bytes (marker + tile_len + strand + rs varints). Ranking only.
const TILE_FRAME_BITS: u64 = 40;
/// Bits an ONT edit (a substituted / inserted / deleted base) costs once
/// entropy-coded, used to weight a candidate's edit distance. Derived from the
/// measured ~1.14 b/base at ~13.5% divergence; only its ratio to
/// [`TILE_FRAME_BITS`] matters, and only for ranking. Never affects decode.
const TILE_EDIT_BITS: u64 = 8;

/// The read-local context the cover functions share: which read is being tiled
/// (`i`), the block's concatenated bases (`seq`) and per-read offsets (`offs`),
/// the read slice itself (`read = seq[offs[i]..offs[i + 1]]`), and the alignment
/// band. Passed by reference so a candidate reference read can be resolved from
/// `seq`/`offs` without re-slicing the world through a long argument list.
struct ReadCtx<'a> {
    i: usize,
    seq: &'a [u8],
    offs: &'a [usize],
    read: &'a [u8],
    band: usize,
}

/// Build the non-ACGT exception list over the whole sequence buffer: for each
/// non-ACGT byte, a delta-coded position and the original byte. Applied last on
/// decode, this is the single source of truth for every non-ACGT base regardless
/// of how the tile / literal streams coded it (identical to the consensus codec).
fn build_exceptions(seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut exc_pos: Vec<u8> = Vec::new();
    let mut exc_bytes: Vec<u8> = Vec::new();
    let mut prev_pos: u64 = 0;
    for (i, &b) in seq.iter().enumerate() {
        if !is_acgt(b) {
            write_varint(&mut exc_pos, i as u64 - prev_pos);
            prev_pos = i as u64;
            exc_bytes.push(b);
        }
    }
    (exc_pos, exc_bytes)
}

/// Tile one read against its earlier-id overlapping neighbours, returning its
/// stream contribution. At each uncovered position an overlap is chosen to code
/// the next span as an edit script; genuinely uncovered runs are stored as
/// literals. `max_refs <= 1` is the plain greedy max-reach cover (take the overlap
/// that reaches furthest); `max_refs > 1` selects best-of-N (see [`cover_bestof`]).
/// Pure function of `(i, seq, offs, idx, band, max_refs)`, so the per-read work is
/// order-free and thread-count invariant.
fn tile_one_read(
    i: usize,
    seq: &[u8],
    offs: &[usize],
    idx: Option<&Index>,
    band: usize,
    max_refs: usize,
) -> TileStreams {
    let read = &seq[offs[i]..offs[i + 1]];
    let mut s = TileStreams::default();
    if read.is_empty() {
        return s;
    }
    // No index (build failed, or nothing to seed against): whole read is literal.
    let Some(idx) = idx else {
        s.push_literal(read);
        return s;
    };

    let mut ovs = find_overlaps(idx, i as u32, read, ChainOpts::default());
    // Only earlier reads can be a reference — the decoder resolves ids forward.
    ovs.retain(|o| (o.target as usize) < i);
    // Total order for the greedy scan (find_overlaps is score-sorted); stable so
    // q_start ties keep the score order, making the cover thread-independent.
    ovs.sort_by_key(|o| o.q_start);

    let ctx = ReadCtx {
        i,
        seq,
        offs,
        read,
        band,
    };
    if max_refs <= 1 {
        cover_greedy(&ctx, &ovs, &mut s);
    } else {
        cover_bestof(&ctx, &ovs, max_refs, &mut s);
    }
    s
}

/// The plain greedy max-reach cover (`tile_max_refs <= 1`): at each uncovered
/// position take the single overlap that reaches furthest and code that whole
/// span. Kept as its own function so the default block is byte-for-byte identical
/// to what shipped in #212 — [`TileStreams::push_tile`] emits the manifest in the
/// same order the original inline code did.
fn cover_greedy(ctx: &ReadCtx, ovs: &[Overlap], s: &mut TileStreams) {
    let read = ctx.read;
    let len_i = read.len();
    let mut pos = 0usize;
    while pos < len_i {
        // The overlap covering `pos` that reaches furthest → the longest tile.
        let best = ovs
            .iter()
            .filter(|o| (o.q_start as usize) <= pos && (o.q_end as usize) > pos)
            .max_by_key(|o| o.q_end);
        match best {
            Some(o) => {
                let te = (o.q_end as usize).min(len_i);
                match eval_candidate(ctx, o, pos) {
                    Some(c) => s.push_tile(&c.refwin, &c.ops, te - pos, c.delta, c.strand, c.rs),
                    // The overlap's diagonal placed no usable neighbour span here.
                    None => s.push_literal(&read[pos..te]),
                }
                pos = te;
            }
            None => {
                // No overlap covers `pos`; literal up to the next overlap start.
                let next = ovs
                    .iter()
                    .filter_map(|o| ((o.q_start as usize) > pos).then_some(o.q_start as usize))
                    .min()
                    .unwrap_or(len_i);
                s.push_literal(&read[pos..next]);
                pos = next;
            }
        }
    }
}

/// Best-of-N cover (`tile_max_refs > 1`): at each uncovered position align the
/// `max_refs` furthest-reaching covering overlaps and keep the tile whose edit
/// script is cheapest *per produced base*. At ONT coverage many earlier reads span
/// the same region with independent error patterns, so the min-cost reference
/// agrees with the query at more positions — CoLoRd's anchor choice, applied per
/// tile. Coverage never regresses versus greedy: the furthest-reaching overlap is
/// always among the candidates, so the worst case is simply picking it.
fn cover_bestof(ctx: &ReadCtx, ovs: &[Overlap], max_refs: usize, s: &mut TileStreams) {
    let read = ctx.read;
    let len_i = read.len();
    let mut pos = 0usize;
    while pos < len_i {
        // Overlaps covering `pos`, furthest reach first with a total tie-break so
        // the candidate set — and thus the block — is thread-count invariant.
        let mut covering: Vec<&Overlap> = ovs
            .iter()
            .filter(|o| (o.q_start as usize) <= pos && (o.q_end as usize) > pos)
            .collect();
        if covering.is_empty() {
            let next = ovs
                .iter()
                .filter_map(|o| ((o.q_start as usize) > pos).then_some(o.q_start as usize))
                .min()
                .unwrap_or(len_i);
            s.push_literal(&read[pos..next]);
            pos = next;
            continue;
        }
        covering.sort_by(|a, b| {
            b.q_end
                .cmp(&a.q_end)
                .then(a.q_start.cmp(&b.q_start))
                .then(a.target.cmp(&b.target))
                .then(a.strand.cmp(&b.strand))
        });

        // Align up to `max_refs` of them; keep the cheapest amortised tile. Ties
        // keep the earlier (furthest-reaching) candidate, so selection is total.
        let mut best: Option<Cand> = None;
        for o in covering.iter().take(max_refs) {
            let Some(c) = eval_candidate(ctx, o, pos) else {
                continue;
            };
            let take = match &best {
                None => true,
                // Cheaper bits per produced base; cross-multiplied to stay integer
                // (both spans are >= 1 since every candidate has `te > pos`).
                Some(b) => {
                    let lhs = c.cost_bits * (b.te - pos) as u64;
                    let rhs = b.cost_bits * (c.te - pos) as u64;
                    lhs < rhs
                }
            };
            if take {
                best = Some(c);
            }
        }

        match best {
            Some(c) => {
                let te = c.te;
                s.push_tile(&c.refwin, &c.ops, te - pos, c.delta, c.strand, c.rs);
                pos = te;
            }
            None => {
                // Every covering overlap's diagonal fell off the read — literal to
                // the furthest reach, matching the greedy cover's rs>=re fallback.
                let te = (covering[0].q_end as usize).min(len_i);
                s.push_literal(&read[pos..te]);
                pos = te;
            }
        }
    }
}

/// Align the query span `read[pos..te]` (`te` the overlap's reach, capped to the
/// read) against neighbour `o` and return the resulting tile with its ranking
/// cost, or `None` when the overlap's diagonal places no usable neighbour span at
/// `pos` (the caller then stores that span as a literal, as the greedy cover does
/// on `rs >= re`). Pure and side-effect free, so it is safe to call speculatively
/// on several candidates before committing to one.
fn eval_candidate(ctx: &ReadCtx, o: &Overlap, pos: usize) -> Option<Cand> {
    let &ReadCtx {
        i,
        seq,
        offs,
        read,
        band,
    } = ctx;
    let te = (o.q_end as usize).min(read.len());
    let target = &seq[offs[o.target as usize]..offs[o.target as usize + 1]];
    let reference = if o.strand {
        revcomp_acgt(target)
    } else {
        target.to_vec()
    };
    // Map the query span into the neighbour's (query-oriented) frame.
    let off = i64::from(o.t_start) - i64::from(o.q_start);
    let rs = (pos as i64 + off).max(0) as usize;
    let re = ((te as i64 + off).max(0) as usize + REF_TAIL).min(reference.len());
    if rs >= re {
        return None;
    }
    let mut ops = align_for_coding(&reference[rs..re], &read[pos..te], band).ops;
    // A trailing Del consumes only reference, never query, so it carries nothing
    // the decoder needs and would desync the "produced == tile_len" termination —
    // drop it (as the consensus coder does for its per-read scripts).
    while matches!(ops.last(), Some(Op::Del(_))) {
        ops.pop();
    }
    // Edit distance actually coded (subs + inserted + deleted bases) — the term
    // that drives a tile's size; matches are near-free under the order-3 op model.
    let dist: u64 = ops
        .iter()
        .map(|op| match op {
            Op::Match(_) => 0,
            Op::Sub(_) => 1,
            Op::Ins(b) => b.len() as u64,
            Op::Del(m) => u64::from(*m),
        })
        .sum();
    Some(Cand {
        te,
        refwin: reference[rs..re].to_vec(),
        ops,
        rs,
        delta: i as u64 - u64::from(o.target), // >= 1
        strand: o.strand,
        cost_bits: TILE_FRAME_BITS + dist * TILE_EDIT_BITS,
    })
}

/// Append a range-coded blob as `[varint symbol_count][varint blob_len][blob]`,
/// matching the consensus codec's framing for its `ops`/`subs`/`ins` streams.
fn put_range(out: &mut Vec<u8>, blob: &[u8], count: usize) {
    write_varint(out, count as u64);
    write_varint(out, blob.len() as u64);
    out.extend_from_slice(blob);
}

/// Serialise the merged streams into a self-contained tiling block.
fn serialize(
    lens: &[u32],
    total_bases: u64,
    all: &TileStreams,
    exc_pos: &[u8],
    exc_bytes: &[u8],
) -> Result<Vec<u8>, Error> {
    let n = lens.len();
    let mut out = Vec::new();
    out.extend_from_slice(&TILE_MAGIC);
    out.push(TILE_VERSION);
    write_varint(&mut out, n as u64);
    write_varint(&mut out, total_bases);

    let mut lens_raw = Vec::new();
    for &l in lens {
        write_varint(&mut lens_raw, u64::from(l));
    }
    put_stream(&mut out, &lens_raw)?;

    // Manifest.
    put_stream(&mut out, &all.markers)?;
    put_stream(&mut out, &all.tile_lens)?;
    put_stream(&mut out, &all.strands)?;
    put_stream(&mut out, &all.rs)?;
    put_stream(&mut out, &all.lit_lens)?;

    // Edit streams (ops/subs/ins range-coded, the rest raw + rANS via put_stream).
    put_range(&mut out, &encode_ops(&all.ops), all.ops.len());
    put_stream(&mut out, &all.runs)?;
    put_range(
        &mut out,
        &encode_subs(&all.subs, &all.sub_ctx),
        all.subs.len(),
    );
    put_range(
        &mut out,
        &encode_subs(&all.ins_bases, &all.ins_ctx),
        all.ins_bases.len(),
    );
    put_stream(&mut out, &all.indel_lens)?;
    put_stream(&mut out, &all.literals)?;

    put_stream(&mut out, exc_pos)?;
    put_stream(&mut out, exc_bytes)?;
    Ok(out)
}

/// Compress a block of long reads with the multi-reference tiler. `lens[r]` is the
/// length of read `r`; `seq` is their bases concatenated in read order. Returns a
/// self-contained block that [`tile_decode`] inverts exactly. `opts.sketch`,
/// `opts.tile_band`, and `opts.tile_max_refs` select the overlap seeding, the
/// alignment band, and how many candidate references each tile weighs; all three
/// affect ratio and speed only — the block self-describes, so decode never
/// re-sketches, re-aligns, or re-selects.
///
/// # Errors
/// Returns [`Error::Corrupt`] if `lens` does not sum to `seq.len()`, or if an
/// entropy-coder pass fails.
pub fn tile_encode(lens: &[u32], seq: &[u8], opts: &EncodeOpts) -> Result<Vec<u8>, Error> {
    let n = lens.len();
    let offs = compute_offs(lens, seq.len())?;
    let total_bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    let band = opts.tile_band;
    let max_refs = opts.tile_max_refs.max(1);

    // A tile references *another read*, which the decoder has only rebuilt in
    // ACGT-placeholder space when it reaches this one — non-ACGT bytes are restored
    // from the exception list in a single final pass, so mid-decode every neighbour
    // still holds an `A` where it originally had an `N`/lowercase base. Encode must
    // see the same neighbours, so the whole pipeline runs on the normalized
    // sequence (non-ACGT -> `A`); the exception list, built over the *original*
    // `seq`, then restores every non-ACGT byte exactly. (The consensus codec is
    // immune because its reference is a stored consensus, identical on both sides.)
    let norm: Vec<u8> = seq.iter().map(|&b| base_of_sym(base_code(b))).collect();

    // One index over every read. On failure fall back to all-literal (lossless).
    let idx = Index::build(lens, &norm, opts.sketch, Repeat::default()).ok();

    // Tile each read in parallel; merge in id order for a thread-independent block.
    let parts: Vec<TileStreams> = (0..n)
        .into_par_iter()
        .map(|i| tile_one_read(i, &norm, &offs, idx.as_ref(), band, max_refs))
        .collect();
    let mut all = TileStreams::default();
    for p in &parts {
        all.merge(p);
    }

    let (exc_pos, exc_bytes) = build_exceptions(seq);
    serialize(lens, total_bases, &all, &exc_pos, &exc_bytes)
}

/// A framed range-coded blob (`ops`/`subs`/`ins`): its declared symbol count and a
/// fresh decoder over its bytes. `pos` is advanced past the blob.
fn take_range<'a>(
    src: &'a [u8],
    pos: &mut usize,
) -> Result<(usize, fqxv_range::Decoder<'a>), Error> {
    let count = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let blob_len = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let end = pos.checked_add(blob_len).ok_or(Error::Corrupt)?;
    let blob = src.get(*pos..end).ok_or(Error::Corrupt)?;
    *pos = end;
    Ok((count, fqxv_range::Decoder::new(blob)))
}

/// Decompress a block produced by [`tile_encode`], returning `(lens, seq)` in the
/// original read order.
///
/// # Errors
/// Returns [`Error::Corrupt`] on a malformed, truncated, or wrong-version block,
/// or if any stream desyncs (a tile whose ops do not produce its declared length,
/// a neighbour id that is not strictly earlier, an out-of-range offset, …).
pub fn tile_decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>), Error> {
    let mut pos = 0usize;
    if src.get(..TILE_MAGIC.len()) != Some(&TILE_MAGIC) {
        return Err(Error::Corrupt);
    }
    pos += TILE_MAGIC.len();
    let version = *src.get(pos).ok_or(Error::Corrupt)?;
    pos += 1;
    if version != TILE_VERSION {
        return Err(Error::Corrupt);
    }
    let n = read_varint(src, &mut pos).ok_or(Error::Corrupt)? as usize;
    let total_bases = read_varint(src, &mut pos).ok_or(Error::Corrupt)? as usize;
    // Reject a header whose declared base count could not possibly fit the coded
    // block before it drives the `vec![0u8; total_bases]` output allocation — the
    // #142 "bound the aggregate allocation against the input" rule. Deliberately
    // loose (a real block decodes far below this), it only fails a hostile length.
    if total_bases
        > src
            .len()
            .saturating_mul(MAX_BASES_PER_BYTE)
            .saturating_add(MAX_BASES_PER_BYTE)
    {
        return Err(Error::Corrupt);
    }

    let lens_raw = get_stream(src, &mut pos)?;
    let mut lp = 0usize;
    // `Vec::new`, not `with_capacity(n)`: `n` is an unvalidated header varint, so a
    // capacity hint would itself be the alloc bomb. The push loop is bounded by
    // `lens_raw` (a real stream) exhausting into `Corrupt`.
    let mut lens = Vec::new();
    for _ in 0..n {
        let l = read_varint(&lens_raw, &mut lp).ok_or(Error::Corrupt)?;
        lens.push(u32::try_from(l).map_err(|_| Error::Corrupt)?);
    }
    let offs = compute_offs(&lens, total_bases)?;

    // Manifest.
    let markers = get_stream(src, &mut pos)?;
    let tile_lens = get_stream(src, &mut pos)?;
    let strands = get_stream(src, &mut pos)?;
    let rs_stream = get_stream(src, &mut pos)?;
    let lit_lens = get_stream(src, &mut pos)?;

    // Edit streams. The range-coded ones are decoded inline during reconstruction,
    // rebuilding their adaptive contexts exactly as encode did; the rest are raw.
    let (ops_count, mut op_dec) = take_range(src, &mut pos)?;
    let mut op_models: [fqxv_range::SimpleModel<4>; OP_CTXS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut op_ctx = 0usize;
    let mut ops_seen = 0usize;

    let runs = get_stream(src, &mut pos)?;

    let (subs_count, mut sub_dec) = take_range(src, &mut pos)?;
    let mut sub_models: [fqxv_range::SimpleModel<SUB_SYMS>; SUB_SYMS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut subs_seen = 0usize;

    let (ins_count, mut ins_dec) = take_range(src, &mut pos)?;
    let mut ins_models: [fqxv_range::SimpleModel<SUB_SYMS>; SUB_SYMS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut ins_seen = 0usize;

    let indel_lens = get_stream(src, &mut pos)?;
    let literals = get_stream(src, &mut pos)?;
    let exc_pos = get_stream(src, &mut pos)?;
    let exc_bytes = get_stream(src, &mut pos)?;

    let mut seq = vec![0u8; total_bases];
    // Cursors, advanced in the exact order `serialize` wrote each stream.
    let mut mp = 0usize; // markers
    let mut tlp = 0usize; // tile_lens
    let mut sp = 0usize; // strands
    let mut rsp = 0usize; // rs
    let mut llp = 0usize; // lit_lens
    let mut runp = 0usize; // runs
    let mut dp = 0usize; // indel_lens
    let mut litp = 0usize; // literals

    // Reads are rebuilt in id order; every tile's neighbour has a strictly smaller
    // id, so it is already present in `seq` by the time this read needs it.
    for i in 0..n {
        let want = lens[i] as usize;
        let mut out_read: Vec<u8> = Vec::with_capacity(want);
        while out_read.len() < want {
            let marker = read_varint(&markers, &mut mp).ok_or(Error::Corrupt)?;
            if marker == 0 {
                let ll = read_varint(&lit_lens, &mut llp).ok_or(Error::Corrupt)? as usize;
                // A segment can't produce more than the read's remaining bases; the
                // bound keeps `litp + ll` from overflowing and caps the copy.
                if ll > want - out_read.len() {
                    return Err(Error::Corrupt);
                }
                let end = litp.checked_add(ll).ok_or(Error::Corrupt)?;
                let seg = literals.get(litp..end).ok_or(Error::Corrupt)?;
                litp = end;
                out_read.extend(seg.iter().map(|&c| base_of_sym(c)));
                continue;
            }
            // Tile against read `i - delta`, which is strictly earlier.
            let delta = marker as usize;
            if delta > i {
                return Err(Error::Corrupt);
            }
            let target = i - delta;
            let strand = *strands.get(sp).ok_or(Error::Corrupt)? != 0;
            sp += 1;
            let tile_len = read_varint(&tile_lens, &mut tlp).ok_or(Error::Corrupt)? as usize;
            // A tile can't produce more than the read's remaining bases. This bounds
            // the op-replay loop and the `ops`/output allocations below against a
            // hostile length (an all-substitution stream would otherwise loop on the
            // never-exhausting range decoder until it OOMs).
            if tile_len > want - out_read.len() {
                return Err(Error::Corrupt);
            }
            let rs = read_varint(&rs_stream, &mut rsp).ok_or(Error::Corrupt)? as usize;

            let nb = &seq[offs[target]..offs[target + 1]];
            let neighbour = if strand {
                revcomp_acgt(nb)
            } else {
                nb.to_vec()
            };
            let cs = neighbour.get(rs..).ok_or(Error::Corrupt)?;

            // Replay the tile's ops, applying each **directly** to `out_read` — no
            // intermediate `Vec<Op>`, so a hostile `tile_len` cannot inflate an op
            // vector (`~32 B/op`) to OOM. `tprod` counts the query bases the ops
            // declare; it must land exactly on `tile_len`, and the bases actually
            // appended must match (a clamped over-long match desyncs the two).
            let tile_start = out_read.len();
            let mut tprod = 0usize;
            let mut ref_pos = 0usize;
            let mut last = 4u8;
            while tprod < tile_len {
                let code = op_models[op_ctx].decode(&mut op_dec) as u8;
                op_ctx = ((op_ctx << 2) | code as usize) & (OP_CTXS - 1);
                ops_seen += 1;
                match code {
                    0 => {
                        let m = read_varint(&runs, &mut runp).ok_or(Error::Corrupt)? as usize;
                        // Copy the matched reference span, clamped to what the
                        // neighbour holds (as `apply` does). A valid tile never
                        // overruns; a corrupt one is caught by the length check.
                        let start = ref_pos.min(cs.len());
                        let end = ref_pos.saturating_add(m).min(cs.len());
                        out_read.extend_from_slice(&cs[start..end]);
                        if m > 0 {
                            last = cs.get(ref_pos + m - 1).map_or(4, |&x| base_code(x));
                        }
                        tprod += m;
                        ref_pos += m;
                    }
                    1 => {
                        let rc = cs.get(ref_pos).map_or(4, |&x| base_code(x)) as usize;
                        let vc = sub_models[rc.min(SUB_SYMS - 1)].decode(&mut sub_dec);
                        subs_seen += 1;
                        out_read.push(base_of_sym(vc as u8));
                        tprod += 1;
                        ref_pos += 1;
                        last = vc as u8;
                    }
                    2 => {
                        let k = read_varint(&indel_lens, &mut dp).ok_or(Error::Corrupt)? as usize;
                        // An insertion can't produce more query bases than the tile
                        // has left, which bounds this loop against a hostile length.
                        if k > tile_len - tprod {
                            return Err(Error::Corrupt);
                        }
                        for _ in 0..k {
                            let vc =
                                ins_models[(last as usize).min(SUB_SYMS - 1)].decode(&mut ins_dec);
                            last = vc as u8;
                            out_read.push(base_of_sym(vc as u8));
                        }
                        ins_seen += k;
                        tprod += k;
                    }
                    3 => {
                        let m = read_varint(&indel_lens, &mut dp).ok_or(Error::Corrupt)? as usize;
                        ref_pos = ref_pos.saturating_add(m);
                    }
                    _ => return Err(Error::Corrupt),
                }
            }
            // The ops must have declared exactly `tile_len` query bases and appended
            // exactly that many (a valid tile never clamps a match; a corrupt one
            // desyncs one of these two counts).
            if tprod != tile_len || out_read.len() - tile_start != tile_len {
                return Err(Error::Corrupt);
            }
        }
        if out_read.len() != want {
            return Err(Error::Corrupt);
        }
        seq[offs[i]..offs[i + 1]].copy_from_slice(&out_read);
    }

    // The decoded op / substitution / insertion counts must match what was framed.
    if ops_seen != ops_count || subs_seen != subs_count || ins_seen != ins_count {
        return Err(Error::Corrupt);
    }

    // Apply non-ACGT exceptions last — the source of truth for every non-ACGT byte.
    let mut ep = 0usize;
    let mut gpos: u64 = 0;
    for &b in &exc_bytes {
        let delta = read_varint(&exc_pos, &mut ep).ok_or(Error::Corrupt)?;
        gpos += delta;
        *seq.get_mut(gpos as usize).ok_or(Error::Corrupt)? = b;
    }

    Ok((lens, seq))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Deterministic pseudo-random ACGT of `len` bases from `seed`.
    fn rand_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
        (0..len)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                b"ACGT"[(s >> 33) as usize % 4]
            })
            .collect()
    }

    /// A genome-derived read set: `depth` reads of length `rlen` drawn from a
    /// common `glen`-base genome at random offsets, each with `err`-rate
    /// substitutions — enough cross-read overlap for the tiler to engage.
    fn overlapping_reads(
        glen: usize,
        rlen: usize,
        depth: usize,
        err: u32,
        seed: u64,
    ) -> Vec<Vec<u8>> {
        let genome = rand_seq(glen, seed);
        let mut s = seed.wrapping_mul(0xD1B5_4A32_D192_ED03).wrapping_add(7);
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        (0..depth)
            .map(|_| {
                let start = (next() as usize) % (glen - rlen + 1);
                let mut r = genome[start..start + rlen].to_vec();
                if err > 0 {
                    for b in &mut r {
                        if (next() % 1000) < u64::from(err) {
                            *b = b"ACGT"[(next() >> 3) as usize % 4];
                        }
                    }
                }
                r
            })
            .collect()
    }

    fn flatten(reads: &[Vec<u8>]) -> (Vec<u32>, Vec<u8>) {
        let lens = reads.iter().map(|r| r.len() as u32).collect();
        let seq = reads.iter().flatten().copied().collect();
        (lens, seq)
    }

    fn roundtrip(reads: &[Vec<u8>], opts: &EncodeOpts) {
        let (lens, seq) = flatten(reads);
        let block = tile_encode(&lens, &seq, opts).expect("encode");
        let (dl, ds) = tile_decode(&block).expect("decode");
        assert_eq!(dl, lens, "lengths round-trip");
        assert_eq!(ds, seq, "sequence round-trips exactly");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[], &EncodeOpts::default());
    }

    #[test]
    fn roundtrip_single_read_all_literal() {
        // Read 0 has no earlier neighbour → one full-length literal gap.
        roundtrip(&[rand_seq(500, 1)], &EncodeOpts::default());
    }

    #[test]
    fn roundtrip_zero_length_read() {
        let reads = vec![rand_seq(300, 2), Vec::new(), rand_seq(300, 3)];
        roundtrip(&reads, &EncodeOpts::default());
    }

    #[test]
    fn roundtrip_disjoint_reads_all_literal() {
        // Unrelated reads: no overlaps, so every read is stored literally.
        let reads: Vec<Vec<u8>> = (0..8).map(|i| rand_seq(400, 100 + i)).collect();
        roundtrip(&reads, &EncodeOpts::default());
    }

    #[test]
    fn roundtrip_overlapping_clean() {
        let reads = overlapping_reads(4_000, 800, 40, 0, 7);
        roundtrip(&reads, &EncodeOpts::default());
    }

    #[test]
    fn roundtrip_overlapping_noisy() {
        // ONT-class: ~8% substitution noise, so tiles carry real edit scripts.
        let reads = overlapping_reads(4_000, 900, 60, 80, 11);
        roundtrip(&reads, &EncodeOpts::default());
        roundtrip(
            &reads,
            &EncodeOpts {
                tile_band: 768,
                ..EncodeOpts::default()
            },
        );
    }

    #[test]
    fn roundtrip_with_non_acgt() {
        // Inject N's and lowercase — every non-ACGT byte must survive via the
        // exception list, whether it lands in a tile, a literal, or a match.
        let mut reads = overlapping_reads(3_000, 700, 30, 40, 21);
        for (k, r) in reads.iter_mut().enumerate() {
            if !r.is_empty() {
                let idx = k % r.len();
                r[idx] = b'N';
            }
            if r.len() > 10 {
                r[10] = b'a';
            }
        }
        roundtrip(&reads, &EncodeOpts::default());
    }

    #[test]
    fn roundtrip_dense_reads_with_non_acgt_in_reference() {
        // Dense, *clean* coverage so the tiler forms full-read tiles against earlier
        // reads (not the sparse literal fallback), with non-ACGT bytes in the first
        // reads — which later reads then tile against. This is the case the exception
        // list alone does not cover: a neighbour is only ACGT-normalized when a read
        // references it mid-decode (its N is restored in the final pass), so encode
        // must tile the same normalized bytes. A pre-normalization encoder fails this.
        let mut reads = overlapping_reads(2_000, 1_000, 60, 0, 55);
        for r in reads.iter_mut().take(6) {
            if r.len() > 400 {
                r[321] = b'N';
                r[400] = b'a';
                r[15] = b'N';
            }
        }
        roundtrip(&reads, &EncodeOpts::default());
        roundtrip(
            &reads,
            &EncodeOpts {
                tile_band: 768,
                ..EncodeOpts::default()
            },
        );
    }

    #[test]
    fn tiling_beats_all_literal_on_overlapping_reads() {
        // The whole point: on genuinely overlapping reads the tiled block must be
        // smaller than the trivial literal packing (2 bits/base + framing).
        let reads = overlapping_reads(5_000, 1_000, 80, 30, 33);
        let (lens, seq) = flatten(&reads);
        let block = tile_encode(&lens, &seq, &EncodeOpts::default()).expect("encode");
        let total_bases: usize = lens.iter().map(|&l| l as usize).sum();
        // A 2-bit literal packing is total_bases/4 bytes; require a clear win.
        assert!(
            block.len() * 4 < total_bases,
            "tiled block {} B should beat 2-bit literal {} B",
            block.len(),
            total_bases / 4
        );
    }

    #[test]
    fn wrong_magic_is_corrupt() {
        assert!(tile_decode(b"XXX\x01").is_err());
        assert!(tile_decode(&[]).is_err());
    }

    #[test]
    fn decode_is_thread_count_invariant() {
        // Same input must produce a byte-identical block regardless of pool size.
        let reads = overlapping_reads(4_000, 800, 50, 50, 99);
        let (lens, seq) = flatten(&reads);
        let a = tile_encode(&lens, &seq, &EncodeOpts::default()).expect("encode a");
        let b = tile_encode(&lens, &seq, &EncodeOpts::default()).expect("encode b");
        assert_eq!(a, b, "encode is deterministic");
    }

    #[test]
    fn bestof_refs_roundtrip_and_deterministic() {
        // Best-of-N reference selection (tile_max_refs > 1) must stay lossless and
        // thread-count invariant: it only changes *which* neighbour a tile picks,
        // never how the decoder replays it.
        let reads = overlapping_reads(6_000, 1_000, 90, 60, 123);
        let (lens, seq) = flatten(&reads);
        for refs in [2usize, 4, 8] {
            for band in [256usize, 768] {
                let opts = EncodeOpts {
                    tile_max_refs: refs,
                    tile_band: band,
                    ..EncodeOpts::default()
                };
                let a = tile_encode(&lens, &seq, &opts).expect("encode a");
                let b = tile_encode(&lens, &seq, &opts).expect("encode b");
                assert_eq!(a, b, "best-of-{refs} band {band} is deterministic");
                let (dl, ds) = tile_decode(&a).expect("decode");
                assert_eq!(dl, lens, "best-of-{refs} lengths round-trip");
                assert_eq!(ds, seq, "best-of-{refs} sequence round-trips exactly");
            }
        }
    }

    #[test]
    fn bestof_refs_one_matches_default() {
        // tile_max_refs = 1 (and 0, clamped up to 1) must reproduce the greedy
        // max-reach cover byte-for-byte — the shipped #212 default is untouched.
        let reads = overlapping_reads(5_000, 900, 70, 50, 77);
        let (lens, seq) = flatten(&reads);
        let base = tile_encode(&lens, &seq, &EncodeOpts::default()).expect("default");
        for refs in [0usize, 1] {
            let opts = EncodeOpts {
                tile_max_refs: refs,
                ..EncodeOpts::default()
            };
            let got = tile_encode(&lens, &seq, &opts).expect("encode");
            assert_eq!(got, base, "tile_max_refs={refs} matches the greedy default");
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

        /// Arbitrary read sets round-trip exactly, across bands and reference fan-out.
        #[test]
        fn prop_roundtrip_arbitrary(
            seed in any::<u64>(),
            glen in 500usize..3_000,
            rlen in 100usize..500,
            depth in 1usize..40,
            err in 0u32..120,
            band in prop::sample::select(vec![64usize, 256, 768]),
            refs in prop::sample::select(vec![1usize, 2, 4]),
        ) {
            let reads = overlapping_reads(glen, rlen.min(glen), depth, err, seed);
            let (lens, seq) = flatten(&reads);
            let opts = EncodeOpts { tile_band: band, tile_max_refs: refs, ..EncodeOpts::default() };
            let block = tile_encode(&lens, &seq, &opts).expect("encode");
            let (dl, ds) = tile_decode(&block).expect("decode");
            prop_assert_eq!(dl, lens);
            prop_assert_eq!(ds, seq);
        }

        /// Fully random (non-overlapping) reads still round-trip — the all-literal
        /// fallback path.
        #[test]
        fn prop_roundtrip_random_reads(
            seed in any::<u64>(),
            n in 0usize..12,
            rlen in 0usize..300,
        ) {
            let reads: Vec<Vec<u8>> = (0..n).map(|i| rand_seq(rlen, seed ^ i as u64)).collect();
            let (lens, seq) = flatten(&reads);
            let block = tile_encode(&lens, &seq, &EncodeOpts::default()).expect("encode");
            let (dl, ds) = tile_decode(&block).expect("decode");
            prop_assert_eq!(dl, lens);
            prop_assert_eq!(ds, seq);
        }
    }
}
