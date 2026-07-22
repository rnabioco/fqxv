//! The lossless container codec: `encode` reads → a self-contained byte block,
//! `decode` the block → the exact reads back, in the original order.
//!
//! [`encode`] runs the crate's pipeline (subsample → overlaps → layout →
//! consensus → place every read → banded edit script) and, unlike
//! `examples/encode.rs` (which only *measures* the stream sizes), serialises
//! every stream the decoder needs into one framed blob:
//!
//! - the per-contig **reference** consensi (stored raw, entropy-coded);
//! - a **manifest** — for each read, in encode order, its original id and
//!   whether it was coded as an edit script or as a standalone literal — which
//!   is what lets [`decode`] restore the *original* read order from the
//!   reference-grouped streams;
//! - the per-read **flip / placement / op / run / sub / indel** streams;
//! - **literals** for reads that landed on no reference;
//! - a **non-ACGT exception list** `(global base index, byte)`, applied last, so
//!   losslessness never depends on how the aligner coded an `N` or a lowercase
//!   base.
//!
//! Read lengths terminate each read's op run on decode (ops are consumed until
//! the produced length equals `lens[r]`), so the format self-describes lengths
//! the same way [`fqxv_seq`] does.
//!
//! ## Shared reference (whole-file)
//!
//! The reference above is by far the largest single stream, and at high coverage
//! the container splits a file into several blocks that each self-assemble and
//! **re-store the same reference**. [`build_reference`] + [`encode_against`] hoist
//! it out: assemble one [`Reference`] over the whole file, store it **once** (the
//! container writes it as a shared frame), and code every block's reads against
//! that frozen frame with [`encode_against`] (a distinct bitstream version), which
//! is byte-for-byte the standalone block *minus* the reference streams.
//! [`decode_against`] inverts it given the same reference. Placement is per-read
//! against a frozen frame, so a read codes identically regardless of which block
//! it lands in — no block-boundary penalty.

use std::collections::HashSet;

use rayon::prelude::*;

use fqxv_bytes::{read_varint, write_varint};
use fqxv_dna::{base_of_sym, is_acgt, revcomp_acgt};

use crate::anchorgap::{anchor_chain, anchorgap_build};
use crate::{
    Alignment, Anchored, ChainOpts, ConsensusOpts, Error, Index, LayoutOpts, MAX_BASES_PER_BYTE,
    Op, Repeat, Sketch, align_banded, apply, consensus, find_overlaps, layout, place_all,
    wfa_align_opt,
};

/// Above this drift-derived band the banded DP is the faster aligner; at or below
/// it the read is similar enough to its reference that score-proportional WFA
/// wins decisively (≈35× at HiFi's ~0.5% error, where a 13 kb read needs ~60
/// edits). The threshold is the measured WFA-vs-DP crossover (~5% divergence,
/// band ≈ 32); see `examples/wfa_bench.rs`. Above it — ONT-class reads — WFA's
/// O(s²) time and memory both lose, so the DP stays.
const WFA_CODING_MAX_BAND: usize = 32;

/// Align a read to its reference consensus for coding, picking the aligner by the
/// read's expected divergence (its chain-drift `band`): WFA for the low-drift
/// (HiFi-class) reads, the banded DP otherwise. WFA returns a different
/// equal-cost path, so a divergence-gated encoder is NOT byte-identical to the
/// pure-DP one — but both round-trip, and the decoder is agnostic to which
/// produced the edit script. On the rare read whose small band under-estimated
/// its true divergence, WFA caps out and we fall back to the DP for a real
/// alignment rather than take the bounded `[Del, Ins]` rewrite.
pub(crate) fn align_for_coding(refr: &[u8], query: &[u8], band: usize) -> Alignment {
    if band <= WFA_CODING_MAX_BAND {
        // A low-drift read's true edit distance is ~O(band); cap generously above
        // it so a genuinely similar read never caps, while a mis-estimated one
        // bails after O(cap²) work instead of exploring unbounded wavefronts.
        let cap = (band as u32 * 4).max(64);
        if let Some(a) = wfa_align_opt(refr, query, cap) {
            return a;
        }
    }
    align_banded(refr, query, band)
}

/// Format tag for a `fqxv-lroverlap` sequence block.
const MAGIC: [u8; 3] = *b"LRO";
/// Bitstream version for a self-contained block ([`encode`]/[`decode`]): the
/// reference consensi are stored inline. Bump on any layout change (nothing on
/// disk is stable yet).
const VERSION: u8 = 1;
/// Bitstream version for a block coded against a shared [`Reference`]
/// ([`encode_against`]/[`decode_against`]): identical to [`VERSION`] except the
/// reference consensi are NOT stored — the caller supplies them.
const VERSION_SHARED: u8 = 2;

/// Format tag for a standalone shared [`Reference`] frame (see
/// [`Reference::encode`]).
const REF_MAGIC: [u8; 3] = *b"LRR";
/// Bitstream version of a [`Reference`] frame.
const REF_VERSION: u8 = 1;

/// Coverage the layout is fed after subsampling. The layout is excellent at
/// ~40× and starves past it (see `layout.rs`); every read is then placed against
/// the consensus it produces, so the reference is amortised over the whole
/// block. Deriving the stride from this target keeps it out of the caller's API.
const TARGET_LAYOUT_COVERAGE: u64 = 40;

/// Hard ceiling on the layout subsample, independent of the coverage estimate.
///
/// The subsample exists only to build reference consensi — a few tens of ×
/// coverage of each locus is ample and the layout starves past ~40× anyway, so
/// no locus needs more than this many reads. The estimate that sizes the stride
/// (`derive_stride`) collapses on high-redundancy input — an amplicon block is
/// one locus copied ~150k times, whose sequencing errors explode the
/// distinct-minimizer genome estimate, so the estimated coverage floors and the
/// stride falls to 1. Every read then enters the overlap→layout→consensus
/// pipeline: the search runs for minutes and the fragmented layout's placement
/// matrix OOM-kills the block (#139). This cap bounds the subsample regardless.
/// It sits far above what any normal block subsamples to (target coverage of a
/// real genome), so those blocks keep their exact stride and stay byte-identical.
const MAX_LAYOUT_SUBSAMPLE: usize = 32_768;

/// One extra reference base past a read's own length, so the banded aligner has
/// somewhere to place a trailing deletion.
pub(crate) const REF_TAIL: usize = 64;

/// Options for [`encode`]. None affect [`decode`] — the block self-describes.
#[derive(Debug, Clone)]
pub struct EncodeOpts {
    /// Minimizer sketch; pick [`Sketch::hifi`] for low-error long reads,
    /// [`Sketch::ont`] for noisy ones. Affects ratio and speed, not correctness.
    pub sketch: Sketch,
    /// Layout subsample stride. `None` derives it from estimated coverage so the
    /// layout sees ~40×; `Some(1)` disables subsampling.
    pub stride: Option<usize>,
    /// Slack added to each read's own chain drift when sizing its alignment band.
    pub band_margin: usize,
    /// Hard cap on the band, so one pathological chain cannot allocate a
    /// quadratic DP table.
    pub band_cap: usize,
    /// Alignment band half-width for the multi-reference tiler
    /// ([`tile_encode`](crate::tile_encode)). Each read is coded against *raw*
    /// earlier neighbour reads (not a voted consensus), so a tile sits at the
    /// sum of both reads' error (~13.5% on ONT) and needs a wider band than the
    /// consensus coder. 256 is the cheaper measured operating point (1.142
    /// b/base on ONT); 768 trades encode time for ~1.056. Affects ratio and
    /// speed only — the block self-describes, so decode is unaffected.
    pub tile_band: usize,
    /// How many candidate reference reads the tiler tries per tile before
    /// committing ([`tile_encode`](crate::tile_encode)). At high (ONT) coverage
    /// many earlier reads overlap the same span, each with an independent error
    /// pattern; aligning the furthest-reaching `tile_max_refs` of them and
    /// keeping the cheapest edit script picks a reference that happens to agree
    /// with the query at more positions, the way CoLoRd chooses its anchor. `1`
    /// is the plain greedy max-reach cover (one alignment per tile, the cheapest
    /// operating point); higher trades ~linearly more alignment time for ratio.
    /// Affects ratio and speed only — the block self-describes, so decode is
    /// unaffected.
    pub tile_max_refs: usize,
}

impl Default for EncodeOpts {
    fn default() -> Self {
        Self {
            sketch: Sketch::ont(),
            stride: None,
            band_margin: 32,
            band_cap: 2048,
            tile_band: 256,
            tile_max_refs: 1,
        }
    }
}

/// A frozen cross-read reference: the per-contig consensus sequences reads are
/// coded against. Built once by [`build_reference`] (over a block, or the whole
/// file) and shared immutably — it is `Send + Sync` — across [`encode_against`] /
/// [`decode_against`] calls. Storing it once and coding every block against it is
/// what removes the redundant per-block reference copies (issue #168).
#[derive(Debug, Clone, Default)]
pub struct Reference {
    /// Consensus contigs in index order (the order [`place_all`] assigns targets).
    consensi: Vec<Vec<u8>>,
}

impl Reference {
    /// The consensus contigs, in index order.
    #[must_use]
    pub fn consensi(&self) -> &[Vec<u8>] {
        &self.consensi
    }

    /// Number of consensus contigs.
    #[must_use]
    pub fn len(&self) -> usize {
        self.consensi.len()
    }

    /// Whether the reference has no contigs (assembly found no shared locus).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.consensi.is_empty()
    }

    /// Total consensus bases across all contigs.
    #[must_use]
    pub fn total_bases(&self) -> usize {
        self.consensi.iter().map(Vec::len).sum()
    }

    /// Serialise the reference to a self-contained, entropy-coded frame: contig
    /// count, per-contig lengths, and the concatenated consensus bases. The
    /// container stores this once per file (see issue #168); [`Reference::decode`]
    /// inverts it. This is byte-for-byte the reference portion a self-contained
    /// block stores inline, just lifted out with its own magic/version.
    pub fn encode(&self) -> Result<Vec<u8>, Error> {
        let mut out = Vec::new();
        out.extend_from_slice(&REF_MAGIC);
        out.push(REF_VERSION);
        write_varint(&mut out, self.consensi.len() as u64);
        let mut ref_lens_raw = Vec::new();
        let mut ref_bytes = Vec::new();
        for cs in &self.consensi {
            write_varint(&mut ref_lens_raw, cs.len() as u64);
            ref_bytes.extend_from_slice(cs);
        }
        put_stream(&mut out, &ref_lens_raw)?;
        put_stream(&mut out, &ref_bytes)?;
        Ok(out)
    }

    /// Decode a frame produced by [`Reference::encode`].
    pub fn decode(src: &[u8]) -> Result<Self, Error> {
        let mut pos = 0usize;
        if src.get(..REF_MAGIC.len()) != Some(&REF_MAGIC) {
            return Err(Error::Corrupt);
        }
        pos += REF_MAGIC.len();
        let version = *src.get(pos).ok_or(Error::Corrupt)?;
        pos += 1;
        if version != REF_VERSION {
            return Err(Error::Corrupt);
        }
        let n_contigs = read_varint(src, &mut pos).ok_or(Error::Corrupt)? as usize;
        let ref_lens_raw = get_stream(src, &mut pos)?;
        let ref_bytes = get_stream(src, &mut pos)?;
        // `n_contigs` is untrusted, so do NOT pre-reserve it (a crafted count is an
        // allocation bomb); each contig length costs at least one varint byte, so
        // the loop errors out of `ref_lens_raw` long before an absurd count is
        // reached. See the reorder decode path's `reserve_checked` for the same
        // rule.
        let mut consensi: Vec<Vec<u8>> = Vec::new();
        let mut rlp = 0usize;
        let mut roff = 0usize;
        for _ in 0..n_contigs {
            let rl = read_varint(&ref_lens_raw, &mut rlp).ok_or(Error::Corrupt)? as usize;
            let end = roff.checked_add(rl).ok_or(Error::Corrupt)?;
            consensi.push(ref_bytes.get(roff..end).ok_or(Error::Corrupt)?.to_vec());
            roff = end;
        }
        if roff != ref_bytes.len() {
            return Err(Error::Corrupt);
        }
        Ok(Reference { consensi })
    }
}

/// Strict 2-bit base code: `A C G T` → `0 1 2 3`, everything else → `4` (a
/// placeholder the exception list overwrites on decode). Deliberately *not*
/// case-folding — a lowercase base is preserved verbatim through the exceptions.
#[inline]
pub(crate) fn base_code(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// Entropy-code `raw` with the smaller of order-0 / order-1 rANS and frame it as
/// `[varint raw_len][varint enc_len][enc bytes]`. Empty streams store only a
/// zero length. The order tag lives inside the rANS header, so [`get_stream`]
/// needs no discriminator of its own.
pub(crate) fn put_stream(out: &mut Vec<u8>, raw: &[u8]) -> Result<(), Error> {
    write_varint(out, raw.len() as u64);
    if raw.is_empty() {
        return Ok(());
    }
    let e0 = fqxv_rans::encode(raw, fqxv_rans::Order::Zero).map_err(|_| Error::Corrupt)?;
    let e1 = fqxv_rans::encode(raw, fqxv_rans::Order::One).map_err(|_| Error::Corrupt)?;
    let enc = if e1.len() <= e0.len() { e1 } else { e0 };
    write_varint(out, enc.len() as u64);
    out.extend_from_slice(&enc);
    Ok(())
}

/// Substitution alphabet: ACGT (codes 0-3) plus a non-ACGT placeholder (4).
pub(crate) const SUB_SYMS: usize = 5;

/// Range-code the substituted bases conditioned on the reference base each replaced.
/// ONT substitutions are strongly reference-base-dependent, so a per-ref-base
/// adaptive model beats the flat order-0/1 rANS (CoLoRd conditions substitutions on
/// the reference base similarly). The models start uniform and adapt in stream order;
/// [`decode`] rebuilds the same models with the same contexts (the consensus base at
/// each position), so it round-trips exactly.
pub(crate) fn encode_subs(subs: &[u8], sub_ctx: &[u8]) -> Vec<u8> {
    let mut enc = fqxv_range::Encoder::new();
    let mut models: [fqxv_range::SimpleModel<SUB_SYMS>; SUB_SYMS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    for (&v, &c) in subs.iter().zip(sub_ctx) {
        let ctx = (c as usize).min(SUB_SYMS - 1);
        models[ctx].encode(&mut enc, (v as usize).min(SUB_SYMS - 1));
    }
    enc.finish()
}

/// Number of order-3 op-history contexts (last 3 ops × 2 bits).
pub(crate) const OP_CTXS: usize = 64;

/// Range-code the op-type stream (Match/Sub/Ins/Del = 0..3) with an order-3 context
/// over the three preceding op types. Op sequences are strongly self-correlated
/// (Match follows Match; homopolymer indels cluster), so this beats the flat rANS.
/// The rolling context is derived from the ops themselves, so decode rebuilds it
/// identically as it decodes each op during reconstruction.
pub(crate) fn encode_ops(ops: &[u8]) -> Vec<u8> {
    let mut enc = fqxv_range::Encoder::new();
    let mut models: [fqxv_range::SimpleModel<4>; OP_CTXS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut ctx = 0usize;
    for &o in ops {
        let sym = (o as usize).min(3);
        models[ctx].encode(&mut enc, sym);
        ctx = ((ctx << 2) | sym) & (OP_CTXS - 1);
    }
    enc.finish()
}

pub(crate) fn get_stream(src: &[u8], pos: &mut usize) -> Result<Vec<u8>, Error> {
    let raw_len = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    if raw_len == 0 {
        return Ok(Vec::new());
    }
    let enc_len = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let end = pos.checked_add(enc_len).ok_or(Error::Corrupt)?;
    let enc = src.get(*pos..end).ok_or(Error::Corrupt)?;
    *pos = end;
    let out = fqxv_rans::decode_bounded(enc, raw_len).map_err(|_| Error::Corrupt)?;
    if out.len() != raw_len {
        return Err(Error::Corrupt);
    }
    Ok(out)
}

/// A single read's contribution to the seven edit streams, plus its literal
/// bases when it lands on no reference.
#[derive(Default)]
struct Streams {
    ops: Vec<u8>,
    runs: Vec<u8>,
    subs: Vec<u8>,
    /// Parallel to `subs`: the reference (consensus) base code at each substitution,
    /// used as the entropy-coding context (ONT substitutions are ref-base-dependent).
    sub_ctx: Vec<u8>,
    ins_bases: Vec<u8>,
    /// Parallel to `ins_bases`: the preceding read base at each inserted base, the
    /// entropy context (ONT insertions are mostly homopolymer extensions).
    ins_ctx: Vec<u8>,
    indel_lens: Vec<u8>,
    placements: Vec<u8>,
    flips: Vec<u8>,
    literals: Vec<u8>,
    /// Manifest, encode order: original read id of every read in this group.
    read_ids: Vec<u8>,
    /// Manifest, encode order: `0` = edit-script coded, `1` = literal.
    kinds: Vec<u8>,
}

impl Streams {
    fn merge(&mut self, o: &Streams) {
        self.ops.extend_from_slice(&o.ops);
        self.runs.extend_from_slice(&o.runs);
        self.subs.extend_from_slice(&o.subs);
        self.sub_ctx.extend_from_slice(&o.sub_ctx);
        self.ins_bases.extend_from_slice(&o.ins_bases);
        self.ins_ctx.extend_from_slice(&o.ins_ctx);
        self.indel_lens.extend_from_slice(&o.indel_lens);
        self.placements.extend_from_slice(&o.placements);
        self.flips.extend_from_slice(&o.flips);
        self.literals.extend_from_slice(&o.literals);
        self.read_ids.extend_from_slice(&o.read_ids);
        self.kinds.extend_from_slice(&o.kinds);
    }
}

/// Everything [`serialize`] needs beyond the reference and read lengths: the
/// per-read edit streams grouped by contig (plus the trailing orphan group), the
/// group counts, and the non-ACGT exception list. Produced by [`code_against`].
struct Coded {
    /// One count per contig (reads placed on it), then the orphan group count.
    counts: Vec<u64>,
    /// All edit streams, concatenated in group order.
    all: Streams,
    /// Non-ACGT exception positions (delta-coded) and their bytes.
    exc_pos: Vec<u8>,
    exc_bytes: Vec<u8>,
}

/// Estimate the block's fold coverage from its distinct-minimizer count
/// (genome ≈ distinct × (w+1)/2, minimizer density 2/(w+1)) and turn it into a
/// layout stride targeting [`TARGET_LAYOUT_COVERAGE`]. Repeats bias the genome
/// estimate low and the stride high; both are bounded and only affect ratio.
fn derive_stride(sketch: Sketch, reads: &[&[u8]], total_bases: u64) -> usize {
    // The distinct minimizer count over every read — a full serial sketch plus a
    // hash insert per minimizer, which dominated encode setup on deep blocks (a
    // ~1-core stall of several seconds). Each read's minimizers are independent,
    // so fold them into per-worker sets on the rayon pool and union the sets. Only
    // the set *cardinality* feeds the stride, and set union is commutative, so the
    // result is identical for any worker count — the derived stride, and thus the
    // whole encoding, stays byte-identical regardless of thread count.
    let distinct = reads
        .par_iter()
        .fold(HashSet::<u64>::new, |mut s, r| {
            for m in sketch.seeds(r) {
                s.insert(m.hash);
            }
            s
        })
        .reduce(HashSet::<u64>::new, |a, b| {
            // Extend the larger set with the smaller to minimize rehashing.
            let (mut big, small) = if a.len() >= b.len() { (a, b) } else { (b, a) };
            big.extend(small);
            big
        });
    if distinct.is_empty() {
        return 1;
    }
    let genome = (distinct.len() as u64 * (sketch.w as u64 + 1) / 2).max(1);
    let coverage = total_bases / genome;
    ((coverage / TARGET_LAYOUT_COVERAGE) as usize).max(1)
}

/// Cumulative byte offsets into `seq` for reads of the given `lens`, validating
/// that the lengths sum to `seq.len()`.
pub(crate) fn compute_offs(lens: &[u32], seq_len: usize) -> Result<Vec<usize>, Error> {
    let mut offs = Vec::with_capacity(lens.len() + 1);
    let mut acc = 0usize;
    for &l in lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);
    if acc != seq_len {
        return Err(Error::Corrupt);
    }
    Ok(offs)
}

/// Assemble a [`Reference`] from long reads: subsample the layout to a target
/// coverage, find overlaps, lay out contigs, and vote a consensus per multi-read
/// contig. This is the assembly half of [`encode`],
/// hoisted so the container can build one reference over the whole file and share
/// it across blocks (issue #168). Pure function of the input (and thread-count
/// invariant), so the shared reference is deterministic.
pub fn build_reference(lens: &[u32], seq: &[u8], opts: &EncodeOpts) -> Result<Reference, Error> {
    let sketch = opts.sketch;
    let n = lens.len();
    let offs = compute_offs(lens, seq.len())?;
    let read_at = |r: usize| &seq[offs[r]..offs[r + 1]];
    let all_refs: Vec<&[u8]> = (0..n).map(read_at).collect();
    let total_bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();

    // ---- subsample the layout ------------------------------------------
    let stride = match opts.stride {
        Some(s) => s.max(1),
        None => derive_stride(sketch, &all_refs, total_bases),
    };
    // Raise the stride if needed so the subsample cannot exceed the ceiling; a
    // stride of `ceil(n / cap)` yields at most `cap` subsampled reads (#139).
    let stride = stride.max(n.div_ceil(MAX_LAYOUT_SUBSAMPLE.max(1)));
    let sub: Vec<u32> = (0..n as u32).step_by(stride).collect();
    let sub_lens: Vec<u32> = sub.iter().map(|&r| lens[r as usize]).collect();
    let mut sub_seq: Vec<u8> = Vec::new();
    for &r in &sub {
        sub_seq.extend_from_slice(read_at(r as usize));
    }
    let mut sub_offs = Vec::with_capacity(sub.len() + 1);
    let mut acc = 0usize;
    for &l in &sub_lens {
        sub_offs.push(acc);
        acc += l as usize;
    }
    sub_offs.push(acc);
    let sub_read_at = |i: usize| &sub_seq[sub_offs[i]..sub_offs[i + 1]];

    // ---- overlaps → layout → consensus over the subsample --------------
    let idx = Index::build(&sub_lens, &sub_seq, sketch, Repeat::default())?;
    let ovs: Vec<Vec<_>> = (0..sub.len())
        .into_par_iter()
        .map(|r| find_overlaps(&idx, r as u32, sub_read_at(r), ChainOpts::default()))
        .collect();
    let contigs = layout(&sub_lens, &ovs, LayoutOpts::default());

    let consensi: Vec<Vec<u8>> = contigs
        .par_iter()
        .filter(|c| c.reads.len() > 1)
        .filter_map(|c| {
            let mut by_id: Vec<Vec<u8>> = vec![Vec::new(); sub.len()];
            for p in &c.reads {
                let r = sub_read_at(p.read as usize);
                by_id[p.read as usize] = if p.flip { revcomp_acgt(r) } else { r.to_vec() };
            }
            let cons = consensus(
                c,
                &by_id,
                ConsensusOpts {
                    sketch,
                    ..ConsensusOpts::default()
                },
            );
            if cons.is_empty() {
                None
            } else {
                Some(cons.seq)
            }
        })
        .collect();

    Ok(Reference { consensi })
}

/// Place every read on its best-scoring consensus and code it as an edit script
/// (or a literal orphan). This is the coding half of [`encode`], factored so both
/// the self-contained [`encode`] and the shared-reference [`encode_against`] share
/// it verbatim — the only difference between the two block versions is whether the
/// caller stores `consensi` in the block.
fn code_against(
    consensi: &[Vec<u8>],
    lens: &[u32],
    seq: &[u8],
    offs: &[usize],
    opts: &EncodeOpts,
) -> Result<Coded, Error> {
    let sketch = opts.sketch;
    let n = lens.len();
    let read_at = |r: usize| &seq[offs[r]..offs[r + 1]];
    let all_refs: Vec<&[u8]> = (0..n).map(read_at).collect();

    // Anchor-restricted coding (#226) is the default; `FQXV_LRO_NO_ANCHORGAP=1`
    // forces the whole-window banded DP everywhere for A/B (mirrors the tiler's
    // `FQXV_TILE_NO_ANCHORGAP`). Even on the anchor path, a read whose reference
    // shares too few exact k-mers to chain still falls back to the DP — the anchor
    // path never *loses* coverage, it only avoids re-deriving the exact matches a
    // recovered chain already proved. Read once here, not per read.
    let anchorgap = std::env::var_os("FQXV_LRO_NO_ANCHORGAP").is_none();

    // ---- place every read on its best-scoring reference ----------------
    // One combined index over all consensi, queried once per read, rather than
    // one `place_against` pass per consensus (`consensi × reads` overlap searches
    // — tens of millions on an amplicon block, the encode-time sink of #139).
    // `place_all` returns `(consensus index, placement)` per read; a read on no
    // reference stays `None` and codes standalone.
    let best: Vec<Option<(usize, Anchored)>> =
        place_all(consensi, &all_refs, sketch, ChainOpts::default());

    // Group placed reads by reference, ordered by (offset, read) so delta-coded
    // placements are small and the order is total.
    let mut by_ref: Vec<Vec<(u32, u32)>> = vec![Vec::new(); consensi.len()];
    for (r, b) in best.iter().enumerate() {
        if let Some((ci, a)) = b {
            by_ref[*ci].push((a.offset, r as u32));
        }
    }
    for v in &mut by_ref {
        v.sort_unstable();
    }

    // ---- code each read against its reference --------------------------
    let per_contig: Vec<Streams> = consensi
        .par_iter()
        .enumerate()
        .map(|(ci, cs)| {
            let aligned: Vec<Option<(usize, Vec<Op>)>> = by_ref[ci]
                .par_iter()
                .map(|&(_, r)| {
                    let a = best[r as usize].expect("assigned").1;
                    let raw = read_at(r as usize);
                    let read = if a.flip {
                        revcomp_acgt(raw)
                    } else {
                        raw.to_vec()
                    };
                    let start = a.offset as usize;
                    let end = (start + read.len() + REF_TAIL).min(cs.len());
                    if start >= end {
                        return None;
                    }
                    let band = (a.drift as usize + opts.band_margin).min(opts.band_cap);
                    let refwin = &cs[start..end];
                    // Emit the read's shared exact-k-mer anchors against its
                    // consensus as free copy-runs and align only the inter-anchor
                    // gaps + flanks, instead of a whole-window banded DP. On HiFi
                    // the read shares nearly all of its k-mers with a good
                    // consensus, so the DP area collapses from O(read_len * band)
                    // to the sum of the (tiny) gaps. `anchor_chain` returns `None`
                    // for a read that shares too few k-mers to chain, and the DP
                    // then codes it whole. Encoder-only — the decoder replays
                    // whichever `Op` stream this produced.
                    let ops = match anchorgap
                        .then(|| anchor_chain(refwin, &read, sketch))
                        .flatten()
                    {
                        Some(chain) => {
                            anchorgap_build(refwin, &read, band, &chain, sketch.k as u32)
                        }
                        None => align_for_coding(refwin, &read, band).ops,
                    };
                    Some((start, ops))
                })
                .collect();

            let mut s = Streams::default();
            let mut prev_start: i64 = 0;
            for (idx, slot) in aligned.into_iter().enumerate() {
                let r = by_ref[ci][idx].1 as usize;
                write_varint(&mut s.read_ids, r as u64);
                let Some((start, mut ops)) = slot else {
                    s.kinds.push(1);
                    s.literals.extend(read_at(r).iter().map(|&b| base_code(b)));
                    continue;
                };
                // Trailing deletions consume only reference (the `REF_TAIL`
                // slack), never query, so they carry no information the decoder
                // needs — and it terminates a read's ops when the produced length
                // reaches `lens[r]`, so a trailing Del would desync the stream.
                while matches!(ops.last(), Some(Op::Del(_))) {
                    ops.pop();
                }
                s.kinds.push(0);
                s.flips.push(u8::from(best[r].expect("assigned").1.flip));
                let d = start as i64 - prev_start;
                prev_start = start as i64;
                write_varint(&mut s.placements, fqxv_bytes::zigzag(d));
                // Track the reference position (for the substitution's ref base) and
                // the preceding read base (for the insertion's homopolymer context).
                let mut ref_pos = 0usize;
                let mut last: u8 = 4; // 4 = no preceding base (read start)
                for op in &ops {
                    match op {
                        Op::Match(m) => {
                            s.ops.push(0);
                            write_varint(&mut s.runs, u64::from(*m));
                            let mm = *m as usize;
                            if mm > 0 {
                                last = cs
                                    .get(start + ref_pos + mm - 1)
                                    .map_or(4, |&x| base_code(x));
                            }
                            ref_pos += mm;
                        }
                        Op::Sub(b) => {
                            s.ops.push(1);
                            let bc = base_code(*b);
                            s.subs.push(bc);
                            s.sub_ctx
                                .push(cs.get(start + ref_pos).map_or(4, |&x| base_code(x)));
                            last = bc;
                            ref_pos += 1;
                        }
                        Op::Ins(bs) => {
                            s.ops.push(2);
                            write_varint(&mut s.indel_lens, bs.len() as u64);
                            for &b in bs {
                                let bc = base_code(b);
                                s.ins_ctx.push(last);
                                s.ins_bases.push(bc);
                                last = bc;
                            }
                        }
                        Op::Del(m) => {
                            s.ops.push(3);
                            write_varint(&mut s.indel_lens, u64::from(*m));
                            ref_pos += *m as usize;
                        }
                    }
                }
            }
            s
        })
        .collect();

    // Orphans — reads on no reference — as a final literal-only group.
    let mut orphans = Streams::default();
    let mut orphan_count = 0u64;
    for (r, b) in best.iter().enumerate() {
        if b.is_none() {
            write_varint(&mut orphans.read_ids, r as u64);
            orphans.kinds.push(1);
            orphans
                .literals
                .extend(read_at(r).iter().map(|&x| base_code(x)));
            orphan_count += 1;
        }
    }

    // Group counts (encode order): one per contig, then the orphan group.
    let mut counts: Vec<u64> = by_ref.iter().map(|v| v.len() as u64).collect();
    counts.push(orphan_count);

    let mut all = Streams::default();
    for s in &per_contig {
        all.merge(s);
    }
    all.merge(&orphans);

    // Non-ACGT exceptions, in original order: overwrite the placeholder after
    // reconstruction, so correctness never depends on how a base was coded.
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

    Ok(Coded {
        counts,
        all,
        exc_pos,
        exc_bytes,
    })
}

/// Serialise a coded block. `store_ref` is `Some(consensi)` for a self-contained
/// [`VERSION`] block (reference inline) and `None` for a shared-reference
/// [`VERSION_SHARED`] block (the caller holds the reference). Everything else is
/// identical, so the two block versions differ only by the reference streams.
fn serialize(
    version: u8,
    store_ref: Option<&[Vec<u8>]>,
    lens: &[u32],
    total_bases: u64,
    coded: &Coded,
) -> Result<Vec<u8>, Error> {
    let n = lens.len();
    let mut out = Vec::new();
    out.extend_from_slice(&MAGIC);
    out.push(version);
    write_varint(&mut out, n as u64);
    write_varint(&mut out, total_bases);

    let mut lens_raw = Vec::new();
    for &l in lens {
        write_varint(&mut lens_raw, u64::from(l));
    }
    put_stream(&mut out, &lens_raw)?;

    if let Some(consensi) = store_ref {
        write_varint(&mut out, consensi.len() as u64);
        let mut ref_lens_raw = Vec::new();
        let mut ref_bytes = Vec::new();
        for cs in consensi {
            write_varint(&mut ref_lens_raw, cs.len() as u64);
            ref_bytes.extend_from_slice(cs);
        }
        put_stream(&mut out, &ref_lens_raw)?;
        put_stream(&mut out, &ref_bytes)?;
    }

    let mut counts_raw = Vec::new();
    for &c in &coded.counts {
        write_varint(&mut counts_raw, c);
    }
    put_stream(&mut out, &counts_raw)?;

    put_stream(&mut out, &coded.all.read_ids)?;
    put_stream(&mut out, &coded.all.kinds)?;
    put_stream(&mut out, &coded.all.flips)?;
    put_stream(&mut out, &coded.all.placements)?;
    // ops: range-coded with an order-3 op-history context (see `encode_ops`).
    {
        let blob = encode_ops(&coded.all.ops);
        write_varint(&mut out, coded.all.ops.len() as u64);
        write_varint(&mut out, blob.len() as u64);
        out.extend_from_slice(&blob);
    }
    put_stream(&mut out, &coded.all.runs)?;
    // subs: range-coded conditioned on the reference base (see `encode_subs`).
    {
        let blob = encode_subs(&coded.all.subs, &coded.all.sub_ctx);
        write_varint(&mut out, coded.all.subs.len() as u64);
        write_varint(&mut out, blob.len() as u64);
        out.extend_from_slice(&blob);
    }
    // ins_bases: range-coded conditioned on the preceding read base (homopolymer
    // context). Same 5-context base coder as subs, keyed on `ins_ctx`.
    {
        let blob = encode_subs(&coded.all.ins_bases, &coded.all.ins_ctx);
        write_varint(&mut out, coded.all.ins_bases.len() as u64);
        write_varint(&mut out, blob.len() as u64);
        out.extend_from_slice(&blob);
    }
    put_stream(&mut out, &coded.all.indel_lens)?;
    put_stream(&mut out, &coded.all.literals)?;
    put_stream(&mut out, &coded.exc_pos)?;
    put_stream(&mut out, &coded.exc_bytes)?;

    Ok(out)
}

/// Compress a block of long reads. `lens[r]` is the length of read `r`; `seq` is
/// their bases concatenated in order. Returns a self-contained block that
/// [`decode`] inverts exactly — the reference consensi are stored inline. For the
/// shared-reference variant (reference stored once per file) see
/// [`build_reference`] + [`encode_against`].
pub fn encode(lens: &[u32], seq: &[u8], opts: &EncodeOpts) -> Result<Vec<u8>, Error> {
    let offs = compute_offs(lens, seq.len())?;
    let total_bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    let refr = build_reference(lens, seq, opts)?;
    let coded = code_against(&refr.consensi, lens, seq, &offs, opts)?;
    serialize(VERSION, Some(&refr.consensi), lens, total_bases, &coded)
}

/// Compress a block of long reads against a shared, externally-stored
/// [`Reference`]. The returned block omits the reference streams (so the reference
/// is not re-stored per block); [`decode_against`] inverts it given the same
/// reference. `opts` need not match the reference's build options — placement and
/// coding are pure functions of the reference the caller passes in — though the
/// sketch should match so overlaps are found the same way.
pub fn encode_against(
    reference: &Reference,
    lens: &[u32],
    seq: &[u8],
    opts: &EncodeOpts,
) -> Result<Vec<u8>, Error> {
    let offs = compute_offs(lens, seq.len())?;
    let total_bases: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    let coded = code_against(&reference.consensi, lens, seq, &offs, opts)?;
    serialize(VERSION_SHARED, None, lens, total_bases, &coded)
}

/// Parse the shared block header (`MAGIC`, expected version, `n`, `total_bases`)
/// and read the per-read lengths, returning `(n, total_bases, lens, offs)` with
/// `pos` advanced to the first stream after the lengths.
fn read_header(
    src: &[u8],
    pos: &mut usize,
    expect_version: u8,
) -> Result<(usize, usize, Vec<u32>, Vec<usize>), Error> {
    if src.get(..3) != Some(&MAGIC) {
        return Err(Error::Corrupt);
    }
    *pos += 3;
    let version = *src.get(*pos).ok_or(Error::Corrupt)?;
    *pos += 1;
    if version != expect_version {
        return Err(Error::Corrupt);
    }
    let n = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let total_bases = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    // Reject a header whose declared base count could not possibly fit the coded
    // block before it drives the `vec![0u8; total_bases]` output allocation — and,
    // via each read's length, the op-replay buffers — in `reconstruct`. The #142
    // "bound the aggregate allocation against the input" rule the tiler already
    // applies (see `MAX_BASES_PER_BYTE`). Deliberately loose (a real block decodes
    // far below this); it only fails a hostile length. In `decode_against`, `src`
    // excludes the shared reference, so bases/byte runs higher there, but the
    // constant keeps a wide margin over any real block.
    if total_bases
        > src
            .len()
            .saturating_mul(MAX_BASES_PER_BYTE)
            .saturating_add(MAX_BASES_PER_BYTE)
    {
        return Err(Error::Corrupt);
    }

    let lens_raw = get_stream(src, pos)?;
    let mut lp = 0usize;
    // `Vec::new`, not `with_capacity(n)`: `n` is an unvalidated header varint, so a
    // capacity hint would itself be the alloc bomb. The push loop is bounded by
    // `lens_raw` (a real stream) exhausting into `Corrupt`. Mirrors `tile_decode`
    // and `Reference::decode`.
    let mut lens = Vec::new();
    for _ in 0..n {
        let l = read_varint(&lens_raw, &mut lp).ok_or(Error::Corrupt)?;
        lens.push(u32::try_from(l).map_err(|_| Error::Corrupt)?);
    }
    let offs = compute_offs(&lens, total_bases)?;
    Ok((n, total_bases, lens, offs))
}

/// Reconstruct the reads from the edit streams, given the reference contigs
/// `refs`. Shared by [`decode`] (reference read from the block) and
/// [`decode_against`] (reference supplied by the caller). `pos` must point at the
/// group-counts stream.
#[allow(clippy::too_many_arguments)]
fn reconstruct(
    src: &[u8],
    pos: &mut usize,
    refs: &[&[u8]],
    n: usize,
    lens: &[u32],
    offs: &[usize],
    total_bases: usize,
) -> Result<Vec<u8>, Error> {
    let n_contigs = refs.len();

    // Groups + manifest + streams.
    let counts_raw = get_stream(src, pos)?;
    let mut cp = 0usize;
    let mut counts = Vec::with_capacity(n_contigs + 1);
    for _ in 0..n_contigs + 1 {
        counts.push(read_varint(&counts_raw, &mut cp).ok_or(Error::Corrupt)? as usize);
    }

    let read_ids = get_stream(src, pos)?;
    let kinds = get_stream(src, pos)?;
    let flips = get_stream(src, pos)?;
    let placements = get_stream(src, pos)?;
    // ops: range-coded with an order-3 op-history context (see `encode_ops`).
    // Decoded inline during reconstruction; the rolling context is rebuilt from the
    // decoded ops exactly as on encode.
    let ops_count = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let ops_blob_len = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let ops_end = pos.checked_add(ops_blob_len).ok_or(Error::Corrupt)?;
    let ops_blob = src.get(*pos..ops_end).ok_or(Error::Corrupt)?;
    *pos = ops_end;
    let mut op_dec = fqxv_range::Decoder::new(ops_blob);
    let mut op_models: [fqxv_range::SimpleModel<4>; OP_CTXS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut op_ctx = 0usize;
    let mut ops_seen = 0usize;
    let runs = get_stream(src, pos)?;
    // subs: range-coded, conditioned on the reference base (see `encode_subs`).
    // Decoded inline during reconstruction — the consensus base at each substitution
    // selects the model, exactly as on encode.
    let subs_count = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let subs_blob_len = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let subs_end = pos.checked_add(subs_blob_len).ok_or(Error::Corrupt)?;
    let subs_blob = src.get(*pos..subs_end).ok_or(Error::Corrupt)?;
    *pos = subs_end;
    let mut sub_dec = fqxv_range::Decoder::new(subs_blob);
    let mut sub_models: [fqxv_range::SimpleModel<SUB_SYMS>; SUB_SYMS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut subs_seen = 0usize;
    // ins_bases: range-coded, conditioned on the preceding read base. Decoded
    // inline during reconstruction with the same rolling context.
    let ins_count = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let ins_blob_len = read_varint(src, pos).ok_or(Error::Corrupt)? as usize;
    let ins_end = pos.checked_add(ins_blob_len).ok_or(Error::Corrupt)?;
    let ins_blob = src.get(*pos..ins_end).ok_or(Error::Corrupt)?;
    *pos = ins_end;
    let mut ins_dec = fqxv_range::Decoder::new(ins_blob);
    let mut ins_models: [fqxv_range::SimpleModel<SUB_SYMS>; SUB_SYMS] =
        std::array::from_fn(|_| fqxv_range::SimpleModel::new());
    let mut ins_seen = 0usize;
    let indel_lens = get_stream(src, pos)?;
    let literals = get_stream(src, pos)?;
    let exc_pos = get_stream(src, pos)?;
    let exc_bytes = get_stream(src, pos)?;

    let mut seq = vec![0u8; total_bases];
    // Cursors into each stream, advanced in the exact order `serialize` wrote them.
    let mut idp = 0usize; // read_ids
    let mut kp = 0usize; // kinds
    let mut fp = 0usize; // flips
    let mut pp = 0usize; // placements
    let mut rp = 0usize; // runs
    let mut dp = 0usize; // indel_lens
    let mut litp = 0usize; // literals

    let write_read = |seq: &mut [u8], r: usize, bytes: &[u8]| -> Result<(), Error> {
        let want = lens[r] as usize;
        if bytes.len() != want {
            return Err(Error::Corrupt);
        }
        seq[offs[r]..offs[r + 1]].copy_from_slice(bytes);
        Ok(())
    };

    for (ci, &count) in counts.iter().enumerate() {
        let is_orphan_group = ci == n_contigs;
        let mut prev_start: i64 = 0;
        for _ in 0..count {
            let r = read_varint(&read_ids, &mut idp).ok_or(Error::Corrupt)? as usize;
            if r >= n {
                return Err(Error::Corrupt);
            }
            let kind = *kinds.get(kp).ok_or(Error::Corrupt)?;
            kp += 1;
            if kind == 1 {
                // Literal: `lens[r]` base codes, corrected later by exceptions.
                let want = lens[r] as usize;
                let bytes: Vec<u8> = literals
                    .get(litp..litp + want)
                    .ok_or(Error::Corrupt)?
                    .iter()
                    .map(|&c| base_of_sym(c))
                    .collect();
                litp += want;
                write_read(&mut seq, r, &bytes)?;
                continue;
            }
            if is_orphan_group {
                return Err(Error::Corrupt); // orphan group is literal-only
            }
            let flip = *flips.get(fp).ok_or(Error::Corrupt)? != 0;
            fp += 1;
            let d = fqxv_bytes::unzigzag(read_varint(&placements, &mut pp).ok_or(Error::Corrupt)?);
            let start = (prev_start + d) as usize;
            prev_start += d;
            let cs = *refs.get(ci).ok_or(Error::Corrupt)?;
            let cs = cs.get(start..).ok_or(Error::Corrupt)?;

            // Rebuild this read's ops until it has produced `lens[r]` bases. Track
            // the reference position so each substitution is decoded with the same
            // reference-base context the encoder used.
            let want = lens[r] as usize;
            let mut produced = 0usize;
            let mut ref_pos = 0usize;
            let mut last: u8 = 4; // preceding read base, for the insertion context
            let mut read_ops: Vec<Op> = Vec::new();
            while produced < want {
                let code = op_models[op_ctx].decode(&mut op_dec) as u8;
                op_ctx = ((op_ctx << 2) | code as usize) & (OP_CTXS - 1);
                ops_seen += 1;
                match code {
                    0 => {
                        let m = read_varint(&runs, &mut rp).ok_or(Error::Corrupt)? as usize;
                        if m > 0 {
                            last = cs.get(ref_pos + m - 1).map_or(4, |&x| base_code(x));
                        }
                        produced += m;
                        ref_pos += m;
                        read_ops.push(Op::Match(m as u32));
                    }
                    1 => {
                        let rc = cs.get(ref_pos).map_or(4, |&x| base_code(x)) as usize;
                        let vc = sub_models[rc.min(SUB_SYMS - 1)].decode(&mut sub_dec);
                        subs_seen += 1;
                        produced += 1;
                        ref_pos += 1;
                        last = vc as u8;
                        read_ops.push(Op::Sub(base_of_sym(vc as u8)));
                    }
                    2 => {
                        let k = read_varint(&indel_lens, &mut dp).ok_or(Error::Corrupt)? as usize;
                        // A single insertion cannot legitimately produce more bases
                        // than the whole read (`want`), so reject a hostile `k`
                        // before it sizes the allocation below — `ins_dec` is a
                        // range decoder that never exhausts to stop the loop itself.
                        if k > want {
                            return Err(Error::Corrupt);
                        }
                        let mut bs: Vec<u8> = Vec::with_capacity(k);
                        for _ in 0..k {
                            let vc =
                                ins_models[(last as usize).min(SUB_SYMS - 1)].decode(&mut ins_dec);
                            last = vc as u8;
                            bs.push(base_of_sym(vc as u8));
                        }
                        ins_seen += k;
                        produced += k;
                        read_ops.push(Op::Ins(bs));
                    }
                    3 => {
                        let m = read_varint(&indel_lens, &mut dp).ok_or(Error::Corrupt)? as u32;
                        ref_pos += m as usize;
                        read_ops.push(Op::Del(m));
                    }
                    _ => return Err(Error::Corrupt),
                }
            }
            if produced != want {
                return Err(Error::Corrupt);
            }
            let oriented = apply(cs, &read_ops);
            if oriented.len() != want {
                return Err(Error::Corrupt);
            }
            let read = if flip {
                revcomp_acgt(&oriented)
            } else {
                oriented
            };
            write_read(&mut seq, r, &read)?;
        }
    }
    // The decoded op / substitution / insertion counts must match what was framed.
    if subs_seen != subs_count || ops_seen != ops_count || ins_seen != ins_count {
        return Err(Error::Corrupt);
    }

    // Apply non-ACGT exceptions last.
    let mut ep = 0usize;
    let mut gpos: u64 = 0;
    for &b in &exc_bytes {
        let delta = read_varint(&exc_pos, &mut ep).ok_or(Error::Corrupt)?;
        gpos += delta;
        let idx = gpos as usize;
        *seq.get_mut(idx).ok_or(Error::Corrupt)? = b;
        // The first delta is the absolute index; subsequent are gaps. Because
        // the first `prev_pos` was 0 and the first stored value is the index
        // itself, gpos already holds the absolute position.
    }

    Ok(seq)
}

/// Decompress a block produced by [`encode`], returning `(lens, seq)` in the
/// original read order.
pub fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>), Error> {
    let mut pos = 0usize;
    let (n, total_bases, lens, offs) = read_header(src, &mut pos, VERSION)?;

    // References, stored inline in a `VERSION` block.
    let n_contigs = read_varint(src, &mut pos).ok_or(Error::Corrupt)? as usize;
    let ref_lens_raw = get_stream(src, &mut pos)?;
    let ref_bytes = get_stream(src, &mut pos)?;
    let mut rlp = 0usize;
    // `n_contigs` is an untrusted header varint; do not pre-reserve (a crafted
    // count is an allocation bomb). The loop is bounded by `ref_lens_raw`
    // exhausting into `Corrupt`, exactly as `Reference::decode` documents.
    let mut refs: Vec<&[u8]> = Vec::new();
    let mut roff = 0usize;
    for _ in 0..n_contigs {
        let rl = read_varint(&ref_lens_raw, &mut rlp).ok_or(Error::Corrupt)? as usize;
        let end = roff.checked_add(rl).ok_or(Error::Corrupt)?;
        refs.push(ref_bytes.get(roff..end).ok_or(Error::Corrupt)?);
        roff = end;
    }

    let seq = reconstruct(src, &mut pos, &refs, n, &lens, &offs, total_bases)?;
    Ok((lens, seq))
}

/// Decompress a block produced by [`encode_against`], returning `(lens, seq)` in
/// the original read order. `reference` must be the exact same [`Reference`] the
/// block was coded against.
pub fn decode_against(reference: &Reference, src: &[u8]) -> Result<(Vec<u32>, Vec<u8>), Error> {
    let mut pos = 0usize;
    let (n, total_bases, lens, offs) = read_header(src, &mut pos, VERSION_SHARED)?;
    let refs: Vec<&[u8]> = reference.consensi.iter().map(Vec::as_slice).collect();
    let seq = reconstruct(src, &mut pos, &refs, n, &lens, &offs, total_bases)?;
    Ok((lens, seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random ACGT of `len` bases from `seed`.
    fn rand_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
        (0..len)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                b"ACGT"[(s & 3) as usize]
            })
            .collect()
    }

    fn roundtrip(reads: &[Vec<u8>], opts: &EncodeOpts) {
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        let enc = encode(&lens, &seq, opts).expect("encode");
        let (dl, ds) = decode(&enc).expect("decode");
        assert_eq!(dl, lens, "lengths");
        assert_eq!(ds, seq, "sequence must round-trip exactly");
    }

    /// Round-trip through the shared-reference path: build one reference over the
    /// reads, code them against it (no inline reference), and decode with the same
    /// reference. Must reconstruct the reads exactly.
    fn roundtrip_shared(reads: &[Vec<u8>], opts: &EncodeOpts) {
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        let reference = build_reference(&lens, &seq, opts).expect("build_reference");
        // The reference frame round-trips independently.
        let frame = reference.encode().expect("reference encode");
        let reference2 = Reference::decode(&frame).expect("reference decode");
        assert_eq!(
            reference2.consensi(),
            reference.consensi(),
            "reference frame"
        );

        let enc = encode_against(&reference, &lens, &seq, opts).expect("encode_against");
        let (dl, ds) = decode_against(&reference2, &enc).expect("decode_against");
        assert_eq!(dl, lens, "lengths");
        assert_eq!(
            ds, seq,
            "sequence must round-trip exactly against a shared reference"
        );
    }

    fn owned(reads: &[&[u8]]) -> Vec<Vec<u8>> {
        reads.iter().map(|r| r.to_vec()).collect()
    }

    #[test]
    fn decode_rejects_oversized_header_counts_without_aborting() {
        // Decode-path allocation guards (#142 class): a hostile read count must not
        // drive a `Vec::with_capacity` bomb, and a hostile `total_bases` must be
        // rejected before it sizes the `vec![0u8; total_bases]` output zero-fill.
        // Each returns `Corrupt` rather than aborting on the allocation. `decode`
        // and `decode_against` share `read_header`, so the `total_bases` guard
        // covers both.

        // (a) Oversized read count `n`, empty lengths stream, zero total_bases:
        // reaches the (now unreserved) length loop and errors on the exhausted
        // stream instead of pre-reserving `n` elements.
        let mut src = Vec::new();
        src.extend_from_slice(&MAGIC);
        src.push(VERSION);
        write_varint(&mut src, 1u64 << 62);
        write_varint(&mut src, 0);
        put_stream(&mut src, &[]).unwrap();
        assert!(matches!(decode(&src), Err(Error::Corrupt)));

        // (b) `total_bases` far larger than the coded block could hold: rejected up
        // front by the MAX_BASES_PER_BYTE guard.
        let mut src = Vec::new();
        src.extend_from_slice(&MAGIC);
        src.push(VERSION);
        write_varint(&mut src, 1);
        write_varint(&mut src, 1u64 << 50);
        assert!(matches!(decode(&src), Err(Error::Corrupt)));

        // (c) Oversized inline contig count: the reference-collection loop is also
        // unreserved and errors on the empty contig-length stream.
        let mut src = Vec::new();
        src.extend_from_slice(&MAGIC);
        src.push(VERSION);
        write_varint(&mut src, 0); // n reads
        write_varint(&mut src, 0); // total_bases
        put_stream(&mut src, &[]).unwrap(); // lens_raw
        write_varint(&mut src, 1u64 << 62); // absurd n_contigs
        put_stream(&mut src, &[]).unwrap(); // ref_lens_raw
        put_stream(&mut src, &[]).unwrap(); // ref_bytes
        assert!(matches!(decode(&src), Err(Error::Corrupt)));
    }

    #[test]
    fn empty_and_singletons() {
        let opts = EncodeOpts {
            stride: Some(1),
            ..EncodeOpts::default()
        };
        roundtrip(&[], &opts);
        roundtrip(&owned(&[b"ACGTACGTACGT"]), &opts);
        roundtrip(&owned(&[b"", b"A", b"ACGT"]), &opts);
        // The shared path must also degrade gracefully on trivial input (a
        // reference with no contigs → every read is an orphan literal).
        roundtrip_shared(&[], &opts);
        roundtrip_shared(&owned(&[b"ACGTACGTACGT"]), &opts);
        roundtrip_shared(&owned(&[b"", b"A", b"ACGT"]), &opts);
    }

    #[test]
    fn non_acgt_and_lowercase_survive() {
        let opts = EncodeOpts {
            stride: Some(1),
            ..EncodeOpts::default()
        };
        // These form no contig, so they exercise the literal + exception path.
        roundtrip(&owned(&[b"ACGTNNNNacgtACGTRYKM"]), &opts);
        roundtrip(&owned(&[b"NNNN", b"acgtacgt", b"ACGTNacgtN"]), &opts);
        roundtrip_shared(&owned(&[b"ACGTNNNNacgtACGTRYKM"]), &opts);
        roundtrip_shared(&owned(&[b"NNNN", b"acgtacgt", b"ACGTNacgtN"]), &opts);
    }

    /// Reads that tile a genome form a contig and are coded against its
    /// consensus — the edit-script/reference path, plus an unrelated orphan and
    /// embedded non-ACGT bases that must survive it.
    fn tiled_reads(seed: u64) -> Vec<Vec<u8>> {
        let genome = rand_seq(30_000, seed);
        let mut reads: Vec<Vec<u8>> = (0..10)
            .map(|i| {
                let mut s = genome[i * 2000..i * 2000 + 6000].to_vec();
                // ~0.7% substitutions: sequencing error the aligner must code.
                for (j, b) in s.iter_mut().enumerate() {
                    if (i * 7 + j) % 143 == 0 {
                        *b = b"ACGT"[((*b as usize) + 1) & 3];
                    }
                }
                // Every other read comes off the sequencer reverse-complemented.
                if i % 2 == 1 { revcomp_acgt(&s) } else { s }
            })
            .collect();
        // A non-ACGT base and a lowercase base inside a placed read.
        reads[0][100] = b'N';
        reads[2][250] = b'n';
        // An unrelated read: lands on no reference, coded as a literal orphan.
        reads.push(rand_seq(4000, seed ^ 0xdead));
        reads
    }

    #[test]
    fn tiled_contig_roundtrips_through_edit_scripts() {
        let opts = EncodeOpts {
            sketch: Sketch::ont(),
            stride: Some(1),
            ..EncodeOpts::default()
        };
        for seed in [1u64, 88, 12345] {
            roundtrip(&tiled_reads(seed), &opts);
        }
    }

    #[test]
    fn tiled_contig_roundtrips_against_shared_reference() {
        let opts = EncodeOpts {
            sketch: Sketch::ont(),
            stride: Some(1),
            ..EncodeOpts::default()
        };
        for seed in [1u64, 88, 12345] {
            roundtrip_shared(&tiled_reads(seed), &opts);
        }
    }

    /// A reference built over the WHOLE set codes every one of several disjoint
    /// blocks — the container's use case (issue #168): one shared reference,
    /// coded per block, each block decoded against the same reference.
    #[test]
    fn shared_reference_codes_disjoint_blocks() {
        let opts = EncodeOpts {
            sketch: Sketch::ont(),
            stride: Some(1),
            ..EncodeOpts::default()
        };
        let all = tiled_reads(88);
        let lens: Vec<u32> = all.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = all.iter().flat_map(|r| r.iter().copied()).collect();
        let reference = build_reference(&lens, &seq, &opts).expect("build_reference");

        // Split the reads into three contiguous blocks; each is coded against the
        // whole-file reference and must decode to its own reads exactly.
        let bounds = [0usize, 4, 8, all.len()];
        for w in bounds.windows(2) {
            let (s, e) = (w[0], w[1]);
            let blk_lens: Vec<u32> = lens[s..e].to_vec();
            let blk_seq: Vec<u8> = all[s..e].iter().flat_map(|r| r.iter().copied()).collect();
            let enc =
                encode_against(&reference, &blk_lens, &blk_seq, &opts).expect("encode_against");
            let (dl, ds) = decode_against(&reference, &enc).expect("decode_against");
            assert_eq!(dl, blk_lens, "block [{s},{e}) lengths");
            assert_eq!(ds, blk_seq, "block [{s},{e}) sequence");
        }
    }

    proptest::proptest! {
        // A full assemble + encode + decode per case; cap the count so this
        // stays a fast regression net rather than a slow test.
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(64))]

        /// Any pile of arbitrary-byte reads round-trips exactly — the exception
        /// list makes losslessness independent of the alphabet, and short/empty
        /// reads exercise the boundaries. Both the self-contained and the
        /// shared-reference paths must reconstruct exactly.
        #[test]
        fn arbitrary_reads_round_trip(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::prelude::any::<u8>(), 0..40),
                0..24,
            )
        ) {
            let opts = EncodeOpts { stride: Some(1), ..EncodeOpts::default() };
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();

            let enc = encode(&lens, &seq, &opts).expect("encode");
            let (dl, ds) = decode(&enc).expect("decode");
            proptest::prop_assert_eq!(&dl, &lens);
            proptest::prop_assert_eq!(&ds, &seq);

            let reference = build_reference(&lens, &seq, &opts).expect("build_reference");
            let enc2 = encode_against(&reference, &lens, &seq, &opts).expect("encode_against");
            let (dl2, ds2) = decode_against(&reference, &enc2).expect("decode_against");
            proptest::prop_assert_eq!(dl2, lens);
            proptest::prop_assert_eq!(ds2, seq);
        }
    }

    #[test]
    fn output_is_thread_count_invariant() {
        let reads = tiled_reads(88);
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        let opts = EncodeOpts {
            sketch: Sketch::ont(),
            stride: Some(1),
            ..EncodeOpts::default()
        };
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| encode(&lens, &seq, &opts).expect("encode"))
        };
        assert_eq!(
            run(1),
            run(8),
            "encode output must be byte-identical regardless of thread count"
        );
    }

    /// The shared-reference path is likewise thread-count invariant: both the
    /// reference assembly and the per-block coding must be pure functions.
    #[test]
    fn shared_reference_is_thread_count_invariant() {
        let reads = tiled_reads(88);
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        let opts = EncodeOpts {
            sketch: Sketch::ont(),
            stride: Some(1),
            ..EncodeOpts::default()
        };
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let reference = build_reference(&lens, &seq, &opts).expect("build_reference");
                    let frame = reference.encode().expect("reference encode");
                    let block =
                        encode_against(&reference, &lens, &seq, &opts).expect("encode_against");
                    (frame, block)
                })
        };
        assert_eq!(
            run(1),
            run(8),
            "shared-reference output must be byte-identical regardless of thread count"
        );
    }
}
