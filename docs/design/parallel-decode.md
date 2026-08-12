# Parallel-First Decode

Design + measurement round for within-block decode parallelism (issues #273 and
#269, PR #275's follow-up). **Nothing here ships a format change**: this
document specifies the candidate design, and the numbers come from encoder-side
diagnostics (`FQXV_DIAG_QCHUNK`, `FQXV_DIAG_TILECHUNK`) that re-code real blocks
under the candidate layouts without touching the emitted archive. The go/no-go
bars in this document were written **before** the measurements landed.

> **Status: shipped.** The go/no-go issue (#278) exercised the GO on both
> long-read platforms, and the chunked layout below ships as quality modes
> **5** (`MODE_SEQ_BINMIX_CHUNKED`) and **6** (its per-block-quantizer twin,
> mode 5 : 6 :: 3 : 4) — warm-clone at K = 8 by platform default on
> Nanopore/PacBio, `--max` pins serial, reset (variant 0) reserved unimplemented,
> and warm-frozen carried for its O(1)-memory decode. See `fqxv-fqzcomp`'s
> mode docs for the final wire layout; this document remains the design
> record and measurement archive.

## Problem

`.fqxv` decode parallelism is per block. Long-read archives have few blocks
(a 576 MB MinION file is ~3 blocks at the 256 MiB budget), so decode saturated
at ~3 threads. #275 shrank the Nanopore block budget to 64 MiB
(`LONGREAD_BLOCK_SEQ_BYTES`), buying 3.2× full-decode at 8 threads for +1.08%
archive size — but block-shrink has hit its knee: halving again measured +2.7%
for 5.3× at 16 threads, and the per-block ratio cost compounds while
within-block decode stays serial. The serial core is the quality stream: on ONT,
quality is **92% of single-thread decode time and 81% of archive bytes** (#269),
coded by an adaptive context mixer that threads model state across every read in
the block, conditioned on the already-decoded bases.

#273's verdict on within-block parallelism was "dead end without a format
change." This round asks what the format change would actually cost, measured on
real data, so a future mode can be a decision rather than a bet.

### What makes quality chunkable at all

Two structural facts (both now pinned by tests in `fqxv-fqzcomp`):

1. **Every per-read context feature resets at read boundaries**
   (`binmix::encode_reads`, `encode_payload_into`): the recent-quality window,
   previous base, homopolymer run, position, delta counter. Coding reads
   `[a, b)` depends only on the bases of reads `[a, b)` plus cross-read state.
2. **The cross-read state is exactly the model** (tier tables + mixer weights +
   range coder), which is structurally trivial to snapshot — `BinMixer` is three
   `Vec<u16>` tiers plus per-bit weights, now `Clone`. Memory, not structure, is
   the design constraint (see the memory budget below).

## Measured decode decomposition

Full decode vs `--fasta` (no quality stream) at 1 thread, 3 reps, median, output
piped to `wc -c`; archives built by this branch's binary at the default level
(`bench/scripts/qchunk_diag.sh`, results `qchunk_decode_share.tsv`).

| dataset | platform | t_full @1t | t_fasta @1t | quality share | t_full @8t | t_full @16t |
|---|---|---|---|---|---|---|
| DRR205413 | ONT MinION | 64.1 s (9.4 MB/s) | 5.0 s | **92.3%** (re-confirms #269's 92%) | 17.7 s (3.62×) | 17.6 s (3.64× — saturated) |
| ecoli_hifi | PacBio HiFi | 284.2 s (10.9 MB/s) | 6.2 s | **97.8% — first HiFi decomposition** | 66.3 s (4.29×) | 66.0 s (4.31× — saturated) |
| novaseq_4m | NovaSeq (control) | 46.3 s (29.5 MB/s) | 37.8 s | 18.4% | 14.4 s (3.21×) | 14.3 s (3.23×) |

The ONT 8→16-thread step is flat: five blocks (#275's 64 MiB budget) exhaust
block-level parallelism, exactly the ceiling this design targets.

## Measured chunking cost

`FQXV_DIAG_QCHUNK`: every quality block re-coded under
{reset, warm-clone, warm-frozen} × K ∈ {4, 8, 16, 32} × warmup ∈ {total/K,
8 MiB}, charged with per-chunk range-coder flush and an estimated chunk header
(variant byte, segment count, warmup boundary + length, per-chunk byte lengths
as varints). Baseline = the payload bytes the archive actually carries (the
probe runs at the winner-decision point, under the winning mode and quantizer).

Variants:

- **reset** — every chunk starts from a fresh model; K parallel units; each
  decoding chunk owns a **mutable** model.
- **warm-clone** — segment 0 (the warmup prefix) is coded from scratch and
  decoded serially; a snapshot is taken at its end; the remaining reads split
  into K−1 chunks, each starting from a **clone** of the snapshot and adapting
  within the chunk.
- **warm-frozen** — as warm-clone, but the chunks code against the **shared
  read-only** snapshot (`ADAPT = false`): one model for any number of chunks,
  the O(1)-memory variant.

Chunk boundaries are fixed-K, whole-read, equal-cumulative-bases — a pure
function of `(lens, K)` (see Determinism).

Per-dataset totals below (payload-relative delta = Σ(total−baseline)/Σbaseline;
archive-relative divides by the archive size instead; full per-block cells in
`qchunk_diag.tsv`).

### ONT (DRR205413, mode 3, five 64 MiB blocks, k = 65)

Quality payload 164.4 MB of the 203.1 MB archive (**80.9% of archive bytes**).
Cell = payload-relative delta, with the archive-relative delta in parentheses.

| variant | warmup | K=4 | K=8 | K=16 | K=32 |
|---|---|---|---|---|---|
| reset | — | +1.17% (0.95%) | +2.09% (1.69%) | +3.30% (2.67%) | +4.88% (3.95%) |
| warm-clone | total/K | +0.15% (0.12%) | +0.48% (0.39%) | +1.02% (0.83%) | +1.81% (1.47%) |
| warm-clone | 8 MiB | +0.27% (0.22%) | **+0.44% (0.36%)** | +0.55% (0.45%) | +0.62% (0.50%) |
| warm-frozen | total/K | +0.47% (0.38%) | +0.90% (0.73%) | +1.94% (1.57%) | +2.76% (2.24%) |
| warm-frozen | 8 MiB | +0.76% (0.61%) | +0.76% (0.61%) | **+0.76% (0.61%)** | +0.76% (0.61%) |

Two shapes worth naming:

- **warm-frozen with a fixed warmup is K-invariant** (+0.756…0.757% payload at
  every K): once adaptation stops at the snapshot, the loss is "no adaptation
  after the warmup," which does not depend on how many chunks share the frozen
  model — only the per-chunk flush/header grows, and that is noise. Frozen buys
  **unbounded fan-out at a fixed price**.
- **warm-clone degrades gently in K** and beats frozen at every K measured —
  within-chunk adaptation recovers most of what freezing gives up — but its
  memory scales with the number of concurrently-decoding chunks.

### HiFi (ecoli_hifi, mode 4, six 256 MiB blocks, k = 93)

Quality payload 635.4 MB of the 649.5 MB archive (**97.8% of archive bytes** —
this dataset is, byte-wise, a quality stream with a rounding error).

| variant | warmup | K=4 | K=8 | K=16 | K=32 |
|---|---|---|---|---|---|
| reset | — | +0.59% (0.58%) | +1.12% (1.10%) | +1.89% (1.85%) | +2.97% (2.91%) |
| warm-clone | total/K | +0.05% (0.05%) | **+0.19% (0.19%)** | +0.45% (0.44%) | +0.88% (0.86%) |
| warm-clone | 8 MiB | +0.23% (0.22%) | +0.46% (0.45%) | +0.68% (0.66%) | +0.85% (0.84%) |
| warm-frozen | total/K | +5.49% (5.37%) | +6.26% (6.12%) | +6.63% (6.49%) | +9.38% (9.18%) |
| warm-frozen | 8 MiB | +9.59% (9.38%) | +9.59% (9.38%) | +9.59% (9.38%) | +9.59% (9.38%) |

The platform contrast is the finding of this table. On ONT, freezing the model
after a warmup costs 0.76% — ONT quality statistics are near-stationary within
a block. On HiFi the same freeze costs **9.6%**: CCS quality drifts (per-ZMW
pass count, read-length mixes), so continued adaptation is where mode 4's ratio
lives, and any parallel variant must let chunks keep adapting (warm-clone) or
pay double digits. warm-frozen's K-invariance shows again (9.587% at every K) —
the property is real, the price is platform-dependent. Chunking itself is
*cheaper* on HiFi than ONT (reset @K=4: 0.59% vs 1.17%): with 93 symbols over
a slowly-drifting distribution, a cold model re-learns quickly relative to the
268 MB it codes.

### NovaSeq control (novaseq_4m, MODE_POS, four blocks, k = 4)

Quality payload 30.7 MB of the 152.6 MB archive (20.1% — short-read archives
are sequence-dominated, and short-read files already decode block-parallel, so
this is a control, not a target).

| variant | warmup | K=4 | K=8 | K=16 | K=32 |
|---|---|---|---|---|---|
| reset | — | +0.031% | **+0.068%** | +0.134% | +0.252% |
| warm-clone | total/K | +0.001% | +0.004% | +0.010% | +0.025% |
| warm-clone | 8 MiB | +0.004% | +0.007% | +0.011% | +0.015% |
| warm-frozen | total/K | +0.173% | +0.240% | +0.216% | +0.259% |
| warm-frozen | 8 MiB | +0.292% | +0.292% | +0.292% | +0.294% |

(Payload-relative; archive-relative is ~5× smaller.) With four symbols the
`SimpleModel` table warms in almost nothing, so even reset is far under its
0.5% bar — chunking MODE_POS is free whenever anyone wants it, in any variant,
and reset would do (no snapshot machinery: the models are ~3 MiB at this size
class).

### ONT tiler confinement (`FQXV_DIAG_TILECHUNK`)

The ONT **sequence** stream (`SEQ_METHOD_TILE`) is a read-chain: each read's
tiles reference earlier raw reads, so `tile_decode` is serial. Chunk-parallel
decode needs the encoder to confine references to the read's own chunk
(`EncodeOpts::diag_confine_chunks`). Measured cost of that confinement:

| K | confinement cost (of the tile stream) |
|---|---|
| 4 | **+21.8%** (per block: +18.9…+30.9%) |
| 8 | **+41.1%** (per block: +37.0…+49.3%) |

The tile stream is 38.6 MB of the 203.1 MB ONT archive (~19%), so confinement
would cost **+4.1% / +7.8% of the archive** at K = 4 / 8 — an order of
magnitude past any acceptable bar. A read's cheapest reference is usually a
near neighbour, but confinement forces every chunk to restart its reference
pool from zero coverage, so early-chunk reads fall back to literals. This
lever is dead as measured; the ONT sequence stream stays serial.

## Wire sketch: `MODE_SEQ_BINMIX_CHUNKED = 5`

The evolution policy blesses a **new fqzcomp mode byte** for an alternative
encoding of the same stream — exactly how modes 3 (binmix) and 4 (binmix + Q)
shipped. No `required_features` bit, no footer change, no container change: the
quality stream stays one length-prefixed blob per block; only its internal
layout gains segments.

```text
version(=2) | binning | mode=5 | k | syms[k]
[1] variant        (0 reset, 1 warm-clone, 2 warm-frozen)
[v] K              (segment count)
[v] warm_bases     (warmup boundary in bases; 0 for reset)
[v] warm_len       (warmup segment's coded length; absent for reset)
[v] chunk_len × (K-1 | K)   (each parallel chunk's coded length)
lens               (the existing read-length array)
payload            (segments concatenated: [warmup] chunk_0 .. chunk_last)
```

The mode-5 header sits where mode 4's qtable sits (after `syms`, before
`lens`). **Chunk boundaries are never transmitted**: both sides recompute them
from `(lens, K)` (and `warm_bases`, itself derived from the warmup policy) —
the same trick `GlobalReference::encode_blocked` uses. Segment byte lengths ARE
transmitted so the decoder can slice the payload and hand each chunk to a
worker without parsing it.

### Decode pipeline

A sequence-conditioned chunk needs its reads' bases, so the block DAG becomes:

```text
names ──────────────────────────────┐
seq ─┬─ (bases of chunk c ready) ─┬─┴─ assemble block
     │   qual chunk 0 (warmup)    │
     │   qual chunk 1..K-1  (par) │
```

With the sequence stream still serial (ONT tile chain; lroverlap mostly
parallelizable, see Secondary levers), quality chunks start as soon as the
bases exist; the block's critical path drops from `T_seq + T_qual` toward
`T_seq + T_qual·w + T_qual·(1−w)/(K−1)` (see the speedup model).

`needs_sequence` must learn mode 5 — **already done**: as of this branch it
fail-safes to `true` for every non-`MODE_POS` mode under the recognized
version, and `decode_seq` refuses unknown modes with
`Error::UnsupportedMode` ("this archive needs a newer fqxv") instead of a
generic corruption error. Every reader from this branch on therefore refuses a
future mode-5 archive cleanly, with the sequence already in hand.

## Determinism

Byte-determinism (thread-count invariance) is a hard repo invariant. The
chunked design preserves it because nothing about the layout is thread-derived:

- Boundaries are a pure function of `(lens, K)`; K is a pure function of the
  effort level/platform, never of the worker count (pinned by
  `chunk_boundaries_are_a_pure_function_of_lens`).
- Each chunk's bytes are a pure function of (chunk reads, snapshot state), and
  the snapshot is a pure function of the warmup prefix. Encode may code chunks
  in any order or in parallel and concatenate by index — same bytes.
- The integer coders are already platform-deterministic (fixed-point binmix,
  Subbotin range coder); a cloned/frozen snapshot replays identically.
- Decoded-content digests are unaffected: chunking changes how bytes are
  *coded*, never the qualities they decode to.

## Memory budget

Model sizes (this branch, from the diag `k` per dataset):

`BinMixer` ≈ (2^16 + 2^18 + 2^20) tier slots × nnodes × 2 B with
`nnodes = 2^ceil(log2 k)`; the MODE_POS table is 2^18 × (2·NM + 4) B.

| dataset | k | nnodes | model size |
|---|---|---|---|
| ONT DRR205413 | 65 | 128 | **336 MiB** (full-alphabet tree depth) |
| HiFi ecoli_hifi | 93 | 128 | **336 MiB** |
| NovaSeq (MODE_POS, size class NM=4) | 4 | — | **3 MiB** (2^18 × (2·4+4) B) |

Decode-side memory as f(variant, concurrency), per **in-flight block**:

- **reset**: every decoding chunk owns a mutable model → `K × model` while the
  block decodes. At ONT's model size this is the variant that does NOT fit: a
  16-chunk 64 MiB block would want multiple GiB of model alone. Reset is only
  viable where the model is small (MODE_POS size classes, small-k long reads).
- **warm-clone**: warmup decodes with one model; each in-flight chunk clones
  the snapshot → `(active chunks + 1) × model`. Same asymptotics as reset with
  a smaller constant only if the scheduler bounds active chunks.
- **warm-frozen**: one snapshot, shared read-only by every chunk → `O(1)`
  model memory regardless of K. **The only variant whose memory does not scale
  with parallelism.** Its ratio cost vs warm-clone is the price of not
  adapting within chunks — measured above.

## Secondary levers, ranked

1. **lroverlap phases A/D (HiFi seq)** — the consensus codec's cross-read
   entropy state is ~928 B; its 10 rANS streams and per-read edit application
   are independent and could fan out **today with no format change**. Paper
   estimate only this round — and the measured HiFi decode share (quality
   97.8%, everything else 2.2%) says seq fan-out cannot move HiFi wall time
   until quality chunking ships.
2. **Tile confinement (ONT seq)** — measured here (`FQXV_DIAG_TILECHUNK`) at
   +21.8%/+41.1% of a stream that is only ~19% of the ONT archive and ~8% of
   its decode time: **dead as a lever**. Only worth revisiting if quality
   chunking ships, the sequence stream becomes the new serial floor, and a
   cheaper confinement (overlapping chunk prefixes, shared early reads) is
   designed.
3. **Names fan** — ~2–3% of decode; free but low payoff.
4. **Op-chunk reset (lroverlap edit streams)** — smallest models, smallest
   payoff; revisit only if HiFi seq measures > ~20% of decode.

## Out of scope

- **fqxv-seq (order-k)**: 42–579 MB of model state; reset-only chunking would
  be the whole design. Not this round.
- **Footer/container changes**: the quality stream chunks internally; the
  container's three-streams-per-block layout is untouched.
- **Shipping mode 5**: this round produces the numbers and the spec sketch;
  the go/no-go issue decides implementation.

## Go/no-go bars (written before the numbers landed)

- **GO (ONT quality chunking)** iff, at K=8, some variant costs **≤ 1.0% of
  archive bytes** (≤ 1.6% at K=16), AND its decode-memory overhead is ≤ 2× the
  current peak, AND the projected full-decode speedup at 16 threads is ≥ 2×
  over the #275 baseline.
- **GO (HiFi)**: same bars, plus quality must measure **≥ 50% of HiFi decode**
  (below that, chunking quality cannot move the total enough to justify a
  mode).
- **MODE_POS extension**: only if reset @K=8 costs ≤ 0.5% of quality-stream
  bytes — short-read archives already decode block-parallel, so the bar is
  deliberately higher.
- **Warm variants must beat reset** by enough to justify carrying snapshot
  machinery (clone) or the frozen-adaptation loss (frozen); otherwise ship the
  simpler variant that passes.

### Verdicts (measured)

- **ONT: GO.** warm-clone 8 MiB @K=8 = +0.36% of archive (bar: ≤1.0%);
  @K=16 = +0.45% (bar: ≤1.6%); warm-frozen O(1) memory satisfies the memory
  bar outright; projected 4.0–4.4× over the #275 baseline @16t (bar: ≥2×).
  Recommended: **warm-clone (8 MiB warmup) at K=8**, or **warm-frozen** where
  unbounded fan-out or strict memory matters (+0.61%, K-invariant).
- **HiFi: GO** — and it is the bigger prize: quality is 97.8% of decode (bar:
  ≥50%). warm-clone total/K @K=8 = +0.19% of archive (bar: ≤1.0%); projected
  3.7× over the baseline @16t. **warm-frozen is disqualified on HiFi**
  (+5.5–9.6%): the mode-5 header must carry the variant byte precisely so the
  encoder can pick per platform. Memory bar holds only with bounded
  concurrently-active chunks (each clone is 336 MiB at k=93) — the decoder
  must schedule chunks, not spawn all K.
- **Warm beats reset everywhere measured** (ONT K=8: 0.44% vs 2.09%; HiFi K=8:
  0.19% vs 1.12%), so the snapshot machinery earns its place; reset's only
  niche would be MODE_POS (see below).
- **MODE_POS: PASS but low-stakes.** reset @K=8 = 0.068% of quality-stream
  bytes (bar: ≤0.5%); any variant qualifies. Not worth a mode byte on its own —
  ride along only if mode 5 ships for long reads.

## Speedup model

Let `w` = warmup fraction of quality bytes (reset: `w = 0`), `T_seq`/`T_qual`
the per-block serial stream times, and K the segment count:

```text
T_block(chunked) ≈ T_names ∥ (T_seq + w·T_qual + (1−w)·T_qual/(K−1))   (warm)
T_block(reset)   ≈ T_names ∥ (T_seq + T_qual/K)
```

Combined with the measured decode shares and per-block times (biggest block's
critical path CP vs total-work bound T1/threads; wall ≈ max of the two):

**ONT** (T1 = 64.1 s, quality 59.2 s; big block: t_qual ≈ 13.2 s, t_other ≈
1.1 s; w = 8.2 M/67.1 M = 0.122 for the 8 MiB warmup):

- warm-clone 8 MiB, K=8: CP ≈ 1.11 + 0.122·13.19 + 0.878·13.19/7 ≈ **4.4 s** →
  wall@16t ≈ max(64.1/16, 4.4) = 4.4 s ≈ **14.7× vs 1 thread, 4.0× over the
  #275 baseline** (17.6 s), for +0.36% archive.
- warm-frozen 8 MiB, K=16: CP ≈ 3.5 s → wall@16t ≈ 4.0 s (work-bound) ≈
  **16× vs 1 thread, 4.4× over baseline**, for +0.61% archive — and the same
  archive serves any thread count, since frozen's cost is K-invariant.

**HiFi** (T1 = 284.2 s, quality 278.1 s; big block: t_qual ≈ 48.2 s, t_other ≈
1.1 s; perk warmup w = 1/K gives CP ≈ t_other + 2·t_qual/K):

- warm-clone total/K, K=8: CP ≈ 1.1 + 12.1 = 13.2 s → wall@16t ≈
  max(284.2/16, 13.2) = 17.8 s ≈ **16× vs 1 thread, 3.7× over the baseline**
  (66.0 s), for **+0.19% archive**. Higher K buys nothing at 16 threads (the
  total-work bound dominates) but keeps scaling at 32+.

Both projections assume enough scheduler freedom to co-run chunks of different
blocks; they are upper bounds on wall, not throughput promises.

## Migration & compatibility

- Mode 5 is a **mode byte**, not a feature bit: old readers refuse it with
  `UnsupportedMode` ("needs a newer fqxv") — the error path shipped in this
  branch, ahead of any writer. Readers older than this branch fail with a
  generic malformed-stream error; that window is why the error fix ships now.
- The encoder would gate mode 5 exactly like mode 4: code chunked and serial,
  keep the smaller (never-worse by construction), with the platform/effort
  gates deciding where the trial runs.
- `fqxv-python` needs nothing: it wraps the container decode, which sees the
  same `(lens, quals)` from `decode_seq`.
