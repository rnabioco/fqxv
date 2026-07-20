# Long-read (ONT / PacBio) support

Status: **design note + measured results.** fqxv is structurally ready for long
reads (see "What already works" below); this note scopes the two ratio levers —
a cross-read **overlap sequence codec** and **long-read quality tuning** —
against the field (CoLoRd, NanoSpring, ENANO/RENANO).

Where each lever stands:

| Lever | Status |
| --- | --- |
| Quality binning (`--quality-bin ont` / `hifi`) | **shipped** — usable from the CLI |
| Quality base-context | not started; measured headroom is ~nil (see Lever 1) |
| Overlap sequence codec (`fqxv-lroverlap`) | **shipped** — auto-selected for long-read blocks, kept only when it beats the order-k model; codes the HiFi sequence stream ~6× smaller through a real archive. See [Wiring](#wiring-and-the-per-block-coverage-cap). |

The rest of this note is the analysis that set those priorities; it is written
against the pre-`lroverlap` baseline, so the "fqxv seq" rows below are the
within-read order-k model — now the *fallback*, used only when the overlap codec
does not win a block.

## Why long reads are different

Three properties break the assumptions a short-read codec is tuned for:

1. **Ragged lengths.** ONT reads average ~10–20 kb and range from ~200 bp to
   >100 kb in one file; PacBio HiFi ~10–25 kb. There is no "read length" — every
   per-read overhead and every fixed positional model must be length-agnostic.
2. **A different quality model.** ONT quality is comparatively flat per read and
   correlates with *local sequence context* (homopolymers, specific k-mers), not
   with cycle position the way Illumina does. HiFi quality is a narrow band near
   the top of the scale. Both span the full Phred/Sanger range (Q0–Q93).
3. **High, base-level error.** ONT ~5–10 % (older basecallers worse), PacBio CLR
   similar; HiFi <1 %. Cross-read redundancy is still large at typical coverage,
   but overlaps must be found through a noisy channel — exact-anchor tricks that
   work on near-identical short reads do not survive a 10 % error rate.

## What already works (no change needed)

The container and codecs were hardened for long reads already:

- Read lengths are `u32` end to end (`container/*`, `fqxv-seq`, `fqxv-fqzcomp`)
  — a single read can be 4.29 Gbp; no truncation, no per-read overflow.
- Blocks are cut by a **256 MiB raw-sequence byte budget**
  (`MAX_BLOCK_SEQ_BYTES`), not just a read count, so a file of 14 kb reads does
  not collapse into one giant row group — parallelism and random access survive.
- The quality alphabet cap is the full Sanger range (`QMAX = 94`); the model is
  sized to the alphabet actually present, so narrow HiFi Q pays nothing for
  unused levels.
- The short-read **reorder** codec auto-disables on long reads
  (`is_long_read`, mean length > 500 bp) and falls back to the deep-context
  sequence path — its single-anchor minimizer clustering is useless at ONT
  error rates and is correctly skipped.

So the open question is **ratio**, not correctness: on the sequence and quality
streams, how do we close the gap to CoLoRd?

## Where the bytes are (measured)

Measured per-stream, lossless, vs CoLoRd `-q org` (its lossless quality mode).
Two datasets, deliberately spanning the regimes.

**How the split is measured.** fqxv has real separate streams, so its columns come
straight from `fqxv info`. CoLoRd does not, so its split is taken by difference:
`-q org` is the whole lossless archive, `-q none` discards quality (Q1 for every
base), and the quality column is `org - none`. That is **additive by
construction** — non-quality + quality equals the total, to the byte — which is
the point: an earlier version of this table took the total from a *different*
CoLoRd run than the columns beside it, and the rows did not add up (its quality
stream alone exceeded its whole archive).

Note the caveat this method carries: `-q none` still contains **names and
container overhead**, so CoLoRd's non-quality column is an *upper bound* on its
sequence stream, while fqxv's `seq` is the real thing.

For a seq-vs-seq comparison use **M1b's direct measurement of CoLoRd's sequence
stream: 0.0676 bits/base** (13.1M on `ecoli_hifi`), which excludes names — that is
the number the 9.7× gap below is computed from, and it is consistent with the
0.069 upper bound here (the difference is names, ~0.0014 bits/base).

The three sources are what went wrong before: the old table took `seq` from M1b,
`qual` from `org - none`, and `total` from the benchmark harness — three separate
runs, so the rows could not add up and the HiFi quality column exceeded its own
archive. Each column above now comes from the same pair of runs.

**`ecoli_ont`** (DRR205413, 287M bases, mean Q≈11.5 — noisy older-basecaller):

| tool | total | non-quality (seq+names) | qual | non-quality bits/base |
| --- | --- | --- | --- | --- |
| CoLoRd `-q org` | 197.9M | **31.4M** | 166.5M | 0.88 |
| fqxv (binmix qual) | 222.6M | **58.8M** (seq only) | 163.7M | 1.64 |

**`ecoli_hifi`** (SRR11434954 subset, 1.55G bases, mean Q≈27, ~300× — narrow
high-Q, low error):

| tool | total | non-quality (seq+names) | qual | non-quality bits/base |
| --- | --- | --- | --- | --- |
| CoLoRd `-q org` | 697.7M | **13.4M** | 684.3M | 0.069 |
| fqxv (binmix qual) | 768.3M | **126.3M** (seq only) | 641.8M | 0.653 |

Two facts, confirmed on **both** platforms:

1. **Quality now leads on both platforms.** fqxv's binary-decomposition
   context-mixing quality coder codes the HiFi quality stream to **641.8M vs
   CoLoRd's 684.3M** (~6% smaller) and ONT to **163.7M vs 166.5M** (~2% smaller).
   So quality is now a *credit* against CoLoRd on both sets — enough to carry the
   HiFi lossless total ahead of CoLoRd, though on ONT the far larger sequence
   deficit below still dominates. Lever 1 has flipped in fqxv's favor.
2. **The entire lossless gap to CoLoRd is the sequence stream, and it widens
   with coverage/fidelity.** ONT: 1.87× (0.87 vs 1.64 bits/base). HiFi:
   **9.7×** (0.0676 vs 0.653 bits/base). At ~300× HiFi the same genome is read
   hundreds of times; CoLoRd's overlap assembly encodes it once + diffs
   (0.068 bits/base — its published high-coverage regime) while fqxv's
   within-read model re-encodes every copy. **Lever 2 is the priority lever**,
   decisively, and most of all on HiFi.

### Lossy regime — the gap explodes

Measured vs `colord-lossy` (CoLoRd's default lossy quality). fqxv lossy rows run
at *default* sequence level (seq 65.3M ONT / 247.1M HiFi); `-l9` would cut seq to
58.8M / 126.3M but not change the verdict.

| dataset | tool | total | ratio | seq | qual | Δqual mae |
| --- | --- | --- | --- | --- | --- | --- |
| ONT | colord-lossy | **73.2M** | 7.87 | — | — | — |
| ONT | fqxv-binont | 114.6M | 5.03 | 65.3M (57%) | 49.2M (43%) | 3.35 |
| HiFi | colord-lossy | **45.3M** | 65.3 | ~13M | ~32M | — |
| HiFi | fqxv-binhifi | 397.2M | 7.44 | 247.1M (62%) | 149.9M (38%) | 14.33 |

The bin tables work — ONT bins cut fqxv's quality stream 163.7M → 49.2M (3.4×).
But binning removes the stream fqxv is *good* at and leaves the stream it is bad
at: **DNA becomes 62–88 % of the lossy archive**, so the sequence gap swallows
the result (HiFi: 8.8× larger than CoLoRd overall, almost entirely DNA). This is
the strongest argument for Lever 2 — the overlap codec is worth *more* in lossy
mode, which is the mode long-read archives actually ship in.

### Platform tables are not interchangeable

`--quality-bin ont` applied to HiFi data is a fidelity disaster: mae **42.84**,
99.4 % of bases changed, because it folds Q93 into the 26+ bin and destroys
HiFi's max-quality semantics (`binhifi` keeps Q93 exact → mae 14.33). On ONT the
two tables are byte-identical (ONT never reaches Q93). Consider warning when the
selected table mismatches the detected `Platform` (the `warn_redundant_binning`
pattern already exists). CoLoRd's lossy quality is still ~4.7× smaller than
`binhifi` on HiFi — worth a look once Lever 2 lands, but it is not the headline.

Cutpoints must ultimately be set by downstream fidelity (`concordance.sh`), not
raw ratio.

---

## Lever 1 — long-read quality tuning (cheap, high-leverage)

fqxv's quality context (`fqxv-fqzcomp`, `context()`) is
`q1(6b) | q2(4b) | delta(2b) | q3(2b) | pos_bucket(4b)` = 18 bits. Two pieces
are short-read-shaped:

- **`pos_bucket` saturates at ~base 224.** Over >99 % of a 14 kb read, all 4
  position bits are a constant — dead context capacity. ONT quality barely
  depends on absolute position anyway.
- **No surrounding-DNA-base context.** Both ENANO (2 prev-Q + 6 bases) and
  CoLoRd (6 prev-Q + 4 bases) mix neighbouring bases into the quality context
  because ONT quality tracks local sequence. fqxv uses none.

**Proposal.** Add a mode (auto-selected when `is_long_read`) that repurposes the
4 dead position bits for more useful context, gated on the long-read flag so the
short-read context stays byte-identical (preserving determinism and golden
outputs). Two variants, in increasing cost:

- *Local (cheap):* spend the freed bits on **more quality history** — widen `q2`,
  add a fourth previous quality, or a coarse "running average Q" bucket.
  Contained entirely to `fqxv-fqzcomp` + one header bit. `encode`/`decode` today
  take only `(lens, quals)`, so this needs no new inputs.
- *Base-context (the ENANO/CoLoRd trick, moderate):* mix **neighbouring DNA
  bases** into the quality context. This is the bigger win but is **not** local:
  `fqzcomp::encode`/`decode` receive no sequence bytes today, so the container
  must plumb the block's seq stream into the quality codec (both directions).
  A cross-cutting but well-bounded change.

**Binning — shipped.** `QualityBinning::Bin4/Bin8/Bin2` are Illumina-calibrated
(HiSeq / NovaSeq cutpoints, absolute Phred). On HiFi's narrow high-Q band they
collapse almost everything into the top bin. CoLoRd ships platform-specific
tables (ONT 4-level, HiFi 5-level with Q93 kept separate). fqxv now ships the
same two: `QualityBinning::BinOnt` and `BinHifi` (header tags 4 and 5), exposed
as `--quality-bin ont` and `--quality-bin hifi`. There is published evidence that
ONT Q-resolution can be reduced with little downstream effect (Bonito → 4
levels), so this is a safe, large lossy lever — measured at 3.4× on the ONT
quality stream above.

Effort: **small** (both changes are local to `fqxv-fqzcomp` + a header bit).

---

## Lever 2 — cross-read overlap sequence codec (big lever, big effort)

fqxv's sequence codec is a *within-read* order-k context model. It cannot
exploit **cross-read overlap**: at 65× coverage the same genomic locus is read
~65 times, and today each read is coded independently. NanoSpring (0.35–0.65
bits/base) and CoLoRd (down to ~0.01–0.02 bits/base at high coverage) get their
DNA wins entirely from modelling reads against each other.

### How the field does it

- **NanoSpring:** MinHash index (k=23, 60 hashes) finds candidate overlaps →
  minimap2 aligns → consensus graph; each read stored as (start pos on
  consensus, edits, RC flag, read index). DNA only — drops names + quality.
- **CoLoRd:** k-mer similarity graph; each read encoded as an **edit script**
  (anchors, matches, ins, del, subs) against a similar earlier read, with up to
  3 recursion levels; edit tuples entropy-coded with context from the previous
  tuple type + neighbouring DNA base. Reference set decays probabilistically to
  bound memory.

Both are the same shape: **find overlaps through noise → pick a reference →
encode the difference**.

### Proposed fqxv design

A new block-local method (either a `fqxv-lroverlap` crate → `rans`, `seq`, or a
long-read layout mode inside `fqxv-reorder`, which already has consensus /
assembly-window machinery). Per block (bounded by the 256 MiB seq budget, so
memory and parallelism stay bounded and deterministic):

1. **Overlap detection robust to ~10 % error.** (w, k) minimizers (small k so
   each read carries many); index minimizer → (read, pos); find candidate
   overlaps by **shared minimizer chains** (minimap2-style chaining), not single
   anchors. Reverse-complement aware via canonical minimizers.
2. **Layout.** Greedy overlap layout / union-find clustering ordered by chain
   score; deterministic tie-break on read index to keep output thread-order-free.
3. **Encode each read against its best predecessor:** (reference id or consensus
   offset, RC flag, edit script). Edit ops rANS-coded with a CoLoRd-style
   context (previous op type + local base). Bases with no good overlap fall back
   to the existing order-k context model — so the codec degrades gracefully to
   today's behaviour when coverage/redundancy is low.
4. **Container integration:** a new sequence-method tag (alongside the existing
   reorder reference methods), gated behind `is_long_read` and an effort level.

### Invariants to respect

- **Determinism:** minimizer index build and greedy layout must be order-free;
  seed all tie-breaks on read index. (Same discipline as the existing reorder
  codec.)
- **Block-local:** no cross-block references, so blocks stay independently
  decodable and `rayon`-parallel.
- **Graceful fallback:** a read with no confident overlap is coded exactly as
  today; the codec can never do *worse* than the order-k baseline by more than
  its small per-read header.

### Result (measured)

The design above is implemented in the **`fqxv-lroverlap`** crate: minimizers →
chain → overlaps → layout → consensus → `place_against` → edit scripts → rANS.
Measured on `ecoli_hifi` (120k reads, 1547 Mbase, ~300×, 5.16 Mb genome) by its
own harness, `crates/fqxv-lroverlap/examples/encode.rs`:

| path | seq bits/base |
| --- | --- |
| fqxv within-read order-k (the fallback) | 0.653 |
| **`fqxv-lroverlap`, whole-file harness** | **0.067** |
| CoLoRd | 0.068 |
| oracle (true reference, known placement) | 0.040 |

This is **parity with CoLoRd, not a win.** Across minimizer strides 4–14 the
result spans 0.067–0.072 bits/base, and that ~6% spread is sample noise — any
claim finer than "parity" is unsupported by this data. The oracle bound at 0.040
says roughly a further 40% is theoretically on the table; per the crate's own
analysis, that margin is the difference between coding against another erroneous
read (~0.005 edits/base — *both* reads' errors) and coding against a voted
consensus (0.0025).

### Wiring and the per-block coverage cap

`fqxv-lroverlap` is wired into the container. The block
sequence stream carries a leading method byte; long-read blocks (mean length
over 500 bp) code with both the overlap codec and the order-k model and keep the
smaller, so the overlap path never regresses a block. Selection is automatic —
no CLI flag — and the archive round-trips exactly (per-block content digest plus
`compress --verify`) and is byte-identical across thread counts.

**The wired codec runs per block, and that caps its coverage.** Each 256 MiB
block self-assembles its own reference — which is what preserves blocked
parallelism, per-block random access, and thread-determinism — but a block holds
only ~256 MiB of bases, so on a 5 Mb genome one block sees tens of ×, not the
whole file's 300×. The harness's 0.067 codes the *whole file* as one reference;
the container could not without giving up the per-block invariants. With a
per-block reference, a real round-trip-verified archive of the whole `ecoli_hifi`
file (120k reads, 1.55 Gbase, 6 blocks at ~52×) put the sequence stream at 0.107
bits/base — 6.1× smaller than the 0.653 fallback (total archive 4.04×) — but
above CoLoRd's whole-file 0.068, the gap being coverage per reference. That is
what the shared whole-file reference below closes: the same archive now measures
**0.065 bits/base at a 4.73× total, past CoLoRd's 0.068**.

**Shared whole-file reference (implemented).** The gap *is* reference
duplication: each block re-stored the same assembled genome, so a ~300× file kept
~6 copies where 1 is needed. The container now assembles **one** consensus over
the whole file, stores it once in a framed region between the header and the first
block (gated by the `GLOBAL_REFERENCE` feature bit and flag bit5), and codes every
block's reads against that frozen frame (sequence method byte 2,
`fqxv_lroverlap::encode_against`). Because placement is per-read against an
immutable frame, a read codes identically regardless of which block holds it —
**no block-boundary penalty** — so blocks stay 256 MiB, `rayon`-parallel, and
independently decodable given the shared frame. A whole-file never-worse gate
adopts the layout only when `reference frame + Σ chosen sequence` beats the plain
order-k total; otherwise no frame is written and the archive is the plain layout,
so it can only ever shrink. This is the same pattern the reorder path uses for its
global reference, and it removes the redundant reference copies without touching
the near-optimal per-read edit term. Measured on `hifi_40k` (516 Mbase, 2 blocks,
default order-11) the sequence stream drops **0.102 → 0.084 bits/base (−18%)**,
storing the ~5 Mb consensus once (a 1.26 MB frame) instead of per block; the win
widens with block count on deeper files. The compress path must **buffer** the
input for this (the streaming single-end path keeps the per-block method-1
fallback), and random-access single-stream projection of a shared-reference block
fails closed — it has no access to the frame. See issue #168.

### Cost / benefit

The benchmark settles the priority: the **entire** measured lossless gap to
CoLoRd is this stream (58.8M → 31.4M on offer, ~27M / 12% of the whole archive),
while quality is already at parity. So the overlap codec is the lever that
closes the CoLoRd gap — it is the larger project (overlap detection through
noise is the crux; too sparse misses overlaps, too dense explodes the index;
minimap2's chaining is the proven recipe) but it is *the* lever, not a
secondary one. Expect the payoff to grow on modern high-Q ONT and HiFi, where
the field reaches 0.35 bits/base and fqxv's within-read model cannot.

### Suggested sequencing

1. ~~**ONT/HiFi lossy binning tables**~~ (from Lever 1) — **done.** Smallest
   change, and the biggest *absolute* archive shrink since quality is 72–74% of
   the file. Still needs a downstream-fidelity check (`concordance.sh`) to set
   cutpoints honestly; the current cutpoints are CoLoRd's, not ones fqxv
   validated.
2. **Overlap sequence codec** (Lever 2) — **shipped.** The algorithm
   (`fqxv-lroverlap`) is wired into the container behind a sequence-method tag
   with `is_long_read` gating and a full decode path, auto-selected and kept only
   when it beats order-k. Two follow-ups have since landed on top: the consensus
   reference is assembled once for the whole file and stored once rather than
   per block, and ONT seeds with closed syncmers instead of window minimizers.
   What remains is consensus *quality* — the draft consensus still sits above a
   raw read's error rate, which bounds the ONT edit cost.
3. Quality base-context (rest of Lever 1) — **shipped**, and it overshot the
   estimate: the sequence-conditioned, context-mixed coder took HiFi quality
   below CoLoRd rather than the ~nil headroom predicted from this ONT data.

## Benchmark

Long-read runs are rows of the main parallel matrix rather than separate job
files — see `bench/README.md`. Both datasets live in `bench/panels/datasets.tsv`
and are compared against CoLoRd (`-q org`, lossless) plus gzip/zstd/xz:

- **`ecoli_ont`** (DRR205413: 21,140 reads, mean 14.2 kb, max 91.7 kb, mean
  Q≈11.5 — ragged, noisy older-basecaller ONT, ~65× E. coli).
- **`ecoli_hifi`** — the narrow-Q / low-error regime where the DNA lever matters
  most, subsampled from SRR11434954 (E. coli E2348/69 HiFi) to ~120k reads /
  ~1.5 Gbp (~300×) by `bench/slurm/prep.sbatch`.

The lossy binning tables are matrix points too, compared against `colord-lossy`.
Submit the matrix with `bench/scripts/submit_parallel.sh`. The `fqxv-lroverlap`
bits/base figures come from the crate's own harness rather than the matrix:

```bash
cargo run --release -p fqxv-lroverlap --example encode -- reads.fastq [ont|hifi]
```
