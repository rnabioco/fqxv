# fqxv benchmark harness (M0)

Baseline the field so every `fqxv` release can be measured against it. All tools
live in a pixi env; data and results live on `$SCRATCH` (never in the repo).

## Layout

```
bench/
  scripts/   the drivers        run_bench.sh, corpus.sh, edgecases.sh, sra_compare.sh,
                               bam_identity.sh, concordance.sh, fetch.sh, report.py, …
  slurm/     batch job files    bench.sbatch, bench_cell.sbatch, prep.sbatch, merge.sbatch, …
  panels/    input tables       datasets.tsv (ratio matrix), sra_panel.tsv (vs .sra)
  tools/     rustc helpers      fqdigest.rs (round-trip digest), bamcmp.rs (BAM digests),
                               fqsim.rs (synthetic per-platform FASTQ)
  RESULTS.md headline numbers
```

Scripts locate their siblings relative to their own path, so they can be invoked
from anywhere. The batch files resolve the harness root from `SLURM_SUBMIT_DIR`
(Slurm runs a spool copy, so `$0` is not the original path) and tolerate
submission from the repo root, `bench/`, or `bench/slurm/`.

## Environment

The repo has a **single** pixi manifest at the root with two environments:
`default` (docs toolchain) and `bench` (this harness). Everything below runs
under `-e bench`; `pixi` finds the manifest by searching upward, so the commands
work from the repo root or from `bench/`.

```bash
pixi install -e bench   # conda-forge + bioconda: sracha, spring, fqz_comp, colord, zstd, xz, pigz, seqkit, samtools, bwa, bcftools
```

## 1. Fetch data (login node is fine — it's I/O bound)

```bash
pixi run -e bench bash scripts/fetch.sh              # all accessions in datasets.tsv, raw .fastq to $SCRATCH/fqxv/data
pixi run -e bench bash scripts/fetch.sh SRR453566    # just one
```

Datasets span both quality regimes (the thing a quality codec cares about):

| accession | label | platform | quality |
| --- | --- | --- | --- |
| SRR2627175 | ecoli_miseq | MiSeq | full-range |
| ERR10213669 | human_novaseq_exome | NovaSeq 6000 | binned |
| DRR174812 | rnaseq_novaseq | NovaSeq 6000 | binned |
| SRR453566 | rnaseq_fullrange | GAIIx | full-range |

Check the quality alphabet of any file with `pixi run -e bench qstats <file.fastq>`.

## 2. Run the benchmark (compute node via Slurm — NOT the login node)

```bash
sbatch slurm/bench.sbatch                 # one exclusive amilan node, 64 threads
```

Or interactively:

```bash
srun --partition=amilan --qos=normal --nodes=1 --ntasks=1 --cpus-per-task=64 \
     --exclusive --mem=0 --time=02:00:00 --pty bash
pixi run -e bench bash scripts/run_bench.sh
```

Knobs (env vars): `FQXV_THREADS`, `FQXV_INPUT=r1|cat` (R1 only vs R1+R2
concatenated), `FQXV_TOOLS="gzip zstd19 fqz_comp"` (subset), `FQXV_DATA_DIR`,
`FQXV_RESULTS_DIR`.

### Which tools run

By default each dataset runs the tool set for **its platform**, defined in one
place — [`scripts/toolsets.sh`](scripts/toolsets.sh) — and shared by both the
sequential driver (`run_bench.sh`) and the parallel Slurm driver
(`submit_parallel.sh`). Illumina gets the full field matrix including
`fqxv-max`/`fqxv-shuffle` and SPRING's two lossy modes; ONT and PacBio get the
long-read set with CoLoRd and platform-calibrated quality bins.

Platform filtering keeps the matrix meaningful rather than merely large: SPRING
is Illumina-only, CoLoRd long-read-only, and `fqxv-reorder*` auto-disables above
the long-read length threshold (so it would just duplicate the plain `fqxv` rows
on a separate node). `fqz_comp` is deliberately kept in the long-read sets even
though it cannot parse long reads — it records `rt=no`, which distinguishes
"inapplicable" from "never tested".

Add a tool to `toolsets.sh` and every driver picks it up. `FQXV_TOOLS` overrides
the platform sets with one explicit list for every dataset; `FQXV_LR_TOOLS_ONT`
and `FQXV_LR_TOOLS_HIFI` override a single long-read platform.

## 3. Read the table

```bash
pixi run -e bench python scripts/report.py
```

`vs_gzip` is how many times smaller than plain `.fastq.gz` a tool is — the number
that actually matters for an archive. `bits/base` is the size normalized to
sequence length.

`rt` is a **content** round-trip check, not just a record count: the decompressed
output is reduced to a sorted multiset of `name / sequence / quality` tuples and
hashed against the input, so any corrupted base or quality fails it. It is
order-independent (SPRING and `fqxv --order any` reorder reads), and it excludes
the `+` line (fqxv normalizes it — the one documented lossy-by-design deviation).
For **lossy-quality** tools the check verifies the intended lossy output: the
fqxv `--quality-bin` rows are hashed against the input passed through that exact
bin table (a full end-to-end check of names + bases + *binned* quality), while
the SPRING lossy rows — whose internal tables we don't reproduce — are checked on
names + bases only. `det` is fqxv's thread-determinism check: the archive built
with `--threads 1` must be byte-identical to the many-threaded one (a core
invariant).

Every lossy row also reports a **quality-distortion** line (`Δqual`): the mean
absolute error, RMSE, and percentage of bases whose quality changed, measured
against the original qualities (records matched by name, so reordering is fine).
This is the fidelity half of the lossy tradeoff — read it next to the ratio, not
instead of it.

The fqxv rows also print a per-stream breakdown (`names / seq / qual` bytes, from
`fqxv info --tsv`) so you can see which stream to invest in. The matrix runs
several fqxv points — `fqxv` (level 5), `fqxv9` (level 9), `fqxv-reorder`
(`--order any`), `fqxv-max` (`--max`, i.e. `-l 9 --order any` — the advertised
best-ratio preset combining the deepest sequence context *and* read reordering),
and the lossy quality sweep `fqxv-bin8` / `fqxv-bin4`
/ `fqxv-bin2` — plus a `fqxv-paired` self-check that compresses R1+R2 as one
spot-interleaved archive. For a **like-for-like** lossy comparison the matrix also
runs SPRING's own binning: `spring-illbin` (`-q ill_bin`, Illumina 8-level —
compare to `fqxv-bin8`) and `spring-binary` (`-q binary 25 37 15`, 2-level —
compare to `fqxv-bin2`). SPRING is the only field tool with Illumina-comparable
binning; `fqz_comp`/`fqzcomp5` have no such mode, so they run lossless only and
`fqxv-bin4` has no competitor row.

## Baselines

- **gzip** (`pigz -6`) — the baseline of record.
- **zstd -19 --long**, **xz -9** — general-purpose strong baselines.
- **fqz_comp** — FASTQ-specific quality/context coder (htscodecs family).
- **fqzcomp5** — Bonfield's newer quality/context + LZP coder (built from source).
- **spring** — reference-free read-reordering archiver.

### From-source tools

`fqzcomp5` and `PgRC` aren't in bioconda; build them with:

```bash
pixi run -e bench build-tools        # clones + compiles into $SCRATCH/fqxv/tools/bin
```

`fqzcomp5` is a full lossless FASTQ compressor and is in the default tool set.

**PgRC is deliberately *not* in the default tool set.** It is a read/*sequence*
compressor: by default it drops read names and simplifies quality, and its
decompressed output is one sequence per line (not FASTQ), so its ratio is not
comparable to full-FASTQ archivers. It's kept built for sequence-stream
experiments (relevant to the M4 reordering work):

```bash
PgRC -o -Q -t 8 -i reads.fastq arch     # compress (‑o keep order, ‑Q lossless qual)
PgRC -d -t 8 arch                        # -> arch_out (sequences, one per line)
```

Once `fqxv` produces output, it joins the comparable table as another tool.

## fqxv vs the archive the data shipped in (`scripts/sra_compare.sh`)

`.sra` is the format most public data actually ships in, and it is itself a
compressed columnar store (2-bit-packed bases, a quality model, spot/name
columns). For anyone pulling from SRA/ENA the honest question is "is an `.fqxv`
worth keeping instead of the `.sra` I'd otherwise store?" — not just "how does it
compare to gzip". This records, per run, the `.sra` size and the fqxv size at
several operating points, and reports **fqxv/.sra** (< 1.0 = fqxv is smaller).

The `.sra` size comes from `sracha info --format tsv` — **metadata only**, so the
multi-GB archive is never downloaded. That splits the run in two, because the
lookup needs the network and the compression needs a compute node:

```bash
pixi run -e bench bash scripts/sra_compare.sh sizes   # login node: cache .sra sizes
# then inside an srun/sbatch allocation:
pixi run -e bench bash scripts/sra_compare.sh run     # compute: compress + join
# results -> $FQXV_RESULTS_DIR/sra_compare.tsv, rendered by scripts/report.py
```

The panel is `panels/sra_panel.tsv`; `SRA_POINTS` selects the fqxv points
(default `max bin8 bin4 bin2`) and `SRA_INCLUDE_LARGE=1` opts into the big rows.
Compare **lossless fqxv vs `.sra` first** — `.sra` bundles both mates, read names
and spot metadata, and its quality handling is platform-adaptive (binned Illumina
vs full-range), so match the fqxv point to the `.sra`'s own quality regime before
reading anything into the lossy rows. Reference-compressed cSRA runs aren't
reconstructable to raw FASTQ and are flagged/excluded.

## Downstream fidelity — variant-call concordance (`scripts/concordance.sh`)

`scripts/run_bench.sh` reports the *rate* of lossy quality binning (ratio) and its raw
*distortion* (`Δqual`). `scripts/concordance.sh` answers the question that actually
matters for lossy quality: **does binning change variant calls?**

For each dataset with a `reference` genome (the last column of `panels/datasets.tsv`),
it aligns the lossless reads and each `fqxv --quality-bin` output to the
reference, calls variants (`bwa mem` → `bcftools mpileup`/`call`), and reports
SNP and indel **recall** and **precision** of the binned calls against the
lossless baseline. The binned FASTQ is produced by the real codec (compress
`--quality-bin` then decompress), so the concordance is exactly what fqxv would
store. It is **heavy** (alignment + calling) and is *not* part of `scripts/run_bench.sh`.

```bash
# inside an srun/sbatch allocation:
pixi run -e bench bash scripts/concordance.sh              # every dataset with a reference
pixi run -e bench bash scripts/concordance.sh ecoli_miseq  # just one
# results land in $FQXV_RESULTS_DIR/concordance.tsv
```

Read the concordance next to the ratio — a bin that compresses well but drops
SNP precision (coarser binning tends to inflate false-positive calls) is not
free. Only DNA-resequencing datasets carry a reference; RNA-seq needs
splice-aware calling and is left out (`-` in `panels/datasets.tsv`). Extra pixi deps:
`bwa`, `bcftools` (`samtools` was already present).

## Downstream fidelity — BAM round-trip (`scripts/bam_identity.sh`)

Where `scripts/concordance.sh` asks "do variant calls change?", `scripts/bam_identity.sh` asks
the sharper question: **is the aligned BAM itself identical before vs after
fqxv?** It aligns the original reads and each fqxv round-trip (`bwa mem`) and
compares the alignments with `bamcmp` — an `rustc -O` streaming tool (sibling of
`tools/fqdigest.rs`) that emits order-independent multiset digests of `samtools view`
in one pass (no `samtools sort`, no `sort | md5sum`), plus a fast `qualdelta`.
Three digests per BAM: `content` (whole record), `body` (drop QNAME — survives
renaming), `place` (FLAG..SEQ — survives renaming *and* quality binning).

```bash
# inside an srun/sbatch allocation:
sbatch slurm/bam_identity.sbatch              # default dataset (ecoli_miseq)
sbatch slurm/bam_identity.sbatch ecoli_miseq  # any datasets.tsv row with a reference
# results -> $FQXV_RESULTS_DIR/bam_identity.tsv
```

On E. coli MiSeq (2.19 M reads): **lossless is byte-identical at the BAM level**
(content/body/place and the coordinate-sorted file all match); the reorder modes
preserve output order on real SRA data (`order_changed=no`) so they are identical
too; and lossy `--quality-bin` never moves a read (`place` matches) — only QUAL
changes, by a bounded amount (`qualdelta`). A `reorder-forced` control shuffles
the identical read set and realigns: it exposes that **`bwa mem` itself is
order-sensitive** (~1.2% of reads realign differently — deterministic, not
threading, unaffected by `-K`), which is why preserving read order is what makes
a BAM reproducible. Tools: `bwa`, `samtools`, `rustc`.

## Robustness corpus (`scripts/corpus.sh`) — correctness, not ratio

`scripts/run_bench.sh` is a curated *ratio* comparison against the field.
`scripts/corpus.sh` is the opposite: a *correctness net* that throws a random pile of
real SRA runs at fqxv and asserts every archive round-trips and builds
identically regardless of thread count. It's how we shake out the messy long
tail (odd read lengths, wide/narrow quality alphabets, Ns, empty reads, long
names) that four curated datasets never hit. Modeled on sracha-rs'
`validation/random_corpus.sh`.

For each accession, per compression **mode** it runs:

- `compress` → `decompress` → order-independent **content** round-trip
  (name/seq/qual multiset, `+` line excluded — fqxv's one documented deviation;
  reordering-safe for `--order any`);
- `compress --threads 1` → **byte-compare** against the many-threaded archive
  (the determinism invariant).

Modes: `default` (`-l5` lossless), `max` (`-l9 --order any`), and `bin8`
(lossy Illumina 8-level, checked against the *binned-expected* content). Set
`FQXV_MODES` to change. Outcomes: `PASS` / `FAIL_RT` / `FAIL_DET` /
`FAIL_COMPRESS` / `FAIL_DECOMPRESS` / `ERROR_FETCH`.

```bash
pixi run -e bench bash scripts/corpus.sh sample -n 20 -s 42   # random ENA accessions -> $CORPUS_DIR/accessions.txt
pixi run -e bench bash scripts/corpus.sh build-digest         # compile fqdigest (fast round-trip hash; optional)
pixi run -e bench bash scripts/corpus.sh fetch                # sracha get every accession (login node — I/O bound)
pixi run -e bench bash scripts/corpus.sh sbatch               # fan out as a slurm array (COMPUTE) — one accession/task
pixi run -e bench bash scripts/corpus.sh summary              # pass/fail tally + failing accessions
```

`sample` draws random Illumina runs by default (`-p all` for every platform);
seed the shuffle for a reproducible corpus. Everything lands under
`FQXV_CORPUS_DIR` (default `$SCRATCH/fqxv/corpus`): `data/`, per-accession
`logs/`, and `results.tsv`. On a **`FAIL_*`** the offending archive is kept in
`work/<acc>/` for post-mortem (`fqxv info`, re-decompress with `-v`).

`tools/fqdigest.rs` is the round-trip hash: an order-independent multiset digest
(sum of per-record hashes — commutative, so read order doesn't matter) in one
O(n) streaming pass with bounded memory, replacing a slow `awk | sort | md5sum`
(~4× faster, and it sidesteps awk's per-byte quality-binning loop entirely).
`scripts/corpus.sh` uses it when built and falls back to the awk pipeline otherwise, so
results are identical either way. Build it once with `corpus.sh build-digest`
(needs `rustc` on PATH; the harness never builds on the login node otherwise).

## Synthetic input (`tools/fqsim.rs`)

For tests that need FASTQ but not a specific accession — interrupt handling,
edge cases, quick ratio sanity checks — `fqsim` generates it far faster than
fetching, and reproducibly from a seed.

```bash
rustc -O -o "$SCRATCH/fqxv/tools/bin/fqsim" bench/tools/fqsim.rs

fqsim --platform novaseq --reads 1000000 --paired sample   # sample_1/_2.fastq
fqsim --platform ont  --reads 20000 --coverage 25 -o ont.fastq
fqsim --platform hifi --reads 20000 --coverage 25 -o hifi.fastq
```

Platforms (`novaseq`, `hiseq`, `ont`, `hifi`) carry their own read length,
error mix, quality model and name layout: NovaSeq emits 4-level RTA3-binned
quality where HiSeq emits a continuous range, and the long-read profiles raise
the indel rate inside homopolymer runs. Individual knobs (`--len`,
`--sub-rate`, `--ins-rate`, `--del-rate`) override the profile.

The important part is that reads are **sampled from a generated reference
genome** at a chosen `--coverage`, not drawn independently. Independent random
reads share no sequence, so the codecs built to exploit cross-read redundancy
have nothing to find and measure the same as the plain context coder — which
makes such data actively misleading. With genome-backed sampling the redundancy
is real and coverage-tunable: on 400k NovaSeq reads, `--order any` beats
`--order preserve` by 47% at 5× coverage and 31% at 60×.

Output is plain FASTQ; pipe to `bgzip`/`gzip` if you need it compressed. Same
`--seed` gives byte-identical output, so a failing case can be reproduced
without staging a data file.
