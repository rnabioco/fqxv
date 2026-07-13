# fqxv benchmark harness (M0)

Baseline the field so every `fqxv` release can be measured against it. All tools
live in a pixi env; data and results live on `$SCRATCH` (never in the repo).

## Environment

```bash
cd bench
pixi install        # solves conda-forge + bioconda: sracha, spring, fqz_comp, zstd, xz, pigz, seqkit, samtools
```

## 1. Fetch data (login node is fine — it's I/O bound)

```bash
pixi run bash fetch.sh              # all accessions in datasets.tsv, raw .fastq to $SCRATCH/fqxv/data
pixi run bash fetch.sh SRR453566    # just one
```

Datasets span both quality regimes (the thing a quality codec cares about):

| accession | label | platform | quality |
| --- | --- | --- | --- |
| SRR2627175 | ecoli_miseq | MiSeq | full-range |
| ERR10213669 | human_novaseq_exome | NovaSeq 6000 | binned |
| DRR174812 | rnaseq_novaseq | NovaSeq 6000 | binned |
| SRR453566 | rnaseq_fullrange | GAIIx | full-range |

Check the quality alphabet of any file with `pixi run qstats <file.fastq>`.

## 2. Run the benchmark (compute node via Slurm — NOT the login node)

```bash
sbatch bench.sbatch                 # one exclusive amilan node, 64 threads
```

Or interactively:

```bash
srun --partition=amilan --qos=normal --nodes=1 --ntasks=1 --cpus-per-task=64 \
     --exclusive --mem=0 --time=02:00:00 --pty bash
pixi run bash run_bench.sh
```

Knobs (env vars): `FQXV_THREADS`, `FQXV_INPUT=r1|cat` (R1 only vs R1+R2
concatenated), `FQXV_TOOLS="gzip zstd19 fqz_comp"` (subset), `FQXV_DATA_DIR`,
`FQXV_RESULTS_DIR`.

## 3. Read the table

```bash
pixi run python report.py
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
(`--order any`), and the lossy quality sweep `fqxv-bin8` / `fqxv-bin4`
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
pixi run build-tools        # clones + compiles into $SCRATCH/fqxv/tools/bin
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

## Downstream fidelity — variant-call concordance (`concordance.sh`)

`run_bench.sh` reports the *rate* of lossy quality binning (ratio) and its raw
*distortion* (`Δqual`). `concordance.sh` answers the question that actually
matters for lossy quality: **does binning change variant calls?**

For each dataset with a `reference` genome (the last column of `datasets.tsv`),
it aligns the lossless reads and each `fqxv --quality-bin` output to the
reference, calls variants (`bwa mem` → `bcftools mpileup`/`call`), and reports
SNP and indel **recall** and **precision** of the binned calls against the
lossless baseline. The binned FASTQ is produced by the real codec (compress
`--quality-bin` then decompress), so the concordance is exactly what fqxv would
store. It is **heavy** (alignment + calling) and is *not* part of `run_bench.sh`.

```bash
# inside an srun/sbatch allocation:
pixi run bash concordance.sh              # every dataset with a reference
pixi run bash concordance.sh ecoli_miseq  # just one
pixi run python -c 'pass'  # results land in $FQXV_RESULTS_DIR/concordance.tsv
```

Read the concordance next to the ratio — a bin that compresses well but drops
SNP precision (coarser binning tends to inflate false-positive calls) is not
free. Only DNA-resequencing datasets carry a reference; RNA-seq needs
splice-aware calling and is left out (`-` in `datasets.tsv`). Extra pixi deps:
`bwa`, `bcftools` (`samtools` was already present).
