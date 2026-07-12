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
sequence length; `rt` is a round-trip record-count sanity check.

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
