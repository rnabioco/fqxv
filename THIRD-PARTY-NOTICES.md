# Third-party notices

`fqxv` contains clean-room reimplementations of published compression algorithms.
No third-party source code is vendored. The algorithms were implemented from
public specifications and papers; we acknowledge the original authors and the
reference implementations we cross-checked against for correctness.

## CRAM 3.1 codecs (rANS Nx16, fqzcomp quality model, name tokenizer)

- Specification: CRAM codecs specification, `samtools/hts-specs`
  (https://samtools.github.io/hts-specs/CRAMcodecs.pdf).
- Reference C implementation: **htscodecs** (https://github.com/samtools/htscodecs),
  © Genome Research Ltd, BSD 3-Clause. Author: James Bonfield.
- Reference Rust implementation cross-checked for test vectors:
  **noodles-cram** (https://github.com/zaeleus/noodles), © 2018 Michael Macias, MIT.

## rANS

- The rANS entropy-coder design derives from Jarek Duda's asymmetric numeral
  systems and Fabien Giesen's `ryg_rans` (public domain / CC0).
- Range-coder design after Eugene Shelwien (public domain).

## Read reordering (sequence stream)

- PgRC2 (Kowalski & Grabowski, *Bioinformatics* 2025) and SPRING
  (Chandak et al., *Bioinformatics* 2019) — algorithmic references for the
  pseudogenome / read-reordering engine. Reimplemented from the papers.

## Long reads (quality binning, overlap codec)

- CoLoRd (Kokot, Gudyś, Li & Deorowicz, *Nature Methods* 2022, MIT) —
  algorithmic reference for long-read compression. The `--quality-bin ont` and
  `--quality-bin hifi` cutpoints follow its platform-specific quality tables,
  and its edit-script sequence model is the reference for the `fqxv-lroverlap`
  overlap work. Reimplemented from the paper.
- minimap2 (Heng Li, *Bioinformatics* 2018, MIT) — the minimizer-index and
  colinear-chaining design that `fqxv-lroverlap`'s overlap detection follows.
  Reimplemented from the paper; no source translated.
- miniasm (Heng Li, *Bioinformatics* 2016, MIT) — overlap–layout–consensus
  reference used to check the long-read assembly's collapse.

None of the above impose obligations beyond attribution; all are permissive
(BSD 3-Clause / MIT) or public domain. This project is licensed MIT OR Apache-2.0.
