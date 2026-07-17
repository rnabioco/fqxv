# Acknowledgments

`fqxv` builds on a large body of prior work in genomic and general-purpose
compression. Every codec in this project is a **clean-room implementation** from
public specifications and papers — no third-party source code is vendored — but
these projects and their authors made the work possible, and we cross-checked
against several of them to verify correctness.

## Entropy coding

- **htscodecs** — [samtools/htscodecs](https://github.com/samtools/htscodecs),
  by James Bonfield (Genome Research Ltd), and the [CRAM 3.1 codecs
  specification](https://samtools.github.io/hts-specs/CRAMcodecs.pdf). These are
  the reference for our `fqxv-rans` (rANS Nx16) coder, the `fqxv-fqzcomp`
  quality model, and the `fqxv-tokenizer` read-name tokenizer.
- **rANS** — Jarek Duda's work on asymmetric numeral systems, and Fabien
  Giesen's [`ryg_rans`](https://github.com/rygorous/ryg_rans) (public domain /
  CC0), which shaped the interleaved-state design of our rANS coder.
- **Range coding** — Eugene Shelwien's range-coder design (public domain)
  underpins `fqxv-range`.
- **noodles** — [zaeleus/noodles](https://github.com/zaeleus/noodles), by
  Michael Macias — a Rust CRAM implementation we cross-checked test vectors
  against.

## Quality-score compression

- **fqzcomp** — the quality-score context model by James Bonfield that our
  `fqxv-fqzcomp` codec is modeled on.

## Sequence reordering

- **SPRING** — Chandak, Tatwawadi, Ochoa, Hernaez & Weissman,
  *Bioinformatics* 2019.
- **PgRC2** — Kowalski & Grabowski, *Bioinformatics* 2025.

These are the algorithmic references for the read-reordering engine in
`fqxv-reorder`, reimplemented from the papers.

## Long reads

- **CoLoRd** — Kokot, Gudyś, Li & Deorowicz, *Nature Methods* 2022. The
  algorithmic reference for long-read compression, and the tool we measure
  against in [Long-read support](design/longread.md). Its platform-specific
  quality tables are the source of the cutpoints in our `--quality-bin ont` and
  `--quality-bin hifi` bins, and its edit-script approach — code each read
  against a similar earlier read — is the shape the `fqxv-lroverlap` sequence
  work follows.
- **minimap2** — Heng Li, *Bioinformatics* 2018. Its minimizer indexing and
  colinear-anchor chaining are the proven recipe for finding overlaps through a
  noisy channel; `fqxv-lroverlap`'s chainer follows that shape.
- **miniasm** — Heng Li, *Bioinformatics* 2016. The overlap–layout–consensus
  reference we checked our assembly's collapse against.
- **NanoSpring** — approximate-assembly long-read sequence compression; field
  context for the same lever (overlap index → align → consensus graph).

As with everything else here, these are reimplemented from the published papers
— no CoLoRd, minimap2, or miniasm source is vendored or translated.

## Licenses

All of the above are permissive (BSD 3-Clause / MIT) or public domain and impose
no obligations beyond attribution. See
[`THIRD-PARTY-NOTICES.md`](https://github.com/rnabioco/fqxv/blob/main/THIRD-PARTY-NOTICES.md)
in the repository for the full attribution and license details. `fqxv` itself is
dual-licensed **MIT OR Apache-2.0**.
