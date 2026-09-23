# Third-party notices

`fqxv` contains clean-room reimplementations of published compression algorithms.
No third-party source code is vendored — with one deliberate exception for
cryptographic primitives, see [Encryption](#encryption-fqxv-crypt) below. The
algorithms were implemented from public specifications and papers; we
acknowledge the original authors and the reference implementations we
cross-checked against for correctness.

## Encryption (`fqxv-crypt`)

Unlike every other entry in this file, `fqxv`'s optional passphrase encryption
does **not** reimplement its cryptographic primitives from a spec. Rolling
your own AEAD or password-based KDF is a well-known way to build something
that looks correct and isn't; this is a deliberate, explicit exception to the
workspace's otherwise strict clean-room policy (see `docs/design/encryption.md`
and CLAUDE.md's Architecture section, the `fqxv-crypt` bullet).

- **ChaCha20-Poly1305** (Bernstein's ChaCha20 stream cipher; RFC 8439's
  Poly1305-AEAD construction) — via the RustCrypto **`chacha20poly1305`**
  crate (https://github.com/RustCrypto/AEADs, MIT OR Apache-2.0). Used
  directly as a dependency, not reimplemented.
- **Argon2id** (Biryukov, Dinu & Khovratovich; the Password Hashing
  Competition winner; RFC 9106) — via the RustCrypto **`argon2`** crate
  (https://github.com/RustCrypto/password-hashes, MIT OR Apache-2.0). Used
  directly as a dependency, not reimplemented.
- **`getrandom`** (https://github.com/rust-random/getrandom, MIT OR
  Apache-2.0) — a thin OS-CSPRNG syscall wrapper (no algorithm of its own),
  used to draw the per-archive Argon2id salt and AEAD nonce identifier.
- **`zeroize`** (https://github.com/RustCrypto/utils, MIT OR Apache-2.0) —
  wipes derived key material from memory on drop.

See `docs/design/encryption.md` for the on-disk scheme these crates
implement.

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
- NanoSpring (Meng, Chandak et al.) — approximate-assembly long-read sequence
  compression; field context for the same lever (overlap index → align →
  consensus graph). Referenced for context only; not reimplemented from, and no
  source used.

## LZ compression (`fqxv-seq::lzma`)

- LZMA (Igor Pavlov, 7-Zip / LZMA SDK, public domain) — the *design* behind the
  clean-room LZ byte coder used for the long-range sequence path: adaptive
  bit-models, a 12-state machine, context literals with matched-byte prediction, a
  length coder, position-slot + aligned distance coding, and rep0–3 short codes.
  Not a port of liblzma or xz; the bit-models, match finder, parse, and stream
  layout are local, built only on this project's own range coder.

## Sequence alignment (`fqxv-align`)

- Wavefront alignment (Marco-Sola, Moure, Moreto & Espinosa, *Bioinformatics*
  2021) — the recurrences behind `wfa_align` / `wfa_align_opt`, whose work scales
  with alignment score rather than sequence length. Implemented from the paper;
  the reference implementation **WFA2-lib**
  (https://github.com/smarco/WFA2-lib, MIT) was consulted for behavior only, and
  no source was translated.
- Banded Needleman–Wunsch (Needleman & Wunsch 1970; Sellers 1974 edit-distance
  formulation) — `align_banded`, the predictable-cost reference the wavefront
  path is checked against. Textbook dynamic programming, implemented directly.

None of the implemented-from references above impose obligations beyond
attribution; all are permissive (BSD 3-Clause / MIT) or public domain. The
`fqxv-crypt` dependencies above are likewise permissively licensed (MIT OR
Apache-2.0). This project is licensed MIT OR Apache-2.0.
