# Golden archive fixtures

Small `.fqxv` archives written by an earlier release, kept so later releases can
prove they still decode them. `manifest.tsv` records what each one is expected to
decode to (byte length and CRC-32C of the decoded stream); the archives
themselves are the fixture, the expected plaintext is not checked in.

Consumed by `../golden.rs`. Also seeds the container fuzz target's corpus, which
otherwise starts from nothing and has to guess the magic and a valid header CRC
before it reaches any decode path.

## Why these exist

Every other test in the workspace compares a build against itself — the
round-trips, the proptests, the fuzz targets, the thread-determinism
byte-comparisons, the accession corpus (which recompresses from FASTQ each run).
All of them are invariant to a change that moves the encoder and decoder
together, and that is exactly the change that breaks archives already on disk.
Reordering context-model initialization in `fqxv-seq`, altering rANS frequency
normalization, repacking a fqzcomp context, changing a minimizer parameter: each
round-trips perfectly, passes everything, and leaves every existing `.fqxv`
either rejected or misdecoded. These fixtures are the only thing that notices.

## Coverage

One fixture per on-disk layout and per sequence codec, because a fixture whose
input is too thin falls back to the plain order-k path and pins nothing the
`plain_se` fixture already covers:

| fixture | what it pins |
| --- | --- |
| `plain_se` | plain layout, order-k sequence codec — the default path |
| `paired_labels` | `G=2` interleaving, and the header extension region (non-critical tag `0x01`, member labels) |
| `reorder` | the whole-file global-cluster reorder layout, which has its own header and framing |
| `longread` | the long-read cross-read overlap codec (a sequence method byte other than 0) |

`golden.rs` asserts this set is present, so dropping one cannot quietly leave a
layout uncovered while the suite still passes.

## Regenerating

```sh
cargo run --release -p fqxv --example make_fixtures -- crates/fqxv/tests/fixtures
```

**Regenerating is almost always the wrong response to a failing golden test.** A
failure means this build no longer reproduces what an earlier one wrote, which is
the regression the fixtures exist to catch; regenerating makes the red go away
and ships it. Regenerate only when adding a fixture, or when a format change has
been made deliberately and gated per the evolution policy in
`docs/design/container.md` — and then keep the old fixtures as well, since the
whole point is that they must still decode.

Adding a fixture for a newly added layout or codec is encouraged. Adding one for
a *release* is not necessary: these pin on-disk bytes, and the bytes only change
when a codec or the container does.
