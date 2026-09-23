# Native Archive Encryption

> **Status: shipped.** The container format, `fqxv-crypt`, and the `fqxv`
> library's compress/decompress/inspect/verify/random-access paths implement
> the design below. This document is the design record: read it before
> touching the `0x81` header extension, the footer's encrypted shape, or the
> nonce/AAD scheme in `fqxv-crypt`.

## Problem

`.fqxv` archives can carry sensitive genomic data — unpublished results,
patient or participant samples, anything under an embargo — and until now the
format had no confidentiality mechanism. The workaround (pipe the archive
through `age`/`gpg`) works, but it throws away everything the container format
gives you: `fqxv info`/`inspect` can no longer read the header without first
decrypting the whole file, `decompress --recover` has no block structure to
recover into, and the archive stops being a self-identifying `.fqxv` file at
all. This adds a real, native `--encrypt` capability instead: the ciphertext
lives inside the container's own framing, gated so an old reader refuses
cleanly rather than misdecoding it.

## Scope decisions

1. **Native format support**, not an external-tool wrapper. Encryption is a
   real feature of the container — a header feature bit, a critical extension
   record, per-block ciphertext — not something layered on top of the file.
2. **Passphrase-based** key derivation via Argon2id. No recipient/public-key
   mode, no raw keyfile mode, in this release (see [Out of
   scope](#out-of-scope-for-this-release)).
3. **Whole-block encryption, not per-stream.** `fqxv-python`'s `fqxv.remote`
   per-stream column-projection random access is explicitly *not* preserved on
   encrypted archives — a single stream (names, sequence, or quality) cannot
   be fetched or authenticated independently of the block that contains it
   once that block is sealed as one ciphertext. `random_access::Index::read`
   refuses an encrypted archive outright (`Error::EncryptedArchiveNotSupported`)
   rather than exposing a misleading partial projection.
4. **Real Cargo dependencies for the cryptography** — RustCrypto's
   `chacha20poly1305` and `argon2` — rather than a clean-room reimplementation.
   This is a deliberate, explicit exception to the rest of the workspace's
   clean-room codec policy (CLAUDE.md, `THIRD-PARTY-NOTICES.md`): hand-rolling
   AEAD/KDF primitives carries constant-time, nonce-misuse, and side-channel
   risks that the "reimplement from a published spec" methodology used
   everywhere else in this codebase was never built to guarantee. `fqxv-crypt`
   is the one crate in the DAG that takes on this dependency, quarantined so
   no codec crate is aware encryption exists or gains a transitive crypto
   dependency (CLAUDE.md's Architecture section).

## The `fqxv-crypt` crate

A leaf crate wrapping RustCrypto directly: `chacha20poly1305` (AEAD),
`argon2` (KDF, `default-features = false` — the low-level
`Argon2::hash_password_into` API needs only the `alloc` feature, not the
`password-hash` PHC-string machinery, since `fqxv` owns on-disk salt storage
itself), `getrandom` (CSPRNG for the salt/nonce material), and `zeroize`
(wipe key material on drop). See its module doc comment
(`crates/fqxv-crypt/src/lib.rs`) for the full public API
(`KdfParams`, `Passphrase`, `ArchiveKey`, `derive_key`, `random_bytes32`,
`ArchiveCipher`).

`KdfParams::DEFAULT` is RFC 9106 §4's "less memory available" profile — 64
MiB, 3 iterations, 4 lanes — deliberately not the stronger 2 GiB profile: KDF
cost is paid once per archive *open* (not per block, so it never sits in the
per-block `rayon` hot path), but `fqxv` runs everywhere from a laptop to a
memory-constrained CI runner, and 2 GiB just to open one archive is a real
availability cost for marginal extra resistance over a threshold that is
already far beyond a realistic cracking budget here. The chosen params are
recorded per archive (not baked into a build), so a later release can raise
the default without breaking any archive already written.

## On-disk layout

### Header: feature bit + critical extension

`crates/fqxv/src/lib.rs` defines:

```rust
pub const ENCRYPTED: u64 = 1 << 1;
```

in `feature`, folded into `KNOWN_FEATURES`. An archive with this bit set
cannot be decoded at all without a passphrase, so — exactly like
`GLOBAL_REFERENCE` before it — it is a `required_features` bit: an old reader
refuses at `read_header`, before a single block is touched, rather than
attempting to entropy-decode ciphertext as if it were a coded stream. This is
the format's second worked example of [`container.md`'s evolution
policy](container.md#versioning-and-evolution-policy): "capability → feature
bit."

The accompanying metadata rides a new **critical** header extension tag,
`0x81` (`EXT_CRITICAL_BIT | 0x01`, `format.rs::EXT_TAG_ENCRYPTION`) — critical,
unlike the non-critical member-label tag, because there is no safe way to
"decode without it": without the salt/nonce/KDF params, no block in the
archive can be opened. Its 42-byte payload:

```text
[1]  crypt_version        -- 1 = ChaCha20-Poly1305 + Argon2id (the only value
                              this build writes or understands; an unknown
                              value is Error::UnsupportedEncryptionVersion,
                              distinct from an unrecognized *tag*)
[16] salt                 -- Argon2id salt, random per archive
[16] nonce_id             -- random per archive; seeds every block/footer/
                              reference-frame nonce and AAD (see below)
[4]  argon2_m_cost_kib (LE u32)
[4]  argon2_t_cost     (LE u32)
[1]  argon2_p_cost      (u8)
```

`read_header` cross-checks the feature bit against the extension's presence —
one set without the other is a structural `Malformed` error, never silently
tolerated (`format.rs`, right after the unknown-flags check). `format.rs`'s
`setup_encryption` (compress side) and `open_encryption` (decode side) are the
only places that turn a passphrase into an `ArchiveCipher`; the container code
elsewhere only ever holds `Option<&ArchiveCipher>`.

### Block payload: sealed, framing unchanged

The outer block frame (`[4 BLOCK_MAGIC][8 payload_len][4 crc32c][payload]`) is
byte-for-byte identical to the plaintext layout. For an encrypted archive,
`payload` is `ArchiveCipher::seal_block(block_index, plaintext_payload)` —
the same `[24 digests][4 n_reads]([4 len][bytes])×3` bytes any archive
produces, now ciphertext-plus-tag. The frame's own `crc32c` still covers
exactly those on-disk bytes, checked before the (much more expensive) AEAD
open — a cheap, keyless first-pass corruption filter that coexists with AEAD
rather than being superseded by it: AEAD is the *authenticity* check (defeats
a deliberate adversary), the frame CRC is the cheap *accidental-corruption*
filter, and there's no reason to give up the second for having the first.

### Footer: same stride, degraded per-stream index, plus an authentication tag

The footer keeps its existing wire shape and per-group stride (`[8 offset][4
read_count]` + three `[8 offset][4 len][4 crc32c]` `StreamLoc` triples,
`FOOTER_GROUP_BYTES` = 60) — no new distinct encrypted-footer shape. What
changes is what's written into those triples: `block::write_blocks`'s
`whole_block_stream_locs` gives all three `StreamLoc` entries the *same*
`(offset, len, crc32c)`, pointing at the whole sealed block rather than at
three real sub-block ranges. That's the honest value to record once a stream
can't be fetched or authenticated independently of its block — inventing
per-stream offsets inside ciphertext would be actively misleading, and
duplicating the whole-block location keeps the on-disk stride (and every
existing size-exactness check in `parse_footer_body`) unchanged rather than
adding yet another footer shape to validate. `inspect`'s encrypted branch
sums just one of the three (identical) entries per group into
`Info::encrypted_bytes`, rather than fabricating a names/seq/qual split that
cannot be known without decrypting every block.

A second, real addition: a 16-byte **footer authentication tag**
(`ArchiveCipher::mac_footer`/`verify_footer`), inserted between
`whole_file_crc` and `footer_crc` when the archive is encrypted
(`format.rs::footer_crc_tail`). It authenticates the whole footer body — the
row-group offsets, read counts, and `total_reads` — a tail-truncation
attacker would otherwise leave unchecked (see [Nonce/AAD
scheme](#nonceaad-scheme-and-what-it-does-and-does-not-defeat) below).
`footer_crc` covers the tag too, so its own bytes stay bit-rot-protected like
the rest of the body. Checking it needs the passphrase (`verify`'s
`opts.password`), so `verify`/`verify_report` treat it as an *additional*
check layered on the existing keyless ones, explicitly reported as skipped
— not silently omitted — when no password is supplied.

### Reference frame: sealed when present

The plain layout's optional whole-file long-read reference frame
(`FLAG_GLOBAL_REFERENCE`) carries real assembled consensus sequence, so it is
sealed too whenever both it and `feature::ENCRYPTED` are set —
`read_reference_frame`/its writer gain an `Option<&ArchiveCipher>` branch,
using a reserved nonce/AAD index (`GREF_INDEX = u64::MAX - 1`) distinct from
any real block ordinal.

### Reorder layout: out of scope

`FLAG_GLOBAL_REORDER` (`--order any`/`shuffle`/`--max`) is a different,
footer-less, many-small-frames shape (flip bitmap, permutation, name
template, per-block sequence/names/quality frames, trailing output digest)
that would need its own nonce/AAD and tail-truncation design — not a drop-in
reuse of the plain layout's scheme. `encode_reordered` (the one choke point
every reorder entry point passes through) rejects `Params { encrypt:
Some(_), .. }` with a clear `Malformed` error before any I/O; see [Out of
scope](#out-of-scope-for-this-release).

## Nonce/AAD scheme, and what it does and does not defeat

`fqxv-crypt::ArchiveCipher` builds every nonce and AAD internally from the
archive's `nonce_id` and a position index — callers never construct one:

- `nonce(idx) = nonce_id[0..4] ‖ LE64(idx)` (96-bit ChaCha20-Poly1305 nonce).
  A block's `idx` is its 0-based footer ordinal, unique and strictly
  increasing within one archive by construction, so `(key, nonce)` never
  repeats within an archive — the property a stream cipher's nonce must have.
  It's a pure function of position rather than a stored random value, so it
  costs nothing on disk.
- `aad(domain, idx) = domain ‖ nonce_id ‖ LE64(idx)`, with distinct 8-byte
  domain tags for blocks, the footer, and the reference frame. This binds
  three things: domain separation (a block ciphertext can't be replayed as a
  footer tag or reference frame even though all three share one key), archive
  identity (`nonce_id`), and — for blocks — *position*.

**Position-bound AAD is what defeats splicing, reordering, and duplication
for the primary decode path.** `decompress`/`decompress_split` read blocks
strictly sequentially from the start of the block region (`for_each_block_batch`
tracks the ordinal as it reads, exactly the index a block was sealed under);
they never skip past a failure. So for that path, "the block's position in
the byte stream" and "the ordinal it was sealed under" are the same thing by
construction — an attacker who moves, drops, or duplicates a block changes
what a decoder encounters at that stream position, and the recomputed AAD no
longer matches, so authentication fails before any bytes are trusted.

**It does not, on its own, catch dropping only the archive's *trailing*
blocks behind a forged terminator** — every surviving block still sits at its
correct ordinal, so each still authenticates. `mac_footer`/`verify_footer`
(above) closes that specific gap for any reader that checks the footer with a
passphrase.

**Marker-scan recovery (`decompress_recover`'s footer-unreadable fallback) is
refused for encrypted archives — not merely out of scope, but actually
unsafe to attempt naively.** `recover_via_scan` finds blocks by scanning for
the `BLOCK_MAGIC` sync marker, and its block-ordinal counter only advances on
a *successfully CRC-validated* frame; a frame whose CRC fails is silently
absorbed as noise and the counter does not advance past it. For a plaintext
archive that's harmless (the counter is cosmetic there). For an encrypted
archive it is not: if an interior frame's CRC has failed, every block after
it would be opened under an assumed ordinal one (or more) lower than the true
one it was sealed under, turning one genuinely corrupt block into a cascade
of spurious AEAD failures on blocks that were otherwise intact. Rather than
attempt this and silently lose more data than the corruption alone caused,
`decompress_recover` refuses scan recovery outright for an encrypted archive
once the footer is confirmed unreadable, with a message directing the caller
to recompress from source. **Footer-driven recovery is unaffected**: the
footer's row-group index gives each block's *true* ordinal directly
(`footer.groups.iter().enumerate()`), regardless of what corruption exists in
the blocks themselves, so it remains the fully supported recovery path.

### Determinism

Two claims that don't conflict:

1. **Ciphertext is *not* byte-identical across separate `compress` runs of
   identical input and passphrase — deliberately.** `setup_encryption` draws
   a fresh random `salt`/`nonce_id` (one `getrandom` call) on every
   invocation. A fixed scheme would let anyone holding two archives compare
   their ciphertext bytes and learn whether the underlying content matched,
   without ever decrypting either — exactly the confidentiality leak a
   passphrase is supposed to prevent.
2. **Output *is* byte-identical across thread counts within one run** — the
   existing invariant (CLAUDE.md) holds unchanged. `salt`/`nonce_id` are drawn
   once, up front, before the per-block parallel loop runs; each block's
   nonce/AAD is a pure function of that fixed `nonce_id` and the block's
   ordinal, itself established (via `block_ranges`/`FooterIndex`) before any
   thread starts encoding. Sealing happens inside the same per-block `rayon`
   closure that already does codec work (`block::write_blocks`) — encryption
   introduces no new serialization point against the parallel-blocks
   invariant.

## Out of scope for this release

- **Recipient/public-key mode** and **raw keyfile mode** — passphrase-only for
  now (scope decision 2). Both are addable later as a new `crypt_version`/
  extension shape without disturbing archives already written.
- **Per-stream AEAD / preserving `fqxv.remote` column projection** on
  encrypted archives (scope decision 3) — `Index::read` refuses outright.
- **Encrypting the `FLAG_GLOBAL_REORDER` layout** (`--order
  any`/`shuffle`/`--max`) — rejected at `encode_reordered`; needs its own
  nonce/AAD/truncation design given that layout's different shape.
- **CLI-tunable Argon2id cost.** `KdfParams::DEFAULT` is used
  unconditionally; the chosen params still travel per-archive, so a future
  release can both expose a CLI override and raise the shipped default
  without breaking any archive already written.
- **Solving streaming-decode tail-truncation detection.** `fqxv decompress -`
  (piped, non-seekable) never reads the footer and so cannot detect a
  truncated tail regardless of encryption — a pre-existing limitation for
  plaintext archives too, not something this feature regresses or resolves.
- **Passphrase rotation without a full re-compress.** Changing the passphrase
  means re-running `compress --encrypt` with a new passphrase (fresh
  salt/nonce_id/key) — the same as re-encrypting from scratch.

## See also

- [`container.md`](container.md) — the base format's byte layout and
  evolution policy this feature extends.
- `crates/fqxv-crypt/src/lib.rs` — the crate implementing the scheme above.
- `crates/fqxv/src/container/format.rs` — the on-disk wiring
  (`EXT_TAG_ENCRYPTION`, `EncryptionHeader`, `setup_encryption`/
  `open_encryption`, `footer_crc_tail`, `Footer::body_bytes`).
