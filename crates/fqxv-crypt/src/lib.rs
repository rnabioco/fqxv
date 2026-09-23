//! Passphrase-based AEAD envelope for the `fqxv` container format.
//!
//! This crate is the one place in the workspace that takes real, non-clean-room
//! cryptographic dependencies (see `THIRD-PARTY-NOTICES.md`): [Argon2id] for
//! key derivation and [ChaCha20-Poly1305] for authenticated encryption, both
//! via RustCrypto crates rather than reimplemented from a spec. Every other
//! codec crate in the DAG is clean-room by policy; hand-rolling a cipher/KDF
//! is a different risk class (constant-time behavior, nonce-misuse resistance,
//! side channels) that policy was never meant to cover.
//!
//! [Argon2id]: https://www.rfc-editor.org/rfc/rfc9106
//! [ChaCha20-Poly1305]: https://www.rfc-editor.org/rfc/rfc8439
//!
//! # Scheme
//!
//! One passphrase-derived 256-bit key and one random 16-byte `nonce_id` cover
//! a whole archive. Every block, the footer, and the optional whole-file
//! reference frame are sealed independently under that one key, each with its
//! own deterministic nonce/associated-data derived from `nonce_id` and a
//! position index (see [`ArchiveCipher`]) — never a single stream cipher over
//! the whole file, so sealing/opening a block stays an independent, order-free
//! operation (the container's blocks are the unit of `rayon` parallelism).
//!
//! `nonce_id` and the Argon2id salt are freshly random on every archive, even
//! for identical input and passphrase: fixing them would let two archives'
//! ciphertext be compared to learn whether their contents matched, without
//! ever decrypting either — the exact leak a passphrase is meant to prevent.
//!
//! ```
//! use fqxv_crypt::{ArchiveCipher, KdfParams, Passphrase, derive_key, random_bytes32};
//!
//! let material = random_bytes32().unwrap();
//! let (salt, nonce_id) = (
//!     <[u8; 16]>::try_from(&material[..16]).unwrap(),
//!     <[u8; 16]>::try_from(&material[16..]).unwrap(),
//! );
//! let pass = Passphrase::from(b"correct horse battery staple".to_vec());
//! let key = derive_key(&pass, &salt, KdfParams::DEFAULT).unwrap();
//! let cipher = ArchiveCipher::new(key, nonce_id);
//!
//! let ciphertext = cipher.seal_block(0, b"a coded block payload");
//! assert_eq!(cipher.open_block(0, &ciphertext).unwrap(), b"a coded block payload");
//! // Wrong index: this ciphertext was sealed at index 0, not 1.
//! assert!(cipher.open_block(1, &ciphertext).is_err());
//! ```

use argon2::{Algorithm, Argon2, Params as Argon2Params, Version as Argon2Version};
use chacha20poly1305::{
    ChaCha20Poly1305, Key, Nonce, Tag,
    aead::{Aead, AeadInOut, KeyInit, Payload, inout::InOutBuf},
};
use std::fmt;
use thiserror::Error;
use zeroize::Zeroizing;

/// Bytes of the Argon2id salt stored per archive.
pub const SALT_LEN: usize = 16;
/// Bytes of the random per-archive nonce identifier stored per archive.
pub const NONCE_ID_LEN: usize = 16;
/// Bytes of the derived key ([`ArchiveKey`] / ChaCha20-Poly1305 key size).
pub const KEY_LEN: usize = 32;
/// Bytes of a ChaCha20-Poly1305 authentication tag.
pub const TAG_LEN: usize = 16;

/// Bytes of a ChaCha20-Poly1305 nonce. Deterministic per use (see
/// [`ArchiveCipher`]), never random, so it never needs to be stored.
const NONCE_LEN: usize = 12;

/// Domain-separation tag mixed into a block's associated data, so a block
/// ciphertext can never be mistaken for (or replayed as) a footer tag or
/// reference-frame ciphertext even though all three share one key.
const AAD_BLOCK: &[u8; 8] = b"FQXVBLK1";
/// See [`AAD_BLOCK`]; the footer's domain tag.
const AAD_FOOTER: &[u8; 8] = b"FQXVFTR1";
/// See [`AAD_BLOCK`]; the whole-file reference frame's domain tag.
const AAD_GREF: &[u8; 8] = b"FQXVGRF1";

/// Reserved block index for the footer's authentication tag. Real block
/// indices are bounded by the container's per-block read/byte budgets, so
/// this is unreachable by any actual block — [`ArchiveCipher::mac_footer`]/
/// [`ArchiveCipher::verify_footer`] debug-assert against a real block ever
/// colliding with it.
const FOOTER_INDEX: u64 = u64::MAX;
/// Reserved block index for the whole-file reference frame; see
/// [`FOOTER_INDEX`].
const GREF_INDEX: u64 = u64::MAX - 1;

/// Errors from key derivation and authenticated encryption/decryption.
///
/// Decryption failure is deliberately one variant ([`Error::Open`]): AEAD
/// cannot distinguish "wrong passphrase" from "tampered or corrupted
/// ciphertext," and reporting a more specific cause would be misleading, not
/// helpful.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// Argon2id key derivation failed (an invalid parameter combination, not
    /// a runtime condition — validated once, at the source of the params).
    #[error("key derivation failed: {0}")]
    Kdf(String),
    /// Decryption or footer-tag verification failed: wrong passphrase, or the
    /// data was corrupted or tampered with after it was sealed.
    #[error("decryption failed: wrong passphrase, or the data is corrupted or tampered")]
    Open,
    /// The operating system's random number generator was unavailable.
    #[error("system random number generator unavailable: {0}")]
    Rng(String),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Argon2id cost parameters. Stored on disk per archive (not baked into a
/// build), so a later `fqxv` release can raise [`KdfParams::DEFAULT`] without
/// breaking any archive already written — decoding always uses whatever
/// params the archive recorded, never a build-time constant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KdfParams {
    /// Memory cost in kibibytes.
    pub m_cost_kib: u32,
    /// Iteration count.
    pub t_cost: u32,
    /// Degree of parallelism.
    pub p_cost: u8,
}

impl KdfParams {
    /// RFC 9106 §4's "less memory available" recommended parameter set: 64
    /// MiB, 3 iterations, 4 lanes. Paid once per archive *open* (to derive the
    /// key), not per block, so it never sits in the per-block `rayon` hot
    /// path. The stronger 2 GiB profile is deliberately not the default here:
    /// `fqxv` runs everywhere from a laptop to a memory-constrained CI runner,
    /// and 2 GiB just to open one archive is a real availability cost for
    /// marginal extra resistance over 64 MiB, which is already far beyond a
    /// realistic cracking budget for this threat model.
    pub const DEFAULT: Self = Self {
        m_cost_kib: 65_536,
        t_cost: 3,
        p_cost: 4,
    };
}

/// A passphrase, held as raw bytes and zeroized on drop.
///
/// Deliberately has no encoding opinion: a caller decides how a passphrase
/// becomes bytes (UTF-8 for an interactive/env-var/string source, raw file
/// contents for a `--password-file`), and the same source must be used to
/// decrypt as was used to encrypt, since different encodings of "the same"
/// passphrase derive different keys.
pub struct Passphrase(Zeroizing<Vec<u8>>);

impl fmt::Debug for Passphrase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("Passphrase").field(&"<redacted>").finish()
    }
}

impl From<Vec<u8>> for Passphrase {
    fn from(bytes: Vec<u8>) -> Self {
        Passphrase(Zeroizing::new(bytes))
    }
}

impl From<&[u8]> for Passphrase {
    fn from(bytes: &[u8]) -> Self {
        Passphrase(Zeroizing::new(bytes.to_vec()))
    }
}

/// A derived 256-bit archive key, zeroized on drop.
pub struct ArchiveKey(Zeroizing<[u8; KEY_LEN]>);

impl fmt::Debug for ArchiveKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("ArchiveKey").field(&"<redacted>").finish()
    }
}

/// Fill `out` with cryptographically random bytes from the OS CSPRNG. Used to
/// draw a fresh Argon2id salt and archive `nonce_id` (16 bytes each) once per
/// `compress` invocation.
///
/// # Errors
/// Returns [`Error::Rng`] if the OS random source is unavailable.
pub fn random_bytes32() -> Result<[u8; 32]> {
    let mut buf = [0u8; 32];
    getrandom::fill(&mut buf).map_err(|e| Error::Rng(e.to_string()))?;
    Ok(buf)
}

/// Derive a 256-bit archive key from a passphrase and salt via Argon2id.
///
/// # Errors
/// Returns [`Error::Kdf`] if `kdf` describes an invalid Argon2 parameter
/// combination (out of range memory/time/parallelism costs).
pub fn derive_key(pass: &Passphrase, salt: &[u8; SALT_LEN], kdf: KdfParams) -> Result<ArchiveKey> {
    let params = Argon2Params::new(
        kdf.m_cost_kib,
        kdf.t_cost,
        u32::from(kdf.p_cost),
        Some(KEY_LEN),
    )
    .map_err(|e| Error::Kdf(e.to_string()))?;
    let argon2 = Argon2::new(Algorithm::Argon2id, Argon2Version::V0x13, params);
    let mut out = [0u8; KEY_LEN];
    argon2
        .hash_password_into(&pass.0, salt, &mut out)
        .map_err(|e| Error::Kdf(e.to_string()))?;
    Ok(ArchiveKey(Zeroizing::new(out)))
}

/// A key bound to one archive's random `nonce_id`, sealing/opening its
/// blocks, footer, and optional reference frame.
///
/// # Nonce and associated-data construction
///
/// Never exposed to callers — every `seal_*`/`open_*`/`mac_*`/`verify_*`
/// method builds its own nonce and AAD internally from `nonce_id` and a
/// position index:
///
/// - `nonce(idx) = nonce_id[0..4] ‖ LE64(idx)` — a block's `idx` is its
///   0-based footer ordinal: unique and strictly increasing within one
///   archive by construction, so `(key, nonce)` never repeats within an
///   archive (the property a stream cipher's nonce must have). A fully
///   random per-block nonce was not used because it would have to be stored
///   per block; this one is a pure function of position, free on disk.
/// - `aad(domain, idx) = domain ‖ nonce_id ‖ LE64(idx)` binds three things:
///   domain separation (a block can't be replayed as a footer or reference
///   frame), archive identity (`nonce_id`), and *position* — a decoder
///   recomputes the expected AAD from where a block sits in the stream, so
///   splicing, reordering, or duplicating blocks fails authentication
///   immediately, before any bytes are trusted.
///
/// Position-bound AAD does not catch dropping only the *trailing* blocks
/// behind a forged terminator (every surviving block still sits at its
/// correct ordinal) — [`ArchiveCipher::mac_footer`] closes that gap for any
/// reader that checks the footer with a passphrase.
pub struct ArchiveCipher {
    cipher: ChaCha20Poly1305,
    nonce_id: [u8; NONCE_ID_LEN],
}

impl fmt::Debug for ArchiveCipher {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ArchiveCipher")
            .field("nonce_id", &self.nonce_id)
            .finish_non_exhaustive()
    }
}

impl ArchiveCipher {
    /// Bind a derived key to one archive's `nonce_id`.
    #[must_use]
    pub fn new(key: ArchiveKey, nonce_id: [u8; NONCE_ID_LEN]) -> Self {
        let cipher = ChaCha20Poly1305::new(&Key::from(*key.0));
        ArchiveCipher { cipher, nonce_id }
    }

    fn nonce_for(&self, index: u64) -> Nonce {
        let mut n = [0u8; NONCE_LEN];
        n[..4].copy_from_slice(&self.nonce_id[..4]);
        n[4..].copy_from_slice(&index.to_le_bytes());
        Nonce::from(n)
    }

    fn aad_for(&self, domain: &[u8; 8], index: u64) -> [u8; 8 + NONCE_ID_LEN + 8] {
        let mut a = [0u8; 8 + NONCE_ID_LEN + 8];
        a[..8].copy_from_slice(domain);
        a[8..8 + NONCE_ID_LEN].copy_from_slice(&self.nonce_id);
        a[8 + NONCE_ID_LEN..].copy_from_slice(&index.to_le_bytes());
        a
    }

    /// Seal one block's already-coded plaintext payload, bound to its 0-based
    /// footer ordinal `block_index`. Returns ciphertext with the 16-byte
    /// authentication tag appended.
    ///
    /// # Panics
    /// Panics if `block_index` collides with a reserved index ([`FOOTER_INDEX`]/
    /// [`GREF_INDEX`]) — unreachable in practice, since real block counts are
    /// bounded far below `u64::MAX - 1` by the container's per-block budgets.
    #[must_use]
    pub fn seal_block(&self, block_index: u64, plaintext: &[u8]) -> Vec<u8> {
        debug_assert!(
            block_index < GREF_INDEX,
            "block index collides with a reserved index"
        );
        let nonce = self.nonce_for(block_index);
        let aad = self.aad_for(AAD_BLOCK, block_index);
        self.cipher
            .encrypt(
                &nonce,
                Payload {
                    msg: plaintext,
                    aad: &aad,
                },
            )
            .expect("ChaCha20Poly1305 seal cannot fail for in-bounds input")
    }

    /// Open one block's ciphertext at its 0-based footer ordinal
    /// `block_index`. Fails if the passphrase is wrong, the ciphertext was
    /// tampered with, or it was sealed at a different index (defeating
    /// splicing/reordering) — see [`ArchiveCipher`]'s AAD construction.
    ///
    /// # Errors
    /// Returns [`Error::Open`] on authentication failure.
    pub fn open_block(&self, block_index: u64, ciphertext: &[u8]) -> Result<Vec<u8>> {
        let nonce = self.nonce_for(block_index);
        let aad = self.aad_for(AAD_BLOCK, block_index);
        self.cipher
            .decrypt(
                &nonce,
                Payload {
                    msg: ciphertext,
                    aad: &aad,
                },
            )
            .map_err(|_| Error::Open)
    }

    /// Authenticate (not encrypt — the footer must stay directly readable
    /// without a key, for `inspect`/size reporting) the plaintext footer
    /// body. Returns a 16-byte tag binding the *entire* footer body,
    /// including the read/block counts a tail-truncation attacker would
    /// otherwise leave unchecked (see [`ArchiveCipher`]'s docs).
    #[must_use]
    pub fn mac_footer(&self, footer_body: &[u8]) -> [u8; TAG_LEN] {
        let nonce = self.nonce_for(FOOTER_INDEX);
        let aad = self.footer_aad(footer_body);
        let mut empty: [u8; 0] = [];
        let tag: Tag = self
            .cipher
            .encrypt_inout_detached(&nonce, &aad, InOutBuf::from(&mut empty[..]))
            .expect("ChaCha20Poly1305 detached seal cannot fail for in-bounds input");
        tag.into()
    }

    /// Verify a tag produced by [`ArchiveCipher::mac_footer`] over the same
    /// footer body.
    ///
    /// # Errors
    /// Returns [`Error::Open`] on authentication failure.
    pub fn verify_footer(&self, footer_body: &[u8], tag: &[u8; TAG_LEN]) -> Result<()> {
        let nonce = self.nonce_for(FOOTER_INDEX);
        let aad = self.footer_aad(footer_body);
        let mut empty: [u8; 0] = [];
        let tag_arr = Tag::from(*tag);
        self.cipher
            .decrypt_inout_detached(&nonce, &aad, InOutBuf::from(&mut empty[..]), &tag_arr)
            .map_err(|_| Error::Open)
    }

    /// `mac_footer`/`verify_footer`'s associated data: domain tag, archive
    /// identity, then the footer body itself (there is only one footer per
    /// archive, so no position index is needed the way blocks need one).
    fn footer_aad(&self, footer_body: &[u8]) -> Vec<u8> {
        let mut a = Vec::with_capacity(8 + NONCE_ID_LEN + footer_body.len());
        a.extend_from_slice(AAD_FOOTER);
        a.extend_from_slice(&self.nonce_id);
        a.extend_from_slice(footer_body);
        a
    }

    /// Seal the whole-file long-read reference frame (`FLAG_GLOBAL_REFERENCE`)
    /// — it carries real assembled sequence, so it must be sealed whenever
    /// present in an encrypted archive.
    #[must_use]
    pub fn seal_reference_frame(&self, plaintext: &[u8]) -> Vec<u8> {
        let nonce = self.nonce_for(GREF_INDEX);
        let aad = self.aad_for(AAD_GREF, GREF_INDEX);
        self.cipher
            .encrypt(
                &nonce,
                Payload {
                    msg: plaintext,
                    aad: &aad,
                },
            )
            .expect("ChaCha20Poly1305 seal cannot fail for in-bounds input")
    }

    /// Open the whole-file long-read reference frame sealed by
    /// [`ArchiveCipher::seal_reference_frame`].
    ///
    /// # Errors
    /// Returns [`Error::Open`] on authentication failure.
    pub fn open_reference_frame(&self, ciphertext: &[u8]) -> Result<Vec<u8>> {
        let nonce = self.nonce_for(GREF_INDEX);
        let aad = self.aad_for(AAD_GREF, GREF_INDEX);
        self.cipher
            .decrypt(
                &nonce,
                Payload {
                    msg: ciphertext,
                    aad: &aad,
                },
            )
            .map_err(|_| Error::Open)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn cipher_for(
        passphrase: &[u8],
        salt: [u8; SALT_LEN],
        nonce_id: [u8; NONCE_ID_LEN],
    ) -> ArchiveCipher {
        let pass = Passphrase::from(passphrase.to_vec());
        let key = derive_key(&pass, &salt, KdfParams::DEFAULT).unwrap();
        ArchiveCipher::new(key, nonce_id)
    }

    /// One `ArchiveCipher`, derived once and reused across every proptest case
    /// below. Argon2id at [`KdfParams::DEFAULT`] costs real wall time by design
    /// (§ the crate's own docs on why); a proptest run generates on the order of
    /// hundreds of cases; deriving fresh per case turned two property tests into
    /// a multi-minute KDF benchmark instead of an AEAD correctness check, which
    /// is the property actually under test here — key derivation itself has its
    /// own dedicated, non-proptest tests above.
    fn shared_cipher() -> &'static ArchiveCipher {
        static CIPHER: std::sync::OnceLock<ArchiveCipher> = std::sync::OnceLock::new();
        CIPHER.get_or_init(|| cipher_for(b"proptest passphrase", [3; SALT_LEN], [4; NONCE_ID_LEN]))
    }

    #[test]
    fn derive_key_is_deterministic() {
        let pass = Passphrase::from(b"hunter2".to_vec());
        let salt = [7u8; SALT_LEN];
        let k1 = derive_key(&pass, &salt, KdfParams::DEFAULT).unwrap();
        let k2 = derive_key(&pass, &salt, KdfParams::DEFAULT).unwrap();
        assert_eq!(&*k1.0, &*k2.0);
    }

    #[test]
    fn derive_key_is_sensitive_to_salt() {
        let pass = Passphrase::from(b"hunter2".to_vec());
        let k1 = derive_key(&pass, &[1u8; SALT_LEN], KdfParams::DEFAULT).unwrap();
        let k2 = derive_key(&pass, &[2u8; SALT_LEN], KdfParams::DEFAULT).unwrap();
        assert_ne!(&*k1.0, &*k2.0);
    }

    #[test]
    fn seal_open_block_round_trips() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let ct = cipher.seal_block(3, b"hello, block");
        assert_eq!(cipher.open_block(3, &ct).unwrap(), b"hello, block");
    }

    #[test]
    fn seal_open_block_round_trips_empty_plaintext() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let ct = cipher.seal_block(0, b"");
        assert_eq!(cipher.open_block(0, &ct).unwrap(), b"");
    }

    #[test]
    fn open_block_rejects_wrong_key() {
        let a = cipher_for(b"pw-a", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let b = cipher_for(b"pw-b", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let ct = a.seal_block(0, b"secret");
        assert!(matches!(b.open_block(0, &ct), Err(Error::Open)));
    }

    #[test]
    fn open_block_rejects_wrong_nonce_id() {
        let a = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let b = cipher_for(b"pw", [1; SALT_LEN], [9; NONCE_ID_LEN]);
        let ct = a.seal_block(0, b"secret");
        assert!(matches!(b.open_block(0, &ct), Err(Error::Open)));
    }

    #[test]
    fn open_block_rejects_wrong_index() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let ct = cipher.seal_block(3, b"secret");
        assert!(matches!(cipher.open_block(5, &ct), Err(Error::Open)));
    }

    #[test]
    fn open_block_rejects_tampered_ciphertext() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let mut ct = cipher.seal_block(0, b"secret payload");
        let last = ct.len() - 1;
        ct[last] ^= 0xff;
        assert!(matches!(cipher.open_block(0, &ct), Err(Error::Open)));
    }

    #[test]
    fn mac_verify_footer_round_trips() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let body = b"a footer body worth authenticating";
        let tag = cipher.mac_footer(body);
        assert!(cipher.verify_footer(body, &tag).is_ok());
    }

    #[test]
    fn verify_footer_rejects_tampered_body() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let mut body = b"a footer body worth authenticating".to_vec();
        let tag = cipher.mac_footer(&body);
        body[0] ^= 0xff;
        assert!(matches!(
            cipher.verify_footer(&body, &tag),
            Err(Error::Open)
        ));
    }

    #[test]
    fn seal_open_reference_frame_round_trips() {
        let cipher = cipher_for(b"pw", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let ct = cipher.seal_reference_frame(b"a whole-file reference");
        assert_eq!(
            cipher.open_reference_frame(&ct).unwrap(),
            b"a whole-file reference"
        );
    }

    #[test]
    fn open_reference_frame_rejects_wrong_key() {
        let a = cipher_for(b"pw-a", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let b = cipher_for(b"pw-b", [1; SALT_LEN], [2; NONCE_ID_LEN]);
        let ct = a.seal_reference_frame(b"reference bytes");
        assert!(matches!(b.open_reference_frame(&ct), Err(Error::Open)));
    }

    #[test]
    fn random_bytes32_is_not_constant() {
        // Not a statistical test — just a sanity check that two draws differ,
        // which would fail if this were ever wired to a fixed buffer by mistake.
        let a = random_bytes32().unwrap();
        let b = random_bytes32().unwrap();
        assert_ne!(a, b);
    }

    proptest! {
        #[test]
        fn seal_open_round_trips_arbitrary_plaintext(
            index in any::<u64>().prop_filter("below reserved indices", |i| *i < GREF_INDEX),
            plaintext in proptest::collection::vec(any::<u8>(), 0..2048),
        ) {
            let cipher = shared_cipher();
            let ct = cipher.seal_block(index, &plaintext);
            prop_assert_eq!(cipher.open_block(index, &ct).unwrap(), plaintext);
        }

        #[test]
        fn tampering_any_ciphertext_byte_is_caught(
            plaintext in proptest::collection::vec(any::<u8>(), 1..256),
            flip_at in 0usize..256,
        ) {
            let cipher = shared_cipher();
            let mut ct = cipher.seal_block(0, &plaintext);
            let idx = flip_at % ct.len();
            ct[idx] ^= 0x01;
            prop_assert!(cipher.open_block(0, &ct).is_err());
        }
    }
}
