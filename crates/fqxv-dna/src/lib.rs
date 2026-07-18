//! Shared nucleotide primitives for the fqxv sequence codecs.
//!
//! The 2-bit `ACGT` encodings and the reverse-complement mapping were previously
//! copy-pasted — subtly, and not always identically — into `fqxv-seq`,
//! `fqxv-reorder`, and `fqxv-lroverlap`. This leaf crate is the single source of
//! truth for them, so the variants can no longer drift apart silently.
//!
//! Two distinctions the codecs genuinely rely on are made explicit here rather
//! than hidden inside each caller:
//!
//! * **Case sensitivity of the 2-bit code.** [`code_strict`] treats *only*
//!   uppercase `A/C/G/T` as bases (lowercase becomes the non-ACGT sentinel and is
//!   exception-coded byte-for-byte), which is what a byte-lossless packer needs.
//!   [`code_fold`] additionally folds lowercase onto the same code, which is what
//!   a clustering/minimizer front-end wants (it only needs k-mer identity, never
//!   byte-exactness).
//! * **Case handling of reverse complement.** [`revcomp`] complements `A/C/G/T`
//!   in *both* cases; [`revcomp_acgt`] complements only the uppercase four and
//!   passes every other byte (including lowercase) through untouched. Callers
//!   that pair `revcomp` with case-folded coding use the former; callers that
//!   keep lowercase as verbatim exceptions use the latter.

/// The four 2-bit symbols in code order: `0 -> A`, `1 -> C`, `2 -> G`, `3 -> T`.
///
/// Index it directly with a symbol known to be `0..4`; use [`base_of_sym`] when
/// the input may be an out-of-range placeholder.
pub const SYM2BASE: [u8; 4] = *b"ACGT";

/// The non-ACGT sentinel returned by [`code_strict`] / [`code_fold`] and stored
/// in [`BASE_LUT`] for every byte that is not a canonical base.
pub const NON_ACGT: u8 = 255;

/// Lookup table: uppercase `A/C/G/T` map to `0..4`, every other byte (lowercase,
/// `N`, IUPAC ambiguity codes, anything) maps to [`NON_ACGT`]. This is the
/// case-sensitive, byte-lossless mapping — see [`code_strict`].
pub const BASE_LUT: [u8; 256] = base_lut();

const fn base_lut() -> [u8; 256] {
    let mut t = [NON_ACGT; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
}

/// Case-*sensitive* 2-bit code: uppercase `A/C/G/T` become `0..4`; every other
/// byte becomes [`NON_ACGT`]. Lowercase is deliberately *not* folded, so a packer
/// can round-trip it verbatim through an exception list. Equivalent to indexing
/// [`BASE_LUT`].
#[inline]
#[must_use]
pub fn code_strict(b: u8) -> u8 {
    BASE_LUT[b as usize]
}

/// Case-*insensitive* 2-bit code: `A/a`, `C/c`, `G/g`, `T/t` all fold onto
/// `0..4`; every other byte becomes [`NON_ACGT`]. Use this only where k-mer
/// identity matters and the original bytes are recovered elsewhere (clustering,
/// minimizers), never where the exact byte must be restored.
#[inline]
#[must_use]
pub fn code_fold(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => NON_ACGT,
    }
}

/// Inverse of the 2-bit code for `sym` in `0..4`; any other value maps to `A`
/// (the placeholder an exception list is expected to overwrite).
#[inline]
#[must_use]
pub fn base_of_sym(sym: u8) -> u8 {
    match sym {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'A',
    }
}

/// True iff `b` is one of the uppercase canonical bases `A/C/G/T`.
#[inline]
#[must_use]
pub fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Reverse complement, complementing `A/C/G/T` in *both* cases (`a -> t`, …) and
/// passing every other byte (including `N`) through unchanged, reversed.
#[inline]
#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; seq.len()];
    revcomp_into(seq, &mut out);
    out
}

/// Reverse-complement `seq` into `dst`, allocating nothing. Same mapping as
/// [`revcomp`] — which delegates here, so the two cannot drift apart.
///
/// # Panics
/// If `dst.len() != seq.len()`.
#[inline]
pub fn revcomp_into(seq: &[u8], dst: &mut [u8]) {
    assert_eq!(dst.len(), seq.len(), "revcomp_into: length mismatch");
    for (d, &b) in dst.iter_mut().zip(seq.iter().rev()) {
        *d = match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        };
    }
}

/// Reverse complement of an uppercase-only sequence: complements `A/C/G/T` and
/// passes every other byte — *including lowercase* — through unchanged, reversed.
///
/// This is the mapping for codecs that keep non-uppercase bytes as verbatim
/// exceptions; unlike [`revcomp`] it does not touch lowercase `a/c/g/t`.
#[inline]
#[must_use]
pub fn revcomp_acgt(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn sym_roundtrip() {
        for sym in 0u8..4 {
            assert_eq!(code_strict(SYM2BASE[sym as usize]), sym);
            assert_eq!(base_of_sym(sym), SYM2BASE[sym as usize]);
        }
        assert_eq!(base_of_sym(4), b'A');
        assert_eq!(base_of_sym(255), b'A');
    }

    #[test]
    fn code_strict_is_case_sensitive() {
        assert_eq!(code_strict(b'A'), 0);
        assert_eq!(code_strict(b'a'), NON_ACGT);
        assert_eq!(code_strict(b'N'), NON_ACGT);
    }

    #[test]
    fn code_fold_folds_case() {
        assert_eq!(code_fold(b'a'), 0);
        assert_eq!(code_fold(b'T'), 3);
        assert_eq!(code_fold(b't'), 3);
        assert_eq!(code_fold(b'N'), NON_ACGT);
    }

    #[test]
    fn base_lut_matches_code_strict() {
        for b in 0u8..=255 {
            assert_eq!(BASE_LUT[b as usize], code_strict(b));
        }
    }

    proptest! {
        #[test]
        fn revcomp_is_involutive(bases in proptest::collection::vec(0u8..=255, 0..256)) {
            prop_assert_eq!(revcomp(&revcomp(&bases)), bases.clone());
        }

        #[test]
        fn revcomp_acgt_is_involutive(bases in proptest::collection::vec(0u8..=255, 0..256)) {
            prop_assert_eq!(revcomp_acgt(&revcomp_acgt(&bases)), bases.clone());
        }

        #[test]
        fn revcomp_into_matches_revcomp(bases in proptest::collection::vec(0u8..=255, 0..256)) {
            let mut dst = vec![0u8; bases.len()];
            revcomp_into(&bases, &mut dst);
            prop_assert_eq!(dst, revcomp(&bases));
        }

        #[test]
        fn revcomp_variants_agree_on_uppercase(
            bases in proptest::collection::vec(prop::sample::select(&b"ACGTN"[..]), 0..256)
        ) {
            // With no lowercase present the two revcomp variants must coincide.
            prop_assert_eq!(revcomp(&bases), revcomp_acgt(&bases));
        }
    }
}
