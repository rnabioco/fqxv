//! CRC-32C (Castagnoli) integrity checksums for the container format.
//!
//! Every coded payload carries a CRC-32C so a single flipped bit is detected
//! and localized to one block rather than silently decoded into wrong bases or
//! quality scores. CRC-32C is the same checksum BGZF/BAM and CRAM use, and it
//! rides the SSE4.2 `crc32` instruction on hardware (a software table is used
//! here as the portable reference — a per-block CRC is negligible next to the
//! entropy coding either way).
//!
//! Reflected polynomial `0x82F63B78`, init/xor-out `0xFFFF_FFFF` — the standard
//! CRC-32C parameters (check value `0xE306_9283` over `b"123456789"`).

use std::io::{self, Write};

/// Reflected CRC-32C generator polynomial (Castagnoli).
const POLY: u32 = 0x82F6_3B78;

/// Byte-wise lookup table, generated at compile time from [`POLY`].
const TABLE: [u32; 256] = {
    let mut table = [0u32; 256];
    let mut i = 0usize;
    while i < 256 {
        let mut crc = i as u32;
        let mut bit = 0;
        while bit < 8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ POLY
            } else {
                crc >> 1
            };
            bit += 1;
        }
        table[i] = crc;
        i += 1;
    }
    table
};

/// Incremental CRC-32C accumulator. Holds the running (un-finalized) state so it
/// can be fed in pieces; [`Hasher::finalize`] applies the final xor.
#[derive(Debug, Clone)]
pub(crate) struct Hasher {
    state: u32,
}

impl Hasher {
    /// A fresh hasher over the empty input.
    pub(crate) fn new() -> Self {
        Hasher { state: !0 }
    }

    /// Fold `bytes` into the running checksum.
    pub(crate) fn update(&mut self, bytes: &[u8]) {
        let mut state = self.state;
        for &b in bytes {
            state = TABLE[((state ^ b as u32) & 0xFF) as usize] ^ (state >> 8);
        }
        self.state = state;
    }

    /// The CRC-32C of everything fed so far.
    pub(crate) fn finalize(&self) -> u32 {
        self.state ^ !0
    }
}

impl Default for Hasher {
    fn default() -> Self {
        Self::new()
    }
}

/// CRC-32C of a single contiguous buffer.
pub(crate) fn crc32c(bytes: &[u8]) -> u32 {
    let mut h = Hasher::new();
    h.update(bytes);
    h.finalize()
}

/// Number of bits in a CRC-32 register — the dimension of the GF(2) operator
/// matrices used by [`crc32c_combine`].
const GF2_DIM: usize = 32;

/// Multiply the GF(2) column-vector `vec` by the bit-matrix `mat` (each `mat[i]`
/// is a column of the operator). Used to apply a precomputed "append N zero
/// bytes" operator to a finalized CRC.
fn gf2_matrix_times(mat: &[u32; GF2_DIM], mut vec: u32) -> u32 {
    let mut sum = 0u32;
    let mut i = 0;
    while vec != 0 {
        if vec & 1 != 0 {
            sum ^= mat[i];
        }
        vec >>= 1;
        i += 1;
    }
    sum
}

/// Square the operator `mat` into `square` (compose it with itself), doubling the
/// number of appended zero bits it represents.
fn gf2_matrix_square(square: &mut [u32; GF2_DIM], mat: &[u32; GF2_DIM]) {
    for n in 0..GF2_DIM {
        square[n] = gf2_matrix_times(mat, mat[n]);
    }
}

/// Combine two CRC-32C checksums: given `crc_a = crc32c(A)`, `crc_b = crc32c(B)`,
/// and `len_b = B.len()`, return `crc32c(A ++ B)` without rescanning either half.
///
/// This is zlib's `crc32_combine` operator specialized to the Castagnoli
/// reflected polynomial. It lets a whole-buffer CRC be reassembled from
/// independently-checksummed chunks in any grouping, which is what the parallel
/// [`verify`](crate::verify) path uses to keep its check byte-for-byte identical
/// to the serial single-pass CRC. Operates on finalized CRC values.
pub(crate) fn crc32c_combine(crc_a: u32, crc_b: u32, len_b: u64) -> u32 {
    if len_b == 0 {
        return crc_a;
    }
    // `odd` holds the operator for an odd power-of-two run of zero bits, `even`
    // the next even one; they leapfrog by squaring as we walk the bits of len_b.
    let mut odd = [0u32; GF2_DIM];
    let mut even = [0u32; GF2_DIM];

    // Operator for a single zero bit: the reflected generator in column 0, then
    // the shifted identity in the remaining columns.
    odd[0] = POLY;
    let mut row = 1u32;
    for col in odd.iter_mut().skip(1) {
        *col = row;
        row <<= 1;
    }

    gf2_matrix_square(&mut even, &odd); // two zero bits
    gf2_matrix_square(&mut odd, &even); // four zero bits

    let mut crc = crc_a;
    let mut len = len_b;
    loop {
        // `even` now represents twice `odd`'s zeros; on the first pass that is one
        // whole zero *byte*. Apply it when this bit of the byte count is set.
        gf2_matrix_square(&mut even, &odd);
        if len & 1 != 0 {
            crc = gf2_matrix_times(&even, crc);
        }
        len >>= 1;
        if len == 0 {
            break;
        }
        gf2_matrix_square(&mut odd, &even);
        if len & 1 != 0 {
            crc = gf2_matrix_times(&odd, crc);
        }
        len >>= 1;
        if len == 0 {
            break;
        }
    }
    crc ^ crc_b
}

/// A `Write` adapter that tees every byte through a [`Hasher`], so the caller can
/// read back a whole-stream CRC-32C without a second pass over the output.
pub(crate) struct CrcWriter<W> {
    inner: W,
    hasher: Hasher,
}

impl<W: Write> CrcWriter<W> {
    /// Wrap `inner`, checksumming everything written to it.
    pub(crate) fn new(inner: W) -> Self {
        CrcWriter {
            inner,
            hasher: Hasher::new(),
        }
    }

    /// CRC-32C of every byte written through this adapter so far.
    pub(crate) fn crc(&self) -> u32 {
        self.hasher.finalize()
    }
}

impl<W: Write> Write for CrcWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        // Only the bytes the sink actually accepted are part of the stream.
        self.hasher.update(&buf[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_answers() {
        // Standard CRC-32C check vectors.
        assert_eq!(crc32c(b""), 0x0000_0000);
        assert_eq!(crc32c(b"123456789"), 0xE306_9283);
    }

    #[test]
    fn incremental_matches_oneshot() {
        let data = b"the quick brown fox jumps over the lazy dog";
        let mut h = Hasher::new();
        h.update(&data[..10]);
        h.update(&data[10..]);
        assert_eq!(h.finalize(), crc32c(data));
    }

    #[test]
    fn crc_writer_tees_and_forwards() {
        let mut sink = Vec::new();
        {
            let mut w = CrcWriter::new(&mut sink);
            w.write_all(b"123456789").unwrap();
            assert_eq!(w.crc(), 0xE306_9283);
        }
        assert_eq!(sink, b"123456789");
    }

    #[test]
    fn detects_single_bit_flip() {
        let a = crc32c(b"ACGTACGTACGT");
        let b = crc32c(b"ACGTACGTACGA");
        assert_ne!(a, b);
    }

    #[test]
    fn combine_reassembles_split_buffer() {
        let data = b"the quick brown fox jumps over the lazy dog";
        for split in 0..=data.len() {
            let (a, b) = data.split_at(split);
            let combined = crc32c_combine(crc32c(a), crc32c(b), b.len() as u64);
            assert_eq!(combined, crc32c(data), "split at {split}");
        }
    }

    #[test]
    fn combine_with_empty_tail_is_identity() {
        let a = crc32c(b"123456789");
        assert_eq!(crc32c_combine(a, crc32c(b""), 0), a);
    }

    proptest::proptest! {
        // Combining chunk CRCs must equal the CRC of the concatenation, for any
        // split point — this is what lets `verify` checksum chunks in parallel.
        #[test]
        fn combine_matches_concatenation(a: Vec<u8>, b: Vec<u8>) {
            let mut cat = a.clone();
            cat.extend_from_slice(&b);
            let combined = crc32c_combine(crc32c(&a), crc32c(&b), b.len() as u64);
            proptest::prop_assert_eq!(combined, crc32c(&cat));
        }
    }
}
