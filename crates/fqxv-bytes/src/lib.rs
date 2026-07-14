//! Shared byte-serialization primitives for the fqxv codec crates.
//!
//! These are the exact LEB128 varint and zig-zag encodings that the `fqxv-seq`,
//! `fqxv-reorder`, `fqxv-fqzcomp`, and `fqxv-tokenizer` crates use on disk. They
//! were previously copy-pasted, byte-identical, into each crate; this leaf crate
//! is the single source of truth. The byte layout must stay stable — every codec
//! stream that stores a varint reads it back with [`read_varint`].

/// Append `v` to `out` as an unsigned LEB128 varint.
///
/// The value is emitted 7 bits per byte, little-endian, with the high bit of
/// each byte set to signal that another byte follows (and clear on the last
/// byte). A zero value encodes as a single `0x00`.
pub fn write_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            break;
        }
        out.push(byte | 0x80);
    }
}

/// Decode an unsigned LEB128 varint from `src` starting at `*pos`, advancing
/// `*pos` past the bytes consumed.
///
/// Returns `None` on truncation (the encoding runs off the end of `src`) or on
/// an over-long (more than 64-bit) encoding. This is the inverse of
/// [`write_varint`].
pub fn read_varint(src: &[u8], pos: &mut usize) -> Option<u64> {
    let mut v = 0u64;
    let mut shift = 0u32;
    loop {
        let byte = *src.get(*pos)?;
        *pos += 1;
        v |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Some(v);
        }
        shift += 7;
        if shift >= 64 {
            return None;
        }
    }
}

/// Map a signed delta to an unsigned, varint-friendly value (zig-zag encoding).
///
/// Small-magnitude values of either sign map to small unsigned values, so they
/// stay cheap under [`write_varint`]. Inverse of [`unzigzag`].
pub fn zigzag(d: i64) -> u64 {
    ((d << 1) ^ (d >> 63)) as u64
}

/// Recover the signed delta from its zig-zag encoding. Inverse of [`zigzag`].
pub fn unzigzag(z: u64) -> i64 {
    ((z >> 1) as i64) ^ -((z & 1) as i64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn explicit_varint_vectors() {
        let cases: &[(u64, &[u8])] = &[
            (0, &[0x00]),
            (1, &[0x01]),
            (127, &[0x7f]),
            (128, &[0x80, 0x01]),
            (300, &[0xac, 0x02]),
            (16384, &[0x80, 0x80, 0x01]),
            (
                u64::MAX,
                &[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x01],
            ),
        ];
        for &(v, bytes) in cases {
            let mut out = Vec::new();
            write_varint(&mut out, v);
            assert_eq!(out, bytes, "encoding of {v}");
            let mut pos = 0;
            assert_eq!(read_varint(&out, &mut pos), Some(v));
            assert_eq!(pos, bytes.len());
        }
    }

    #[test]
    fn read_varint_truncated_is_none() {
        // Continuation bit set but no following byte.
        assert_eq!(read_varint(&[0x80], &mut 0), None);
        // Empty input.
        assert_eq!(read_varint(&[], &mut 0), None);
    }

    #[test]
    fn read_varint_overlong_is_none() {
        // Ten continuation bytes push shift to 70 (>= 64) before terminating.
        let overlong = [0x80u8; 10];
        assert_eq!(read_varint(&overlong, &mut 0), None);
    }

    #[test]
    fn explicit_zigzag_vectors() {
        let cases: &[(i64, u64)] = &[(0, 0), (-1, 1), (1, 2), (-2, 3), (2, 4)];
        for &(d, z) in cases {
            assert_eq!(zigzag(d), z, "zigzag({d})");
            assert_eq!(unzigzag(z), d, "unzigzag({z})");
        }
    }

    proptest! {
        #[test]
        fn varint_round_trips(v: u64) {
            let mut out = Vec::new();
            write_varint(&mut out, v);
            let mut pos = 0;
            prop_assert_eq!(read_varint(&out, &mut pos), Some(v));
            prop_assert_eq!(pos, out.len());
        }

        #[test]
        fn varint_round_trips_with_trailing_bytes(v: u64, tail: Vec<u8>) {
            // Decoding must stop exactly at the end of the varint, leaving the
            // trailing bytes untouched (this is how the codecs pack streams).
            let mut buf = Vec::new();
            write_varint(&mut buf, v);
            let varint_len = buf.len();
            buf.extend_from_slice(&tail);
            let mut pos = 0;
            prop_assert_eq!(read_varint(&buf, &mut pos), Some(v));
            prop_assert_eq!(pos, varint_len);
        }

        #[test]
        fn zigzag_round_trips(x: i64) {
            prop_assert_eq!(unzigzag(zigzag(x)), x);
        }
    }
}
