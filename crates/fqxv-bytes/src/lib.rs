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

use core::marker::PhantomData;

/// Error constructors a codec crate supplies so it can use [`Reader`].
///
/// [`Reader`] is generic over the caller's error type so each codec crate keeps
/// its own `Error` (and `Result` alias) while sharing the byte-level bounds
/// checks. A crate implements this trait once for its `Error`, aliases
/// `Reader<'a, MyError>`, and every `?` in its decoder then yields `MyError` —
/// with no per-crate reader struct to maintain.
pub trait ReaderError {
    /// The stream ended before a requested byte or slice was available.
    fn truncated() -> Self;
    /// A varint ran past 64 bits (an over-long encoding).
    fn bad_varint() -> Self;
    /// A length prefix was too large to allocate for.
    fn oversized() -> Self;
}

/// A forward-only, bounds-checked cursor over an in-memory byte slice.
///
/// This is the single source of truth for the truncation and overflow guards
/// that the codec crates' decoders used to copy-paste into private
/// `Cursor`/`ByteReader` structs. Every read advances the internal position;
/// failures map to the caller's error type through [`ReaderError`].
pub struct Reader<'a, E> {
    buf: &'a [u8],
    pos: usize,
    _err: PhantomData<fn() -> E>,
}

impl<E> core::fmt::Debug for Reader<'_, E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Reader")
            .field("pos", &self.pos)
            .field("len", &self.buf.len())
            .finish()
    }
}

impl<'a, E: ReaderError> Reader<'a, E> {
    /// Start a reader at the beginning of `buf`.
    pub fn new(buf: &'a [u8]) -> Self {
        Reader {
            buf,
            pos: 0,
            _err: PhantomData,
        }
    }

    /// Read one byte and advance.
    pub fn u8(&mut self) -> Result<u8, E> {
        let b = *self.buf.get(self.pos).ok_or_else(E::truncated)?;
        self.pos += 1;
        Ok(b)
    }

    /// Read a little-endian `u32` and advance 4 bytes.
    pub fn u32(&mut self) -> Result<u32, E> {
        let s = self.take(4)?;
        Ok(u32::from_le_bytes(s.try_into().unwrap()))
    }

    /// Read a little-endian `u64` and advance 8 bytes.
    pub fn u64(&mut self) -> Result<u64, E> {
        let s = self.take(8)?;
        Ok(u64::from_le_bytes(s.try_into().unwrap()))
    }

    /// Read an unsigned LEB128 varint (see [`read_varint`]) and advance.
    pub fn varint(&mut self) -> Result<u64, E> {
        read_varint(self.buf, &mut self.pos).ok_or_else(E::bad_varint)
    }

    /// Borrow the next `n` bytes and advance past them.
    ///
    /// `n` is treated as untrusted: the position add is checked so a hostile
    /// length yields [`ReaderError::truncated`] rather than wrapping.
    pub fn take(&mut self, n: usize) -> Result<&'a [u8], E> {
        let end = self.pos.checked_add(n).ok_or_else(E::truncated)?;
        let s = self.buf.get(self.pos..end).ok_or_else(E::truncated)?;
        self.pos = end;
        Ok(s)
    }

    /// Read a varint length prefix, then borrow that many bytes.
    pub fn take_stream(&mut self) -> Result<&'a [u8], E> {
        let n = self.varint()? as usize;
        self.take(n)
    }

    /// Read a little-endian `u32` length prefix, then borrow that many bytes.
    pub fn slice_u32(&mut self) -> Result<&'a [u8], E> {
        let n = self.u32()? as usize;
        self.take(n)
    }

    /// Borrow all not-yet-consumed bytes without advancing.
    pub fn rest(&self) -> &'a [u8] {
        &self.buf[self.pos..]
    }
}

/// Serialize a read-length array with a fixed-length fast path.
///
/// Writes the count, a flag byte for whether all lengths are equal, then either
/// the single shared length or every length — each as a varint. Inverse of
/// [`read_lens`]. This is the exact on-disk encoding shared by the `fqxv-seq`
/// and `fqxv-fqzcomp` decoders; the byte layout must stay stable.
pub fn write_lens(out: &mut Vec<u8>, lens: &[u32]) {
    write_varint(out, lens.len() as u64);
    let fixed = lens.first().is_some_and(|&f| lens.iter().all(|&l| l == f));
    out.push(u8::from(fixed));
    if fixed {
        if let Some(&f) = lens.first() {
            write_varint(out, u64::from(f));
        }
    } else {
        for &l in lens {
            write_varint(out, u64::from(l));
        }
    }
}

/// Deserialize a read-length array written by [`write_lens`].
pub fn read_lens<E: ReaderError>(r: &mut Reader<'_, E>) -> Result<Vec<u32>, E> {
    let n = r.varint()? as usize;
    let fixed = r.u8()? != 0;
    let mut lens = Vec::new();
    if fixed {
        // The fixed path allocates all `n` entries up front regardless of how
        // many input bytes remain, so an untrusted `n` must not abort the
        // process on a hostile allocation — turn it into a clean error.
        if n > 0 {
            let f = r.varint()? as u32;
            lens.try_reserve_exact(n).map_err(|_| E::oversized())?;
            lens.resize(n, f);
        }
    } else {
        // Self-limiting: each length is a varint consuming >= 1 input byte, so
        // `n` is bounded by the remaining input. Cap the speculative reserve.
        lens.try_reserve(n.min(1 << 20))
            .map_err(|_| E::oversized())?;
        for _ in 0..n {
            lens.push(r.varint()? as u32);
        }
    }
    Ok(lens)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Minimal error type so the tests can exercise the generic [`Reader`].
    #[derive(Debug, PartialEq, Eq)]
    enum TestErr {
        Truncated,
        BadVarint,
        Oversized,
    }

    impl ReaderError for TestErr {
        fn truncated() -> Self {
            TestErr::Truncated
        }
        fn bad_varint() -> Self {
            TestErr::BadVarint
        }
        fn oversized() -> Self {
            TestErr::Oversized
        }
    }

    #[test]
    fn reader_reads_and_advances() {
        let mut buf = Vec::new();
        buf.push(0xAB); // u8
        write_varint(&mut buf, 300); // varint
        buf.extend_from_slice(&7u32.to_le_bytes()); // u32
        buf.extend_from_slice(&9u64.to_le_bytes()); // u64
        buf.extend_from_slice(b"hi"); // take(2)

        let mut r = Reader::<TestErr>::new(&buf);
        assert_eq!(r.u8(), Ok(0xAB));
        assert_eq!(r.varint(), Ok(300));
        assert_eq!(r.u32(), Ok(7));
        assert_eq!(r.u64(), Ok(9));
        assert_eq!(r.take(2), Ok(&b"hi"[..]));
        assert_eq!(r.rest(), b"");
    }

    #[test]
    fn reader_truncation_is_an_error() {
        let mut r = Reader::<TestErr>::new(&[0x01]);
        assert_eq!(r.u8(), Ok(1));
        assert_eq!(r.u8(), Err(TestErr::Truncated));
        // A length that overflows `usize` must error, not wrap.
        let mut r = Reader::<TestErr>::new(&[]);
        assert_eq!(r.take(usize::MAX), Err(TestErr::Truncated));
    }

    #[test]
    fn slice_prefixes_round_trip() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&3u32.to_le_bytes());
        buf.extend_from_slice(b"abc");
        write_varint(&mut buf, 2);
        buf.extend_from_slice(b"de");

        let mut r = Reader::<TestErr>::new(&buf);
        assert_eq!(r.slice_u32(), Ok(&b"abc"[..]));
        assert_eq!(r.take_stream(), Ok(&b"de"[..]));
    }

    #[test]
    fn lens_round_trip_fixed_and_variable() {
        for lens in [vec![], vec![100, 100, 100], vec![5, 9, 5, 2]] {
            let mut out = Vec::new();
            write_lens(&mut out, &lens);
            let mut r = Reader::<TestErr>::new(&out);
            assert_eq!(read_lens(&mut r), Ok(lens.clone()));
            assert_eq!(r.rest(), b"", "consumed all bytes for {lens:?}");
        }
    }

    proptest! {
        #[test]
        fn lens_round_trip_prop(lens: Vec<u32>) {
            let mut out = Vec::new();
            write_lens(&mut out, &lens);
            let mut r = Reader::<TestErr>::new(&out);
            prop_assert_eq!(read_lens(&mut r), Ok(lens));
        }
    }

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
