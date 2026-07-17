//! SPRING-faithful reference coder: **2-bit-pack the ACGT consensus, then run a
//! byte-domain compressor on the packed bytes.**
//!
//! This is exactly what SPRING does (`pack_compress_seq` in its `encoder.cpp`:
//! 4 bases → 1 byte, then BSC on the packed stream), and on real data it beats
//! fqxv's order-k context model on the raw bases by ~7% — because the packing is
//! a hard 2 bits/base floor and a byte-domain LZ/BWT then captures the long-range
//! near-duplicate-contig repeats the context model cannot see. The order-k model
//! and a BWT/LZMA run directly on the *raw* 4-symbol bases both underperform this;
//! the 2-bit packing is the missing front-end.
//!
//! The consensus is pure ACGT by construction (plurality vote), with only a
//! handful of `N`/other bytes in practice; those are packed as `A` and restored
//! from a tiny exception list, so the codec stays byte-exact. Packing shrinks the
//! input 4× (≈22 MB for an 89 M-base reference), so the whole reference fits one
//! LZ window / BWT block — full long-range reach, and fast.
//!
//! Entropy backend: the clean-room LZMA core ([`super::reflzma`]), which on the
//! packed bytes is in its natural byte-oriented regime (LZ/BSC on packed win here;
//! see the measurements in the reference-coder notes). Gated never-worse by the
//! container.

use fqxv_bytes::{read_varint, write_varint};

use crate::{reflzma, Error, Result};

const SYM2BASE: [u8; 4] = *b"ACGT";

/// Strict, case-SENSITIVE base code: only uppercase `A/C/G/T` pack to `0..4`;
/// everything else (lowercase, `N`, IUPAC) is `255` and exception-coded, so the
/// unpack restores the exact byte. (The shared [`code`](crate::code) is
/// case-insensitive and would upper-case lowercase bases — not byte-exact.)
#[inline]
fn pack_code(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 255,
    }
}

/// 2-bit-pack `seq` (4 bases/byte, little-endian within a byte). Non-ACGT bytes
/// are packed as `A` (code 0) and recorded as `(index, byte)` exceptions for
/// verbatim restore. Returns `(packed, exceptions)`.
fn pack(seq: &[u8]) -> (Vec<u8>, Vec<(usize, u8)>) {
    let mut packed = Vec::with_capacity(seq.len() / 4 + 1);
    let mut exceptions = Vec::new();
    let mut cur = 0u8;
    let mut nb = 0u8;
    for (i, &b) in seq.iter().enumerate() {
        let c = pack_code(b);
        let two = if c < 4 {
            c
        } else {
            exceptions.push((i, b));
            0
        };
        cur |= two << (2 * nb);
        nb += 1;
        if nb == 4 {
            packed.push(cur);
            cur = 0;
            nb = 0;
        }
    }
    if nb > 0 {
        packed.push(cur);
    }
    (packed, exceptions)
}

/// Inverse of [`pack`]: unpack `total` bases from `packed`, then restore
/// exceptions. `packed` must hold at least `total.div_ceil(4)` bytes.
fn unpack(packed: &[u8], total: usize, exceptions: &[(usize, u8)]) -> Result<Vec<u8>> {
    if packed.len() < total.div_ceil(4) {
        return Err(Error::Malformed("refpack: packed stream too short"));
    }
    // `total` is bounded by `packed.len() * 4` (checked above), but reserve
    // fallibly anyway so an over-large value can never abort on the allocation.
    let mut seq = Vec::new();
    seq.try_reserve_exact(total)
        .map_err(|_| Error::Malformed("refpack: output too large to allocate"))?;
    for i in 0..total {
        let byte = packed[i / 4];
        let c = ((byte >> (2 * (i % 4))) & 3) as usize;
        seq.push(SYM2BASE[c]);
    }
    for &(pos, b) in exceptions {
        *seq.get_mut(pos)
            .ok_or(Error::Malformed("refpack: exception out of range"))? = b;
    }
    Ok(seq)
}

/// Encode a whole reference (`lens` + concatenated `seq`) SPRING-style: 2-bit-pack
/// the bases and LZMA the packed stream; the per-contig lengths are LZMA'd too
/// (varints, ~0.5 MB raw on a 377 k-contig reference, so worth coding). Frame:
/// `[varint n_contigs][varint lens_raw_len][varint lens_coded_len][lens_coded]`
/// `[varint total_bases][varint n_exc][exc: dpos varint + byte]...[lzma(packed)]`.
/// (LZMA beat BWT decisively on the packed stream, so there is a single backend.)
pub(crate) fn encode(lens: &[u32], seq: &[u8]) -> Result<Vec<u8>> {
    let sum: usize = lens.iter().map(|&l| l as usize).sum();
    if sum != seq.len() {
        return Err(Error::Malformed("refpack: length/seq disagreement"));
    }
    let (packed, exceptions) = pack(seq);

    // Contig lengths as varints, then LZMA'd (concurrent with the bases).
    let mut lens_raw = Vec::with_capacity(lens.len() * 2);
    for &l in lens {
        write_varint(&mut lens_raw, u64::from(l));
    }
    let (packed_coded, lens_coded) = rayon::join(
        || reflzma::lzma_encode(&packed),
        || reflzma::lzma_encode(&lens_raw),
    );

    let mut out =
        Vec::with_capacity(packed_coded.len() + lens_coded.len() + exceptions.len() * 3 + 32);
    write_varint(&mut out, lens.len() as u64);
    write_varint(&mut out, lens_raw.len() as u64);
    write_varint(&mut out, lens_coded.len() as u64);
    out.extend_from_slice(&lens_coded);
    write_varint(&mut out, seq.len() as u64);
    write_varint(&mut out, exceptions.len() as u64);
    let mut prev = 0usize;
    for &(pos, b) in &exceptions {
        write_varint(&mut out, (pos - prev) as u64);
        out.push(b);
        prev = pos;
    }
    out.extend_from_slice(&packed_coded);
    Ok(out)
}

/// Inverse of [`encode`]: returns `(lens, seq)`.
pub(crate) fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut p = 0usize;
    let nc =
        read_varint(src, &mut p).ok_or(Error::Malformed("refpack: bad contig count"))? as usize;
    // Bomb guard scale for all declared lengths. (`nc` is NOT bounded by
    // `src.len()` — the lengths are compressed — but parsing `nc` varints from the
    // decompressed `lens_raw` self-limits, and the `lens` reserve is capped below.)
    let max_plausible = src.len().saturating_mul(1 << 18).saturating_add(1 << 18);
    if nc > max_plausible {
        return Err(Error::Malformed("refpack: implausible contig count"));
    }
    // Compressed contig lengths.
    let lens_raw_len =
        read_varint(src, &mut p).ok_or(Error::Malformed("refpack: bad lens raw len"))? as usize;
    let lens_coded_len =
        read_varint(src, &mut p).ok_or(Error::Malformed("refpack: bad lens coded len"))? as usize;
    if lens_raw_len > max_plausible {
        return Err(Error::Malformed("refpack: lens length exceeds capacity"));
    }
    // `lens_coded_len` is an untrusted varint, so `p + lens_coded_len` overflows.
    // Release builds disable overflow checks, so it wrapped to a value below `p`,
    // `get` saw an inverted range and returned None, and the error was right for
    // the wrong reason — but the fuzz build has debug assertions and panicked
    // here ("attempt to add with overflow", found by the `reorder` target once the
    // OOM above stopped masking it). Add checked, like every other offset in this
    // function already does.
    let end = p
        .checked_add(lens_coded_len)
        .ok_or(Error::Malformed("refpack: lens stream length overflows"))?;
    let lens_coded = src
        .get(p..end)
        .ok_or(Error::Malformed("refpack: truncated lens stream"))?;
    p = end;
    let lens_raw = reflzma::lzma_decode(lens_coded, lens_raw_len)?;
    let mut lens = Vec::with_capacity(nc.min(1 << 20));
    let mut lp = 0usize;
    for _ in 0..nc {
        lens.push(
            read_varint(&lens_raw, &mut lp).ok_or(Error::Malformed("refpack: bad contig len"))?
                as u32,
        );
    }
    let total = read_varint(src, &mut p).ok_or(Error::Malformed("refpack: bad total"))? as usize;
    if total > max_plausible {
        return Err(Error::Malformed(
            "refpack: declared length exceeds capacity",
        ));
    }
    let n_exc =
        read_varint(src, &mut p).ok_or(Error::Malformed("refpack: bad exc count"))? as usize;
    if n_exc > total {
        return Err(Error::Malformed("refpack: too many exceptions"));
    }
    let mut exceptions = Vec::with_capacity(n_exc.min(1 << 20));
    let mut pos = 0usize;
    for _ in 0..n_exc {
        let delta =
            read_varint(src, &mut p).ok_or(Error::Malformed("refpack: bad exc pos"))? as usize;
        pos = pos
            .checked_add(delta)
            .filter(|&x| x < total)
            .ok_or(Error::Malformed("refpack: exception position out of range"))?;
        let b = *src
            .get(p)
            .ok_or(Error::Malformed("refpack: truncated exc byte"))?;
        p += 1;
        exceptions.push((pos, b));
    }
    let packed = reflzma::lzma_decode(&src[p..], total.div_ceil(4))?;
    let seq = unpack(&packed, total, &exceptions)?;
    let s: usize = lens.iter().map(|&l| l as usize).sum();
    if s != seq.len() {
        return Err(Error::Malformed("refpack: length/seq disagreement"));
    }
    Ok((lens, seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(lens: &[u32], seq: &[u8]) {
        let coded = encode(lens, seq).expect("encode");
        let (out_lens, out_seq) = decode(&coded).expect("decode");
        assert_eq!(out_lens, lens, "lens");
        assert_eq!(out_seq, seq, "seq");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[], b"");
    }

    #[test]
    fn roundtrip_short_and_tail() {
        // Lengths that exercise every packed-tail remainder (0..3 bases over).
        roundtrip(&[1], b"A");
        roundtrip(&[2], b"AC");
        roundtrip(&[3], b"ACG");
        roundtrip(&[4], b"ACGT");
        roundtrip(&[5], b"ACGTA");
    }

    #[test]
    fn roundtrip_with_exceptions() {
        roundtrip(&[12], b"ACGTNNRYACGT");
    }

    #[test]
    fn roundtrip_repeats_compress() {
        let unit: &[u8] = b"ACGTACGGTATTCCGGAACCTTGGACGTACGGTATTCCGGAACCTTGG";
        let mut seq = Vec::new();
        let mut lens = Vec::new();
        for _ in 0..400 {
            seq.extend_from_slice(unit);
            lens.push(unit.len() as u32);
        }
        roundtrip(&lens, &seq);
        let coded = encode(&lens, &seq).unwrap();
        assert!(
            coded.len() < seq.len() / 16,
            "2-bit+LZMA should crush repeats: {} -> {}",
            seq.len(),
            coded.len()
        );
    }

    #[test]
    fn decode_rejects_garbage_without_panic() {
        for seed in 0u8..64 {
            let junk: Vec<u8> = (0..256u16).map(|i| (i as u8) ^ seed).collect();
            let _ = decode(&junk);
        }
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(b"ACGTNacgtRYKM".to_vec()), 0..90),
                0..70),
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            let coded = encode(&lens, &seq).expect("encode");
            let (out_lens, out_seq) = decode(&coded).expect("decode");
            proptest::prop_assert_eq!(out_lens, lens);
            proptest::prop_assert_eq!(out_seq, seq);
        }

        #[test]
        fn decode_never_panics(bytes in proptest::collection::vec(0u8..=255, 0..400)) {
            let _ = decode(&bytes);
        }
    }
}
