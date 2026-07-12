//! Reference-free read reordering by canonical-minimizer clustering.
//!
//! Each read is reduced to its minimum canonical k-mer (the smaller of a k-mer
//! and its reverse complement, over the whole read). Sorting reads by that key
//! (then by oriented sequence) brings exact and reverse-complement duplicates —
//! and near-duplicates sharing a minimizer — next to each other, so the order-k
//! sequence model in `fqxv-seq` sees long runs of near-identical reads. This is
//! the cross-read redundancy lever (PgRC/SPRING-class), the piece that plain
//! per-read context modeling can't reach.
//!
//! [`plan`] returns the emission [`Plan::order`] and per-read [`Plan::flip`]
//! (whether a read is stored reverse-complemented so its minimizer is in
//! canonical orientation). The caller reorders names/sequence/quality by `order`,
//! reverse-complements + reverses the quality of flipped reads, and stores the
//! permutation to restore the original order.
//!
//! ```
//! use fqxv_reorder::{plan, revcomp};
//! // read 1 is the reverse complement of read 0 — they should cluster.
//! let a = b"ACGTTTGACCGATT";
//! let ra = revcomp(a);
//! let mut seq = a.to_vec();
//! seq.extend_from_slice(&ra);
//! let lens = [a.len() as u32, ra.len() as u32];
//! let p = plan(&lens, &seq, 7);
//! assert_eq!(p.order.len(), 2);
//! // exactly one of the mates is flipped so both store identically.
//! assert_ne!(p.flip[0], p.flip[1]);
//! ```

use rayon::prelude::*;
use thiserror::Error;

/// Errors returned by the reordering engine.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed clustered stream: {0}")]
    Malformed(&'static str),
    /// Underlying rANS coder failure.
    #[error(transparent)]
    Rans(#[from] fqxv_rans::Error),
    /// Underlying sequence codec failure.
    #[error(transparent)]
    Seq(#[from] fqxv_seq::Error),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// A reordering: the order to emit reads, and which reads are stored
/// reverse-complemented.
#[derive(Debug, Clone, Default)]
pub struct Plan {
    /// Original read indices in the order they should be emitted.
    pub order: Vec<u32>,
    /// Per original-read flag: store the read reverse-complemented.
    pub flip: Vec<bool>,
}

/// Default minimizer k-mer length.
pub const DEFAULT_K: usize = 15;

/// byte -> 2-bit code, 255 for non-ACGT.
#[inline]
fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 255,
    }
}

/// Reverse complement (complementing A/C/G/T in either case; other bytes,
/// including `N`, are passed through), reversed.
#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// Minimum canonical k-mer of `read` and whether that minimizer sits on the
/// reverse strand (i.e. the read should be flipped to canonicalize it).
fn min_canonical(read: &[u8], k: usize) -> (u64, bool) {
    if read.len() < k || k == 0 || k > 32 {
        return (u64::MAX, false);
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k as u64 - 1);
    let (mut fwd, mut rc, mut valid) = (0u64, 0u64, 0usize);
    let (mut best, mut best_flip) = (u64::MAX, false);
    for &b in read {
        let c = code(b);
        if c == 255 {
            fwd = 0;
            rc = 0;
            valid = 0;
            continue;
        }
        let c = u64::from(c);
        fwd = ((fwd << 2) | c) & mask;
        rc = ((rc >> 2) | ((3 - c) << shift)) & mask;
        valid += 1;
        if valid >= k {
            let (canon, is_rc) = if fwd <= rc { (fwd, false) } else { (rc, true) };
            if canon < best {
                best = canon;
                best_flip = is_rc;
            }
        }
    }
    (best, best_flip)
}

/// Build a clustering [`Plan`] for the reads in `seq` (lengths in `lens`).
///
/// `k` is the minimizer length (see [`DEFAULT_K`]); it is clamped to `1..=32`.
#[must_use]
pub fn plan(lens: &[u32], seq: &[u8], k: usize) -> Plan {
    let k = k.clamp(1, 32);
    let n = lens.len();

    // Byte offset of each read (so the key build can run in parallel).
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // (canonical minimizer, oriented sequence, original index, flip). Building
    // each key is independent, so it runs across cores.
    let mut keys: Vec<(u64, Vec<u8>, u32, bool)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let read = &seq[offs[i]..offs[i + 1]];
            let (canon, flip) = min_canonical(read, k);
            let oriented = if flip { revcomp(read) } else { read.to_vec() };
            (canon, oriented, i as u32, flip)
        })
        .collect();

    // Parallel sort. The final tie-break on original index makes the comparator
    // a TOTAL order, so duplicate reads (equal minimizer + sequence) get one
    // deterministic ordering regardless of thread count — the byte-identical
    // output invariant. Clustering quality is unaffected (dups stay adjacent).
    keys.par_sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
    });

    let order = keys.iter().map(|key| key.2).collect();
    let mut flip = vec![false; n];
    for key in &keys {
        flip[key.2 as usize] = key.3;
    }
    Plan { order, flip }
}

// --- clustered differential codec -------------------------------------------

const OP_MATCH: u8 = 0;
const OP_DELTA: u8 = 1;
const OP_LITERAL: u8 = 2;

/// Differentially encode reads that are already in clustered order (e.g. via
/// [`plan`]).
///
/// Each read is coded relative to the previous read: `MATCH` (identical),
/// `DELTA` (same length, at most ¼ mismatches — stored as position + byte), or
/// `LITERAL` (coded with the [`fqxv_seq`] context model). Byte-exact, so it
/// preserves `N` and any other bytes. `seq_order` is the literal model order.
pub fn encode_clustered(reads: &[&[u8]], seq_order: usize) -> Result<Vec<u8>> {
    let mut ops = Vec::with_capacity(reads.len());
    let (mut nmis, mut pos, mut subs) = (Vec::new(), Vec::new(), Vec::new());
    let (mut lit_seq, mut lit_lens): (Vec<u8>, Vec<u32>) = (Vec::new(), Vec::new());

    for (i, &cur) in reads.iter().enumerate() {
        // The first read has no predecessor, so it is always a literal.
        let prev: &[u8] = if i > 0 { reads[i - 1] } else { &[] };
        let mism = (i > 0 && cur.len() == prev.len()).then(|| {
            (0..cur.len())
                .filter(|&j| cur[j] != prev[j])
                .collect::<Vec<_>>()
        });
        if i > 0 && cur == prev {
            ops.push(OP_MATCH);
        } else if let Some(mism) = mism.filter(|m| m.len() <= cur.len() / 4) {
            ops.push(OP_DELTA);
            write_varint(&mut nmis, mism.len() as u64);
            let mut last = 0usize;
            for &m in &mism {
                write_varint(&mut pos, (m - last) as u64);
                last = m;
                subs.push(cur[m]);
            }
        } else {
            ops.push(OP_LITERAL);
            lit_seq.extend_from_slice(cur);
            lit_lens.push(cur.len() as u32);
        }
    }

    let ops_c = fqxv_rans::encode(&ops, fqxv_rans::Order::One)?;
    let nmis_c = fqxv_rans::encode(&nmis, fqxv_rans::Order::Zero)?;
    let pos_c = fqxv_rans::encode(&pos, fqxv_rans::Order::Zero)?;
    let subs_c = fqxv_rans::encode(&subs, fqxv_rans::Order::One)?;
    let lit_c = fqxv_seq::encode(&lit_lens, &lit_seq, seq_order)?;

    let mut out = Vec::new();
    out.push(0u8); // version
    write_varint(&mut out, reads.len() as u64);
    for s in [&ops_c, &nmis_c, &pos_c, &subs_c, &lit_c] {
        write_varint(&mut out, s.len() as u64);
        out.extend_from_slice(s);
    }
    Ok(out)
}

/// Decode a stream produced by [`encode_clustered`], returning the reads in the
/// same (clustered) order.
pub fn decode_clustered(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != 0 {
        return Err(Error::Malformed("unsupported version"));
    }
    let n = r.varint()? as usize;
    let ops = fqxv_rans::decode(r.take_stream()?)?;
    let nmis = fqxv_rans::decode(r.take_stream()?)?;
    let pos = fqxv_rans::decode(r.take_stream()?)?;
    let subs = fqxv_rans::decode(r.take_stream()?)?;
    let (lit_lens, lit_seq) = fqxv_seq::decode(r.take_stream()?)?;

    let mut c_nmis = Cursor::new(&nmis);
    let mut c_pos = Cursor::new(&pos);
    let (mut subs_pos, mut lit_pos, mut lit_idx) = (0usize, 0usize, 0usize);
    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n.min(1 << 22));

    for i in 0..n {
        let op = *ops.get(i).ok_or(Error::Malformed("op underrun"))?;
        let read = match op {
            OP_MATCH => reads
                .last()
                .ok_or(Error::Malformed("MATCH with no previous"))?
                .clone(),
            OP_DELTA => {
                let mut read = reads
                    .last()
                    .ok_or(Error::Malformed("DELTA with no previous"))?
                    .clone();
                let m = c_nmis.varint()? as usize;
                let mut p = 0usize;
                for _ in 0..m {
                    p += c_pos.varint()? as usize;
                    let b = *subs
                        .get(subs_pos)
                        .ok_or(Error::Malformed("subs underrun"))?;
                    subs_pos += 1;
                    *read
                        .get_mut(p)
                        .ok_or(Error::Malformed("delta position out of range"))? = b;
                }
                read
            }
            OP_LITERAL => {
                let l = *lit_lens
                    .get(lit_idx)
                    .ok_or(Error::Malformed("lit len underrun"))? as usize;
                lit_idx += 1;
                let bytes = lit_seq
                    .get(lit_pos..lit_pos + l)
                    .ok_or(Error::Malformed("lit data underrun"))?
                    .to_vec();
                lit_pos += l;
                bytes
            }
            _ => return Err(Error::Malformed("unknown op")),
        };
        reads.push(read);
    }
    Ok(reads)
}

fn write_varint(out: &mut Vec<u8>, mut v: u64) {
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

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Cursor { buf, pos: 0 }
    }
    fn u8(&mut self) -> Result<u8> {
        let b = *self
            .buf
            .get(self.pos)
            .ok_or(Error::Malformed("truncated"))?;
        self.pos += 1;
        Ok(b)
    }
    fn varint(&mut self) -> Result<u64> {
        let (mut v, mut shift) = (0u64, 0u32);
        loop {
            let b = self.u8()?;
            v |= u64::from(b & 0x7f) << shift;
            if b & 0x80 == 0 {
                return Ok(v);
            }
            shift += 7;
            if shift >= 64 {
                return Err(Error::Malformed("varint too long"));
            }
        }
    }
    fn take_stream(&mut self) -> Result<&'a [u8]> {
        let n = self.varint()? as usize;
        let end = self.pos + n;
        let s = self
            .buf
            .get(self.pos..end)
            .ok_or(Error::Malformed("stream truncated"))?;
        self.pos = end;
        Ok(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_basic() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
        assert_eq!(revcomp(b"AACG"), b"CGTT");
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
    }

    #[test]
    fn order_is_a_permutation() {
        let reads: Vec<&[u8]> = vec![b"ACGTACGT", b"TTTTGGGG", b"ACGTACGT", b"CCCCAAAA"];
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.concat();
        let p = plan(&lens, &seq, 5);
        let mut sorted = p.order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3]);
        assert_eq!(p.flip.len(), 4);
    }

    #[test]
    fn duplicates_cluster_adjacently() {
        // Two copies of read A and two of read B, shuffled; A's and B's should
        // each become adjacent after planning.
        let a: &[u8] = b"ACGTTTGACCGATTGCA";
        let b: &[u8] = b"GGGGCCCCAAAATTTTG";
        let reads = [a, b, a, b, a];
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.concat();
        let p = plan(&lens, &seq, 9);
        // Map emitted order back to which template (a=0/b=1) each read is.
        let tmpl: Vec<u8> = p
            .order
            .iter()
            .map(|&i| u8::from(reads[i as usize] == b))
            .collect();
        // Count adjacency changes; clustered => few transitions (<= 1).
        let transitions = tmpl.windows(2).filter(|w| w[0] != w[1]).count();
        assert!(transitions <= 1, "reads not clustered: {tmpl:?}");
    }

    #[test]
    fn revcomp_duplicate_flips_to_match() {
        let a = b"ACGTTTGACCGATTGCA";
        let ra = revcomp(a);
        let seq: Vec<u8> = a.iter().chain(ra.iter()).copied().collect();
        let lens = [a.len() as u32, ra.len() as u32];
        let p = plan(&lens, &seq, 9);
        // The two mates share a canonical minimizer; exactly one is flipped so
        // both are stored in the same orientation.
        assert_ne!(p.flip[0], p.flip[1]);
        // After orienting, the stored sequences are identical.
        let s0 = if p.flip[0] { revcomp(a) } else { a.to_vec() };
        let s1 = if p.flip[1] { revcomp(&ra) } else { ra.clone() };
        assert_eq!(s0, s1);
    }

    #[test]
    fn handles_short_and_n_reads() {
        let reads: Vec<&[u8]> = vec![b"AC", b"NNNNNNN", b"", b"ACGTACGTAC"];
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.concat();
        let p = plan(&lens, &seq, 5);
        assert_eq!(p.order.len(), 4);
    }

    #[test]
    fn clustered_roundtrip() {
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGTACGT", // literal (first)
            b"ACGTACGTACGT", // match
            b"ACGTAAGTACGT", // delta (1 mismatch)
            b"ACGTACGTNCGT", // delta including an N
            b"TTTTGGGGCCCC", // literal
            b"",             // empty read
            b"",             // match (empty == empty)
        ];
        let enc = encode_clustered(&reads, 4).expect("encode");
        let dec = decode_clustered(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn clustered_empty() {
        let enc = encode_clustered(&[], 4).unwrap();
        assert!(decode_clustered(&enc).unwrap().is_empty());
    }

    #[test]
    fn clustered_deduplicates() {
        // Many copies of one read collapse to MATCH ops; the encoded size is
        // dominated by fixed per-stream overhead, not the read count.
        let read = b"ACGTTTGACCGATTGCAACGTTTGACCGATTGCA";
        let reads: Vec<&[u8]> = vec![&read[..]; 10_000];
        let raw = read.len() * reads.len();
        let enc = encode_clustered(&reads, 6).unwrap();
        assert!(
            enc.len() < raw / 50,
            "expected heavy dedup, got {} for {raw} raw",
            enc.len()
        );
        assert_eq!(decode_clustered(&enc).unwrap().len(), 10_000);
    }

    proptest::proptest! {
        #[test]
        fn clustered_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..30),
                0..50)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let enc = encode_clustered(&refs, 4).expect("encode");
            let dec = decode_clustered(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, reads);
        }

        #[test]
        fn plan_is_valid_permutation(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..40),
                0..50),
            k in 1usize..=20,
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            let p = plan(&lens, &seq, k);
            let mut sorted = p.order.clone();
            sorted.sort_unstable();
            let expect: Vec<u32> = (0..reads.len() as u32).collect();
            proptest::prop_assert_eq!(sorted, expect);
            proptest::prop_assert_eq!(p.flip.len(), reads.len());
        }
    }
}
