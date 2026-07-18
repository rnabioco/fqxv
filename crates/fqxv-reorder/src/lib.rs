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

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use fqxv_bytes::ReaderError;
pub(crate) use fqxv_bytes::{unzigzag, write_varint, zigzag};
pub(crate) use fqxv_dna::code_fold;
use thiserror::Error;

// Reverse complement lives in `fqxv-dna` now; re-exported so the long-standing
// `fqxv_reorder::revcomp` / `revcomp_into` paths (used by the container and
// examples) keep working.
pub use fqxv_dna::{revcomp, revcomp_into};

mod reflzma;
mod refpack;

/// A minimal integer hasher for the assembly maps. Their keys are already
/// well-mixed — 2-bit-packed k-mers and dense contig ids — so a single
/// multiplicative (Fibonacci) mix beats SipHash on the ~10^8 inserts/probes the
/// global assembler drives, and it is the throughput bottleneck of the v4 encode
/// path. Byte-output-preserving: these maps are only ever probed by key (never
/// iterated), and callers sort any candidate set deterministically, so the hash
/// choice cannot change the encoded stream.
#[derive(Default)]
pub(crate) struct IntHasher(u64);

impl Hasher for IntHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write_u64(&mut self, v: u64) {
        self.0 = v.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
    #[inline]
    fn write_u32(&mut self, v: u32) {
        self.0 = u64::from(v).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        // Only u32/u64 keys are hashed in this crate; keep a correct fallback so
        // the impl is total regardless of future key types.
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(b);
        }
        self.0 = self.0.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

/// A `HashMap` over integer keys using [`IntHasher`].
pub(crate) type IntMap<K, V> = HashMap<K, V, BuildHasherDefault<IntHasher>>;

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

mod clustered;
mod column;
mod global;
mod merge;
mod plan;
mod rescue;

pub use clustered::*;
pub(crate) use column::*;
pub use global::*;
pub use merge::*;
pub use plan::*;
pub use rescue::*;

/// Shared byte cursor specialized to this crate's [`Error`].
pub(crate) type Cursor<'a> = fqxv_bytes::Reader<'a, Error>;

/// The largest sequence, in bytes, this codec will reconstruct from one block.
///
/// Every rANS stream in a block is bounded by what the block could legitimately
/// hold, and this is that ceiling. It mirrors the container's per-block sequence
/// budget (`MAX_BLOCK_SEQ_BYTES`, 256 MiB) — a block never carries more, so no
/// honest stream decodes past it.
///
/// Duplicating the container's constant rather than taking it as a parameter is
/// deliberate: this crate sits *below* the container in the DAG, so it cannot
/// read it, and every public decode here is reachable from untrusted bytes (the
/// `reorder` fuzz target calls [`decode_clustered_auto`] on arbitrary input). If
/// the container's budget ever grows past this, decoding fails loudly here rather
/// than silently allowing a larger allocation.
pub const MAX_DECODED_BASES: usize = 256 << 20;

/// The largest read count this codec will accept in one block.
///
/// `n` arrives as a varint from an untrusted header and every per-read bound is
/// derived from it, so it must be sanity-checked before it multiplies. 2^26 reads
/// is far past any block the container writes and still keeps `n * 10` well under
/// a gigabyte.
pub const MAX_DECODED_READS: usize = 1 << 26;

/// Reject an untrusted read count before it is used to size anything.
pub(crate) fn check_n(n: usize) -> Result<usize> {
    if n > MAX_DECODED_READS {
        return Err(Error::Malformed("read count exceeds decode limit"));
    }
    Ok(n)
}

/// Bound for a stream holding one varint per read: 10 bytes is the longest LEB128
/// encoding of a `u64`, so no honest stream of `n` varints exceeds this.
pub(crate) fn per_read_varints(n: usize) -> usize {
    n.saturating_mul(10)
}

impl ReaderError for Error {
    fn truncated() -> Self {
        Error::Malformed("truncated")
    }
    fn bad_varint() -> Self {
        Error::Malformed("varint too long")
    }
    fn oversized() -> Self {
        Error::Malformed("length count too large to allocate")
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
    fn alloc_read_rejects_huge_len() {
        // The per-read length guard used by the clustered decoders: a corrupt
        // length must error rather than abort on a huge infallible allocation.
        assert!(alloc_read(usize::MAX).is_err());
        assert_eq!(alloc_read(8).unwrap().len(), 8);
    }

    #[test]
    fn decode_blocked_rejects_huge_block_count() {
        // A corrupt reference block count must fail the reservation, not abort.
        let mut src = Vec::new();
        write_varint(&mut src, 1u64 << 62); // absurd block count
        assert!(GlobalReference::decode_blocked(&src).is_err());
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
        let enc = encode_clustered(&reads, &vec![0u32; reads.len()], 4).expect("encode");
        let dec = decode_clustered(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn clustered_empty() {
        let enc = encode_clustered(&[], &[], 4).unwrap();
        assert!(decode_clustered(&enc).unwrap().is_empty());
    }

    /// A crafted archive must not be able to dictate an allocation.
    ///
    /// Each rANS stream states its own decoded length as a `u64` in its header,
    /// and unbounded `decode` simply believes it: claim 2 GiB and 2 GiB is
    /// allocated before a byte is validated. The container CRCs every payload
    /// first, which is why this was survivable — but a CRC detects accidents, not
    /// intent, and recomputing a valid one over a hostile stream is free. So the
    /// bound, not the CRC, is what closes this.
    ///
    /// Rewrites the declared length of the FIRST stream (`ops`), which is exactly
    /// one byte per read and so is bounded by `n`.
    #[test]
    fn a_crafted_output_length_cannot_force_an_allocation() {
        let reads: Vec<&[u8]> = vec![b"ACGTACGTACGT", b"ACGTACGTACGT", b"ACGTAAGTACGT"];
        let mut enc = encode_clustered(&reads, &vec![0u32; reads.len()], 4).expect("encode");
        assert!(decode_clustered(&enc).is_ok(), "fixture must decode clean");

        // [version u8][varint n][varint comp_len][rANS stream ...], and the rANS
        // header carries the decoded length at its bytes 1..9.
        let mut pos = 1usize;
        fqxv_bytes::read_varint(&enc, &mut pos).expect("n");
        let comp_len = fqxv_bytes::read_varint(&enc, &mut pos).expect("comp_len") as usize;
        assert!(
            comp_len >= 9 && pos + 9 <= enc.len(),
            "stream header present"
        );
        enc[pos + 1..pos + 9].copy_from_slice(&(1u64 << 31).to_le_bytes()); // claim 2 GiB

        // Assert WHICH error, not merely that there is one. `is_err()` alone does
        // not test this: unbounded `decode` also returns Err here — it allocates
        // and zeroes the 2 GiB first, then fails on input exhaustion. Measured,
        // this test passes either way and takes 10.9s unbounded against 0.38s
        // bounded, the runtime being the only tell. "Returned an error" and
        // "refused to allocate" are different claims, and only the second one is
        // the fix.
        let err = decode_clustered(&enc).expect_err("2 GiB claim must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("exceeds caller bound"),
            "must be refused BY THE BOUND, before allocating; got: {msg}"
        );
    }

    /// The read count is the seed for every per-read bound, and it arrives as an
    /// untrusted varint — so it has to be checked before it multiplies.
    #[test]
    fn an_absurd_read_count_is_refused() {
        let reads: Vec<&[u8]> = vec![b"ACGTACGTACGT"];
        let enc = encode_clustered(&reads, &[0u32], 4).expect("encode");
        let mut bomb = vec![enc[0]];
        fqxv_bytes::write_varint(&mut bomb, u64::MAX / 2);
        bomb.extend_from_slice(&enc[2..]);
        let err = decode_clustered(&bomb).expect_err("absurd read count must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("read count exceeds decode limit"),
            "must be refused by the count check, before it seeds any bound; got: {msg}"
        );
    }

    #[test]
    fn clustered_deduplicates() {
        // Many copies of one read collapse to MATCH ops; the encoded size is
        // dominated by fixed per-stream overhead, not the read count.
        let read = b"ACGTTTGACCGATTGCAACGTTTGACCGATTGCA";
        let reads: Vec<&[u8]> = vec![&read[..]; 10_000];
        let raw = read.len() * reads.len();
        let enc = encode_clustered(&reads, &vec![0u32; reads.len()], 6).unwrap();
        assert!(
            enc.len() < raw / 50,
            "expected heavy dedup, got {} for {raw} raw",
            enc.len()
        );
        assert_eq!(decode_clustered(&enc).unwrap().len(), 10_000);
    }

    #[test]
    fn clustered_shift_overlap_roundtrips() {
        // Overlapping windows of a reference share sequence at a shift, which
        // should trigger SHIFT ops (and round-trip exactly).
        let reference = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGCATCGATCGTAGCTAGCAT";
        let win = 30usize;
        let (mut lens, mut seq) = (Vec::new(), Vec::new());
        for start in 0..=(reference.len() - win) {
            seq.extend_from_slice(&reference[start..start + win]);
            lens.push(win as u32);
        }
        let p = plan(&lens, &seq, DEFAULT_K);
        let mut offs = vec![0usize];
        for &l in &lens {
            offs.push(offs.last().unwrap() + l as usize);
        }
        let cl: Vec<Vec<u8>> = p
            .order
            .iter()
            .map(|&oi| {
                let oi = oi as usize;
                let s = &seq[offs[oi]..offs[oi + 1]];
                if p.flip[oi] {
                    revcomp(s)
                } else {
                    s.to_vec()
                }
            })
            .collect();
        let refs: Vec<&[u8]> = cl.iter().map(Vec::as_slice).collect();
        let anchors: Vec<u32> = p.order.iter().map(|&oi| p.anchor[oi as usize]).collect();
        let enc = encode_clustered(&refs, &anchors, 8).unwrap();
        let expect: Vec<Vec<u8>> = refs.iter().map(|r| r.to_vec()).collect();
        assert_eq!(decode_clustered(&enc).unwrap(), expect);
    }

    #[test]
    fn rescue_roundtrip() {
        // The same mix the v2 codec round-trips must also round-trip under rescue.
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGTACGT",
            b"ACGTACGTACGT", // match
            b"ACGTAAGTACGT", // 1 mismatch
            b"ACGTACGTNCGT", // mismatch incl. N
            b"TTTTGGGGCCCC", // literal
            b"",             // empty
            b"",             // match (empty == empty)
        ];
        let enc = encode_clustered_rescue(&reads, &vec![0u32; reads.len()], 4).expect("encode");
        let dec = decode_clustered_rescue(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn rescue_attaches_to_earlier_contig() {
        // Two unrelated references interleaved: A, B, A, B, A. After B seeds the
        // "current" contig, the v2 codec would strand the next A as a LITERAL;
        // the rescue index lets it re-attach to the earlier A contig. This asserts
        // the round-trip; op_stats_rescue reports the recovered literal.
        let a = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGC";
        let b = b"TTAGGCCATTACAGGTACCATGACATTGGACATTACAGGTTCAAGT";
        let reads: Vec<&[u8]> = vec![&a[..], &b[..], &a[..], &b[..], &a[..]];
        let anchors = vec![0u32; reads.len()];
        let enc = encode_clustered_rescue(&reads, &anchors, 6).expect("encode");
        let dec = decode_clustered_rescue(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
        // The third A (index 2) should be rescued onto A's contig, not stranded:
        // only the first A and the first B are literals.
        let st = op_stats_rescue(&reads, &anchors);
        assert_eq!(st.literals, 2, "expected only the two seeds to be literals");
    }

    /// Encode a whole read set with the two-pass v4 codec as a single block and
    /// decode it against the frozen reference; must be byte-exact.
    fn v4_roundtrip_one_block(reads: &[&[u8]], anchors: &[u32]) -> Vec<Vec<u8>> {
        let (reference, places) = assemble_global(reads, anchors);
        // Reference must serialize/deserialize losslessly too.
        let ref_bytes = reference.encode(6, 0, 0).expect("ref encode");
        let reference = GlobalReference::decode(&ref_bytes).expect("ref decode");
        let enc = encode_global_block(reads, &places, &reference).expect("v4 encode");
        // Determinism: re-encoding the same input is byte-identical.
        assert_eq!(
            encode_global_block(reads, &places, &reference).expect("v4 again"),
            enc,
            "v4 encode not deterministic"
        );
        decode_global_block(&enc, &reference).expect("v4 decode")
    }

    /// Assemble, overlap-merge the reference, then encode/decode every read
    /// against the MERGED reference — the placement remap must stay byte-exact.
    fn v4_merged_roundtrip(reads: &[&[u8]], anchors: &[u32]) -> (usize, usize, Vec<Vec<u8>>) {
        let (reference, places) = assemble_global(reads, anchors);
        let before = reference.n_contigs();
        let (merged, mplaces) = merge_reference(reads, &reference, &places);
        let after = merged.n_contigs();
        // Merged reference must serialize/reload and every read must round-trip.
        let ref_bytes = merged.encode(6, 0, 0).expect("ref encode");
        let merged = GlobalReference::decode(&ref_bytes).expect("ref decode");
        let enc = encode_global_block(reads, &mplaces, &merged).expect("v4 encode");
        let dec = decode_global_block(&enc, &merged).expect("v4 decode");
        (before, after, dec)
    }

    #[test]
    fn merge_reference_roundtrips_and_shrinks() {
        // Overlapping windows of one reference fragment into several contigs under
        // the greedy pass; overlap-merge should chain them into fewer contigs and
        // still round-trip byte-exactly.
        let reference_seq = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGCATCGATCGTAGCTAGCATTTACAGGTACCATGACATTGG";
        let win = 40usize;
        let (mut lens, mut seq) = (Vec::new(), Vec::new());
        for start in (0..=(reference_seq.len() - win)).step_by(3) {
            seq.extend_from_slice(&reference_seq[start..start + win]);
            lens.push(win as u32);
        }
        let p = plan(&lens, &seq, DEFAULT_K);
        let mut offs = vec![0usize];
        for &l in &lens {
            offs.push(offs.last().unwrap() + l as usize);
        }
        let cl: Vec<Vec<u8>> = p
            .order
            .iter()
            .map(|&oi| {
                let oi = oi as usize;
                let s = &seq[offs[oi]..offs[oi + 1]];
                if p.flip[oi] {
                    revcomp(s)
                } else {
                    s.to_vec()
                }
            })
            .collect();
        let refs: Vec<&[u8]> = cl.iter().map(Vec::as_slice).collect();
        let anchors: Vec<u32> = p.order.iter().map(|&oi| p.anchor[oi as usize]).collect();
        let (before, after, dec) = v4_merged_roundtrip(&refs, &anchors);
        let expect: Vec<Vec<u8>> = refs.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect, "merged reference must round-trip");
        assert!(after <= before, "merge must not increase contig count");
    }

    #[test]
    fn global_roundtrip() {
        let reads: Vec<&[u8]> = vec![
            b"ACGTACGTACGT",
            b"ACGTACGTACGT", // match
            b"ACGTAAGTACGT", // 1 mismatch
            b"ACGTACGTNCGT", // mismatch incl. N
            b"TTTTGGGGCCCC", // seeds a second contig
            b"",             // empty
            b"",             // match (empty == empty)
        ];
        let dec = v4_roundtrip_one_block(&reads, &vec![0u32; reads.len()]);
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn global_attaches_to_earlier_contig() {
        // The A,B,A,B,A interleave that strands the third A as a LITERAL under
        // v2: v4's global reference has one A contig every A read maps onto, so
        // there are exactly two contigs (one A, one B) — the reference dedups.
        let a = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGC";
        let b = b"TTAGGCCATTACAGGTACCATGACATTGGACATTACAGGTTCAAGT";
        let reads: Vec<&[u8]> = vec![&a[..], &b[..], &a[..], &b[..], &a[..]];
        let anchors = vec![0u32; reads.len()];
        let (reference, _places) = assemble_global(&reads, &anchors);
        assert_eq!(
            reference.n_contigs(),
            2,
            "reference should hold two contigs"
        );
        let dec = v4_roundtrip_one_block(&reads, &anchors);
        let expect: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
        assert_eq!(dec, expect);
    }

    #[test]
    fn global_multi_block_shares_reference() {
        // Assemble globally, then code in several small blocks against the one
        // frozen reference (the container's parallel-block shape). Reads at block
        // boundaries can't be block-local MATCH, so this exercises the `(ci, off)`
        // fallback placement for every read.
        let reference_seq = b"ACGTTGCAACCGGTTACGTAGCTAGCATCGATCGATCGTAGCATGCATCGATCGTAGCTAGCAT";
        let win = 30usize;
        let (mut lens, mut seq) = (Vec::new(), Vec::new());
        for start in 0..=(reference_seq.len() - win) {
            seq.extend_from_slice(&reference_seq[start..start + win]);
            lens.push(win as u32);
        }
        let p = plan(&lens, &seq, DEFAULT_K);
        let mut offs = vec![0usize];
        for &l in &lens {
            offs.push(offs.last().unwrap() + l as usize);
        }
        let cl: Vec<Vec<u8>> = p
            .order
            .iter()
            .map(|&oi| {
                let oi = oi as usize;
                let s = &seq[offs[oi]..offs[oi + 1]];
                if p.flip[oi] {
                    revcomp(s)
                } else {
                    s.to_vec()
                }
            })
            .collect();
        let refs: Vec<&[u8]> = cl.iter().map(Vec::as_slice).collect();
        let anchors: Vec<u32> = p.order.iter().map(|&oi| p.anchor[oi as usize]).collect();

        let (reference, places) = assemble_global(&refs, &anchors);
        let ref_bytes = reference.encode(8, 0, 0).expect("ref encode");
        let reference = GlobalReference::decode(&ref_bytes).expect("ref decode");
        // Cut into 4-read blocks and code/decode each against the shared ref.
        let mut got: Vec<Vec<u8>> = Vec::new();
        let mut s = 0usize;
        while s < refs.len() {
            let e = (s + 4).min(refs.len());
            let enc = encode_global_block(&refs[s..e], &places[s..e], &reference).expect("enc");
            got.extend(decode_global_block(&enc, &reference).expect("dec"));
            s = e;
        }
        let expect: Vec<Vec<u8>> = refs.iter().map(|r| r.to_vec()).collect();
        assert_eq!(got, expect);
    }

    proptest::proptest! {
        #[test]
        fn rescue_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..30),
                0..50)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let enc = encode_clustered_rescue(&refs, &vec![0u32; refs.len()], 4).expect("encode");
            let dec = decode_clustered_rescue(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, reads);
        }

        // Harden the rescue codec against the way the container actually drives it:
        // longer reads, a realistic `seq_order`, and NON-uniform anchors (the
        // anchors above are all zero, so the anchor-implied offset path — the v2
        // fast candidate inside `Assembler::place` — went largely unexercised).
        // Also pins encode determinism, which the container's determinism invariant
        // relies on now that rescue is on by default.
        #[test]
        fn rescue_roundtrip_varied_anchors(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..80)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let anchors: Vec<u32> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| ((r.len() * 7 + i * 3) % 41) as u32)
                .collect();
            for order in [1usize, 8] {
                let v3 = encode_clustered_rescue(&refs, &anchors, order).expect("v3 encode");
                proptest::prop_assert_eq!(
                    decode_clustered_rescue(&v3).expect("v3 decode"),
                    reads.clone()
                );
                // Deterministic: re-encoding the same input yields identical bytes.
                proptest::prop_assert_eq!(
                    encode_clustered_rescue(&refs, &anchors, order).expect("v3 again"),
                    v3
                );
                // v2 must also round-trip under the same non-uniform anchors.
                let v2 = encode_clustered(&refs, &anchors, order).expect("v2 encode");
                proptest::prop_assert_eq!(
                    decode_clustered(&v2).expect("v2 decode"),
                    reads.clone()
                );
            }
        }

        #[test]
        fn clustered_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..30),
                0..50)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let enc = encode_clustered(&refs, &vec![0u32; refs.len()], 4).expect("encode");
            let dec = decode_clustered(&enc).expect("decode");
            proptest::prop_assert_eq!(dec, reads);
        }

        // v4 two-pass: assemble globally, serialize+reload the reference, then
        // code the whole set as one block and as several small blocks (the
        // parallel-decode shape) — both against the shared frozen reference. Uses
        // non-uniform anchors so the anchor-implied candidate path is exercised.
        #[test]
        fn global_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..80)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let anchors: Vec<u32> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| ((r.len() * 7 + i * 3) % 41) as u32)
                .collect();
            let (reference, places) = assemble_global(&refs, &anchors);
            let ref_bytes = reference.encode(4, 0, 0).expect("ref encode");
            proptest::prop_assert_eq!(
                reference.encode(4, 0, 0).expect("ref again"), ref_bytes.clone(),
                "reference encode not deterministic"
            );
            let reference = GlobalReference::decode(&ref_bytes).expect("ref decode");
            // Single block.
            let enc = encode_global_block(&refs, &places, &reference).expect("v4 encode");
            proptest::prop_assert_eq!(
                decode_global_block(&enc, &reference).expect("v4 decode"),
                reads.clone()
            );
            // Several small blocks sharing the reference.
            let mut got: Vec<Vec<u8>> = Vec::new();
            let mut s = 0usize;
            while s < refs.len() {
                let e = (s + 7).min(refs.len());
                let b = encode_global_block(&refs[s..e], &places[s..e], &reference).expect("enc");
                got.extend(decode_global_block(&b, &reference).expect("dec"));
                s = e;
            }
            proptest::prop_assert_eq!(got, reads);
        }

        // The overlap-merge refinement must preserve the v4 round-trip: assemble,
        // merge the reference (remapping placements), then encode/decode every read
        // against the merged reference. Also pins merge determinism and that it
        // never grows the contig set.
        #[test]
        fn merge_reference_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..80)
        ) {
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            let anchors: Vec<u32> = reads
                .iter()
                .enumerate()
                .map(|(i, r)| ((r.len() * 5 + i * 7) % 37) as u32)
                .collect();
            let (reference, places) = assemble_global(&refs, &anchors);
            let (merged, mplaces) = merge_reference(&refs, &reference, &places);
            proptest::prop_assert!(merged.n_contigs() <= reference.n_contigs());
            // Merge is deterministic.
            let (merged2, _) = merge_reference(&refs, &reference, &places);
            proptest::prop_assert_eq!(merged.total_bases(), merged2.total_bases());
            let ref_bytes = merged.encode(4, 0, 0).expect("ref encode");
            let merged = GlobalReference::decode(&ref_bytes).expect("ref decode");
            let enc = encode_global_block(&refs, &mplaces, &merged).expect("v4 encode");
            proptest::prop_assert_eq!(
                decode_global_block(&enc, &merged).expect("v4 decode"),
                reads.clone()
            );
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
