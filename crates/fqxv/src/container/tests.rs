use super::verify::verify_whole_file_crc;
use super::*;

const SAMPLE: &[u8] = b"\
@SRR1.1 INST:1:FC:1:1101:1000:2000 length=8\n\
ACGTACGT\n\
+SRR1.1 INST:1:FC:1:1101:1000:2000 length=8\n\
IIIIFFF#\n\
@SRR1.2 INST:1:FC:1:1101:1005:2050 length=8\n\
NNGGCCTA\n\
+\n\
###IIIFF\n";

fn compress_bytes(input: &[u8], params: Params) -> Vec<u8> {
    let mut out = Vec::new();
    compress(input, &mut out, params).expect("compress");
    out
}

/// A truncated FASTQ must be an error in BOTH `--order` modes, not an abort in
/// one of them.
///
/// The reorder path derives its byte offsets from `lens` alone and then indexes
/// `seq`/`qual` with them, so a quality line cut short made every slice out of
/// bounds — and the release profile is `panic = "abort"`, so `--order any`
/// core-dumped on a file `--order preserve` rejected cleanly. Truncation is the
/// most ordinary corruption there is (an interrupted download, a broken pipe),
/// and the two modes disagreeing about whether a file is valid is the bug.
///
/// Asserts they agree, rather than merely that each returns something: an
/// `is_err()` on the reorder path alone would have passed while it aborted, since
/// an abort is not an `Err`.
#[test]
fn a_truncated_fastq_errors_in_both_order_modes() {
    // Last record's quality line is cut short: lens sum past the qual buffer.
    let truncated = b"@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2\nTTGGCCAATTGG\n+\nFFF";

    let plain = compress(&truncated[..], &mut Vec::new(), Params::default());
    let reordered = compress(
        &truncated[..],
        &mut Vec::new(),
        Params {
            reorder: true,
            ..Params::default()
        },
    );

    assert!(plain.is_err(), "plain path must reject a truncated record");
    assert!(
        reordered.is_err(),
        "reorder path must reject it too — it used to abort here"
    );
    assert_eq!(
        plain.unwrap_err().to_string(),
        reordered.unwrap_err().to_string(),
        "both modes must reject the same input with the same message"
    );
}

/// A per-record sequence/quality length mismatch must be rejected at parse time,
/// even when a second record mis-compensates so the *block totals* match.
///
/// The block-level check only compares total quality bytes against the summed
/// read lengths, so `(seq 4, qual 2)` followed by `(seq 2, qual 4)` netted out and
/// slipped through — then decode sliced both streams with the same `lens` and
/// silently handed r1 two of r2's quality bytes. Worse, `--verify` re-derived its
/// digests from that same wrong `lens` and reported success. The fix rejects the
/// first mis-sized record up front; this asserts it does, in both order modes.
#[test]
fn compensating_seq_qual_mismatch_is_rejected() {
    // Totals are 6 seq / 6 qual, but neither record has seq_len == qual_len.
    let input = b"@r1\nAAAA\n+\nII\n@r2\nCC\n+\nGGGG\n";

    let plain = compress(&input[..], &mut Vec::new(), Params::default());
    let reordered = compress(
        &input[..],
        &mut Vec::new(),
        Params {
            reorder: true,
            ..Params::default()
        },
    );

    assert!(
        matches!(plain, Err(Error::RecordLengthMismatch { seq: 4, qual: 2 })),
        "plain path must reject the first mis-sized record, got {plain:?}"
    );
    assert!(
        reordered.is_err(),
        "reorder path must reject it too, not silently misalign quality"
    );
}

#[test]
fn roundtrip_normalizes_plus() {
    let archive = compress_bytes(SAMPLE, Params::default());
    let mut fastq = Vec::new();
    decompress(&archive[..], &mut fastq, 1).expect("decompress");
    let expected = b"\
@SRR1.1 INST:1:FC:1:1101:1000:2000 length=8\n\
ACGTACGT\n\
+\n\
IIIIFFF#\n\
@SRR1.2 INST:1:FC:1:1101:1005:2050 length=8\n\
NNGGCCTA\n\
+\n\
###IIIFF\n";
    assert_eq!(fastq, expected);
}

/// #49: the raw definition line must round-trip byte-exactly — a trailing
/// separator or a tab must not be dropped or rewritten to a space.
const SEP_EDGE: &[u8] = b"\
@a:b:c:d \nACGTACGTACGT\n+\nIIIIIIIIIIII\n\
@id\tdesc\nTTGGCCAATTGG\n+\nFFFFFFFFFFFF\n\
@plain desc\nAACCGGTTAACC\n+\nJJJJJJJJJJJJ\n\
@bare\nACGTACGTACGT\n+\nKKKKKKKKKKKK\n";

#[test]
fn preserves_trailing_and_tab_separators_plain() {
    let archive = compress_bytes(SEP_EDGE, Params::default());
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(
        out, SEP_EDGE,
        "plain path must preserve headers byte-exactly"
    );
}

#[test]
fn preserves_trailing_and_tab_separators_reorder() {
    let archive = compress_bytes(
        SEP_EDGE,
        Params {
            reorder: true,
            ..Params::default()
        },
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    // Single-end reorder may change read order, so compare the record set.
    assert_eq!(
        record_set(&out),
        record_set(SEP_EDGE),
        "reorder path must preserve headers byte-exactly"
    );
}

#[test]
fn preserves_separators_multi_file() {
    let r1: &[u8] = b"@a:b \nACGT\n+\nIIII\n@c\td\nTTTT\n+\nFFFF\n";
    let r2: &[u8] = b"@a:b \nGGGG\n+\nJJJJ\n@c\td\nCCCC\n+\nKKKK\n";
    let readers: Vec<Box<dyn io::Read + Send>> = vec![Box::new(r1), Box::new(r2)];
    let mut archive = Vec::new();
    compress_multi(readers, &mut archive, Params::default()).unwrap();
    let (mut m1, mut m2) = (Vec::new(), Vec::new());
    decompress_split(&archive[..], &mut [&mut m1, &mut m2], 1).unwrap();
    assert_eq!(m1, r1, "mate 1 header not byte-exact");
    assert_eq!(m2, r2, "mate 2 header not byte-exact");
}

#[test]
fn classify_header_reads_platform_from_name_grammar() {
    // Illumina Casava 1.8 name + description.
    assert_eq!(
        classify_header(b"M01234:12:000-ABC:1:1101:1234:5678 1:N:0:ATCACG"),
        Platform::Illumina
    );
    // Older Illumina with #index/mate.
    assert_eq!(
        classify_header(b"HWUSI:2:3:4:5#ATCACG/1"),
        Platform::Illumina
    );
    // Nanopore: UUID name, and the runid= description tag alone.
    assert_eq!(
        classify_header(b"1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d runid=x read=1 ch=100"),
        Platform::Nanopore
    );
    assert_eq!(
        classify_header(b"anything runid=deadbeef ch=42"),
        Platform::Nanopore
    );
    // PacBio movie/zmw/ccs.
    assert_eq!(
        classify_header(b"m64011_190228_190319/1001/ccs"),
        Platform::PacBio
    );
    // MGI/BGI V-prefixed flowcell.
    assert_eq!(
        classify_header(b"V300026399L1C001R0010000001/1"),
        Platform::MgiBgi
    );
    // Bare names match nothing.
    assert_eq!(classify_header(b"read_42"), Platform::Unknown);
    assert_eq!(classify_header(b"SRR1.1"), Platform::Unknown);
}

#[test]
fn platform_is_detected_stored_and_reported() {
    let ont = b"\
@1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d runid=x read=1 ch=100\n\
ACGTACGT\n+\nIIIIFFF#\n";
    let archive = compress_bytes(ont, Params::default());
    assert_eq!(
        inspect(std::io::Cursor::new(&archive)).unwrap().platform,
        Platform::Nanopore
    );
    // peek reads it from the header flags too.
    assert_eq!(peek(&archive[..]).unwrap().platform, Platform::Nanopore);
}

#[test]
fn platform_override_forces_recorded_value() {
    // Bare names would auto-detect Unknown; the override wins.
    let params = Params {
        platform: Some(Platform::PacBio),
        ..Params::default()
    };
    let archive = compress_bytes(SAMPLE, params);
    assert_eq!(
        inspect(std::io::Cursor::new(&archive)).unwrap().platform,
        Platform::PacBio
    );
}

#[test]
fn platform_survives_paired_and_reorder_paths() {
    // Paired input flows through compress_multi's streaming drive path.
    let mate = |m: &str| {
        format!("@M01234:12:000-ABC:1:1101:1000:2000 {m}:N:0:ATCACG\nACGT\n+\nIIII\n").into_bytes()
    };
    let (r1, r2) = (mate("1"), mate("2"));
    let mut archive = Vec::new();
    let readers: Vec<Box<dyn Read + Send>> =
        vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
    compress_multi(readers, &mut archive, Params::default()).unwrap();
    assert_eq!(peek(&archive[..]).unwrap().platform, Platform::Illumina);

    // Reorder (global-cluster) path stores it in its own header.
    let ont = b"\
@1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d runid=x read=1 ch=100\n\
ACGTACGTACGT\n+\nIIIIFFF#IIII\n";
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let archive = compress_bytes(ont, params);
    assert_eq!(peek(&archive[..]).unwrap().platform, Platform::Nanopore);
}

#[test]
fn illumina_reorder_without_reference_roundtrips() {
    // Regression: the platform tag used to live in flags bits 5-7, where
    // `Platform::Illumina` (code 1) produced exactly 0x20 == FLAG_GLOBAL_REFERENCE.
    // An Illumina archive in the global-reorder layout that did not adopt a
    // reference (`use_reference` is false whenever the reference does not
    // strictly win — e.g. on non-redundant data like this) still advertised one,
    // and decode died on `Corrupt { what: "reorder global reference" }`. The tag
    // now has its own header byte; the platform must never imply a flag.
    let mut input = Vec::new();
    let mut state: u32 = 12345;
    for i in 0..64 {
        let mut seq = Vec::new();
        for _ in 0..60 {
            // Deterministic LCG: non-redundant bases so no reference wins.
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            seq.push(b"ACGT"[(state >> 16) as usize % 4]);
        }
        input.extend_from_slice(
            format!("@M01234:12:000-ABC:1:1101:{i}:2000 1:N:0:ATCACG\n").as_bytes(),
        );
        input.extend_from_slice(&seq);
        input.extend_from_slice(b"\n+\n");
        input.extend_from_slice(&[b'I'; 60]);
        input.push(b'\n');
    }
    // keep_order so the exact input is reconstructible; the global-reorder
    // layout (and thus the flags/platform header path) is exercised either way.
    let params = Params {
        reorder: true,
        keep_order: true,
        ..Params::default()
    };
    let archive = compress_bytes(&input, params);
    // The platform must survive, and the archive must decode.
    assert_eq!(peek(&archive[..]).unwrap().platform, Platform::Illumina);
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).expect("decompress Illumina+reorder archive");
    assert_eq!(out, input, "Illumina + reorder must round-trip");
}

fn make_reads(tag: &str, n: usize) -> Vec<u8> {
    let mut v = Vec::new();
    for i in 0..n {
        v.extend_from_slice(format!("@r.{i} {tag}\nACGT\n+\nIIII\n").as_bytes());
    }
    v
}

#[test]
fn content_stats_and_metadata() {
    // Two reads of different lengths with hand-countable content.
    let input = b"@r0 x\nACGTACGT\n+\nIIIIIIII\n@r1 x\nACGTN\n+\n!!!!!\n";
    let archive = compress_bytes(
        input,
        Params {
            threads: 1,
            ..Params::default()
        },
    );

    // Metadata is read from the header + footer without decoding.
    let info = inspect(io::Cursor::new(&archive)).unwrap();
    assert_eq!(
        info.format_version,
        (u16::from(FORMAT_MAJOR) << 8) | u16::from(FORMAT_MINOR)
    );
    assert!(info.whole_file_crc.is_some(), "plain layout stores the CRC");
    assert_eq!(info.reads, 2);

    // Content stats require a full decode.
    let cs = content_stats(&archive[..], 1).unwrap();
    assert_eq!(cs.reads, 2);
    assert_eq!(cs.bases, 13);
    assert_eq!((cs.min_len, cs.max_len), (5, 8));
    assert!(!cs.fixed_length());
    assert_eq!((cs.a, cs.c, cs.g, cs.t, cs.n, cs.other), (3, 3, 3, 3, 1, 0));
    assert_eq!(cs.gc_fraction(), Some(0.5));
    assert_eq!(cs.mean_len(), Some(6.5));
    // 'I' is Phred 40, '!' is Phred 0.
    assert_eq!(cs.qual_hist[40], 8);
    assert_eq!(cs.qual_hist[0], 5);
    assert_eq!(cs.qual_sum, 8 * 40);
    assert!((cs.mean_quality().unwrap() - 320.0 / 13.0).abs() < 1e-9);

    // Thread count must not change the summary (decode is deterministic).
    assert_eq!(content_stats(&archive[..], 4).unwrap(), cs);
}

#[test]
fn paired_roundtrip_splits() {
    let r1 = make_reads("a", 3);
    let r2 = make_reads("b", 3);
    let mut archive = Vec::new();
    let readers: Vec<Box<dyn Read + Send>> =
        vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
    let s = compress_multi(readers, &mut archive, Params::default()).unwrap();
    assert_eq!(s.reads, 6);
    assert_eq!(peek(&archive[..]).unwrap().group_size, 2);

    let (mut o1, mut o2) = (Vec::new(), Vec::new());
    {
        let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
        decompress_split(&archive[..], &mut outs, 1).unwrap();
    }
    assert_eq!(o1, r1);
    assert_eq!(o2, r2);
}

#[test]
fn single_cell_four_way_roundtrip() {
    // 10x-style: R1, R2, I1, I2.
    let files: Vec<Vec<u8>> = ["R1", "R2", "I1", "I2"]
        .iter()
        .map(|t| make_reads(t, 5))
        .collect();
    let mut archive = Vec::new();
    let readers: Vec<Box<dyn Read + Send>> = files
        .iter()
        .map(|f| Box::new(&f[..]) as Box<dyn Read + Send>)
        .collect();
    compress_multi(readers, &mut archive, Params::default()).unwrap();
    assert_eq!(peek(&archive[..]).unwrap().group_size, 4);

    let mut outs: Vec<Vec<u8>> = vec![Vec::new(); 4];
    decompress_split(&archive[..], &mut outs, 1).unwrap();
    assert_eq!(outs, files);
}

#[test]
fn paired_split_spans_multiple_blocks() {
    // De-interleaving is a block-local `i % g`, correct only because every block
    // starts on member 0. A single-block archive can't exercise that; force many
    // blocks with a tiny target so the split must re-anchor per block.
    let r1 = make_reads("a", 25);
    let r2 = make_reads("b", 25);
    let mut archive = Vec::new();
    let readers: Vec<Box<dyn Read + Send>> =
        vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
    compress_multi(
        readers,
        &mut archive,
        Params {
            block_reads: 4,
            ..Params::default()
        },
    )
    .unwrap();
    assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
    assert!(
        inspect(io::Cursor::new(&archive[..])).unwrap().blocks > 1,
        "test must span multiple blocks to cover per-block member re-anchoring"
    );

    let (mut o1, mut o2) = (Vec::new(), Vec::new());
    {
        let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
        decompress_split(&archive[..], &mut outs, 1).unwrap();
    }
    assert_eq!(o1, r1, "R1 must reassemble across block boundaries");
    assert_eq!(o2, r2, "R2 must reassemble across block boundaries");
}

#[test]
fn grouped_block_rejects_partial_spot() {
    // The whole-spots invariant (a block's read count is a multiple of g) is
    // enforced at encode but nowhere recorded on disk, so `decode_block_group`
    // guards it: a block whose count is not a multiple of g would otherwise
    // silently misroute the trailing partial spot. Build a valid 3-read block
    // (a plain single-end archive) and feed it to the grouped splitter as g = 2.
    let archive = compress_bytes(&make_reads("x", 3), Params::default());
    let len_off = HEADER_LEN + BLOCK_MAGIC.len();
    let payload_len =
        u64::from_le_bytes(archive[len_off..len_off + 8].try_into().unwrap()) as usize;
    let start = HEADER_LEN + FRAME_HEAD_LEN;
    let payload = &archive[start..start + payload_len];

    // g that divides the count still decodes; g = 2 does not divide 3 and errors.
    assert!(decode_block_group(payload, 1, None).is_ok());
    assert!(decode_block_group(payload, 3, None).is_ok());
    let err = decode_block_group(payload, 2, None).unwrap_err();
    assert!(
        matches!(err, Error::Malformed(_)),
        "partial-spot block must be rejected, got {err:?}"
    );
}

#[test]
fn grouped_archive_streams_interleaved() {
    let r1 = b"@r.1 a\nACGT\n+\nIIII\n";
    let r2 = b"@r.1 b\nGGGG\n+\n####\n";
    let mut archive = Vec::new();
    let readers: Vec<Box<dyn Read + Send>> =
        vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
    compress_multi(readers, &mut archive, Params::default()).unwrap();
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, b"@r.1 a\nACGT\n+\nIIII\n@r.1 b\nGGGG\n+\n####\n");
}

// Two paired spots, mates interleaved on one stream with /1 /2 names.
const INTERLEAVED: &[u8] = b"\
@s1/1\nAAAA\n+\nIIII\n\
@s1/2\nTTTT\n+\nFFFF\n\
@s2/1\nCCCC\n+\nIIII\n\
@s2/2\nGGGG\n+\nFFFF\n";

#[test]
fn interleaved_stream_forces_pairing_and_splits() {
    let mut archive = Vec::new();
    let s = compress_interleaved(INTERLEAVED, &mut archive, Params::default(), 2).unwrap();
    assert_eq!(s.reads, 4);
    assert_eq!(s.group_size, 2);
    assert_eq!(peek(&archive[..]).unwrap().group_size, 2);

    let (mut o1, mut o2) = (Vec::new(), Vec::new());
    {
        let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
        decompress_split(&archive[..], &mut outs, 1).unwrap();
    }
    assert_eq!(o1, b"@s1/1\nAAAA\n+\nIIII\n@s2/1\nCCCC\n+\nIIII\n");
    assert_eq!(o2, b"@s1/2\nTTTT\n+\nFFFF\n@s2/2\nGGGG\n+\nFFFF\n");
}

#[test]
fn interleaved_odd_count_errors() {
    let mut truncated = INTERLEAVED.to_vec();
    truncated.extend_from_slice(b"@s3/1\nACGT\n+\nIIII\n"); // dangling mate
    let err = compress_interleaved(&truncated[..], &mut Vec::new(), Params::default(), 2);
    assert!(matches!(err, Err(Error::Malformed(_))));
}

#[test]
fn auto_detects_interleaved_pairing() {
    let mut archive = Vec::new();
    let s = compress_auto(INTERLEAVED, &mut archive, Params::default()).unwrap();
    assert_eq!(
        s.group_size, 2,
        "paired /1 /2 names should auto-detect as paired"
    );
    assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
}

#[test]
fn auto_leaves_single_end_ungrouped() {
    // Distinct, unpaired names must not be mistaken for mates.
    let single = make_reads("x", 6);
    let mut archive = Vec::new();
    let s = compress_auto(&single[..], &mut archive, Params::default()).unwrap();
    assert_eq!(s.group_size, 1);

    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, single);
}

#[test]
fn unequal_mate_counts_error() {
    let r1 = make_reads("a", 2);
    let r2 = make_reads("b", 1);
    let readers: Vec<Box<dyn Read + Send>> =
        vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
    let err = compress_multi(readers, &mut Vec::new(), Params::default());
    assert!(matches!(err, Err(Error::Malformed(_))));
}

#[test]
fn split_count_mismatch_errors() {
    let archive = compress_bytes(SAMPLE, Params::default()); // group_size 1
    let mut outs: Vec<Vec<u8>> = vec![Vec::new(); 2];
    let err = decompress_split(&archive[..], &mut outs, 1);
    assert!(matches!(err, Err(Error::Malformed(_))));
}

#[test]
fn inspect_reports_streams() {
    let archive = compress_bytes(SAMPLE, Params::default());
    let info = inspect(io::Cursor::new(&archive[..])).expect("inspect");
    assert_eq!(info.reads, 2);
    assert_eq!(info.blocks, 1);
    assert_eq!(info.group_size, 1);
    assert!(info.plus_normalized);
    assert!(info.names_bytes > 0 && info.seq_bytes > 0 && info.qual_bytes > 0);
}

// Concatenate every Nth record line (line 4 = quality, line 2 = sequence)
// across a FASTQ byte stream, in order.
fn record_line(fastq: &[u8], which: usize) -> Vec<u8> {
    fastq
        .split(|&b| b == b'\n')
        .enumerate()
        .filter(|(i, l)| i % 4 == which && !l.is_empty())
        .flat_map(|(_, l)| l.iter().copied())
        .collect()
}

#[test]
fn lossy_binning_roundtrips_and_reports_tag() {
    for (bin, tag) in [
        (QualityBinning::Bin8, 1u8),
        (QualityBinning::Bin4, 2),
        (QualityBinning::Bin2, 3),
    ] {
        let params = Params {
            quality_binning: bin,
            ..Params::default()
        };
        let archive = compress_bytes(SAMPLE, params);

        // The header tag round-trips through inspect.
        assert_eq!(
            inspect(io::Cursor::new(&archive[..]))
                .expect("inspect")
                .quality_binning,
            tag,
            "info tag for {bin:?}"
        );

        let mut fastq = Vec::new();
        decompress(&archive[..], &mut fastq, 1).expect("decompress");

        // Lossy contract: recovered qualities equal the input qualities passed
        // through the same bin table; bases survive exactly.
        let want: Vec<u8> = record_line(SAMPLE, 3)
            .iter()
            .map(|&b| bin.apply(b))
            .collect();
        assert_eq!(record_line(&fastq, 3), want, "binned qualities for {bin:?}");
        assert_eq!(
            record_line(&fastq, 1),
            record_line(SAMPLE, 1),
            "bases must be exact for {bin:?}"
        );
    }
}

#[test]
fn reorder_with_lossy_binning_roundtrips() {
    // Regression: the reorder whole-output digest must be folded over the
    // *stored* (post-binning) quality, not the original input. Folding the
    // original made a `--order any --quality-bin` archive fail its own
    // output-digest check on decode — it recovered the data but returned
    // Err(Corrupt { "reorder output digest" }). Use quality that actually
    // shifts under binning so a digest-over-original would mismatch, and a
    // revcomp pair so clustering exercises the un-flip path.
    let a = b"ACGTTTGACCGATTGCAACGT";
    let ra = fqxv_reorder::revcomp(a);
    let ql: Vec<u8> = (0..a.len()).map(|i| b'!' + (i as u8 * 3 % 40)).collect();
    let read = |i: u32| -> Vec<u8> {
        let s = if i.is_multiple_of(2) {
            a.to_vec()
        } else {
            ra.clone()
        };
        let mut rec = format!("@read.{i}\n").into_bytes();
        rec.extend_from_slice(&s);
        rec.extend_from_slice(b"\n+\n");
        rec
    };
    let mut input = Vec::new();
    for i in 0..40u32 {
        let mut rec = read(i);
        rec.extend_from_slice(&ql);
        rec.push(b'\n');
        input.extend_from_slice(&rec);
    }
    for bin in [
        QualityBinning::Bin8,
        QualityBinning::Bin4,
        QualityBinning::Bin2,
    ] {
        let params = Params {
            reorder: true,
            quality_binning: bin,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        assert_eq!(
            archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
            FLAG_GLOBAL_REORDER
        );
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1)
            .unwrap_or_else(|e| panic!("reorder + {bin:?} decode failed: {e:?}"));
        // Order-independent (reorder permutes reads): recovered records equal
        // the input with quality passed through the bin table.
        let binned: Vec<u8> = ql.iter().map(|&q| bin.apply(q)).collect();
        let mut want: Vec<Vec<u8>> = (0..40u32)
            .map(|i| {
                let mut rec = read(i);
                rec.extend_from_slice(&binned);
                rec
            })
            .collect();
        want.sort();
        assert_eq!(record_set(&out), want, "content for reorder + {bin:?}");
    }
}

#[test]
fn lossless_default_reports_zero_tag() {
    let archive = compress_bytes(SAMPLE, Params::default());
    assert_eq!(
        inspect(io::Cursor::new(&archive[..]))
            .unwrap()
            .quality_binning,
        0
    );
}

fn dup_rich_input(keep_order_marker: char) -> Vec<u8> {
    // Duplicate-rich single-end reads, including a reverse-complement pair so
    // clustering flips a read (exercises the un-flip path).
    let a = b"ACGTTTGACCGATTGCAACGT";
    let ra = fqxv_reorder::revcomp(a);
    let mut input = Vec::new();
    for i in 0..40u32 {
        let s = match i % 3 {
            0 => a.to_vec(),
            1 => ra.clone(),
            _ => b"TTTTGGGGCCCCAAAATTTTG".to_vec(),
        };
        input.extend_from_slice(format!("@read.{i} {keep_order_marker}\n").as_bytes());
        input.extend_from_slice(&s);
        input.extend_from_slice(format!("\n+\n{}\n", "I".repeat(s.len())).as_bytes());
    }
    input
}

fn record_set(fastq: &[u8]) -> Vec<Vec<u8>> {
    let lines: Vec<&[u8]> = fastq.split(|&b| b == b'\n').collect();
    let mut recs: Vec<Vec<u8>> = lines
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| c.join(&b"\n"[..]))
        .collect();
    recs.sort();
    recs
}

#[test]
fn reorder_keep_order_is_byte_exact() {
    let input = dup_rich_input('d');
    let params = Params {
        reorder: true,
        keep_order: true,
        rescue: false,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    assert_eq!(
        inspect(io::Cursor::new(&archive[..])).unwrap().group_size,
        1
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, input, "reorder --keep-order must be byte-exact");
}

#[test]
fn long_reads_skip_reorder() {
    // Requesting reorder on long-read data must auto-fall-back to the
    // non-reorder layout (smaller and far cheaper there): the global-reorder
    // flag stays clear, and the archive still round-trips byte-exact.
    let seq: Vec<u8> = (0..800u32).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
    let qual = vec![b'I'; seq.len()];
    let mut input = Vec::new();
    for i in 0..30u32 {
        write_record(&mut input, format!("read.{i}").as_bytes(), &seq, &qual);
    }
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    assert_eq!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
        0,
        "reorder must be skipped (flag clear) for long-read data"
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, input, "long-read fallback must be byte-exact");
}

/// The leading sequence-stream method byte of the archive's first block
/// (`SEQ_METHOD_ORDERK` or `SEQ_METHOD_OVERLAP`). Parses the plain-layout frame
/// and block payload directly — the same offsets [`decode_block_parts`] uses — so
/// a test can assert *which* sequence codec the container kept, rather than
/// assuming it from the input shape.
fn first_block_seq_method(archive: &[u8]) -> u8 {
    // The shared-reference layout writes a framed reference (`[4 len][4 crc][bytes]`)
    // between the header and the first block, so skip it when present.
    let mut block_start = HEADER_LEN;
    if archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REFERENCE != 0 {
        let ref_len =
            u32::from_le_bytes(archive[HEADER_LEN..HEADER_LEN + 4].try_into().unwrap()) as usize;
        block_start += 4 + CRC_LEN + ref_len;
    }
    let len_off = block_start + BLOCK_MAGIC.len();
    let plen = u64::from_le_bytes(archive[len_off..len_off + 8].try_into().unwrap()) as usize;
    let payload = &archive[block_start + FRAME_HEAD_LEN..block_start + FRAME_HEAD_LEN + plen];
    // payload: [24 stream digests][4 n_reads][4 names_len][names][4 seq_len][seq…]
    let names_len = u32::from_le_bytes(payload[28..32].try_into().unwrap()) as usize;
    payload[32 + names_len + 4] // first byte of the seq stream is the method tag
}

#[test]
fn long_read_path_tries_overlap_but_keeps_order_k_roundtrip_and_determinism() {
    // 600 bp reads clear the long-read gate (`is_long_read`, mean > 500), so the
    // container *tries* the overlap codec and order-k and keeps the smaller — and
    // at this scale order-k wins, so the block is stored order-k. (Overlap only
    // beats order-k at realistic multi-kb reads and high depth; on a tiny genome
    // the assembly/consensus overhead loses.) This test guards that long-read
    // *selection* path: it round-trips byte-exact and is thread-count invariant.
    // The overlap *decode* dispatch is covered separately by
    // `container_decodes_overlap_tagged_sequence_stream`.
    let genome: Vec<u8> = (0..1500u32)
        .map(|i| b"ACGT"[((i.wrapping_mul(2_654_435_761) >> 13) & 3) as usize])
        .collect();
    let mut input = Vec::new();
    for i in 0..24u32 {
        let start = (i as usize * 37) % (genome.len() - 600);
        let mut s = genome[start..start + 600].to_vec();
        // A substitution the aligner must code, and an N, so the exception path
        // is exercised end to end through the container.
        s[100] = b"ACGT"[((i as usize) + 1) & 3];
        s[400] = b'N';
        let qual = vec![b'I'; s.len()];
        write_record(&mut input, format!("read.{i}").as_bytes(), &s, &qual);
    }

    let mut archives = Vec::new();
    for threads in [1usize, 4] {
        let params = Params {
            threads,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        archives.push(archive);
    }
    assert_eq!(
        first_block_seq_method(&archives[0]),
        SEQ_METHOD_ORDERK,
        "order-k must win the size race at this scale; if overlap ever wins here, \
         this test's premise (and name) is stale"
    );
    assert_eq!(
        archives[0], archives[1],
        "long-read selection output must not vary by thread count"
    );
    let mut out = Vec::new();
    decompress(&archives[0][..], &mut out, 4).unwrap();
    assert_eq!(out, input, "long-read archive must round-trip byte-exact");
}

#[test]
fn container_decodes_overlap_tagged_sequence_stream() {
    // The container's sequence path dispatches on a leading method byte; the
    // overlap branch (`SEQ_METHOD_OVERLAP` -> `fqxv_lroverlap::decode`) is not
    // reached by the compress path on synthetic data, because order-k wins the
    // size race there. Cover that decode dispatch directly: build a genuine
    // overlap-coded stream, tag it, and round-trip it through the container's
    // `decode_sequence_stream`. Long reads tiling a genome so the codec forms
    // real overlaps and does substantive assembly/consensus work.
    let genome: Vec<u8> = (0..3000u32)
        .map(|i| b"ACGT"[((i.wrapping_mul(2_654_435_761) >> 13) & 3) as usize])
        .collect();
    let (read_len, step) = (1000usize, 100usize);
    let n = (genome.len() - read_len) / step;
    let mut lens = Vec::with_capacity(n);
    let mut seq = Vec::new();
    for i in 0..n {
        let mut s = genome[i * step..i * step + read_len].to_vec();
        s[50] = b"ACGT"[(i + 1) & 3]; // a substitution the aligner must code
        lens.push(s.len() as u32);
        seq.extend_from_slice(&s);
    }

    let coded = fqxv_lroverlap::encode(&lens, &seq, &fqxv_lroverlap::EncodeOpts::default())
        .expect("overlap encode");
    let mut stream = vec![SEQ_METHOD_OVERLAP];
    stream.extend_from_slice(&coded);

    let (dlens, dseq) = decode_sequence_stream(&stream, None).expect("overlap dispatch decode");
    assert_eq!(dlens, lens, "overlap decode must restore per-read lengths");
    assert_eq!(
        dseq, seq,
        "overlap decode must restore the bases byte-exact"
    );
}

/// Deep-tiling long-read fixture: `n` reads of ~`read_len` bp tiling a genome at a
/// short step (high coverage, so the overlap codec forms a contig), each carrying
/// ~1% random substitution error — the sequencing noise that pollutes the order-k
/// context model but is voted out by the consensus, which is what makes the shared
/// reference win the size race (a clean synthetic genome is too easy for order-k).
/// One `N` per read exercises the non-ACGT exception path. Returns FASTQ.
/// `err_period` sets the per-base substitution rate (`1/err_period`): 300 is
/// HiFi-like (~0.3%), a small value is ONT-like and noisy enough that reads stop
/// collapsing onto a shared consensus.
fn deep_longread_fastq(
    n: usize,
    read_len: usize,
    step: usize,
    genome_len: usize,
    err_period: u64,
) -> Vec<u8> {
    // splitmix64, used for both the genome and the per-read error, so the genome is
    // near-uniform ACGT (order-0 ≈ 2 bits/base) — a biased genome would let even a
    // low-order model compress it below what the reference-coded stream achieves.
    let rng = |x: &mut u64| {
        *x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = *x;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    let genome: Vec<u8> = {
        let mut g = 0x5eed_1234u64;
        (0..genome_len)
            .map(|_| b"ACGT"[(rng(&mut g) % 4) as usize])
            .collect()
    };
    let mut input = Vec::new();
    for i in 0..n {
        let start = (i * step) % (genome.len() - read_len);
        let mut s = genome[start..start + read_len].to_vec();
        let mut state = 0x1234_5678u64.wrapping_add(i as u64);
        for b in s.iter_mut() {
            if rng(&mut state) % err_period == 0 {
                *b = b"ACGT"[(rng(&mut state) % 4) as usize];
            }
        }
        s[read_len / 2] = b'N'; // exercise the non-ACGT exception path
        let qual = vec![b'I'; s.len()];
        write_record(&mut input, format!("read.{i}").as_bytes(), &s, &qual);
    }
    input
}

#[test]
fn shared_reference_gate_measures_against_the_plain_layout() {
    // Issue #184, pinned with the ONT numbers that actually regressed. The gate must
    // compare the shared-reference layout against the plain layout it falls back to
    // — per block, the smaller of the overlap codec and order-k — not against
    // order-k alone.
    //
    // These constants are the measured `ecoli_ont` (DRR205413) totals: a 4.37 MB
    // frame buys a sequence stream only 1.58 MB smaller than the plain layout's, so
    // adopting it inflates the archive by ~2.8 MB. Crucially it still beats order-k
    // (~1.8 b/base on ONT), which is exactly why the old bar let it through.
    let (frame, shared, plain, order_k) = (4_369_101, 48_292_740, 49_868_432, 67_768_843);

    assert!(
        frame + shared > plain,
        "fixture must describe a net loss against the plain layout"
    );
    assert!(
        frame + shared < order_k,
        "fixture must still clear the old order-k bar, or it would not be a regression"
    );
    assert!(
        !adopt_shared_reference(frame, shared, plain),
        "a reference that loses to the per-block overlap layout must be rejected"
    );

    // HiFi inverts the tradeoff: a compact consensus (1.36 MB frame) buys a 7.20 MB
    // smaller sequence stream, so there the reference must still be adopted.
    let (hifi_frame, hifi_plain, hifi_saving) = (1_356_344, 100_000_000, 7_201_416);
    assert!(
        adopt_shared_reference(hifi_frame, hifi_plain - hifi_saving, hifi_plain),
        "a reference that pays for itself must still be adopted"
    );

    // A tie is a loss: an equal-size archive plus a frame is pure overhead.
    let tie = 1_000;
    assert!(!adopt_shared_reference(0, tie, tie));
}

#[test]
fn longread_routing_roundtrips_and_never_loses_to_the_plain_layout() {
    // Broad guard on the long-read routing: whichever layout the whole-file gate
    // picks must decode byte-exact and must not produce a larger archive than the
    // plain layout.
    //
    // This does NOT pin issue #184 — verified by mutation: forcing the gate to
    // always adopt still passes here. Reproducing that regression needs the
    // whole-file assembly to collapse *worse* than the per-block assemblies, which
    // only happens at real coverage and read counts; on a fixture this small the
    // frame is cheap and always pays off. The gate arithmetic is pinned by
    // `shared_reference_gate_measures_against_the_plain_layout` instead, and the
    // end-to-end behaviour by the benchmark suite.
    for (label, err_period) in [("hifi-like", 300u64), ("ont-like", 12u64)] {
        let input = deep_longread_fastq(120, 2000, 20, 4000, err_period);
        // Order-0, small blocks: same fixture regime as the round-trip test above,
        // where the order-k model is weak enough not to memorize the tiny genome.
        let params = Params {
            threads: 4,
            seq_order: 0,
            block_reads: 60,
            ..Params::default()
        };

        let mut routed = Vec::new();
        compress_auto(&input[..], &mut routed, params).expect("compress");
        let mut plain = Vec::new();
        compress_buffered_plain(&input, &mut plain, params, 1).expect("plain layout");

        assert!(
            routed.len() <= plain.len(),
            "{label}: the shared-reference path must never produce a larger archive \
             than the plain layout it falls back to ({} vs {} bytes)",
            routed.len(),
            plain.len(),
        );

        // Whichever layout the gate chose must still decode byte-exact.
        let mut out = Vec::new();
        decompress(&routed[..], &mut out, 4).expect("decompress");
        assert_eq!(out, input, "{label}: archive must round-trip byte-exact");
    }
}

#[test]
fn shared_reference_longread_roundtrips_across_blocks() {
    // The whole-file shared-reference layout (issue #168): deep-tiling long reads
    // split across several blocks are coded against ONE reference stored between the
    // header and the first block. This guards the whole wiring — reference frame
    // written and read back, footer offsets that start past it, and every block's
    // sequence decoded against the shared frame — with a byte-exact round trip,
    // thread-count determinism, and a full verify.
    let input = deep_longread_fastq(120, 2000, 20, 4000, 300);

    // Order-0 keeps the order-k baseline weak enough (~2 bits/base) that the
    // reference-coded stream wins the size race on this compact synthetic fixture;
    // at default order-11 the model would memorize such a tiny genome outright,
    // which is why the plain long-read tests see order-k win. Small blocks force
    // multiple row groups over the one shared reference, exercising the cross-block
    // frame wiring.
    let mk = |threads: usize| Params {
        threads,
        seq_order: 0,
        block_reads: 60,
        ..Params::default()
    };
    let mut archive = Vec::new();
    let stats = compress_auto(&input[..], &mut archive, mk(4)).expect("compress");
    assert!(stats.blocks > 1, "fixture must span multiple blocks");

    // The reference must have been adopted: the header flag and feature bit are set.
    assert_ne!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REFERENCE,
        0,
        "shared reference must be adopted for deep-tiling long reads"
    );
    let features = u64::from_le_bytes(archive[6..14].try_into().unwrap());
    assert_ne!(
        features & crate::feature::GLOBAL_REFERENCE,
        0,
        "the reference feature bit must be set so pre-feature readers refuse the archive"
    );

    // At least one block actually used the reference-coded method (else the frame
    // would be wasted overhead the whole-file gate should have rejected).
    assert_eq!(
        first_block_seq_method(&archive),
        SEQ_METHOD_OVERLAP_REF,
        "the first deep block should code against the shared reference"
    );

    // Byte-exact round trip and full structural + content verification.
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 4).expect("decompress");
    assert_eq!(
        out, input,
        "shared-reference archive must round-trip byte-exact"
    );
    verify_roundtrip(io::Cursor::new(&archive), 4).expect("shared-reference archive must verify");

    // Thread-count determinism: the reference build and per-block coding are pure.
    let mut archive1 = Vec::new();
    compress_auto(&input[..], &mut archive1, mk(1)).expect("compress 1-thread");
    assert_eq!(
        archive, archive1,
        "shared-reference output must be byte-identical regardless of thread count"
    );
}

#[test]
fn sketch_for_platform_picks_hifi_only_for_pacbio() {
    // PacBio's low error rate earns the sparse HiFi sketch at either coverage;
    // every other platform (including a mis- or undetected one) falls back to the
    // dense ONT sketch, which also works on HiFi and never misses ONT overlaps.
    for ctx in [SeedContext::WholeFile, SeedContext::PerBlock] {
        assert_eq!(
            sketch_for(Platform::PacBio, ctx),
            fqxv_lroverlap::Sketch::hifi()
        );
        for p in [
            Platform::Nanopore,
            Platform::Unknown,
            Platform::Illumina,
            Platform::MgiBgi,
        ] {
            let s = sketch_for(p, ctx);
            let ont = fqxv_lroverlap::Sketch::ont();
            assert_eq!((s.w, s.k), (ont.w, ont.k), "{p:?}/{ctx:?} keeps ONT (w, k)");
        }
    }
}

#[test]
fn ont_seeding_scheme_follows_the_index_coverage() {
    // Issue #184 follow-up. Syncmers conserve anchors better at ONT error rates
    // but only pay off at depth, so the two long-read indexes want opposite
    // schemes: the whole-file reference sees full coverage (syncmers, 1.280 vs
    // 1.416 b/base measured on ecoli_ont), one block sees a fraction of it
    // (minimizers, 1.243 vs 1.517). Using syncmers for both — the state before
    // this split — cost 2.79 MB on that archive.
    for p in [Platform::Nanopore, Platform::Unknown, Platform::MgiBgi] {
        assert_eq!(
            sketch_for(p, SeedContext::WholeFile).scheme,
            fqxv_lroverlap::SeedScheme::Syncmer,
            "{p:?}: the whole-file reference sees full coverage"
        );
        assert_eq!(
            sketch_for(p, SeedContext::PerBlock).scheme,
            fqxv_lroverlap::SeedScheme::Minimizer,
            "{p:?}: a per-block index sees only its share of the coverage"
        );
    }
    // Anchor density is 2/(w + 1) and specificity follows k, so the split changes
    // only WHICH positions are selected — not how many, and not the index cost.
    let (whole, block) = (
        sketch_for(Platform::Nanopore, SeedContext::WholeFile),
        sketch_for(Platform::Nanopore, SeedContext::PerBlock),
    );
    assert_eq!((whole.w, whole.k), (block.w, block.k));
    // PacBio is deliberately unsplit: at <1% error nearly every k-mer survives,
    // so minimizers are already near-optimal at both coverages.
    assert_eq!(
        sketch_for(Platform::PacBio, SeedContext::WholeFile),
        sketch_for(Platform::PacBio, SeedContext::PerBlock),
    );
}

#[test]
fn pacbio_long_reads_use_hifi_sketch_and_roundtrip() {
    // Long PacBio reads route the overlap codec through the HiFi sketch (via
    // `sketch_for`) — a path the ONT default never exercised before. Guard that
    // the container detects PacBio, drives the HiFi-sketched overlap encode, and
    // round-trips byte-exact whichever codec ends up smaller. 1000 bp reads
    // tiling a 3 kb genome at high depth give the overlap codec real overlaps to
    // work with, the same shape as `container_decodes_overlap_tagged_sequence_stream`.
    let genome: Vec<u8> = (0..3000u32)
        .map(|i| b"ACGT"[((i.wrapping_mul(2_654_435_761) >> 13) & 3) as usize])
        .collect();
    let (read_len, step) = (1000usize, 100usize);
    let n = (genome.len() - read_len) / step;
    let mut input = Vec::new();
    for i in 0..n {
        let mut s = genome[i * step..i * step + read_len].to_vec();
        s[50] = b"ACGT"[(i + 1) & 3]; // a substitution the aligner must code
        let qual = vec![b'~'; s.len()];
        // PacBio `movie/zmw/ccs` read name so the platform detector fires.
        write_record(&mut input, format!("m64012/{i}/ccs").as_bytes(), &s, &qual);
    }

    let mut archive = Vec::new();
    compress(&input[..], &mut archive, Params::default()).unwrap();
    assert_eq!(
        peek(&archive[..]).unwrap().platform,
        Platform::PacBio,
        "container must detect PacBio from movie/zmw/ccs names"
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 4).unwrap();
    assert_eq!(
        out, input,
        "HiFi-sketched long-read archive must round-trip byte-exact"
    );
}

#[test]
fn reorder_free_preserves_records_as_a_set() {
    let input = dup_rich_input('e');
    let params = Params {
        reorder: true,
        keep_order: false,
        rescue: false,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(record_set(&out), record_set(&input));
}

#[test]
fn reorder_rescue_preserves_records_as_a_set() {
    // The literal-rescue sequence codec must round-trip through the container
    // (decode auto-detects the version byte). Multi-thread to exercise the
    // parallel per-block encode path.
    let input = dup_rich_input('e');
    for threads in [1usize, 4] {
        let params = Params {
            reorder: true,
            keep_order: false,
            rescue: true,
            threads,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(record_set(&out), record_set(&input), "threads={threads}");
    }
}

/// FASTQ of `n` overlapping windows (length `win`) of a fixed pseudo-random
/// reference, emitted in a SHUFFLED order with header text from `name(i)`.
/// The windows share minimizers so clustering re-groups them; because file
/// order is shuffled, clustered order differs from file order — so a
/// positional counter in the name scrambles under clustering (the case where
/// keep-order pays off). Bare `+`, so a keep-order archive round-trips
/// byte-for-byte.
fn windowed_input(name: impl Fn(usize) -> String, n: usize, win: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut x = 0x1234_5678u32;
    let mut lcg = || {
        x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        x
    };
    let mut refseq = Vec::with_capacity(n + win);
    for _ in 0..n + win {
        refseq.push(bases[((lcg() >> 16) & 3) as usize]);
    }
    // Window starts, Fisher-Yates shuffled so file order != clustered order.
    let mut starts: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        starts.swap(i, lcg() as usize % (i + 1));
    }
    let mut v = Vec::new();
    for i in 0..n {
        v.extend_from_slice(name(i).as_bytes());
        v.push(b'\n');
        let s = starts[i];
        v.extend_from_slice(&refseq[s..s + win]);
        v.extend_from_slice(b"\n+\n");
        v.extend(std::iter::repeat_n(b'I', win));
        v.push(b'\n');
    }
    v
}

/// Sorted multiset of sequence lines (record line 1) — the content a
/// reorder-lossy mode must preserve even as it renumbers names.
fn seq_set(fastq: &[u8]) -> Vec<Vec<u8>> {
    let lines: Vec<&[u8]> = fastq.split(|&b| b == b'\n').collect();
    let mut s: Vec<Vec<u8>> = lines
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| c[1].to_vec())
        .collect();
    s.sort();
    s
}

#[test]
fn discard_order_renumbers_and_preserves_content() {
    // Counter-named, reorder-inducing input. Discard-order regenerates the
    // names as a fresh 1..n counter in OUTPUT order (reorder-lossy for names)
    // while preserving the sequence content exactly.
    let input = windowed_input(|i| format!("@read.{} {}", i + 1, i + 1), 3000, 40);
    for threads in [1usize, 4] {
        let params = Params {
            reorder: true,
            regenerate_names: true,
            threads,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        // Discard-order is a non-keep-order layout with regenerated names.
        assert!(!peek(&archive[..]).unwrap().keep_order, "threads={threads}");
        // inspect must skip the name-template frame and report correctly.
        let info = inspect(io::Cursor::new(&archive[..])).unwrap();
        assert!(info.regenerated_names, "threads={threads}");
        assert_eq!(info.reads, 3000, "threads={threads}");

        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();

        let lines: Vec<&[u8]> = out.split(|&b| b == b'\n').collect();
        let recs: Vec<&[&[u8]]> = lines.chunks(4).filter(|c| c.len() == 4).collect();
        assert_eq!(recs.len(), 3000, "threads={threads}");
        // Names regenerated sequentially in output order.
        for (k, c) in recs.iter().enumerate() {
            assert_eq!(
                c[0],
                format!("@read.{} {}", k + 1, k + 1).as_bytes(),
                "name at output {k} (threads={threads})"
            );
        }
        // Sequence multiset preserved exactly.
        assert_eq!(seq_set(&out), seq_set(&input), "threads={threads}");
    }
}

#[test]
fn discard_order_renumbers_non_counter_names() {
    // Illumina tile/x/y names aren't a per-read counter (x/y vary
    // non-monotonically), so they can't be reproduced from position. Renumber
    // mode still engages: the reads are relabeled with a fresh 1..n counter (no
    // name stream, no permutation) rather than falling back to a lossless
    // order-preserving layout. Sequences are preserved as a set; original names
    // are intentionally discarded.
    let input = windowed_input(
        |i| {
            format!(
                "@INST:1:FC:1:1101:{}:{}",
                1000 + (i * 7) % 500,
                2000 + (i * 13) % 500
            )
        },
        2000,
        40,
    );
    let params = Params {
        reorder: true,
        regenerate_names: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    // No permutation stored, names reported as regenerated.
    assert!(!peek(&archive[..]).unwrap().keep_order);
    assert!(
        inspect(io::Cursor::new(&archive[..]))
            .unwrap()
            .regenerated_names
    );

    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    // Sequences preserved exactly as a set (identity discarded).
    assert_eq!(seq_set(&out), seq_set(&input));
    // Names are a fresh 1..n counter in output order.
    let lines: Vec<&[u8]> = out.split(|&b| b == b'\n').collect();
    let recs: Vec<&[&[u8]]> = lines.chunks(4).filter(|c| c.len() == 4).collect();
    assert_eq!(recs.len(), 2000);
    for (k, c) in recs.iter().enumerate() {
        assert_eq!(c[0], format!("@{}", k + 1).as_bytes(), "name at output {k}");
    }
}

#[test]
fn adaptive_keeps_order_for_counter_names() {
    // Counter-style names (the `.N N` pattern) delta-code to almost nothing in
    // original order, so the permutation is cheaper than the scrambled-counter
    // clustered-order stream: adaptive should keep order — and then restore it
    // byte-for-byte.
    let input = windowed_input(|i| format!("@read.{i} {i}"), 2000, 40);
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    assert!(
        peek(&archive[..]).unwrap().keep_order,
        "counter names should trigger keep_order"
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, input, "keep_order must restore original order exactly");
}

#[test]
fn adaptive_drops_order_for_random_names() {
    // Avalanched (splitmix64) names look i.i.d., so they carry no order
    // structure: original- and clustered-order coding cost the same and the
    // permutation is pure overhead — adaptive should NOT keep order.
    let splitmix = |i: usize| -> u64 {
        let mut z = (i as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    };
    let input = windowed_input(|i| format!("@{:016x}", splitmix(i)), 2000, 40);
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    assert!(
        !peek(&archive[..]).unwrap().keep_order,
        "random names should not keep order (permutation is pure overhead)"
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(record_set(&out), record_set(&input));
}

/// Duplicate-rich reads for one mate of a paired set (`n` spots), sharing
/// sequences across spots so clustering has real work — including a
/// reverse-complement so the flip path is exercised.
fn dup_rich_mate(mate: u8, n: u32) -> Vec<u8> {
    let a = b"ACGTTTGACCGATTGCAACGT";
    let ra = fqxv_reorder::revcomp(a);
    let mut v = Vec::new();
    for i in 0..n {
        let s = match (i + mate as u32) % 3 {
            0 => a.to_vec(),
            1 => ra.clone(),
            _ => b"TTTTGGGGCCCCAAAATTTTG".to_vec(),
        };
        v.extend_from_slice(format!("@spot.{i}/{mate}\n").as_bytes());
        v.extend_from_slice(&s);
        v.extend_from_slice(format!("\n+\n{}\n", "I".repeat(s.len())).as_bytes());
    }
    v
}

fn paired_readers<'a>(r1: &'a [u8], r2: &'a [u8]) -> Vec<Box<dyn Read + Send + 'a>> {
    vec![Box::new(r1) as Box<dyn Read + Send>, Box::new(r2)]
}

#[test]
fn reorder_paired_preserves_order_and_splits() {
    let r1 = dup_rich_mate(1, 30);
    let r2 = dup_rich_mate(2, 30);
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress_multi(paired_readers(&r1, &r2), &mut archive, params).unwrap();

    // Grouped reorder records the real group size and stays a reorder archive.
    assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
    assert!(inspect(io::Cursor::new(&archive[..])).unwrap().reordered);

    // Interleaved decode matches a plain (non-reorder) archive byte-for-byte,
    // i.e. the permutation fully restored the original spot interleaving.
    let mut plain = Vec::new();
    compress_multi(paired_readers(&r1, &r2), &mut plain, Params::default()).unwrap();
    let (mut expected, mut got) = (Vec::new(), Vec::new());
    decompress(&plain[..], &mut expected, 1).unwrap();
    decompress(&archive[..], &mut got, 1).unwrap();
    assert_eq!(got, expected, "grouped reorder must restore spot order");

    // Split decode reconstructs each mate file exactly.
    let (mut o1, mut o2) = (Vec::new(), Vec::new());
    {
        let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
        decompress_split(&archive[..], &mut outs, 1).unwrap();
    }
    assert_eq!(o1, r1);
    assert_eq!(o2, r2);
}

#[test]
fn reorder_paired_is_thread_count_deterministic() {
    let r1 = dup_rich_mate(1, 50);
    let r2 = dup_rich_mate(2, 50);
    let mut archives = Vec::new();
    for threads in [1usize, 4] {
        let params = Params {
            reorder: true,
            threads,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress_multi(paired_readers(&r1, &r2), &mut archive, params).unwrap();
        archives.push(archive);
    }
    assert_eq!(
        archives[0], archives[1],
        "reorder output must not vary by threads"
    );
}

#[test]
fn reorder_interleaved_single_stream_splits() {
    // One already-interleaved stream (as `sracha get -Z` emits), reordered.
    let mut stream = Vec::new();
    for i in 0..20u32 {
        for mate in 1..=2u8 {
            let s = if (i + mate as u32).is_multiple_of(2) {
                "ACGTACGTAC"
            } else {
                "TTTTGGGGCC"
            };
            stream.extend_from_slice(format!("@spot.{i}/{mate}\n{s}\n+\nIIIIIIIIII\n").as_bytes());
        }
    }
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress_interleaved(&stream[..], &mut archive, params, 2).unwrap();
    assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, stream, "interleaved reorder must be byte-exact");
}

#[test]
fn empty_input() {
    let archive = compress_bytes(b"", Params::default());
    let mut fastq = Vec::new();
    let stats = decompress(&archive[..], &mut fastq, 1).unwrap();
    assert_eq!(stats.reads, 0);
    assert!(fastq.is_empty());
}

#[test]
fn bad_magic() {
    let err = decompress(&b"not an fqxv file at all"[..], &mut Vec::new(), 1);
    assert!(matches!(err, Err(Error::BadMagic)));
}

// --- v1 footer: index, determinism --------------------------------------

#[test]
fn block_ranges_cuts_on_reads_bytes_and_spots() {
    // 10 reads × 100 bp, one parse chunk.
    let mut fq = Vec::new();
    for i in 0..10 {
        fq.extend_from_slice(
            format!("@r{i}\n{}\n+\n{}\n", "A".repeat(100), "I".repeat(100)).as_bytes(),
        );
    }
    let chunks = vec![parse_chunk(&fq, 0, fq.len()).unwrap()];

    // Byte budget binds first: cut every 3 reads (300 B ≥ 250 B).
    assert_eq!(
        block_ranges(&chunks, 1000, 250, 1),
        vec![(0, 3), (3, 6), (6, 9), (9, 10)]
    );
    // Read budget binds first.
    assert_eq!(
        block_ranges(&chunks, 4, usize::MAX, 1),
        vec![(0, 4), (4, 8), (8, 10)]
    );
    // A tiny byte budget still only cuts on whole spots (g = 2).
    assert_eq!(
        block_ranges(&chunks, 1000, 1, 2),
        vec![(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)]
    );
}

#[test]
fn archive_is_deterministic_across_threads() {
    let input = make_reads("y", 500);
    let mk = |threads| {
        compress_bytes(
            &input,
            Params {
                block_reads: 32,
                threads,
                ..Params::default()
            },
        )
    };
    // Byte-identical (header, blocks, and footer offsets) regardless of pool.
    assert_eq!(mk(1), mk(4));
}

#[test]
fn inspect_falls_back_without_trailer() {
    let archive = compress_bytes(SAMPLE, Params::default());
    // Drop the trailing "FQXF" magic so the footer can't be located — a
    // partial download loses the EOF trailer this way. inspect must fall back
    // to a forward scan and still report the intact blocks rather than error.
    let truncated = &archive[..archive.len() - 4];
    let info = inspect(io::Cursor::new(truncated)).expect("fallback scan");
    assert_eq!(info.reads, 2);
    assert_eq!(info.blocks, 1);
}

#[test]
fn ragged_lengths_roundtrip_multiblock() {
    // Mixed read lengths (10..=310 bp) exercise variable-length framing; a
    // small block target spreads them over several row groups.
    let mut input = Vec::new();
    for i in 0..30usize {
        let len = 10 + (i % 7) * 50;
        let seq: String = "ACGT".chars().cycle().take(len).collect();
        input.extend_from_slice(format!("@read.{i}\n{seq}\n+\n{}\n", "I".repeat(len)).as_bytes());
    }
    let params = Params {
        block_reads: 5,
        ..Params::default()
    };
    let archive = compress_bytes(&input, params);
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, input, "ragged variable-length reads must round-trip");
}

#[test]
fn nanopore_wide_quality_roundtrip_deterministic() {
    // The nanopore profile: long, variable-length reads with N bases and the
    // full Sanger quality range (Phred 0..=93) — which previously tripped the
    // 64-symbol quality cap. A small block target spreads the long reads over
    // several row groups. Output must be byte-exact and identical regardless
    // of thread count.
    let mut input = Vec::new();
    let mut st = 0x9e37_79b9u32;
    for i in 0..12usize {
        let len = 800 + (i % 5) * 900; // 800..=4400 bp
        let mut seq = String::with_capacity(len);
        let mut qual = String::with_capacity(len);
        for _ in 0..len {
            st ^= st << 13;
            st ^= st >> 17;
            st ^= st << 5;
            let r = st % 20;
            seq.push(if r == 0 {
                'N'
            } else {
                b"ACGT"[(r % 4) as usize] as char
            });
            st ^= st << 13;
            st ^= st >> 17;
            st ^= st << 5;
            qual.push((b'!' + (st % 94) as u8) as char); // '!'..='~'
        }
        input.extend_from_slice(format!("@read.{i} ch={i}\n{seq}\n+\n{qual}\n").as_bytes());
    }
    let base = Params {
        block_reads: 4,
        ..Params::default()
    };
    let a1 = compress_bytes(&input, Params { threads: 1, ..base });
    let a4 = compress_bytes(&input, Params { threads: 4, ..base });
    assert_eq!(
        a1, a4,
        "nanopore archive must be deterministic across threads"
    );
    let mut out = Vec::new();
    decompress(&a1[..], &mut out, 1).unwrap();
    assert_eq!(
        out, input,
        "wide-quality long reads must round-trip byte-exact"
    );
}

// --- integrity: CRC detection, recovery, truncation --------------------

/// Build a multi-block archive of `n` uniform reads (small block target).
fn multiblock_archive(n: usize, block_reads: usize) -> Vec<u8> {
    let input = make_reads("x", n);
    compress_bytes(
        &input,
        Params {
            block_reads,
            ..Params::default()
        },
    )
}

#[test]
fn verify_accepts_intact_archive() {
    let archive = multiblock_archive(40, 8);
    verify(io::Cursor::new(&archive), 1).expect("intact archive verifies");
}

#[test]
fn verify_rejects_payload_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    // Header(10) + [8 len][4 crc]; the first payload byte is at offset 22.
    archive[HEADER_LEN + FRAME_HEAD_LEN] ^= 0x01;
    let err = verify(io::Cursor::new(&archive), 1).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn parallel_whole_file_crc_matches_serial() {
    // Exceed CHUNK (1 MiB) and the batch size so several parallel batches run,
    // then confirm the combined result is byte-identical to a single-pass CRC.
    let data: Vec<u8> = (0..5_000_000u32)
        .map(|i| (i.wrapping_mul(2_654_435_761) >> 13) as u8)
        .collect();
    let got = verify_whole_file_crc(&mut io::Cursor::new(&data), data.len() as u64, 0).unwrap();
    assert_eq!(got, crc32c(&data), "full buffer");
    // A partial covered_len must hash only that prefix (not a chunk boundary).
    let partial = 3_000_001usize;
    let got = verify_whole_file_crc(&mut io::Cursor::new(&data), partial as u64, 0).unwrap();
    assert_eq!(got, crc32c(&data[..partial]), "partial prefix");
}

#[test]
fn verify_rejects_footer_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    // Flip a byte inside the footer body (just before the EOF trailer).
    let i = archive.len() - TRAILER_LEN - 1;
    archive[i] ^= 0x01;
    let err = verify(io::Cursor::new(&archive), 1).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn verify_roundtrip_accepts_intact_archive() {
    // Returns the decoded read count so a caller can cross-check what it wrote.
    let archive = multiblock_archive(40, 8);
    let reads = verify_roundtrip(io::Cursor::new(&archive), 1).expect("intact archive round-trips");
    assert_eq!(reads, 40);
}

#[test]
fn verify_roundtrip_honors_thread_budget() {
    // The `threads` argument must reach both the whole-file CRC and the decode —
    // the plumbing that makes a CLI `--verify` respect the command's `--threads`
    // instead of grabbing every core. Output is thread-count-independent (the
    // determinism invariant), so a 1-worker and a many-worker verify must agree.
    let archive = multiblock_archive(40, 8);
    let single = verify_roundtrip(io::Cursor::new(&archive), 1).expect("1-thread verify");
    let many = verify_roundtrip(io::Cursor::new(&archive), 8).expect("many-thread verify");
    assert_eq!(
        single, many,
        "verify read count must not depend on thread budget"
    );
    assert_eq!(single, 40);
}

#[test]
fn verify_roundtrip_counts_empty_archive() {
    // An empty input still writes a valid header/terminator/footer; verify must
    // handle the zero-read case without erroring.
    let archive = compress_bytes(b"", Params::default());
    let reads = verify_roundtrip(io::Cursor::new(&archive), 1).expect("empty archive round-trips");
    assert_eq!(reads, 0);
}

#[test]
fn verify_roundtrip_rejects_payload_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    archive[HEADER_LEN + FRAME_HEAD_LEN] ^= 0x01;
    let err = verify_roundtrip(io::Cursor::new(&archive), 1).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn verify_roundtrip_accepts_reorder_archive() {
    // The globally-clustered layout has no whole-file CRC; verify_roundtrip must
    // decode it (frame CRCs + output digest) and still report the read count.
    let input = dup_rich_input('q');
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            ..Params::default()
        },
    )
    .unwrap();
    assert_eq!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER,
        "test archive must use the reorder layout"
    );
    let reads =
        verify_roundtrip(io::Cursor::new(&archive), 1).expect("intact reorder archive round-trips");
    assert_eq!(reads, 40);
}

#[test]
fn verify_roundtrip_rejects_reorder_corruption() {
    let input = dup_rich_input('q');
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            ..Params::default()
        },
    )
    .unwrap();
    // Flip a content byte in the trailing whole-output digest frame (the last
    // thing in this layout); its frame CRC must reject it. Deliberately not a
    // framing length prefix, which is a separate (allocation-guard) path.
    *archive.last_mut().unwrap() ^= 0xFF;
    assert!(
        verify_roundtrip(io::Cursor::new(&archive), 1).is_err(),
        "corrupt reorder archive must not round-trip"
    );
}

#[test]
fn verify_roundtrip_rejects_corrupt_reorder_read_count() {
    // Regression: a corrupt read-count prefix (the u64 immediately past the
    // header) must fail *gracefully* rather than abort the process on a
    // multi-terabyte allocation. The reorder decode's capacity hints used to be
    // infallible `Vec::with_capacity(n)`.
    let input = dup_rich_input('q');
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            ..Params::default()
        },
    )
    .unwrap();
    assert_eq!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER
    );
    // The read count is a little-endian u64 at `HEADER_LEN`; set its top byte to
    // blow the value past any real dataset.
    archive[HEADER_LEN + 7] ^= 0xFF;
    let err = verify_roundtrip(io::Cursor::new(&archive), 1).unwrap_err();
    assert!(matches!(err, Error::Malformed(_)), "got {err:?}");
}

/// Write `bytes` to a fresh temp file and return an open handle plus its path.
/// The name is unique per process *and* per call so it is safe whether tests
/// run as separate processes (nextest) or as threads (`cargo test`).
fn temp_archive(bytes: &[u8]) -> (File, std::path::PathBuf) {
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQ: AtomicU64 = AtomicU64::new(0);
    let n = SEQ.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("fqxv-quick-{}-{n}.fqxv", std::process::id()));
    std::fs::write(&path, bytes).unwrap();
    (File::open(&path).unwrap(), path)
}

#[test]
fn verify_quick_accepts_intact_archive() {
    let archive = multiblock_archive(40, 8);
    let (file, path) = temp_archive(&archive);
    let result = verify_quick(&file, 1);
    std::fs::remove_file(&path).ok();
    result.expect("intact archive passes quick verify");
}

#[test]
fn verify_quick_rejects_payload_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    archive[HEADER_LEN + FRAME_HEAD_LEN] ^= 0x01;
    let (file, path) = temp_archive(&archive);
    let err = verify_quick(&file, 1).unwrap_err();
    std::fs::remove_file(&path).ok();
    // The per-block check localizes the failure to the offending block.
    assert!(
        matches!(&err, Error::Corrupt { what } if what.starts_with("block")),
        "got {err:?}"
    );
}

#[test]
fn verify_quick_falls_back_for_reorder_layout() {
    // The globally-clustered layout has no per-block footer index, so quick
    // verify must transparently run the full decode-driven check.
    let input = dup_rich_input('q');
    let params = Params {
        reorder: true,
        ..Params::default()
    };
    let mut archive = Vec::new();
    compress(&input[..], &mut archive, params).unwrap();
    // Header flags byte sits at offset 8 ([4]magic [2]ver [1]order [1]binning).
    assert_eq!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER,
        "test archive must use the reorder layout to exercise the fallback"
    );
    let (file, path) = temp_archive(&archive);
    let result = verify_quick(&file, 1);
    std::fs::remove_file(&path).ok();
    result.expect("intact reorder archive passes quick verify via fallback");
}

#[test]
fn verify_report_intact_lists_passing_checks() {
    let archive = multiblock_archive(40, 8);
    let (file, path) = temp_archive(&archive);
    let report = verify_report(&file, false, 1).expect("readable archive");
    std::fs::remove_file(&path).ok();
    assert!(report.passed());
    assert!(report.failed_blocks.is_empty());
    let names: Vec<&str> = report.checks.iter().map(|c| c.name.as_str()).collect();
    assert_eq!(names, ["header", "footer", "block CRCs", "whole-file CRC"]);
    assert!(report.blocks_total >= 2, "expected several small blocks");
}

#[test]
fn verify_report_localizes_corrupt_block() {
    let archive = multiblock_archive(40, 8);
    let footer = read_footer(&mut io::Cursor::new(&archive)).unwrap();
    assert!(footer.groups.len() >= 3, "need multiple blocks to localize");
    // Corrupt the payload of the second block (index 1): its footer offset
    // points at the frame head ([4 marker][8 len][4 crc]), so the payload follows.
    let mut archive = archive;
    archive[footer.groups[1].0 as usize + FRAME_HEAD_LEN] ^= 0xFF;
    let (file, path) = temp_archive(&archive);
    let report = verify_report(&file, false, 1).expect("still structurally readable");
    std::fs::remove_file(&path).ok();

    assert!(!report.passed());
    assert_eq!(report.failed_blocks, vec![1]);
    let blocks = report
        .checks
        .iter()
        .find(|c| c.name == "block CRCs")
        .unwrap();
    assert!(!blocks.ok);
    assert!(
        blocks.detail.contains("failed: 1"),
        "detail: {}",
        blocks.detail
    );
}

#[test]
fn verify_report_quick_skips_whole_file_crc() {
    let archive = multiblock_archive(40, 8);
    let (file, path) = temp_archive(&archive);
    let report = verify_report(&file, true, 1).expect("readable archive");
    std::fs::remove_file(&path).ok();
    assert!(report.passed());
    // Quick mode stops at the per-block CRCs; no whole-file digest row.
    let names: Vec<&str> = report.checks.iter().map(|c| c.name.as_str()).collect();
    assert_eq!(names, ["header", "footer", "block CRCs"]);
}

#[test]
fn decompress_detects_block_corruption() {
    let mut archive = multiblock_archive(40, 8);
    archive[HEADER_LEN + FRAME_HEAD_LEN] ^= 0xFF;
    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn stream_digests_localize_and_pin_boundaries() {
    let d = |names: &[&[u8]], lens: &[u32], seq: &[u8], qual: &[u8]| {
        stream_digests(names.len(), names.iter().copied(), lens, seq, qual)
    };
    let base = d(&[b"r1", b"r2"], &[3, 3], b"ACGTTT", b"IIIFFF");

    // Localization: a change in one stream perturbs ONLY that stream's digest,
    // so a mismatch names the stream that regressed.
    let name_changed = d(&[b"r1", b"rX"], &[3, 3], b"ACGTTT", b"IIIFFF");
    assert_ne!(
        base.names, name_changed.names,
        "name change -> names digest"
    );
    assert_eq!(base.seq, name_changed.seq, "name change leaves seq digest");
    assert_eq!(
        base.qual, name_changed.qual,
        "name change leaves qual digest"
    );

    let seq_changed = d(&[b"r1", b"r2"], &[3, 3], b"ACGTTA", b"IIIFFF");
    assert_ne!(base.seq, seq_changed.seq, "seq change -> seq digest");
    assert_eq!(
        base.names, seq_changed.names,
        "seq change leaves names digest"
    );
    assert_eq!(base.qual, seq_changed.qual, "seq change leaves qual digest");

    let qual_changed = d(&[b"r1", b"r2"], &[3, 3], b"ACGTTT", b"IIIFF#");
    assert_ne!(base.qual, qual_changed.qual, "qual change -> qual digest");
    assert_eq!(
        base.names, qual_changed.names,
        "qual change leaves names digest"
    );
    assert_eq!(base.seq, qual_changed.seq, "qual change leaves seq digest");

    // Boundary pinning per stream: the same concatenated bytes split differently
    // between reads must not collide (a byte "sliding" across a read boundary is
    // exactly the silent-corruption shape the length folds catch).
    let a = d(&[b"AB", b"C"], &[1, 2], b"G", b"I");
    let b = d(&[b"A", b"BC"], &[2, 1], b"G", b"I");
    assert_ne!(a.names, b.names, "names boundary");
    let c = d(&[b"n"], &[1, 2], b"ACG", b"III");
    let e = d(&[b"n"], &[2, 1], b"ACG", b"III");
    assert_ne!(c.seq, e.seq, "seq boundary");
    assert_ne!(c.qual, e.qual, "qual boundary");
}

#[test]
fn decompress_detects_and_localizes_content_digest_mismatch() {
    // The failure mode CRC cannot see: frame CRC intact, but the decoded content
    // does not match a stored digest (a codec round-trip bug). Simulate by
    // flipping a byte inside one of the three per-stream digests and repairing the
    // frame CRC, so only the content-digest check can reject it — and the reported
    // stream must match the digest that was corrupted.
    for (digest_idx, stream) in [(0usize, "names"), (1, "sequence"), (2, "quality")] {
        let mut archive = multiblock_archive(20, 64); // block_reads > n => one block
        let payload_start = HEADER_LEN + FRAME_HEAD_LEN;
        archive[payload_start + digest_idx * DIGEST_LEN] ^= 0xFF; // that stream's digest
        let len_at = HEADER_LEN + BLOCK_MAGIC.len();
        let len = u64::from_le_bytes(archive[len_at..len_at + 8].try_into().unwrap()) as usize;
        let repaired = crc32c(&archive[payload_start..payload_start + len]);
        archive[len_at + 8..payload_start].copy_from_slice(&repaired.to_le_bytes());

        let mut out = Vec::new();
        let err = decompress(&archive[..], &mut out, 1).unwrap_err();
        let expected = format!("block {stream} digest");
        assert!(
            matches!(&err, Error::Corrupt { what } if *what == expected),
            "flipping the {stream} digest should report {expected:?}, got {err:?}"
        );
    }
}

#[test]
fn content_digest_accepts_lossy_binning_roundtrip() {
    // The encode-side digest is over the POST-binning quality, so a lossy
    // archive must decode without a false digest failure (guards the scoping:
    // the digest checks the round-trip, not the lossy transform).
    let input = make_reads("x", 30);
    let archive = compress_bytes(
        &input,
        Params {
            quality_binning: QualityBinning::Bin4,
            ..Params::default()
        },
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1)
        .expect("lossy archive round-trips without a false content-digest failure");
}

#[test]
fn decompress_rejects_header_bit_flip() {
    // The header CRC catches a flipped field byte (here the lossy binning tag)
    // that would otherwise silently change how the archive is interpreted.
    let mut archive = multiblock_archive(20, 64);
    archive[HDR_OFF_BINNING] ^= 0x02; // quality-binning tag, inside the CRC'd header fields
    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(
        matches!(&err, Error::Corrupt { what } if what == "header"),
        "got {err:?}"
    );
}

#[test]
fn reorder_header_is_crc_protected() {
    // The reorder layout previously left its header (incl. the lossy binning
    // tag and flags) covered by no checksum; the header CRC now covers it too.
    let input = dup_rich_input('q');
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            ..Params::default()
        },
    )
    .unwrap();
    assert_eq!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER
    );
    archive[HDR_OFF_BINNING] ^= 0x02; // binning tag in the reorder layout's header
    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(
        matches!(&err, Error::Corrupt { what } if what == "header"),
        "got {err:?}"
    );
}

#[test]
fn decompress_detects_reorder_output_digest_mismatch() {
    // Reorder analog of the per-block content-digest test: repair the trailing
    // digest frame's CRC after corrupting the stored digest, so only the
    // whole-output content check can reject the (otherwise valid) archive.
    let input = dup_rich_input('q');
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            ..Params::default()
        },
    )
    .unwrap();
    assert_eq!(
        archive[HDR_OFF_FLAGS] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER
    );
    // Trailing frame is [4 len=8][4 crc][8 digest] at the very end (no footer).
    let len = archive.len();
    let dig_start = len - DIGEST_LEN;
    let crc_start = dig_start - CRC_LEN;
    archive[len - 1] ^= 0xFF; // flip a stored-digest byte
    let repaired = crc32c(&archive[dig_start..len]);
    archive[crc_start..dig_start].copy_from_slice(&repaired.to_le_bytes());
    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(
        matches!(&err, Error::Corrupt { what } if what.contains("output digest")),
        "got {err:?}"
    );
}

#[test]
fn oversized_block_length_is_rejected_not_allocated() {
    let mut archive = multiblock_archive(40, 8);
    // Overwrite the first block's [8] length (past the marker) with a hostile
    // value; the reader must reject it up front instead of allocating exabytes.
    let len_at = HEADER_LEN + BLOCK_MAGIC.len();
    archive[len_at..len_at + 8].copy_from_slice(&u64::MAX.to_le_bytes());
    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(matches!(err, Error::Malformed(_)), "got {err:?}");
}

#[test]
fn recover_skips_corrupt_block_and_keeps_the_rest() {
    let mut archive = multiblock_archive(40, 8);
    let footer = read_footer(&mut io::Cursor::new(&archive)).unwrap();
    assert!(
        footer.groups.len() >= 3,
        "need several blocks for this test"
    );
    let (off1, rc1) = footer.groups[1];
    // Corrupt one byte in block 1's payload (past its frame head).
    archive[off1 as usize + FRAME_HEAD_LEN] ^= 0xFF;

    let mut out = Vec::new();
    let rec = decompress_recover(io::Cursor::new(&archive), &mut out, 1).unwrap();
    assert_eq!(rec.blocks_skipped, 1);
    assert_eq!(rec.reads_lost, u64::from(rc1));
    assert_eq!(rec.stats.reads, footer.total_reads - u64::from(rc1));
    assert_eq!(
        rec.blocks_recovered,
        footer.groups.len() as u64 - 1,
        "every other block recovered"
    );
    // Output is valid FASTQ for the recovered reads (4 lines each).
    assert_eq!(
        out.iter().filter(|&&b| b == b'\n').count() as u64,
        rec.stats.reads * 4
    );
}

/// Footer-independent recovery: a lost tail (no footer) still recovers every
/// whole block via the sync-marker scan.
#[test]
fn recover_scans_when_footer_is_lost() {
    let full = multiblock_archive(40, 8);
    let footer = read_footer(&mut io::Cursor::new(&full)).unwrap();
    assert!(footer.groups.len() >= 3, "need several blocks");
    // Drop the footer + trailer at the terminator boundary — a truncated tail.
    let n = full.len();
    let footer_offset =
        u64::from_le_bytes(full[n - TRAILER_LEN..n - 4].try_into().unwrap()) as usize;
    let truncated = full[..footer_offset - (BLOCK_MAGIC.len() + 8)].to_vec();

    let mut out = Vec::new();
    let rec = decompress_recover(io::Cursor::new(&truncated), &mut out, 1).unwrap();
    assert_eq!(rec.blocks_recovered, footer.groups.len() as u64);
    assert_eq!(rec.stats.reads, footer.total_reads);
    let mut full_out = Vec::new();
    decompress(&full[..], &mut full_out, 1).unwrap();
    assert_eq!(out, full_out, "scan recovery yields the same reads");
}

/// Footer lost AND a middle block corrupt: the scan resynchronizes to the next
/// marker, skipping the bad block and recovering the rest — where streaming
/// `decompress` would stop at the bad block and lose everything after it.
#[test]
fn recover_scans_past_a_corrupt_block_without_footer() {
    let full = multiblock_archive(40, 8);
    let footer = read_footer(&mut io::Cursor::new(&full)).unwrap();
    assert!(footer.groups.len() >= 3, "need several blocks");
    let n = full.len();
    let footer_offset =
        u64::from_le_bytes(full[n - TRAILER_LEN..n - 4].try_into().unwrap()) as usize;
    let mut archive = full[..footer_offset - (BLOCK_MAGIC.len() + 8)].to_vec();
    // Corrupt a payload byte of the middle block (its marker stays intact).
    let (off1, rc1) = footer.groups[1];
    archive[off1 as usize + FRAME_HEAD_LEN] ^= 0xFF;

    let mut out = Vec::new();
    let rec = decompress_recover(io::Cursor::new(&archive), &mut out, 1).unwrap();
    assert_eq!(
        rec.blocks_recovered,
        footer.groups.len() as u64 - 1,
        "every block but the corrupt one recovered"
    );
    assert_eq!(rec.stats.reads, footer.total_reads - u64::from(rc1));
    assert_eq!(
        out.iter().filter(|&&b| b == b'\n').count() as u64,
        rec.stats.reads * 4
    );
}

#[test]
fn truncated_at_block_boundary_streams_prefix() {
    let full = multiblock_archive(40, 8);
    // The trailer's back-pointer gives footer_offset; the terminator frame
    // (marker + [8] 0) sits just before it, so subtracting its size lands on a
    // clean block boundary.
    let n = full.len();
    let footer_offset =
        u64::from_le_bytes(full[n - TRAILER_LEN..n - 4].try_into().unwrap()) as usize;
    let truncated = &full[..footer_offset - (BLOCK_MAGIC.len() + 8)];

    // Streaming decode reads every whole block, then stops at the clean EOF.
    let mut out_trunc = Vec::new();
    decompress(truncated, &mut out_trunc, 1).expect("prefix decodes");
    let mut out_full = Vec::new();
    decompress(&full[..], &mut out_full, 1).unwrap();
    assert_eq!(
        out_trunc, out_full,
        "boundary-truncated file yields all blocks"
    );
}

#[test]
fn reorder_archive_detects_frame_corruption() {
    let input = make_reads("y", 200);
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            keep_order: true,
            rescue: false,
            block_reads: 64,
            ..Params::default()
        },
    )
    .unwrap();
    // Flip a byte well past the header, inside a framed payload. The frame
    // CRC (or a downstream consistency check) must catch it — never a silent
    // wrong decode.
    let mid = archive.len() / 2;
    archive[mid] ^= 0xFF;
    assert!(verify(io::Cursor::new(&archive), 1).is_err());
}

// --- random access: per-stream column projection (Gap 1 + Gap 3) -------------

/// One record's `(header, sequence, quality)` bytes.
type Record = (Vec<u8>, Vec<u8>, Vec<u8>);
/// A per-read column of byte slices (names, sequences, or qualities).
type Column = Vec<Vec<u8>>;

/// FASTQ with per-read varied names, sequence lengths, and quality, plus the
/// matching [`Record`] triples for checking a projection reconstructs exactly
/// what was stored.
fn varied_records(n: usize) -> (Vec<u8>, Vec<Record>) {
    let bases = *b"ACGT";
    let mut bytes = Vec::new();
    let mut recs = Vec::with_capacity(n);
    for i in 0..n {
        let len = 20 + (i % 30);
        let seq: Vec<u8> = (0..len).map(|j| bases[(i + j) % 4]).collect();
        let qual: Vec<u8> = (0..len).map(|j| b'!' + ((i + j) % 40) as u8).collect();
        let header = format!("read.{i} sample:{i}").into_bytes();
        bytes.push(b'@');
        bytes.extend_from_slice(&header);
        bytes.push(b'\n');
        bytes.extend_from_slice(&seq);
        bytes.extend_from_slice(b"\n+\n");
        bytes.extend_from_slice(&qual);
        bytes.push(b'\n');
        recs.push((header, seq, qual));
    }
    (bytes, recs)
}

/// Slice a group's concatenated stream buffer back into per-read pieces.
fn split_by_lens(buf: &[u8], lens: &[u32]) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(lens.len());
    let mut off = 0usize;
    for &l in lens {
        out.push(buf[off..off + l as usize].to_vec());
        off += l as usize;
    }
    out
}

/// The footer's per-stream index lets each column be fetched and decoded in
/// isolation, and every projected stream reconstructs exactly what was stored.
#[test]
fn index_projects_every_stream_and_roundtrips() {
    let (input, recs) = varied_records(130);
    // Small groups so the index has several rows to project across.
    let archive = compress_bytes(
        &input,
        Params {
            block_reads: 16,
            ..Params::default()
        },
    );
    let index = Index::read(io::Cursor::new(&archive)).unwrap();
    assert_eq!(index.total_reads(), recs.len() as u64);
    assert!(index.groups().len() >= 4, "expected several row groups");

    let (mut names, mut seqs, mut quals): (Column, Column, Column) =
        (Vec::new(), Vec::new(), Vec::new());
    for g in 0..index.groups().len() {
        // Names.
        let r = &index.byte_ranges(&[g], Stream::Names).unwrap()[0];
        let coded = &archive[r.start as usize..r.end as usize];
        index.verify_stream(g, Stream::Names, coded).unwrap();
        names.extend(decode_names(coded).unwrap());
        // Sequence.
        let r = &index.byte_ranges(&[g], Stream::Sequence).unwrap()[0];
        let coded = &archive[r.start as usize..r.end as usize];
        index.verify_stream(g, Stream::Sequence, coded).unwrap();
        let (lens, seq) = decode_sequence(coded).unwrap();
        seqs.extend(split_by_lens(&seq, &lens));
        // Quality.
        let r = &index.byte_ranges(&[g], Stream::Quality).unwrap()[0];
        let coded = &archive[r.start as usize..r.end as usize];
        index.verify_stream(g, Stream::Quality, coded).unwrap();
        let (lens, qual) = decode_quality(coded).unwrap();
        quals.extend(split_by_lens(&qual, &lens));
    }

    let exp_names: Vec<Vec<u8>> = recs.iter().map(|(h, _, _)| h.clone()).collect();
    let exp_seqs: Vec<Vec<u8>> = recs.iter().map(|(_, s, _)| s.clone()).collect();
    let exp_quals: Vec<Vec<u8>> = recs.iter().map(|(_, _, q)| q.clone()).collect();
    assert_eq!(names, exp_names, "projected names");
    assert_eq!(seqs, exp_seqs, "projected sequences");
    assert_eq!(quals, exp_quals, "projected qualities");
}

/// A whole fetched block payload decodes via the public IO-free entry point.
#[test]
fn index_decodes_a_whole_fetched_block() {
    let (input, recs) = varied_records(40);
    let archive = compress_bytes(
        &input,
        Params {
            block_reads: 16,
            ..Params::default()
        },
    );
    let index = Index::read(io::Cursor::new(&archive)).unwrap();
    let g = &index.groups()[0];
    // Fetch the block frame at its offset, check the frame CRC, decode the payload.
    let payload = read_block(&mut io::Cursor::new(&archive[g.block_offset as usize..]), 0)
        .unwrap()
        .unwrap();
    let block = decode_block_contents(&payload).unwrap();
    assert_eq!(block.names.len(), g.read_count as usize);
    assert_eq!(block.names[0], recs[0].0, "first read name");
    let first_len = block.lengths[0] as usize;
    assert_eq!(
        &block.sequence[..first_len],
        &recs[0].1[..],
        "first read seq"
    );
}

/// #142 F4: a block frame claiming the maximum payload but truncated must error
/// as `Truncated` after reading only the bytes present — not zero-fill the full
/// 2 GB claim before the short read is discovered.
#[test]
fn read_block_truncated_body_does_not_alloc_the_claim() {
    let mut frame = Vec::new();
    frame.extend_from_slice(&BLOCK_MAGIC);
    frame.extend_from_slice(&MAX_BLOCK_PAYLOAD.to_le_bytes()); // len = the cap
    frame.extend_from_slice(&0u32.to_le_bytes()); // crc
    frame.extend_from_slice(b"short"); // body far below the claim
    assert!(matches!(
        read_block(&mut io::Cursor::new(&frame), 0),
        Err(Error::Truncated)
    ));
}

/// A frame length past `MAX_BLOCK_PAYLOAD` is rejected outright.
#[test]
fn read_block_rejects_over_cap_length() {
    let mut frame = Vec::new();
    frame.extend_from_slice(&BLOCK_MAGIC);
    frame.extend_from_slice(&(MAX_BLOCK_PAYLOAD + 1).to_le_bytes());
    frame.extend_from_slice(&0u32.to_le_bytes());
    assert!(matches!(
        read_block(&mut io::Cursor::new(&frame), 0),
        Err(Error::Malformed(_))
    ));
}

/// A bit flipped in a projected stream is caught by the index's per-stream CRC —
/// the guarantee the joint block content digest can't provide on a column fetch.
#[test]
fn index_stream_crc_detects_corruption() {
    let (input, _) = varied_records(40);
    let archive = compress_bytes(
        &input,
        Params {
            block_reads: 16,
            ..Params::default()
        },
    );
    let index = Index::read(io::Cursor::new(&archive)).unwrap();
    let r = &index.byte_ranges(&[0], Stream::Sequence).unwrap()[0];
    let mut coded = archive[r.start as usize..r.end as usize].to_vec();
    coded[0] ^= 0xFF;
    let err = index
        .verify_stream(0, Stream::Sequence, &coded)
        .unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

/// `from_suffix` parses the index out of a fetched tail, and asks for a longer
/// tail (with the exact length) when the first fetch was too short.
#[test]
fn index_from_suffix_handles_short_and_full_tails() {
    let (input, recs) = varied_records(200);
    let archive = compress_bytes(
        &input,
        Params {
            block_reads: 16,
            ..Params::default()
        },
    );
    let file_len = archive.len() as u64;
    let full = Index::read(io::Cursor::new(&archive)).unwrap();

    // The whole file as a suffix parses to the same index.
    match Index::from_suffix(&archive, file_len).unwrap() {
        SuffixParse::Parsed(idx) => {
            assert_eq!(idx.total_reads(), recs.len() as u64);
            assert_eq!(idx.groups().len(), full.groups().len());
            assert_eq!(idx.whole_file_crc(), full.whole_file_crc());
        }
        SuffixParse::NeedAtLeast(_) => panic!("full file should parse"),
    }

    // A tiny tail can't reach the footer; it reports the exact length to refetch,
    // and that refetch then parses.
    let tiny = 16usize.min(archive.len());
    let need = match Index::from_suffix(&archive[archive.len() - tiny..], file_len).unwrap() {
        SuffixParse::NeedAtLeast(n) => n,
        SuffixParse::Parsed(_) => panic!("tiny tail should be insufficient"),
    };
    assert!(need as usize > tiny && need <= file_len);
    let refetched = &archive[archive.len() - need as usize..];
    match Index::from_suffix(refetched, file_len).unwrap() {
        SuffixParse::Parsed(idx) => {
            assert_eq!(idx.groups().len(), full.groups().len());
        }
        SuffixParse::NeedAtLeast(_) => panic!("exact refetch should parse"),
    }
}

/// The per-stream sizes `inspect` reports equal the sum of the index's per-stream
/// lengths — the same footer data, read two ways.
#[test]
fn index_stream_sizes_match_inspect() {
    let (input, _) = varied_records(90);
    let archive = compress_bytes(
        &input,
        Params {
            block_reads: 16,
            ..Params::default()
        },
    );
    let index = Index::read(io::Cursor::new(&archive)).unwrap();
    let sum = |s: Stream| -> u64 {
        index
            .groups()
            .iter()
            .map(|g| {
                let r = g.stream_range(s);
                r.end - r.start
            })
            .sum()
    };
    let info = inspect(io::Cursor::new(&archive)).unwrap();
    assert_eq!(info.names_bytes, sum(Stream::Names));
    assert_eq!(info.seq_bytes, sum(Stream::Sequence));
    assert_eq!(info.qual_bytes, sum(Stream::Quality));
}

/// The footer-less reorder layout has no row-group index, so `Index::read`
/// rejects it rather than pretending random access is possible.
#[test]
fn index_rejects_reorder_layout() {
    let input = make_reads("z", 200);
    let mut archive = Vec::new();
    compress(
        &input[..],
        &mut archive,
        Params {
            reorder: true,
            keep_order: true,
            rescue: false,
            block_reads: 64,
            ..Params::default()
        },
    )
    .unwrap();
    let err = Index::read(io::Cursor::new(&archive)).unwrap_err();
    assert!(matches!(err, Error::Malformed(_)), "got {err:?}");
}
