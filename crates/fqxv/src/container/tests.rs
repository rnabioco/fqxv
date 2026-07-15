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
    assert_eq!(info.format_version, FORMAT_VERSION);
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
        assert_eq!(archive[8] & FLAG_GLOBAL_REORDER, FLAG_GLOBAL_REORDER);
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
        archive[8] & FLAG_GLOBAL_REORDER,
        0,
        "reorder must be skipped (flag clear) for long-read data"
    );
    let mut out = Vec::new();
    decompress(&archive[..], &mut out, 1).unwrap();
    assert_eq!(out, input, "long-read fallback must be byte-exact");
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
    verify(io::Cursor::new(&archive)).expect("intact archive verifies");
}

#[test]
fn verify_rejects_payload_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    // Header(10) + [8 len][4 crc]; the first payload byte is at offset 22.
    archive[HEADER_LEN + 8 + CRC_LEN] ^= 0x01;
    let err = verify(io::Cursor::new(&archive)).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn parallel_whole_file_crc_matches_serial() {
    // Exceed CHUNK (1 MiB) and the batch size so several parallel batches run,
    // then confirm the combined result is byte-identical to a single-pass CRC.
    let data: Vec<u8> = (0..5_000_000u32)
        .map(|i| (i.wrapping_mul(2_654_435_761) >> 13) as u8)
        .collect();
    let got = verify_whole_file_crc(&mut io::Cursor::new(&data), data.len() as u64).unwrap();
    assert_eq!(got, crc32c(&data), "full buffer");
    // A partial covered_len must hash only that prefix (not a chunk boundary).
    let partial = 3_000_001usize;
    let got = verify_whole_file_crc(&mut io::Cursor::new(&data), partial as u64).unwrap();
    assert_eq!(got, crc32c(&data[..partial]), "partial prefix");
}

#[test]
fn verify_rejects_footer_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    // Flip a byte inside the footer body (just before the EOF trailer).
    let i = archive.len() - TRAILER_LEN - 1;
    archive[i] ^= 0x01;
    let err = verify(io::Cursor::new(&archive)).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn verify_roundtrip_accepts_intact_archive() {
    // Returns the decoded read count so a caller can cross-check what it wrote.
    let archive = multiblock_archive(40, 8);
    let reads = verify_roundtrip(io::Cursor::new(&archive)).expect("intact archive round-trips");
    assert_eq!(reads, 40);
}

#[test]
fn verify_roundtrip_counts_empty_archive() {
    // An empty input still writes a valid header/terminator/footer; verify must
    // handle the zero-read case without erroring.
    let archive = compress_bytes(b"", Params::default());
    let reads = verify_roundtrip(io::Cursor::new(&archive)).expect("empty archive round-trips");
    assert_eq!(reads, 0);
}

#[test]
fn verify_roundtrip_rejects_payload_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    archive[HEADER_LEN + 8 + CRC_LEN] ^= 0x01;
    let err = verify_roundtrip(io::Cursor::new(&archive)).unwrap_err();
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
        archive[8] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER,
        "test archive must use the reorder layout"
    );
    let reads =
        verify_roundtrip(io::Cursor::new(&archive)).expect("intact reorder archive round-trips");
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
        verify_roundtrip(io::Cursor::new(&archive)).is_err(),
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
    assert_eq!(archive[8] & FLAG_GLOBAL_REORDER, FLAG_GLOBAL_REORDER);
    // The read count is a little-endian u64 at `HEADER_LEN`; set its top byte to
    // blow the value past any real dataset.
    archive[HEADER_LEN + 7] ^= 0xFF;
    let err = verify_roundtrip(io::Cursor::new(&archive)).unwrap_err();
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
    let result = verify_quick(&file);
    std::fs::remove_file(&path).ok();
    result.expect("intact archive passes quick verify");
}

#[test]
fn verify_quick_rejects_payload_bit_flip() {
    let mut archive = multiblock_archive(40, 8);
    archive[HEADER_LEN + 8 + CRC_LEN] ^= 0x01;
    let (file, path) = temp_archive(&archive);
    let err = verify_quick(&file).unwrap_err();
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
        archive[8] & FLAG_GLOBAL_REORDER,
        FLAG_GLOBAL_REORDER,
        "test archive must use the reorder layout to exercise the fallback"
    );
    let (file, path) = temp_archive(&archive);
    let result = verify_quick(&file);
    std::fs::remove_file(&path).ok();
    result.expect("intact reorder archive passes quick verify via fallback");
}

#[test]
fn verify_report_intact_lists_passing_checks() {
    let archive = multiblock_archive(40, 8);
    let (file, path) = temp_archive(&archive);
    let report = verify_report(&file, false).expect("readable archive");
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
    // points at the [8 len][4 crc] frame header, so the payload follows.
    let mut archive = archive;
    archive[footer.groups[1].0 as usize + 8 + CRC_LEN] ^= 0xFF;
    let (file, path) = temp_archive(&archive);
    let report = verify_report(&file, false).expect("still structurally readable");
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
    let report = verify_report(&file, true).expect("readable archive");
    std::fs::remove_file(&path).ok();
    assert!(report.passed());
    // Quick mode stops at the per-block CRCs; no whole-file digest row.
    let names: Vec<&str> = report.checks.iter().map(|c| c.name.as_str()).collect();
    assert_eq!(names, ["header", "footer", "block CRCs"]);
}

#[test]
fn decompress_detects_block_corruption() {
    let mut archive = multiblock_archive(40, 8);
    archive[HEADER_LEN + 8 + CRC_LEN] ^= 0xFF;
    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
}

#[test]
fn content_digest_distinguishes_streams_and_boundaries() {
    let d = |names: &[&[u8]], lens: &[u32], seq: &[u8], qual: &[u8]| {
        content_digest(names.len(), names.iter().copied(), lens, seq, qual)
    };
    let base = d(&[b"r1", b"r2"], &[3, 3], b"ACGTTT", b"IIIFFF");
    // Sensitive to each of the three decoded streams.
    assert_ne!(
        base,
        d(&[b"r1", b"rX"], &[3, 3], b"ACGTTT", b"IIIFFF"),
        "name"
    );
    assert_ne!(
        base,
        d(&[b"r1", b"r2"], &[3, 3], b"ACGTTA", b"IIIFFF"),
        "seq"
    );
    assert_ne!(
        base,
        d(&[b"r1", b"r2"], &[3, 3], b"ACGTTT", b"IIIFF#"),
        "qual"
    );
    // Boundary pinning: the same concatenated bytes split differently between
    // name and sequence must not collide (a byte "sliding" across a stream
    // boundary is exactly the silent-corruption shape the length folds catch).
    assert_ne!(
        d(&[b"AB"], &[2], b"CD", b"II"),
        d(&[b"ABC"], &[1], b"D", b"I"),
        "boundary"
    );
}

#[test]
fn decompress_detects_content_digest_mismatch() {
    // The failure mode CRC cannot see: frame CRC intact, but the decoded
    // content does not match the stored digest (a codec round-trip bug).
    // Simulate by flipping a byte in the payload's leading digest and repairing
    // the frame CRC, so only the content-digest check can reject it.
    let mut archive = multiblock_archive(20, 64); // block_reads > n => one block
    let payload_start = HEADER_LEN + 8 + CRC_LEN;
    archive[payload_start] ^= 0xFF; // first byte of the content digest
    let len = u64::from_le_bytes(archive[HEADER_LEN..HEADER_LEN + 8].try_into().unwrap()) as usize;
    let repaired = crc32c(&archive[payload_start..payload_start + len]);
    archive[HEADER_LEN + 8..payload_start].copy_from_slice(&repaired.to_le_bytes());

    let mut out = Vec::new();
    let err = decompress(&archive[..], &mut out, 1).unwrap_err();
    assert!(
        matches!(&err, Error::Corrupt { what } if what.contains("content digest")),
        "got {err:?}"
    );
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
    archive[7] ^= 0x02; // quality-binning tag, inside the CRC'd header fields
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
    assert_eq!(archive[8] & FLAG_GLOBAL_REORDER, FLAG_GLOBAL_REORDER);
    archive[7] ^= 0x02; // binning tag in the reorder layout's header
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
    assert_eq!(archive[8] & FLAG_GLOBAL_REORDER, FLAG_GLOBAL_REORDER);
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
    // Overwrite the first block's [8] length with a hostile value; the reader
    // must reject it up front instead of trying to allocate exabytes.
    archive[HEADER_LEN..HEADER_LEN + 8].copy_from_slice(&u64::MAX.to_le_bytes());
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
    // Corrupt one byte in block 1's payload (past its [8 len][4 crc]).
    archive[off1 as usize + 8 + CRC_LEN] ^= 0xFF;

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

#[test]
fn truncated_at_block_boundary_streams_prefix() {
    let full = multiblock_archive(40, 8);
    // The trailer's back-pointer gives footer_offset; the 8-byte terminator
    // sits just before it, so footer_offset - 8 is a clean block boundary.
    let n = full.len();
    let footer_offset =
        u64::from_le_bytes(full[n - TRAILER_LEN..n - 4].try_into().unwrap()) as usize;
    let truncated = &full[..footer_offset - 8];

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
    assert!(verify(io::Cursor::new(&archive)).is_err());
}
