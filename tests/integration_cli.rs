use std::{fs, process::Command};

use rust_htslib::bam::{Format, Header, HeaderView, Writer, header::HeaderRecord, record::Record};
use tempfile::tempdir;

fn write_fixture_bam(path: &std::path::Path) {
    let mut header = Header::new();
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &"chr1")
            .push_tag(b"LN", &100_u32),
    );
    header.push_record(
        HeaderRecord::new(b"RG")
            .push_tag(b"ID", &"rg1")
            .push_tag(b"SM", &"sample_from_rg"),
    );
    let header_view = HeaderView::from_header(&header);

    let mut writer = Writer::from_path(path, &header, Format::Bam).expect("create fixture BAM");
    let sam_lines: [&[u8]; 4] = [
        b"read1\t0\tchr1\t1\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0\tMD:Z:5\tRG:Z:rg1",
        b"read2\t0\tchr1\t2\t60\t1S4M\t*\t0\t0\tTACGA\t!IIII\tNM:i:1\tMD:Z:2T1\tRG:Z:rg1",
        b"read3\t1024\tchr1\t3\t60\t5M\t*\t0\t0\tGGGGG\tIIIII\tNM:i:0\tMD:Z:5\tRG:Z:rg1",
        b"read4\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAA\tIIIII",
    ];

    for line in sam_lines {
        let record = Record::from_sam(&header_view, line).expect("valid SAM line");
        writer.write(&record).expect("write BAM record");
    }
}

fn resolve_binary_path() -> std::path::PathBuf {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_blammo-qc") {
        return path.into();
    }
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_blammo_qc") {
        return path.into();
    }

    // Fallback for environments that don't export CARGO_BIN_EXE_* during tests.
    let mut path = std::env::current_exe().expect("current exe");
    path.pop(); // test-binary filename
    path.pop(); // deps/
    path.push("blammo-qc");
    path
}

#[test]
fn cli_generates_expected_json_and_html_for_fixture_bam() {
    let tmp = tempdir().expect("create temp dir");
    let bam_path = tmp.path().join("fixture.bam");
    let json_path = tmp.path().join("report.json");
    let html_path = tmp.path().join("report.html");
    write_fixture_bam(&bam_path);

    let bin = resolve_binary_path();
    assert!(
        bin.exists(),
        "binary path does not exist: {}",
        bin.display()
    );

    let output = Command::new(&bin)
        .args([
            "--output-json",
            &json_path.to_string_lossy(),
            "--output-html",
            &html_path.to_string_lossy(),
            "--threads",
            "1",
            &bam_path.to_string_lossy(),
        ])
        .output()
        .expect("run blammo-qc");

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let json_text = fs::read_to_string(&json_path).expect("read JSON report");
    let json: serde_json::Value = serde_json::from_str(&json_text).expect("parse JSON report");

    assert_eq!(json["samples"].as_array().expect("samples array").len(), 1);
    let sample = &json["samples"][0];
    assert_eq!(sample["sample_id"], "sample_from_rg");

    assert_eq!(sample["counts"]["total_records_seen"], 4);
    assert_eq!(sample["counts"]["unmapped_reads"], 1);
    assert_eq!(sample["counts"]["primary_mapped_reads_used"], 2);
    assert_eq!(sample["counts"]["duplicate_reads_excluded"], 1);

    assert_eq!(sample["soft_clips"]["total_soft_clipped_bases"], 1);
    assert_eq!(sample["soft_clips"]["reads_with_soft_clips"], 1);

    assert_eq!(sample["mismatches"]["nm_sum"], 1);
    assert_eq!(sample["mismatches"]["nm_histogram"]["0"], 1);
    assert_eq!(sample["mismatches"]["nm_histogram"]["1"], 1);
    assert_eq!(sample["mismatches"]["by_base_quality_bin"]["40+"], 1);
    assert_eq!(sample["mismatches"]["by_base_quality_bin"]["30-39"], 0);

    assert_eq!(sample["read_length"]["min"], 5);
    assert_eq!(sample["read_length"]["max"], 5);
    assert_eq!(sample["read_length"]["mean"], 5.0);
    assert_eq!(sample["read_length"]["median"], 5.0);
    assert_eq!(sample["read_length"]["histogram"]["5"], 2);

    let genome_hist = sample["depth_distribution"]["genome_wide_histogram"]
        .as_array()
        .expect("depth histogram array");
    assert_eq!(genome_hist.len(), 32_768);
    assert_eq!(genome_hist[1], 1);
    assert_eq!(genome_hist[2], 4);
    assert_eq!(sample["depth_distribution"]["eligible_bases_count"], 5);
    assert_eq!(
        sample["depth_distribution"]["per_chromosome_histograms"]["chr1"][1],
        1
    );
    assert_eq!(
        sample["depth_distribution"]["per_chromosome_histograms"]["chr1"][2],
        4
    );

    assert_eq!(
        sample["warnings"].as_array().expect("warnings array").len(),
        0
    );

    let html = fs::read_to_string(&html_path).expect("read HTML report");
    assert!(html.contains("Blammo QC Report"));
    assert!(html.contains("plotly"));
}
