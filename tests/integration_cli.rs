use std::{fs, process::Command};

use rust_htslib::bam::{Format, Header, HeaderView, Writer, header::HeaderRecord, record::Record};
use serde_json::json;
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
        b"read1\t0\tchr1\t1\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0\tMD:Z:5\tRG:Z:rg1\tZX:i:7",
        b"read2\t0\tchr1\t2\t60\t1S4M\t*\t0\t0\tTACGA\t!IIII\tNM:i:1\tMD:Z:2T1\tRG:Z:rg1\tZX:i:7",
        b"read3\t1024\tchr1\t3\t60\t5M\t*\t0\t0\tGGGGG\tIIIII\tNM:i:0\tMD:Z:5\tRG:Z:rg1\tZX:i:9",
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

#[test]
fn cli_tag_bar_counts_integer_values_and_adds_metric_views() {
    let tmp = tempdir().expect("create temp dir");
    let bam_path = tmp.path().join("fixture.bam");
    let json_path = tmp.path().join("report.json");
    let html_path = tmp.path().join("report.html");
    write_fixture_bam(&bam_path);

    let bin = resolve_binary_path();
    let output = Command::new(&bin)
        .args([
            "--output-json",
            &json_path.to_string_lossy(),
            "--output-html",
            &html_path.to_string_lossy(),
            "--threads",
            "1",
            "--tag-bar",
            "NM",
            "--tag-bar",
            "ZX",
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

    assert_eq!(json["settings"]["tag_bars"], json!(["NM", "ZX"]));
    let sample = &json["samples"][0];
    assert_eq!(sample["tag_value_counts"]["NM"]["0"], 1);
    assert_eq!(sample["tag_value_counts"]["NM"]["1"], 1);
    assert_eq!(sample["tag_value_counts"]["ZX"]["7"], 2);
    assert!(sample["tag_value_counts"]["ZX"]["9"].is_null());

    let html = fs::read_to_string(&html_path).expect("read HTML report");
    assert!(html.contains("NM tag values"));
    assert!(html.contains("ZX tag values"));
    assert!(html.contains("NM tag integer values by sample"));
}

#[test]
fn cli_tag_line_counts_integer_values_and_adds_line_metric_views() {
    let tmp = tempdir().expect("create temp dir");
    let bam_path = tmp.path().join("fixture.bam");
    let json_path = tmp.path().join("report.json");
    let html_path = tmp.path().join("report.html");
    write_fixture_bam(&bam_path);

    let bin = resolve_binary_path();
    let output = Command::new(&bin)
        .args([
            "--output-json",
            &json_path.to_string_lossy(),
            "--output-html",
            &html_path.to_string_lossy(),
            "--threads",
            "1",
            "--tag-line",
            "ZX",
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

    assert_eq!(json["settings"]["tag_lines"], json!(["ZX"]));
    let sample = &json["samples"][0];
    assert_eq!(sample["tag_value_counts"]["ZX"]["7"], 2);
    assert!(sample["tag_value_counts"]["ZX"]["9"].is_null());

    let html = fs::read_to_string(&html_path).expect("read HTML report");
    assert!(html.contains("ZX tag values (line)"));
    assert!(html.contains("ZX tag integer values by sample (line)"));
}

#[test]
fn cli_writes_downstream_data_json_for_fixture_bam() {
    let tmp = tempdir().expect("create temp dir");
    let bam_path = tmp.path().join("fixture.bam");
    let data_json_path = tmp.path().join("nested").join("report.data.json");
    let html_path = tmp.path().join("report.html");
    write_fixture_bam(&bam_path);

    let bin = resolve_binary_path();
    let output = Command::new(&bin)
        .args([
            "--output-data-json",
            &data_json_path.to_string_lossy(),
            "--output-html",
            &html_path.to_string_lossy(),
            "--threads",
            "1",
            "--tag-bar",
            "NM",
            "--tag-line",
            "ZX",
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

    let json_text = fs::read_to_string(&data_json_path).expect("read data JSON report");
    let json: serde_json::Value = serde_json::from_str(&json_text).expect("parse data JSON report");

    assert_eq!(json["data"]["sample_count"], 1);
    assert_eq!(json["data"]["samples"][0]["sample_id"], "sample_from_rg");
    assert_eq!(json["data"]["depth_views"][0]["key"], "genome_wide");
    assert_eq!(json["data"]["depth_views"][0]["kind"], "genome");
    assert_eq!(json["data"]["depth_views"][1]["key"], "chr1");
    assert!(
        !json["data"]["depth_views"][0]
            .as_object()
            .expect("genome depth view")
            .contains_key("omitted_contig_count")
    );
    assert!(
        !json["data"]["depth_views"][1]
            .as_object()
            .expect("contig depth view")
            .contains_key("omitted_contig_count")
    );
    assert_eq!(
        json["data"]["sample_depths"][0]["histograms"]["genome_wide"][1],
        1
    );
    assert_eq!(
        json["data"]["sample_depths"][0]["histograms"]["genome_wide"][2],
        4
    );
    assert_eq!(
        json["data"]["sample_depths"][0]["histograms"]["genome_wide"]
            .as_array()
            .expect("genome histogram")
            .len(),
        3
    );
    assert_eq!(json["data"]["sample_depths"][0]["histograms"]["chr1"][1], 1);

    let metrics = &json["data"]["metrics_table"][0];
    assert_eq!(metrics["sample_id"], "sample_from_rg");
    assert_eq!(metrics["median_depth"], 2);
    assert_eq!(metrics["mode_depth"], 2);
    assert_eq!(metrics["total_records_seen"], 4);
    assert_eq!(metrics["primary_mapped_reads_used"], 2);
    assert_eq!(metrics["primary_mapped_reads_percent"], 50.0);
    assert_eq!(metrics["unmapped_reads"], 1);
    assert_eq!(metrics["unmapped_reads_percent"], 25.0);
    assert_eq!(metrics["total_soft_clipped_bases"], 1);
    assert_eq!(metrics["soft_clipped_bases_per_million_bases"], 100000.0);
    assert_eq!(metrics["reads_with_soft_clips"], 1);
    assert_eq!(metrics["reads_with_soft_clips_percent"], 50.0);
    assert_eq!(metrics["nm_sum"], 1);
    assert_eq!(metrics["nm_per_million_bases"], 100000.0);
    assert_eq!(metrics["read_length_min"], 5);
    assert_eq!(metrics["read_length_p10"], 5);
    assert_eq!(metrics["read_length_median"], 5.0);
    assert_eq!(metrics["read_length_mean"], 5.0);
    assert_eq!(metrics["read_length_p90"], 5);
    assert_eq!(metrics["read_length_max"], 5);
    assert_eq!(metrics["warnings_count"], 0);
    assert_eq!(metrics["bq_mismatch_bins"]["40+"], 1);

    let tag_metrics = json["data"]["tag_metrics"]
        .as_array()
        .expect("tag metrics array");
    assert_eq!(tag_metrics.len(), 2);
    assert_eq!(tag_metrics[0]["tag"], "NM");
    assert_eq!(tag_metrics[0]["view"], "bar");
    assert_eq!(tag_metrics[0]["values"], json!([0, 1]));
    assert_eq!(tag_metrics[0]["sample_counts"][0]["counts"]["0"], 1);
    assert_eq!(tag_metrics[0]["sample_counts"][0]["counts"]["1"], 1);
    assert_eq!(tag_metrics[1]["tag"], "ZX");
    assert_eq!(tag_metrics[1]["view"], "line");
    assert_eq!(tag_metrics[1]["values"], json!([7]));
    assert_eq!(tag_metrics[1]["sample_counts"][0]["counts"]["7"], 2);
}

#[test]
fn cli_skips_json_output_when_flag_not_provided() {
    let tmp = tempdir().expect("create temp dir");
    let bam_path = tmp.path().join("fixture.bam");
    let html_path = tmp.path().join("report.html");
    let implicit_json_path = tmp.path().join("qc.json");
    let implicit_data_json_path = tmp.path().join("data.json");
    write_fixture_bam(&bam_path);

    let bin = resolve_binary_path();
    assert!(
        bin.exists(),
        "binary path does not exist: {}",
        bin.display()
    );

    let output = Command::new(&bin)
        .current_dir(tmp.path())
        .args([
            "--output-html",
            "report.html",
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

    assert!(html_path.exists(), "expected HTML report to be written");
    assert!(
        !implicit_json_path.exists(),
        "did not expect default JSON output at {}",
        implicit_json_path.display()
    );
    assert!(
        !implicit_data_json_path.exists(),
        "did not expect default data JSON output at {}",
        implicit_data_json_path.display()
    );
}
