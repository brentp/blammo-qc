use std::collections::BTreeMap;

use chrono::{DateTime, Utc};
use serde::Serialize;

pub const DEPTH_HIST_BINS: usize = 2 * 16_384;
pub const READ_LENGTH_HIST_BINS: usize = 200_000;
pub const BQ_BIN_LABELS: [&str; 5] = ["0-9", "10-19", "20-29", "30-39", "40+"];

#[derive(Debug, Clone, Serialize)]
pub struct QcReport {
    pub tool_version: String,
    pub run_timestamp: DateTime<Utc>,
    pub settings: RunSettings,
    pub samples: Vec<SampleQc>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunSettings {
    pub min_base_quality: u8,
    pub min_mapping_quality: u8,
    pub threads: usize,
    pub reference: Option<String>,
    pub depth_scope: String,
    pub plot_max_contigs: usize,
    pub tag_bars: Vec<String>,
    pub tag_lines: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SampleQc {
    pub sample_id: String,
    pub input_path: String,
    pub counts: ReadCounts,
    pub soft_clips: SoftClipSummary,
    pub mismatches: MismatchSummary,
    pub read_length: ReadLengthSummary,
    pub depth_distribution: DepthSummary,
    pub tag_value_counts: BTreeMap<String, BTreeMap<i64, u64>>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct ReadCounts {
    pub total_records_seen: u64,
    pub unmapped_reads: u64,
    pub primary_mapped_reads_used: u64,
    pub duplicate_reads_excluded: u64,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct SoftClipSummary {
    pub total_soft_clipped_bases: u64,
    pub reads_with_soft_clips: u64,
    pub clipped_length_histogram: BTreeMap<u32, u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MismatchSummary {
    pub nm_sum: u64,
    pub nm_histogram: BTreeMap<u32, u64>,
    pub by_base_quality_bin: BTreeMap<String, u64>,
}

impl Default for MismatchSummary {
    fn default() -> Self {
        let by_base_quality_bin = BQ_BIN_LABELS
            .iter()
            .map(|label| ((*label).to_string(), 0_u64))
            .collect();
        Self {
            nm_sum: 0,
            nm_histogram: BTreeMap::new(),
            by_base_quality_bin,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct ReadLengthSummary {
    pub min: u32,
    pub max: u32,
    pub mean: f64,
    pub median: f64,
    pub p10: u32,
    pub p90: u32,
    pub histogram: BTreeMap<u32, u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DepthSummary {
    pub genome_wide_histogram: Vec<u64>,
    pub per_chromosome_histograms: BTreeMap<String, Vec<u64>>,
    pub overflow_bin_index: usize,
    pub eligible_bases_count: u64,
    pub depth_scope: String,
    pub reference_bases_total: u64,
    pub covered_bases_total: u64,
    pub zero_depth_bases_total: u64,
}

impl Default for DepthSummary {
    fn default() -> Self {
        Self {
            genome_wide_histogram: vec![0; DEPTH_HIST_BINS],
            per_chromosome_histograms: BTreeMap::new(),
            overflow_bin_index: DEPTH_HIST_BINS - 1,
            eligible_bases_count: 0,
            depth_scope: "covered_bases".to_string(),
            reference_bases_total: 0,
            covered_bases_total: 0,
            zero_depth_bases_total: 0,
        }
    }
}

pub fn base_quality_bin(quality: u8) -> usize {
    match quality {
        0..=9 => 0,
        10..=19 => 1,
        20..=29 => 2,
        30..=39 => 3,
        _ => 4,
    }
}
