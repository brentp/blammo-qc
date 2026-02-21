use std::{
    cmp::Reverse,
    collections::{BTreeMap, BinaryHeap},
    path::Path,
};

use anyhow::{Context, Result};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar},
    Read, Record,
};

use crate::{
    cli::{is_cram_path, DepthScope, SampleInput},
    model::{
        base_quality_bin, DepthSummary, MismatchSummary, ReadCounts, ReadLengthSummary, SampleQc,
        SoftClipSummary, BQ_BIN_LABELS, DEPTH_HIST_BINS, READ_LENGTH_HIST_BINS,
    },
};

#[derive(Debug, Clone)]
pub struct ProcessingOptions<'a> {
    pub reference: Option<&'a Path>,
    pub min_base_quality: u8,
    pub min_mapping_quality: u8,
    pub include_duplicates: bool,
    pub depth_scope: DepthScope,
}

#[derive(Debug, Clone, Copy)]
struct Interval {
    start: u32,
    end: u32,
}

#[derive(Debug, Default)]
struct MismatchBins {
    by_bq: [u64; 5],
    skipped_offsets: u64,
}

#[derive(Debug, Default)]
struct WarningCounters {
    missing_nm: u64,
    missing_md: u64,
    malformed_md: u64,
    md_offsets_without_query_base: u64,
}

#[derive(Debug)]
struct ReadLengthAccumulator {
    bins: Vec<u64>,
    count: u64,
    sum: u64,
    min: u32,
    max: u32,
}

impl ReadLengthAccumulator {
    fn new() -> Self {
        Self {
            bins: vec![0_u64; READ_LENGTH_HIST_BINS],
            count: 0,
            sum: 0,
            min: u32::MAX,
            max: 0,
        }
    }

    fn push(&mut self, read_len: u32) {
        let index = (read_len as usize).min(READ_LENGTH_HIST_BINS - 1);
        self.bins[index] += 1;
        self.count += 1;
        self.sum += read_len as u64;
        self.min = self.min.min(read_len);
        self.max = self.max.max(read_len);
    }

    fn into_summary(self) -> ReadLengthSummary {
        if self.count == 0 {
            return ReadLengthSummary::default();
        }

        let mut histogram = BTreeMap::new();
        for (len, count) in self.bins.iter().enumerate() {
            if *count > 0 {
                histogram.insert(len as u32, *count);
            }
        }

        let median = if self.count % 2 == 1 {
            value_at_rank(&self.bins, self.count / 2) as f64
        } else {
            let left = value_at_rank(&self.bins, (self.count / 2) - 1) as f64;
            let right = value_at_rank(&self.bins, self.count / 2) as f64;
            (left + right) / 2.0
        };

        ReadLengthSummary {
            min: self.min,
            max: self.max,
            mean: self.sum as f64 / self.count as f64,
            median,
            p10: percentile_from_bins(&self.bins, self.count, 0.10),
            p90: percentile_from_bins(&self.bins, self.count, 0.90),
            histogram,
        }
    }
}

#[derive(Debug)]
struct DepthCollector {
    target_names: Vec<String>,
    target_lengths: Vec<u64>,
    depth_scope: DepthScope,
    summary: DepthSummary,
    processed_tids: Vec<bool>,
    current_tid: Option<usize>,
    current_pos: u32,
    current_covered_bases: u64,
    current_histogram: Vec<u64>,
    active_ends: BinaryHeap<Reverse<u32>>,
}

impl DepthCollector {
    fn new(target_names: Vec<String>, target_lengths: Vec<u64>, depth_scope: DepthScope) -> Self {
        let mut summary = DepthSummary::default();
        summary.depth_scope = depth_scope.as_str().to_string();
        let processed_tids = vec![false; target_names.len()];
        Self {
            target_names,
            target_lengths,
            depth_scope,
            summary,
            processed_tids,
            current_tid: None,
            current_pos: 0,
            current_covered_bases: 0,
            current_histogram: Vec::new(),
            active_ends: BinaryHeap::new(),
        }
    }

    fn push_intervals(&mut self, tid: usize, intervals: &[Interval]) -> Result<()> {
        self.enter_tid(tid)?;
        for interval in intervals {
            self.push_interval(*interval)?;
        }
        Ok(())
    }

    fn finish(mut self) -> DepthSummary {
        self.finalize_current_tid();

        if matches!(self.depth_scope, DepthScope::ReferenceBases) {
            for tid in 0..self.target_names.len() {
                if self.processed_tids[tid] {
                    continue;
                }
                let reference_bases = self.target_lengths.get(tid).copied().unwrap_or(0);
                if reference_bases == 0 {
                    continue;
                }
                self.summary.reference_bases_total += reference_bases;
                self.summary.genome_wide_histogram[0] += reference_bases;
                let contig_name = self
                    .target_names
                    .get(tid)
                    .cloned()
                    .unwrap_or_else(|| format!("tid_{tid}"));
                if let Some(canonical_name) = canonical_depth_chromosome(&contig_name) {
                    let mut histogram = vec![0_u64; DEPTH_HIST_BINS];
                    histogram[0] = reference_bases;
                    merge_per_chrom_histogram(
                        &mut self.summary.per_chromosome_histograms,
                        canonical_name,
                        histogram,
                    );
                }
            }
        }

        self.summary.zero_depth_bases_total = self.summary.genome_wide_histogram[0];
        self.summary
    }

    fn enter_tid(&mut self, tid: usize) -> Result<()> {
        if tid >= self.target_names.len() {
            return Ok(());
        }
        match self.current_tid {
            None => {
                self.start_tid(tid);
                Ok(())
            }
            Some(current) if current == tid => Ok(()),
            Some(current) if tid > current => {
                self.finalize_current_tid();
                self.start_tid(tid);
                Ok(())
            }
            Some(current) => anyhow::bail!(
                "depth streaming requires coordinate-sorted input; encountered tid {} after tid {}",
                tid,
                current
            ),
        }
    }

    fn start_tid(&mut self, tid: usize) {
        self.current_tid = Some(tid);
        self.processed_tids[tid] = true;
        self.current_pos = 0;
        self.current_covered_bases = 0;
        self.current_histogram = vec![0_u64; DEPTH_HIST_BINS];
        self.active_ends.clear();
    }

    fn push_interval(&mut self, interval: Interval) -> Result<()> {
        if interval.end <= interval.start {
            return Ok(());
        }
        if interval.start < self.current_pos {
            anyhow::bail!(
                "depth streaming requires coordinate-sorted input within contig; saw start {} after {}",
                interval.start,
                self.current_pos
            );
        }
        self.advance_to(interval.start);
        self.active_ends.push(Reverse(interval.end));
        Ok(())
    }

    fn advance_to(&mut self, target: u32) {
        while self.current_pos < target {
            self.pop_ended();
            let next_boundary = self
                .active_ends
                .peek()
                .map(|Reverse(end)| (*end).min(target))
                .unwrap_or(target);
            if next_boundary <= self.current_pos {
                break;
            }
            let span = (next_boundary - self.current_pos) as u64;
            let depth = self.active_ends.len();
            if depth > 0 {
                self.add_span(depth, span);
            } else if matches!(self.depth_scope, DepthScope::ReferenceBases) {
                self.add_span(0, span);
            }
            self.current_pos = next_boundary;
        }
    }

    fn pop_ended(&mut self) {
        while let Some(Reverse(end)) = self.active_ends.peek() {
            if *end <= self.current_pos {
                self.active_ends.pop();
            } else {
                break;
            }
        }
    }

    fn drain_active_depth(&mut self) {
        loop {
            self.pop_ended();
            let Some(next_end) = self.active_ends.peek().map(|Reverse(end)| *end) else {
                break;
            };
            if next_end <= self.current_pos {
                continue;
            }
            let span = (next_end - self.current_pos) as u64;
            self.add_span(self.active_ends.len(), span);
            self.current_pos = next_end;
        }
        self.pop_ended();
    }

    fn add_span(&mut self, depth: usize, span: u64) {
        if span == 0 {
            return;
        }
        let bin = depth.min(DEPTH_HIST_BINS - 1);
        self.summary.genome_wide_histogram[bin] += span;
        if !self.current_histogram.is_empty() {
            self.current_histogram[bin] += span;
        }
        if depth > 0 {
            self.summary.eligible_bases_count += span;
            self.summary.covered_bases_total += span;
            self.current_covered_bases += span;
        }
    }

    fn finalize_current_tid(&mut self) {
        let Some(tid) = self.current_tid.take() else {
            return;
        };
        self.drain_active_depth();

        let reference_bases = self.target_lengths.get(tid).copied().unwrap_or(0);
        if matches!(self.depth_scope, DepthScope::ReferenceBases) {
            let consumed = self.current_pos as u64;
            if reference_bases > consumed {
                self.add_span(0, reference_bases - consumed);
            }
            self.summary.reference_bases_total += reference_bases.max(self.current_covered_bases);
        } else {
            self.summary.reference_bases_total += self.current_covered_bases;
        }

        let include_contig = match self.depth_scope {
            DepthScope::CoveredBases => self.current_covered_bases > 0,
            DepthScope::ReferenceBases => reference_bases > 0 || self.current_covered_bases > 0,
        };
        if include_contig {
            let contig_name = self
                .target_names
                .get(tid)
                .cloned()
                .unwrap_or_else(|| format!("tid_{tid}"));
            if let Some(canonical_name) = canonical_depth_chromosome(&contig_name) {
                merge_per_chrom_histogram(
                    &mut self.summary.per_chromosome_histograms,
                    canonical_name,
                    std::mem::take(&mut self.current_histogram),
                );
            }
        } else {
            self.current_histogram.clear();
        }

        self.current_pos = 0;
        self.current_covered_bases = 0;
        self.active_ends.clear();
    }
}

pub fn process_sample(input: &SampleInput, opts: &ProcessingOptions<'_>) -> Result<SampleQc> {
    let mut reader = bam::Reader::from_path(&input.path)
        .with_context(|| format!("failed to open alignment file {}", input.path.display()))?;

    if is_cram_path(&input.path) {
        let reference = opts
            .reference
            .context("CRAM input requires --reference but none was provided")?;
        reader
            .set_reference(reference)
            .with_context(|| format!("failed to set reference {}", reference.display()))?;
    }
    reader
        .set_threads(2)
        .context("failed to set alignment reader threads to 2")?;

    let target_names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();
    let target_lengths: Vec<u64> = (0..target_names.len())
        .map(|tid| reader.header().target_len(tid as u32).unwrap_or(0))
        .collect();
    let mut depth_collector = DepthCollector::new(target_names, target_lengths, opts.depth_scope);

    let mut counts = ReadCounts::default();
    let mut soft_clips = SoftClipSummary::default();
    let mut mismatch_summary = MismatchSummary::default();
    let mut mismatch_bins = MismatchBins::default();
    let mut warning_counts = WarningCounters::default();

    let mut read_length_acc = ReadLengthAccumulator::new();
    let mut depth_interval_buf: Vec<Interval> = Vec::new();

    for record_result in reader.records() {
        let record = record_result.with_context(|| {
            format!(
                "error while reading records from {}",
                input.path.to_string_lossy()
            )
        })?;
        counts.total_records_seen += 1;

        let classification = classify_record(&record, opts.include_duplicates);
        if classification.is_unmapped {
            counts.unmapped_reads += 1;
        }
        if classification.duplicate_excluded {
            counts.duplicate_reads_excluded += 1;
            continue;
        }
        if !classification.include_in_read_metrics {
            continue;
        }
        counts.primary_mapped_reads_used += 1;

        update_soft_clip_metrics(&record, &mut soft_clips);

        let read_len = record.seq_len() as u32;
        read_length_acc.push(read_len);

        if let Some(nm) = extract_nm_value(&record) {
            mismatch_summary.nm_sum += nm as u64;
            *mismatch_summary.nm_histogram.entry(nm).or_insert(0) += 1;
        } else {
            warning_counts.missing_nm += 1;
        }

        match extract_md_string(&record) {
            Some(md) => match mismatch_bins_from_md(&record, md) {
                Ok(md_bins) => {
                    for (i, count) in md_bins.by_bq.iter().enumerate() {
                        mismatch_bins.by_bq[i] += count;
                    }
                    warning_counts.md_offsets_without_query_base += md_bins.skipped_offsets;
                }
                Err(_) => {
                    warning_counts.malformed_md += 1;
                }
            },
            None => {
                warning_counts.missing_md += 1;
            }
        }

        if record.mapq() < opts.min_mapping_quality {
            continue;
        }
        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let tid = tid as usize;
        depth_intervals_from_record(&record, opts.min_base_quality, &mut depth_interval_buf);
        depth_collector
            .push_intervals(tid, &depth_interval_buf)
            .with_context(|| {
                format!(
                    "failed depth streaming for {} (input must be coordinate-sorted)",
                    input.path.display()
                )
            })?;
    }

    for (index, label) in BQ_BIN_LABELS.iter().enumerate() {
        mismatch_summary
            .by_base_quality_bin
            .insert((*label).to_string(), mismatch_bins.by_bq[index]);
    }

    let read_length = read_length_acc.into_summary();
    let depth_distribution = depth_collector.finish();
    let warnings = build_warnings(warning_counts);

    Ok(SampleQc {
        sample_id: input.sample_id.clone(),
        input_path: input.path.to_string_lossy().to_string(),
        counts,
        soft_clips,
        mismatches: mismatch_summary,
        read_length,
        depth_distribution,
        warnings,
    })
}

#[derive(Debug, Clone, Copy)]
struct RecordClassification {
    is_unmapped: bool,
    duplicate_excluded: bool,
    include_in_read_metrics: bool,
}

fn classify_record(record: &Record, include_duplicates: bool) -> RecordClassification {
    let is_unmapped = record.is_unmapped();
    let is_primary_mapped = !is_unmapped && !record.is_secondary() && !record.is_supplementary();
    let duplicate_excluded = is_primary_mapped && record.is_duplicate() && !include_duplicates;
    let include_in_read_metrics = is_primary_mapped && !duplicate_excluded;
    RecordClassification {
        is_unmapped,
        duplicate_excluded,
        include_in_read_metrics,
    }
}

fn update_soft_clip_metrics(record: &Record, soft_clips: &mut SoftClipSummary) {
    let mut clipped_bases: u32 = 0;
    for op in record.cigar().iter() {
        if let Cigar::SoftClip(len) = *op {
            clipped_bases += len;
        }
    }

    if clipped_bases > 0 {
        soft_clips.total_soft_clipped_bases += clipped_bases as u64;
        soft_clips.reads_with_soft_clips += 1;
        *soft_clips
            .clipped_length_histogram
            .entry(clipped_bases)
            .or_insert(0) += 1;
    }
}

fn extract_nm_value(record: &Record) -> Option<u32> {
    match record.aux(b"NM").ok()? {
        Aux::I8(v) if v >= 0 => Some(v as u32),
        Aux::U8(v) => Some(v as u32),
        Aux::I16(v) if v >= 0 => Some(v as u32),
        Aux::U16(v) => Some(v as u32),
        Aux::I32(v) if v >= 0 => Some(v as u32),
        Aux::U32(v) => Some(v),
        _ => None,
    }
}

fn extract_md_string(record: &Record) -> Option<&str> {
    match record.aux(b"MD").ok()? {
        Aux::String(value) => Some(value),
        _ => None,
    }
}

fn mismatch_bins_from_md(record: &Record, md: &str) -> Result<MismatchBins> {
    let mismatch_offsets = parse_md_mismatch_offsets(md)?;
    let quals = record.qual();

    let mut bins = MismatchBins::default();
    if mismatch_offsets.is_empty() {
        return Ok(bins);
    }

    let mut mismatch_index = 0_usize;
    let mut ref_offset = 0_u32;
    let mut query_pos = 0_usize;

    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let op_end = ref_offset.saturating_add(len);
                while mismatch_index < mismatch_offsets.len()
                    && mismatch_offsets[mismatch_index] < op_end
                {
                    let rel = (mismatch_offsets[mismatch_index] - ref_offset) as usize;
                    let qpos = query_pos + rel;
                    if qpos < quals.len() {
                        bins.by_bq[base_quality_bin(quals[qpos])] += 1;
                    } else {
                        bins.skipped_offsets += 1;
                    }
                    mismatch_index += 1;
                }
                ref_offset = op_end;
                query_pos += len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                let op_end = ref_offset.saturating_add(len);
                while mismatch_index < mismatch_offsets.len()
                    && mismatch_offsets[mismatch_index] < op_end
                {
                    bins.skipped_offsets += 1;
                    mismatch_index += 1;
                }
                ref_offset = op_end;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                query_pos += len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    bins.skipped_offsets += (mismatch_offsets.len() - mismatch_index) as u64;

    Ok(bins)
}

fn parse_md_mismatch_offsets(md: &str) -> Result<Vec<u32>> {
    let bytes = md.as_bytes();
    let mut offsets = Vec::new();
    let mut i = 0_usize;
    let mut ref_offset = 0_u32;

    while i < bytes.len() {
        let current = bytes[i];
        if current.is_ascii_digit() {
            let mut value = 0_u32;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                value = value
                    .checked_mul(10)
                    .and_then(|v| v.checked_add((bytes[i] - b'0') as u32))
                    .context("overflow while parsing MD numeric run")?;
                i += 1;
            }
            ref_offset = ref_offset
                .checked_add(value)
                .context("overflow while advancing MD reference offset")?;
            continue;
        }

        if current == b'^' {
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                ref_offset = ref_offset
                    .checked_add(1)
                    .context("overflow while parsing MD deletion block")?;
                i += 1;
            }
            continue;
        }

        if current.is_ascii_alphabetic() {
            offsets.push(ref_offset);
            ref_offset = ref_offset
                .checked_add(1)
                .context("overflow while parsing MD mismatch offset")?;
            i += 1;
            continue;
        }

        anyhow::bail!("unexpected character '{}' in MD tag", current as char);
    }

    Ok(offsets)
}

fn depth_intervals_from_record(
    record: &Record,
    min_base_quality: u8,
    intervals: &mut Vec<Interval>,
) {
    intervals.clear();
    if record.pos() < 0 {
        return;
    }

    let mut ref_pos = record.pos() as u32;
    let mut query_pos: usize = 0;
    let qualities = record.qual();
    let mut open_start: Option<u32> = None;

    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for offset in 0..len as usize {
                    let qpos = query_pos + offset;
                    if qpos >= qualities.len() {
                        break;
                    }
                    let curr_ref = ref_pos + offset as u32;
                    if qualities[qpos] >= min_base_quality {
                        if open_start.is_none() {
                            open_start = Some(curr_ref);
                        }
                    } else if let Some(start) = open_start.take() {
                        if start < curr_ref {
                            intervals.push(Interval {
                                start,
                                end: curr_ref,
                            });
                        }
                    }
                }
                ref_pos += len;
                query_pos += len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if let Some(start) = open_start.take() {
                    if start < ref_pos {
                        intervals.push(Interval {
                            start,
                            end: ref_pos,
                        });
                    }
                }
                ref_pos += len;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                query_pos += len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    if let Some(start) = open_start {
        if start < ref_pos {
            intervals.push(Interval {
                start,
                end: ref_pos,
            });
        }
    }
}

fn percentile_from_bins(bins: &[u64], total_count: u64, fraction: f64) -> u32 {
    if total_count == 0 {
        return 0;
    }
    let rank = ((total_count - 1) as f64 * fraction).round() as u64;
    value_at_rank(bins, rank)
}

fn value_at_rank(bins: &[u64], rank: u64) -> u32 {
    let mut cumulative = 0_u64;
    for (index, count) in bins.iter().enumerate() {
        cumulative += *count;
        if cumulative > rank {
            return index as u32;
        }
    }
    (bins.len().saturating_sub(1)) as u32
}

fn merge_per_chrom_histogram(
    per_chrom: &mut BTreeMap<String, Vec<u64>>,
    chromosome: String,
    histogram: Vec<u64>,
) {
    if let Some(existing) = per_chrom.get_mut(&chromosome) {
        for (i, value) in histogram.into_iter().enumerate().take(DEPTH_HIST_BINS) {
            existing[i] += value;
        }
    } else {
        per_chrom.insert(chromosome, histogram);
    }
}

fn canonical_depth_chromosome(name: &str) -> Option<String> {
    let raw = name.trim();
    let without_chr = raw.strip_prefix("chr").unwrap_or(raw);

    if let Ok(num) = without_chr.parse::<u8>() {
        if (1..=22).contains(&num) {
            return Some(format!("chr{num}"));
        }
        return None;
    }

    if without_chr.eq_ignore_ascii_case("x") {
        return Some("chrX".to_string());
    }
    if without_chr.eq_ignore_ascii_case("y") {
        return Some("chrY".to_string());
    }
    None
}

fn build_warnings(counts: WarningCounters) -> Vec<String> {
    let mut warnings = Vec::new();
    if counts.missing_nm > 0 {
        warnings.push(format!(
            "NM tag missing on {} primary mapped reads; NM totals skip those reads.",
            counts.missing_nm
        ));
    }
    if counts.missing_md > 0 {
        warnings.push(format!(
            "MD tag missing on {} primary mapped reads; mismatch-by-base-quality bins skip those reads.",
            counts.missing_md
        ));
    }
    if counts.malformed_md > 0 {
        warnings.push(format!(
            "Malformed MD tag detected on {} reads; mismatch-by-base-quality bins skip those reads.",
            counts.malformed_md
        ));
    }
    if counts.md_offsets_without_query_base > 0 {
        warnings.push(format!(
            "{} MD mismatch positions could not be mapped to query bases and were skipped.",
            counts.md_offsets_without_query_base
        ));
    }
    warnings
}

#[cfg(test)]
mod tests {
    use super::{classify_record, parse_md_mismatch_offsets, DepthCollector, Interval};
    use crate::cli::DepthScope;

    #[test]
    fn md_parser_handles_mismatches_and_deletions() {
        let offsets = parse_md_mismatch_offsets("10A5^CC3T0G").unwrap();
        assert_eq!(offsets, vec![10, 21, 22]);
    }

    #[test]
    fn streaming_depth_collector_tracks_depth_without_interval_buffering() {
        let mut collector = DepthCollector::new(
            vec!["chr1".to_string()],
            vec![100],
            DepthScope::CoveredBases,
        );
        collector
            .push_intervals(
                0,
                &[
                    Interval { start: 10, end: 20 },
                    Interval { start: 10, end: 20 },
                    Interval { start: 10, end: 20 },
                ],
            )
            .unwrap();
        let summary = collector.finish();
        assert_eq!(summary.covered_bases_total, 10);
        assert_eq!(summary.genome_wide_histogram[3], 10);
    }

    #[test]
    fn depth_collector_keeps_only_chr1_to_22_x_y() {
        let mut collector = DepthCollector::new(
            vec!["chr1".to_string(), "chrM".to_string(), "1".to_string()],
            vec![100, 100, 100],
            DepthScope::CoveredBases,
        );
        collector
            .push_intervals(0, &[Interval { start: 10, end: 20 }])
            .unwrap();
        collector
            .push_intervals(1, &[Interval { start: 10, end: 20 }])
            .unwrap();
        collector
            .push_intervals(2, &[Interval { start: 10, end: 20 }])
            .unwrap();
        let summary = collector.finish();
        assert_eq!(summary.per_chromosome_histograms.len(), 1);
        assert!(summary.per_chromosome_histograms.contains_key("chr1"));
        assert_eq!(summary.per_chromosome_histograms["chr1"][1], 20);
    }

    #[test]
    fn record_classification_respects_duplicate_policy() {
        let mut record = rust_htslib::bam::Record::new();
        record.set_flags(0); // primary mapped
        let no_dups = classify_record(&record, false);
        assert!(no_dups.include_in_read_metrics);

        record.set_flags(0x400); // duplicate
        let filtered = classify_record(&record, false);
        assert!(filtered.duplicate_excluded);
        assert!(!filtered.include_in_read_metrics);

        let included = classify_record(&record, true);
        assert!(included.include_in_read_metrics);
    }
}
