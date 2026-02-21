use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::BufWriter,
    path::Path,
};

use anyhow::{Context, Result};
use plotly::{
    common::{Mode, Visible},
    layout::{
        update_menu::{Button, ButtonMethod, UpdateMenu, UpdateMenuDirection, UpdateMenuType},
        Axis, BarMode,
    },
    Bar, Configuration, Layout, Plot, Scatter,
};
use serde_json::json;

use crate::model::{QcReport, SampleQc, BQ_BIN_LABELS, DEPTH_HIST_BINS};

pub fn write_json_report(report: &QcReport, path: &Path) -> Result<()> {
    let filtered = filter_json_chromosomes(report);
    let file = File::create(path)
        .with_context(|| format!("failed to create JSON output {}", path.display()))?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &filtered)
        .with_context(|| format!("failed to write JSON output {}", path.display()))?;
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum TraceView {
    SoftClips,
    Unmapped,
    MismatchByBaseQuality,
    ReadLength,
    Depth { chromosome: String },
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ViewSelection {
    SoftClips,
    Unmapped,
    MismatchByBaseQuality,
    ReadLength,
    Depth { chromosome: String },
}

pub fn write_html_report(report: &QcReport, path: &Path, plot_max_contigs: usize) -> Result<()> {
    let mut plot = Plot::new();
    let mut trace_views: Vec<TraceView> = Vec::new();
    let sample_ids: Vec<String> = report
        .samples
        .iter()
        .map(|sample| sample.sample_id.clone())
        .collect();

    let soft_clip_values: Vec<u64> = report
        .samples
        .iter()
        .map(|sample| sample.soft_clips.total_soft_clipped_bases)
        .collect();
    plot.add_trace(
        Bar::new(sample_ids.clone(), soft_clip_values)
            .name("Soft-clipped bases")
            .visible(Visible::True),
    );
    trace_views.push(TraceView::SoftClips);

    let unmapped_values: Vec<u64> = report
        .samples
        .iter()
        .map(|sample| sample.counts.unmapped_reads)
        .collect();
    plot.add_trace(
        Bar::new(sample_ids.clone(), unmapped_values)
            .name("Unmapped reads")
            .visible(Visible::False),
    );
    trace_views.push(TraceView::Unmapped);

    for sample in &report.samples {
        let y: Vec<u64> = BQ_BIN_LABELS
            .iter()
            .map(|label| sample_count(sample, label))
            .collect();
        let x: Vec<&str> = BQ_BIN_LABELS.to_vec();
        plot.add_trace(
            Scatter::new(x, y)
                .mode(Mode::LinesMarkers)
                .name(format!("{} mismatch BQ", sample.sample_id))
                .visible(Visible::False),
        );
        trace_views.push(TraceView::MismatchByBaseQuality);
    }

    for sample in &report.samples {
        let (x, y) = histogram_points(&sample.read_length.histogram);
        plot.add_trace(
            Scatter::new(x, y)
                .mode(Mode::Lines)
                .name(format!("{} read length", sample.sample_id))
                .visible(Visible::False),
        );
        trace_views.push(TraceView::ReadLength);
    }

    let (selected_contigs, omitted_contig_count) = select_depth_contigs(report, plot_max_contigs);
    let selected_contig_set: BTreeSet<String> = selected_contigs.iter().cloned().collect();
    let other_histograms: Vec<Vec<u64>> = if omitted_contig_count > 0 {
        report
            .samples
            .iter()
            .map(|sample| combine_other_contig_histogram(sample, &selected_contig_set))
            .collect()
    } else {
        Vec::new()
    };

    let mut depth_views = vec!["__genome__".to_string()];
    depth_views.extend(selected_contigs.iter().cloned());
    if omitted_contig_count > 0 {
        depth_views.push("__other__".to_string());
    }

    for chromosome in &depth_views {
        for (sample_index, sample) in report.samples.iter().enumerate() {
            let histogram = if chromosome == "__genome__" {
                sample.depth_distribution.genome_wide_histogram.as_slice()
            } else if chromosome == "__other__" {
                other_histograms
                    .get(sample_index)
                    .map(Vec::as_slice)
                    .unwrap_or(&[])
            } else {
                sample
                    .depth_distribution
                    .per_chromosome_histograms
                    .get(chromosome)
                    .map(Vec::as_slice)
                    .unwrap_or(&[])
            };
            let (x, y) = depth_points(histogram);
            let label = if chromosome == "__genome__" {
                format!("{} depth (genome)", sample.sample_id)
            } else if chromosome == "__other__" {
                format!("{} depth (other contigs)", sample.sample_id)
            } else {
                format!("{} depth ({chromosome})", sample.sample_id)
            };
            plot.add_trace(
                Scatter::new(x, y)
                    .mode(Mode::Lines)
                    .name(label)
                    .visible(Visible::False),
            );
            trace_views.push(TraceView::Depth {
                chromosome: chromosome.clone(),
            });
        }
    }

    let mut view_options = vec![
        (
            "Soft clips".to_string(),
            ViewSelection::SoftClips,
            "Soft-clipped bases by sample".to_string(),
            "Sample".to_string(),
            "Soft-clipped bases".to_string(),
        ),
        (
            "Unmapped reads".to_string(),
            ViewSelection::Unmapped,
            "Unmapped reads by sample".to_string(),
            "Sample".to_string(),
            "Reads".to_string(),
        ),
        (
            "Mismatches by base quality".to_string(),
            ViewSelection::MismatchByBaseQuality,
            "Mismatch counts partitioned by base quality".to_string(),
            "Base-quality bin".to_string(),
            "Mismatches".to_string(),
        ),
        (
            "Read length distribution".to_string(),
            ViewSelection::ReadLength,
            "Read-length distribution".to_string(),
            "Read length".to_string(),
            "Reads".to_string(),
        ),
    ];

    view_options.push((
        "Depth distribution (genome)".to_string(),
        ViewSelection::Depth {
            chromosome: "__genome__".to_string(),
        },
        format!(
            "Depth distribution (genome-wide, {})",
            depth_scope_summary(report)
        ),
        "Depth".to_string(),
        depth_y_axis_label(report),
    ));
    for chromosome in depth_views
        .iter()
        .filter(|v| v.as_str() != "__genome__" && v.as_str() != "__other__")
    {
        view_options.push((
            format!("Depth distribution ({chromosome})"),
            ViewSelection::Depth {
                chromosome: chromosome.clone(),
            },
            format!(
                "Depth distribution ({chromosome}, {})",
                depth_scope_summary(report)
            ),
            "Depth".to_string(),
            depth_y_axis_label(report),
        ));
    }
    if omitted_contig_count > 0 {
        view_options.push((
            format!("Depth distribution (Other {omitted_contig_count} contigs)"),
            ViewSelection::Depth {
                chromosome: "__other__".to_string(),
            },
            format!(
                "Depth distribution (other contigs, {}, n={omitted_contig_count})",
                depth_scope_summary(report)
            ),
            "Depth".to_string(),
            depth_y_axis_label(report),
        ));
    }

    let metric_buttons: Vec<Button> = view_options
        .iter()
        .map(|(label, selection, title, x_label, y_label)| {
            let visible: Vec<bool> = trace_views
                .iter()
                .map(|trace_view| trace_is_visible(trace_view, selection))
                .collect();
            Button::new()
                .label(label.clone())
                .method(ButtonMethod::Update)
                .args(json!([
                    { "visible": visible },
                    {
                        "title": { "text": title },
                        "xaxis": { "title": { "text": x_label } },
                        "yaxis": { "title": { "text": y_label } }
                    }
                ]))
        })
        .collect();

    let y_scale_buttons = vec![
        Button::new()
            .label("Y: Linear")
            .method(ButtonMethod::Relayout)
            .args(json!([{ "yaxis.type": "linear" }])),
        Button::new()
            .label("Y: Log")
            .method(ButtonMethod::Relayout)
            .args(json!([{ "yaxis.type": "log" }])),
    ];

    let metric_menu = UpdateMenu::new()
        .ty(UpdateMenuType::Dropdown)
        .buttons(metric_buttons)
        .active(0)
        .x(0.0)
        .y(1.20);
    let y_scale_menu = UpdateMenu::new()
        .ty(UpdateMenuType::Buttons)
        .direction(UpdateMenuDirection::Right)
        .buttons(y_scale_buttons)
        .active(0)
        .x(0.55)
        .y(1.20);

    let layout = Layout::new()
        .title("Blammo QC report")
        .bar_mode(BarMode::Group)
        .x_axis(Axis::new().title("Sample"))
        .y_axis(Axis::new().title("Soft-clipped bases"))
        .update_menus(vec![metric_menu, y_scale_menu]);
    plot.set_layout(layout);
    plot.set_configuration(Configuration::new().responsive(true));
    plot.write_html(path);
    Ok(())
}

fn sample_count(sample: &SampleQc, label: &str) -> u64 {
    sample
        .mismatches
        .by_base_quality_bin
        .get(label)
        .copied()
        .unwrap_or(0)
}

fn trace_is_visible(trace_view: &TraceView, selection: &ViewSelection) -> bool {
    match (trace_view, selection) {
        (TraceView::SoftClips, ViewSelection::SoftClips) => true,
        (TraceView::Unmapped, ViewSelection::Unmapped) => true,
        (TraceView::MismatchByBaseQuality, ViewSelection::MismatchByBaseQuality) => true,
        (TraceView::ReadLength, ViewSelection::ReadLength) => true,
        (
            TraceView::Depth {
                chromosome: trace_chromosome,
            },
            ViewSelection::Depth {
                chromosome: selected_chromosome,
            },
        ) => trace_chromosome == selected_chromosome,
        _ => false,
    }
}

fn histogram_points(hist: &std::collections::BTreeMap<u32, u64>) -> (Vec<u32>, Vec<u64>) {
    let mut x = Vec::with_capacity(hist.len());
    let mut y = Vec::with_capacity(hist.len());
    for (read_len, count) in hist {
        x.push(*read_len);
        y.push(*count);
    }
    (x, y)
}

fn depth_points(hist: &[u64]) -> (Vec<u32>, Vec<u64>) {
    if hist.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let upper = usize::min(hist.len(), DEPTH_HIST_BINS);
    let max_nonzero = hist
        .iter()
        .take(upper)
        .enumerate()
        .skip(1)
        .rfind(|(_, count)| **count > 0)
        .map(|(depth, _)| depth);
    let Some(max_nonzero) = max_nonzero else {
        return (Vec::new(), Vec::new());
    };

    let mut x = Vec::with_capacity(max_nonzero);
    let mut y = Vec::with_capacity(max_nonzero);
    for depth in 1..=max_nonzero {
        x.push(depth as u32);
        y.push(hist[depth]);
    }
    (x, y)
}

fn select_depth_contigs(report: &QcReport, plot_max_contigs: usize) -> (Vec<String>, usize) {
    let mut totals: BTreeMap<String, u64> = BTreeMap::new();
    for sample in &report.samples {
        for (chromosome, histogram) in &sample.depth_distribution.per_chromosome_histograms {
            let score = histogram.iter().sum::<u64>();
            *totals.entry(chromosome.clone()).or_insert(0) += score;
        }
    }

    let mut ranked: Vec<(String, u64)> = totals.into_iter().collect();
    ranked.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let selected = ranked
        .iter()
        .take(plot_max_contigs)
        .map(|(chromosome, _)| chromosome.clone())
        .collect::<Vec<_>>();
    let omitted = ranked.len().saturating_sub(plot_max_contigs);
    (selected, omitted)
}

fn combine_other_contig_histogram(sample: &SampleQc, selected: &BTreeSet<String>) -> Vec<u64> {
    let mut combined = vec![0_u64; DEPTH_HIST_BINS];
    for (chromosome, histogram) in &sample.depth_distribution.per_chromosome_histograms {
        if selected.contains(chromosome) {
            continue;
        }
        for (i, value) in histogram.iter().enumerate().take(DEPTH_HIST_BINS) {
            combined[i] += value;
        }
    }
    combined
}

fn depth_scope_summary(report: &QcReport) -> &'static str {
    if report.settings.depth_scope == "reference_bases" {
        "includes depth 0"
    } else {
        "covered bases only"
    }
}

fn depth_y_axis_label(report: &QcReport) -> String {
    if report.settings.depth_scope == "reference_bases" {
        "Bases (including depth 0)".to_string()
    } else {
        "Covered bases".to_string()
    }
}

fn filter_json_chromosomes(report: &QcReport) -> QcReport {
    let mut filtered_samples = Vec::with_capacity(report.samples.len());
    for sample in &report.samples {
        let mut canonical: BTreeMap<String, Vec<u64>> = BTreeMap::new();
        for (chromosome, histogram) in &sample.depth_distribution.per_chromosome_histograms {
            let Some(canonical_name) = canonicalize_json_chromosome(chromosome) else {
                continue;
            };
            let entry = canonical
                .entry(canonical_name)
                .or_insert_with(|| vec![0_u64; DEPTH_HIST_BINS]);
            for (i, value) in histogram.iter().enumerate().take(DEPTH_HIST_BINS) {
                entry[i] += value;
            }
        }
        filtered_samples.push(SampleQc {
            sample_id: sample.sample_id.clone(),
            input_path: sample.input_path.clone(),
            counts: sample.counts.clone(),
            soft_clips: sample.soft_clips.clone(),
            mismatches: sample.mismatches.clone(),
            read_length: sample.read_length.clone(),
            depth_distribution: crate::model::DepthSummary {
                genome_wide_histogram: sample.depth_distribution.genome_wide_histogram.clone(),
                per_chromosome_histograms: canonical,
                overflow_bin_index: sample.depth_distribution.overflow_bin_index,
                eligible_bases_count: sample.depth_distribution.eligible_bases_count,
                depth_scope: sample.depth_distribution.depth_scope.clone(),
                reference_bases_total: sample.depth_distribution.reference_bases_total,
                covered_bases_total: sample.depth_distribution.covered_bases_total,
                zero_depth_bases_total: sample.depth_distribution.zero_depth_bases_total,
            },
            warnings: sample.warnings.clone(),
        });
    }
    QcReport {
        tool_version: report.tool_version.clone(),
        run_timestamp: report.run_timestamp,
        settings: report.settings.clone(),
        samples: filtered_samples,
    }
}

fn canonicalize_json_chromosome(name: &str) -> Option<String> {
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

#[cfg(test)]
mod tests {
    use super::{
        canonicalize_json_chromosome, depth_points, trace_is_visible, TraceView, ViewSelection,
    };

    #[test]
    fn depth_trace_visibility_matches_selected_chromosome() {
        let trace = TraceView::Depth {
            chromosome: "chr1".to_string(),
        };
        let selected_chr1 = ViewSelection::Depth {
            chromosome: "chr1".to_string(),
        };
        let selected_chr2 = ViewSelection::Depth {
            chromosome: "chr2".to_string(),
        };
        assert!(trace_is_visible(&trace, &selected_chr1));
        assert!(!trace_is_visible(&trace, &selected_chr2));
    }

    #[test]
    fn depth_points_preserves_zero_gaps() {
        let mut hist = vec![0_u64; 8];
        hist[1] = 5;
        hist[3] = 2;
        let (x, y) = depth_points(&hist);
        assert_eq!(x, vec![1, 2, 3]);
        assert_eq!(y, vec![5, 0, 2]);
    }

    #[test]
    fn canonical_chromosome_filter_allows_only_autosomes_x_y() {
        assert_eq!(
            canonicalize_json_chromosome("chr1"),
            Some("chr1".to_string())
        );
        assert_eq!(
            canonicalize_json_chromosome("22"),
            Some("chr22".to_string())
        );
        assert_eq!(canonicalize_json_chromosome("X"), Some("chrX".to_string()));
        assert_eq!(
            canonicalize_json_chromosome("chrY"),
            Some("chrY".to_string())
        );
        assert_eq!(canonicalize_json_chromosome("chrM"), None);
        assert_eq!(canonicalize_json_chromosome("GL000207.1"), None);
    }
}
