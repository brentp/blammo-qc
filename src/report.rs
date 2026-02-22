use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result};
use plotly::{
    Bar, Configuration, Layout, Plot, Scatter,
    common::{Anchor, DashType, Font, Line, Marker, Mode, Orientation, Title, Visible},
    configuration::{
        DisplayModeBar, DoubleClick, ImageButtonFormats, ModeBarButtonName, ToImageButtonOptions,
    },
    layout::{
        Axis, AxisType, BarMode, DragMode, HoverMode, Legend, Margin, ModeBar, RangeMode,
        SpikeMode,
        update_menu::{Button, ButtonMethod, UpdateMenu, UpdateMenuType},
    },
};
use serde_json::json;

use crate::model::{BQ_BIN_LABELS, DEPTH_HIST_BINS, QcReport, SampleQc};

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
enum MetricTraceView {
    MappedReads,
    SoftClips,
    Unmapped,
    MismatchByBaseQuality,
    ReadLength,
    TagBar(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MetricSelection {
    MappedReads,
    SoftClips,
    Unmapped,
    MismatchByBaseQuality,
    ReadLength,
    TagBar(String),
}

pub fn write_html_report(report: &QcReport, path: &Path, plot_max_contigs: usize) -> Result<()> {
    let sample_ids: Vec<String> = report
        .samples
        .iter()
        .map(|sample| sample.sample_id.clone())
        .collect();

    let (mut selected_contigs, omitted_contig_count) =
        select_depth_contigs(report, plot_max_contigs);
    sort_depth_chromosomes(&mut selected_contigs);
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

    let mut depth_plot = Plot::new();
    let mut depth_trace_chromosomes: Vec<String> = Vec::new();
    let mut depth_x_max_by_chromosome: BTreeMap<String, u32> = BTreeMap::new();

    for chromosome in &depth_views {
        let mut combined_histogram = vec![0_u64; DEPTH_HIST_BINS];
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
            for (depth_bin, count) in histogram.iter().enumerate().take(DEPTH_HIST_BINS) {
                combined_histogram[depth_bin] += *count;
            }
            let (x, y) = depth_points(histogram);
            let label = if chromosome == "__genome__" {
                format!("{} depth (genome)", sample.sample_id)
            } else if chromosome == "__other__" {
                format!("{} depth (other contigs)", sample.sample_id)
            } else {
                format!("{} depth ({chromosome})", sample.sample_id)
            };
            depth_plot.add_trace(
                Scatter::new(x, y)
                    .mode(Mode::Lines)
                    .line(Line::new().width(2.0))
                    .hover_template("Depth: %{x}<br>Bases: %{y}<extra></extra>")
                    .name(label)
                    .visible(if chromosome == "__genome__" {
                        Visible::True
                    } else {
                        Visible::False
                    }),
            );
            depth_trace_chromosomes.push(chromosome.clone());
        }
        depth_x_max_by_chromosome.insert(
            chromosome.clone(),
            depth_x_max_from_histogram(&combined_histogram),
        );
    }

    let mut depth_options = vec![(
        "Genome-wide".to_string(),
        "__genome__".to_string(),
        format!(
            "Depth distribution (genome-wide, {})",
            depth_scope_summary(report)
        ),
        depth_x_max_by_chromosome
            .get("__genome__")
            .copied()
            .unwrap_or(1),
    )];
    for chromosome in depth_views
        .iter()
        .filter(|v| v.as_str() != "__genome__" && v.as_str() != "__other__")
    {
        depth_options.push((
            chromosome.clone(),
            chromosome.clone(),
            format!(
                "Depth distribution ({chromosome}, {})",
                depth_scope_summary(report)
            ),
            depth_x_max_by_chromosome
                .get(chromosome)
                .copied()
                .unwrap_or(1),
        ));
    }
    if omitted_contig_count > 0 {
        depth_options.push((
            format!("Other ({omitted_contig_count} contigs)"),
            "__other__".to_string(),
            format!(
                "Depth distribution (other contigs, {}, n={omitted_contig_count})",
                depth_scope_summary(report)
            ),
            depth_x_max_by_chromosome
                .get("__other__")
                .copied()
                .unwrap_or(1),
        ));
    }

    let depth_default_x_max = depth_x_max_by_chromosome
        .get("__genome__")
        .copied()
        .unwrap_or(1);
    let depth_y_axis = depth_y_axis_label(report);
    let depth_buttons: Vec<Button> = depth_options
        .iter()
        .map(|(label, selected_chromosome, title, x_max)| {
            let visible: Vec<bool> = depth_trace_chromosomes
                .iter()
                .map(|trace_chromosome| trace_chromosome == selected_chromosome)
                .collect();
            Button::new()
                .label(label.clone())
                .method(ButtonMethod::Update)
                .args(json!([
                    { "visible": visible },
                    {
                        "title": { "text": title },
                        "xaxis": {
                            "title": { "text": "Depth" },
                            "type": "linear",
                            "autorange": false,
                            "range": [0, x_max],
                            "rangemode": "nonnegative",
                            "fixedrange": false
                        },
                        "yaxis": { "title": { "text": depth_y_axis }, "fixedrange": true },
                        "dragmode": "zoom"
                    }
                ]))
        })
        .collect();

    let depth_menu = UpdateMenu::new()
        .ty(UpdateMenuType::Dropdown)
        .buttons(depth_buttons)
        .background_color("#FFFFFF")
        .border_color("#CBD5E1")
        .border_width(1)
        .font(Font::new().size(12).color("#0F172A"))
        .active(0)
        .x(0.01)
        .y(1.22);

    let depth_layout = Layout::new()
        .title(
            Title::with_text(format!(
                "Depth distribution (genome-wide, {})",
                depth_scope_summary(report)
            ))
            .font(
                Font::new()
                    .family("IBM Plex Sans, Arial, sans-serif")
                    .size(22)
                    .color("#0F172A"),
            ),
        )
        .height(680)
        .drag_mode(DragMode::Zoom)
        .font(
            Font::new()
                .family("IBM Plex Sans, Arial, sans-serif")
                .size(13)
                .color("#0F172A"),
        )
        .paper_background_color("#F8FAFC")
        .plot_background_color("#FFFFFF")
        .colorway(vec![
            "#0EA5E9", "#F97316", "#22C55E", "#A855F7", "#EF4444", "#14B8A6", "#F59E0B", "#3B82F6",
        ])
        .hover_mode(HoverMode::XUnified)
        .margin(Margin::new().top(130).right(30).bottom(70).left(70))
        .legend(
            Legend::new()
                .orientation(Orientation::Vertical)
                .x(0.01)
                .x_anchor(Anchor::Left)
                .y(0.99)
                .y_anchor(Anchor::Top),
        )
        .mode_bar(
            ModeBar::new()
                .background_color("#FFFFFF")
                .color("#334155")
                .active_color("#0EA5E9"),
        )
        .x_axis(
            Axis::new()
                .title("Depth")
                .type_(AxisType::Linear)
                .auto_range(false)
                .range(vec![0.0, depth_default_x_max as f64])
                .range_mode(RangeMode::NonNegative)
                .auto_margin(true)
                .line_color("#CBD5E1")
                .grid_color("#E2E8F0")
                .show_spikes(true)
                .spike_mode(SpikeMode::Across)
                .spike_color("#94A3B8")
                .spike_dash(DashType::Dot)
                .zero_line(false),
        )
        .y_axis(
            Axis::new()
                .title(depth_y_axis.clone())
                .fixed_range(true)
                .auto_margin(true)
                .line_color("#CBD5E1")
                .grid_color("#E2E8F0")
                .show_spikes(true)
                .spike_mode(SpikeMode::Across)
                .spike_color("#94A3B8")
                .spike_dash(DashType::Dot)
                .zero_line(false),
        )
        .update_menus(vec![depth_menu]);
    depth_plot.set_layout(depth_layout);
    depth_plot.set_configuration(
        Configuration::new()
            .responsive(true)
            .scroll_zoom(true)
            .display_logo(false)
            .display_mode_bar(DisplayModeBar::True)
            .double_click(DoubleClick::ResetAutoSize)
            .mode_bar_buttons_to_remove(vec![
                ModeBarButtonName::SendDataToCloud,
                ModeBarButtonName::Pan2d,
            ])
            .to_image_button_options(
                ToImageButtonOptions::new()
                    .format(ImageButtonFormats::Png)
                    .filename("blammo-depth-distribution")
                    .height(900)
                    .width(1500)
                    .scale(2),
            ),
    );

    let mut metrics_plot = Plot::new();
    let mut metric_trace_views: Vec<MetricTraceView> = Vec::new();

    let mapped_values: Vec<u64> = report
        .samples
        .iter()
        .map(|sample| sample.counts.primary_mapped_reads_used)
        .collect();
    metrics_plot.add_trace(
        Bar::new(sample_ids.clone(), mapped_values)
            .name("Mapped reads")
            .hover_template("Sample: %{x}<br>Mapped reads: %{y}<extra></extra>")
            .visible(Visible::False),
    );
    metric_trace_views.push(MetricTraceView::MappedReads);

    let soft_clip_values: Vec<u64> = report
        .samples
        .iter()
        .map(|sample| sample.soft_clips.total_soft_clipped_bases)
        .collect();
    metrics_plot.add_trace(
        Bar::new(sample_ids.clone(), soft_clip_values)
            .name("Soft-clipped bases")
            .hover_template("Sample: %{x}<br>Soft-clipped bases: %{y}<extra></extra>")
            .visible(Visible::False),
    );
    metric_trace_views.push(MetricTraceView::SoftClips);

    let unmapped_values: Vec<u64> = report
        .samples
        .iter()
        .map(|sample| sample.counts.unmapped_reads)
        .collect();
    metrics_plot.add_trace(
        Bar::new(sample_ids.clone(), unmapped_values)
            .name("Unmapped reads")
            .hover_template("Sample: %{x}<br>Unmapped reads: %{y}<extra></extra>")
            .visible(Visible::False),
    );
    metric_trace_views.push(MetricTraceView::Unmapped);

    for sample in &report.samples {
        let y: Vec<u64> = BQ_BIN_LABELS
            .iter()
            .map(|label| sample_count(sample, label))
            .collect();
        let x: Vec<&str> = BQ_BIN_LABELS.to_vec();
        metrics_plot.add_trace(
            Scatter::new(x, y)
                .mode(Mode::LinesMarkers)
                .line(Line::new().width(2.0))
                .marker(Marker::new().size(7))
                .hover_template("BQ bin: %{x}<br>Mismatches: %{y}<extra></extra>")
                .name(format!("{} mismatch BQ", sample.sample_id))
                .visible(Visible::False),
        );
        metric_trace_views.push(MetricTraceView::MismatchByBaseQuality);
    }

    for sample in &report.samples {
        let (x, y) = read_length_points(&sample.read_length.histogram);
        metrics_plot.add_trace(
            Scatter::new(x, y)
                .mode(Mode::Lines)
                .line(Line::new().width(2.0))
                .hover_template("Read length: %{x}<br>Reads: %{y:.2f}<extra></extra>")
                .name(format!("{} read length", sample.sample_id))
                .visible(Visible::True),
        );
        metric_trace_views.push(MetricTraceView::ReadLength);
    }
    let read_length_x_max = read_length_x_max_from_histograms(
        report
            .samples
            .iter()
            .map(|sample| &sample.read_length.histogram),
    ) as f64;
    let mut tag_metric_options: Vec<(
        String,
        MetricSelection,
        String,
        String,
        String,
        String,
        serde_json::Value,
        bool,
        bool,
        bool,
        Option<f64>,
    )> = Vec::new();
    for tag in &report.settings.tag_bars {
        let mut tag_values = BTreeSet::new();
        for sample in &report.samples {
            if let Some(counts) = sample.tag_value_counts.get(tag) {
                tag_values.extend(counts.keys().copied());
            }
        }
        let tag_values: Vec<i64> = tag_values.into_iter().collect();
        let x_values: Vec<String> = tag_values.iter().map(|value| value.to_string()).collect();
        for sample in &report.samples {
            let y_values: Vec<u64> = tag_values
                .iter()
                .map(|value| {
                    sample
                        .tag_value_counts
                        .get(tag)
                        .and_then(|counts| counts.get(value).copied())
                        .unwrap_or(0)
                })
                .collect();
            metrics_plot.add_trace(
                Bar::new(x_values.clone(), y_values)
                    .name(sample.sample_id.clone())
                    .hover_template(format!(
                        "Sample: {}<br>{tag} value: %{{x}}<br>Reads: %{{y}}<extra></extra>",
                        sample.sample_id
                    ))
                    .visible(Visible::False),
            );
            metric_trace_views.push(MetricTraceView::TagBar(tag.clone()));
        }
        tag_metric_options.push((
            format!("{tag} tag values"),
            MetricSelection::TagBar(tag.clone()),
            format!("{tag} tag integer values by sample"),
            format!("{tag} value"),
            "Reads".to_string(),
            "category".to_string(),
            json!(x_values),
            true,
            true,
            false,
            None,
        ));
    }

    let mut metric_options = vec![
        (
            "Mapped reads".to_string(),
            MetricSelection::MappedReads,
            "Mapped reads by sample".to_string(),
            "Sample".to_string(),
            "Mapped reads".to_string(),
            "category".to_string(),
            json!(sample_ids.clone()),
            true,
            true,
            false,
            None,
        ),
        (
            "Soft clips".to_string(),
            MetricSelection::SoftClips,
            "Soft-clipped bases by sample".to_string(),
            "Sample".to_string(),
            "Soft-clipped bases".to_string(),
            "category".to_string(),
            json!(sample_ids.clone()),
            true,
            true,
            false,
            None,
        ),
        (
            "Unmapped reads".to_string(),
            MetricSelection::Unmapped,
            "Unmapped reads by sample".to_string(),
            "Sample".to_string(),
            "Reads".to_string(),
            "category".to_string(),
            json!(sample_ids.clone()),
            true,
            true,
            false,
            None,
        ),
        (
            "Mismatches by base quality".to_string(),
            MetricSelection::MismatchByBaseQuality,
            "Mismatch counts partitioned by base quality".to_string(),
            "Base-quality bin".to_string(),
            "Mismatches".to_string(),
            "category".to_string(),
            json!(BQ_BIN_LABELS),
            false,
            false,
            false,
            None,
        ),
        (
            "Read length distribution".to_string(),
            MetricSelection::ReadLength,
            "Read-length distribution".to_string(),
            "Read length".to_string(),
            "Reads".to_string(),
            "linear".to_string(),
            json!(null),
            false,
            true,
            true,
            Some(read_length_x_max),
        ),
    ];
    metric_options.extend(tag_metric_options);
    let metric_buttons: Vec<Button> = metric_options
        .iter()
        .map(
            |(
                label,
                selection,
                title,
                x_label,
                y_label,
                x_axis_type,
                category_array,
                x_fixed,
                y_fixed,
                x_nonnegative,
                x_max,
            )| {
                let visible: Vec<bool> = metric_trace_views
                    .iter()
                    .map(|trace_view| metric_trace_is_visible(trace_view, selection))
                    .collect();
                let category_order = if category_array.is_null() {
                    serde_json::Value::Null
                } else {
                    json!("array")
                };
                let x_range_mode = if *x_nonnegative {
                    json!("nonnegative")
                } else {
                    serde_json::Value::Null
                };
                let x_autorange = x_max.is_none();
                let x_range = x_max
                    .map(|max| json!([0.0, max]))
                    .unwrap_or(serde_json::Value::Null);
                let drag_mode = if *x_fixed && *y_fixed { "pan" } else { "zoom" };
                Button::new()
                    .label(label.clone())
                    .method(ButtonMethod::Update)
                    .args(json!([
                        { "visible": visible },
                        {
                            "title": { "text": title },
                            "xaxis": {
                                "title": { "text": x_label },
                                "type": x_axis_type,
                                "categoryorder": category_order,
                                "categoryarray": category_array,
                                "autorange": x_autorange,
                                "range": x_range,
                                "rangemode": x_range_mode,
                                "fixedrange": x_fixed
                            },
                            "yaxis": {
                                "title": { "text": y_label },
                                "fixedrange": y_fixed
                            },
                            "dragmode": drag_mode
                        }
                    ]))
            },
        )
        .collect();

    let metrics_menu = UpdateMenu::new()
        .ty(UpdateMenuType::Dropdown)
        .buttons(metric_buttons)
        .background_color("#FFFFFF")
        .border_color("#CBD5E1")
        .border_width(1)
        .font(Font::new().size(12).color("#0F172A"))
        .active(4)
        .x(0.01)
        .y(1.22);

    let metrics_layout = Layout::new()
        .title(
            Title::with_text("Read-length distribution").font(
                Font::new()
                    .family("IBM Plex Sans, Arial, sans-serif")
                    .size(22)
                    .color("#0F172A"),
            ),
        )
        .height(680)
        .bar_mode(BarMode::Group)
        .drag_mode(DragMode::Zoom)
        .font(
            Font::new()
                .family("IBM Plex Sans, Arial, sans-serif")
                .size(13)
                .color("#0F172A"),
        )
        .paper_background_color("#F8FAFC")
        .plot_background_color("#FFFFFF")
        .colorway(vec![
            "#0EA5E9", "#F97316", "#22C55E", "#A855F7", "#EF4444", "#14B8A6", "#F59E0B", "#3B82F6",
        ])
        .hover_mode(HoverMode::XUnified)
        .margin(Margin::new().top(130).right(30).bottom(70).left(70))
        .legend(
            Legend::new()
                .orientation(Orientation::Vertical)
                .x(0.99)
                .x_anchor(Anchor::Right)
                .y(0.99)
                .y_anchor(Anchor::Top),
        )
        .mode_bar(
            ModeBar::new()
                .background_color("#FFFFFF")
                .color("#334155")
                .active_color("#0EA5E9"),
        )
        .x_axis(
            Axis::new()
                .title("Read length")
                .type_(AxisType::Linear)
                .auto_range(false)
                .range(vec![0.0, read_length_x_max])
                .range_mode(RangeMode::NonNegative)
                .fixed_range(false)
                .auto_margin(true)
                .line_color("#CBD5E1")
                .grid_color("#E2E8F0")
                .show_spikes(true)
                .spike_mode(SpikeMode::Across)
                .spike_color("#94A3B8")
                .spike_dash(DashType::Dot)
                .zero_line(false),
        )
        .y_axis(
            Axis::new()
                .title("Reads")
                .fixed_range(true)
                .auto_margin(true)
                .line_color("#CBD5E1")
                .grid_color("#E2E8F0")
                .show_spikes(true)
                .spike_mode(SpikeMode::Across)
                .spike_color("#94A3B8")
                .spike_dash(DashType::Dot)
                .zero_line(false),
        )
        .update_menus(vec![metrics_menu]);
    metrics_plot.set_layout(metrics_layout);
    metrics_plot.set_configuration(
        Configuration::new()
            .responsive(true)
            .scroll_zoom(false)
            .display_logo(false)
            .display_mode_bar(DisplayModeBar::True)
            .double_click(DoubleClick::ResetAutoSize)
            .mode_bar_buttons_to_remove(vec![
                ModeBarButtonName::SendDataToCloud,
                ModeBarButtonName::Zoom2d,
                ModeBarButtonName::Pan2d,
                ModeBarButtonName::Select2d,
                ModeBarButtonName::Lasso2d,
                ModeBarButtonName::ZoomIn2d,
                ModeBarButtonName::ZoomOut2d,
                ModeBarButtonName::AutoScale2d,
                ModeBarButtonName::ResetScale2d,
            ])
            .to_image_button_options(
                ToImageButtonOptions::new()
                    .format(ImageButtonFormats::Png)
                    .filename("blammo-qc-metrics")
                    .height(900)
                    .width(1500)
                    .scale(2),
            ),
    );

    let html = compose_two_plot_html(report, &depth_plot, &metrics_plot);
    let file = File::create(path)
        .with_context(|| format!("failed to create HTML output {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(html.as_bytes())
        .with_context(|| format!("failed to write HTML output {}", path.display()))?;
    writer
        .flush()
        .with_context(|| format!("failed to flush HTML output {}", path.display()))?;
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

fn metric_trace_is_visible(trace_view: &MetricTraceView, selection: &MetricSelection) -> bool {
    match (trace_view, selection) {
        (MetricTraceView::MappedReads, MetricSelection::MappedReads) => true,
        (MetricTraceView::SoftClips, MetricSelection::SoftClips) => true,
        (MetricTraceView::Unmapped, MetricSelection::Unmapped) => true,
        (MetricTraceView::MismatchByBaseQuality, MetricSelection::MismatchByBaseQuality) => true,
        (MetricTraceView::ReadLength, MetricSelection::ReadLength) => true,
        (MetricTraceView::TagBar(trace_tag), MetricSelection::TagBar(selected_tag)) => {
            trace_tag == selected_tag
        }
        _ => false,
    }
}

fn compose_two_plot_html(report: &QcReport, depth_plot: &Plot, metrics_plot: &Plot) -> String {
    let asset_bundle = extract_plotly_asset_bundle(&depth_plot.to_html()).unwrap_or_else(|| {
        r#"<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
<script src="https://cdn.plot.ly/plotly-2.12.1.min.js"></script>"#
            .to_string()
    });
    let depth_inline = depth_plot.to_inline_html(Some("depth-distribution-plot"));
    let metrics_inline = metrics_plot.to_inline_html(Some("metrics-plot"));
    let metrics_table = build_non_depth_metrics_table(report);
    let sample_count = report.samples.len();
    let tool_version = escape_html(&report.tool_version);
    let generated_at = report.run_timestamp.to_rfc3339();
    let depth_scope = if report.settings.depth_scope == "reference_bases" {
        "Reference bases (includes depth 0)"
    } else {
        "Covered bases only"
    };
    format!(
        r#"<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>Blammo QC Report</title>
    <style>
        :root {{
            --bg-top: #f4f7fb;
            --bg-bottom: #e8eef7;
            --ink: #0f172a;
            --muted: #475569;
            --card: rgba(255, 255, 255, 0.84);
            --stroke: rgba(148, 163, 184, 0.28);
            --accent-a: #0ea5e9;
            --accent-b: #1d4ed8;
            --shadow: 0 14px 40px rgba(15, 23, 42, 0.12);
            --radius-xl: 20px;
            --radius-md: 12px;
        }}
        * {{
            box-sizing: border-box;
        }}
        body {{
            margin: 0;
            color: var(--ink);
            font-family: "IBM Plex Sans", Arial, sans-serif;
            background:
                radial-gradient(900px 420px at 92% -6%, rgba(14, 165, 233, 0.18), transparent 60%),
                radial-gradient(820px 440px at -8% 10%, rgba(29, 78, 216, 0.14), transparent 58%),
                linear-gradient(180deg, var(--bg-top), var(--bg-bottom));
        }}
        .report-layout {{
            max-width: 1620px;
            margin: 0 auto;
            padding: 28px 24px 36px 24px;
        }}
        .hero {{
            margin-bottom: 24px;
            border-radius: var(--radius-xl);
            border: 1px solid var(--stroke);
            background:
                linear-gradient(135deg, rgba(255, 255, 255, 0.93), rgba(241, 245, 249, 0.8));
            box-shadow: var(--shadow);
            padding: 22px 22px 18px 22px;
        }}
        .hero h1 {{
            margin: 0;
            font-size: 28px;
            line-height: 1.15;
            letter-spacing: -0.02em;
        }}
        .hero h1 .title-link {{
            color: inherit;
            text-decoration: none;
            border-bottom: 2px solid rgba(14, 165, 233, 0.35);
        }}
        .hero h1 .title-link:hover {{
            border-bottom-color: rgba(29, 78, 216, 0.65);
        }}
        .hero p {{
            margin: 8px 0 0 0;
            color: var(--muted);
            font-size: 14px;
        }}
        .meta-grid {{
            margin-top: 16px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 10px;
        }}
        .meta-pill {{
            border: 1px solid rgba(148, 163, 184, 0.32);
            border-radius: 10px;
            background: rgba(248, 250, 252, 0.85);
            padding: 10px 12px;
        }}
        .meta-pill .label {{
            display: block;
            color: #64748b;
            font-size: 11px;
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }}
        .meta-pill .value {{
            display: block;
            margin-top: 4px;
            font-size: 14px;
            font-weight: 600;
            color: #0b1220;
        }}
        .plot-section {{
            margin-bottom: 24px;
            border-radius: var(--radius-xl);
            border: 1px solid var(--stroke);
            background: var(--card);
            box-shadow: var(--shadow);
            backdrop-filter: blur(6px);
            overflow: hidden;
        }}
        .section-head {{
            padding: 18px 22px 8px 22px;
        }}
        .section-head h2 {{
            margin: 0;
            font-size: 21px;
            line-height: 1.3;
            font-weight: 600;
        }}
        .section-head p {{
            margin: 6px 0 0 0;
            font-size: 13px;
            color: var(--muted);
        }}
        .section-body {{
            padding: 0 14px 14px 14px;
        }}
        .section-body .plotly-graph-div {{
            border-radius: var(--radius-md);
        }}
        .metrics-table-wrap {{
            margin-top: 20px;
            border: 1px solid #e2e8f0;
            border-radius: 10px;
            background: rgba(255, 255, 255, 0.96);
            overflow-x: auto;
        }}
        .metrics-table {{
            width: 100%;
            min-width: 1320px;
            border-collapse: collapse;
            font-size: 12px;
        }}
        .metrics-table thead th {{
            position: sticky;
            top: 0;
            z-index: 2;
        }}
        .metrics-table th,
        .metrics-table td {{
            padding: 9px 10px;
            border-bottom: 1px solid #e2e8f0;
            text-align: right;
            white-space: nowrap;
        }}
        .metrics-table th {{
            text-align: right;
            background: #eef2f7;
            font-weight: 600;
            color: #1e293b;
        }}
        .metrics-table th:first-child,
        .metrics-table td:first-child {{
            text-align: left;
            position: sticky;
            left: 0;
            background: #f8fafc;
            z-index: 1;
        }}
        .metrics-table th:first-child {{
            z-index: 3;
        }}
        .metrics-table tbody tr:hover {{
            background: #f8fafc;
        }}
        .metrics-table tbody tr:hover td:first-child {{
            background: #f1f5f9;
        }}
        @media (max-width: 760px) {{
            .report-layout {{
                padding: 16px 10px 24px 10px;
            }}
            .hero {{
                padding: 14px;
            }}
            .hero h1 {{
                font-size: 23px;
            }}
            .section-head {{
                padding: 14px 14px 8px 14px;
            }}
            .section-head h2 {{
                font-size: 18px;
            }}
            .section-body {{
                padding: 0 8px 10px 8px;
            }}
        }}
    </style>
    {asset_bundle}
</head>
<body>
    <main class="report-layout">
        <header class="hero">
            <h1><a class="title-link" href="https://github.com/brentp/blammo-qc" target="_blank" rel="noopener noreferrer">Blammo QC Report</a></h1>
            <p>Sequencing QC summary across {sample_count} sample(s).</p>
            <div class="meta-grid">
                <div class="meta-pill">
                    <span class="label">Tool Version</span>
                    <span class="value">{tool_version}</span>
                </div>
                <div class="meta-pill">
                    <span class="label">Generated</span>
                    <span class="value">{generated_at}</span>
                </div>
                <div class="meta-pill">
                    <span class="label">Depth Scope</span>
                    <span class="value">{depth_scope}</span>
                </div>
            </div>
        </header>
        <section class="plot-section">
            <div class="section-head">
                <h2>Depth Distribution</h2>
                <p>Defaults to genome-wide depth. Use the selector to switch chromosomes.</p>
            </div>
            <div class="section-body">
                {depth_inline}
            </div>
        </section>
        <section class="plot-section">
            <div class="section-head">
                <h2>Other QC Metrics</h2>
                <p>Defaults to mapped reads; switch other non-depth metrics from the dropdown.</p>
            </div>
            <div class="section-body">
                {metrics_table}
                {metrics_inline}
            </div>
        </section>
    </main>
</body>
</html>
"#
    )
}

fn extract_plotly_asset_bundle(full_html: &str) -> Option<String> {
    let start_marker =
        r#"<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>"#;
    let end_marker = r#"<div id="plotly-html-element""#;
    let start = full_html.find(start_marker)?;
    let end = full_html.find(end_marker)?;
    if start >= end {
        return None;
    }
    Some(full_html[start..end].trim().to_string())
}

fn build_non_depth_metrics_table(report: &QcReport) -> String {
    let mut html = String::new();
    html.push_str("<div class=\"metrics-table-wrap\"><table class=\"metrics-table\"><thead><tr>");
    html.push_str("<th>Sample</th>");
    html.push_str("<th>Median depth</th>");
    html.push_str("<th>Mode depth</th>");
    html.push_str("<th>Total reads seen</th>");
    html.push_str("<th>Mapped reads</th>");
    html.push_str("<th>Unmapped reads</th>");
    html.push_str("<th>Soft-clipped bases</th>");
    html.push_str("<th>Reads with soft clips</th>");
    html.push_str("<th>NM sum</th>");
    html.push_str("<th>Read length min</th>");
    html.push_str("<th>Read length p10</th>");
    html.push_str("<th>Read length median</th>");
    html.push_str("<th>Read length mean</th>");
    html.push_str("<th>Read length p90</th>");
    html.push_str("<th>Read length max</th>");
    html.push_str("<th>Warnings</th>");
    html.push_str("<th>BQ 0-9 mismatches</th>");
    html.push_str("<th>BQ 10-19 mismatches</th>");
    html.push_str("<th>BQ 20-29 mismatches</th>");
    html.push_str("<th>BQ 30-39 mismatches</th>");
    html.push_str("<th>BQ 40+ mismatches</th>");
    html.push_str("</tr></thead><tbody>");

    for sample in &report.samples {
        let median_depth =
            depth_median_from_histogram(&sample.depth_distribution.genome_wide_histogram);
        let mode_depth =
            depth_mode_from_histogram(&sample.depth_distribution.genome_wide_histogram);
        let total_reads = sample.counts.total_records_seen;
        html.push_str("<tr>");
        html.push_str(&format!("<td>{}</td>", escape_html(&sample.sample_id)));
        html.push_str(&format!("<td>{}</td>", median_depth));
        html.push_str(&format!("<td>{}</td>", mode_depth));
        html.push_str(&format!("<td>{}</td>", total_reads));
        html.push_str(&format!(
            "<td>{}</td>",
            format_count_with_percent(sample.counts.primary_mapped_reads_used, total_reads)
        ));
        html.push_str(&format!(
            "<td>{}</td>",
            format_count_with_percent(sample.counts.unmapped_reads, total_reads)
        ));
        html.push_str(&format!(
            "<td>{}</td>",
            sample.soft_clips.total_soft_clipped_bases
        ));
        html.push_str(&format!(
            "<td>{}</td>",
            format_count_with_percent(
                sample.soft_clips.reads_with_soft_clips,
                sample.counts.primary_mapped_reads_used
            )
        ));
        html.push_str(&format!("<td>{}</td>", sample.mismatches.nm_sum));
        html.push_str(&format!("<td>{}</td>", sample.read_length.min));
        html.push_str(&format!("<td>{}</td>", sample.read_length.p10));
        html.push_str(&format!("<td>{:.0}</td>", sample.read_length.median));
        html.push_str(&format!("<td>{:.0}</td>", sample.read_length.mean));
        html.push_str(&format!("<td>{}</td>", sample.read_length.p90));
        html.push_str(&format!("<td>{}</td>", sample.read_length.max));
        html.push_str(&format!("<td>{}</td>", sample.warnings.len()));
        html.push_str(&format!("<td>{}</td>", sample_count(sample, "0-9")));
        html.push_str(&format!("<td>{}</td>", sample_count(sample, "10-19")));
        html.push_str(&format!("<td>{}</td>", sample_count(sample, "20-29")));
        html.push_str(&format!("<td>{}</td>", sample_count(sample, "30-39")));
        html.push_str(&format!("<td>{}</td>", sample_count(sample, "40+")));
        html.push_str("</tr>");
    }

    html.push_str("</tbody></table></div>");
    html
}

fn escape_html(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(ch),
        }
    }
    escaped
}

fn format_count_with_percent(count: u64, total: u64) -> String {
    let percent = if total == 0 {
        0.0
    } else {
        (count as f64 * 100.0) / total as f64
    };
    format!("{count} ({percent:.2}%)")
}

fn depth_median_from_histogram(hist: &[u64]) -> u32 {
    let total: u64 = hist.iter().take(DEPTH_HIST_BINS).sum();
    if total == 0 {
        return 0;
    }
    let rank = (total - 1) / 2;
    let mut cumulative = 0_u64;
    for (depth, count) in hist.iter().take(DEPTH_HIST_BINS).enumerate() {
        cumulative += *count;
        if cumulative > rank {
            return depth as u32;
        }
    }
    0
}

fn depth_mode_from_histogram(hist: &[u64]) -> u32 {
    let mut mode_depth = 0_u32;
    let mut mode_count = 0_u64;
    for (depth, count) in hist.iter().take(DEPTH_HIST_BINS).enumerate() {
        if *count > mode_count {
            mode_count = *count;
            mode_depth = depth as u32;
        }
    }
    mode_depth
}

fn sort_depth_chromosomes(chromosomes: &mut [String]) {
    chromosomes.sort_by_key(|chromosome| chromosome_sort_key(chromosome));
}

fn chromosome_sort_key(name: &str) -> (u8, u16, String) {
    let trimmed = name.trim();
    let normalized = trimmed
        .strip_prefix("chr")
        .or_else(|| trimmed.strip_prefix("CHR"))
        .unwrap_or(trimmed);

    if let Ok(num) = normalized.parse::<u16>() {
        if (1..=22).contains(&num) {
            return (0, num, String::new());
        }
    }
    if normalized.eq_ignore_ascii_case("x") {
        return (0, 23, String::new());
    }
    if normalized.eq_ignore_ascii_case("y") {
        return (0, 24, String::new());
    }

    (1, u16::MAX, trimmed.to_ascii_lowercase())
}

fn read_length_points(hist: &std::collections::BTreeMap<u32, u64>) -> (Vec<u32>, Vec<f64>) {
    let mut x = Vec::with_capacity(hist.len());
    let mut y = Vec::with_capacity(hist.len());
    for (read_len, count) in hist {
        x.push(*read_len);
        if *read_len > 1_000 {
            y.push(*count as f64 / 50.0);
        } else {
            y.push(*count as f64);
        }
    }
    (x, y)
}

fn read_length_x_max_from_histograms<'a, I>(histograms: I) -> u32
where
    I: IntoIterator<Item = &'a BTreeMap<u32, u64>>,
{
    let mut max_with_at_least_ten = 0_u32;
    let mut max_nonzero = 0_u32;
    for histogram in histograms {
        for (read_len, count) in histogram {
            if *count > 0 {
                max_nonzero = max_nonzero.max(*read_len);
            }
            if *count >= 10 {
                max_with_at_least_ten = max_with_at_least_ten.max(*read_len);
            }
        }
    }
    if max_with_at_least_ten > 0 {
        max_with_at_least_ten
    } else if max_nonzero > 0 {
        max_nonzero
    } else {
        1
    }
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

fn depth_x_max_from_histogram(hist: &[u64]) -> u32 {
    if hist.is_empty() {
        return 1;
    }

    let upper = usize::min(hist.len(), DEPTH_HIST_BINS);
    let total_bases: u64 = hist.iter().take(upper).skip(1).sum();
    if total_bases == 0 {
        return 1;
    }
    let threshold = (total_bases + 1_999) / 2_000; // 0.05% of covered bases, rounded up.

    let Some(mut bin_index) = hist
        .iter()
        .take(upper)
        .enumerate()
        .skip(1)
        .rfind(|(_, count)| **count > 0)
        .map(|(depth, _)| depth)
    else {
        return 1;
    };

    while bin_index > 1 && hist[bin_index] < threshold {
        bin_index -= 1;
    }

    bin_index as u32
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
            tag_value_counts: sample.tag_value_counts.clone(),
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
    use std::collections::BTreeMap;

    use super::{
        MetricSelection, MetricTraceView, canonicalize_json_chromosome,
        depth_median_from_histogram, depth_mode_from_histogram, depth_points,
        depth_x_max_from_histogram, format_count_with_percent, metric_trace_is_visible,
        read_length_points, read_length_x_max_from_histograms, sort_depth_chromosomes,
    };

    #[test]
    fn metric_trace_visibility_matches_selected_metric() {
        assert!(metric_trace_is_visible(
            &MetricTraceView::MappedReads,
            &MetricSelection::MappedReads
        ));
        assert!(!metric_trace_is_visible(
            &MetricTraceView::MappedReads,
            &MetricSelection::Unmapped
        ));
        assert!(metric_trace_is_visible(
            &MetricTraceView::ReadLength,
            &MetricSelection::ReadLength
        ));
        assert!(metric_trace_is_visible(
            &MetricTraceView::TagBar("NM".to_string()),
            &MetricSelection::TagBar("NM".to_string())
        ));
        assert!(!metric_trace_is_visible(
            &MetricTraceView::TagBar("NM".to_string()),
            &MetricSelection::TagBar("AS".to_string())
        ));
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
    fn depth_x_max_keeps_bins_until_point_zero_five_percent_threshold() {
        let mut hist = vec![0_u64; 16];
        hist[1] = 3_000;
        hist[2] = 700;
        hist[3] = 200;
        hist[4] = 60;
        hist[5] = 30;
        hist[6] = 2;
        hist[7] = 1;
        assert_eq!(depth_x_max_from_histogram(&hist), 6);
    }

    #[test]
    fn read_length_points_normalize_wide_bins() {
        let mut hist = BTreeMap::new();
        hist.insert(999, 20);
        hist.insert(1000, 10);
        hist.insert(1050, 100);
        let (x, y) = read_length_points(&hist);
        assert_eq!(x, vec![999, 1000, 1050]);
        assert_eq!(y, vec![20.0, 10.0, 2.0]);
    }

    #[test]
    fn read_length_x_max_uses_last_bin_with_at_least_ten_reads() {
        let mut hist_a = BTreeMap::new();
        hist_a.insert(900, 10);
        hist_a.insert(1500, 9);
        let mut hist_b = BTreeMap::new();
        hist_b.insert(1200, 12);
        hist_b.insert(2000, 1);
        let x_max = read_length_x_max_from_histograms([&hist_a, &hist_b]);
        assert_eq!(x_max, 1200);
    }

    #[test]
    fn format_count_with_percent_uses_total_reads() {
        assert_eq!(format_count_with_percent(9923, 10_000), "9923 (99.23%)");
        assert_eq!(format_count_with_percent(0, 0), "0 (0.00%)");
    }

    #[test]
    fn depth_summary_stats_use_histogram() {
        let mut hist = vec![0_u64; 8];
        hist[0] = 1;
        hist[2] = 3;
        hist[5] = 2;
        assert_eq!(depth_mode_from_histogram(&hist), 2);
        assert_eq!(depth_median_from_histogram(&hist), 2);
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

    #[test]
    fn depth_chromosome_sort_is_human_readable() {
        let mut chromosomes = vec![
            "chr10".to_string(),
            "chr2".to_string(),
            "chrX".to_string(),
            "chr1".to_string(),
            "GL000192.1".to_string(),
            "chrY".to_string(),
        ];
        sort_depth_chromosomes(&mut chromosomes);
        assert_eq!(
            chromosomes,
            vec!["chr1", "chr2", "chr10", "chrX", "chrY", "GL000192.1"]
        );
    }
}
