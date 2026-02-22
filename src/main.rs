mod cli;
mod metrics;
mod model;
mod report;

use anyhow::{Context, Result};
use chrono::Utc;
use clap::Parser;
use rayon::prelude::*;

use crate::{
    cli::Cli,
    metrics::{ProcessingOptions, process_sample},
    model::{QcReport, RunSettings},
    report::{write_html_report, write_json_report},
};

fn main() {
    init_logging();

    if let Err(err) = run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}

fn init_logging() {
    let _ = env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or("blammo_qc=info"),
    )
    .try_init();
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    let config = cli.into_config()?;

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build()
        .context("failed to build rayon thread pool")?;

    let processing_opts = ProcessingOptions {
        reference: config.reference.as_deref(),
        min_base_quality: config.min_base_quality,
        min_mapping_quality: config.min_mapping_quality,
        depth_scope: config.depth_scope,
        tag_bars: &config.tag_bars,
    };

    let sample_results = worker_pool.install(|| {
        config
            .inputs
            .par_iter()
            .map(|input| process_sample(input, &processing_opts))
            .collect::<Vec<_>>()
    });
    let mut samples = Vec::with_capacity(sample_results.len());
    for sample_result in sample_results {
        samples.push(sample_result?);
    }

    let report = QcReport {
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        run_timestamp: Utc::now(),
        settings: RunSettings {
            min_base_quality: config.min_base_quality,
            min_mapping_quality: config.min_mapping_quality,
            threads: config.threads,
            reference: config
                .reference
                .as_ref()
                .map(|path| path.to_string_lossy().to_string()),
            depth_scope: config.depth_scope.as_str().to_string(),
            plot_max_contigs: config.plot_max_contigs,
            tag_bars: config.tag_bars.clone(),
        },
        samples,
    };

    if let Some(output_json) = &config.output_json {
        write_json_report(&report, output_json)?;
    }
    write_html_report(&report, &config.output_html, config.plot_max_contigs)?;
    Ok(())
}
