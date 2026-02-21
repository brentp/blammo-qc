use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};

#[derive(Debug, Parser)]
#[command(
    name = "blammo-qc",
    version,
    about = "Compute BAM/CRAM QC metrics and generate JSON + Plotly HTML reports."
)]
pub struct Cli {
    /// Input BAM/CRAM files (one or more).
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output JSON report path.
    #[arg(long)]
    pub output_json: PathBuf,

    /// Output HTML report path. Defaults to <output-json-stem>.html.
    #[arg(long)]
    pub output_html: Option<PathBuf>,

    /// Reference FASTA path. Required for CRAM input.
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Minimum base quality used in depth distribution.
    #[arg(long, default_value_t = 20)]
    pub min_base_quality: u8,

    /// Minimum mapping quality used in depth distribution.
    #[arg(long, default_value_t = 20)]
    pub min_mapping_quality: u8,

    /// Number of worker threads (defaults to logical CPUs).
    #[arg(long)]
    pub threads: Option<usize>,

    /// Optional sample name overrides, in the same order as inputs.
    #[arg(long = "sample-name")]
    pub sample_names: Vec<String>,

    /// Include duplicate-marked reads (default behavior excludes duplicates).
    #[arg(long)]
    pub include_duplicates: bool,

    /// Depth histogram basis: covered bases only, or full reference space (includes depth=0).
    #[arg(long, value_enum, default_value_t = DepthScope::CoveredBases)]
    pub depth_scope: DepthScope,

    /// Maximum number of per-contig depth views in HTML (remaining contigs are grouped as "Other").
    #[arg(long, default_value_t = 50)]
    pub plot_max_contigs: usize,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum DepthScope {
    CoveredBases,
    ReferenceBases,
}

impl DepthScope {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CoveredBases => "covered_bases",
            Self::ReferenceBases => "reference_bases",
        }
    }
}

#[derive(Debug, Clone)]
pub struct SampleInput {
    pub path: PathBuf,
    pub sample_id: String,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub inputs: Vec<SampleInput>,
    pub output_json: PathBuf,
    pub output_html: PathBuf,
    pub reference: Option<PathBuf>,
    pub min_base_quality: u8,
    pub min_mapping_quality: u8,
    pub threads: usize,
    pub include_duplicates: bool,
    pub depth_scope: DepthScope,
    pub plot_max_contigs: usize,
}

impl Cli {
    pub fn into_config(self) -> Result<Config> {
        if !self.sample_names.is_empty() && self.sample_names.len() != self.inputs.len() {
            bail!(
                "--sample-name count ({}) must match number of inputs ({})",
                self.sample_names.len(),
                self.inputs.len()
            );
        }

        let has_cram = self.inputs.iter().any(|path| is_cram_path(path));
        if has_cram && self.reference.is_none() {
            bail!("CRAM input detected; provide --reference <FASTA>");
        }

        if let Some(reference) = &self.reference {
            if !reference.exists() {
                bail!("reference FASTA not found: {}", reference.display());
            }
        }

        let threads = self.threads.unwrap_or_else(num_cpus::get).max(1);
        let output_html = self
            .output_html
            .unwrap_or_else(|| default_html_path(&self.output_json));

        if let Some(parent) = self.output_json.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent).with_context(|| {
                    format!(
                        "failed to create output JSON directory {}",
                        parent.display()
                    )
                })?;
            }
        }
        if let Some(parent) = output_html.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent).with_context(|| {
                    format!(
                        "failed to create output HTML directory {}",
                        parent.display()
                    )
                })?;
            }
        }

        let inputs = self
            .inputs
            .iter()
            .enumerate()
            .map(|(i, path)| SampleInput {
                path: path.clone(),
                sample_id: self
                    .sample_names
                    .get(i)
                    .cloned()
                    .unwrap_or_else(|| infer_sample_name(path)),
            })
            .collect();

        Ok(Config {
            inputs,
            output_json: self.output_json,
            output_html,
            reference: self.reference,
            min_base_quality: self.min_base_quality,
            min_mapping_quality: self.min_mapping_quality,
            threads,
            include_duplicates: self.include_duplicates,
            depth_scope: self.depth_scope,
            plot_max_contigs: self.plot_max_contigs,
        })
    }
}

pub fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .map(|ext| ext.eq_ignore_ascii_case("cram"))
        .unwrap_or(false)
}

fn default_html_path(output_json: &Path) -> PathBuf {
    let mut html = output_json.to_path_buf();
    html.set_extension("html");
    html
}

fn infer_sample_name(path: &Path) -> String {
    path.file_stem()
        .map(|stem| stem.to_string_lossy().to_string())
        .filter(|name| !name.is_empty())
        .unwrap_or_else(|| path.to_string_lossy().to_string())
}

#[cfg(test)]
mod tests {
    use super::default_html_path;

    #[test]
    fn derives_default_html_path() {
        let path = default_html_path(std::path::Path::new("out/report.json"));
        assert_eq!(path.to_string_lossy(), "out/report.html");
    }
}
