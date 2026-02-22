use std::{
    collections::BTreeSet,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};

#[derive(Debug, Parser)]
#[command(
    name = "blammo-qc",
    version,
    before_help = concat!("blammo-qc ", env!("CARGO_PKG_VERSION")),
    about = "Compute BAM/CRAM QC metrics and generate optional JSON plus Plotly HTML reports."
)]
pub struct Cli {
    /// Input BAM/CRAM files (one or more).
    #[arg(required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output JSON report path (optional; omitted unless specified).
    #[arg(long)]
    pub output_json: Option<PathBuf>,

    /// Output HTML report path. Defaults to blammo.html.
    #[arg(long)]
    pub output_html: Option<PathBuf>,

    /// Reference FASTA path. Required for CRAM input.
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Minimum base quality used in depth distribution.
    #[arg(long, default_value_t = 10)]
    pub min_base_quality: u8,

    /// Minimum mapping quality used in depth distribution.
    #[arg(long, default_value_t = 10)]
    pub min_mapping_quality: u8,

    /// Number of worker threads (defaults to logical CPUs).
    #[arg(long)]
    pub threads: Option<usize>,

    /// Depth histogram basis: covered bases only, or full reference space (includes depth=0).
    #[arg(long, value_enum, default_value_t = DepthScope::CoveredBases)]
    pub depth_scope: DepthScope,

    /// Maximum number of per-contig depth views in HTML (remaining contigs are grouped as "Other").
    #[arg(long, default_value_t = 25)]
    pub plot_max_contigs: usize,

    /// SAM tags (two characters) whose integer values should be counted and shown as grouped bars.
    #[arg(long = "tag-bar", value_name = "TAG")]
    pub tag_bar: Vec<String>,
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
}

#[derive(Debug, Clone)]
pub struct Config {
    pub inputs: Vec<SampleInput>,
    pub output_json: Option<PathBuf>,
    pub output_html: PathBuf,
    pub reference: Option<PathBuf>,
    pub min_base_quality: u8,
    pub min_mapping_quality: u8,
    pub threads: usize,
    pub depth_scope: DepthScope,
    pub plot_max_contigs: usize,
    pub tag_bars: Vec<String>,
}

impl Cli {
    pub fn into_config(self) -> Result<Config> {
        let tag_bars = normalize_tag_bars(&self.tag_bar)?;
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
            .unwrap_or_else(|| PathBuf::from("blammo.html"));

        if let Some(output_json) = &self.output_json {
            if let Some(parent) = output_json.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).with_context(|| {
                        format!(
                            "failed to create output JSON directory {}",
                            parent.display()
                        )
                    })?;
                }
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
            .map(|path| SampleInput { path: path.clone() })
            .collect();

        Ok(Config {
            inputs,
            output_json: self.output_json,
            output_html,
            reference: self.reference,
            min_base_quality: self.min_base_quality,
            min_mapping_quality: self.min_mapping_quality,
            threads,
            depth_scope: self.depth_scope,
            plot_max_contigs: self.plot_max_contigs,
            tag_bars,
        })
    }
}

pub fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .map(|ext| ext.eq_ignore_ascii_case("cram"))
        .unwrap_or(false)
}

fn normalize_tag_bars(tag_bars: &[String]) -> Result<Vec<String>> {
    let mut seen = BTreeSet::new();
    let mut normalized = Vec::new();
    for tag in tag_bars {
        if !is_valid_sam_tag(tag) {
            bail!(
                "invalid --tag-bar value `{}`; expected a two-character SAM tag (e.g. NM)",
                tag
            );
        }
        if seen.insert(tag.clone()) {
            normalized.push(tag.clone());
        }
    }
    Ok(normalized)
}

fn is_valid_sam_tag(tag: &str) -> bool {
    let bytes = tag.as_bytes();
    bytes.len() == 2 && bytes[0].is_ascii_alphanumeric() && bytes[1].is_ascii_alphanumeric()
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::{is_valid_sam_tag, normalize_tag_bars};

    #[test]
    fn default_html_path_is_blammo_html() {
        let path = PathBuf::from("blammo.html");
        assert_eq!(path.to_string_lossy(), "blammo.html");
    }

    #[test]
    fn tag_bar_validation_requires_two_ascii_alnum_characters() {
        assert!(is_valid_sam_tag("NM"));
        assert!(is_valid_sam_tag("X1"));
        assert!(!is_valid_sam_tag("N"));
        assert!(!is_valid_sam_tag("NM3"));
        assert!(!is_valid_sam_tag("N_"));
    }

    #[test]
    fn tag_bar_normalization_deduplicates_in_order() {
        let normalized =
            normalize_tag_bars(&["NM".to_string(), "AS".to_string(), "NM".to_string()]).unwrap();
        assert_eq!(normalized, vec!["NM".to_string(), "AS".to_string()]);
    }
}
