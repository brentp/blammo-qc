# blammo-qc

`blammo-qc` is a Rust CLI for BAM/CRAM QC that writes:
- an optional JSON metrics report (only when `--output-json` is provided)
- an optional downstream data JSON report (only when `--output-data-json` is provided)
- a self-contained interactive Plotly HTML report

It does a full depth pass of each bam/cram.

It is written entirely with ~~vibe coding~~ *agentic-engineering*, written 
with a combination of GPT-5.2 and GPT-5.4 via codex and Opus 4.5 via claude-code`

## Example usage

```bash
blammo-qc \
  --reference ref.fa \
  --tag-bar np \ # (make a bar plot grouped-by sample of the `np` integer tag)
  --tag-line ZX \ # (make a line plot for each sample of the `ZX` tag)
  --output-html blammo.cohort.html \
  sample1.bam sample2.cram 
```

## Example output

![Example output 1](images/blammo1.png)

![Example output 2](images/blammo2.png)

![Example output 3](images/blammo3.png)

<details>
<summary>Downstream data JSON shape</summary>

`--output-data-json` writes a structured JSON file for downstream tools. It includes extracted and derived report data, but not Plotly traces, Plotly layout, HTML, CSS, JavaScript, or browser assets.

Depth histograms are written only through the highest nonzero depth bin for that sample and view. For example, if the highest observed depth is `2`, the histogram is `[depth_0_count, depth_1_count, depth_2_count]`, not the full internal 32,768-bin array.

Tag metrics keep a global `values` list for ordering and emit sparse per-sample `counts` maps containing only nonzero counts.

```json
{
  "tool_version": "0.1.9",
  "run_timestamp": "2026-04-24T12:34:56Z",
  "settings": {
    "min_base_quality": 10,
    "min_mapping_quality": 10,
    "threads": 1,
    "reference": null,
    "depth_scope": "covered_bases",
    "plot_max_contigs": 25,
    "tag_bars": ["NM"],
    "tag_lines": ["ZX"]
  },
  "data": {
    "sample_count": 1,
    "samples": [
      {
        "sample_id": "sample_from_rg",
        "input_path": "/path/to/fixture.bam"
      }
    ],
    "depth_views": [
      {
        "key": "genome_wide",
        "label": "Genome-wide",
        "kind": "genome",
        "x_max": 2
      },
      {
        "key": "chr1",
        "label": "chr1",
        "kind": "contig",
        "x_max": 2
      }
    ],
    "sample_depths": [
      {
        "sample_id": "sample_from_rg",
        "histograms": {
          "genome_wide": [0, 1, 4],
          "chr1": [0, 1, 4]
        }
      }
    ],
    "metrics_table": [
      {
        "sample_id": "sample_from_rg",
        "median_depth": 2,
        "mode_depth": 2,
        "total_records_seen": 4,
        "primary_mapped_reads_used": 2,
        "primary_mapped_reads_percent": 50.0,
        "unmapped_reads": 1,
        "unmapped_reads_percent": 25.0,
        "total_soft_clipped_bases": 1,
        "soft_clipped_bases_per_million_bases": 100000.0,
        "reads_with_soft_clips": 1,
        "reads_with_soft_clips_percent": 50.0,
        "nm_sum": 1,
        "nm_per_million_bases": 100000.0,
        "read_length_min": 5,
        "read_length_p10": 5,
        "read_length_median": 5.0,
        "read_length_mean": 5.0,
        "read_length_p90": 5,
        "read_length_max": 5,
        "warnings_count": 0,
        "bq_mismatch_bins": {
          "0-9": 0,
          "10-19": 0,
          "20-29": 0,
          "30-39": 0,
          "40+": 1
        }
      }
    ],
    "tag_metrics": [
      {
        "tag": "NM",
        "view": "bar",
        "values": [0, 1],
        "sample_counts": [
          {
            "sample_id": "sample_from_rg",
            "counts": {
              "0": 1,
              "1": 1
            }
          }
        ]
      },
      {
        "tag": "ZX",
        "view": "line",
        "values": [7],
        "sample_counts": [
          {
            "sample_id": "sample_from_rg",
            "counts": {
              "7": 2
            }
          }
        ]
      }
    ]
  }
}
```

</details>
