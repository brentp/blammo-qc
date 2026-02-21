# blammo-qc

`blammo-qc` is a Rust CLI for BAM/CRAM QC that writes:
- an optional JSON metrics report (only when `--output-json` is provided)
- a self-contained interactive Plotly HTML report

## Example usage

```bash
blammo-qc \
  --reference ref.fa \
  --output-html blammo.cohort.html \
  sample1.bam sample2.cram 
```
