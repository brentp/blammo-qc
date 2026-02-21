# blammo-qc

`blammo-qc` is a Rust CLI for BAM/CRAM QC that writes:
- a JSON metrics report
- a self-contained interactive Plotly HTML report

## Example usage

```bash
cargo run -- \
  --output-json qc.json \
  --depth-scope reference-bases \
  --plot-max-contigs 30 \
  sample1.bam sample2.cram --reference ref.fa
```

