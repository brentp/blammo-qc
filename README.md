# blammo-qc

`blammo-qc` is a Rust CLI for BAM/CRAM QC that writes:
- an optional JSON metrics report (only when `--output-json` is provided)
- a self-contained interactive Plotly HTML report

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

![Example image](https://private-user-images.githubusercontent.com/1739/553228223-02fe052e-d1c6-460e-96c0-d2aea1793aa5.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NzE3ODg0MTgsIm5iZiI6MTc3MTc4ODExOCwicGF0aCI6Ii8xNzM5LzU1MzIyODIyMy0wMmZlMDUyZS1kMWM2LTQ2MGUtOTZjMC1kMmFlYTE3OTNhYTUucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQVZDT0RZTFNBNTNQUUs0WkElMkYyMDI2MDIyMiUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyNjAyMjJUMTkyMTU4WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9MjRkODFhZTY0MTQ5YzExOTg0NGExZWZiNDFkYTZlNTY0MDk3ZWIyMDA5YzM3MmY0NDZmMDZkYWZkYTIyYWNjZiZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QifQ.xbfZRJ3R9E_pnCQGVxhHkabMnyPpJd2JZ-HhHtPmFG0)

![Example image 2](https://private-user-images.githubusercontent.com/1739/553228225-24a90162-9128-4d39-adb0-8543ffccdf14.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NzE3ODg0MTgsIm5iZiI6MTc3MTc4ODExOCwicGF0aCI6Ii8xNzM5LzU1MzIyODIyNS0yNGE5MDE2Mi05MTI4LTRkMzktYWRiMC04NTQzZmZjY2RmMTQucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQVZDT0RZTFNBNTNQUUs0WkElMkYyMDI2MDIyMiUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyNjAyMjJUMTkyMTU4WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9ZjNkYTk1NWRjNDRkYzVlMjI4NzhjZTM4MTU5Mzg3Zjg5MmU4YmMwYzRhZGQ0MmI4MTgzZjYwNDE5M2FjMGQzNiZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QifQ.NlbO7i-5VEhATjQbuuzv12lDc5UHIR2NT5_4_ziOwnQ)
