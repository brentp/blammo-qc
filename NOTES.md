# Notes

## Read-length histogram behavior

- Includes primary mapped, non-duplicate reads only.
- Each included read contributes one count using `seq_len`.
- JSON histogram values are raw read counts.
- Bins are 1 bp up to 1000 bp.
- Bins are 50 bp above 1000 bp, keyed by lower bound (`1000`, `1050`, `1100`, ...).
- Internal accumulation is capped at 200,000 bp; lengths `>= 200000` go to the overflow bin.

