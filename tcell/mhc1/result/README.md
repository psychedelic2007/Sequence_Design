
# MHC-I per-allele results

Per-allele **residue-level T-cell scores** generated from [`../`](../) input CSVs by [`../../TCell_Calculation.py`](../../TCell_Calculation.py).

Each file is named like the input allele (e.g. `hla_a_0101.csv`) with columns:

- `Residue` — 1-based index along the 223-aa RBD
- `Score` — weighted aggregation of `1 / %Rank_EL` over overlapping peptides

Use these for per-allele inspection; [`../../mhc1_averaged_scores.csv`](../../mhc1_averaged_scores.csv) is the pooled input for master sequence generation.
