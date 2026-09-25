
# MHC class II inputs

Per-allele **NetMHCIIpan** (or compatible) prediction tables for **HLA class II** molecules (DR, DQ, DP loci).

## Files

Fifty-four CSV/TXT allele tables (e.g. `drb1_0401.csv`, `hla_dqa10102_dqb10604.csv`). Required columns for scoring: `MHC`, `Peptide`, `%Rank_EL`.

## Processing

Run [`../TCell_Calculation.py`](../TCell_Calculation.py) with `mhc_class="II"` to write per-allele residue scores under [`result/`](result/).

## Downstream

Averaged into [`../mhc2_averaged_scores.csv`](../mhc2_averaged_scores.csv) and [`../combined_mhc_scores.csv`](../combined_mhc_scores.csv).
