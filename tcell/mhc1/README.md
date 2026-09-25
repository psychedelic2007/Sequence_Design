
# MHC class I inputs

Per-allele **NetMHCpan** (or compatible) prediction tables for **HLA class I** molecules on the Wuhan RBD sequence.

## Files

Seventeen CSV files, e.g. `hla_a_0101.csv`, `hla_b_0702.csv`, `hla_c_0602.csv`. Each row is a predicted peptide binder with `%Rank_EL` (lower rank = stronger predicted binding).

## Processing

Run [`../TCell_Calculation.py`](../TCell_Calculation.py) with `mhc_class="I"` to produce residue scores in [`result/`](result/) (one output CSV per input allele).

## Downstream

Per-residue scores are averaged into [`../mhc1_averaged_scores.csv`](../mhc1_averaged_scores.csv) for master sequence design.
