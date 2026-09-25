
# T-cell epitope (MHC) scoring

Residue-level **T-cell epitope** scores derived from MHC-I and MHC-II peptide binding predictions (NetMHCpan-style `%Rank_EL` output).

## Main script

[`TCell_Calculation.py`](TCell_Calculation.py):

- Reads allele CSV files from [`mhc1/`](mhc1/) and [`mhc2/`](mhc2/).
- Maps peptides to RBD residue indices (223 residues) using position weights for 8-mer (MHC-I) and 18-mer (MHC-II) peptides.
- Writes per-allele outputs under `mhc1/result/` and `mhc2/result/` (configure paths in the script; default uses `mhc1/results` — align with your `result` folder name).

## Aggregated outputs (repository root of `tcell/`)

| File | Description |
|------|-------------|
| [`mhc1_averaged_scores.csv`](mhc1_averaged_scores.csv) | Mean MHC-I score per RBD residue |
| [`mhc2_averaged_scores.csv`](mhc2_averaged_scores.csv) | Mean MHC-II score per RBD residue |
| [`combined_mhc_scores.csv`](combined_mhc_scores.csv) | Normalized MHC-I, MHC-II, and combined score per position |

These CSVs feed [`../master_sequence_generation.py`](../master_sequence_generation.py) (high T-cell scores favor **retaining** residues for immune recognition).

## Input format

Each allele file must include columns: `MHC`, `Peptide`, `%Rank_EL`. Example: [`mhc1/hla_a_0101.csv`](mhc1/hla_a_0101.csv).

## Subfolders

- [`mhc1/`](mhc1/) — Class I alleles (HLA-A, HLA-B, HLA-C)
- [`mhc2/`](mhc2/) — Class II alleles (HLA-DR, HLA-DQ, HLA-DP)
