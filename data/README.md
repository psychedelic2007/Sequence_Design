# Variant sequence data

This directory holds **one folder per SARS-CoV-2 lineage** used in the Sequence_Design pipeline. Each variant follows the same general processing stages; file names use a lineage-specific prefix (see each subfolder’s README).

## Lineage folders

| Folder | Lineage (common name) | File prefix |
|--------|------------------------|-------------|
| [`alpha/`](alpha/) | Alpha (B.1.1.7) | `alpha_` |
| [`beta/`](beta/) | Beta (B.1.351) | `beta_` |
| [`gamma/`](gamma/) | Gamma (P.1) | `gamma_` |
| [`delta/`](delta/) | Delta (B.1.617.2) | `delta_` |
| [`omicron/`](omicron/) | Omicron (B.1.1.529) | `omicron_` |
| [`ba1/`](ba1/) | BA.1 | `ba1_` |
| [`ba2/`](ba2/) | BA.2 | `ba2_` |
| [`ba4/`](ba4/) | BA.4 | `ba4_` |
| [`ba5/`](ba5/) | BA.5 | `ba5_` |
| [`bq1.1/`](bq1.1/) | BQ.1.1 | `bq11_` |
| [`jn1/`](jn1/) | JN.1 | `jn1_` |
| [`kp2.3/`](kp2.3/) | KP.2.3 | `kp23_` |
| [`kp.3.1.1/`](kp.3.1.1/) | KP.3.1.1 | `kp311_` |
| [`xec/`](xec/) | XEC | `xec_` |

## Typical files in each variant folder

| Pattern | Description |
|---------|-------------|
| `{prefix}_raw.fasta` or `{prefix}_raw.zip` | Unprocessed sequences from surveillance / GISAID-style dumps |
| `{prefix}_preprocessed.fasta` | Full spike after cleaning and deduplication |
| `{prefix}_preprocessed_rbd.fasta` | RBD-only FASTA (223 aa) |
| `{prefix}_preprocessed_rbd_cleaned.fasta` | RBD with empty entries removed |
| `{prefix}_preprocessed_rbd_cleaned_align.fasta` | Aligned RBD for comparative analysis (when present) |
| `{prefix}_rbd_mutation_analysis_mutation_stats.csv` | Per-position mutation statistics vs. reference |
| `{prefix}_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Alternate stats export (includes `most_common_aa` for master sequence filling) |
| `{prefix}_rbd_mutation_analysis_detailed_report.txt` | Human-readable mutation summary |

Not every variant includes every intermediate file (for example, Alpha uses `alpha_preprocessed_rbd_align.fasta`; BA.4 may omit a cleaned-align file).

## How these files are produced

1. Preprocessing and RBD extraction: [`../Preprocessing/`](../Preprocessing/)
2. Mutation analysis: [`../mutational_analysis.py`](../mutational_analysis.py) with variant RBD FASTA and Wuhan reference

Pooled analysis can combine preprocessed full-length files via [`../Preprocessing/combine_preprocessed_file.py`](../Preprocessing/combine_preprocessed_file.py).
