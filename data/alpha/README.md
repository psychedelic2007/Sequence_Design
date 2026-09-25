
# Alpha (B.1.1.7)

SARS-CoV-2 **Alpha (B.1.1.7)** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `alpha_`  
**Folder:** [`data/alpha/`](.) in the repository root.

## Raw input

`alpha_raw.zip` (extract before preprocessing)

## Processing pipeline

1. Preprocess full spike → `alpha_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `alpha_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `alpha_preprocessed_rbd_cleaned.fasta`, `alpha_preprocessed_rbd_cleaned_align.fasta` when present
- `alpha_preprocessed_rbd_align.fasta` — aligned RBD (this variant uses `_rbd_align` instead of `_rbd_cleaned_align`)
4. Mutational analysis → `alpha_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `alpha_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `alpha_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `alpha_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

No `*_preprocessed_rbd_cleaned.fasta` in this folder; alignment file is named `alpha_preprocessed_rbd_align.fasta`.

See the [data index](../README.md) for cross-variant naming conventions.
