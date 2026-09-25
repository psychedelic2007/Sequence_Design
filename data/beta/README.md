
# Beta (B.1.351)

SARS-CoV-2 **Beta (B.1.351)** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `beta_`  
**Folder:** [`data/beta/`](.) in the repository root.

## Raw input

`beta_raw.fasta`

## Processing pipeline

1. Preprocess full spike → `beta_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `beta_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `beta_preprocessed_rbd_cleaned.fasta`, `beta_preprocessed_rbd_cleaned_align.fasta` when present
4. Mutational analysis → `beta_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `beta_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `beta_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `beta_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

Standard file set for this lineage; see table above.

See the [data index](../README.md) for cross-variant naming conventions.
