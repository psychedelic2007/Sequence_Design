
# JN.1

SARS-CoV-2 **JN.1** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `jn1_`  
**Folder:** [`data/jn1/`](.) in the repository root.

## Raw input

`jn1_raw.fasta`

## Processing pipeline

1. Preprocess full spike → `jn1_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `jn1_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `jn1_preprocessed_rbd_cleaned.fasta`, `jn1_preprocessed_rbd_cleaned_align.fasta` when present
4. Mutational analysis → `jn1_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `jn1_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `jn1_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `jn1_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

Standard file set for this lineage; see table above.

See the [data index](../README.md) for cross-variant naming conventions.
