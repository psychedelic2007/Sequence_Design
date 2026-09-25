
# Delta (B.1.617.2)

SARS-CoV-2 **Delta (B.1.617.2)** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `delta_`  
**Folder:** [`data/delta/`](.) in the repository root.

## Raw input

`delta_raw.fasta`

## Processing pipeline

1. Preprocess full spike → `delta_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `delta_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `delta_preprocessed_rbd_cleaned.fasta`, `delta_preprocessed_rbd_cleaned_align.fasta` when present
4. Mutational analysis → `delta_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `delta_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `delta_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `delta_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

Standard file set for this lineage; see table above.

See the [data index](../README.md) for cross-variant naming conventions.
