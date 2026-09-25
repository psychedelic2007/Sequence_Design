
# Gamma (P.1)

SARS-CoV-2 **Gamma (P.1)** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `gamma_`  
**Folder:** [`data/gamma/`](.) in the repository root.

## Raw input

`gamma_raw.fasta`

## Processing pipeline

1. Preprocess full spike → `gamma_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `gamma_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `gamma_preprocessed_rbd_cleaned.fasta`, `gamma_preprocessed_rbd_cleaned_align.fasta` when present
4. Mutational analysis → `gamma_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `gamma_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `gamma_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `gamma_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

Standard file set for this lineage; see table above.

See the [data index](../README.md) for cross-variant naming conventions.
