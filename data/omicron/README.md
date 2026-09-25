
# Omicron (B.1.1.529)

SARS-CoV-2 **Omicron (B.1.1.529)** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `omicron_`  
**Folder:** [`data/omicron/`](.) in the repository root.

## Raw input

`omicron_raw.fasta`

## Processing pipeline

1. Preprocess full spike → `omicron_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `omicron_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `omicron_preprocessed_rbd_cleaned.fasta`, `omicron_preprocessed_rbd_cleaned_align.fasta` when present
4. Mutational analysis → `omicron_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `omicron_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `omicron_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `omicron_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

Standard file set for this lineage; see table above.

See the [data index](../README.md) for cross-variant naming conventions.
