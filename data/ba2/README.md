
# BA.2 (Omicron sublineage)

SARS-CoV-2 **BA.2 (Omicron sublineage)** spike sequences and RBD-focused mutation analysis for the Sequence_Design pipeline.

**File prefix:** `ba2_`  
**Folder:** [`data/ba2/`](.) in the repository root.

## Raw input

`ba2_raw.zip` (extract before preprocessing)

## Processing pipeline

1. Preprocess full spike → `ba2_preprocessed.fasta` ([`../../Preprocessing/preprocessing.py`](../../Preprocessing/preprocessing.py))
2. Extract RBD (223 aa) → `ba2_preprocessed_rbd.fasta` ([`../../Preprocessing/extracting_rbd_sequence.py`](../../Preprocessing/extracting_rbd_sequence.py))
3. Clean / align → `ba2_preprocessed_rbd_cleaned.fasta`, `ba2_preprocessed_rbd_cleaned_align.fasta` when present
4. Mutational analysis → `ba2_rbd_mutation_analysis_*` ([`../../mutational_analysis.py`](../../mutational_analysis.py))

## Analysis outputs

| File | Description |
|------|-------------|
| `ba2_rbd_mutation_analysis_mutation_stats.csv` | Position, reference AA, mutation frequency, conservation, statistics |
| `ba2_rbd_mutation_analysis_12feb2025_mutation_stats.csv` | Extended stats including `most_common_aa` (used by [`../../clean_master_sequence.py`](../../clean_master_sequence.py)) |
| `ba2_rbd_mutation_analysis_detailed_report.txt` | Text summary of top mutated positions |

## Notes

Standard file set for this lineage; see table above.

See the [data index](../README.md) for cross-variant naming conventions.
