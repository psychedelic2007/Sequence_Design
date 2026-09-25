# Preprocessing

Scripts to prepare SARS-CoV-2 spike sequences before RBD extraction, pooling, and mutational analysis. Run from this directory or adjust paths in each script’s `__main__` block.

## Scripts

| File | Purpose |
|------|---------|
| [`preprocessing.py`](preprocessing.py) | Read a raw FASTA file; remove sequences containing gaps or ambiguous letters (`X`, `J`, `B`, `O`, `U`, `Z`); deduplicate; write `{variant}_preprocessed.fasta`. Optional one-hot encoding for ML. |
| [`extracting_rbd_sequence.py`](extracting_rbd_sequence.py) | Extract the RBD from full spike using flanking motifs **start** `RVQP` and **end** `CVNF`; write `{variant}_preprocessed_rbd.fasta`. |
| [`remove_empty_entries.py`](remove_empty_entries.py) | Optional: strip empty FASTA records (e.g. into `{variant}_preprocessed_rbd_cleaned.fasta`). |
| [`combine_preprocessed_file.py`](combine_preprocessed_file.py) | Concatenate multiple variant `*_preprocessed.fasta` files into one combined FASTA for cross-lineage analysis. |

## Recommended order

1. `preprocessing.py` on `{prefix}_raw.fasta` (or extracted sequences from `{prefix}_raw.zip`).
2. `extracting_rbd_sequence.py` on the preprocessed full-length file.
3. `remove_empty_entries.py` if needed, then external alignment → `*_preprocessed_rbd_cleaned_align.fasta` (produced outside this folder).
4. `combine_preprocessed_file.py` when building a multi-variant dataset for `mutational_analysis.py` at the repo root.

## Output naming

Outputs are written next to inputs under [`../data/<variant>/`](../data/). Prefix examples: `alpha_`, `ba1_`, `bq11_`, `kp23_`, `kp311_`, `xec_`.

## Dependencies

Biopython and NumPy (`preprocessing.py` only).
