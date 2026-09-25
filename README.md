# Sequence_Design

Computational pipeline for designing a SARS-CoV-2 spike **receptor-binding domain (RBD)** master sequence by combining variant sequence data, mutational statistics, immune epitope scores (B-cell, T-cell / MHC), antibody escape probabilities, and structural docking results.

## Repository layout

| Path | Description |
|------|-------------|
| [`Preprocessing/`](Preprocessing/) | Clean full-length spike FASTA files, extract RBD, combine variants |
| [`data/`](data/) | Per-variant raw sequences, preprocessed FASTA, and RBD mutation analyses |
| [`tcell/`](tcell/) | MHC-I/II binding predictions and residue-level T-cell epitope scores |
| [`docking/`](docking/) | HADDOCK (or similar) cluster structures and interaction analysis plots |
| [`mutational_analysis.py`](mutational_analysis.py) | Position-wise mutation statistics vs. a reference RBD |
| [`master_sequence_generation.py`](master_sequence_generation.py) | Weighted scoring to mark positions for modification (`X`) or retention |
| [`clean_master_sequence.py`](clean_master_sequence.py) | Fill `X` placeholders using the most common amino acid per position |
| [`master_seq_verification.py`](master_seq_verification.py) | Apply lineage-specific residue replacements for validation |

## Typical workflow

1. **Preprocess** variant FASTA files (`Preprocessing/preprocessing.py`): remove gaps/ambiguous residues and duplicates.
2. **Extract RBD** (`Preprocessing/extracting_rbd_sequence.py`): slice residues using `RVQP` … `CVNF` patterns (223 aa RBD).
3. **Optional cleanup** (`Preprocessing/remove_empty_entries.py`): drop empty records; align cleaned RBD if needed.
4. **Combine variants** (`Preprocessing/combine_preprocessed_file.py`): merge lineage preprocessed files for pooled analysis.
5. **Mutational analysis** (`mutational_analysis.py`): compare aligned RBD sequences to Wuhan reference; export CSV, heatmaps, and reports (see each [`data/<variant>/`](data/) folder).
6. **T-cell scoring** (`tcell/TCell_Calculation.py`): convert NetMHCpan-style `%Rank_EL` tables into per-residue scores; use [`tcell/mhc1_averaged_scores.csv`](tcell/mhc1_averaged_scores.csv) and [`tcell/mhc2_averaged_scores.csv`](tcell/mhc2_averaged_scores.csv) in master sequence generation.
7. **Master sequence** (`master_sequence_generation.py`): integrate mutation, escape, B-cell, and T-cell inputs (update CSV paths in the script). Run [`clean_master_sequence.py`](clean_master_sequence.py) to resolve `X` to consensus residues.
8. **Docking** ([`docking/`](docking/)): compare wild-type vs. master-sequence RBD–ACE2 (or partner) structural ensembles.

## Reference RBD

Scripts assume a **223-residue** Wuhan-Hu-1 RBD segment (SPIKE positions **319–541**). The reference amino acid string used in `master_sequence_generation.py` begins with `RVQPTESI…`.

## Dependencies

Python 3 with:

- `numpy`, `pandas`, `matplotlib`, `seaborn`, `scipy`, `statsmodels`
- `Biopython` (`Bio`)

Install example:

```bash
pip install numpy pandas matplotlib seaborn scipy statsmodels biopython
```

## Data variants

Fourteen SARS-CoV-2 lineages are under [`data/`](data/): Alpha, Beta, Gamma, Delta, Omicron, BA.1, BA.2, BA.4, BA.5, BQ.1.1, JN.1, KP.2.3, KP.3.1.1, and XEC. Each subdirectory documents its files and naming prefix.

## Notes

- Several scripts use **hard-coded paths**; edit `input_file` / `output_file` variables before running on your machine.
- `master_sequence_generation.py` is a template: ensure normalization helpers and input CSV paths are configured for your run.
- Raw sequences for some lineages are shipped as `.zip` archives; extract before preprocessing.

## License

Add your license file here if not already present in the GitHub repository settings.
