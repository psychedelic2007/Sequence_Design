
# Master-sequence RBD docking

Docked structural ensemble for the **designed master RBD** sequence (after scoring and optional `X` → consensus replacement in [`../../clean_master_sequence.py`](../../clean_master_sequence.py)).

## Files

- **`cluster1_1.pdb` … `cluster8_4.pdb`** — Thirty-two PDB structures: eight clusters with four models each.
- **`plots/`** — Same HTML metric plots as wild-type (interface RMSD, FCC, electrostatic/vdw/desolvation summaries).

## Comparison

Benchmark against [`../Wild_type/`](../Wild_type/) using matching plot types and cluster spreads. Large degradations in FCC or increased i-RMSD may flag positions worth revisiting in [`../../master_sequence_generation.py`](../../master_sequence_generation.py).

See also [`../README.md`](../README.md).
