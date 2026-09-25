
# Structural docking

Structural models and analysis for comparing **wild-type** and **master-sequence** RBD complexes (typically RBD–ACE2 or similar), organized by construct.

## Subdirectories

| Directory | Contents |
|-----------|----------|
| [`Wild_type/`](Wild_type/) | Docking clusters for the reference / wild-type RBD sequence |
| [`Master_sequence/`](Master_sequence/) | Docking clusters for the designed master RBD sequence |

Each construct folder contains:

- **`cluster{N}_{M}.pdb`** — representative structures from HADDOCK (or compatible) clustering (Wild_type: 5 clusters × 4 members; Master_sequence: 8 clusters × 4 members).
- **`plots/`** — Interactive HTML plots for interface metrics (e.g. i-RMSD, FCC, electrostatics, desolvation, H-bonds). Open in a web browser.
- **Screenshots** (Wild_type only) — static figures captured from the docking analysis UI.

## Using these results

Use cluster PDBs for visual inspection in PyMOL, ChimeraX, or VMD. Compare Wild_type vs. Master_sequence ensembles to assess whether designed changes preserve binding geometry and favorable interaction energy profiles suggested by the boxplot/scatter HTML files in each `plots/` folder.

## Related code

Sequence design and scoring live at the repository root ([`master_sequence_generation.py`](../master_sequence_generation.py), [`master_seq_verification.py`](../master_seq_verification.py)). This folder stores **outputs** only; docking jobs are run externally.
