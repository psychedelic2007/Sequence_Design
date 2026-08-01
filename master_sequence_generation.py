import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

mutation_df = pd.read_csv("mutation_freq_combined.csv") 
escape_df = pd.read_csv("escape_mutation.csv")
bcell_df = pd.read_csv("bepipred.csv")
data5 = pd.read_csv("t_cell/mhc1_averaged_scores.csv")
data6 = pd.read_csv("t_cell/mhc2_averaged_scores.csv")

mutation_values = mutation_df["mutation_frequency"].values 
escape_values = escape_df["Escape"].values  
bcell_values = bcell_df["Score"].values 
tcell_mhc1_values = data5["Average_Score"].values 
tcell_mhc2_values = data6["Average_Score"].values

tcell_values = 0.5 * tcell_mhc1_values + 0.5 * tcell_mhc2_values

print(f"Length of mutation: {len(mutation_values)}")
print(f"Length of escape: {len(escape_values)}")
print(f"Length of bcell: {len(bcell_values)}")
print(f"Length of tcell MHC-I: {len(tcell_mhc1_values)}")
print(f"Length of tcell MHC-II: {len(tcell_mhc2_values)}")
print(f"Length of tcell combined: {len(tcell_values)}")


# Min-max normalization
norm_mutation = (mutation_values - np.nanmin(mutation_values)) / (np.nanmax(mutation_values) - np.nanmin(mutation_values))
norm_escape = (escape_values - np.nanmin(escape_values)) / (np.nanmax(escape_values) - np.nanmin(escape_values))
norm_bcell = (bcell_values - np.nanmin(bcell_values)) / (np.nanmax(bcell_values) - np.nanmin(bcell_values))
norm_tcell = (tcell_values - np.nanmin(tcell_values)) / (np.nanmax(tcell_values) - np.nanmin(tcell_values))
norm_tcell_mhc1 = (tcell_mhc1_values - np.nanmin(tcell_mhc1_values)) / (np.nanmax(tcell_mhc1_values) - np.nanmin(tcell_mhc1_values))
norm_tcell_mhc2 = (tcell_mhc2_values - np.nanmin(tcell_mhc2_values)) / (np.nanmax(tcell_mhc2_values) - np.nanmin(tcell_mhc2_values))

norm_escape_full = np.full(223, np.nan)
norm_escape_full[12:213] = norm_escape

print(f"\nFull array lengths after padding:")
print(f"Mutation: {len(norm_mutation)}")
print(f"Escape (padded): {len(norm_escape_full)}")
print(f"B-cell: {len(norm_bcell)}")
print(f"T-cell combined: {len(norm_tcell)}")
print(f"T-cell MHC-I: {len(norm_tcell_mhc1)}")
print(f"T-cell MHC-II: {len(norm_tcell_mhc2)}")


position_labels = np.arange(319, 542)  # 319, 320, ..., 541
print(f"\nPosition labels: {len(position_labels)} positions (319-541)")

# Thresholds
mutation_threshold = 0.5
bcell_threshold = 0.7
tcell_threshold = 0.7
escape_threshold = 0.6

# For all 223 positions
mutation_flagged_mask = norm_mutation > mutation_threshold
bcell_flagged_mask = norm_bcell > bcell_threshold
tcell_flagged_mask = norm_tcell > tcell_threshold
tcell_mhc1_flagged_mask = norm_tcell_mhc1 > tcell_threshold
tcell_mhc2_flagged_mask = norm_tcell_mhc2 > tcell_threshold

# Get position numbers for flagged residues
mutation_flagged_positions = position_labels[mutation_flagged_mask]
bcell_flagged_positions = position_labels[bcell_flagged_mask]
tcell_flagged_positions = position_labels[tcell_flagged_mask]
tcell_mhc1_flagged_positions = position_labels[tcell_mhc1_flagged_mask]
tcell_mhc2_flagged_positions = position_labels[tcell_mhc2_flagged_mask]

# For escape, only consider positions 331-541 (indices 12:213)
escape_flagged_mask = norm_escape_full[12:213] > escape_threshold
escape_flagged_positions = position_labels[12:213][escape_flagged_mask]

print(f"\n=== THRESHOLD ANALYSIS ===")
print(f"Mutation probability threshold: {mutation_threshold}")
print(f"B-cell and T-cell threshold: {tcell_threshold}")
print(f"Escape threshold: {escape_threshold}")

print(f"\nPositions above thresholds:")
print(f"Mutation frequency > {mutation_threshold}: {len(mutation_flagged_positions)} positions")
print(f"B-cell epitope > {bcell_threshold}: {len(bcell_flagged_positions)} positions")
print(f"T-cell combined > {tcell_threshold}: {len(tcell_flagged_positions)} positions")
print(f"MHC-I > {tcell_threshold}: {len(tcell_mhc1_flagged_positions)} positions")
print(f"MHC-II > {tcell_threshold}: {len(tcell_mhc2_flagged_positions)} positions")
print(f"Escape > {escape_threshold}: {len(escape_flagged_positions)} positions")


# Hard constraint: Epitope-protected positions (B-cell OR T-cell >= 0.9) are NEVER modified
epitope_protected = bcell_flagged_mask | tcell_flagged_mask
escape_mask_full = np.full(223, False)
escape_mask_full[12:213] = norm_escape_full[12:213] > escape_threshold
modify_candidates = (mutation_flagged_mask | escape_mask_full) & ~epitope_protected
positions_to_modify = position_labels[modify_candidates]

print(f"\n=== MODIFICATION DECISION ===")
print(f"Epitope-protected positions (B-cell OR T-cell >= {bcell_threshold}): {np.sum(epitope_protected)}")
print(f"Positions flagged for modification: {len(positions_to_modify)}")
print(f"Positions retained: {223 - len(positions_to_modify)}")
print(f"\nPositions to modify (residue numbers): {sorted(positions_to_modify.tolist())}")


# WT RBD sequence (Wuhan-Hu-1)
wt_rbd = "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNFASFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
print(f"\nWT RBD length: {len(wt_rbd)}")

master_seq_list = list(wt_rbd)

for i, pos in enumerate(position_labels):
    if pos in positions_to_modify:
        master_seq_list[i] = 'X'  # Mark for modification
    # else: keep WT residue (already in master_seq_list)

master_sequence = ''.join(master_seq_list)

print(f"\nMaster Sequence:")
print(master_sequence)
print(f"\nModified positions (residue numbers): {sorted(positions_to_modify.tolist())}")

# Save master sequence
with open("master_sequence_template.txt", "w") as f:
    f.write(f"Master sequence (X = positions to modify):\n")
    f.write(master_sequence)
    f.write(f"\n\nPositions to modify: {sorted(positions_to_modify.tolist())}")
print("\n✓ Master sequence saved to 'master_sequence_template.txt'")

results_df = pd.DataFrame({
    'Position': position_labels,
    'Mutation_Norm': norm_mutation,
    'Escape_Norm': norm_escape_full,
    'Bcell_Norm': norm_bcell,
    'Tcell_Combined_Norm': norm_tcell,
    'Tcell_MHC1_Norm': norm_tcell_mhc1,
    'Tcell_MHC2_Norm': norm_tcell_mhc2,
    'Epitope_Protected': epitope_protected,
    'Flagged_For_Modify': modify_candidates,
    'WT_Residue': list(wt_rbd),
    'Master_Template': master_seq_list
})

results_df.to_csv("master_sequence_analysis.csv", index=False)
print("✓ Detailed analysis saved to 'master_sequence_analysis.csv'")
