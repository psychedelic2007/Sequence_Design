import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

bcell_threshold = 0.7      # Above 0.7 means strong B-cell epitope (retain)
tcell_threshold = 0.7      # Above 0.7 means strong T-cell epitope (retain)
mutation_threshold = 0.5   # Above 0.5 means high mutation frequency
escape_threshold = 0.6     # Above 0.6 means high escape probability

# (Preference order: escape > mutation > T-cell > B-cell)
w_mutation = 0.25
w_escape   = 0.35
w_tcell    = 0.2
w_bcell    = 0.2

# For positions before escape data is available, re-scale weights.
w_mutation_no_escape = 0.4
w_tcell_no_escape    = 0.3
w_bcell_no_escape    = 0.3

ref_seq = ("RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF")
master_seq_list = list(ref_seq)
n_total = len(master_seq_list)

offset = 331 - 319  # equals 12

mutation_df = pd.read_csv("path_to/rbd_mutation_analysis_12feb2025_mutation_stats.csv")
escape_df = pd.read_csv("path_to/escape_scores.csv")  # From Bloom et al
bcell_df = pd.read_csv("path_to/bcell_predictions.csv")  # From BEpipred
mhc1_df = pd.read_csv("t_cell/mhc1_averaged_scores.csv")
mhc2_df = pd.read_csv("t_cell/mhc2_averaged_scores.csv")

mutation_values = mutation_df["most_common_aa"].values 
escape_values = escape_df["escape"].values  
bcell_values = bcell_df["bcell_score"].values 
tcell_mhc1_values = mhc1_df["score"].values  
tcell_mhc2_values = mhc2_df["score"].values 

report_rows = []
modified_positions = []

for i in range(n_total):
    ref_residue = master_seq_list[i]
    mut_score = normalised_mutation_frq[i]
    bcell_score = normalised_bcell[i]
    tcell_score = normalised_tcell[i]

    if i < offset:
        esc_score = np.nan
        current_w_mutation = w_mutation_no_escape
        current_w_tcell = w_tcell_no_escape
        current_w_bcell = w_bcell_no_escape
        escape_above = False
    else:
        esc_score = normalised_escape[i - offset]
        current_w_mutation = w_mutation
        current_w_tcell = w_tcell
        current_w_bcell = w_bcell
        current_w_escape = w_escape
        escape_above = not np.isnan(esc_score)

    mutation_above = mut_score >= mutation_threshold
    bcell_above = bcell_score >= bcell_threshold
    tcell_above = tcell_score >= tcell_threshold

    if i < offset:  # Positions 319-330 (no escape data)
        weighted_score = (current_w_mutation * mut_score) - (current_w_tcell * tcell_score) - (current_w_bcell * bcell_score)
    else:  # Positions 331-541 (with escape data)
        weighted_score = (current_w_mutation * mut_score) + (current_w_escape * esc_score) - (current_w_tcell * tcell_score) - (current_w_bcell * bcell_score)

    if bcell_above or tcell_above:
        decision = "Retain"
        explanation = "High B-cell or T-cell epitope score; retained for immune recognition."
    elif mutation_above or escape_above:
        decision = "Modify"
        explanation = "Modification: High mutation and/or escape probability (and epitope scores not high)."
        master_seq_list[i] = "X"
        modified_positions.append(i)
    else:
        decision = "Retain"
        explanation = "None of the thresholds met; retained by default."

    report_rows.append({
        "Position": i + 319,
        "ReferenceResidue": ref_residue,
        "MutationFrequency": mut_score,
        "EscapeProbability": esc_score if i >= offset else None,
        "BcellEpitopeScore": bcell_score,
        "TcellEpitopeScore": tcell_score,
        "WeightedScore": weighted_score,
        "Decision": decision,
        "Explanation": explanation
    })

report_df = pd.DataFrame(report_rows)
report_df.to_csv("detailed_report.csv", index=False)

master_seq_str = "".join(master_seq_list)
print("Master Sequence:")
print(master_seq_str)
print("\nModified positions (residue numbers):", [pos + 319 for pos in modified_positions])

positions = np.arange(n_total) + 319
plt.figure(figsize=(12, 6))
plt.scatter(positions, normalised_mutation_frq, label="Mutation Frequency", color="blue", marker="o")
plt.scatter(positions[offset:], normalised_escape, label="Escape Probability", color="red", marker="s")
plt.scatter(positions, normalised_bcell, label="B-cell Epitope Score", color="magenta", marker="^")
plt.scatter(positions, normalised_tcell, label="T-cell Epitope Score", color="green", marker="d")

for pos in modified_positions:
    plt.axvline(x=pos + 319, color="gray", linestyle="--", alpha=0.5)

plt.xlabel("Residue Position (Actual Number)")
plt.ylabel("Normalized Score")
plt.title("Overlay of Mutation, Escape, B-cell, and T-cell Epitope Scores\n(Gray lines indicate modified positions)")
plt.legend()
plt.tight_layout()
plt.savefig("master_sequence.png", dpi=300)
plt.show()
