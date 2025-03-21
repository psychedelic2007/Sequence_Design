import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data1 = pd.read_csv("mutation_freq_combined.csv")
data2 = pd.read_csv("escape_mutation.csv")
data3 = pd.read_csv("bepipred.csv")
data4 = pd.read_csv("t_cell/combined_mhc_scores.csv")
data5 = pd.read_csv("t_cell/mhc1_averaged_scores.csv")
data6 = pd.read_csv("t_cell/mhc2_averaged_scores.csv")

normalised_mutation_frq = data1["mutation_frequency"]
normalised_escape = (data2["Escape"] - data2["Escape"].min()) / (data2["Escape"].max() - data2["Escape"].min())
normalised_bcell = (data3["Score"] - data3["Score"].min()) / (data3["Score"].max() - data3["Score"].min())
normalised_tcell = (data4["Combined_Score"] - data4["Combined_Score"].min())/ (data4["Combined_Score"].max() - data4["Combined_Score"].min())
normalised_tcell1 = (data5["Average_Score"] - data5["Average_Score"].min())/ (data5["Average_Score"].max() - data5["Average_Score"].min())
normalised_tcell2 = (data6["Average_Score"] - data6["Average_Score"].min())/ (data6["Average_Score"].max() - data6["Average_Score"].min())

print("Lenght of mutation is::", len(normalised_mutation_frq))
print("Lenght of escape is::", len(normalised_escape))
print("Lenght of bcell is::", len(normalised_bcell))
print("Lenght of tcell is::", len(normalised_tcell))

x = np.arange(201)  # Assuming 201 mutation positions
x1 = np.arange(223)
new_labels = np.arange(331, 532)  # Corresponding actual positions
labels = np.arange(319,541)

fig, axes = plt.subplots(4, 1, figsize=(10, 6), sharex=True)
axes[0].bar(x, normalised_mutation_frq[22:], color='blue', alpha=0.7)
axes[0].set_ylabel("Mutation Probability")
axes[0].set_title("Mutation Probability vs Position")
axes[1].bar(x, normalised_escape, color='red', alpha=0.7)
axes[1].set_ylabel("Escape Probability")
axes[1].set_title("Escape Probability vs Position")
axes[2].bar(x, normalised_bcell[22:], color="magenta", alpha=0.7)
axes[2].set_ylabel("BCell Epitope")
axes[2].set_title("BCell Epitope vs Position")
axes[3].bar(x, normalised_tcell[22:], color="green", alpha=0.7)
axes[3].set_ylabel("TCell Epitope")
axes[3].set_title("TCell Epitope vs Position")
axes[3].set_xticks(x[::10])
axes[3].set_xticklabels(new_labels[::10], rotation=45)
plt.xlabel("Position (RBD Region)")
plt.tight_layout()
plt.savefig("mutation_escape_bcell_tcell.png", dpi=300)
plt.show()

fig, axes = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
axes[0].bar(x1, normalised_tcell1, color='blue', alpha=0.7)
axes[0].set_ylabel("MHC-I Score")
axes[0].set_title("MHC-I vs Position")
axes[1].bar(x1, normalised_tcell2, color='red', alpha=0.7)
axes[1].set_ylabel("MHC-II Score")
axes[1].set_title("MHC-II vs Position")
axes[2].bar(x1, normalised_tcell, color="magenta", alpha=0.7)
axes[2].set_ylabel("TCell Epitope Score")
axes[2].set_title("TCell Epitope vs Position")
axes[2].set_xticks(x1[::10])
axes[2].set_xticklabels(labels[::10], rotation=45)
plt.xlabel("Position (RBD Region)")
plt.tight_layout()
plt.savefig("t_cell/tcell_combined.png", dpi=300)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- Define thresholds for decision-making ---
bcell_threshold = 0.7      # Above 0.7 means strong B-cell epitope (retain)
tcell_threshold = 0.7      # Above 0.7 means strong T-cell epitope (retain)
mutation_threshold = 0.5   # Above 0.5 means high mutation frequency
escape_threshold = 0.6     # Above 0.6 means high escape probability

# --- Define weight factors for positions with available escape data ---
# (Preference order: escape > mutation > T-cell > B-cell)
w_mutation = 0.25
w_escape   = 0.35
w_tcell    = 0.2
w_bcell    = 0.2

# For positions before escape data is available, re-scale weights.
w_mutation_no_escape = 0.4
w_tcell_no_escape    = 0.3
w_bcell_no_escape    = 0.3

# --- Prepare a reference sequence ---
ref_seq = ("RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF")

master_seq_list = list(ref_seq)
n_total = len(master_seq_list)

# Define the offset: escape data starts at residue 331 whereas RBD starts at 319.
offset = 331 - 319  # equals 12

# --- Generate dummy data (replace with actual data) ---
np.random.seed(42)
normalised_mutation_frq = np.random.rand(n_total)
normalised_escape = np.random.rand(n_total - offset)
normalised_bcell = np.random.rand(n_total)
normalised_tcell = np.random.rand(n_total)

# --- Create a detailed report list ---
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

    if i < offset:
        weighted_score = (current_w_mutation * mut_score) + (current_w_tcell * tcell_score) + (current_w_bcell * bcell_score)
    else:
        weighted_score = (current_w_mutation * mut_score) + (current_w_escape * esc_score) + (current_w_tcell * tcell_score) + (current_w_bcell * bcell_score)

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

# --- Visualization ---
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
