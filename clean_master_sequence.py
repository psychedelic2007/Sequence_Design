import pandas as pd

df = pd.read_csv("rbd_mutation_analysis_12feb2025_mutation_stats.csv")
most_common_residues = df["most_common_aa"].tolist()

# Provided master sequence (with 'X' as placeholders)
master_seq = "RVXPTESIXXFPNITXLXPFXEVFXXXRFXSVXAXNXKXIXNXXAXXXVXYNXASFXXFKCXGVSPTXLNDXCFTNXXADSXXIRGXXVRQIXXGQTGKXADYXXXLXXDFXGXVIXWXXNNXXSXVXGXXXYXYXXFRKSNLKXXXRDISXXIYXAGXXPXXGXXXXNXXFXXXXXGXQPTXXVGYXXXRXXVLXFEXXXXXAXXXXXXXSTXLXXXXCXXX"

# Convert the master sequence to a list to allow modifications.
seq_list = list(master_seq)

# Replace each 'X' with the corresponding residue from the CSV file.
for i in range(len(seq_list)):
    if seq_list[i] == "X":
        # Positions are 1-indexed in the CSV. Here, i is 0-indexed.
        # Replace with the most common residue for this position if available.
        if i < len(most_common_residues):
            seq_list[i] = most_common_residues[i]
        else:
            # If CSV doesn't have a value for this position, keep "X".
            pass

# Join the list back into the final sequence string.
final_master_seq = "".join(seq_list)
print("Final Master Sequence:")
print(final_master_seq)
