import pandas as pd

df = pd.read_csv("rbd_mutation_analysis_12feb2025_mutation_stats.csv")
most_common_residues = df["most_common_aa"].tolist()

# master sequence (with 'X' as placeholders)
master_seq = "RVXPTESIXXFPNITXLXPFXEVFXXXRFXSVXAXNXKXIXNXXAXXXVXYNXASFXXFKCXGVSPTXLNDXCFTNXXADSXXIRGXXVRQIXXGQTGKXADYXXXLXXDFXGXVIXWXXNNXXSXVXGXXXYXYXXFRKSNLKXXXRDISXXIYXAGXXPXXGXXXXNXXFXXXXXGXQPTXXVGYXXXRXXVLXFEXXXXXAXXXXXXXSTXLXXXXCXXX"

seq_list = list(master_seq)

# Replace 'X' with the corresponding residue from the CSV file.
for i in range(len(seq_list)):
    if seq_list[i] == "X":
        if i < len(most_common_residues):
            seq_list[i] = most_common_residues[i]
        else:
            pass

final_master_seq = "".join(seq_list)
print("Final Master Sequence:")
print(final_master_seq)
