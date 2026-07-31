def replace_residues(sequence, positions, replacements):
    seq_list = list(sequence)
    if len(positions) != len(replacements):
        raise ValueError("Number of positions must match number of replacements")

    for i, pos in enumerate(positions):
        zero_indexed_pos = pos - 1
        if zero_indexed_pos < 0 or zero_indexed_pos >= len(sequence):
            raise ValueError(f"Position {pos} is out of bounds for sequence of length {len(sequence)}")

        seq_list[zero_indexed_pos] = replacements[i]
    return ''.join(seq_list)

# master sequence
master_sequence = "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNFASFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"

# Positions to replace (1-indexed) (add 318 to get original rbd position)
positions = [21, 28, 38, 50, 53, 55, 57, 58, 85, 87, 90, 99, 122, 126, 127, 128, 132, 134, 137, 138, 142, 159, 160, 163, 166, 168, 172, 175,
             178, 180, 183, 187]

# Replacement residues (These are characterstic mutations for different variants)
replacements = "DTTILPFAKNSNKTHSDWSLKNKKKVSRSRYH"

modified_sequence = replace_residues(master_sequence, positions, replacements)

print("Original sequence length:", len(master_sequence))
print("Modified sequence length:", len(modified_sequence))
print("\nModified sequence:")
print(modified_sequence)

print("\nVerification of changes:")
print("Position | Original | New")
print("-" * 30)
for i, pos in enumerate(positions):
    original_residue = master_sequence[pos-1]
    new_residue = modified_sequence[pos-1]
    print(f"{pos:8d} | {original_residue:8s} | {new_residue}")
