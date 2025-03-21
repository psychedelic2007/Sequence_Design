def replace_residues(sequence, positions, replacements):
    """
    Replace residues at specific positions in a sequence with new residues.
    
    Args:
        sequence (str): The original sequence
        positions (list): List of 1-indexed positions to replace
        replacements (str): String of replacement residues in order
    
    Returns:
        str: The modified sequence
    """
    # Convert sequence to list for easier manipulation
    seq_list = list(sequence)
    
    # Make sure positions and replacements have the same length
    if len(positions) != len(replacements):
        raise ValueError("Number of positions must match number of replacements")
    
    # Replace each residue (adjusting for 0-indexing)
    for i, pos in enumerate(positions):
        # Convert from 1-indexed to 0-indexed
        zero_indexed_pos = pos - 1
        
        # Check if position is valid
        if zero_indexed_pos < 0 or zero_indexed_pos >= len(sequence):
            raise ValueError(f"Position {pos} is out of bounds for sequence of length {len(sequence)}")
            
        # Replace the residue
        seq_list[zero_indexed_pos] = replacements[i]
    
    # Convert back to string
    return ''.join(seq_list)

# Your master sequence
master_sequence = "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNFASFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"

# Positions to replace (1-indexed) (add 318 to get original rbd position)
positions = [21, 28, 38, 50, 53, 55, 57, 58, 85, 87, 90, 99, 122, 126, 127, 128, 132, 134, 137, 138, 142, 159, 160, 163, 166, 168, 172, 175,
             178, 180, 183, 187]

# Replacement residues (These are characterstic mutations for different variants)
replacements = "DTTILPFAKNSNKTHSDWSLKNKKKVSRSRYH"

# Perform the replacements
modified_sequence = replace_residues(master_sequence, positions, replacements)

# Print results
print("Original sequence length:", len(master_sequence))
print("Modified sequence length:", len(modified_sequence))
print("\nModified sequence:")
print(modified_sequence)

# Verify the changes by showing original and new residues at each position
print("\nVerification of changes:")
print("Position | Original | New")
print("-" * 30)
for i, pos in enumerate(positions):
    original_residue = master_sequence[pos-1]
    new_residue = modified_sequence[pos-1]
    print(f"{pos:8d} | {original_residue:8s} | {new_residue}")
