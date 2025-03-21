#This is an optional code. If your preprocessed fasta files
#contain gaps or empty sequences, use this code to remove them

from Bio import SeqIO

def remove_empty_entries(input_fasta, output_fasta):
    """Removes entries with no sequence from the given FASTA file."""
    # List to hold sequences that are not empty
    valid_records = []

    # Read the input FASTA file
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            # Check if sequence is not empty (not just whitespace)
            if record.seq.strip():  # Strip to remove whitespace
                valid_records.append(record)

    # Write the valid (non-empty) sequences to the output file
    with open(output_fasta, "w") as outfile:
        SeqIO.write(valid_records, outfile, "fasta")

    print(f"Removed empty sequences. Cleaned FASTA saved to {output_fasta}")

remove_empty_entries("xec/xec_preprocessed_rbd.fasta", "xec/xec_preprocessed_rbd_cleaned.fasta")
