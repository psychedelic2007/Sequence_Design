import numpy as np
from Bio import SeqIO
import os

def one_hot_encode_sequence(sequence, valid_aas="ACDEFGHIKLMNPQRSTVWY"):
    """
    One-hot encode a protein sequence, ignoring gaps and invalid characters.
    - sequence: Input sequence (string).
    - valid_aas: Valid amino acids to include (others are ignored).
    Returns: One-hot encoded matrix (L x 20).
    """
    # Create a mapping from amino acids to indices
    aa_to_index = {aa: i for i, aa in enumerate(valid_aas)}
    encoded = np.zeros((len(sequence), len(valid_aas)))
    
    # Fill matrix
    for i, aa in enumerate(sequence):
        if aa in aa_to_index:  # Only encode valid amino acids
            encoded[i, aa_to_index[aa]] = 1
    
    return encoded

def preprocess_protein_sequences(records):
    """Preprocess protein sequences by removing gaps, unwanted characters, and duplicates."""
    # Defining unwanted characters
    unwanted_chars = set('XJBOUZ-')
    
    print("Removing gaps and sequences containing unwanted characters (X, J, B, O, U, Z)...")
    cleaned_records = []
    for record in records:
        sequence = str(record.seq)
        # Checking if any unwanted character exists in the sequence
        if not any(char in sequence for char in unwanted_chars):
            # Creating a new SeqRecord with the cleaned sequence
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            new_record = SeqRecord(
                Seq(sequence),
                id=record.id,
                description=record.description
            )
            cleaned_records.append(new_record)
    
    print(f"Number of sequences after removing unwanted characters: {len(cleaned_records)}")
    
    # Removing duplicate sequences
    print("Removing duplicate sequences...")
    sequences = []
    unique_records = []
    for record in cleaned_records:
        sequence = str(record.seq)
        if sequence not in sequences:
            unique_records.append(record)
            sequences.append(sequence)
    print(f"Final number of sequences after removing duplicates: {len(unique_records)}")
    
    return unique_records

def main(input_fasta, output_fasta):
    """Main function to process the FASTA file."""
    print(f"Reading sequences from {input_fasta}...")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Initial number of sequences: {len(records)}")
    
    filtered_records = preprocess_protein_sequences(records)
    
    # Writing to new fasta file
    print(f"Writing filtered sequences to {output_fasta}...")
    SeqIO.write(filtered_records, output_fasta, "fasta")
    print("Processing completed!")
    
    return filtered_records

if __name__ == "__main__":
    input_file = "xec_raw.fasta"  # Replace with your input file path
    output_file = "xec_preprocessed.fasta"  # Replace with your desired output file path
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found!")
    else:
        processed_records = main(input_file, output_file)
        
        # Optionally, create one-hot encodings for the processed sequences
        print("Creating one-hot encodings...")
        encodings = [one_hot_encode_sequence(str(record.seq)) for record in processed_records]
        print(f"Created one-hot encodings for {len(encodings)} sequences")
