import numpy as np
from Bio import SeqIO
import os

def one_hot_encode_sequence(sequence, valid_aas="ACDEFGHIKLMNPQRSTVWY"):
    aa_to_index = {aa: i for i, aa in enumerate(valid_aas)}
    encoded = np.zeros((len(sequence), len(valid_aas)))

    for i, aa in enumerate(sequence):
        if aa in aa_to_index:
            encoded[i, aa_to_index[aa]] = 1
    
    return encoded

def preprocess_protein_sequences(records):
    unwanted_chars = set('XJBOUZ-')
    
    print("Removing gaps and sequences containing unwanted characters (X, J, B, O, U, Z)...")
    cleaned_records = []
    for record in records:
        sequence = str(record.seq)
        if not any(char in sequence for char in unwanted_chars):
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            new_record = SeqRecord(
                Seq(sequence),
                id=record.id,
                description=record.description
            )
            cleaned_records.append(new_record)
    
    print(f"Number of sequences after removing unwanted characters: {len(cleaned_records)}")
    
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
    print(f"Reading sequences from {input_fasta}...")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Initial number of sequences: {len(records)}")
    
    filtered_records = preprocess_protein_sequences(records)
    
    print(f"Writing filtered sequences to {output_fasta}...")
    SeqIO.write(filtered_records, output_fasta, "fasta")
    print("Processing completed!")
    
    return filtered_records

if __name__ == "__main__":
    input_file = "xec_raw.fasta"
    output_file = "xec_preprocessed.fasta"
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found!")
    else:
        processed_records = main(input_file, output_file)
        
        print("Creating one-hot encodings...")
        encodings = [one_hot_encode_sequence(str(record.seq)) for record in processed_records]
        print(f"Created one-hot encodings for {len(encodings)} sequences")
