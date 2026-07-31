#This is an optional code. If your preprocessed fasta files
#contain gaps or empty sequences, use this code to remove them

from Bio import SeqIO

def remove_empty_entries(input_fasta, output_fasta):
    valid_records = []

    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.seq.strip():
                valid_records.append(record)

    with open(output_fasta, "w") as outfile:
        SeqIO.write(valid_records, outfile, "fasta")

    print(f"Removed empty sequences. Cleaned FASTA saved to {output_fasta}")

remove_empty_entries("xec/xec_preprocessed_rbd.fasta", "xec/xec_preprocessed_rbd_cleaned.fasta")
