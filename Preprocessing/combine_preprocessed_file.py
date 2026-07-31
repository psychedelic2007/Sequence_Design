from Bio import SeqIO
import os

def count_sequences(fasta_file):
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

def combine_fasta_files(input_files, output_file):
    file_stats = {}
    combined_records = []
    
    print("\nProcessing files:")
    print("-" * 50)
    for fasta_file in input_files:
        if os.path.exists(fasta_file):
            sequence_count = count_sequences(fasta_file)
            file_stats[fasta_file] = sequence_count

            records = list(SeqIO.parse(fasta_file, "fasta"))
            combined_records.extend(records)
            
            print(f"{os.path.basename(fasta_file)}: {sequence_count:,} sequences")
        else:
            print(f"Warning: File not found - {fasta_file}")
    
    if combined_records:
        SeqIO.write(combined_records, output_file, "fasta")
        
        print("\nSummary Statistics:")
        print("-" * 50)
        for file_name, count in file_stats.items():
            print(f"{os.path.basename(file_name):<30} {count:>10,} sequences")
        print("-" * 50)
        print(f"{'Total combined sequences:':<30} {len(combined_records):>10,}")
        print(f"\nCombined sequences written to: {output_file}")
    else:
        print("No sequences were combined. Check input files.")

if __name__ == "__main__":
    input_files = [
        "new_results_13kseqs/variants/alpha/alpha_preprocessed.fasta",
        "new_results_13kseqs/variants/beta/beta_preprocessed.fasta",
        "new_results_13kseqs/variants/gamma/gamma_preprocessed.fasta",
        "new_results_13kseqs/variants/delta/delta_preprocessed.fasta",
        "new_results_13kseqs/variants/omicron/omicron_preprocessed.fasta",
        "new_results_13kseqs/variants/ba1/ba1_preprocessed.fasta",
        "new_results_13kseqs/variants/ba2/ba2_preprocessed.fasta",
        "new_results_13kseqs/variants/ba4/ba4_preprocessed.fasta",
        "new_results_13kseqs/variants/ba5/ba5_preprocessed.fasta",
        "new_results_13kseqs/variants/xec/xec_preprocessed.fasta",
        "new_results_13kseqs/variants/bq11/bq11_preprocessed.fasta",
        "new_results_13kseqs/variants/jn1/jn1_preprocessed.fasta",
        "new_results_13kseqs/variants/kp23/kp23_preprocessed.fasta",
        "new_results_13kseqs/variants/kp311/kp311_preprocessed.fasta"
    ]

    output_file = "new_results_13kseqs/variants/combined_variants.fasta"
    
    # Combine files and get statistics
    combine_fasta_files(input_files, output_file)
