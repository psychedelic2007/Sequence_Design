import pandas as pd
import numpy as np
import os
import re

def get_mhc_weights(peptide_length, mhc_class):
    if mhc_class == "I":
        if peptide_length == 8:
            return [0.2, 1.0, 0.5, 0.3, 0.3, 0.5, 0.2, 1.0]  # MHC-I 8-mer
    elif mhc_class == "II":
        if peptide_length == 18:
            return [0.1, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 0.7, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1]  # MHC-II 18-mer
    raise ValueError(f"No weights for {mhc_class} {peptide_length}-mer")

def determine_mhc_class(mhc_name):
    if mhc_name.startswith(('HLA-DR', 'HLA-DQ', 'HLA-DP')):
        return "II"
    else:
        return "I" 

def extract_start_position(peptide, protein_sequence):
    match = protein_sequence.find(peptide)
    if match != -1:
        return match + 1
    return None

def calculate_residue_scores(peptides_df, protein_length, protein_sequence=None):
    residue_scores = np.zeros(protein_length)
    weight_sums = np.zeros(protein_length)  # To normalize based on weights
    
    skipped_count = 0
    processed_count = 0
    
    for _, row in peptides_df.iterrows():
        try:
            peptide = str(row['Peptide'])
            rank_el = float(row['%Rank_EL'])

            mhc_allele = str(row['MHC'])
            mhc_class = determine_mhc_class(mhc_allele)

            start_pos = None
            if 'start_pos' in peptides_df.columns:
                start_pos = int(row['start_pos'])
            elif protein_sequence is not None:
                start_pos = extract_start_position(peptide, protein_sequence)
            
            if start_pos is None:
                skipped_count += 1
                continue
                
            pep_len = len(peptide)

            if rank_el <= 0 or start_pos < 1 or (start_pos + pep_len - 1) > protein_length:
                skipped_count += 1
                continue

            try:
                weights = get_mhc_weights(pep_len, mhc_class)
            except ValueError:
                skipped_count += 1
                continue
                
            total_weight = sum(weights)
            score = 1 / rank_el  #Lower %Rank_EL=stronger binder

            for i, weight in enumerate(weights):
                residue_idx = start_pos + i - 1
                if residue_idx >= protein_length:
                    continue
                residue_scores[residue_idx] += (weight / total_weight) * score
                weight_sums[residue_idx] += weight
                
            processed_count += 1
            
        except (ValueError, KeyError, TypeError) as e:
            skipped_count += 1
            continue

    residue_scores = np.divide(residue_scores, weight_sums, where=weight_sums != 0)
    print(f"Processed {processed_count} peptides, skipped {skipped_count} peptides")
    return residue_scores

def process_mhc_data(mhc_class, input_dir, output_dir, protein_length, protein_sequence=None):
    os.makedirs(output_dir, exist_ok=True)
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    
    for csv_file in csv_files:
        try:
            file_path = os.path.join(input_dir, csv_file)
            peptides_df = pd.read_csv(file_path)

            required_cols = ['MHC', 'Peptide', '%Rank_EL']
            if not all(col in peptides_df.columns for col in required_cols):
                print(f"Skipping {csv_file}: missing required columns")
                continue

            if mhc_class == "I":
                peptides_df = peptides_df[peptides_df['MHC'].apply(determine_mhc_class) == "I"]
            elif mhc_class == "II":
                peptides_df = peptides_df[peptides_df['MHC'].apply(determine_mhc_class) == "II"]

            print(f"Processing {csv_file} for MHC class {mhc_class}...")
            tcell_scores = calculate_residue_scores(peptides_df, protein_length, protein_sequence)
            tcell_scores_df = pd.DataFrame({
                "Residue_Position": np.arange(1, protein_length + 1),
                "Tcell_Score": tcell_scores
            })

            output_filename = f"tcell_scores_{os.path.splitext(csv_file)[0]}_MHC{mhc_class}.csv"
            output_path = os.path.join(output_dir, output_filename)
            tcell_scores_df.to_csv(output_path, index=False)
            print(f"Processed {csv_file} and saved results to {output_path}")
            
        except Exception as e:
            print(f"Error processing {csv_file}: {str(e)}")

if __name__ == "__main__":
    rbd_length = 223
    protein_sequence = None

    mhc1_input_dir = "mhc1"
    mhc1_output_dir = "mhc1/results"
    
    mhc2_input_dir = "mhc2"
    mhc2_output_dir = "mhc2/results"

    print("Processing MHC-I data...")
    process_mhc_data("I", mhc1_input_dir, mhc1_output_dir, rbd_length, protein_sequence)

    print("Processing MHC-II data...")
    process_mhc_data("II", mhc2_input_dir, mhc2_output_dir, rbd_length, protein_sequence)
    
    print("All processing complete!")
