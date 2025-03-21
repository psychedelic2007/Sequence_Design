from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.stats.multitest as multi

def analyze_mutations(records, reference_sequence):
    """
    Enhanced mutation analysis with conservation scores and detailed mutation tracking
    """
    position_frequencies = defaultdict(lambda: defaultdict(int))
    mutation_data = defaultdict(list)
    total_sequences = len(records)
    
    # Track all mutations and their frequencies
    for record in records:
        sequence = str(record.seq).replace('-', '')
        for pos in range(min(len(sequence), len(reference_sequence))):
            aa = sequence[pos]
            ref_aa = reference_sequence[pos]
            position_frequencies[pos][aa] += 1
            
            if aa != ref_aa:
                mutation_data[pos].append({
                    'original': ref_aa,
                    'mutation': aa,
                    'sequence_id': record.id
                })
    
    # Calculate conservation scores
    conservation_scores = {}
    for pos in position_frequencies:
        total = sum(position_frequencies[pos].values())
        # Shannon entropy calculation
        entropy = 0
        for aa_count in position_frequencies[pos].values():
            if aa_count > 0:
                freq = aa_count / total
                entropy -= freq * np.log2(freq)
        conservation_scores[pos] = 1 - (entropy / np.log2(20))  # Normalize by max possible entropy
    
    return position_frequencies, mutation_data, conservation_scores

def calculate_mutation_statistics(position_frequencies, mutation_data, total_sequences, conservation_scores, reference_sequence):
    """
    Enhanced statistical analysis including significance testing
    """
    mutation_stats = []
    
    for pos in sorted(position_frequencies.keys()):
        freq_dict = position_frequencies[pos]
        mutations = mutation_data.get(pos, [])
        total_mutations = len(mutations)
        
        # Count frequency of each mutation
        mutation_counts = defaultdict(int)
        for mut in mutations:
            mutation_counts[mut['mutation']] += 1
            
        # Sort mutations by frequency
        sorted_mutations = sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True)
        mutation_string = ', '.join(f"{aa}({count})" for aa, count in sorted_mutations)
        
        # Calculate statistical significance
        # Chi-square test against uniform distribution
        if total_mutations > 0:
            observed = list(mutation_counts.values())
            expected = [total_mutations/len(mutation_counts)] * len(mutation_counts)
            chi2, p_value = stats.chisquare(observed, expected)
        else:
            chi2, p_value = 0, 1.0
            
        stats_dict = {
            'position': pos + 1,
            'reference_aa': reference_sequence[pos],
            'total_mutations': total_mutations,
            'mutation_frequency': total_mutations / total_sequences,
            'mutations': mutation_string,
            'conservation_score': conservation_scores[pos],
            'chi_square': chi2,
            'p_value': p_value
        }
        mutation_stats.append(stats_dict)
    
    # Convert to DataFrame and adjust p-values for multiple testing
    stats_df = pd.DataFrame(mutation_stats)
    stats_df['adjusted_p_value'] = multi.multipletests(stats_df['p_value'], method='fdr_bh')[1]
    
    return stats_df

def plot_enhanced_heatmap(position_frequencies, reference_sequence, output_file):
    """
    Create an enhanced heatmap with better visualization
    """
    all_aas = "ACDEFGHIKLMNPQRSTVWY"
    positions = sorted(position_frequencies.keys())
    
    # Create frequency matrix
    matrix_data = []
    for pos in positions:
        total = sum(position_frequencies[pos].values())
        row = [position_frequencies[pos].get(aa, 0) / total * 100 for aa in all_aas]
        matrix_data.append(row)
    
    matrix = np.array(matrix_data)
    
    # Plot setup
    plt.figure(figsize=(20, 12))
    
    # Create heatmap with improved aesthetics
    sns.heatmap(matrix,
                xticklabels=list(all_aas),
                yticklabels=[f"{i+1}\n{reference_sequence[i]}" for i in positions],
                cmap='viridis',
                cbar_kws={'label': 'Frequency (%)'})
    
    plt.xlabel('Mutated Amino Acids', fontsize=12)
    plt.ylabel('Position (Reference AA)', fontsize=12)
    plt.title('Mutation Frequency Heatmap', fontsize=14, pad=20)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def plot_mutation_frequency_distribution(stats_df, output_file):
    """
    Plot distribution of mutation frequencies
    """
    plt.figure(figsize=(15, 8))
    
    # Create main frequency plot
    sns.barplot(data=stats_df.sort_values('mutation_frequency', ascending=False),
                x='position',
                y='mutation_frequency',
                color='skyblue')
    
    plt.xlabel('Position', fontsize=12)
    plt.ylabel('Mutation Frequency', fontsize=12)
    plt.title('Distribution of Mutation Frequencies Across Positions', fontsize=14)
    plt.xticks(rotation=45)
    
    # Add significance markers
    significant_positions = stats_df[stats_df['adjusted_p_value'] < 0.05]['position']
    for pos in significant_positions:
        plt.plot(stats_df[stats_df['position'] == pos].index, 
                stats_df[stats_df['position'] == pos]['mutation_frequency'],
                'r*', markersize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_detailed_report(stats_df, mutation_data, reference_sequence, total_sequences, output_file):
    """
    Generate comprehensive analysis report
    """
    with open(output_file, 'w') as f:
        f.write("=== Comprehensive Mutation Analysis Report ===\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total sequences analyzed: {total_sequences}\n\n")
        
        f.write("=== Statistical Summary ===\n")
        f.write(f"Average mutation frequency: {stats_df['mutation_frequency'].mean():.3f}\n")
        f.write(f"Median mutation frequency: {stats_df['mutation_frequency'].median():.3f}\n")
        
        f.write("=== Top Mutation Positions (Ordered by Frequency) ===\n")
        top_positions = stats_df.sort_values('mutation_frequency', ascending=False)
        
        for _, row in top_positions.iterrows():
            if row['mutation_frequency'] > 0:
                f.write(f"\nPosition {int(row['position'])} "
                       f"(Reference: {row['reference_aa']}):\n")
                f.write(f"Mutation Frequency: {row['mutation_frequency']:.3f}\n")
                f.write(f"Conservation Score: {row['conservation_score']:.3f}\n")
                f.write(f"Mutations: {row['mutations']}\n")

def main(fasta_file, reference_sequence_file, output_prefix):
    """
    Enhanced main function with additional analysis
    """
    # Read sequences
    print("Reading sequences...")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    reference_sequence = str(next(SeqIO.parse(reference_sequence_file, "fasta")).seq)
    
    # Analyze mutations
    print("Analyzing mutations...")
    position_frequencies, mutation_data, conservation_scores = analyze_mutations(records, reference_sequence)
    
    # Calculate statistics
    print("Calculating statistics...")
    stats_df = calculate_mutation_statistics(position_frequencies, mutation_data, len(records), 
                                          conservation_scores, reference_sequence)
    
    # Save results
    print("Saving results...")
    stats_df.to_csv(f"{output_prefix}_mutation_stats.csv", index=False)
    
    # Generate visualizations
    print("Creating visualizations...")
    plot_enhanced_heatmap(position_frequencies, reference_sequence, 
                         f"{output_prefix}_mutation_heatmap.png")
    plot_mutation_frequency_distribution(stats_df,
                                      f"{output_prefix}_mutation_frequency.png")
    
    # Generate detailed report
    print("Generating detailed report...")
    generate_detailed_report(stats_df, mutation_data, reference_sequence,
                           len(records), f"{output_prefix}_detailed_report.txt")
    
    print("Analysis complete!")
    return stats_df, position_frequencies, mutation_data

if __name__ == "__main__":
    fasta_file = "combined_preprocessed_rbd.fasta"
    reference_file = "wuhan_rbd.fasta"
    output_prefix = "rbd_mutation_analysis"
    
    stats_df, frequencies, mutations = main(fasta_file, reference_file, output_prefix)
