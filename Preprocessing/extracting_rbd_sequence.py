from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def find_rbd_boundaries(sequence, start_pattern, end_pattern):
    """
    Find the start and end positions of RBD based on sequence patterns.
    
    Args:
        sequence: Full protein sequence string
        start_pattern: String of amino acids marking RBD start
        end_pattern: String of amino acids marking RBD end
    
    Returns:
        Tuple of (start_index, end_index) or None if patterns not found
    """
    try:
        start_idx = sequence.index(start_pattern)
        end_idx = sequence.index(end_pattern) + len(end_pattern)
        return start_idx, end_idx
    except ValueError:
        return None

def extract_rbd_region(records, start_pattern, end_pattern):
    """
    Extract the RBD region from protein sequences using pattern matching.
    
    Args:
        records: List of SeqRecord objects
        start_pattern: String of amino acids marking RBD start
        end_pattern: String of amino acids marking RBD end
    
    Returns:
        Tuple of (list of SeqRecord objects containing RBD region, number of skipped sequences)
    """
    rbd_records = []
    skipped = 0
    pattern_not_found = 0
    
    print(f"\nExtracting RBD region using patterns:")
    print(f"Start pattern: {start_pattern}")
    print(f"End pattern: {end_pattern}")
    
    for record in records:
        sequence = str(record.seq)
        boundaries = find_rbd_boundaries(sequence, start_pattern, end_pattern)
        
        if boundaries:
            start_idx, end_idx = boundaries
            rbd_sequence = sequence[start_idx:end_idx]
            
            # Create new record with RBD sequence
            rbd_record = SeqRecord(
                Seq(rbd_sequence),
                id=record.id,
                description=f"{record.description} | RBD region (pattern-based extraction)"
            )
            rbd_records.append(rbd_record)
        else:
            pattern_not_found += 1
            skipped += 1
    
    print(f"\nPatterns not found in {pattern_not_found} sequences")
    return rbd_records, skipped

def main(input_fasta, output_fasta, start_pattern, end_pattern):
    """Main function to process the FASTA file and extract RBD regions."""
    # Read input sequences
    print(f"Reading sequences from {input_fasta}...")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total_sequences = len(records)
    print(f"Total sequences read: {total_sequences:,}")
    
    # Extract RBD regions
    rbd_records, skipped = extract_rbd_region(records, start_pattern, end_pattern)
    
    # Write RBD sequences to new file
    SeqIO.write(rbd_records, output_fasta, "fasta")
    
    # Print statistics
    print("\nProcessing Statistics:")
    print("-" * 50)
    print(f"Input sequences:     {total_sequences:,}")
    print(f"Sequences processed: {len(rbd_records):,}")
    print(f"Sequences skipped:   {skipped:,}")
    print(f"\nRBD sequences written to: {output_fasta}")

if __name__ == "__main__":
    input_file = "xec/xec_preprocessed.fasta"
    output_file = "xec/xec_preprocessed_rbd.fasta"

    START_PATTERN = "RVQP" 
    END_PATTERN = "CVNF"
    
    # Process the sequences
    main(input_file, output_file, START_PATTERN, END_PATTERN)
