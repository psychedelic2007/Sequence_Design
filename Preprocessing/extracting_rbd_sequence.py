from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def find_rbd_boundaries(sequence, start_pattern, end_pattern):
    try:
        start_idx = sequence.index(start_pattern)
        end_idx = sequence.index(end_pattern) + len(end_pattern)
        return start_idx, end_idx
    except ValueError:
        return None

def extract_rbd_region(records, start_pattern, end_pattern):
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
    print(f"Reading sequences from {input_fasta}...")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total_sequences = len(records)
    print(f"Total sequences read: {total_sequences:,}")
    
    rbd_records, skipped = extract_rbd_region(records, start_pattern, end_pattern)
    
    SeqIO.write(rbd_records, output_fasta, "fasta")
    
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

    main(input_file, output_file, START_PATTERN, END_PATTERN)
