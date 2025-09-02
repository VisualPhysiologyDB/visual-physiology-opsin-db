"""
Generates reciprocal mutants between two aligned DNA sequences and writes them
to a new FASTA file.

This script identifies differing amino acid positions between two sequences
(seq2 and seq3) in an alignment and creates mutants by swapping these
amino acids. A third sequence (seq1) is used as a positional reference for
naming the mutations.

The script can be run from the command line or imported as a module.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Tuple

def _get_ref_pos_label(aln_pos: int, aligned_to_ungapped_ref: Dict[int, int]) -> str:
    """
    Generates a position label for a mutation based on a reference sequence.

    If the alignment position corresponds to a residue in the reference, it returns
    the 1-based position. If it's a gap, it finds the nearest preceding residue
    and appends 'gap'.

    Args:
        aln_pos: The 0-based index in the alignment.
        aligned_to_ungapped_ref: A mapping from aligned position to the ungapped
                                 index of the reference sequence.

    Returns:
        A string label for the position (e.g., "127" or "126gap").
    """
    if aln_pos in aligned_to_ungapped_ref:
        return str(aligned_to_ungapped_ref[aln_pos] + 1)  # Return 1-based position
    
    # If the position is a gap in the reference, backtrack to find the nearest
    # non-gap residue to create a label like "126gap".
    for back_pos in range(aln_pos - 1, -1, -1):
        if back_pos in aligned_to_ungapped_ref:
            return f"{aligned_to_ungapped_ref[back_pos] + 1}gap"
            
    return "1gap"  # Fallback if no preceding residue is found

def generate_reciprocal_mutants(records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Takes a list of three aligned sequences and generates reciprocal mutants.

    It compares the second and third sequences, finds differences, and creates
    mutants for each. The first sequence is used as a reference for naming.

    Args:
        records: A list of three Biopython SeqRecord objects from an alignment.

    Returns:
        A list of SeqRecord objects containing the three original (ungapped)
        sequences and all generated mutants.
        
    Raises:
        AssertionError: If the input does not contain exactly three sequences.
    """
    assert len(records) == 3, "Input must contain exactly three aligned sequences."

    # Unpack records and get aligned sequences as strings
    ref_record, seq2_record, seq3_record = records
    aligned_ref = str(ref_record.seq)
    aligned_2 = str(seq2_record.seq)
    aligned_3 = str(seq3_record.seq)

    # Create ungapped versions of the sequences
    ungapped_ref_seq = aligned_ref.replace("-", "")
    ungapped_seq_2 = aligned_2.replace("-", "")
    ungapped_seq_3 = aligned_3.replace("-", "")

    # --- Build mappings from aligned position to ungapped index ---
    # This is crucial for correctly placing the mutation in the final sequence.
    aligned_to_ungapped_ref: Dict[int, int] = {}
    aligned_to_ungapped_2: Dict[int, int] = {}
    aligned_to_ungapped_3: Dict[int, int] = {}
    
    ungapped_idx_ref, ungapped_idx_2, ungapped_idx_3 = 0, 0, 0

    for i in range(len(aligned_ref)):
        if aligned_ref[i] != "-":
            aligned_to_ungapped_ref[i] = ungapped_idx_ref
            ungapped_idx_ref += 1
        if aligned_2[i] != "-":
            aligned_to_ungapped_2[i] = ungapped_idx_2
            ungapped_idx_2 += 1
        if aligned_3[i] != "-":
            aligned_to_ungapped_3[i] = ungapped_idx_3
            ungapped_idx_3 += 1

    # --- Find differences between seq2 and seq3 ---
    # We only care about positions where neither sequence has a gap.
    diffs: List[Tuple[int, str, str]] = []
    for i, (aa_2, aa_3) in enumerate(zip(aligned_2, aligned_3)):
        if aa_2 != "-" and aa_3 != "-" and aa_2 != aa_3:
            diffs.append((i, aa_2, aa_3))  # (alignment_index, AA_in_seq2, AA_in_seq3)

    # --- Generate mutant SeqRecord objects ---
    mutants: List[SeqRecord] = []

    # Mutate sequence 2 to be like sequence 3 at different positions
    for aln_pos, aa_2, aa_3 in diffs:
        ungapped_pos = aligned_to_ungapped_2[aln_pos]
        mutant_seq_list = list(ungapped_seq_2)
        mutant_seq_list[ungapped_pos] = aa_3
        
        pos_label = _get_ref_pos_label(aln_pos, aligned_to_ungapped_ref)
        mutant_id = f"{seq2_record.id}_{aa_2}{pos_label}{aa_3}"
        
        mutant_record = SeqRecord(Seq("".join(mutant_seq_list)), id=mutant_id, description='')
        mutants.append(mutant_record)

    # Mutate sequence 3 to be like sequence 2 at different positions
    for aln_pos, aa_2, aa_3 in diffs:
        ungapped_pos = aligned_to_ungapped_3[aln_pos]
        mutant_seq_list = list(ungapped_seq_3)
        mutant_seq_list[ungapped_pos] = aa_2
        
        pos_label = _get_ref_pos_label(aln_pos, aligned_to_ungapped_ref)
        mutant_id = f"{seq3_record.id}_{aa_3}{pos_label}{aa_2}"
        
        mutant_record = SeqRecord(Seq("".join(mutant_seq_list)), id=mutant_id, description='')
        mutants.append(mutant_record)

    # --- Combine original sequences with mutants for final output ---
    originals = [
        SeqRecord(Seq(ungapped_ref_seq), id=ref_record.id, description=''),
        SeqRecord(Seq(ungapped_seq_2), id=seq2_record.id, description=''),
        SeqRecord(Seq(ungapped_seq_3), id=seq3_record.id, description=''),
    ]

    return originals + mutants

def main():
    """Main function to run the script from the command line."""
    parser = argparse.ArgumentParser(
        description="Generates reciprocal mutants from an aligned FASTA file containing three sequences."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input aligned FASTA file (e.g., 'test3_seqs_ALN.fasta')."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path for the output FASTA file (e.g., 'all_mutants.fasta')."
    )
    args = parser.parse_args()

    try:
        # Load aligned sequences from the input file
        records = list(SeqIO.parse(args.input_file, "fasta"))
        
        # Generate the mutants using the core function
        all_sequences_to_write = generate_reciprocal_mutants(records)
        
        # Write the results to the output file
        SeqIO.write(all_sequences_to_write, args.output_file, "fasta")
        
        num_mutants = len(all_sequences_to_write) - 3
        print(f"✅ Success! Generated {num_mutants} mutants and included 3 original sequences.")
        print(f"Output written to {args.output_file}")

    except FileNotFoundError:
        print(f"❌ Error: Input file not found at '{args.input_file}'")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()