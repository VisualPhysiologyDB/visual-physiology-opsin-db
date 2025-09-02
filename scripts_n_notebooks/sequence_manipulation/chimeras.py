"""
A script for generating chimeric protein sequences.

This script constructs a new protein sequence by stitching together segments from
one or more source proteins. The final chimera can also have specific point
mutations applied.

The script requires a special format to define the chimera:
  Acc1_Start1_End1-Acc2_Start2_End2[Mutation1,Mutation2,...]

Where:
- Segments are defined by 'Accession_Start_End' and separated by '-'.
- 'Accession' can be an NCBI accession or a predefined ancestral name.
- 'Start' and 'End' are 1-based coordinates relative to the reference sequence.
- An optional block of point mutations can be added at the end in brackets,
  separated by commas (e.g., [A123G,F45S]).

Command-Line Usage Example:
  python chimera_refactored.py \\
      --chimera "AncSW1_1_150-AncSW2_151_348[C203A,F204Y]" \\
      --output_file my_chimera.fasta \\
      --reference_accession NM_001014890.2 \\
      --email "your.email@example.com"
"""

import re
import argparse
from Bio import Entrez, SeqIO
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from Bio.Align import substitution_matrices

# --- Global Configuration ---
ENTREZ_EMAIL = 'your.email@example.com'

def fetch_sequence_from_nucleotide(accession):
    """
    Fetches a protein sequence from a nucleotide accession or a predefined dict.

    This function first checks a dictionary of hardcoded sequences. If not found,
    it queries the NCBI nucleotide database, retrieves the GenBank record, and
    parses the translated protein sequence from the 'CDS' feature.

    Args:
        accession (str): The accession name or number. Can also be 'manual'.

    Returns:
        str: The fetched amino acid sequence.
    """
    if accession.lower() == "manual":
        return input("Manual Sequence Request Detected! \nEnter Sequence Here: ")

    try:
        print(f"Fetching '{accession}' from NCBI Nucleotide DB...")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "gb")
        handle.close()
        # Find the protein translation in the features
        for feature in record.features:
            if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                return feature.qualifiers['translation'][0]
        # If no CDS with translation is found
        raise ValueError(f"No CDS translation found for {accession}")
    except Exception as e:
        print(f"Error fetching or parsing {accession}: {e}")
        manual_seq = input(f"Please enter the sequence for '{accession}' manually: ")
        return manual_seq


def parse_chimera_string(chimera_string):
    """
    Parses the complex chimera definition string into structured data.

    Args:
        chimera_string (str): The string defining the chimera.
                              e.g., "Acc1_1_100-Acc2_101_200[A50G]"

    Returns:
        tuple: A tuple containing:
               - list: A list of dicts, where each dict defines a segment.
               - list: A list of point mutations to apply.
    """
    # Regex to separate the main segments from the optional mutation block
    main_match = re.match(r'([^\[\]]+)(\[(.+)\])?', chimera_string)
    if not main_match:
        raise ValueError("Invalid chimera string format.")

    segments_part = main_match.group(1)
    mutations_part = main_match.group(3)

    # Parse segments
    segments = []
    segment_strings = segments_part.split('-')
    for seg_str in segment_strings:
        parts = seg_str.strip().split('_')
        if len(parts) != 3:
            raise ValueError(f"Invalid segment format: '{seg_str}'")
        segments.append({
            'accession': parts[0],
            'start': int(parts[1]),
            'end': int(parts[2])
        })

    # Parse mutations
    mutations = []
    if mutations_part:
        mutations = [m.strip() for m in mutations_part.split(',')]

    return segments, mutations


def apply_point_mutations(sequence, mutations, reference_protein):
    """
    Applies a list of point mutations to a sequence.

    The sequence is first aligned to the reference to ensure correct positioning.

    Args:
        sequence (str): The protein sequence to mutate.
        mutations (list): A list of mutations (e.g., ['A123G']).
        reference_protein (skbio.Protein): The reference for alignment.

    Returns:
        str: The final mutated sequence.
    """
    print("Applying point mutations to the chimera...")
    protein_to_mutate = Protein(sequence)
    substitution_matrix = substitution_matrices.load("BLOSUM62")

    # Align the new chimera to the reference to apply mutations correctly
    alignment, _, _ = global_pairwise_align_protein(
        reference_protein, protein_to_mutate,
        gap_open_penalty=11, gap_extend_penalty=1,
        substitution_matrix=substitution_matrix
    )
    aligned_ref = str(alignment[0])
    aligned_chimera = str(alignment[1])
    mutated_seq_list = list(aligned_chimera)

    for mutation in mutations:
        match = re.match(r'([A-Z])(\d+)([A-Z])', mutation, re.IGNORECASE)
        if not match:
            print(f"  - WARNING: Skipping invalid mutation format '{mutation}'")
            continue

        original_aa, pos, new_aa = match.groups()
        pos = int(pos) - 1  # 0-based index

        # Find the equivalent position in the aligned sequence
        ref_pos_count = 0
        aligned_idx = -1
        for i, char in enumerate(aligned_ref):
            if char != '-':
                if ref_pos_count == pos:
                    aligned_idx = i
                    break
                ref_pos_count += 1

        if aligned_idx == -1:
            print(f"  - WARNING: Position {pos+1} not found in reference. Skipping.")
            continue

        if aligned_chimera[aligned_idx].upper() != original_aa.upper():
            print(f"  - WARNING: Mismatch at {pos+1}. Expected '{original_aa}', "
                  f"found '{aligned_chimera[aligned_idx]}'. Skipping.")
            continue

        mutated_seq_list[aligned_idx] = new_aa
        print(f"  - Applied mutation {mutation}")

    return "".join(mutated_seq_list).replace('-', '')


def create_chimera(segments, point_mutations, reference_accession):
    """
    Constructs the chimera by fetching, aligning, and slicing segments.

    Args:
        segments (list): A list of segment definition dictionaries.
        point_mutations (list): A list of point mutations to apply.
        reference_accession (str): Accession for the reference sequence.

    Returns:
        str: The final amino acid sequence of the chimera.
    """
    print(f"Fetching reference sequence '{reference_accession}'...")
    reference_protein = Protein(fetch_sequence_from_nucleotide(reference_accession))
    substitution_matrix = substitution_matrices.load("BLOSUM62")
    chimera_pieces = []

    print("\nProcessing chimera segments...")
    for seg in segments:
        acc, start, end = seg['accession'], seg['start'], seg['end']
        print(f"  - Segment from {acc}, region {start}-{end}")

        # Fetch the sequence for the current segment
        wt_protein = Protein(fetch_sequence_from_nucleotide(acc))

        # Align it to the reference to get coordinates right
        alignment, _, _ = global_pairwise_align_protein(
            reference_protein, wt_protein,
            gap_open_penalty=11, gap_extend_penalty=1,
            substitution_matrix=substitution_matrix
        )
        aligned_ref = str(alignment[0])
        aligned_wt = str(alignment[1])

        # Find the slice indices in the *aligned* sequence
        ref_pos_count = 0
        slice_start_idx, slice_end_idx = -1, -1

        for i, char in enumerate(aligned_ref):
            if char != '-':
                ref_pos_count += 1
                if ref_pos_count == start:
                    slice_start_idx = i
                if ref_pos_count == end:
                    slice_end_idx = i
                    break # Found the end of the slice

        if slice_start_idx == -1 or slice_end_idx == -1:
            print(f"    - WARNING: Could not find range {start}-{end} in alignment. Skipping segment.")
            continue

        # Extract the segment and add to our list
        segment_sequence = aligned_wt[slice_start_idx : slice_end_idx + 1]
        chimera_pieces.append(segment_sequence)

    # Join the pieces and remove any alignment gaps
    raw_chimera_seq = "".join(chimera_pieces).replace('-', '')

    # Apply point mutations if they exist
    if point_mutations:
        final_sequence = apply_point_mutations(
            raw_chimera_seq, point_mutations, reference_protein
        )
    else:
        final_sequence = raw_chimera_seq

    return final_sequence


def main():
    """Main function to parse arguments and run the chimera workflow."""
    parser = argparse.ArgumentParser(
        description='A script for generating chimeric proteins.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-co', '--chimera',
        required=True,
        help="Chimeric protein definition string.\n"
             "Format: Acc1_Start1_End1-Acc2_Start2_End2[Mutation1,etc...]"
    )
    parser.add_argument(
        '-o', '--output_file',
        required=True,
        help='Path to the output FASTA file.'
    )
    parser.add_argument(
        '-ra', '--reference_accession',
        default='NM_001014890.2', # Bos taurus rhodopsin DNA
        help='Reference NUCLEOTIDE accession for sequence numbering.\n'
             'Default: NM_001014890.2 (Bos taurus rh1).'
    )
    parser.add_argument(
        '--email',
        default=ENTREZ_EMAIL,
        help='Your email address for NCBI Entrez queries.'
    )
    args = parser.parse_args()

    # --- Setup ---
    if args.email != 'your.email@example.com':
        Entrez.email = args.email
    else:
        print("Warning: Using default Entrez email. Please provide your own with the --email flag.")

    # --- Core Logic ---
    try:
        # 1. Parse the input string
        segments, point_mutations = parse_chimera_string(args.chimera)

        # 2. Build the chimera
        final_chimera_sequence = create_chimera(
            segments, point_mutations, args.reference_accession
        )

        # 3. Write to output file
        if final_chimera_sequence:
            chimera_name = re.sub(r'[\[\]]', '_', args.chimera).replace(',', '_')
            with open(args.output_file, 'w') as f:
                f.write(f'>{chimera_name}\n')
                f.write(f'{final_chimera_sequence}\n')
            print(f"\nSuccessfully generated chimera '{chimera_name}'.")
            print(f"Output saved to '{args.output_file}'.")
        else:
            print("\nCould not generate chimera sequence due to previous warnings.")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == '__main__':
    main()
