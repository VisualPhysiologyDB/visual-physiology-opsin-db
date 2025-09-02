"""
A script for performing in-silico mutagenesis on protein sequences.

This script can generate mutated protein sequences based on a set of specified
mutations against a wild-type sequence. It can fetch sequences from NCBI or use
a predefined set of ancestral sequences.

The main functionalities are:
1.  Generate all combinatorial mutants from a list of mutations against a
    wild-type accession.
2.  Take a predefined list of mutant accessions (in the format
    'ACCESSION_M1,M2,...') and generate the corresponding mutated protein
    sequences.

The script aligns the target sequence to a reference sequence (e.g., Bos taurus)
to ensure mutation positions are correct, even with insertions or deletions.

Command-Line Usage Examples:

1.  Generate all combinations from a list of mutations:
    python mutagenesis_refactored.py \\
        --wt_accession AncBovine \\
        --mutations "A116S,S119A,G121A" \\
        --output_file combined_mutants.fasta \\
        --reference_accession NM_001014890.2

2.  Generate sequences for specific mutants from a file:
    python mutagenesis_refactored.py \\
        --mutant_file my_mutants.txt \\
        --output_file specific_mutants.fasta

3.  Generate a sequence for a single mutant:
    python mutagenesis_refactored.py \\
        --mutant_accession "AncBovine_A116S,G121A" \\
        --output_file single_mutant.fasta
"""

import re
import argparse
import itertools
import pandas as pd
from Bio import Entrez, SeqIO
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from Bio.Align import substitution_matrices

# --- Global Configuration ---

# It's good practice for the user of the script to provide their email
# for NCBI Entrez.
ENTREZ_EMAIL = 'your.email@example.com'
Entrez.email = ENTREZ_EMAIL
CACHED_SEQUENCES = {}



def fetch_protein_sequence(accession, sequence_type='target'):
    """
    Fetches a protein sequence from various sources.

    Order of operations:
    1. Checks for manually provided sequence string.
    2. Checks the hardcoded ancestral sequence dictionary.
    3. Queries NCBI Entrez for the accession number.
    4. If all else fails, prompts the user for a manual sequence string.

    Args:
        accession (str): The accession name, number, or 'manual'.
        sequence_type (str): The type of sequence being fetched ('reference'
                             or 'target'), used for manual input prompt.

    Returns:
        str: The fetched amino acid sequence.
    """
    if accession.lower() == "manual":
        return input(f"Enter {sequence_type.capitalize()} Sequence: ")

    try:
        print(f"Fetching '{accession}' from NCBI...")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "gb")
        handle.close()
        for i,feature in enumerate(record.features):
            if feature.type=='CDS':
                aa = feature.qualifiers['translation'][0]
        
        return aa
    except Exception as e:
        print(f"Warning: Could not fetch accession '{accession}' from NCBI. Error: {e}")
        manual_seq = input(f"Please enter the {sequence_type} sequence for '{accession}' manually: ")
        return manual_seq


def get_mutant_combinations(wt_accession, mutations_list):
    """
    Generates a list of mutant accession strings from all combinations of mutations.

    For a wt_accession 'WT' and mutations ['A1B', 'C2D'], it will generate:
    ['WT_A1B', 'WT_C2D', 'WT_A1B,C2D']

    Args:
        wt_accession (str): The accession name of the wild-type sequence.
        mutations_list (list): A list of mutations as strings (e.g., ['A123G', 'F45S']).

    Returns:
        list: A list of formatted mutant accession strings.
    """
    mutant_accessions = []
    for i in range(1, len(mutations_list) + 1):
        # Get all combinations of mutations of length i
        combinations = itertools.combinations(mutations_list, i)
        for combo in combinations:
            # Format the accession string: e.g., "MyProtein_A1B,C2D"
            mut_str = ",".join(combo)
            mutant_accessions.append(f"{wt_accession}_{mut_str}")
    return mutant_accessions


def parse_mutant_accession(mutant_accession):
    """
    Parses a mutant accession string into its components.

    Args:
        mutant_accession (str): The string to parse (e.g., "MyProtein_A1B,C2D").

    Returns:
        tuple: A tuple containing (wild_type_accession, list_of_mutations).
               Returns (mutant_accession, []) if no mutations are found.
    """
    mutant_accession = mutant_accession.strip()
    count = mutant_accession.count('_')
    
    if count == 1 and len(mutant_accession.split('_')[0]) > 2:
        parts = mutant_accession.split('_', 1)
        wt_accession = parts[0]
        mutations_str = parts[1]
        # Handle single or multiple mutations
        mutations = mutations_str.split(',')
        return wt_accession, mutations
    elif count == 2:
        parts = mutant_accession.split('_', 2)
        wt_accession = f'{parts[0]}_{parts[1]}'
        mutations_str = parts[2]
        # Handle single or multiple mutations
        mutations = mutations_str.split(',')
        return wt_accession, mutations
    else:
        # It's a wild-type sequence with no mutations
        return mutant_accession, []


def get_mutant_seqs(
    mutant_accessions,
    output_file,
    reference_accession,
    output_format='fasta',
    allow_wt=True,
    email='your.email@example.com'
):
    """
    Generates mutated sequences and writes them to an output file.

    For each mutant accession string, this function will:
    1. Parse the wild-type accession and the required mutations.
    2. Fetch the wild-type and reference protein sequences.
    3. Align the WT to the reference to correctly map mutation sites.
    4. Apply the mutations.
    5. Write the final sequence to the output file.

    Args:
        mutant_accessions (list): List of mutant accession strings to process.
        output_file (str): Path to the output file.
        reference_accession (str): Accession of the sequence for numbering.
        output_format (str): 'fasta' or 'tsv'.
        allow_wt (bool): If True, allows processing of sequences with no mutations.
    """
    if email != 'your.email@example.com':
        Entrez.email = email
    CURRENT_WT = ''
    
    if reference_accession not in CACHED_SEQUENCES.keys():
        print(f"Fetching reference sequence '{reference_accession}'...")
        CACHED_SEQUENCES[reference_accession] = Protein(fetch_protein_sequence(reference_accession, 'reference'))

    reference_protein = CACHED_SEQUENCES[reference_accession]

    # Prepare output file
    if output_format == 'tsv':
        with open(output_file, 'w') as f:
            f.write('Accession\tSequence\n')

    print("\nProcessing mutations...")
    for mutant_acc in mutant_accessions:
        wt_acc, mutations = parse_mutant_accession(mutant_acc)

        if not mutations:
            if allow_wt:
                print(f"Processing Wild-Type: {wt_acc}")
                if wt_acc not in CACHED_SEQUENCES.keys():
                    CACHED_SEQUENCES[wt_acc] = Protein(fetch_protein_sequence(wt_acc, 'target'))
                wt_seq = CACHED_SEQUENCES[wt_acc]
                final_seq = wt_seq
                final_acc = wt_acc
            else:
                print(f"Skipping Wild-Type sequence '{wt_acc}' as per settings.")
                continue
        else:
            print(f"Processing: {mutant_acc}")
            if wt_acc not in CACHED_SEQUENCES.keys():
                CACHED_SEQUENCES[wt_acc] = Protein(fetch_protein_sequence(wt_acc, 'target'))
            
            wt_protein = CACHED_SEQUENCES[wt_acc]

            # Align the WT to the reference to get the correct numbering
            if wt_protein != CURRENT_WT:
                substitution_matrix = substitution_matrices.load("BLOSUM62")
                alignment, _, _ = global_pairwise_align_protein(
                    reference_protein, wt_protein,
                    gap_open_penalty=11,
                    gap_extend_penalty=1,
                    substitution_matrix=substitution_matrix
                )
                
                aligned_ref = str(alignment[0])
                aligned_wt = str(alignment[1])
                CURRENT_WT = wt_protein
            
            mutated_seq_list = list(aligned_wt)

            for mutation in mutations:
                # Regex to parse mutations like 'A123G'
                match = re.match(r'([A-Z])(\d+)([A-Z])', mutation, re.IGNORECASE)
                if not match:
                    print(f"  - WARNING: Skipping invalid mutation format '{mutation}'")
                    continue
                
                original_aa, pos, new_aa = match.groups()
                pos = int(pos) - 1  # Convert to 0-based index

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
                    print(f"  - WARNING: Position {pos+1} not found in reference alignment. Skipping.")
                    continue
                
                # Check if the original AA matches the sequence at that position
                if aligned_wt[aligned_idx].upper() != original_aa.upper():
                    print(f"  - WARNING: Mismatch at position {pos+1}. "
                          f"Expected '{original_aa}', found '{aligned_wt[aligned_idx]}'. "
                          "This can happen if the reference and target are very divergent. Skipping mutation.")
                    continue
                
                # Apply mutation
                mutated_seq_list[aligned_idx] = new_aa

            # Create final sequence by removing gaps
            final_seq = "".join(mutated_seq_list).replace('-', '')
            final_acc = mutant_acc

        # Write to output file
        with open(output_file, 'a+') as f:
            if output_format == 'fasta':
                f.write(f'>{final_acc}\n{final_seq}\n')
            else:
                f.write(f'{final_acc}\t{final_seq}\n')
        print(f"  -> Saved as {final_acc}")

    print(f"\nProcessing complete. Output saved to '{output_file}'.")


def main():
    """Main function to parse arguments and run the mutagenesis workflow."""
    parser = argparse.ArgumentParser(
        description='A script for in-silico site-directed mutagenesis.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
            This script allows two main modes of operation:
            1. Generate all combinations of mutants from a base wild-type sequence.
            (Use --wt_accession and --mutations)
            2. Generate sequences for a pre-defined list of mutants.
            (Use --mutant_file or --mutant_accession)
        """
    )

    # --- Input Modes (Mutually Exclusive) ---
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--mutant_file',
        help='Path to a text file containing mutant accessions, one per line.'
    )
    input_group.add_argument(
        '--mutant_accession',
        help="A single mutant accession string (e.g., 'MyProtein_A123G,F45S')."
    )
    input_group.add_argument(
        '--wt_accession',
        help='Wild-type accession to generate mutants from (used with --mutations).'
    )

    # --- Arguments for Combination Generation ---
    parser.add_argument(
        '--mutations',
        help='Comma-separated list of mutations for combinatorial generation (e.g., "A123G,F45S").\nRequired if --wt_accession is used.'
    )

    # --- General Configuration ---
    parser.add_argument(
        '-o', '--output_file',
        required=True,
        help='Path to the output file for the generated sequences.'
    )
    parser.add_argument(
        '-ra', '--reference_accession',
        default='NP_001014890.1', # Bos taurus rhodopsin protein
        help='Reference accession for sequence numbering. Default: NP_001014890.1 (Bos taurus rh1).'
    )
    parser.add_argument(
        '--output_format',
        choices=['fasta', 'tsv'],
        default='fasta',
        help="Format for the output file. Default: 'fasta'."
    )
    parser.add_argument(
        '--no_wt',
        action='store_true',
        help='Flag to prevent processing of wild-type sequences (those with no mutations).'
    )
    parser.add_argument(
        '--email',
        default=ENTREZ_EMAIL,
        help='Your email address for NCBI Entrez queries.'
    )

    args = parser.parse_args()

    # --- Argument Validation ---
    if args.email != 'your.email@example.com':
        Entrez.email = args.email
    else:
        print("Warning: Using default Entrez email. Please provide your own with the --email flag.")

    if args.wt_accession and not args.mutations:
        parser.error('--mutations is required when using --wt_accession.')

    mutant_list = []

    # --- Workflow Selection ---
    if args.wt_accession:
        # Mode 1: Generate combinations
        print("Mode: Generating mutant combinations.")
        mutations_list = [m.strip() for m in args.mutations.split(',')]
        mutant_list = get_mutant_combinations(args.wt_accession, mutations_list)
        # Also add the WT to the list if desired
        if not args.no_wt:
            mutant_list.insert(0, args.wt_accession)

    elif args.mutant_file:
        # Mode 2a: Process mutants from a file
        print(f"Mode: Processing mutants from file '{args.mutant_file}'.")
        with open(args.mutant_file, 'r') as f:
            mutant_list = [line.strip() for line in f if line.strip()]

    elif args.mutant_accession:
        # Mode 2b: Process a single mutant
        print(f"Mode: Processing single mutant '{args.mutant_accession}'.")
        mutant_list = [args.mutant_accession]

    # --- Execute Core Logic ---
    if mutant_list:
        get_mutant_seqs(
            mutant_accessions=mutant_list,
            output_file=args.output_file,
            reference_accession=args.reference_accession,
            output_format=args.output_format,
            allow_wt=not args.no_wt
        )
    else:
        print("No mutants to process. Exiting.")


if __name__ == '__main__':
    main()