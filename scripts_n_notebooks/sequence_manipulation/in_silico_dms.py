"""
A script for in-silico Deep Mutational Scanning (DMS).

This script generates a complete single-mutant library for a given set of
sites on a wild-type protein. For each specified site, it creates 20 new
sequences, one for each possible amino acid substitution.

The script aligns the target (wild-type) sequence to a reference sequence
(e.g., Bos taurus rhodopsin) to ensure mutation positions are mapped correctly,
even if the target sequence has insertions or deletions relative to the reference.

Command-Line Usage Example:

  python dms_scanner.py \\
      --wt_accession "AncRho1" \\
      --sites "S121,A185,G203" \\
      --output_file dms_library.fasta \\
      --reference_accession "NM_001014890.2" \\
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
# The 20 standard amino acids for substitution
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def fetch_protein_sequence(accession, sequence_type='target'):
    """
    Fetches a protein sequence from a nucleotide accession or a predefined dict.

    This function queries the NCBI nucleotide database, retrieves the GenBank
    record, and parses the translated protein sequence from the 'CDS' feature.

    Args:
        accession (str): The accession name or number. Can also be 'manual'.
        sequence_type (str): The type of sequence being fetched ('reference'
                             or 'target'), used for manual input prompt.

    Returns:
        str: The fetched amino acid sequence.
    """
    if accession.lower() == "manual":
        return input(f"Enter {sequence_type.capitalize()} Sequence: ")

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


def generate_dms_library(
    wt_accession,
    sites,
    output_file,
    reference_accession,
    email='your.email@example.com'):
    """
    Generates a deep mutational scanning library for specified sites.

    For each site, this function creates 20 mutants, one for each amino acid.

    Args:
        wt_accession (str): Accession of the wild-type sequence to mutate.
        sites (list): A list of sites to scan (e.g., ['S121', 'A185']).
        output_file (str): Path to the output FASTA file.
        reference_accession (str): Accession of the sequence for numbering.
    """
    if email != 'your.email@example.com':
        Entrez.email = email
    # --- 1. Fetch Sequences and Align ---
    print(f"Fetching reference sequence '{reference_accession}'...")
    reference_protein = Protein(fetch_protein_sequence(reference_accession, 'reference'))

    print(f"Fetching wild-type sequence '{wt_accession}'...")
    wt_protein = Protein(fetch_protein_sequence(wt_accession, 'target'))

    print("Aligning WT to reference to map coordinates...")
    substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignment, _, _ = global_pairwise_align_protein(
        reference_protein, wt_protein,
        gap_open_penalty=11,
        gap_extend_penalty=1,
        substitution_matrix=substitution_matrix
    )
    aligned_ref = str(alignment[0])
    aligned_wt = str(alignment[1])

    # Prepare the output file
    with open(output_file, 'w') as f:
        f.write('') # Create an empty file to overwrite previous runs

    print("\n--- Generating Mutant Library ---")
    # --- 2. Iterate Through Each Site for Scanning ---
    for site in sites:
        match = re.match(r'([A-Z])(\d+)', site, re.IGNORECASE)
        if not match:
            print(f"WARNING: Skipping invalid site format '{site}'. Must be like 'A123'.")
            continue

        expected_aa, pos = match.groups()
        pos = int(pos) - 1  # Convert to 0-based index

        # --- 3. Find and Validate the Site in the Alignment ---
        ref_pos_count = 0
        aligned_idx = -1
        for i, char in enumerate(aligned_ref):
            if char != '-': # Count only non-gap positions in the reference
                if ref_pos_count == pos:
                    aligned_idx = i
                    break
                ref_pos_count += 1

        if aligned_idx == -1:
            print(f"WARNING: Position {pos + 1} not found in reference alignment. Skipping site '{site}'.")
            continue

        actual_aa = aligned_wt[aligned_idx]
        if actual_aa.upper() != expected_aa.upper():
            print(f"WARNING: Mismatch at site {site}. Expected '{expected_aa}' based on input, "
                  f"but found '{actual_aa}' in the sequence. Skipping this site.")
            continue

        print(f"Scanning site: {site} (Found '{actual_aa}' at position {pos + 1})")
        
        # --- 4. Generate all 20 Mutations for the Validated Site ---
                
        # Find the corresponding index in the *unaligned* WT sequence
        # by counting non-gap characters up to the aligned index.
        unaligned_idx = len(aligned_wt[:aligned_idx+1].replace('-', '')) - 1

        for new_aa in AMINO_ACIDS:
            if actual_aa == new_aa:
                pass
            else:
                wt_protein = str(wt_protein)
                mutant_seq_list = list(wt_protein)
                mutant_seq_list[unaligned_idx] = new_aa
                mutant_sequence = "".join(mutant_seq_list)

                # Create a descriptive name for the FASTA header
                mutation_name = f"{actual_aa}{pos + 1}{new_aa}"
                fasta_header = f">{wt_accession}_{mutation_name}"

                # Append the new mutant to the output file
                with open(output_file, 'a') as f:
                    f.write(f"{fasta_header}\n")
                    f.write(f"{mutant_sequence}\n")
        
        print(f"  -> Generated all mutants for site {site} and saved to '{output_file}'.")

    print("\nDeep mutational scanning library generation complete.")


def main():
    """Main function to parse arguments and run the DMS workflow."""
    parser = argparse.ArgumentParser(
        description='A script for in-silico Deep Mutational Scanning (DMS).',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--wt_accession',
        required=True,
        help='Accession ID for the wild-type sequence to be mutated.'
    )
    parser.add_argument(
        '--sites',
        required=True,
        help="Comma-separated list of sites to scan (e.g., 'S121,A185,G203').\n"
             "The format is OriginalAminoAcid-Position."
    )
    parser.add_argument(
        '-o', '--output_file',
        required=True,
        help='Path to the output FASTA file for the generated library.'
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

    sites_list = [s.strip() for s in args.sites.split(',')]

    # --- Execute Core Logic ---
    generate_dms_library(
        wt_accession=args.wt_accession,
        sites=sites_list,
        output_file=args.output_file,
        reference_accession=args.reference_accession
    )


if __name__ == '__main__':
    main()
