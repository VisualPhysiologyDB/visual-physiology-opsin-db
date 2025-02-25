import pandas as pd
import os

def translate_candidate_sts(report_dir, ref_seq_name, reference_seq):
    """
    Translates candidate sites from a reference sequence to bovine standard equivalents 
    and annotates the importance report with true positions, TMD regions, and amino acids.

    Args:
        report_dir (str): Directory containing the 'importance_report.csv' file.
        ref_seq_name (str): Name of the reference sequence (e.g., 'bovine').
        reference_seq (pd.Series): Pandas Series representing the reference sequence. 
                                     Values are expected to be amino acid characters or 'nan' for gaps.

    Returns:
        pd.DataFrame: Modified DataFrame with 'true_position', 'TMD', and 'amino_acid' columns added.
    """
    transmembrane_domain = ''
    site_counter = 0
    gap_count = 0
    true_positions = []
    amino_acids = []
    transmembrane_domains = []

    tmd_ranges = {
        'bovine': {
            'N-Termina': range(3, 37),
            '1': range(37, 62),
            '2': range(74, 96),
            '3': range(111, 133),
            '4': range(153, 174),
            '5': range(203, 225),
            '6': range(253, 275),
            '7': range(287, 309),
            'CT/EC': None
        },
        'other': {
            'N-Termina': range(3, 34),
            '1': range(34, 59),
            '2': range(71, 97),
            '3': range(110, 132),
            '4': range(152, 173),
            '5': range(200, 225),
            '6': range(262, 284),
            '7': range(294, 315),
            'CT/EC': None
        }
    }
    reference_tmd_ranges = tmd_ranges.get(ref_seq_name, tmd_ranges['other'])

    df = pd.read_csv(f'{report_dir}\importance_report.csv')

    for index, seq_value in reference_seq.items(): # Iterate through series items (index and value)
        seq_value_str = str(seq_value)

        if seq_value_str == 'nan':
            gap_count += 1
            site_counter += 1
            true_positions.append('NA')
            amino_acids.append('-')
            transmembrane_domains.append('NA')
        else:
            site_counter += 1
            translated_site = site_counter - gap_count # Calculate translated site by subtracting gap count from current site count.
                                                        # This aligns the site number to the reference sequence without gaps.
            transmembrane_domain = 'CT/EC' # Default TMD value if no range matches.
            for tmd_name, site_range in reference_tmd_ranges.items():
                if tmd_name == 'CT/EC':
                    continue
                if site_range and translated_site in site_range:
                    transmembrane_domain = tmd_name
                    break

            true_positions.append(str(translated_site))
            amino_acids.append(seq_value_str)
            transmembrane_domains.append(transmembrane_domain)

    true_positions.pop()
    transmembrane_domains.pop()
    amino_acids.pop()
    df['true_position'] = true_positions
    df['TMD'] = transmembrane_domains
    df['amino_acid'] = amino_acids
    df.to_csv(path_or_buf=os.path.join(report_dir, 'importance_report.csv'), index='Feature', mode="w")

    return df