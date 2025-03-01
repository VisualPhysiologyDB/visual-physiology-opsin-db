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
    if ref_seq_name not in ['bovine', 'squid']:
        raise ValueError("ref_seq_name must be either 'bovine' or 'squid'")

    try:
        df = pd.read_csv(os.path.join(report_dir, 'importance_report.csv'))
    except FileNotFoundError:
        raise FileNotFoundError(f"importance_report.csv not found in {report_dir}")
    
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
    df.to_csv(path_or_buf=os.path.join(report_dir, 'translated_importance_report.csv'), index='Feature', mode="w")

    return df

import pandas as pd
import os

def translate_candidate_sts_aa_props(report_dir, reference_seq, ref_seq_name):
    """Translates candidate STS positions to their equivalent positions in a reference sequence (bovine or squid).

    This function takes an importance report (CSV) and a reference sequence, 
    calculates the corresponding position in the reference sequence,
    and determines the transmembrane domain (TMD) based on pre-defined ranges.

    Args:
        report_dir (str): The directory containing the 'importance_report.csv' file.
        reference_seq (pandas.Series): A Series containing the reference sequence (e.g., bovine or squid).  
                                        Each element represents an amino acid or a gap ('-').
        ref_seq_name (str): The name of the reference sequence ('bovine' or 'squid').  This determines
                            the TMD ranges used.

    Returns:
        pandas.DataFrame: The updated 'importance_report.csv' DataFrame with added columns:
                          'true_position', 'TMD', and 'amino_acid'.  The DataFrame is also
                          saved back to the same CSV file.
    Raises:
        ValueError: If ref_seq_name is not 'bovine' or 'squid'.
        FileNotFoundError: If 'importance_report.csv' is not found in report_dir.

    """

    # Input Validation and Setup
    if (ref_seq_name != 'bovine' and ref_seq_name != 'squid'):
        raise ValueError("ref_seq_name must be either 'bovine' or 'squid'")

    if ref_seq_name == 'bovine':
        tmd_ranges = {
            'N-Termina': (3, 37),
            '1': (37, 62),
            '2': (74, 96),
            '3': (111, 133),
            '4': (153, 174),
            '5': (203, 225),
            '6': (253, 275),
            '7': (287, 309),
            'CT/EC': (None, None)  # Special case for CT/EC
        }
    elif ref_seq_name == 'squid':
        tmd_ranges = {
            'N-Termina': (3, 34),
            '1': (34, 59),
            '2': (71, 97),
            '3': (110, 132),
            '4': (152, 173),
            '5': (200, 225),
            '6': (262, 284),
            '7': (294, 315),
            'CT/EC': (None, None)  # Special case for CT/EC
        }
        
    try:
        df = pd.read_csv(os.path.join(report_dir, 'importance_report.csv'))
    except FileNotFoundError:
        raise FileNotFoundError(f"importance_report.csv not found in {report_dir}")

    # --- Feature Position Extraction ---
    indx_list = df['feature'].tolist()
    pos_list = [int(x.split('_')[0]) for x in indx_list]  # More concise list comprehension
    df['feat_pos'] = pos_list
    df = df.sort_values('feat_pos', ascending=True).reset_index(drop=True)

    # --- Translation Logic ---
    true_pos = []
    aa = []
    tmd = []
    gaps = 0
    k = 1  # k is the index in the *reference* sequence

    for rows in reference_seq.values:
        rows = str(rows)  # Handle potential mixed data types consistently

        if rows == 'nan' or rows == 'NaN' or rows == '-' or rows is None:
            gaps += 1
            true_pos.append('NA')
            aa.append('-')
            tmd.append('NA')
        else:
            trans_site = k - gaps  # Calculate translated site
            # Determine TMD based on ref_seq_name

            tm = 'CT/EC' # Default value
            for region, (start, end) in tmd_ranges.items():
                if start is not None and end is not None and trans_site in range(start, end):
                    tm = region
                    break  # Important: exit inner loop once a match is found

            true_pos.append(str(trans_site))
            aa.append(rows)
            tmd.append(tm)

        k += 1  # Increment *after* processing the current position

    # --- Merge Translated Data ---
    #  Create a lookup dictionary for efficient merging
    lookup_dict = {pos: (tp, t, a) for pos, tp, t, a in zip(range(1, len(true_pos) + 1), true_pos, tmd, aa)}

    exp_true_pos = []
    exp_tmd = []
    exp_aa = []
    for feat in df['feat_pos']:
        if feat in lookup_dict:
            tp, t, a = lookup_dict[feat]
            exp_true_pos.append(tp)
            exp_tmd.append(t)
            exp_aa.append(a)
        else:
            # Handle the case where feat is not in lookup_dict.
            # This could happen if there's a mismatch in lengths.
            #  Append 'NA' or some other placeholder.
            exp_true_pos.append('NA')
            exp_tmd.append('NA')
            exp_aa.append('NA')



    df['true_position'] = exp_true_pos
    df['TMD'] = exp_tmd
    df['amino_acid'] = exp_aa
    df = df.drop(columns='feat_pos')  # Corrected line: remove feat_pos column
    df.to_csv(os.path.join(report_dir, 'translated_importance_report.csv'), index='Feature', mode="w")
    return df