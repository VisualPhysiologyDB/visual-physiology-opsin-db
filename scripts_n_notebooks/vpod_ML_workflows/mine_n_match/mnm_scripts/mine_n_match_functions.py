import os
import re
import datetime
import time
import json
import copy
import pandas as pd
import numpy as np
import io
from Bio import Entrez, SeqIO
from pygbif import species
from tqdm import tqdm

# Our custom modules
from mnm_scripts.ncbi_functions import get_prots_from_acc
from mnm_scripts.taxonomy import get_sp_taxon_dict

# API key for NCBI 
api_key = "1efb120056e1cea873ba8d85d6692abd5d09"


def clean_lambda_max(df, lambda_max_column):
    """
    Cleans a DataFrame column containing lambda max values.
    """
    new_rows = []

    for index, row in df.iterrows():
        lambda_max_str = str(row[lambda_max_column])
        other_data = row.drop(lambda_max_column)

        if pd.isna(lambda_max_str) or lambda_max_str.lower() == 'nan':
            new_row = other_data.to_dict()
            new_row[lambda_max_column] = np.nan
            new_rows.append(new_row)
            continue

        if ',' in lambda_max_str:
            values = lambda_max_str.split(',')
            for value in values:
                new_row = other_data.to_dict()
                new_row[lambda_max_column] = clean_single_value(value)
                new_rows.append(new_row)
        elif '-' in lambda_max_str or '–' in lambda_max_str:
            new_row = other_data.to_dict()
            new_row[lambda_max_column] = clean_single_value(lambda_max_str)
            new_rows.append(new_row)
        else:
            new_row = other_data.to_dict()
            new_row[lambda_max_column] = float(lambda_max_str)
            new_rows.append(new_row)

    return pd.DataFrame(new_rows)


def clean_single_value(value_str):
    """
    Cleans a single lambda max value string by averaging ranges or returning floats.
    """
    numbers = re.findall(r"[\d\-]+", value_str)
    
    if not numbers:
        print(f'Not a number: {value_str}')
        return np.nan
    
    extracted_numbers = []
    for num in numbers:
        if '-' in num or '–' in num:
            separator = '-' if '-' in num else '–'
            try:
                start, end = map(float, num.split(separator))
                extracted_numbers.extend([start, end])
            except ValueError:
                print(f'Not a number: {num}')
                return np.nan
        else:
            try:
                extracted_numbers.append(float(num))
            except ValueError:
                print(f'Not a number: {num}')
                return np.nan
                
    return np.mean(extracted_numbers) if extracted_numbers else np.nan


def merge_accessory_dbs(df_list, report_dir):
    """
    Merges a list of DataFrames, retaining specified columns and using index columns as foreign keys.
    """
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    merged_df = pd.DataFrame()
    
    for df in df_list:
        if 'Accession' not in df.columns:
            df['Accession'] = pd.NA  

        index_name = df.index.name 
        temp_df = df[['Full_Species', 'LambdaMax', 'Accession']].copy()
        temp_df = temp_df.reset_index(drop=True)
        temp_df[index_name] = df.index.astype(str).to_list()

        merged_df = pd.concat([merged_df, temp_df], ignore_index=True) if not merged_df.empty else temp_df

    processed_df = merged_df.sort_values(by=['Accession'], ascending=[False])
    
    accession_col = processed_df['Accession'].fillna('').astype(str)
    accession_col = accession_col.str.replace(r'â€“|–|,Â', '-', regex=True)
    accession_col = accession_col.str.replace(r'\(EST\)|\s', '', regex=True)
    processed_df['Accession'] = accession_col.replace({'': pd.NA}) 
    
    processed_df.index.name = 'comp_db_id'
    processed_df.to_csv(f'{report_dir}/vpod_comp_acc_dbs_{dt_label}.csv', index=True)
    processed_df.dropna(subset=['LambdaMax'], inplace=True)
    
    cleaned_df = clean_lambda_max(processed_df.copy(), 'LambdaMax') 
    cleaned_df.dropna(subset=['LambdaMax'], inplace=True)
    cleaned_df = cleaned_df.drop_duplicates(subset=['Full_Species', 'LambdaMax'], keep='first')
    cleaned_df.reset_index(drop=True, inplace=True)
    cleaned_df.index.name = 'comp_db_id'

    cleaned_df_name = f'{report_dir}/VPOD_in_vivo_1.0_{dt_label}.csv'
    cleaned_df.to_csv(cleaned_df_name, index=True)
    return cleaned_df, cleaned_df_name


def create_matched_df(report_dir, mnm_merged_df, source_data, prediction_to_use, err_filter):
    """
    Creates a matched DataFrame mapping predicted values to measured values via Greedy Bipartite Matching.
    Calculates candidate counts and allows fallback mappings within the error threshold.
    """
    unique_species = mnm_merged_df['Full_Species'].unique()
    matched_results = []
    candidate_counts = {}

    unmatched_species_count = 0

    for species in unique_species:
        species_preds = mnm_merged_df[mnm_merged_df['Full_Species'] == species]
        species_meas = source_data[source_data['Full_Species'] == species]

        if species_meas.empty:
            unmatched_species_count += 1
            print(f'\nThis species is missing a match: {species}\n')
            continue

        pairs = []
        
        # 1. Create all possible Sequence <-> Measurement pairings
        for p_idx, p_row in species_preds.iterrows():
            p_val = p_row[prediction_to_use]
            acc = p_row.get('Accession', p_idx)
            
            # Extract identity for tie breaking
            identity = p_row.get('%Identity_Nearest_VPOD_Sequence', 0)
            if isinstance(identity, list) and len(identity) > 0: identity = identity[0]
            try: identity = float(identity)
            except (ValueError, TypeError): identity = 0.0
            
            for m_idx, m_row in species_meas.iterrows():
                m_val = m_row['LambdaMax']
                diff = abs(p_val - m_val)
                
                # comp_db_id is likely the index in source_data
                comp_id = m_row.name if 'comp_db_id' not in m_row else m_row['comp_db_id']
                
                pairs.append({
                    'Accession': acc,
                    'comp_db_id': comp_id,
                    'prediction_value': p_val,
                    'LambdaMax': m_val,
                    'abs_diff': diff,
                    'Full_Species': species,
                    'Identity_Tiebreaker': identity
                })

        # 2. Count candidates that fall within the err_filter threshold for each measurement
        for m_idx, m_row in species_meas.iterrows():
            comp_id = m_row.name if 'comp_db_id' not in m_row else m_row['comp_db_id']
            valid_preds = [p for p in pairs if p['comp_db_id'] == comp_id and p['abs_diff'] <= err_filter]
            candidate_counts[comp_id] = len(valid_preds)

        # 3. Greedy Bipartite Assignment
        # Sort best matches first (smallest distance, highest identity)
        pairs.sort(key=lambda x: (x['abs_diff'], -x['Identity_Tiebreaker']))

        assigned_accessions = set()
        assigned_measurements = set()

        for pair in pairs:
            # Only assign if neither the Accession nor the Measurement has been taken
            if pair['Accession'] not in assigned_accessions and pair['comp_db_id'] not in assigned_measurements:
                assigned_accessions.add(pair['Accession'])
                assigned_measurements.add(pair['comp_db_id'])
                matched_results.append(pair)

    print(f'There were {unmatched_species_count} unmatched species')
    matched_df = pd.DataFrame(matched_results)
    if matched_df.empty: return matched_df

    # 4. Clean Data Merging (Replacing iterative lookups)
    matched_df['candidate_opsin_count'] = matched_df['comp_db_id'].map(candidate_counts)
    
    cols_to_bring = ['%Identity_Nearest_VPOD_Sequence', 'Gene_Description', 'Protein', 'Genus', 'Species', 'Phylum', 'Subphylum', 'Class']
    available_cols = [c for c in cols_to_bring if c in mnm_merged_df.columns]
    
    mnm_subset = mnm_merged_df[['Accession', 'Full_Species'] + available_cols].drop_duplicates(subset=['Accession', 'Full_Species'])
    final_matched_df = pd.merge(matched_df, mnm_subset, on=['Accession', 'Full_Species'], how='left')

    # Drop temporary tiebreaker
    if 'Identity_Tiebreaker' in final_matched_df.columns:
        final_matched_df.drop(columns=['Identity_Tiebreaker'], inplace=True)

    # Standardize column order
    ordered_cols = ['Accession', 'Phylum', 'Subphylum', 'Class', 'Genus', 'Species', 'Full_Species', 
                    '%Identity_Nearest_VPOD_Sequence', 'prediction_value', 'LambdaMax', 'abs_diff', 
                    'candidate_opsin_count', 'comp_db_id', 'Protein', 'Gene_Description', 'Notes']
    ordered_cols = [c for c in ordered_cols if c in final_matched_df.columns]
    final_matched_df = final_matched_df.reindex(columns=ordered_cols)

    final_matched_df.to_csv(f'./{report_dir}/uncleaned_matched_mnm_df.csv', index=False)
    return final_matched_df


# --- Refactored Pipeline Helpers ---

def _prepare_source_data(source_file, report_dir):
    """Loads and standardizes naming conventions for source accessory data."""
    path = source_file if os.path.isfile(source_file) else f'./{report_dir}/{source_file}'
    df = pd.read_csv(path, index_col=0)
    
    if 'Genus' in df.columns and 'Species' in df.columns and 'Full_Species' not in df.columns:
        df['Full_Species'] = df['Genus'] + '_' + df['Species']
    elif 'Full_Species' in df.columns and ('Genus' not in df.columns or 'Species' not in df.columns):
        # Fix encoding errors dynamically
        clean_sp = df['Full_Species'].str.replace('Â', ' ', regex=False).str.strip()
        splits = clean_sp.str.split(n=1, expand=True)
        df['Genus'] = splits[0]
        df['Species'] = splits[1].fillna(splits[0])
        
    return df

def _load_data(file_val, report_dir, sep=','):
    """Loads data standardizing paths between local and report_dir."""
    if isinstance(file_val, pd.DataFrame): return file_val.copy()
    path = file_val if os.path.isfile(file_val) else f'./{report_dir}/{file_val}'
    return pd.read_csv(path, sep=sep)

def _handle_redundant_proteins(df):
    """Flags items resolving to the exact same sequence but holding different Accessions."""
    result_list = []
    
    for protein, group in df.groupby('Protein'):
        if group['Accession'].nunique() > 1:
            group = group.sort_values('abs_diff')
            first_row = group.iloc[[0]].copy()
            
            comp_db_id_list = group.iloc[1:]['comp_db_id'].tolist()
            comp_db_id_str = ", ".join(map(str, comp_db_id_list))
            first_row['Notes'] = "Entry with different Accession from different species but redundant sequence found. Comp_db_id: " + comp_db_id_str
            
            result_list.append(first_row)
        else:
            result_list.append(group)
            
    return pd.concat(result_list).reset_index(drop=True)

def _append_accessory_data(filtered_df, source_data, email, report_dir):
    """Appends un-matched known sequences from the accessory database with taxonomy processing."""
    sd_2 = source_data.copy()
    if "Accession" not in sd_2.columns: return filtered_df.copy()
    
    sd_2 = sd_2.dropna(subset=["Accession"]).reset_index()
    if sd_2.empty: return filtered_df.copy()
    
    # Filter out dashed Accessions & entries already existing in our results
    sd_2 = sd_2[~sd_2["Accession"].str.contains(r"-|–|â€“|,", regex=True)].reset_index(drop=True)
    sd_2_filtered = sd_2[~sd_2["Accession"].isin(filtered_df["Accession"])].reset_index(drop=True)
    if sd_2_filtered.empty: return filtered_df.copy()

    sd_2_filtered['Protein'] = get_prots_from_acc(sd_2_filtered['Accession'].tolist())
    sd_2_filtered['Genus'] = sd_2_filtered['Full_Species'].apply(lambda x: str(x).split(' ')[0])
    sd_2_filtered['Species'] = sd_2_filtered['Full_Species'].apply(lambda x: ' '.join(str(x).split(' ')[1:]))

    # Protect against sequence redundancy
    sd_2_filtered = sd_2_filtered[~sd_2_filtered["Protein"].isin(filtered_df["Protein"])].reset_index(drop=True)
    if sd_2_filtered.empty: return filtered_df.copy()

    # Taxonomy fetching
    taxon_file = './data_sources/taxonomy/ncbi_taxon_dict.json'
    existing_taxon_dict = {}
    
    if os.path.isfile(taxon_file):
        try:
            with open(taxon_file, 'r') as f: existing_taxon_dict = json.load(f)
            print('Existing Taxon Dictionary Found! One Moment While We Update It...\n')
        except Exception as e:
            print(f"Error reading taxon file: {e}")

    species_list = sd_2_filtered['Full_Species'].tolist()
    species_taxon_dict = get_sp_taxon_dict(species_list=species_list, email=email, taxon_file=taxon_file, sp_taxon_dict=copy.deepcopy(existing_taxon_dict))
    
    if len(species_taxon_dict) >= len(existing_taxon_dict):
        os.makedirs(os.path.dirname(taxon_file), exist_ok=True)
        try:
            with open(taxon_file, 'w') as f: json.dump(species_taxon_dict, f, indent=4)
        except Exception as e:
            print(f"Error saving taxon file: {e}")

    # Map taxons
    sd_2_filtered["Phylum"] = sd_2_filtered['Full_Species'].map(lambda x: species_taxon_dict.get(x, {}).get("Phylum"))
    sd_2_filtered["Subphylum"] = sd_2_filtered['Full_Species'].map(lambda x: species_taxon_dict.get(x, {}).get("Subphylum"))
    sd_2_filtered["Class"] = sd_2_filtered['Full_Species'].map(lambda x: species_taxon_dict.get(x, {}).get("Class"))

    # If 'index' column exists from reset_index() but should be comp_db_id:
    if 'index' in sd_2_filtered.columns and 'comp_db_id' not in sd_2_filtered.columns:
        sd_2_filtered = sd_2_filtered.rename(columns={'index': 'comp_db_id'})

    sd_2_filtered["Notes"] = "Known sequence specified in accessory database"
    sd_2_filtered["candidate_opsin_count"] = 1 # By default, exactly 1 sequence mapped for known ones.

    cols_to_keep = ["comp_db_id", "Accession", "Phylum", "Subphylum", "Class", "Genus", "Species", "LambdaMax", "Protein", "Notes", "candidate_opsin_count"]
    cols_to_keep = [c for c in cols_to_keep if c in sd_2_filtered.columns]
    sd_2_filtered = sd_2_filtered[cols_to_keep]

    return pd.concat([filtered_df, sd_2_filtered]).reset_index(drop=True)


def mine_n_match(email, report_dir, source_file, ncbi_q_file, optics_pred_file, out='unnamed', err_filter=15, per_identity_minimum=20, prediction_to_use='Prediction_Medians', identity_filter=True):
    """
    Mines and matches opsin sequences based on prediction logic.
    """
    os.makedirs(report_dir, exist_ok=True)
    
    # 1. Prepare base data
    source_data = _prepare_source_data(source_file, report_dir)
    ncbi_data = _load_data(ncbi_q_file, report_dir)
    pred_data = _load_data(optics_pred_file, report_dir, sep='\t')
        
    # 2. Merge and clean NCBI + Predictions
    mnm_merged_df = pd.merge(ncbi_data, pred_data, left_index=True, right_index=True)
    mnm_merged_df.drop_duplicates(subset=['Full_Species', 'Protein'], keep='first', inplace=True)
    
    mnm_merged_df = mnm_merged_df[mnm_merged_df['%Identity_Nearest_VPOD_Sequence'] != 'blastp unsuccessful'].copy()
    mnm_merged_df['%Identity_Nearest_VPOD_Sequence'] = pd.to_numeric(mnm_merged_df['%Identity_Nearest_VPOD_Sequence'], errors='coerce')
    mnm_merged_df = mnm_merged_df[mnm_merged_df['%Identity_Nearest_VPOD_Sequence'] > per_identity_minimum]
    mnm_merged_df.reset_index(inplace=True, drop=True)
    mnm_merged_df.to_csv(f"./{report_dir}/mnm_merged_df.csv", index=True)

    # 3. Create matched df (Now using the Greedy assignment with candidate counts)
    matched_df = create_matched_df(report_dir, mnm_merged_df, source_data, prediction_to_use, err_filter)
    if matched_df.empty:
        print("No matches could be found. Returning empty dataframe.")
        return matched_df

    matched_df.index.name = 'mnm_id'
    matched_df.to_csv(f"./{report_dir}/mnm_on_{out}_results_id_filtered.csv", index=True)
    
    # 4. Strict Error filtering logic
    final_err_filtered_df = matched_df[matched_df['abs_diff'] <= err_filter].copy()
    if identity_filter:
        final_err_filtered_df = final_err_filtered_df[final_err_filtered_df['%Identity_Nearest_VPOD_Sequence'] != 100]
        
    final_err_filtered_df.reset_index(inplace=True, drop=True)
    final_err_filtered_df.index.name = 'mnm_id'
    final_err_filtered_df.to_csv(f"./{report_dir}/mnm_on_{out}_results_err_filtered.csv", index=True)
    
    # 5. Pipeline modifications to remove redundancies and append standard accessories
    redundant_filter_df = _handle_redundant_proteins(final_err_filtered_df)
    final_df = _append_accessory_data(redundant_filter_df, source_data, email, report_dir)
            
    # 6. Final Clean & Save
    final_df.drop_duplicates(subset=['comp_db_id'], keep='first', inplace=True)
    final_df = final_df.reset_index(drop=True)
    final_df.index.name = 'mnm_id'
    final_df.to_csv(f"./{report_dir}/mnm_on_{out}_results_fully_filtered.csv", index=True)

    return final_df


