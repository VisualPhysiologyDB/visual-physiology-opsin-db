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

    Args:
        df: The pandas DataFrame.
        lambda_max_column: The name of the column containing lambda max values.

    Returns:
        A new DataFrame with cleaned lambda max values.
    """

    new_rows = []

    for index, row in df.iterrows():
        lambda_max_str = str(row[lambda_max_column])  # Convert to string to handle potential NaN
        other_data = row.drop(lambda_max_column)

        if pd.isna(lambda_max_str) or lambda_max_str.lower() == 'nan' :  # Handle NaN or "nan" string values
            new_row = other_data.to_dict()
            new_row[lambda_max_column] = np.nan
            new_rows.append(new_row)
            continue


        if ',' in lambda_max_str:
            # Multiple values separated by commas
            values = lambda_max_str.split(',')
            for value in values:
                new_row = other_data.to_dict()
                # Remove non-numerical characters and take average if necessary
                cleaned_value = clean_single_value(value)
                new_row[lambda_max_column] = cleaned_value
                new_rows.append(new_row)

        elif ('-' in lambda_max_str) or ('–' in lambda_max_str):
            # Range of values
            new_row = other_data.to_dict()
            # Remove non-numerical characters and take average
            new_row[lambda_max_column] = clean_single_value(lambda_max_str)
            new_rows.append(new_row)
        else:
            # Single value (potentially with non-numerical characters)
            new_row = other_data.to_dict()
            #new_row[lambda_max_column] = clean_single_value(lambda_max_str)
            new_row[lambda_max_column] = float(lambda_max_str)
            new_rows.append(new_row)

    new_df = pd.DataFrame(new_rows)
    return new_df

def clean_single_value(value_str):
    """
    Cleans a single lambda max value string.
    
    Args:
        value_str: string of lambda max value

    Returns:
        Cleaned numerical lambda max value or np.nan if unable to extract
    """
    
    # Extract numbers and hyphens
    numbers = re.findall(r"[\d\-]+", value_str)
    
    if not numbers:
        print(f'Not a number: {value_str}')
        return np.nan  # Return NaN if no numbers are found
    
    # Handle multiple numbers
    if len(numbers) > 1:
        extracted_numbers = []
        for num in numbers:
            # If hyphen is present, split the values and convert to numbers
            if ('-' in num):
                try:
                    start, end = map(float, num.split('-'))
                    extracted_numbers.extend([start, end])
                except ValueError:
                    print(f'Not a number: {num}')
                    return np.nan  # Return NaN if there's a problem with splitting or converting to float
            elif ('–' in num):
                try:
                    start, end = map(float, num.split('–'))
                    extracted_numbers.extend([start, end])
                except ValueError:
                    print(f'Not a number: {num}')
                    return np.nan  # Return NaN if there's a problem with splitting or converting to float
            else:
                try:
                    extracted_numbers.append(float(num))
                except ValueError:
                    print(f'Not a number: {num}')
                    return np.nan # Return NaN if value cannot be converted to float
        return np.mean(extracted_numbers)
    
    # Handle single range or single number
    elif len(numbers) == 1:
        if ('-' in numbers[0]):
            try:
                start, end = map(float, numbers[0].split('-'))
                return np.mean([start, end])
            except ValueError:
                print(f'Not a number: {num}')
                return np.nan
        elif ('–' in numbers[0]):
            try:
                start, end = map(float, numbers[0].split('–'))
                return np.mean([start, end])
            except ValueError:
                print(f'Not a number: {num}')
                return np.nan     
        else:
            try:
                return float(numbers[0])
            except ValueError:
                print(f'Not a number: {num}')
                return np.nan

def merge_accessory_dbs(df_list, report_dir):
    """
    Merges a list of DataFrames, retaining specified columns and using index 
    columns as foreign keys.

    Args:
      df_list: A list of pandas DataFrames.

    Returns:
      A merged pandas DataFrame.
    """
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    merged_df = pd.DataFrame()  # Initialize an empty DataFrame
    
    for df in df_list:
      # If 'Accession' is missing, create it with NaN values
      if 'Accession' not in df.columns:
          df['Accession'] = pd.NA  

      # Extract the name of the index column
      index_name = df.index.name 
      # Select desired columns, including the index as a new column
      temp_df = df[['Full_Species', 'LambdaMax', 'Accession']].copy()
      temp_df = temp_df.reset_index(drop=True)
      temp_df[index_name] = df.index.astype(str).to_list()

      if merged_df.empty:
          merged_df = temp_df
      else:
          merged_df = pd.concat([merged_df, temp_df], ignore_index=True) 

    #merged_df.drop_duplicates(subset=['Full_Species', 'LambdaMax'], keep='first', inplace=True)
    # Sort the DataFrame by 'Accession' (descending) and then 'obs_id' (ascending)
    processed_df = merged_df.sort_values(by=['Accession'], ascending=[False])
    
    
    # Ensure the 'Accession' column is a string type to avoid errors with non-string data
    # We use .fillna('') to handle NaNs, then replace empty strings back to NA later.
    accession_col = processed_df['Accession'].fillna('').astype(str)

    # 1. Use regex to normalize various dash formats and encoding errors to a standard hyphen '-'
    # This pattern looks for 'â€“', '–' (en-dash), or ',Â'
    accession_col = accession_col.str.replace(r'â€“|–|,Â', '-', regex=True)

    # 2. Use regex to remove specific unwanted text like '(EST)' and all whitespace characters (\s)
    accession_col = accession_col.str.replace(r'\(EST\)|\s', '', regex=True)

    # Assign the cleaned series back to the DataFrame
    processed_df['Accession'] = accession_col.replace({'': pd.NA}) # Put NAs back
    
    #processed_df.reset_index(drop=True, inplace=True)
    processed_df.index.name = 'comp_db_id'
    #processed_df['LambdaMax'] = processed_df['LambdaMax'].str.replace('–','-').replace(' ','')
    processed_df.to_csv(f'{report_dir}/vpod_comp_acc_dbs_{dt_label}.csv', index=True)
    processed_df.dropna(subset=['LambdaMax'], inplace=True)
    
    #Clean the proccessed df column containing lambda max values of rows with mutliple entries or ranges
    cleaned_df = clean_lambda_max(processed_df.copy(), 'LambdaMax')  # Use df.copy() to avoid modifying the original DataFrame
    cleaned_df.dropna(subset=['LambdaMax'], inplace=True)
    # Drop duplicate rows based on 'Full_Species' and 'LambdaMax', keeping the first occurrence
    cleaned_df = cleaned_df.drop_duplicates(subset=['Full_Species', 'LambdaMax'], keep='first')
    cleaned_df.reset_index(drop=True, inplace=True)
    cleaned_df.index.name = 'comp_db_id'

    #save the final clean, merged df
    cleaned_df_name = f'{report_dir}/VPOD_in_vivo_1.0_{dt_label}.csv'
    cleaned_df.to_csv(cleaned_df_name, index=True)
    return cleaned_df, cleaned_df_name


def create_matched_df(report_dir, mnm_merged_df, source_data, prediction_to_use):
    
    """Creates a matched DataFrame by comparing predicted and measured values for opsin sequences.

    This function takes a DataFrame of merged NCBI and OPTICS prediction data (`mnm_merged_df`)
    and a source DataFrame (`source_data`) containing measured LambdaMax values for opsin sequences.
    It iterates through unique species in the `mnm_merged_df`, finds the closest measured LambdaMax 
    value for each predicted value, and creates a new DataFrame (`matched_df`) containing the 
    matched results.

    The matching process involves:
    1. Identifying unique species in the merged prediction data.
    2. Iterating through each species and filtering the prediction and source data accordingly.
    3. For each prediction, calculating the absolute difference between the predicted value 
       and all measured LambdaMax values for that species.
    4. Selecting the measurement with the minimum absolute difference as the closest match.
    5. Handling cases where multiple predictions for the same species and Accession might exist,
       keeping only the prediction with the smallest absolute difference to the measured value.
    6. Adding relevant information from the `mnm_merged_df` to the `matched_df`, including
       %Identity_Nearest_VPOD_Sequence, Gene_Description, Protein sequence, Genus, Species, 
       Phylum, Subphylum, and Class.
    7. Filtering out duplicate entries based on 'comp_db_id' and 'Full_Species', prioritizing 
       entries with the lowest 'abs_diff' and highest '%Identity_Nearest_VPOD_Sequence' when duplicates are found.

    Args:
        report_dir (str): Path to the directory where output files will be saved. (Used for saving intermediate output)
        mnm_merged_df (pandas.DataFrame): DataFrame containing merged NCBI and OPTICS prediction data. 
                                          Must have columns 'Full_Species', 'Accession', 'Prediction_Medians', and others that will be copied to matched_df.
        source_data (pandas.DataFrame): DataFrame containing measured LambdaMax values. Must have 
                                       columns 'Full_Species' and 'LambdaMax'. Also a 'comp_db_id' column that will be used as index later.

    Returns:
        pandas.DataFrame: A DataFrame (`matched_df`) containing the matched prediction and 
                          measurement data, along with additional information from `mnm_merged_df`.
    """    
    # Get unique species from predictions
    unique_species = list(mnm_merged_df['Full_Species'].unique())
    # Initialize a list to store the matched results
    matched_results = []
    i = 0
    # Iterate through each species
    for species in unique_species:
        scp_df_copy = source_data.copy()
        # Filter predictions and measurements for the current species
        species_predictions = mnm_merged_df[mnm_merged_df['Full_Species'] == species]
        species_measurements = scp_df_copy[scp_df_copy['Full_Species'] == species]
        if species_measurements.shape[0] == 0:
            i+=1
            print(f'\nThis species is missing a match: {species}\n')
        else:
            #print(species)
            # Iterate through each prediction for the current species
            for _, prediction_row in species_predictions.iterrows():
                prediction_value = prediction_row[prediction_to_use]
                try:
                    accession = prediction_row['Accession']
                except:
                    accession = prediction_row.name

                # Calculate absolute differences between the prediction and all measurements for the species
                species_measurements.loc[:, 'abs_diff'] = (species_measurements['LambdaMax'] - prediction_value).abs()
                # Find the closest measurement (handling ties)
                closest_measurement_row = species_measurements.sort_values('abs_diff',ascending=True).iloc[0]
                min_diff = closest_measurement_row['abs_diff']
                
                existing_match_index = next((i for i, match in enumerate(matched_results) if match['Accession'] == accession), None)

                if (existing_match_index is not None): 
                    if (min_diff < matched_results[existing_match_index]['abs_diff']):
                        # If it exists, compare the absolute differences and keep the better match
                        matched_results[existing_match_index].update({
                            #'Accession': accession,
                            #'Full_Species': species,
                            'prediction_value': prediction_value,
                            'LambdaMax': closest_measurement_row['LambdaMax'],
                            'abs_diff': min_diff,
                            'comp_db_id': closest_measurement_row.name
                        })
                        #print('Updating match dataframe, better match found.')
                else:
                    # If it doesn't exist, add the new match
                    matched_results.append({
                        'Accession': accession,
                        'Full_Species': species,
                        'prediction_value': prediction_value,
                        'LambdaMax': closest_measurement_row['LambdaMax'],
                        'abs_diff': min_diff,
                        'comp_db_id': closest_measurement_row.name
                    })
                
    # Create a new dataframe from the matched results
    matched_df = pd.DataFrame(matched_results)
    #print(f'This is the raw matched dataframe:\n{matched_df.head(5)}')
    print(f'There were {i} unmatched species')
    
    #Match %Idenitities to Accessions
    iden_list = []
    prot_des_list = []
    aa_seq_list = []
    genus_list = []
    species_list = []
    phylum_list = []
    subphylum_list = []
    class_list = []
    
    #mnm_merged_df.set_index(['Accession'], inplace=True)
    mnm_merged_df.set_index(['Accession','Full_Species'], inplace=True)

    for _, d in matched_df.iterrows():
        acc = d['Accession']
        genus_species = d['Full_Species']
        try:
            iden_list.append(mnm_merged_df.loc[acc,genus_species]['%Identity_Nearest_VPOD_Sequence'][0])
            prot_des_list.append(mnm_merged_df.loc[acc,genus_species]['Gene_Description'][0])
            aa_seq_list.append(mnm_merged_df.loc[acc,genus_species]['Protein'][0])
            genus_list.append(mnm_merged_df.loc[acc,genus_species]['Genus'][0])
            species_list.append(mnm_merged_df.loc[acc,genus_species]['Species'][0])
            try:
                phylum_list.append(mnm_merged_df.loc[acc,genus_species]['Phylum'][0])
                subphylum_list.append(mnm_merged_df.loc[acc,genus_species]['Subphylum'][0])
                class_list.append(mnm_merged_df.loc[acc,genus_species]['Class'][0])
            except:
                pass
        except:
            iden_list.append(mnm_merged_df.loc[acc,genus_species]['%Identity_Nearest_VPOD_Sequence'])
            prot_des_list.append(mnm_merged_df.loc[acc,genus_species]['Gene_Description'])
            aa_seq_list.append(mnm_merged_df.loc[acc,genus_species]['Protein'])
            genus_list.append(mnm_merged_df.loc[acc,genus_species]['Genus'])
            species_list.append(mnm_merged_df.loc[acc,genus_species]['Species'])
            try:
                phylum_list.append(mnm_merged_df.loc[acc,genus_species]['Phylum'])
                subphylum_list.append(mnm_merged_df.loc[acc,genus_species]['Subphylum'])
                class_list.append(mnm_merged_df.loc[acc,genus_species]['Class'])
            except:
                pass
        
    matched_df['%Identity_Nearest_VPOD_Sequence'] = iden_list
    matched_df['Gene_Description'] = prot_des_list
    matched_df['Protein'] = aa_seq_list
    matched_df['Genus'] = genus_list
    matched_df['Species'] = species_list
    matched_df['Phylum'] = phylum_list
    matched_df['Subphylum'] = subphylum_list
    matched_df['Class'] = class_list
    matched_df = matched_df.reindex(columns=['Accession','Phylum','Subphylum','Class','Genus','Species','Full_Species','%Identity_Nearest_VPOD_Sequence','prediction_value','LambdaMax','abs_diff','comp_db_id','Protein','Gene_Description','Notes'])
    matched_df.to_csv(f'./{report_dir}/uncleaned_matched_mnm_df.csv')
    
    #matched_df.sort_values(['abs_diff', 'Accession', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False, False])
    #matched_df.drop_duplicates('Accession', keep='first',inplace=True,ignore_index=True)
    #matched_df.reset_index(drop=True, inplace=True)
    #matched_df.index.name = 'mnm_id'
    #matched_df.to_csv('test_final.csv')
    
    # Group by 'comp_db_id' and count unique accessions
    grouped_counts = matched_df.groupby(['comp_db_id','Full_Species'])['Accession'].nunique()
    # Filter groups with more than one unique accession
    duplicate_comp_db_id_groups = grouped_counts[grouped_counts > 1]
    #print(duplicate_comp_db_id_groups)
    if not duplicate_comp_db_id_groups.empty:
        filtered_results = []
        # Iterate through each mnm_id with duplicates
        for comp_db_id,fs in duplicate_comp_db_id_groups.index:
            #print(f'{comp_db_id} , {fs}')
            # Filter rows with the current mnm_id
            duplicates = matched_df[(matched_df['comp_db_id'] == comp_db_id) & (matched_df['Full_Species'] == fs)]            #print(duplicates)
            # Sort by abs_diff (ascending) and then percent_identity (descending)
            try:
                duplicates = duplicates.sort_values(['abs_diff', 'Accession'], ascending=[True, False])
            except:
                print(f'Duplicates error: {duplicates}')
                raise Exception
            try:
                duplicates = duplicates.sort_values(['abs_diff', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False])
            except:
                print(f'Duplicates error: {duplicates}')
                raise Exception
            filtered_results.append(duplicates.iloc[0])

        # Combine filtered results with non-duplicate rows
        #comp_db_ids = duplicate_comp_db_id_groups.index.get_level_values('comp_db_id')
        #fsns = duplicate_comp_db_id_groups.index.get_level_values('Full_Species')
        non_duplicates = matched_df[~matched_df.set_index(['comp_db_id', 'Full_Species']).index.isin(duplicate_comp_db_id_groups.index)]        
        final_filtered_df = pd.concat([pd.DataFrame(filtered_results), non_duplicates], ignore_index=True)
        #final_filtered_df = final_filtered_df.sort_values(['abs_diff', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False])
        final_filtered_df.reset_index(drop=True, inplace=True)
        final_filtered_df.index.name = 'mnm_id'
    else:
        final_filtered_df = matched_df 
    
    return(final_filtered_df)


def mine_n_match(email, report_dir, source_file, ncbi_q_file, optics_pred_file, out='unnamed', err_filter=15, per_identity_minimum=20, prediction_to_use='Prediction_Medians'):
    
    """Mines and matches opsin sequences based on various criteria and prediction data.

    This function performs a multi-step process to identify and filter opsin sequences 
    based on provided data sources, including a source file with known opsin data, 
    NCBI query results, and prediction data from the OPTICS tool.

    The process involves:
    1. Loading and preprocessing data from the source file, NCBI query file, and OPTICS prediction file.
    2. Merging the NCBI and OPTICS data, filtering by identity to the nearest VPOD sequence, 
       and creating a matched dataframe based on source data.
    3. Filtering the matched dataframe by an error threshold (abs_diff) and removing entries 
       with 100% identity to the nearest VPOD sequence.
    4. Handling redundant sequences with different accession numbers by selecting the entry 
       with the lowest abs_diff and adding a note about the redundancy.
    5. Incorporating additional sequences from the source data that were not found in the 
       initial matching process.
    6. Retrieving and adding taxonomic information (Phylum, Subphylum, Class) for 
       the added sequences.
    7. Generating a final dataframe with filtered and annotated opsin data, including 
       comp_db_id, Accession, taxonomic information, LambdaMax, Protein, and Notes.

    Args:
        email (str): User's email address, used for NCBI Entrez queries.
        report_dir (str): Path to the directory where output files will be saved.
        source_file (str): Path to the source data file (CSV or TSV format). This file should contain a "Full_Species" column or "Genus" and "Species" columns, and an "Accession" column that will be filtered to include only entries without "-" or "–".
        ncbi_q_file (str): Path to the NCBI query results file (CSV format).
        optics_pred_file (str): Path to the OPTICS prediction results file (TSV format).
        out (str, optional): Base name for output files. Defaults to 'unnamed'.
        err_filter (int, optional): Error threshold (abs_diff) for filtering results. Defaults to 15.

    Returns:
        pandas.DataFrame: A DataFrame containing the final filtered and annotated opsin data. 
        The DataFrame is also saved as a CSV file in the report_dir.
    """
    try:
        source_data = pd.read_csv(source_file, index_col = 0)
    except:
        source_file = f'./{report_dir}/{source_file}'
        source_data = pd.read_csv(source_file, index_col = 0)   
        
    if 'Genus' in source_data.columns and 'Species' in source_data.columns and 'Full_Species' not in source_data.columns:
        source_data['Full_Species'] = source_data['Genus'] + '_' + source_data['Species']
    elif 'Genus' not in source_data.columns and 'Species' not in source_data.columns and 'Full_Species' in source_data.columns:
        gn_list = []
        sp_list = []
        og_sp_list = source_data['Full_Species'].to_list()
        for sp in og_sp_list:
            try:
            #    print(sp)
                genus = sp.split(' ', 1)[0]
                species = sp.split(' ', 1)[1]
                gn_list.append(genus)
                sp_list.append(species)
            except:
                try:
                    genus = sp.split('Â', 1)[0]
                    species = sp.split('Â', 1)[1]
                    gn_list.append(genus)
                    sp_list.append(species)
                except:
                    try:
                        genus = sp.split('Â', 1)[0]
                        species = sp.split('Â', 1)[0]
                        gn_list.append(genus)
                        sp_list.append(species)
                    except:
                        raise Exception(f'Problematic species name: {sp}')
        source_data['Genus'] = gn_list
        source_data['Species'] = sp_list

    if isinstance(ncbi_q_file, pd.DataFrame):
        ncbi_data = ncbi_q_file.copy()
    else:
        try:    
            ncbi_data = pd.read_csv(ncbi_q_file)
        except:
            ncbi_q_file = f'./{report_dir}/{ncbi_q_file}'
            ncbi_data = pd.read_csv(ncbi_q_file)
    try:
        pred_data = pd.read_csv(optics_pred_file, sep = '\t')
    except:
        optics_pred_file = f'./{report_dir}/{optics_pred_file}'
        pred_data = pd.read_csv(optics_pred_file, sep = '\t')
        
    mnm_merged_df = pd.merge(ncbi_data, pred_data, left_index=True, right_index=True)
    mnm_merged_df.drop_duplicates(subset=['Full_Species', 'Protein'],  keep='first', inplace=True)
    mnm_merged_df = mnm_merged_df[(mnm_merged_df['%Identity_Nearest_VPOD_Sequence'] != 'blastp unsuccessful') & (mnm_merged_df['%Identity_Nearest_VPOD_Sequence']>per_identity_minimum)]
    mnm_merged_df.reset_index(inplace=True, drop=True)
    #mnm_merged_df.index.name = 'comp_db_id'
    mnm_merged_df.to_csv(path_or_buf=f"./{report_dir}/mnm_merged_df.csv", index=True)

    matched_df = create_matched_df(report_dir=report_dir, mnm_merged_df=mnm_merged_df, source_data=source_data, prediction_to_use=prediction_to_use)
    matched_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_results_id_filtered.csv", index=True)
    #note that prediction values from optics are taken from the mean of the bootstrap predictions
    final_err_filtered_df = matched_df[matched_df['abs_diff'] <= err_filter]
    final_err_filtered_df = final_err_filtered_df[final_err_filtered_df['%Identity_Nearest_VPOD_Sequence'] != 100]
    final_err_filtered_df.reset_index(inplace=True, drop=True)
    final_err_filtered_df.index.name = 'mnm_id'
    final_err_filtered_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_results_err_filtered.csv", index=True)
    
    # Initialize a dictionary with keys as the `Protein` column value and values as a list of the `Accession` column value
    protein_accession_dict = final_err_filtered_df.groupby("Protein")['Accession'].apply(list).to_dict()
    # Initialize a list to store the results
    result_list = []
    # Iterate through all the proteins in the dataset
    for protein, accessions in protein_accession_dict.items():
        # Filter the dataframe for rows where the protein is the same
        temp_df = final_err_filtered_df[final_err_filtered_df["Protein"] == protein].reset_index(drop = True)

        # Check if there are more than 1 unique accession number for the same protein
        if len(set(accessions)) > 1:
            # Sort the dataframe by `abs_diff`
            temp_df = temp_df.sort_values("abs_diff")

            # Take the first row
            first_row = temp_df.iloc[[0]]

            # Take the `comp_db_id` of the other rows
            comp_db_id_list = temp_df.iloc[1:]['comp_db_id'].to_list()

            # If there are multiple `comp_db_id`, join them together
            if len(comp_db_id_list) > 1:
                comp_db_id_str = ", ".join(str(comp_db_id) for comp_db_id in comp_db_id_list)
            else:
                comp_db_id_str = str(comp_db_id_list[0])

            # Add the `comp_db_id` to the `Notes` column
            first_row["Notes"] = "Entry with different Accession from different species but redundant sequence found. Comp_db_id: " + comp_db_id_str

            # Append the first row to the result list
            result_list.append(first_row)
        else:
            # Append the row to the result list
            result_list.append(temp_df)

    # Concatenate all the series into a single dataframe
    redundant_filter_df = pd.concat(result_list).reset_index(drop = True)
    redundant_filter_df.index.name = 'mnm_id'
    
    # Add in any sequences with accessions from the lmax compendium (aka - 'source_data')
    sd_2 = source_data.copy()
    sd_2= sd_2[~sd_2["Accession"].isna()].reset_index()
    
    # Filter the Accession column to include only those without '-' or '–'
    sd_2 = sd_2[~sd_2["Accession"].str.contains("-|–|â€“|,")].reset_index(drop=True)
    #sd_2 = sd_2[~sd_2["Accession"].str.contains("-|–|")].reset_index(drop = True)

    # Create a list of `Accession` in `result_df`
    accession_list = final_err_filtered_df["Accession"].to_list()

    # Filter `sd_2` to include only those whose `Accession` is not in `accession_list`
    sd_2_filtered = sd_2[~sd_2["Accession"].isin(accession_list)].reset_index(drop = True)

    prot_list = get_prots_from_acc(sd_2_filtered['Accession'].to_list())
    sd_2_filtered['Protein'] = prot_list

    # Split the `Full_Species` column by the first space (' ') into `Genus` and `Species`
    g_list = [x.split(' ')[0] for x in sd_2_filtered['Full_Species']]
    sp_list = [' '.join(x.split(' ')[1:]) for x in sd_2_filtered['Full_Species']]
    sp_list = [' '.join(x.split(' ')[0:]) for x in sp_list]

    sd_2_filtered["Genus"] = g_list
    sd_2_filtered["Species"] = sp_list

    # Filter to make sure that no other sequences here are repeats of sequences in the main dataframe
    existing_proteins = redundant_filter_df["Protein"].to_list()
    sd_2_filtered = sd_2_filtered[~sd_2_filtered["Protein"].isin(existing_proteins)].reset_index(drop = True)
    
    ## Filter to make sure that no other sequences here are repeats of sequences in the main dataframe
    #existing_proteins = sd_2_filtered["Protein"].to_list()
    #existing_ids = sd_2_filtered["comp_db_id"].to_list()
    #sd_2_filtered = sd_2_filtered[~sd_2_filtered["Protein"].isin(existing_proteins)].reset_index(drop = True)
    
    #redundant_filter_df = redundant_filter_df[(~redundant_filter_df["Protein"].isin(existing_proteins)) & (~redundant_filter_df["comp_db_id"].isin(existing_ids))].reset_index(drop = True)
    #redundant_filter_df.index.name = 'mnm_id'

    # Need to get the taxon data for the added data
    # Call taxa dictionary function here
    species_list = sd_2_filtered['Full_Species'].to_list()
    
    taxon_file = './data_sources/taxonomy/ncbi_taxon_dict.json'
    if os.path.isfile(taxon_file):
        try:
            with open(taxon_file, 'r') as f:
                existing_taxon_dict = json.load(f)
        except FileNotFoundError:
            print(f"Error: File '{taxon_file}' not found or can't be loaded\n")
        print('Existing Taxon Dictionary Found! One Moment While We Update It...\n')
        species_taxon_dict = get_sp_taxon_dict(species_list = species_list, email = email, taxon_file = taxon_file, sp_taxon_dict = copy.deepcopy(existing_taxon_dict))
    else:
        existing_taxon_dict = {}
        species_taxon_dict = get_sp_taxon_dict(species_list = species_list, email = email, taxon_file = taxon_file)
    
    # Save the taxon dictionary if it doesn't yet exist or if it has been updated since being loaded 
    if (list(species_taxon_dict.keys()) >= list(existing_taxon_dict.keys())):  
        #print('Saving Updated Dictionary') 
        try:
            with open(taxon_file, 'w') as f:
                json.dump(species_taxon_dict, f, indent=4)  # indent for pretty formatting
        except FileNotFoundError:
            print(f"Error: File '{taxon_file}' can't be saved...\n")

    Phylum = []
    Subphylum = []
    Class = []
    for sp in species_list:    
        Phylum.append(species_taxon_dict[sp]["Phylum"])
        Subphylum.append(species_taxon_dict[sp]["Subphylum"])
        Class.append(species_taxon_dict[sp]["Class"])
    sd_2_filtered["Phylum"] = Phylum
    sd_2_filtered['Subphylum'] = Subphylum
    sd_2_filtered['Class'] = Class

    # Take the `comp_db_id`, `Genus`, `Species`, `Accession`, and `LambdaMax` from the filtered dataframe and rename the `LambdaMax` column to `closest_measurement`
    sd_2_filtered = sd_2_filtered[["comp_db_id", "Accession", "Phylum", "Subphylum", "Class", "Genus", "Species",  "LambdaMax", "Protein"]]

    # Add a column named `Notes` with the text 'Known sequence specified in accessory database'
    sd_2_filtered["Notes"] = "Known sequence specified in accessory database"

    # Concatenate this series with the `result_df` dataframe
    final_df = pd.concat([redundant_filter_df, sd_2_filtered]).reset_index(drop = True)
    final_df.drop_duplicates(subset=['comp_db_id'],  keep='first', inplace=True)
    final_df = final_df.reset_index(drop = True)
    final_df.index.name = 'mnm_id'
    final_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_results_fully_filtered.csv", index=True)

    return(final_df)
    


