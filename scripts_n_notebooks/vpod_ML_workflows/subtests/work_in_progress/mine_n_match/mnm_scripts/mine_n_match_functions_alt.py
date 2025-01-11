import os
import re
import datetime
import time
import pandas as pd
from Bio import Entrez, SeqIO
import progressbar


# create a function that return the collection of all NCBI fetch result
def ncbi_fetch(email, term, ncbi_db = "nuccore", rettype = "gb", format = "genbank"):
  Entrez.email = email    # Always tell NCBI who you are
  handle = Entrez.esearch(db=ncbi_db,
                        term=term, 
                        api_key = "1efb120056e1cea873ba8d85d6692abd5d09") # using api key allows 10 fetch per second
  record = Entrez.read(handle)
  #print(record)
  full_res = int(record["Count"])
  
  handle_full = Entrez.esearch(db=ncbi_db,
                        term=term, 
                        api_key = "1efb120056e1cea873ba8d85d6692abd5d09",
                        retmax = full_res + 1)
  record_full = Entrez.read(handle_full)
  #print(record_full)
  q_list = record_full["IdList"]

# create a list for all the sequence fetched from NCBI
  NCBI_seq =[]
  for x in q_list:
    # fetch result from previous search
    fet = Entrez.efetch(db=ncbi_db, id=x, rettype=rettype)
    seq = SeqIO.read(fet, format)
    fet.close()
    
    NCBI_seq.append(seq)

    time.sleep(0.5)
  
  return NCBI_seq

def ncbi_query_to_df(query_list):
    
    # create empty lists
    Accession = []
    DNA = []
    Genus = []
    Species = []
    gene_des = []
    version = []
    Protein = []
    full_sp_names = []
    
    # loop through the result list obtained from the NCBI search
    # may take over 10 minutes
    for query in query_list:
        for seq in query:
            # get genus nd speceis name seperately
            spe_name = seq.annotations["organism"]
            g_s_name = spe_name.split()

            # get and append protein sequence
            if seq.features:
                for feature in seq.features:
                    if feature.type == "CDS":
                        if "translation" in feature.qualifiers.keys():
                            pro_seq = feature.qualifiers['translation'][0]
                        

            # attached them to lists
            Accession.append(str(seq.name))
            Genus.append(str(g_s_name[0]))
            Species.append(str(g_s_name[1]))
            full_sp_names.append(str(g_s_name[0]) + ' ' + str(g_s_name[1]))
            gene_des.append(str(seq.description))
            version.append(str(seq.id))
            Protein.append(str(pro_seq))
    # create a dataframe for the information
    ncbi_q_df = pd.DataFrame(
        {'Accession': version,
        'Genus': Genus,
        'Species': Species,
        'Full_Species': full_sp_names,
        'Protein': Protein,
        'Gene_Description': gene_des
        })
    # This could be where the duplicate issue arises
    ncbi_q_df.drop_duplicates(subset=['Full_Species', 'Protein'],  keep='first', inplace=True)
    return ncbi_q_df
    


def ncbi_fetch_opsins(email, job_label='unnamed', out='unnamed', species_list=None):
    
    query_list = []
    
    # make a progress bar for tracking query progression. Based on length of the species list
    i=0
    with progressbar.ProgressBar(max_value=len(species_list),style='BouncingBar') as bar:
        for species in species_list:
            NCBI_seq = ncbi_fetch(email=email, 
                            term = f"{species}[Organism] AND (opsin[Title] OR rhodopsin[Title] OR Opn[Title] OR rh1[Title] OR rh2[Title] OR Rh1[Title] OR Rh2[Title]) NOT partial[Title] NOT voucher[All Fields] NOT kinase[All Fields] NOT similar[Title] NOT homolog[Title]")
            query_list.append(NCBI_seq)
            bar.update(i)
            i+=1
        bar.finish()
    try:
        dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        report_dir = f'mnm_on_{job_label}_{dt_label}'
        os.makedirs(report_dir)
        print('NCBI Queries Complete!\nNow Extracting and Formatting Results For DataFrame...\n')
        ncbi_q_df = ncbi_query_to_df(query_list=query_list)
        
        if out == 'unnamed':
            out = job_label
        ncbi_q_df.to_csv(path_or_buf=f"./{report_dir}/{out.replace(' ','_')}_ncbi_q_data.csv", index=False)
        print('DataFrame Formatted and Saved to CSV file for future use :)\n')
        fasta_file = f'{report_dir}/mined_{out.replace(" ","_")}_seqs.fasta'
        with open(fasta_file, 'w') as f:
            for id, seq in zip(ncbi_q_df['Accession'], ncbi_q_df['Protein']):
                f.write(f'>{id}\n{seq}\n')
        print('FASTA File Saved...\n')
        
        # Cleaning raw ncbi query df to keep only species entries that match our input sp list - returns the 'clean' df
        # Will also return a dataframe of the 'potential' hits - since it could be that the species that don't match is due to species mispelling or synonymous names
        ncbi_sp_hits = list(set(ncbi_query_df['Full_Species'].to_list()))
        #len(ncbi_sp_hits)
        intersection = list(set(ncbi_sp_hits) & set(species_list))
        #len(intersection)
        sp_no_hits  = list(set(species_list).symmetric_difference(intersection))
        #len(sp_no_hits)
        if len(sp_no_hits) > 0:
            print('Saving txt file with names of species that retrieved no results for opsins...\n')
            no_sp_hits_file = f'{report_dir}/species_w_no_hits.txt'
            with open(no_sp_hits_file, 'w') as f:
                for sp in sp_no_hits:
                    f.write(f'{sp}\n')
        
        sp_rnd_hits  = list(set(ncbi_sp_hits).symmetric_difference(intersection))
        if len(sp_no_hits) > 0:
            print('Saving txt file with names of species that retrieved results for opsins but are NOT in submitted species list...\n')
            #len(sp_rnd_hits)
            rnd_sp_hits_file = f'{report_dir}/potnetial_species_hits.txt'
            with open(rnd_sp_hits_file, 'w') as f:
                for sp in sp_rnd_hits:
                    f.write(f'{sp}\n')
            ncbi_query_df_cleaned = ncbi_query_df[~ncbi_query_df['Full_Species'].isin(sp_rnd_hits)]
            ncbi_query_df_potential_hits = ncbi_query_df[ncbi_query_df['Full_Species'].isin(sp_rnd_hits)]
            #ncbi_query_df_cleaned.shape
            #ncbi_query_df_potential_hits.shape
            print('Saving and returning cleaned dataframe with only species entries from species list...\n')
            ncbi_query_df_cleaned.to_csv(path_or_buf=f'{report_dir}/mnm_on_all_dbs_ncbi_q_data_cleaned.csv', index=False)
            print('Saving another dataframe with species that retrieved results for opsins but are NOT in submitted species list for further examination...\n')
            ncbi_query_df_potential_hits.to_csv(path_or_buf=f'{report_dir}/mnm_on_all_dbs_ncbi_q_potential_hits.csv', index=False)
            
            fasta_file = f'{report_dir}/mined_{out.replace(" ","_")}_seqs_cleaned.fasta'
            with open(fasta_file, 'w') as f:
                for id, seq in zip(ncbi_query_df_cleaned['Accession'], ncbi_query_df_cleaned['Protein']):
                    f.write(f'>{id}\n{seq}\n')
            print('Clean FASTA File Saved...\n')
            return(ncbi_query_df_cleaned)
        
        return(ncbi_q_df)
    
    except:
        print('Error occured when trying to save dataframe of NCBI query data!\nReturning query list instead...')
        return(query_list)


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

        elif '-' in lambda_max_str:
            # Range of values
            new_row = other_data.to_dict()
            # Remove non-numerical characters and take average
            new_row[lambda_max_column] = clean_single_value(lambda_max_str)
            new_rows.append(new_row)
        else:
            # Single value (potentially with non-numerical characters)
            new_row = other_data.to_dict()
            new_row[lambda_max_column] = clean_single_value(lambda_max_str)
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
        return np.nan  # Return NaN if no numbers are found
    
    # Handle multiple numbers
    if len(numbers) > 1:
        extracted_numbers = []
        for num in numbers:
            # If hyphen is present, split the values and convert to numbers
            if '-' in num:
                try:
                    start, end = map(float, num.split('-'))
                    extracted_numbers.extend([start, end])
                except ValueError:
                    return np.nan  # Return NaN if there's a problem with splitting or converting to float
            else:
                try:
                    extracted_numbers.append(float(num))
                except ValueError:
                    return np.nan # Return NaN if value cannot be converted to float
        return np.mean(extracted_numbers)
    
    # Handle single range or single number
    elif len(numbers) == 1:
        if '-' in numbers[0]:
            try:
                start, end = map(float, numbers[0].split('-'))
                return np.mean([start, end])
            except ValueError:
                return np.nan
        else:
            try:
                return float(numbers[0])
            except ValueError:
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
      #temp_df['LambdaMax'] = temp_df['LambdaMax'].astype(str) 
      #print(temp_df.head())

      # Rename the index column to its unique name
      #df = df.rename(columns={index_name: index_name})

      # Perform the merge (outer join to keep all data)
      if merged_df.empty:
          merged_df = temp_df
      else:
          merged_df = pd.concat([merged_df, temp_df]) 
          
    merged_df.index.name = 'comp_db_id'
    #merged_df.drop_duplicates(subset=['Full_Species', 'LambdaMax'], keep='first', inplace=True)
    # Sort the DataFrame by 'Accession' (descending) and then 'obs_id' (ascending)
    processed_df = merged_df.sort_values(by=['Accession', 'comp_db_id'], ascending=[False, True])
    # Drop duplicate rows based on 'Full_Species' and 'LambdaMax', keeping the first occurrence
    processed_df = processed_df.drop_duplicates(subset=['Full_Species', 'LambdaMax'], keep='first')
    processed_df.dropna(subset=['LambdaMax'], inplace=True)
    #re-sort by index and then reset index
    #processed_df = processed_df.sort_values('comp_db_id')
    processed_df = processed_df.reset_index(drop=True)
    processed_df.index.name = 'comp_db_id'
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    processed_df.to_csv(f'{report_dir}/vpod_comp_acc_dbs_{dt_label}.csv', index=True)

    #Clean the proccessed df column containing lambda max values of rows with mutliple entries or ranges
    cleaned_df = clean_lambda_max(processed_df.copy(), 'LambdaMax')  # Use df.copy() to avoid modifying the original DataFrame
    cleaned_df.index.name = 'comp_db_id'
    #save the final clean, merged df
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    cleaned_df.to_csv(f'{report_dir}/clean_vpod_comp_acc_dbs_{dt_label}.csv', index=True)
    
    return cleaned_df


import pandas as pd
import numpy as np
from Bio import Entrez
import time
# Provide your email to NCBI
Entrez.email = "sethfrazer@ucsb.edu"  # **IMPORTANT: Replace with your actual email**

def get_protein_sequence(accession):
    """
    Fetches the amino acid sequence from NCBI using the accession number.
    """
    try:
        # Introduce a delay to avoid overloading the server
        time.sleep(0.5)  # Adjust the delay as needed
        
        print(f"Fetching sequence for accession: '{accession}'")  # Debug print

        # Check if the accession contains extra characters (e.g., brackets)
        if any(c in accession for c in ['[', ']', '(', ')', '{', '}']):
            print(f"Skipping accession {accession} due to invalid characters.")
            return None
        
        # Check if accession is empty or contains only whitespace
        if not accession or accession.isspace():
            print(f"Skipping empty or whitespace accession: '{accession}'")
            return None

        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "gb")
        handle.close()
         # Extract sequence from FASTA record
        for i,feature in enumerate(record.features):
            if feature.type=='CDS':
                sequence = feature.qualifiers['translation'][0]
        return sequence

    except Exception as e:
        print(f"Error fetching sequence for accession '{accession}': {e}")  # Modified error message
        return None

def create_matched_df(mnm_merged_df, source_data):
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
                prediction_value = prediction_row['Prediction_Means']
                
                # Handle cases where accession is not in the columns
                try:
                    accession = prediction_row['Accession']
                except:
                    accession = prediction_row.name

                # Calculate absolute differences between the prediction and all measurements for the species
                species_measurements.loc[:, 'abs_diff'] = (species_measurements['LambdaMax'] - prediction_value).abs()
                
                # Handle Accession precedence
                
                closest_measurement_row = species_measurements.sort_values('abs_diff').iloc[0]
                source_accession = closest_measurement_row['Accession']
                
                #source_supersedes = False # Flag to indicate if source accession supersedes
                
                #if pd.notna(source_accession) and not any(c in source_accession for c in ['-', '–']):
                  # Source file Accession is valid and doesn't contain '-' or '–'
                  
                  #Only change accession if it's not already in mnm_merged_df
                #  if source_accession not in mnm_merged_df.index:
                #    accession = source_accession
                #    source_supersedes = True
                    
                #  closest_measurement_row = species_measurements[species_measurements['Accession'] == accession].iloc[0]
                
                min_diff = closest_measurement_row['abs_diff']
                
                existing_match_index = next((i for i, match in enumerate(matched_results) if match['Accession'] == accession), None)
                
                #if source_supersedes:
                    # Fetch protein sequence from NCBI
                #    protein_sequence = get_protein_sequence(accession)

                if (existing_match_index is not None):
                    if (min_diff < matched_results[existing_match_index]['abs_diff']):
                        # If it exists, compare the absolute differences and keep the better match
                        matched_results[existing_match_index].update({
                            'prediction_value': prediction_value,
                            'closest_measurement': closest_measurement_row['LambdaMax'],
                            'abs_diff': min_diff,
                            'comp_db_id': closest_measurement_row.name
                        })
                        if source_supersedes:
                            matched_results[existing_match_index]['Protein'] = protein_sequence
                else:
                    # If it doesn't exist, add the new match
                    
                    match_data = {
                        'Accession': accession,
                        'Full_Species': species,
                        'prediction_value': prediction_value,
                        'closest_measurement': closest_measurement_row['LambdaMax'],
                        'abs_diff': min_diff,
                        'comp_db_id': closest_measurement_row.name
                    }
                    
                    #if source_supersedes:
                    #    match_data['Protein'] = protein_sequence
                    
                    matched_results.append(match_data)

    # Create a new dataframe from the matched results
    matched_df = pd.DataFrame(matched_results)
    print(f'There were {i} unmatched species')
    
    #Match %Idenitities to Accessions
    iden_list = []
    prot_des_list = []
    aa_seq_list = []
    genus_list = []
    species_list = []
    mnm_merged_df.set_index('Accession', inplace=True)
    for _, d in matched_df.iterrows():
        acc = d['Accession']
        iden_list.append(mnm_merged_df.loc[acc]['%Identity_Nearest_VPOD_Sequence'])
        prot_des_list.append(mnm_merged_df.loc[acc]['Gene_Description'])
        
        # Only fetch sequence if not already present due to source_supersede
        #if pd.isna(d['Protein']):
        aa_seq_list.append(mnm_merged_df.loc[acc]['Protein'])
        #else:
        #    aa_seq_list.append(d['Protein'])
        
        genus_list.append(mnm_merged_df.loc[acc]['Genus'])
        species_list.append(mnm_merged_df.loc[acc]['Species'])
    matched_df['%Identity_Nearest_VPOD_Sequence'] = iden_list
    matched_df['Gene_Description'] = prot_des_list
    matched_df['Protein'] = aa_seq_list
    matched_df['Genus'] = genus_list
    matched_df['Species'] = species_list
    matched_df = matched_df.reindex(columns=['Accession','Genus','Species','%Identity_Nearest_VPOD_Sequence','prediction_value','closest_measurement','abs_diff','comp_db_id','Protein','Gene_Description','Notes'])
    
    # Group by 'comp_db_id' and count unique accessions
    grouped_counts = matched_df.groupby('comp_db_id')['Accession'].nunique()

    # Filter groups with more than one unique accession
    duplicate_comp_db_id_groups = grouped_counts[grouped_counts > 1]

    if not duplicate_comp_db_id_groups.empty:
        filtered_results = []

        # Iterate through each mnm_id with duplicates
        for comp_db_id in duplicate_comp_db_id_groups.index:
            # Filter rows with the current mnm_id
            duplicates = matched_df[matched_df['comp_db_id'] == comp_db_id]
            duplicates = pd.DataFrame(duplicates)
            duplicates = duplicates.sort_values(['abs_diff', 'Accession', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False, False])
            filtered_results.append(duplicates.iloc[0])
        
        # Combine filtered results with non-duplicate rows
        non_duplicates = matched_df[~matched_df['comp_db_id'].isin(duplicate_comp_db_id_groups.index)]
        final_filtered_df = pd.concat([pd.DataFrame(filtered_results), non_duplicates], ignore_index=True)
        final_filtered_df = final_filtered_df.sort_values(['abs_diff', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False])
        final_filtered_df.reset_index(drop=True, inplace=True)
        final_filtered_df.index.name = 'mnm_id'
    else:
        final_filtered_df = matched_df 
    
    return(final_filtered_df)

def mine_n_match(report_dir, source_file, ncbi_q_file, optics_pred_file, out='unnamed', err_filter = 15):

    try:
        source_data = pd.read_csv(source_file, index_col = 0)
    except:
        source_file = f'./{report_dir}/{source_file}'
        source_data = pd.read_csv(source_file, index_col = 0)
    
    try:
        source_data_indx_name = str(source_data.index.name)
    except:
        pass
    
    if 'Genus' in source_data.columns and 'Species' in source_data.columns and 'Full_Species' not in source_data.columns:
        source_data['Full_Species'] = source_data['Genus'] + '_' + source_data['Species']
    elif 'Genus' not in source_data.columns and 'Species' not in source_data.columns and 'Full_Species' in source_data.columns:
        gn_list = []
        sp_list = []
        og_sp_list = source_data['Full_Species'].to_list()
        #print(og_sp_list)
#        try:
        for sp in og_sp_list:
        #    print(sp)
            gn_list.append(sp.split(' ')[0])
            sp_list.append(sp.split(' ')[1:])
        source_data['Genus'] = gn_list
        source_data['Species'] = sp_list
#        except:
#            for sp in og_sp_list:
#                gn_list.append(sp.split('_')[0])
#                sp_list.append(sp.split('_')[1])
#            source_data['Genus'] = gn_list
#            source_data['Species'] = sp_list
            
    try:
        source_data.index.name = source_data_indx_name
    except:
        source_data.index.name = 'comp_db_id'

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
    mnm_merged_df = mnm_merged_df[mnm_merged_df['%Identity_Nearest_VPOD_Sequence'] != 'blastp unsuccessful']
    matched_df = create_matched_df(mnm_merged_df=mnm_merged_df, source_data=source_data)
    matched_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_results_id_filtered.csv", index=True)
    #note that prediction values from optics are taken from the mean of the bootstrap predictions
    final_err_filtered_df = matched_df[matched_df['abs_diff'] <= err_filter]
    #final_err_filtered_df = final_err_filtered_df[final_err_filtered_df['%Identity_Nearest_VPOD_Sequence'] != 'blastp unsuccessful']
    final_err_filtered_df = final_err_filtered_df[final_err_filtered_df['%Identity_Nearest_VPOD_Sequence'] != 100.000]
    final_err_filtered_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_final_results_err_filtered.csv", index=True)

    return(final_err_filtered_df)
    
def post_process_matching(report_dir, mnm_file):
    try:
        mnm_data = pd.read_csv(mnm_file)
    except:
        mnm_file = f'./{report_dir}/{mnm_file}'
        mnm_data = pd.read_csv(mnm_file)

    # Sort the dataframe by `abs_diff` in ascending order
    mnm_data = mnm_data.sort_values('abs_diff')

    # Drop duplicate `Accession` values, keeping only the first (lowest) `abs_diff` value
    mnm_data_unique = mnm_data.drop_duplicates(subset='Accession', keep='first')
    mnm_data_unique.reset_index(inplace=True, drop=True)
    mnm_data_unique.index.name = 'comp_db_id'    
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    mnm_data_unique.to_csv(path_or_buf=f".{report_dir}/mine_n_match_curated_{dt_label}.csv", index=True)
    
    mnm_duplicates = mnm_data[mnm_data.duplicated(subset=['Accession'], keep=False)]
    mnm_duplicates.to_csv(path_or_buf=f"./mine_n_match_duplicates_{dt_label}.csv", index=True)
    
    return(mnm_data_unique)


import requests
from requests.adapters import HTTPAdapter, Retry
import re
import random


def get_prots_from_acc(acc_list):
    
    protein_list = []
    for acc in acc_list:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "gb")
            handle.close()
            for i,feature in enumerate(record.features):
                if feature.type=='CDS':
                    aa = feature.qualifiers['translation'][0]
        except:
            aa = input(f"Accession# {accession} - Not Recognized by NCBI\nPlease Enter Target Sequence to Continue: ")
        
        protein_list.append(aa)
        
    return(protein_list)


def get_random_protein_fasta(query, max_outputs=2000, desired_seqs=5000, out=None):
    # Define the UniProt API endpoint
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    # Loop until all results are fetched
    # Define the parameters for the API request
    params = {
        "query": query,
        "format": "fasta",
        "size": 500,  # Fetch in batches of 500
    }

    # Send the API request
    ex_response = requests.get(base_url, params=params)

    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    all_sequences = []

    url = ex_response.url
    interactions = {}
    for batch, total in get_batch(url, session, re_next_link):
        sequences = batch.text.split(">")[1:]
            #primaryAccession, interactsWith = line.split('\t')
            #interactions[primaryAccession] = len(interactsWith.split(';')) if interactsWith else 0

            # Add the sequences to the list
        all_sequences.extend(sequences)
        
        if len(all_sequences) > desired_seqs:
            break

        print(f'{len(all_sequences)} / {total}')
        
    unique_sequences = list(set(all_sequences))
    
    # Randomly select 1000 sequences from the unique list
    if len(unique_sequences) >= max_outputs:
        random_sequences = random.sample(unique_sequences, desired_seqs)
    else:
        random_sequences = unique_sequences  # Take all unique sequences if less than 1000

    # Save the sequences to a file
    with open(out, "w") as f:
        
        #f.write(">" + ">".join(random_sequences)
        for entry in random_sequences:
            name,sequence = entry.split('\n',1)
            f.write(f'>{name.strip().replace("|","_").replace(" ","_").replace("/","")[0:30]}\n')
            f.write(f'{sequence}')
            
    return(random_sequences)

def get_next_link(headers, re_next_link):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url, session, re_next_link):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers, re_next_link)
        
        
def get_species_synonyms(species_name):
    """Fetches synonyms for a given species name from NCBI Taxonomy."""
    try:
        search_handle = Entrez.esearch(db="taxonomy", term=species_name, usehistory="y")
        search_results = Entrez.read(search_handle)
        query_key = search_results["QueryKey"]
        webenv = search_results["WebEnv"]

        fetch_handle = Entrez.esummary(db="taxonomy", query_key=query_key, webenv=webenv)
        fetch_results = Entrez.read(fetch_handle)

        synonyms = []
        for result in fetch_results:
            if 'Synonym' in result:
                synonyms.extend(result['Synonym'].split('; '))
        return synonyms
    except Exception as e:
        print(f"Error fetching synonyms for {species_name}: {e}")
        return []