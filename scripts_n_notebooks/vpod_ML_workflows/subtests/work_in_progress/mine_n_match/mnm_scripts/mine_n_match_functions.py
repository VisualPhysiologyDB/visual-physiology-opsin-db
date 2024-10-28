import os
import re
import datetime
import time
import pandas as pd
from Bio import Entrez, SeqIO
from progress.bar import ShadyBar


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
  print(record_full)
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
    bar = ShadyBar('Processing Sequences', max=len(query_list), charset='ascii')
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
        bar.next()
    # create a dataframe for the information
    ncbi_q_df = pd.DataFrame(
        {'Accession': version,
        'Genus': Genus,
        'Species': Species,
        'Full_Species': full_sp_names,
        'Protein': Protein,
        'Gene_Description': gene_des
        })
    ncbi_q_df.drop_duplicates(subset=['Full_Species', 'Protein'],  keep='first')
    bar.finish()
    return ncbi_q_df
    


def ncbi_fetch_opsins(email, job_label='unnamed', out='unnamed', species_list=None):
    
    query_list = []
    
    # make a progress bar for tracking query progression. Based on length of the species list
    bar = ShadyBar('Processing Sequences', max=len(species_list), charset='ascii')
    for species in species_list:
        NCBI_seq = ncbi_fetch(email=email, 
                        term = f"{species}[ORGN] AND (opsin[All Fields] AND complete[All Fields] AND cds[All Fields] NOT voucher)")
        query_list.append(NCBI_seq)
        bar.next()
    bar.finish()
    print('NCBI Queries Complete!\nNow Extracting and Formatting Results For DataFrame...\n')
    ncbi_q_df = ncbi_query_to_df(query_list=query_list)
    
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    report_dir = f'mnm_on_{job_label}_{dt_label}'
    os.makedirs(report_dir)
    if out == 'unnamed':
        out = job_label
    ncbi_q_df.to_csv(path_or_buf=f"./{report_dir}/{out.replace(' ','_')}_ncbi_q_data.csv", index=False)
    print('DataFrame Formatted and Saved to CSV file for future use :)\n')
    fasta_file = f'{report_dir}/mined_{out.replace(' ','_')}_seqs.fasta'
    with open(fasta_file, 'w') as f:
        for id, seq in zip(ncbi_q_df['Accession'], ncbi_q_df['Protein']):
            f.write(f'>{id}\n{seq}\n')
    print('FASTA File Saved...\n')

    return(ncbi_q_df)

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
            try:
                accession = prediction_row['Accession']
            except:
                accession = prediction_row.name

            # Calculate absolute differences between the prediction and all measurements for the species
            species_measurements.loc[:, 'abs_diff'] = (species_measurements['LambdaMax'] - prediction_value).abs()
            # Find the closest measurement (handling ties)
            closest_measurement_row = species_measurements.sort_values('abs_diff').iloc[0]
            min_diff = closest_measurement_row['abs_diff']
            
            existing_match_index = next((i for i, match in enumerate(matched_results) if match['Accession'] == accession), None)

            if (existing_match_index is not None): 
                if (min_diff < matched_results[existing_match_index]['abs_diff']):
                    # If it exists, compare the absolute differences and keep the better match
                    matched_results[existing_match_index].update({
                        'prediction_value': prediction_value,
                        'closest_measurement': closest_measurement_row['LambdaMax'],
                        'abs_diff': min_diff,
                        'max_id': closest_measurement_row.name
                    })
            else:
                # If it doesn't exist, add the new match
                matched_results.append({
                    'Accession': accession,
                    'Full_Species': species,
                    'prediction_value': prediction_value,
                    'closest_measurement': closest_measurement_row['LambdaMax'],
                    'abs_diff': min_diff,
                    'max_id': closest_measurement_row.name
                })
            
    # Create a new dataframe from the matched results
    matched_df = pd.DataFrame(matched_results)
    print(f'There were {i} unmatched species')
    
    #Match %Idenitities to Accessions
    iden_list = []
    prot_des_list = []
    aa_seq_list = []
    genus_list = []
    species_list = []
    for _, d in matched_df.iterrows():
        acc = d['Accession']
        iden_list.append(mnm_merged_df.loc[acc]['%Identity_Nearest_VPOD_Sequence'])
        prot_des_list.append(mnm_merged_df.loc[acc]['Gene_Description'])
        aa_seq_list.append(mnm_merged_df.loc[acc]['Protein'])
        genus_list.append(mnm_merged_df.loc[acc]['Genus'])
        species_list.append(mnm_merged_df.loc[acc]['Species'])
    matched_df['%Identity_Nearest_VPOD_Sequence'] = iden_list
    matched_df['Gene_Description'] = prot_des_list
    matched_df['Protein'] = aa_seq_list
    matched_df['Genus'] = genus_list
    matched_df['Species'] = species_list
    matched_df = matched_df.reindex(columns=['Accession','Genus','Species','%Identity_Nearest_VPOD_Sequence','prediction_value','closest_measurement','abs_diff','max_id','Protein','Gene_Description','Notes'])
    
    # Group by 'max_id' and count unique accessions
    grouped_counts = matched_df.groupby('max_id')['Accession'].nunique()

    # Filter groups with more than one unique accession
    duplicate_max_id_groups = grouped_counts[grouped_counts > 1]

    if not duplicate_max_id_groups.empty:
        filtered_results = []

        # Iterate through each max_id with duplicates
        for max_id in duplicate_max_id_groups.index:
            # Filter rows with the current max_id
            duplicates = matched_df[matched_df['max_id'] == max_id]

            # Sort by abs_diff (ascending) and then percent_identity (descending)
            duplicates = duplicates.sort_values(['abs_diff', 'Accession'], ascending=[True, False])
            duplicates = duplicates.sort_values(['abs_diff', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False])

            # Keep only the first row (lowest abs_diff, highest percent_identity if tied)
            filtered_results.append(duplicates.iloc[0])

        # Combine filtered results with non-duplicate rows
        non_duplicates = matched_df[~matched_df['max_id'].isin(duplicate_max_id_groups.index)]
        final_filtered_df = pd.concat([pd.DataFrame(filtered_results), non_duplicates], ignore_index=True)
        final_filtered_df = final_filtered_df.sort_values(['abs_diff', '%Identity_Nearest_VPOD_Sequence'], ascending=[True, False])
        final_filtered_df.reset_index(drop=True, inplace=True)
        final_filtered_df.index.name = 'mnm_id'
    else:
        final_filtered_df = matched_df 
    
    return(final_filtered_df)
    

def mine_n_match(report_dir, source_file, ncbi_q_file, optics_pred_file, out='unnamed', err_filter = 15):

    try:
        source_data = pd.read_csv(source_file)
    except:
        source_file = f'./{report_dir}/{source_file}'
        source_data = pd.read_csv(source_file)    
        
    if 'Genus' in source_data.columns and 'Species' in source_data.columns and 'Full_Species' not in source_data.columns:
        source_data['Full_Species'] = source_data['Genus'] + '_' + source_data['Species']
    elif 'Genus' not in source_data.columns and 'Species' not in source_data.columns and 'Full_Species' in source_data.columns:
        gn_list = []
        sp_list = []
        og_sp_list = source_data['Full_Species'].to_list()
        try:
            for sp in og_sp_list:
                gn_list.append(sp.split(' ')[0])
                sp_list.append(sp.split(' ')[1])
            source_data['Genus'] = gn_list
            source_data['Species'] = sp_list
        except:
            for sp in og_sp_list:
                gn_list.append(sp.split('_')[0])
                sp_list.append(sp.split('_')[1])
            source_data['Genus'] = gn_list
            source_data['Species'] = sp_list
    source_data.index.name = 'mnm_id'

    try:
        ncbi_data = pd.read_csv(ncbi_q_file)
    except:
        ncbi_q_file = f'./{report_dir}/{ncbi_q_file}'
        ncbi_data = pd.read_csv(ncbi_q_file)
    try:
        pred_data = pd.read_csv(optics_pred_file)
    except:
        optics_pred_file = f'./{report_dir}/{optics_pred_file}'
        pred_data = pd.read_csv(optics_pred_file)
        
    mnm_merged_df = pd.merge(ncbi_data, pred_data, left_index=True, right_index=True)
    matched_df = create_matched_df(mnm_merged_df=mnm_merged_df, source_data=source_data)
    matched_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_results_id_filtered.csv", index=True)
    #note that prediction values from optics are taken from the mean of the bootstrap predictions
    final_err_filtered_df = matched_df[matched_df['abs_diff'] <= err_filter]
    final_err_filtered_df = final_err_filtered_df[final_err_filtered_df['%Identity_Nearest_VPOD_Sequence'] != 'blastp unsuccessful']
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
    mnm_data_unique.index.name = 'mnm_id'    
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    mnm_data_unique.to_csv(path_or_buf=f".{report_dir}/mine_n_match_curated_{dt_label}.csv", index=True)
    
    mnm_duplicates = mnm_data[mnm_data.duplicated(subset=['Accession'], keep=False)]
    mnm_duplicates.to_csv(path_or_buf=f"./mine_n_match_duplicates_{dt_label}.csv", index=True)
    
    return(mnm_data_unique)