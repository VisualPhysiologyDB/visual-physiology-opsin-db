import os
import re
import datetime
import time
import json
import copy
import pandas as pd
import numpy as np
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

# create a list for all the entries fetched from NCBI
  NCBI_seq =[]
  for x in q_list:
    # fetch result from previous search
    fet = Entrez.efetch(db=ncbi_db, id=x, rettype=rettype)
    seq = SeqIO.read(fet, format)
    fet.close()
    
    NCBI_seq.append(seq)

    time.sleep(0.25)
  
  return NCBI_seq

def get_species_taxonomy(species_name, email):
    """Fetches synonyms for a given species name from NCBI Taxonomy."""

    Entrez.email = email    # Always tell NCBI who you are

    # Search for the species
    handle = Entrez.esearch(db="taxonomy", term=species_name, api_key = "1efb120056e1cea873ba8d85d6692abd5d09")
    record = Entrez.read(handle)
    taxonomy = {}
    
    try:
        tax_id = record["IdList"][0]
        # Fetch the taxonomy record
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        record = Entrez.read(handle)[0]
        #print(record)
        # Extract synonyms
        synonyms = []
        if "OtherNames" in record:
            for name in record["OtherNames"]["Synonym"]:
                if name not in synonyms:
                    if '(' in name:
                        if ',' in name:
                            rename = name.split('(')[0].strip()
                            if (rename not in synonyms) and (rename != species_name): 
                                synonyms.append(rename)
                        else:
                            rename = name.split('(')[1].strip()
                            rename = rename.split(')')[0].strip() + ' ' + rename.split(')')[1].strip()
                            if (rename not in synonyms) and (rename != species_name):
                                synonyms.append(rename)
                            rename2 = name.split('(')[0].strip() + rename.split(' ')[1]
                            if (rename2 not in synonyms) and (rename2 != species_name):
                                synonyms.append(rename2)
                    else:
                        if ',' in name:
                            rename = name.split(',')[0]
                            rename = rename.split(' ')[0] + ' ' + rename.split(' ')[1]
                            if (rename not in synonyms) and (rename != species_name):
                                synonyms.append(rename)     
                        else:
                            synonyms.append(name)
                else:
                    pass
        for lineage in record["LineageEx"]:
            rank = lineage["Rank"]
            if rank == "phylum": 
                taxonomy["Phylum"] = lineage["ScientificName"]
            elif rank == "subphylum":
                taxonomy["Subphylum"] = lineage["ScientificName"]
            elif rank == "class":
                taxonomy["Class"] = lineage["ScientificName"]
    except:
        synonyms = []
      
    # Asign species synonyms  
    taxonomy['Synonyms'] = synonyms
    
    # Ensure all desired ranks are present
    for rank in ["Phylum", "Subphylum", "Class"]:
        if rank not in taxonomy:
            taxonomy[rank] = "Unknown"  # Or None if you prefer

    handle.close()

    return taxonomy

def get_sp_taxon_dict(species_list, email, taxon_file, sp_taxon_dict={}):
    #sp_taxon_dict = {}
    all_keys = sp_taxon_dict.keys()

    for species in species_list:
        queried = False
        if species in all_keys:
            queried = True
        else: 
            for x in range(10):
                try:
                    if queried == False:
                        taxonomy = get_species_taxonomy(species, email)
                        #print(f"{species} synonyms: {synonyms}")
                        sp_taxon_dict[species] = taxonomy
                        queried = True
                    else:
                        pass
                except:
                    time.sleep(1)
            if queried == False:
                with open(taxon_file, 'w') as f:
                    json.dump(sp_taxon_dict, f, indent=4)  # indent for pretty formatting
                raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')

    return sp_taxon_dict

def ncbi_query_to_df(query_list, species_list, species_taxon_dict, email):
    
    # create empty lists
    Accession = []
    #DNA = []
    Phylum = []
    Subphylum = []
    Class = []
    Genus = []
    Species = []
    gene_des = []
    version = []
    Protein = []
    full_sp_names = []
    Sp_syn_used = []
    
    # loop through the result list obtained from the NCBI search
    for query, sp in zip(query_list, species_list):
        # Get genus and speceis name seperately
        for seq in query:
            g_s_name = sp.split(' ', 1)
            # Search the dictionary of synonymous species names to see 
            # if this is the primary name or synonym.
            entry_spe_name = seq.annotations["organism"]
            
            l1 = len(entry_spe_name)
            l2 = len(sp)
            
            temp = entry_spe_name.split(' ')
            temp2 = sp.split(' ')
            if ((entry_spe_name == sp) or (entry_spe_name[:l1-1] == sp) or (entry_spe_name == sp[:l2-1])):
                found_w_synonym = False
                Phylum.append(species_taxon_dict[sp]["Phylum"])
                Subphylum.append(species_taxon_dict[sp]["Subphylum"])
                Class.append(species_taxon_dict[sp]["Class"])
            elif len(temp) == 3 or len(temp2) == 3:
                if (len(temp) == 3) and ((sp == str(temp[0]+' '+temp[1])) or (sp == str(temp[0]+' '+temp[2]))):
                    found_w_synonym = True
                    Phylum.append(species_taxon_dict[sp]["Phylum"])
                    Subphylum.append(species_taxon_dict[sp]["Subphylum"])
                    Class.append(species_taxon_dict[sp]["Class"])
                elif (len(temp2) == 3) and ((entry_spe_name == str(temp2[0]+' '+temp2[1])) or (entry_spe_name == str(temp2[0]+' '+temp2[2]))):
                    found_w_synonym = True
                    Phylum.append(species_taxon_dict[sp]["Phylum"])
                    Subphylum.append(species_taxon_dict[sp]["Subphylum"])
                    Class.append(species_taxon_dict[sp]["Class"])
                else:
                    g_s_name = entry_spe_name.split(' ', 1)
                    found_w_synonym = False
                    # If the species for this entry is not in the species list or synonyms dict then we will fetch the phylogeny now
                    queried = False
                    for x in range(50):
                        if queried == False:
                            try:        
                                temp_taxon = get_species_taxonomy(entry_spe_name, email)
                                queried = True
                            except:
                                pass
                        else:
                            pass
                    if queried == False:
                        raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')
                    Phylum.append(temp_taxon["Phylum"])
                    Subphylum.append(temp_taxon["Subphylum"])
                    Class.append(temp_taxon["Class"])   
            else:
                if entry_spe_name in species_taxon_dict[sp]["Synonyms"]:
                    found_w_synonym = True
                    Phylum.append(species_taxon_dict[sp]["Phylum"])
                    Subphylum.append(species_taxon_dict[sp]["Subphylum"])
                    Class.append(species_taxon_dict[sp]["Class"])
                else:
                    g_s_name = entry_spe_name.split(' ', 1)
                    found_w_synonym = False
                    # If the species for this entry is not in the species list or synonyms dict then we will fetch the phylogeny now
                    queried = False
                    for x in range(50):
                        if queried == False:
                            try:        
                                temp_taxon = get_species_taxonomy(entry_spe_name, email)
                                queried = True
                            except:
                                pass
                        else:
                            pass
                    if queried == False:
                        raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')
                    Phylum.append(temp_taxon["Phylum"])
                    Subphylum.append(temp_taxon["Subphylum"])
                    Class.append(temp_taxon["Class"])
            if found_w_synonym == True:
                Sp_syn_used.append(entry_spe_name)
            else:
                Sp_syn_used.append('NA')

            # get and append protein sequence
            if seq.features:
                for feature in seq.features:
                    if feature.type == "CDS":
                        if "translation" in feature.qualifiers.keys():
                            pro_seq = feature.qualifiers['translation'][0]
                        
            # Append all meta data to corresponding lists
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
        'Phylum': Phylum,
        'Subphylum': Subphylum,
        'Class': Class,
        'Genus': Genus,
        'Species': Species,
        'Full_Species': full_sp_names,
        'Protein': Protein,
        'Gene_Description': gene_des,
        'Species_Synonym_Used': Sp_syn_used
        })

    # Drop duplicates where the species names and protein sequences are the same...
    ncbi_q_df.drop_duplicates(subset=['Full_Species', 'Protein'],  keep='first', inplace=True)
    ncbi_q_df = ncbi_q_df.reset_index(drop=True)
    return ncbi_q_df
    


def ncbi_mining(email, job_label='unnamed', out='unnamed', species_list=None, query):
    
    print('Creating Job Directory\n')
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    report_dir = f'mnm_data/mnm_on_{job_label}_{dt_label}'
    os.makedirs(report_dir)
    
    print('Saving Species Query List to Text\n')
    with open(f"{report_dir}/species_queried.txt", "w") as f:
        for sp in species_list:
            f.write(str(sp) + "\n")
    
    # Create a taxonomic dictionary, including species synonyms, for all the species in the unique species list  
    print('Constructing Taxon Dictionary, Including Species Synonyms\n')
    # Check to see if an existing taxonomy dictionary already exists to save time.
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
    if (list(species_taxon_dict.keys()) > list(existing_taxon_dict.keys())):         #print('Saving Updated Dictionary') 
        try:
            with open(taxon_file, 'w') as f:
                json.dump(species_taxon_dict, f, indent=4)  # indent for pretty formatting
        except FileNotFoundError:
            print(f"Error: File '{taxon_file}' can't be saved...\n")

    print('Taxon Dictionary Complete!\n')
    
    # List to append query responses to
    query_list = []
    # make a progress bar for tracking query progression. Based on length of the species list
    print('Starting Queries to NCBI for Opsin Sequences\n')
    i=0
    with progressbar.ProgressBar(max_value=len(species_list),style='BouncingBar') as bar:
        for species in species_list:
            try:
                temp = species.split(' ')
            except:
                raise Exception(f'Species Name Causing Error: {species}')
            if len(species_taxon_dict[species]['Synonyms']) > 0:
                sp_for_query = f'("{species}"[Organism]'
                if len(temp) == 3:
                    sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism]'
                    sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism]'

                for syn in species_taxon_dict[species]['Synonyms']:
                    sp_for_query+= f' OR "{syn}"[Organism]'
                sp_for_query+=')'
            elif len(temp) == 3:
                sp_for_query = f'("{species}"[Organism]'
                sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism]'
                sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism]'
                sp_for_query+=')'
            else:
                sp_for_query = f'"{species}"[Organism]'
            #print(sp_for_query)
            queried = False
            for x in range(50):
                if queried == False:
                    try:            
                        NCBI_seq = ncbi_fetch(email=email, 
                                        term = f"{sp_for_query} AND {query}")
                        queried = True
                    except:
                        time.sleep(1)
                        pass
            if queried == False:
                print('Uh-oh, species query failed.\nThis is likely a back-end issue with the querying process\nIf this message continues to appear, please manually interrupt the qury process and restart...')
                    
            query_list.append(NCBI_seq)
            bar.update(i)
            i+=1
        bar.finish()
        
    print('NCBI Queries Complete!\nNow Extracting and Formatting Results For DataFrame...\n')
    #maybe add the syn-species dictionary to the function below
    ncbi_query_df = ncbi_query_to_df(query_list=query_list, species_list=species_list, species_taxon_dict=species_taxon_dict, email=email)
    
    if out == 'unnamed':
        out = job_label
    ncbi_query_df.to_csv(path_or_buf=f"./{report_dir}/{out.replace(' ','_')}_ncbi_q_data.csv", index=False)
    print('DataFrame Formatted and Saved to CSV file for future use :)\n')
    fasta_file = f'{report_dir}/mined_{out.replace(" ","_")}_seqs.fasta'
    with open(fasta_file, 'w') as f:
        for id, seq in zip(ncbi_query_df['Accession'], ncbi_query_df['Protein']):
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
    if len(sp_rnd_hits) > 0:
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
        ncbi_query_df_cleaned.to_csv(path_or_buf=f'{report_dir}/{out}_ncbi_q_data_cleaned}.csv', index=False)
        print('Saving another dataframe with species that retrieved results for opsins but are NOT in submitted species list for further examination...\n')
        ncbi_query_df_potential_hits.to_csv(path_or_buf=f'{report_dir}/{out}_ncbi_q_potential_hits.csv', index=False)
        
        fasta_file = f'{report_dir}/mined_{out.replace(" ","_")}_seqs_cleaned.fasta'
        with open(fasta_file, 'w') as f:
            for id, seq in zip(ncbi_query_df_cleaned['Accession'], ncbi_query_df_cleaned['Protein']):
                f.write(f'>{id}\n{seq}\n')
        print('Clean FASTA File Saved...\n')
        return(ncbi_query_df_cleaned, report_dir)
        
    return(ncbi_query_df, report_dir)
    


import requests
from requests.adapters import HTTPAdapter, Retry
import re
import random


def get_prots_from_acc(acc_list):
    
    protein_list = []
    for accession in acc_list:
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
        search_handle = Entrez.esearch(db="taxonomy", term=species_name, api_key = "1efb120056e1cea873ba8d85d6692abd5d09", usehistory="y")
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
    
    
def fasta_to_dataframe(fasta_file):
  """
  Converts a fasta file into a pandas DataFrame, extracting species name, 
  opsin type, accession, and sequence from the header line and sequence data.

  Args:
    fasta_file: Path to the fasta file.

  Returns:
    A pandas DataFrame with columns: 'species_name', 'opsin_type', 
    'accession', and 'sequence'.
  """

  data = []
  with open(fasta_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
      try:
        header = record.description.split("---")
        species_name = header[0].replace("_", " ")
        opsin_type = header[1].split("-",1)[0]
        accession = header[1].split("-",1)[1]
      except:
        header = record.description.split("_")        
        if str(record.description).count('_') > 1:
            try:
              accession = header[0] + header[1]
              species_name = header[2].split('--')[0].replace('-',' ')
              opsin_type = header[2].split('--')[1]
            except:
              try:
                accession = header[0] + header[1]
                species_name = header[2] + header[3].split('-')[0]
                opsin_type = header[3].split('-')[1]
              except:
                try:
                  accession = header[0] + header[1]
                  species_name = header[2].split('-')[0] + ' ' + header[2].split('-')[1]
                  opsin_type = header[2].split('-')[2]
                except:
                  print(header)
                  raise Exception('Header with unique formating seems to have caused an error')
        else:
            try:
              header = record.description.split("_")
              accession = header[0]
              species_name = header[1].split('--')[0].replace('-',' ')
              opsin_type = header[1].split('--')[1]
            except:
              try:
                accession = header[0]
                species_name = header[1] + header[2].split('-')[0]
                opsin_type = header[2].split('-')[1]
              except:
                try:
                  accession = header[0]
                  species_name = header[1].split('-')[0] + ' ' + header[1].split('-')[1]
                  opsin_type = header[1].split('-')[2]
                except:
                  print(header)
                  raise Exception('Header with unique formating seems to have caused an error')
            
      aln_sequence = str(record.seq).replace(" ", "").replace("\n", "").replace("=", "")
      sequence = str(aln_sequence).replace('-','')
      seq_len = len(sequence)
      data.append([species_name, opsin_type, accession, aln_sequence, sequence, seq_len])
        
  df = pd.DataFrame(data, columns=['species_name', 'opsin_type', 'accession', 'aln_sequence','sequence', 'seq_length'])
  return df