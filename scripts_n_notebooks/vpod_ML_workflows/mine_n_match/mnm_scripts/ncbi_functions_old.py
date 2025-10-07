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
from tqdm import tqdm

api_key = "1efb120056e1cea873ba8d85d6692abd5d09"


def ncbi_fetch_alt(email, term, ncbi_db="nuccore", rettype="gb", format="genbank"):
    """Fetches sequences from NCBI's databases using Entrez.

    This function queries an NCBI database (default: "nuccore") using a provided search term and
    retrieves all corresponding sequence records. It uses the provided email for identification
    and an API key for increased request limits. It also caches results in a JSON file to avoid
    redundant NCBI queries.

    Args:
        email (str): Your email address. Required by NCBI for Entrez queries.
        term (str): The search term to query the NCBI database.
        ncbi_db (str, optional): The NCBI database to search. Defaults to "nuccore".
        rettype (str, optional): The retrieval type for efetch. Defaults to "gb".
        format (str, optional): The format of the sequence records to retrieve. Defaults to "genbank".
        
    Returns:
        list: A list of Biopython SeqRecord objects representing the fetched sequence records.
                
    Notes:
        - An API key should be set using `api_key = 'your_api_key'` in the main code before calling this function for increased requests per second.
        - The function automatically pauses for 0.25 seconds between fetches to avoid overloading 
          the NCBI servers, adhering to the rate limit even with an API key.
        - The function returns all results, even if the number of results exceeds the default `retmax` limit of Entrez.esearch.

    """
    Entrez.email = email  # Always tell NCBI who you are

    # Create queries directory if it doesn't exist
    queries_dir = "data_sources/queries"
    if not os.path.exists(queries_dir):
        os.makedirs(queries_dir)

    cache_file = os.path.join(queries_dir, "cached_ncbi_queries.json")

    # Load cache from file if it exists
    try:
        with open(cache_file, "r") as f:
            cache = json.load(f)
    except FileNotFoundError:
        cache = {}

    handle = Entrez.esearch(db=ncbi_db, term=term, api_key=api_key)
    record = Entrez.read(handle)
    full_res = int(record["Count"])

    handle_full = Entrez.esearch(db=ncbi_db, term=term, api_key=api_key, retmax=full_res + 1)
    record_full = Entrez.read(handle_full)
    q_list = record_full["IdList"]
    #print(f'Here is the query ID return list:\n{q_list}\n')
    NCBI_seq = []
    for x in q_list:
        if x in cache:
            # Fetch from cache
            #print(f"Fetching record {x} from cache.")
            seq_record = SeqIO.read(io.StringIO(cache[x]), format)
            NCBI_seq.append(seq_record)
        else:
            # Fetch from NCBI
            #print(f"Fetching record {x} from NCBI.")
            fet = Entrez.efetch(db=ncbi_db, id=x, rettype=rettype)
            seq = SeqIO.read(fet, format)
            fet.close()

            NCBI_seq.append(seq)
            
            # Store in cache as a string
            with io.StringIO() as string_handle:
                SeqIO.write(seq, string_handle, format)
                cache[x] = string_handle.getvalue()

            time.sleep(0.25)

    # Save cache to file
    with open(cache_file, "w") as f:
        json.dump(cache, f, indent=4)

    return NCBI_seq

def ncbi_query_to_df(query_list, species_list, species_taxon_dict, email):
    
    """Converts a list of NCBI query results into a Pandas DataFrame.

    This function takes a list of Biopython SeqRecord objects (the results of NCBI queries) 
    and extracts relevant information such as accession numbers, taxonomic information, 
    gene descriptions, and protein sequences. It then organizes this data into a Pandas DataFrame.

    Args:
        query_list (list): A list of lists, where each inner list contains Biopython SeqRecord objects 
                           returned from NCBI queries for a specific species.
        species_list (list): A list of species names (strings) that were used in the NCBI queries.
        species_taxon_dict (dict): A dictionary containing taxonomic information for each species in 
                                   `species_list`, including synonyms. The keys are species names, 
                                   and the values are dictionaries with keys like "Phylum", "Subphylum", 
                                   "Class", and "Synonyms".
        email (str): Your email address, used for querying NCBI in the case of missing taxonomy.

    Returns:
        pandas.DataFrame: A DataFrame containing the extracted information from the NCBI query results. 
                          The DataFrame has the following columns:
                            - Accession: The NCBI accession number (version).
                            - Phylum: The phylum of the species.
                            - Subphylum: The subphylum of the species.
                            - Class: The class of the species.
                            - Genus: The genus of the species.
                            - Species: The species name (without the genus).
                            - Full_Species: The full species name (genus + species).
                            - Protein: The amino acid sequence of the protein.
                            - Gene_Description: A description of the gene.
                            - Species_Synonym_Used: The synonymous species name used in the query if a synonym was used, 
                                                    otherwise 'NA'.

    Raises:
        Exception: If a species query fails during the taxonomic lookup for an unknown species. This exception
                   is raised after 50 retries. This is highly unlikely, and likely indicates a back-end problem.
                   If you encounter this exception please save any work and restart the script.

    Notes:
      - The function handles cases where the species name in the NCBI record might be a synonym 
        of the name used in the original query.
      - It attempts to look up the taxonomic classification (Phylum, Subphylum, Class) 
        in the provided `species_taxon_dict`.
      - If a species is not found in `species_taxon_dict` and is not a synonym of any included species, it will attempt to fetch the taxonomic information from NCBI directly using the `get_species_taxonomy` function and the provided email. It will retry this process up to 50 times for each species if the initial attempt fails.
      - Duplicate entries (based on species name and protein sequence) are removed, keeping only the first occurrence.
    """
    
    # create empty lists
    Accession = []
    dna = []
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

            if seq.seq:
                try:
                    if len(str(seq.seq)) < 10000:
                        dna_seq = str(seq.seq).strip().replace('\n','').replace(',','').replace('\t','')
                    else:
                        dna_seq = ""
                except:
                    dna_seq = ""
            else:
                dna_seq = ""

            # get and append protein sequence
            if seq.features:
                for feature in seq.features:
                    if feature.type == "CDS":
                        if "translation" in feature.qualifiers.keys():
                            pro_seq = feature.qualifiers['translation'][0]
                        
            # Append all meta data to corresponding lists
            Accession.append(str(seq.name))
            
            try:
                Species.append(str(g_s_name[1]))
                Genus.append(str(g_s_name[0]))
                full_sp_names.append(str(g_s_name[0]) + ' ' + str(g_s_name[1]))
            except:
                try:
                    genus = sp.split('Â', 1)[0]
                    species = sp.split('Â', 1)[1]
                    Species.append(species)
                    Genus.append(genus)
                except:
                    Genus.append(str(g_s_name[0]))
                    Species.append('')
                    full_sp_names.append(str(g_s_name[0]))
                    print(f'This species, {sp}, is not a complete species name!')

                
            gene_des.append(str(seq.description))
            version.append(str(seq.id))
            try:
                dna.append(str(dna_seq))
            except:
                print(f'DNA-Seq for {seq.id} is undefined...\n')
                dna.append('')
            try:
                Protein.append(str(pro_seq))
            except:
                print(f'Protein-Seq for {seq.id} is undefined...\n')
                Protein.append('')
            
    # create a dataframe for the information
    ncbi_q_df = pd.DataFrame(
        {'Accession': version,
        'Phylum': Phylum,
        'Subphylum': Subphylum,
        'Class': Class,
        'Genus': Genus,
        'Species': Species,
        'Full_Species': full_sp_names,
        'DNA': dna,
        'Protein': Protein,
        'Gene_Description': gene_des,
        'Species_Synonym_Used': Sp_syn_used
        })

    # Drop duplicates where the species names and protein sequences are the same...
    ncbi_q_df.drop_duplicates(subset=['Full_Species', 'Protein'],  keep='first', inplace=True)
    ncbi_q_df = ncbi_q_df.reset_index(drop=True)
    return ncbi_q_df
    


def ncbi_fetch_opsins(email, job_label='unnamed', out='unnamed', species_list=None, filter_len_lower=300, filter_len_upper=700):
    
    """Fetches opsin sequences from NCBI for a given list of species.

    This function performs the following steps:
    1. Creates a directory to store results, labeled with the job name and timestamp.
    2. Saves the list of queried species to a text file.
    3. Constructs a taxonomic dictionary for the species, including synonyms, 
       using an existing dictionary if available to improve efficiency.
    4. Queries NCBI's Nucleotide database for opsin sequences for each species,
       using the species name and its synonyms, as well as search terms related to opsin proteins.
    5. Handles query failures with a retry mechanism (up to 50 attempts with 1 second delays).
    6. Parses the results into a Pandas DataFrame.
    7. Saves the DataFrame to a CSV file.
    8. Saves the retrieved protein sequences to a FASTA file.
    9. Filters the results to include only species present in the original input list.
    10. Saves a list of species with no hits to a text file.
    11. Saves a list of species with hits that were not in the original list (potential hits due to synonyms or misspellings) to a separate file and dataframe.
    12. Returns the cleaned DataFrame (containing only species from the input list) and the path to the report directory.

    Args:
        email (str): Your email address. Required for NCBI Entrez queries.
        job_label (str, optional): A label for the job. Used in directory and file names. Defaults to 'unnamed'.
        out (str, optional): The base name for output files. Defaults to the `job_label`.
        species_list (list): A list of species names (strings) to query.

    Returns:
        tuple: A tuple containing:
            - pandas.DataFrame: A DataFrame containing the cleaned NCBI query results (only species from the input list). 
                                If no species from the list had hits, returns the original unfiltered dataframe.
            - str: The path to the directory where results are saved.

    Raises:
        Exception: If a species name causes issues in the query generation process. This is highly unlikely,
                   and would only be seen if a species name is not formatted as a string.
    """
    
    
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
    if (list(species_taxon_dict.keys()) >= list(existing_taxon_dict.keys())):         
        #print('Saving Updated Dictionary') 
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
    #with progressbar.ProgressBar(max_value=len(species_list),style='BouncingBar') as bar:
    for species in tqdm(species_list, 
                   desc="Processing species queries",
                   colour="#CF9FFF",
                   bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                   dynamic_ncols=True,
                   ascii=" ▖▝▗▚▞█"):           

        # Build standard query for each species, at a time
        sp_for_query = build_species_query(species_taxon_dict, species)
        complete_ncbi_query = build_visual_opsin_query(sp_for_query)
        
        queried = False
        for x in range(50):
            if queried == False:
                try:            
                    NCBI_seq = ncbi_fetch_alt(email=email, term = complete_ncbi_query)
                    queried = True
                except:
                    time.sleep(0.2)
                    pass
        if queried == False:
            print('Uh-oh, species query failed.\nThis is likely a back-end issue with the querying process\nIf this message continues to appear, please manually interrupt the qury process and restart...')
                
        query_list.append(NCBI_seq)
        
    print('NCBI Queries Complete!\nNow Extracting and Formatting Results For DataFrame...\n')
    ncbi_query_df = ncbi_query_to_df(query_list=query_list, species_list=species_list, species_taxon_dict=species_taxon_dict, email=email)
    ncbi_query_df['Prot_Len'] = ncbi_query_df['Protein'].str.len()
    ncbi_query_df = ncbi_query_df[(ncbi_query_df['Prot_Len'] >= filter_len_lower) & (ncbi_query_df['Prot_Len'] <= filter_len_upper)]

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
        rnd_sp_hits_file = f'{report_dir}/potential_species_hits.txt'
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
        return(ncbi_query_df_cleaned, report_dir)
        
    return(ncbi_query_df, report_dir)

