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
import progressbar

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


from pygbif import species
def correct_species_name(species_name):
    
    """Attempts to find a corrected species name and its associated higher taxa using the GBIF backbone taxonomy.

    This function uses the `pygbif` library to query the Global Biodiversity Information Facility (GBIF) 
    backbone taxonomy. It attempts to find a match for the input `species_name` and returns the 
    corrected species name (if different) along with a higher taxonomic group (family, order, or genus).
    The higher taxonomic group can be used to find related taxa from NCBI if a species name is not recognized.

    Args:
        species_name (str): The species name to check.

    Returns:
        tuple: A tuple containing two strings:
            - The corrected species name (if found), or "No match found" otherwise.
            - The associated higher taxa (family, order, or genus) of the corrected species name (if found), or "No match found" otherwise.

    Notes:
        - The function retries the GBIF query up to 10 times if it fails initially.
        - The function prioritizes returning the family, then the order, then the genus if a species is found.
        - If no match is found in the GBIF backbone, or if there's an error during the query, 
          it returns ("No match found", "No match found").
    """
    queried = False
    try:
        for x in range(10):
            if queried == False:
                result = species.name_backbone(name=species_name, rank='species', verbose=True)
                queried = True
    except:
        pass
    
    try:
        if 'species' in result and 'family' in result:
            return result['species'], result['family']
        elif ('species' in result) and ('order' in result):
            return result['species'], result['order']
        elif ('species' in result):
            return result['species'], result['genus']
        else:
            return "No match found", "No match found"
    except:
        return "No match found", "No match found"



def get_species_taxonomy(species_name, email, record_alt=False, use_higher_taxa=False, higher_taxa=None):
    
    """Fetches taxonomic information, including synonyms, for a given species or higher taxon from NCBI Taxonomy.

    This function retrieves taxonomic details from the NCBI Taxonomy database using the Biopython Entrez module.
    It can retrieve synonyms for a species or information about a higher taxon if `use_higher_taxa` is set to True.
    It also attempts to determine the Phylum, Subphylum and Class of the species or taxon.

    Args:
        species_name (str): The scientific name of the species to look up.
        email (str): Your email address, required by NCBI Entrez.
        record_alt (bool, optional): If True, includes the input `species_name` in the synonyms list. Defaults to False.
        use_higher_taxa (bool, optional): If True, searches for a higher taxon instead of a species. Defaults to False.
        higher_taxa (str, optional): The name of the higher taxon to search for when `use_higher_taxa` is True.

    Returns:
        dict: A dictionary containing the taxonomic information, including:
            - Synonyms (list): A list of synonymous names for the species.
            - Phylum (str): The phylum of the species/taxon.
            - Subphylum (str): The subphylum of the species/taxon.
            - Class (str): The class of the species/taxon.
            - If a rank is not found, it will be assigned "Unknown".
    """
    Entrez.email = email    # Always tell NCBI who you are
    queried = False
    synonyms = []

    for x in range(10):
        if queried == False:
            try:
                if use_higher_taxa == False:
                    # Search for the species
                    handle = Entrez.esearch(db="taxonomy", term=species_name, api_key = api_key)
                else:
                    handle = Entrez.esearch(db="taxonomy", term=higher_taxa, api_key = api_key)

                record = Entrez.read(handle)
                taxonomy = {}
            
                tax_id = record["IdList"][0]
                # Fetch the taxonomy record
                handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                record = Entrez.read(handle)[0]
                #if the code makes it this far, then there was a successful query
                queried = True
                #print(record)
                taxonomy["TaxId"] = record["TaxId"]
                if species_name != record["ScientificName"]:
                    synonyms.append(record["ScientificName"])

                # Extract synonyms
                if ("OtherNames" in record) and (use_higher_taxa == False):
                    for name in record["OtherNames"]["Synonym"]:
                        if name not in synonyms and name != species_name:

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
                pass
            
    if record_alt == True:
        synonyms.append(species_name)
        
    # Asign species synonyms  
    taxonomy['Synonyms'] = synonyms
    
    # Ensure all desired ranks are present
    for rank in ["Phylum", "Subphylum", "Class"]:
        if rank not in taxonomy:
            taxonomy[rank] = "Unknown"
            
    handle.close()
    
    return taxonomy

def get_sp_taxon_dict(species_list, email, taxon_file, sp_taxon_dict={}):
    
    """Builds a dictionary of taxonomic information for a list of species.

    This function takes a list of species names and constructs a dictionary where keys are 
    species names and values are dictionaries containing taxonomic information 
    (e.g., Phylum, Subphylum, Class, Order, Synonyms) for each species. It uses the 
    `get_species_taxonomy` function to retrieve the information from NCBI, and can also
    leverage the `correct_species_name` function to handle potential misspellings or outdated names.

    Args:
        species_list (list): A list of species names (strings).
        email (str): Your email address, required for NCBI Entrez queries.
        taxon_file (str): The file path to save the resulting taxonomy dictionary (JSON format). This is used in error handling to save progress.
        sp_taxon_dict (dict, optional): An existing taxonomy dictionary to update. Defaults to an empty dictionary.

    Returns:
        dict: A dictionary where keys are species names and values are dictionaries 
              containing taxonomic information for each species.

    Raises:
        Exception: If there is an issue querying NCBI for taxonomic information. In this case the
                   current progress is saved to the taxon_file, and the user is prompted to restart.
                   This is implemented to avoid issues with querying NCBI's servers.

    Notes:
        - If a species is already present in `sp_taxon_dict`, its information is not fetched again,
        unless the existing entry is missing crucial taxonomic levels, in which case it will try to look
        up the species again using the `correct_species_name` function if possible.
        - The function attempts to correct species names using the `correct_species_name` function if the initial
          taxonomy lookup fails or returns incomplete information (missing Phylum, Order, and Class).
        - The `correct_species_name` function tries to find a corrected name or higher taxonomic
          group from the Global Biodiversity Information Facility (GBIF).
        - Intermediate results are saved to `taxon_file` in case of errors during the process. This allows
          you to resume the process from where it left off if an error occurs, rather than starting from scratch.
    """
    
    all_keys = sp_taxon_dict.keys()

    for species in species_list:
        #if (species in all_keys) and (sp_taxon_dict[species]['Phylum']!="Unknown"):
        if (species in all_keys):
            pass
        else: 
            try:
                taxonomy = get_species_taxonomy(species, email)
                #print(f"{species} synonyms: {synonyms}")
                sp_taxon_dict[species] = taxonomy
                if (len(sp_taxon_dict[species]['Synonyms']) == 0) and (sp_taxon_dict[species]['Phylum']=="Unknown") and (sp_taxon_dict[species]['Subphylum']=="Unknown") and (sp_taxon_dict[species]['Class']=="Unknown"):
                    try:
                        corrected_name, higher_taxa_name = correct_species_name(species)
                        if (corrected_name != "No match found") and (corrected_name != species):
                            taxonomy = get_species_taxonomy(corrected_name, email, record_alt=True)
                            sp_taxon_dict[species] = taxonomy

                            if (sp_taxon_dict[species]['Phylum']=="Unknown") and (sp_taxon_dict[species]['Subphylum']=="Unknown") and (sp_taxon_dict[species]['Class']=="Unknown"):
                                try:
                                    taxonomy = get_species_taxonomy(corrected_name, email, record_alt=True, use_higher_taxa=True, higher_taxa=higher_taxa_name)
                                    sp_taxon_dict[species] = taxonomy
                                    #print(f"Corrected Species' ({species}) Higher Taxa Name Used to Find Taxa!\n")
                                except:
                                    sp_taxon_dict[species] = taxonomy
                                    #print(f'Finding Corrected Species Name for {species} Failed Due to Error Using Higher Taxa as Substitute.\n')  
                            else:
                                pass
                                #print(f'Corrected Species Name for {species} Used to Find Taxa!\n')
                        else:
                            pass
                            #print(f'Finding Corrected Species Name for {species} Failed.\n') 
                    except:
                        print(f'Finding Corrected Species Name for {species} Failed Due to Error.\n')  
            except:
                with open(taxon_file, 'w') as f:
                    json.dump(sp_taxon_dict, f, indent=4)  # indent for pretty formatting
                raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')

    return sp_taxon_dict

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
            Genus.append(str(g_s_name[0]))
            Species.append(str(g_s_name[1]))
            full_sp_names.append(str(g_s_name[0]) + ' ' + str(g_s_name[1]))
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
        try:
            temp = species.split(' ')
        except:
            raise Exception(f'Species Name Causing Error: {species}')
        
        if len(species_taxon_dict[species]['Synonyms']) > 0:
            sp_for_query = f'("{species}"[Organism] OR "{species}"[Title]'
            if (len(temp) == 3) and ('(' not in species) and ('.' not in species):
                sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism] OR "{temp[0]} {temp[1]}"[Title]'
                sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism]  OR "{temp[0]} {temp[2]}"[Title]'
                for syn in species_taxon_dict[species]['Synonyms']:
                    if (syn != species) and (f"{temp[0]} {temp[1]}" != syn) and (f"{temp[0]} {temp[2]}" != syn):
                        sp_for_query+= f' OR "{syn}"[Organism] OR {syn}"[Title]'
                        sp_for_query+=')'
            else:
                for syn in species_taxon_dict[species]['Synonyms']:
                    if (syn != species):
                        sp_for_query+= f' OR "{syn}"[Organism] OR {syn}"[Title]'
                        sp_for_query+=')'
                        
        elif (len(temp) == 3) and ('(' not in species) and ('.' not in species):
            sp_for_query = f'("{species}"[Organism]'
            sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism] OR "{temp[0]} {temp[1]}"[Title]'
            sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism] OR "{temp[0]} {temp[2]}"[Title]'
            sp_for_query+=')'
            
        else:
            sp_for_query = f'"{species}"[Organism]'
        #print(sp_for_query)
        queried = False
        for x in range(50):
            if queried == False:
                try:            
                    NCBI_seq = ncbi_fetch_alt(email=email, 
                                    term = f'{sp_for_query} AND (opsin[Title] OR rhodopsin[Title] OR OPN[Title] OR rh1[Title] OR rh2[Title] OR Rh1[Title] OR Rh2[Title]) NOT partial[Title] NOT voucher[All Fields] NOT kinase[All Fields] NOT kinase-like[All Fields] NOT similar[Title] NOT homolog[Title] NOT enhancer[Title]')
                                    #term = f'{sp_for_query} AND (opsin[Title] OR rhodopsin[Title] OR OPN[Title] OR rh1[Title] OR rh2[Title] OR Rh1[Title] OR Rh2[Title]) NOT partial[Title] NOT voucher[All Fields] NOT kinase[All Fields] NOT kinase-like[All Fields] NOT similar[Title] NOT homolog[Title] NOT opsin-like[Title] NOT enhancer[Title]')

                    queried = True
                except:
                    time.sleep(1)
                    pass
        if queried == False:
            print('Uh-oh, species query failed.\nThis is likely a back-end issue with the querying process\nIf this message continues to appear, please manually interrupt the qury process and restart...')
                
        query_list.append(NCBI_seq)
        #bar.update(i)
        #i+=1
        #time.sleep(0.5)

    #bar.finish()
        
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
    processed_df['Accession'] = processed_df['Accession'].str.replace('â€“','-').replace('–','-').replace(',Â','-').replace(' ','').replace('-','-')
    processed_df.reset_index(drop=True, inplace=True)
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
    cleaned_df_name = f'{report_dir}/clean_vpod_comp_acc_dbs_{dt_label}.csv'
    cleaned_df.to_csv(cleaned_df_name, index=True)
    return cleaned_df, cleaned_df_name


def create_matched_df(report_dir, mnm_merged_df, source_data):
    
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
                prediction_value = prediction_row['Prediction_Medians']
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


def mine_n_match(email, report_dir, source_file, ncbi_q_file, optics_pred_file, out='unnamed', err_filter=15):
    
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
                    species = sp.split('Â', 1)[0]
                    gn_list.append(genus)
                    sp_list.append(species)
                except:
                    raise Exception(f'Problematic species name: {sp}')
        source_data['Genus'] = gn_list
        source_data['Species'] = sp_list

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
    mnm_merged_df = mnm_merged_df[mnm_merged_df['%Identity_Nearest_VPOD_Sequence'] != 'blastp unsuccessful']
    mnm_merged_df.reset_index(inplace=True, drop=True)
    #mnm_merged_df.index.name = 'comp_db_id'
    mnm_merged_df.to_csv(path_or_buf=f"./{report_dir}/mnm_merged_df.csv", index=True)

    matched_df = create_matched_df(report_dir=report_dir, mnm_merged_df=mnm_merged_df, source_data=source_data)
    matched_df.to_csv(path_or_buf=f"./{report_dir}/mnm_on_{out}_results_id_filtered.csv", index=True)
    #note that prediction values from optics are taken from the mean of the bootstrap predictions
    final_err_filtered_df = matched_df[matched_df['abs_diff'] <= err_filter]
    final_err_filtered_df = final_err_filtered_df[final_err_filtered_df['%Identity_Nearest_VPOD_Sequence'] != 100.000]
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