import os
import re
import datetime
import time
import json
import copy
import pandas as pd
import numpy as np
import io
from multiprocessing import Manager
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez, SeqIO
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

from mnm_scripts.taxonomy import get_species_taxonomy, get_sp_taxon_dict

# --- Configuration ---
# The API key allows for up to 10 requests per second.
Entrez.api_key = "1efb120056e1cea873ba8d85d6692abd5d09"

def load_ncbi_cache():
    queries_dir = "data_sources/queries"
    if not os.path.exists(queries_dir):
        os.makedirs(queries_dir)

    cache_file = os.path.join(queries_dir, "cached_ncbi_queries_protein.json")
    try:
        with open(cache_file, "r") as f:
            cache = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f'Cannot read cached file, here is the error:\n{e}')
        cache = {}
    return(cache, cache_file)

def load_taxon_dict(taxon_file):
    print('Constructing Taxon Dictionary, Including Species Synonyms\n')
    os.makedirs(os.path.dirname(taxon_file), exist_ok=True)
    existing_taxon_dict = {}
    if os.path.isfile(taxon_file):
        try:
            with open(taxon_file, 'r') as f:
                existing_taxon_dict = json.load(f)
            print('Existing Taxon Dictionary Found! Checking if we need to update it...\n')
        except (FileNotFoundError, json.JSONDecodeError):
            print(f"Could not load existing taxon file at '{taxon_file}'. Starting fresh.\n")
    return(existing_taxon_dict)


def ncbi_fetch_alt(term, ncbi_db="protein", rettype="gb", format="genbank", cache = {}):
    """
    Fetches sequences from NCBI's databases using an efficient batch method.
    
    This function queries an NCBI database (default: "nuccore") using a provided search term and
    retrieves all corresponding sequence records. It uses the provided email for identification
    and an API key for increased request limits. It also caches results in a JSON file to avoid
    redundant NCBI queries.

    This  function first uses Entrez.esearch to get a list of IDs, then
    fetches the corresponding records in large batches. This dramatically
    reduces the number of API calls and speeds up data retrieval. It also
    caches results to avoid redundant downloads.

    Args:
        term (str): The search term to query the NCBI database.
        ncbi_db (str, optional): The NCBI database to search. Defaults to "protein".
        rettype (str, optional): The retrieval type for efetch. Defaults to "gb".
        format (str, optional): The record format. Defaults to "genbank".

    Returns:
        list: A list of Biopython SeqRecord objects.
    """
    # Step 1: Get the list of IDs for the search term
    handle = Entrez.esearch(db=ncbi_db, term=term, retmax='100000') # Get a large number of IDs
    search_results = Entrez.read(handle)
    handle.close()
    
    id_list = search_results["IdList"]
    count = int(search_results["Count"])
    
    if count == 0:
        return []

    NCBI_seq = []
    
    # Check cache first for all IDs
    non_cached_ids = [seq_id for seq_id in id_list if seq_id not in cache]
    for seq_id in id_list:
        if seq_id in cache:
            # Load from cache
            try:
                seq_record = SeqIO.read(io.StringIO(cache[seq_id]), format)
                NCBI_seq.append(seq_record)
            except Exception:
                # Handle potential format error in cached data
                if seq_id in non_cached_ids: non_cached_ids.remove(seq_id)
                pass


    # Step 2: Fetch non-cached records in batches
    batch_size = 500
    for start in (range(0, len(non_cached_ids), batch_size)):
        end = min(len(non_cached_ids), start + batch_size)
        
        try:
            fetch_handle = Entrez.efetch(
                db=ncbi_db,
                id=",".join(non_cached_ids[start:end]),
                rettype=rettype,
                retmode="text",
                api_key=Entrez.api_key,
            )
            
            records = list(SeqIO.parse(fetch_handle, format))
            for record in records:
                NCBI_seq.append(record)
                # Store new record in cache
                with io.StringIO() as string_handle:
                    SeqIO.write(record, string_handle, format)
                    # Use record.id which often includes version
                    cache_id = record.id.split('.')[0]
                    cache[cache_id] = string_handle.getvalue()
            fetch_handle.close()
        except Exception as e:
            print(f"Error fetching or parsing batch: {e}")
            continue

    return NCBI_seq


def ncbi_query_to_df(query_list, species_list, species_taxon_dict, email, taxon_file):
    """
    Converts a list of NCBI query results into a Pandas DataFrame.
    
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

    Notes:
      - The function handles cases where the species name in the NCBI record might be a synonym 
        of the name used in the original query.
      - It attempts to look up the taxonomic classification (Phylum, Subphylum, Class) 
        in the provided `species_taxon_dict`.
      - If a species is not found in `species_taxon_dict` and is not a synonym of any included species, it will attempt to fetch the taxonomic information from NCBI directly using the `get_species_taxonomy` function and the provided email. It will retry this process up to 50 times for each species if the initial attempt fails.
      - Duplicate entries (based on species name and protein sequence) are removed, keeping only the first occurrence.

    Returns:
        pandas.DataFrame: A DataFrame containing the extracted information from the NCBI query results. 
    """
    Accession, Phylum, Subphylum, Class, Genus, Species = [], [], [], [], [], []
    gene_des, Protein, full_sp_names, Sp_syn_used = [], [], [], []
    
    for i in tqdm(range(len(query_list)), 
            desc="Formatting species queries",
            colour="#CF9FFF",
            bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
            dynamic_ncols=True,
            ascii=" ▖▝▗▚▞█"):
        
        query = query_list[i]
        sp = species_list[i]
        for seq in query:
            g_s_name = sp.split(' ', 1)
            entry_spe_name = seq.annotations.get("organism", sp)
            found_w_synonym = False

            if (entry_spe_name not in species_taxon_dict.keys()) and (entry_spe_name not in species_taxon_dict.get(sp, {}).get("Synonyms", [])):
                temp_sp_list = [entry_spe_name]
                species_taxon_dict = get_sp_taxon_dict(temp_sp_list, email, taxon_file, sp_taxon_dict=species_taxon_dict)
                

            # --- Taxonomy and Species Name Logic ---
            if (entry_spe_name == sp) or (entry_spe_name in species_taxon_dict.get(sp, {}).get("Synonyms", [])) or (sp in species_taxon_dict.get(entry_spe_name, {}).get("Synonyms", [])) or (f'{g_s_name[0]} sp.' == entry_spe_name) or (sp == f'{entry_spe_name.split(" ")[0]} sp.'):
                found_w_synonym = entry_spe_name != sp
            else:
                g_s_name = entry_spe_name.split(' ', 1)
                
            tax_info = species_taxon_dict.get(sp, {})
            Phylum.append(tax_info.get("Phylum", "Unknown"))
            Subphylum.append(tax_info.get("Subphylum", "Unknown"))
            Class.append(tax_info.get("Class", "Unknown"))
            Sp_syn_used.append(entry_spe_name if found_w_synonym else 'NA')
            
            pro_seq = str(seq.seq) if seq.seq else ""

            # Append all meta data to corresponding lists
            Accession.append(str(seq.id))
            Protein.append(pro_seq)
            gene_des.append(str(seq.description))
            
            if len(g_s_name) > 1:
                Genus.append(str(g_s_name[0]))
                Species.append(' '.join(g_s_name[1:]))
                full_sp_names.append(f"{g_s_name[0]} {g_s_name[1]}")
            else:
                Genus.append(str(g_s_name[0]))
                Species.append('sp.')
                full_sp_names.append(f"{g_s_name[0]} sp.")
            
    ncbi_q_df = pd.DataFrame({
        'Accession': Accession, 'Phylum': Phylum, 'Subphylum': Subphylum,
        'Class': Class, 'Genus': Genus, 'Species': Species,
        'Full_Species': full_sp_names, 'Protein': Protein,
        'Gene_Description': gene_des, 'Species_Synonym_Used': Sp_syn_used
    })

    ncbi_q_df.drop_duplicates(subset=['Full_Species', 'Protein'], keep='first', inplace=True)
    return ncbi_q_df.reset_index(drop=True)
         
def _fetch_for_species(species, species_taxon_dict, ncbi_cache):
    """
    Helper function to build and run a query for a single species.
    This function is designed to be called by a thread pool executor.
    """
    try:
        sp_for_query = build_species_query(species_taxon_dict, species)
        complete_ncbi_query = build_visual_opsin_query(sp_for_query)
        return ncbi_fetch_alt(term=complete_ncbi_query, cache=ncbi_cache)
    except Exception as e:
        # Print an error message and return an empty list if a query fails
        print(f"An error occurred while processing species '{species}': {e}")
        return []

def ncbi_fetch_opsins(email, species_list, job_label='unnamed', out='unnamed', filter_len_lower=300, filter_len_upper=500):
    """
    Fetches opsin sequences from NCBI Protein DB.
    """
    Entrez.email = email   # Always tell NCBI who you are
    
    print('Creating Job Directory\n')
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    report_dir = f'mnm_data/mnm_on_{job_label}_{dt_label}'
    os.makedirs(report_dir, exist_ok=True)
    
    print('Saving Species Query List to Text\n')
    with open(f"{report_dir}/species_queried.txt", "w") as f:
        for sp in species_list:
            f.write(str(sp) + "\n")
    
    # load taxon dict
    taxon_file = './data_sources/taxonomy/ncbi_taxon_dict.json'
    existing_taxon_dict = load_taxon_dict(taxon_file=taxon_file) 
    # update taxon dict if neccessary
    species_taxon_dict = get_sp_taxon_dict(email=email,species_list=species_list, taxon_file=taxon_file, sp_taxon_dict=copy.deepcopy(existing_taxon_dict))
    
    if species_taxon_dict.keys() > existing_taxon_dict.keys():
        print('Updated Taxon Dictionary\n') 
    else:
        print('No need to update! Carry on :) \n')


    print('Starting Queries to NCBI for Opsin Sequences\n')
    
    results_dict = {}
    # load ncbi_cache
    ncbi_cache, cache_file = load_ncbi_cache()
    manager = Manager()
    mp_ncbi_cache_dict = manager.dict(ncbi_cache)  # Initialize a multiprocess available cached dict with cached data

    # Use ThreadPoolExecutor for I/O-bound tasks. 
    # max_workers=10 is a safe number to avoid overwhelming the NCBI API.
    with ThreadPoolExecutor(max_workers=10) as executor:
        # Map each future to its corresponding species name for tracking
        future_to_species = {
            executor.submit(_fetch_for_species, species, species_taxon_dict, mp_ncbi_cache_dict): species 
            for species in species_list
        }
        
        # Use tqdm to create a progress bar as futures complete
        progress_bar_opts = {
            "total": len(species_list),
            "desc": "Processing species queries",
            "colour": "#CF9FFF",
            "bar_format": "{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
            "dynamic_ncols": True,
            "ascii": " ▖▝▗▚▞█"
        }
        for future in tqdm(as_completed(future_to_species), **progress_bar_opts):
            species = future_to_species[future]
            try:
                results_dict[species] = future.result()
            except Exception as exc:
                print(f'\nAn exception occurred for {species}: {exc}')
                results_dict[species] = []

    updated_ncbi_cache_dict = dict(mp_ncbi_cache_dict) 
    # Save the taxon dictionary if it doesn't yet exist or if it has been updated since being loaded 
    if list(updated_ncbi_cache_dict.keys()) > list(ncbi_cache.keys()):  
        #print('Saving Updated Dictionary') 
        try:
            with open(cache_file, 'w') as f:
                json.dump(updated_ncbi_cache_dict, f, indent=4)  # indent for pretty formatting
        except FileNotFoundError:
            print(f"Error: Cached prediction file can't be saved...\n")
            

    # Reconstruct the list of results in the original order of species_list
    # This ensures the downstream DataFrame creation function works correctly
    query_list = [results_dict[species] for species in species_list]
        
    print('NCBI Queries Complete!\nNow Extracting and Formatting Results For DataFrame...\n')
    ncbi_query_df = ncbi_query_to_df(query_list=query_list, species_list=species_list, species_taxon_dict=species_taxon_dict, email=email, taxon_file=taxon_file)
    ncbi_query_df['Prot_Len'] = ncbi_query_df['Protein'].str.len()
    ncbi_query_df = ncbi_query_df[(ncbi_query_df['Prot_Len'] >= filter_len_lower) & (ncbi_query_df['Prot_Len'] <= filter_len_upper)]

    if out == 'unnamed':
        out = job_label
    csv_path = f"./{report_dir}/{out.replace(' ','_')}_ncbi_q_data.csv"
    ncbi_query_df.to_csv(path_or_buf=csv_path, index=False)
    print(f'DataFrame Formatted and Saved to {csv_path}\n')
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
        
        
    print("Query workflow complete.")

    return(ncbi_query_df, report_dir)    

def build_species_query(species_taxon_dict: dict, species: str) -> str:
    """
    Constructs a the species portion of the Entrez query string for visual opsins
    for a given species.

    Args:
        species_name: The scientific name of the species to search for
                      (e.g., "Danio rerio").

    Returns:
        A formatted string ready to be used as the 'term' parameter in an
        Entrez API call.
    """
    try:
        temp = species.strip().split(' ')
    except:
        raise Exception(f'Species Name Causing Error: {species}')
    
    if len(species_taxon_dict[species]['Synonyms']) > 0:
        sp_for_query = f'("{species}"[Organism] OR "{species}"[Title]'
        if (len(temp) == 3) and ('(' not in species) and ('.' not in species):
            sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism] OR "{temp[0]} {temp[1]}"[Title]'
            sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism]  OR "{temp[0]} {temp[2]}"[Title]'
            for syn in species_taxon_dict[species]['Synonyms']:
                if (syn != species) and (f"{temp[0]} {temp[1]}" != syn) and (f"{temp[0]} {temp[2]}" != syn) and syn.count(' ') > 0:
                    sp_for_query+= f' OR "{syn}"[Organism] OR {syn}"[Title]'
        else:
            for syn in species_taxon_dict[species]['Synonyms']:
                if (syn != species):
                    sp_for_query+= f' OR "{syn}"[Organism] OR "{syn}"[Title]'
        
        sp_for_query+=')'
                    
    elif (len(temp) == 3) and ('(' not in species) and ('.' not in species):
        sp_for_query = f'("{species}"[Organism] OR "{species}"[Title]'
        sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism] OR "{temp[0]} {temp[1]}"[Title]'
        sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism] OR "{temp[0]} {temp[2]}"[Title]'
        sp_for_query+=')'
        
    else:
        sp_for_query = f'"{species}"[Organism]'
    #print(sp_for_query)    
    
    return(sp_for_query)
         
# --- Helper function to format query parts ---
def format_query_part(terms, field_tag):
    # Add quotes around multi-word terms and apply the field tag
    formatted_terms = [f'"{term}"{field_tag}' if ' ' in term else f'{term}{field_tag}' for term in terms]
    return f'({" OR ".join(formatted_terms)})'

def format_exclusion_part(terms, field_tag):
    # Add quotes and tags, then join with " NOT "
    formatted_terms = [f'"{term}"{field_tag}' if ' ' in term else f'{term}{field_tag}' for term in terms]
    return " NOT ".join(formatted_terms)

def build_visual_opsin_query(sp_for_query: str) -> str:
    """
    Constructs a detailed NCBI Entrez query string for visual opsins
    for a given species.

    Args:
        species_name: The scientific name of the species to search for
                      (e.g., "Danio rerio").

    Returns:
        A formatted string ready to be used as the 'term' parameter in an
        Entrez API call.
    """

    # 1. Broad search terms to find the general family
    # We'll search for these in [All Fields] to cast a wide net initially
    base_terms_all_fields = [
        "opsin",
        "rhodopsin",
        "rhabdomeric*opsin",
        "c*opsin",
        "r*opsin"
    ]

    # 2. Specific INCLUSION terms for visual opsins, targeted to the [Title]
    # Grouping by opsin type makes this much easier to manage
    visual_opsin_terms_title = [
        # RH1/Rod opsin
        "RH1*", "rho*", "rod*",
        # SWS1/UV opsin
        "OPN1SW*", "SWS1*", "blue*sensitive*", "uv*sensitive", "ultra*violet*sensitive*",
        # SWS2/Blue opsin
        "OPN1S2*", "SWS2*", "short*wave*sensitive*",
        # RH2/Green opsin
        "RH2*",
        # LWS/Red opsin
        "OPN1LW*", "LWS*", "long*wave*sensitive*", "red*sensitive*",
        # MWS/Green opsin
        "OPN1MW*", "MWS*", "medium*wave*sensitive*", "green*sensitive*",
        # General visual terms
        "cone*opsin*", "visual*pigment*", "photoreceptor*opsin*", "opsin","rhodopsin", "rhabdomeric*opsin", "c*opsin", "r*opsin"
    ]

    # 3. EXCLUSION terms for non-visual opsins, targeted to the [Title]
    non_visual_opsin_terms_title = [
        "melanopsin*", "encephalopsin*", "tmt*opsin*", "neuropsin*", "pinopsin*",
        "parietopsin*", "parapinopsin*", "peropsin*", "rgr opsin*", "va*opsin*",
        "vertebrate*ancient*", "exorh*", "extra*ocular*", "non*visual*"
        # Grouped OPN numbers for non-visual opsins
        "opn*3*", "opn*4*", "opn*5*", "opn*6*", "opn*7*", "opn*8*", "opn*9*",
        "tmt*1", "tmt*2", "tmt*3", "tmt*4*", "tmt*5*", "tmt*6*", "tmt*7*",
        # Grouped "Opsin X" formats
        "opsin*3", "opsin*4*", "opsin*5", "opsin*6", "opsin*7", "opsin*8", "opsin*9"
    ]

    # 4. EXCLUSION terms for quality control (fragments, predictions, etc.)
    quality_control_terms_title = [
        "partial", "fragment", "incomplete", "putative", 
        #"similar", "homolog",
        "psuedogene", "transcript variant", "exon", "intron", "segment"
    ]

    # 5. EXCLUSION terms for other related genes or regulatory elements
    other_related_terms_title = [
        "kinase", "GTPase", "transducin", "arrestin", "enhancer",
        "promoter", "regulatory"
    ]

    # --- Assemble the final query ---
    query_parts = []

    # Part 1: Species and Base Search
    query_parts.append(f'{sp_for_query}')
    query_parts.append(f'AND {format_query_part(base_terms_all_fields, "[All Fields]")}')

    # Part 2: Specific Visual Opsin Inclusions
    query_parts.append(f'AND {format_query_part(visual_opsin_terms_title, "[Title]")}')

    # Part 3: Exclusions
    # We create one large exclusion block
    all_exclusions = (non_visual_opsin_terms_title +
                    quality_control_terms_title +
                    other_related_terms_title)
    query_parts.append(f'NOT {format_exclusion_part(all_exclusions, "[Title]")}')

    # Join all major parts into the final query string
    final_query = ' '.join(query_parts)
    # for now we are keeping this hard-coded since the way shown above would produce bugs occasionally
    #final_query = f'''({sp_for_query} AND ("opsin"[All Fields] OR "rhodopsin"[All Fields] OR “rhabdomeric opsin"[All Fields] OR “c-opsin”[All Fields] OR “r-opsin”[All Fields]) AND ("RH1"[Title] OR "rh1"[Title] OR "rho"[Title] OR “rh”[Title] OR "RH"[Title] OR "Rho"[Title] OR "RHO"[Title] OR "rod"[Title] OR "RH2"[Title] OR "rh2"[Title] OR "OPN1SW"[Title] OR "opn1sw"[Title] OR "SWS1"[Title] OR "sws1"[Title] OR "blue-sensitive"[Title] OR "blue sensitive"[Title] OR "green-sensitive"[Title] OR "green sensitive"[Title] OR "ultra-violet-sensitive"[Title] "ultra-violet sensitive"[Title] OR "OPN1S2"[Title] OR "opn1s2"[Title] OR "SWS2"[Title] OR "sws2"[Title] OR "short-wave-sensitive"[Title] OR "short-wave-sensitive"[Title] OR "OPN1LW"[Title] OR "opn1lw"[Title] OR "LWS"[Title] OR "lws"[Title] OR "long-wave-sensitive"[Title] OR "long-wave sensitive"[Title] OR "long-wavelength sensitive"[Title] OR "red-sensitive"[Title] OR "red sensitive"[Title] OR "red-wavelength"[Title] OR "OPN1MW"[Title] OR "opn1mw"[Title] OR "MWS"[Title] OR "mws"[Title] OR "medium-wave sensitive"[Title] OR "cone opsin"[Title] OR "visual pigment"[Title] OR "photoreceptor opsin"[Title] OR "opsin"[Title] or "rhodopsin"[Title]))) NOT "melanopsin"[Title] NOT "OPN4"[Title] NOT "opn4"[Title] NOT "non-visual"[Title] NOT "nonvisual"[Title] NOT "encephalopsin"[Title] NOT "TMT opsin"[Title] NOT "OPN3"[Title] NOT "opn3"[Title] NOT "neuropsin"[Title] NOT "OPN5"[Title] NOT "opn5"[Title] NOT "pinopsin"[Title] NOT "pinopsin-like"[All-Fields] NOT "parietopsin"[Title] NOT "parapinopsin"[Title] NOT "peropsin"[Title] NOT "peropsin-like"[All-Fields] NOT "RGR opsin"[Title] NOT "VA opsin"[Title] NOT "vertebrate ancient"[Title] NOT "partial"[Title] NOT "fragment"[Title] NOT "exon"[Title] NOT "intron"[Title] NOT "segment"[Title] NOT "kinase"[Title] NOT "GTPase"[Title] NOT "transducin"[Title] NOT "arrestin"[Title] NOT "enhancer"[Title] NOT "promoter"[Title] NOT "regulatory"[Title] NOT "similar"[Title] NOT "homolog"[Title] NOT "putative"[Title] NOT “incomplete”[Title] NOT "transcript variant"[Title] NOT "extra-ocular rhodopsin"[Title] NOT "exorh"[Title] NOT "psuedogene"[Title] NOT "Opsin 4"[Title] NOT "Opsin 5"[Title] NOT "Opsin 6"[Title] NOT "Opsin 7"[Title] NOT "Opsin 8"[Title] NOT "Opsin 9"[Title] NOT "extra-ocular"[Title] NOT "opsin 4"[Title] NOT "opsin 5"[Title] NOT "opsin 6"[Title] NOT "opsin 7"[Title] NOT "opsin 8"[Title] NOT "opsin 9"[Title]'''
    
    #print(f'This is the query being used:\n{final_query}')

    return final_query

from requests.adapters import HTTPAdapter, Retry
import re


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



import requests
from requests.adapters import HTTPAdapter, Retry
import re
import random

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
        