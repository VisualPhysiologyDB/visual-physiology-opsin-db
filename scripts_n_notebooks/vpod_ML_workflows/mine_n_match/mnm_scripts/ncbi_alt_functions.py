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
import requests
from requests.adapters import HTTPAdapter, Retry

from mnm_scripts.mine_n_match_functions import get_species_taxonomy, get_sp_taxon_dict

# It's good practice to set your email and API key at the start.
# The API key allows for up to 10 requests per second.
Entrez.email = "your.email@example.com"  # Please replace with your email
Entrez.api_key = "1efb120056e1cea873ba8d85d6692abd5d09"

def ncbi_fetch_alt_v2(term, ncbi_db="protein", rettype="gb", format="genbank"):
    """
    Fetches sequences from NCBI's databases using an efficient batch method.

    This improved function uses Entrez.epost to upload a list of IDs to the NCBI
    history server, then fetches the corresponding records in large batches. This
    dramatically reduces the number of API calls and speeds up data retrieval.

    Args:
        term (str): The search term to query the NCBI database.
        ncbi_db (str, optional): The NCBI database to search. Defaults to "protein".
        rettype (str, optional): The retrieval type for efetch. Defaults to "gb".
        format (str, optional): The record format. Defaults to "genbank".

    Returns:
        list: A list of Biopython SeqRecord objects.
    """
    queries_dir = "data_sources/queries"
    if not os.path.exists(queries_dir):
        os.makedirs(queries_dir)

    cache_file = os.path.join(queries_dir, "cached_ncbi_queries.json")
    try:
        with open(cache_file, "r") as f:
            cache = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        cache = {}

    # Step 1: Get the list of IDs for the search term
    handle = Entrez.esearch(db=ncbi_db, term=term, usehistory="y")
    search_results = Entrez.read(handle)
    handle.close()
    
    id_list = search_results["IdList"]
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    if count == 0:
        return []

    NCBI_seq = []
    
    # Check cache first for all IDs
    non_cached_ids = [seq_id for seq_id in id_list if seq_id not in cache]
    for seq_id in id_list:
        if seq_id in cache:
            # Load from cache
            seq_record = SeqIO.read(io.StringIO(cache[seq_id]), format)
            NCBI_seq.append(seq_record)

    # Step 2: Fetch non-cached records in batches
    batch_size = 500
    for start in tqdm(range(0, len(non_cached_ids), batch_size), desc="Fetching NCBI Records"):
        end = min(len(non_cached_ids), start + batch_size)
        
        fetch_handle = Entrez.efetch(
            db=ncbi_db,
            id=",".join(non_cached_ids[start:end]),
            rettype=rettype,
            retmode="text",
            api_key=Entrez.api_key,
        )
        
        try:
            records = list(SeqIO.parse(fetch_handle, format))
            for record in records:
                NCBI_seq.append(record)
                # Store new record in cache
                with io.StringIO() as string_handle:
                    SeqIO.write(record, string_handle, format)
                    # Use record.id which often includes version
                    cache[record.id.split('.')[0]] = string_handle.getvalue()
        except Exception as e:
            print(f"Error parsing records: {e}")
        finally:
            fetch_handle.close()

    # Save updated cache to file
    with open(cache_file, "w") as f:
        json.dump(cache, f, indent=4)

    return NCBI_seq


def build_species_query(species_taxon_dict: dict, species: str) -> str:
    """
    Constructs the species portion of the Entrez query string.
    (This function is from your original script and works well.)
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
        sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism] OR "{temp[0]} {temp[2]}"[Title])'
    else:
        sp_for_query = f'"{species}"[Organism]'
    
    return sp_for_query

def build_visual_opsin_query_v2(sp_for_query: str) -> str:
    """
    Constructs a broader NCBI Entrez query for visual opsins to maximize results.
    This query is optimized for the NCBI 'protein' database.

    Args:
        sp_for_query: The species-specific part of the query string.

    Returns:
        A formatted Entrez query string.
    """
    # Broader search terms for visual opsins
    visual_terms = [
        "opsin", "rhodopsin", "photopigment", "visual pigment",
        "RH1", "RH2", "SWS1", "SWS2", "LWS", "MWS",
        "cone opsin", "rod opsin", "photoreceptor"
    ]
    visual_query = " OR ".join([f'"{term}"[Title]' for term in visual_terms])

    # Essential exclusions for non-visual or problematic proteins
    exclusion_terms = [
        "melanopsin", "peropsin", "neuropsin", "encephalopsin", "pinopsin",
        "parietopsin", "VA opsin", "non-visual", "nonvisual", "RGR opsin", "tmt opsin",
        "kinase", "arrestin", "GTPase", "transducin", "regulator", "interactor",
        "partial", "putative", "hypothetical", "predicted", "uncharacterized", "like"
    ]
    exclusion_query = " NOT ".join([f'"{term}"[All Fields]' for term in exclusion_terms])

    # Assemble the final query
    final_query = f'({sp_for_query}) AND ({visual_query}) NOT ({exclusion_query})'
    return final_query


def ncbi_query_to_df_v2(query_list, species_list, species_taxon_dict):
    """
    Converts a list of NCBI query results into a Pandas DataFrame.
    This version is adapted to handle records from the 'protein' database.
    (Based on your original function with modifications.)
    """
    # This function is long and largely similar to your original.
    # The key change is how the protein sequence is retrieved.
    # For a protein record, seq.seq is the protein sequence directly.
    # For a nuccore record, we must get it from the CDS feature translation.

    Accession, Phylum, Subphylum, Class, Genus, Species = [], [], [], [], [], []
    Protein, full_sp_names, gene_des = [], [], []
    
    for query, sp in zip(query_list, species_list):
        for seq in query:
            # --- Taxonomy and Species Name Logic (adapted from your script) ---
            g_s_name = seq.annotations.get("organism", sp).split(' ', 1)
            if len(g_s_name) < 2: g_s_name.append("") # Avoid index error
                
            entry_spe_name = seq.annotations.get("organism", "")
            
            tax_info = species_taxon_dict.get(sp)
            if tax_info:
                Phylum.append(tax_info.get("Phylum", "Unknown"))
                Subphylum.append(tax_info.get("Subphylum", "Unknown"))
                Class.append(tax_info.get("Class", "Unknown"))
            else:
                # Fallback for species not in the initial dict
                Phylum.append("Unknown")
                Subphylum.append("Unknown")
                Class.append("Unknown")

            # --- Protein Sequence Extraction ---
            pro_seq = ""
            # If the record is from the protein DB, seq.seq is the protein
            if str(seq.seq).isalpha():
                 pro_seq = str(seq.seq)
            # If from nuccore, parse features like your original function
            elif seq.features:
                for feature in seq.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        pro_seq = feature.qualifiers['translation'][0]
                        break # Found the protein, no need to check other features

            # Append data to lists
            Accession.append(seq.id)
            Genus.append(g_s_name[0])
            Species.append(g_s_name[1])
            full_sp_names.append(seq.annotations.get("organism", sp))
            gene_des.append(seq.description)
            Protein.append(pro_seq)
            
    df = pd.DataFrame({
        'Accession': Accession,
        'Phylum': Phylum, 'Subphylum': Subphylum, 'Class': Class,
        'Genus': Genus, 'Species': Species, 'Full_Species': full_sp_names,
        'Protein': Protein, 'Gene_Description': gene_des,
    })

    df.drop_duplicates(subset=['Full_Species', 'Protein'], keep='first', inplace=True)
    return df.reset_index(drop=True)


def ncbi_fetch_opsins_v2(species_list, job_label='unnamed', out='unnamed', filter_len_lower=300, filter_len_upper=700):
    """
    Main orchestration function to fetch opsin sequences from NCBI.
    Uses the new, more efficient fetching and query-building functions.
    (This function is a simplified version of your original, focusing on the core task.)
    """
    print("Starting NCBI Opsin Fetch...")
    
    # NOTE: You will need the `get_sp_taxon_dict` function from your original script
    # to build this dictionary. For this example, we assume it exists.
    # taxon_dict = get_sp_taxon_dict(species_list, Entrez.email, ...)
    
    # For demonstration, we'll create a placeholder. Replace with your actual function call.
    print("Building taxon dictionary (placeholder)...")
    species_taxon_dict = {sp: {"Synonyms": []} for sp in species_list}

    query_list = []
    print("Querying NCBI for opsin sequences...")
    for species in tqdm(species_list, desc="Querying Species"):
        sp_for_query = build_species_query(species_taxon_dict, species)
        complete_ncbi_query = build_visual_opsin_query_v2(sp_for_query)
        
        # Using the new, faster fetch function
        ncbi_seqs = ncbi_fetch_alt_v2(term=complete_ncbi_query, ncbi_db="protein")
        query_list.append(ncbi_seqs)

    print("\nFormatting results into DataFrame...")
    ncbi_query_df = ncbi_query_to_df_v2(query_list, species_list, species_taxon_dict)
    
    # Length filtering
    ncbi_query_df['Prot_Len'] = ncbi_query_df['Protein'].str.len()
    ncbi_query_df = ncbi_query_df[
        (ncbi_query_df['Prot_Len'] >= filter_len_lower) & 
        (ncbi_query_df['Prot_Len'] <= filter_len_upper)
    ].drop(columns=['Prot_Len'])

    print(f"NCBI fetch complete. Found {len(ncbi_query_df)} sequences.")
    return ncbi_query_df


def query_uniprot_for_opsins(species_list: list, reviewed=True):
    """
    Queries UniProt for visual opsins for a list of species in a single request.

    Args:
        species_list (list): A list of species names.
        reviewed (bool): If True, searches only for reviewed (Swiss-Prot) entries.

    Returns:
        pandas.DataFrame: A DataFrame with opsin data, formatted like the NCBI output.
    """
    print("Querying UniProt for opsin sequences...")
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    
    # Build query for all species
    species_query = " OR ".join([f'(organism_name:"{sp}")' for sp in species_list])
    
    opsin_query = (
        '(protein_name:opsin OR protein_name:rhodopsin OR gene:opn* OR gene:rh*)'
        ' NOT (name:melanopsin OR name:peropsin OR name:neuropsin OR name:"non-visual")'
    )
    
    full_query = f"({species_query}) AND ({opsin_query})"
    if reviewed:
        full_query += " AND (reviewed:true)"

    params = {
        "query": full_query,
        "format": "tsv",
        "fields": "accession,protein_name,sequence,organism_name,lineage(PHYLUM),lineage(CLASS)",
        "compressed": "true",
    }
    
    try:
        response = requests.get(base_url, params=params, stream=True)
        response.raise_for_status()
        
        # Use pandas to read the gzipped TSV response directly
        df = pd.read_csv(response.raw, compression='gzip', sep='\t')
        
        # Rename columns to match the desired output format
        df.rename(columns={
            "Entry": "Accession",
            "Protein names": "Gene_Description",
            "Sequence": "Protein",
            "Organism": "Full_Species",
            "Lineage (PHYLUM)": "Phylum",
            "Lineage (CLASS)": "Class"
        }, inplace=True)
        
        # Add missing columns
        df['Subphylum'] = "Unknown" # UniProt doesn't easily provide subphylum
        df[['Genus', 'Species']] = df['Full_Species'].str.split(' ', n=1, expand=True)

        # Reorder columns to match
        final_cols = ['Accession', 'Phylum', 'Subphylum', 'Class', 'Genus', 'Species', 
                        'Full_Species', 'Protein', 'Gene_Description']
        df = df[final_cols]
        
        print(f"UniProt query successful. Found {len(df)} sequences.")
        return df

    except requests.exceptions.RequestException as e:
        print(f"An error occurred with UniProt query: {e}")
        print("Response content:", response.text)
        return pd.DataFrame() # Return empty dataframe on error

# --- Example Usage ---
if __name__ == '__main__':
    # List of species to query
    my_species = [
        "Astyanax mexicanus",
        "Danio rerio",
        "Homo sapiens",
        "Gallus gallus",
        "Mus musculus"
    ]

    # --- Option 1: Use the improved NCBI workflow ---
    # Note: You'll need your 'get_sp_taxon_dict' and other helper functions for this to run fully.
    # ncbi_results_df = ncbi_fetch_opsins_v2(species_list=my_species, job_label="test_run")
    # if not ncbi_results_df.empty:
    #     print("\n--- NCBI Results ---")
    #     print(ncbi_results_df.head())
    #     # ncbi_results_df.to_csv("ncbi_opsin_results.csv", index=False)

    # --- Option 2: Use the fast UniProt workflow ---
    uniprot_results_df = query_uniprot_for_opsins(species_list=my_species)
    if not uniprot_results_df.empty:
        print("\n--- UniProt Results ---")
        print(uniprot_results_df.head())
        uniprot_results_df.to_csv("uniprot_opsin_results.csv", index=False)