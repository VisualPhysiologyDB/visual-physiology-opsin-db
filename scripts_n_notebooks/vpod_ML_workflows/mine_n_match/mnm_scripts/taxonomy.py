import json
from Bio import Entrez
from pygbif import species

Entrez.api_key = "1efb120056e1cea873ba8d85d6692abd5d09"

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
    synonyms = []
    taxonomy = {}
    try:
        if use_higher_taxa == False:
            # Search for the species
            handle = Entrez.esearch(db="taxonomy", term=species_name)
        else:
            handle = Entrez.esearch(db="taxonomy", term=higher_taxa)

        record = Entrez.read(handle)
    
        tax_id = record["IdList"][0]
        # Fetch the taxonomy record
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        record = Entrez.read(handle)[0]
        #if the code makes it this far, then there was a successful query
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

    with open(taxon_file, 'w') as f:
        json.dump(sp_taxon_dict, f, indent=4)  # indent for pretty formatting
    return sp_taxon_dict