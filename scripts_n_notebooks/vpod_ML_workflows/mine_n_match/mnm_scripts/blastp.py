import subprocess
import os
import tempfile
import pandas as pd
import json
import io
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --- Cache Management Functions ---

def get_cache_path(database_path, wrk_dir="."):
    """Generates a unique cache filename based on the BLAST database name."""
    db_name = Path(database_path).stem
    cache_dir = Path(wrk_dir)/ "data_sources/queries"
    cache_dir.mkdir(exist_ok=True)
    return cache_dir / f"blast_cache_{db_name}.json"

def load_blast_cache(cache_path):
    """Loads the BLAST results cache from a JSON file."""
    if cache_path.exists():
        try:
            with open(cache_path, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError:
            print(f"Warning: Cache file at {cache_path} is corrupted. Starting fresh.")
            return {}
    return {}

def save_blast_cache(cache, cache_path):
    """Saves the BLAST results cache to a JSON file."""
    with open(cache_path, 'w') as f:
        json.dump(cache, f, indent=4)

# --- Core Bioinformatic Functions ---

def run_blast_on_the_fly(blast_type, query_file, db_path, evalue, outfmt_str):
    """
    Performs a BLAST search and captures the standard output.
    This is preferred for caching as it avoids writing intermediate files.
    """
    print(f"INFO: Running {blast_type} for new sequences...")
    try:
        blast_cmd = [
            blast_type,
            "-query", query_file,
            "-db", db_path,
            "-evalue", str(evalue),
            "-outfmt", outfmt_str
        ]
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
        return result.stdout
    except FileNotFoundError:
        print(f"ERROR: '{blast_type}' command not found. Make sure BLAST+ is installed and in your system's PATH.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {blast_type} failed with the following error:\n{e.stderr}")
        return None

def filter_best_blast_hits(blast_output_string, column_names):
    """
    Parses BLAST output string and filters for the single best hit per query.
    The best hit is determined by the lowest e-value, with bitscore as a tie-breaker.
    """
    if not blast_output_string or not blast_output_string.strip():
        print("WARNING: BLAST returned no hits.")
        return pd.DataFrame()

    try:
        blast_results_df = pd.read_csv(io.StringIO(blast_output_string), sep="\t", header=None, names=column_names)
    except Exception as e:
        print(f"ERROR: Failed to parse BLAST output. {e}")
        return pd.DataFrame()

    # Sort by e-value (ascending) and bitscore (descending) to rank hits
    blast_results_df.sort_values(by=['evalue', 'bitscore'], ascending=[True, False], inplace=True)
    
    # Keep only the first occurrence of each query ID (which is now the best hit)
    best_hits_df = blast_results_df.drop_duplicates(subset='qseqid', keep='first').reset_index(drop=True)
    
    return best_hits_df

# --- Main Caching Workflow Function ---
def run_blastp_analysis(blast_type, query_file, db_base_path, output_file, evalue="1e-5", wrk_dir="."):
    """
    Runs a BLAST search with a caching system, correctly handling duplicate sequences and preserving input order.
    
    1. Reads the query FASTA file, preserving the original order of all records.
    2. Identifies unique sequences that are not in the cache.
    3. Runs BLAST only on these new, unique sequences.
    4. Updates the cache with the new results (best hit per sequence).
    5. Reconstructs the final results by iterating through the original input records, looking up
       each sequence in the cache, and applying the original record ID.
    """
    # --- 1. Setup ---
    outfmt_str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    column_names = outfmt_str.split()[1:]
    
    cache_path = get_cache_path(db_base_path, wrk_dir)
    blast_cache = load_blast_cache(cache_path)
    
    # --- 2. Read query file and identify uncached sequences ---
    try:
        all_query_records = list(SeqIO.parse(query_file, "fasta"))
    except FileNotFoundError:
        print(f"ERROR: Query file not found at {query_file}")
        return pd.DataFrame()

    uncached_records_map = {} # Use a dict to store only unique new sequences
    for record in all_query_records:
        seq_str = str(record.seq)
        if seq_str not in blast_cache:
            # The key is the sequence, the value is the record. If a sequence appears
            # multiple times before caching, only the first record is used for the BLAST query.
            if seq_str not in uncached_records_map:
                uncached_records_map[seq_str] = record
            
    uncached_records = list(uncached_records_map.values())
    print(f"INFO: Total records in input: {len(all_query_records)}. Unique uncached sequences to process: {len(uncached_records)}.")

    # --- 3. Process uncached sequences and update cache ---
    if uncached_records:
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tfile:
            temp_query_path = tfile.name
            SeqIO.write(uncached_records, tfile, 'fasta')
        
        blast_output = run_blast_on_the_fly(blast_type, temp_query_path, db_base_path, evalue, outfmt_str)
        os.remove(temp_query_path)
        
        if blast_output:
            best_hits_df = filter_best_blast_hits(blast_output, column_names)
            
            if not best_hits_df.empty:
                newly_processed_results = best_hits_df.to_dict('records')
                
                # Create a map of query ID back to sequence string for efficient cache update
                id_to_seq_map = {rec.id: str(rec.seq) for rec in uncached_records}
                for result_dict in newly_processed_results:
                    query_id = result_dict.get('qseqid')
                    original_sequence = id_to_seq_map.get(query_id)
                    if original_sequence:
                        blast_cache[original_sequence] = result_dict
    
    save_blast_cache(blast_cache, cache_path)
    
    # --- 4. Reconstruct final DataFrame in the original input order ---
    final_ordered_results = []
    queries_no_hits = []
    for record in all_query_records:
        seq_str = str(record.seq)
        result_from_cache = blast_cache.get(seq_str)
        
        if result_from_cache:
            # Found a hit in the cache, now apply the correct ID from the input file
            final_result = result_from_cache.copy()
            final_result['qseqid'] = record.id
            final_ordered_results.append(final_result)
        else:
            # This sequence had no BLAST hit, record it as such
            queries_no_hits.append(record.id)

    if queries_no_hits:
        print(f"WARNING: No BLAST hits found for {len(queries_no_hits)} queries: {', '.join(queries_no_hits[:5])}{'...' if len(queries_no_hits) > 5 else ''}")
            
    if not final_ordered_results:
        print("No results to compile (no hits found for any query).")
        return pd.DataFrame()
        
    final_df = pd.DataFrame(final_ordered_results)
    
    # Ensure all columns from the blast output are present, even if no hits were found
    for col in column_names:
        if col not in final_df.columns:
            final_df[col] = None

    # Reorder columns to match the standard BLAST output format
    final_df = final_df[column_names]

    output_path = f'{output_file}_filtered.csv'
    final_df.to_csv(output_path, index=False)
    print(f"\n✅ Analysis complete. Filtered BLAST results saved to '{output_path}'")
    
    return final_df