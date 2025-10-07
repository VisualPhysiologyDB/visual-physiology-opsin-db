from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import os
import pandas as pd

def unalign_fasta(input_file, output_file):
    """
    Removes gap characters ('-') from sequences in a FASTA file.

    Args:
        input_file (str): Path to the input aligned FASTA file.
        output_file (str): Path to the output unaligned FASTA file.
    """
    unaligned_records = []
    
    # Use SeqIO.parse for robust FASTA parsing
    for record in SeqIO.parse(input_file, "fasta"):
        # Convert the sequence to a string and remove gaps
        unaligned_seq_str = str(record.seq).replace("-", "")
        
        # Create a new Seq object from the unaligned string
        record.seq = Seq(unaligned_seq_str)
        
        # Add the modified record to our list
        unaligned_records.append(record)
        
    # Write the unaligned records to the output file
    SeqIO.write(unaligned_records, output_file, "fasta")
    print(f"Successfully created unaligned file at: {output_file}")

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
                  try:
                    accession = str(record.description)
                    if len(header) == 3:
                        if header[1] in ['male','female']:
                            species_name = header[0] + ' sp.'
                            opsin_type = header[1] + ' ' + header[2]
                        else:
                            species_name = header[0] + ' ' + header[1]
                            opsin_type = header[2]
                    else:
                        if header[2] in (['male','female']):
                            species_name = header[0] + ' ' + header[1]
                            opsin_type = header[2] + ' ' + header[3]
                        else:
                            species_name = header[0] + ' ' + header[1]
                            opsin_type = ' '.join(header[2:])
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
                  try:
                    accession = str(record.description)
                    species_name = header[0] + " sp."
                    opsin_type = header[1]
                  except:
                    print(header)
                    raise Exception('Header with unique formating seems to have caused an error')
            
      aln_sequence = str(record.seq).replace(" ", "").replace("\n", "").replace("=", "")
      sequence = str(aln_sequence).replace('-','')
      seq_len = len(sequence)
      data.append([species_name, opsin_type, accession, aln_sequence, sequence, seq_len])
        
  df = pd.DataFrame(data, columns=['species_name', 'opsin_type', 'accession', 'aln_sequence','sequence', 'seq_length'])
  return df

def run_blast_query(blasttyp, QUERY_FILE, DB_BASE_PATH, OUTPUT_FILE, EVALUE="1e-5", OUTFMT="6"):
    """
    Runs a BLAST search and returns the results as a pandas DataFrame.
    DB_BASE_PATH should be the base path for the formatted BLAST database (without .fasta or extensions like .pin).
    """
    
    column_names = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]

    blast_cmd = [
        blasttyp,
        "-query", QUERY_FILE,
        "-db", DB_BASE_PATH,    # Database base path to search against
        "-out", f'{OUTPUT_FILE}.tsv',    # Where to save the raw results
        "-evalue", EVALUE,      # E-value threshold
        "-outfmt", OUTFMT       # Output format (e.g., tabular format 6)
    ]

    #print(f"Running {blasttyp} with query '{os.path.basename(QUERY_FILE)}' against DB '{os.path.basename(DB_BASE_PATH)}'...")
    result = subprocess.run(blast_cmd, capture_output=True, text=True)

    if result.stderr:
        print(f"  Error running {blasttyp}: {result.stderr.strip()}")
        return pd.DataFrame()  # Return empty DataFrame on BLAST error

    #print(f"  {blasttyp} search completed. Raw results intended for: {OUTPUT_FILE}")

    if os.path.exists(f"{OUTPUT_FILE}.tsv") and os.path.getsize(f"{OUTPUT_FILE}.tsv") > 0:
        try:
            blast_df = pd.read_csv(f'{OUTPUT_FILE}.tsv', sep='\t', header=None, names=column_names)
            clean_blast_df = clean_blast_query(blast_df, OUTPUT_FILE)
            #print(f"  Successfully read {len(blast_df)} BLAST hits into DataFrame from {OUTPUT_FILE}.")
            return clean_blast_df
        except pd.errors.EmptyDataError:
            print(f"  Warning: BLAST output file {OUTPUT_FILE} is empty.")
            return pd.DataFrame()
        except Exception as e:
            print(f"  Error reading BLAST output file {OUTPUT_FILE} into DataFrame: {e}")
            return pd.DataFrame()
    else:
        if not os.path.exists(OUTPUT_FILE):
            print(f"  Error: Output file {OUTPUT_FILE} was not created by BLAST.")
        elif os.path.getsize(OUTPUT_FILE) == 0:
            print(f"  Warning: BLAST output file {OUTPUT_FILE} is empty (no hits found or error).")
        return pd.DataFrame()
    
    
def clean_blast_query(blast_results_df, output_file):    
    # Ensure 'pident' is numeric for correct identification of maximum
    blast_results_df['pident'] = pd.to_numeric(blast_results_df['pident'], errors='coerce')
    blast_results_df.dropna(subset=['pident'], inplace=True) # Remove rows if pident wasn't numeric

    # Filter to keep only the hit with the highest 'pident' for each 'qseqid'
    # Using .copy() to avoid SettingWithCopyWarning
    idx = blast_results_df.groupby('qseqid', group_keys=False)['pident'].idxmax()
    best_hits_df = blast_results_df.loc[idx].copy()
    
    best_hits_df.to_csv(f'{output_file}_filtered.tsv', index=False)
    
    return best_hits_df


    
