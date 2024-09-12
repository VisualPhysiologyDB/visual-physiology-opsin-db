# prediction_functions.py
import subprocess 
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
import pandas as pd
import tempfile
import argparse
import os
import csv  # For CSV export


def process_sequence(mafft, sequence, selected_model, alignment_data, gap_threshold = 0.6):

    if sequence == None:
        return Exception('Error: No sequence given')
    if selected_model == None:
        return Exception('Error: No model selected')
    if alignment_data == None:
        return Exception('Error: No alignment data selected')

    temp_seq = "./vpod_scripts/temp/temp_seq.fasta"
    with open(temp_seq, "w") as temp_file:  # Key change
        if '>' in sequence:
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence)
        else:
            sequence = ">placeholder_name\n" + sequence
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence) # Write your data to the file object
        

    new_ali = './vpod_scripts/temp/temp_ali.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    mafft_exe = mafft #change to your own directory for mafft.bat or mafft execution file
    seq_type = 'aa'

    cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
    with open(new_ali, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

    new_seq_test = read_data(new_ali, seq_type = seq_type, is_main=True, gap_threshold=gap_threshold)
    ref_copy = read_data(alignment_data, seq_type = seq_type, is_main=True, gap_threshold=gap_threshold)
    last_seq = ref_copy.shape[0]
    new_seq_test = new_seq_test.iloc[last_seq:].copy()
    #print(new_seq_test)

    # ... (Load the selected model and make a prediction)
    load_top_mod = load_obj(selected_model)
    prediction = load_top_mod.predict(new_seq_test)
    
    return(prediction[0])
 

def process_sequences_from_file(mafft,input_file,output_file,selected_model,alignment_data, gap_threshold = 0.6):
    if input_file == None:
        return Exception('Error: No input file given')
    if output_file == None:
        output_file = 'opsin_predictions.tsv'

    with open(input_file, 'r') as f:
        sequences = []
        names = []
        i = 0
        line_count = 0
        entry = ""
        lines = f.readlines()
        num_lines = len(lines)
        #print(num_lines)

        for line in lines:
            if '>' in line:
                if i == 1:
                    names.append(line.replace('>','').strip())
                    sequences.append(entry)
                    entry = ""
                    entry += line
                    #print(sequences)
                    line_count+=1
                else:
                    names.append(line.replace('>','').strip())
                    entry += line
                    i+=1
                    line_count+=1
            else:
                entry += line
                line_count+=1
                if line_count >= num_lines:
                     sequences.append(entry)
    #print(sequences)
    #print(names)
                     
    predictions = []
    print('Processing predictions...')
    for seq in sequences:
        #print(seq)
        prediction = process_sequence(mafft,seq,selected_model, alignment_data, gap_threshold)  # Process each sequence
        predictions.append(prediction)
    #print(predictions)

    with open(output_file, 'w') as f:
        i = 0
        while i in range(len(names)):
            if i == 0:
                f.write('Names\tPredictions\n')
            f.write(f"{names[i]}\t{predictions[i]}\n")
            print(f"{names[i]}\t{predictions[i]}\n")
            i+=1
    return(print(f'All {i} Predictions Complete'))
