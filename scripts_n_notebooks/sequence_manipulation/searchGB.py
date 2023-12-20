from Bio import Entrez, SeqIO
Entrez.email = "oakley@ucsb.edu"
searchsp = ""

import argparse
import re #regular expressions

ap = argparse.ArgumentParser(description='Converts a 3-column table to LATEX to format species based on columns')
ap.add_argument("-f", "--filein", required=True,
        help="In file in 3-column tab-delimited format")
args = vars(ap.parse_args())
infile = args["filein"]

#read file for data
file1 = open(infile, 'r')
Lines = file1.readlines()

count=0
for line in Lines:
    searchsp = line.strip()

    handle = Entrez.esearch(db="nucleotide", retmax=10, term = searchsp + "[ORGN] rh1", idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    if(int(record["Count"]) == 0) :
        #search again with rhodopsin
        handle = Entrez.esearch(db="nucleotide", retmax=10, term = searchsp + "[ORGN] rhodopsin", idtype="acc")
        record = Entrez.read(handle)
        handle.close()
    if(int(record["Count"]) == 0) :
        print(searchsp + "\tNot found")
    elif(int(record["Count"]) >= 2) :
        print(searchsp + "\tMore than 2 records")
    else :
        #now grab sequence data based on retrieved accession number
        opsacc = record["IdList"][0]
        handle2 = Entrez.efetch(db="nucleotide", id=opsacc, rettype="gb", retmode="text")
        record = SeqIO.read(handle2, "gb")
        handle2.close()
        id = record.id
        dna = record.seq
        organism = record.annotations["organism"]
        orglist = organism.split(" ")
        genus = orglist[0]
        species = orglist[1]

        for i,feature in enumerate(record.features):
             if feature.type=='CDS':
                  aa = feature.qualifiers['translation'][0]

        print(searchsp, end="\t")
        print(genus, end='\t')
        print(species, end='\t')
        print("NCBI", end='\t')
        print(id, end='\t')
      # Comment one of next two lines to print or to not print DNA right now
        print(dna, end='\t')
      #    print(end='\t')
        print(aa, end='\n')
