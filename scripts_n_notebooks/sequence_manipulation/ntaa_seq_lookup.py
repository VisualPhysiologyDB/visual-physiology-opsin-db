import argparse
from ast import And
import re #regular expressions
from Bio import Entrez, SeqIO
Entrez.email = 'oakley@ucsb.edu'
from skbio import TabularMSA
from skbio import Protein

ap = argparse.ArgumentParser(description='Call the Nucleotide and Amino Acid Sequences')

ap.add_argument("-ta", "--target_accession", required=True, 
        help="accession number for sequence to change. This should be the DNA accession and protein will be pulled from translation")

args = vars(ap.parse_args())
accession_aa  = args["target_accession"]
accession_nt = args["target_accession"]


def getAcc(accession_aa):

    handle = Entrez.efetch(db="nucleotide", id=accession_aa, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    for i,feature in enumerate(record.features):
         if feature.type=='CDS':
              aa = feature.qualifiers['translation'][0]         
    return(aa)

amino = Protein(getAcc(accession_aa))

def getNT(accession_nt):
    handle = Entrez.efetch(db="nucleotide", id=accession_nt, rettype='gb', retmode="text")
    record = SeqIO.read(handle, 'gb')
    handle.close()
    nt = record.seq 

    return(nt)

nuc = getNT(accession_nt) 

print("The Nucleotide Seequence Is:")
print(nuc)
print("The Amino Acid Sequence Is:")
print(amino)


