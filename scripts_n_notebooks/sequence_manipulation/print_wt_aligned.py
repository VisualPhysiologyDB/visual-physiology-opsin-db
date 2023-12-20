import argparse
import re #regular expressions
from Bio import Entrez, SeqIO
Entrez.email = 'oakley@ucsb.edu'
from skbio import TabularMSA
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from Bio.Align import substitution_matrices

ap = argparse.ArgumentParser(description='Mutagenesis changes protein sequences with mutations named in the standard way')
ap.add_argument("-m", "--mut", required=False, 
        help="mutation name in the format of XaY where X=old amino acid a=number Y=new amino acid. Must be capital letters")
ap.add_argument("-ra", "--reference_accession", required=False, default = "NM_001014890", 
        help="accession number for reference sequence numbering. This should be the DNA accession and protein will be pulled from translation. Default is Bos taurus rh1")
ap.add_argument("-ta", "--target_accession", required=True, 
        help="accession number for sequence to change. This should be the DNA accession and protein will be pulled from translation")
ap.add_argument("-v", "--verbose", required=False, default="no", 
        help="-v yes will print more information")
args = vars(ap.parse_args())
mutation = args["mut"]
accession  = args["target_accession"]
raccession = args["reference_accession"]
v = args["verbose"]

def getAcc(accession):

    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    for i,feature in enumerate(record.features):
         if feature.type=='CDS':
              aa = feature.qualifiers['translation'][0]
    return(aa)


#Fetch sequences to manipulate and align
wt = Protein(getAcc(accession))
bovine = Protein(getAcc(raccession))


substitution_matrix = substitution_matrices.load("BLOSUM45")

##Simple example for testing
#wt=Protein("ABCDEF")
#bovine=Protein("ABCDEF")

alignment, score, start_end_positions = global_pairwise_align_protein(bovine, wt, substitution_matrix=substitution_matrix)
dic = alignment.to_dict()
aligned_bovine = dic[0]
aligned_wt = dic[1]

print(wt)

print(aligned_wt)

print(aligned_bovine)

print(dic)