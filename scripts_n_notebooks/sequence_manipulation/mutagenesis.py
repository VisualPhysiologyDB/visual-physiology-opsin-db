import argparse
import re #regular expressions
from Bio import Entrez, SeqIO
Entrez.email = 'oakley@ucsb.edu'
from skbio import TabularMSA
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from Bio.Align import substitution_matrices

ap = argparse.ArgumentParser(description='Mutagenesis changes protein sequences with mutations named in the standard way')
ap.add_argument("-m", "--mut", required=True, 
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
switch = int(0)
def getAcc(accession, check):
    
    if accession == "manual":

        if check == 0:
            manual = input("Enter Target Sequence: ")
        else:
            manual = input("Enter Reference Sequence: ")
        return(manual)
    
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    for i,feature in enumerate(record.features):
         if feature.type=='CDS':
              aa = feature.qualifiers['translation'][0]
    return(aa)


#Fetch sequences to manipulate and align
wt = Protein(getAcc(accession, switch))
switch+=1
bovine = Protein(getAcc(raccession, switch))


substitution_matrix = substitution_matrices.load("BLOSUM45")

##Simple example for testing
#wt=Protein("ABCDEF")
#bovine=Protein("ABCDEF")

alignment, score, start_end_positions = global_pairwise_align_protein(bovine, wt, substitution_matrix=substitution_matrix)
dic = alignment.to_dict()
aligned_bovine = dic[0]
aligned_wt = dic[1]
print(aligned_wt)

#Check if single or multiple mutations entered
mutsuite = mutation
if "," in mutation :
    if v == 'yes':
        print("Multiple mutations detected")
    mutations = mutation.split(",")
else :
    if v == 'yes':
        print("Single mutation")
    #declare empty list and append single mutation so later can be used as list if only one mutation
    mutations = []
    mutations.append(mutation)

mutated = aligned_wt #need to copy original sequence and accumulate mutations
for mutation in mutations:
    if v == 'yes':
        print("\n\n" + mutation)
        #print(alignment)

    #check that mutation is formatted properly
    MutationRegex = re.compile(r'([A-Z])(\d+)([A-Z])')
    if re.match(MutationRegex, mutation):
         mu = MutationRegex.search(mutation)
         old = mu.group(1)
         mutsite = int(mu.group(2))
         new = mu.group(3)
    else :
         print("ERROR: Expecting mutation name in the format of XaY where X=old amino acid a=number Y=new amino acid")
         quit()

    #substract one to account for index starting at zero
    mutsite = mutsite-1

    #Count gaps to not count them when finding original site to change
    gaps = aligned_bovine[:mutsite].count('-')

    match = str(mutated[mutsite+gaps:mutsite+gaps+1])

    if v == 'yes':
        refmatch = str(aligned_bovine[mutsite+gaps:mutsite+gaps+1])
        print("AA in reference sequence " + refmatch)
        region = str(aligned_bovine[mutsite+gaps-5:mutsite+gaps+5] )
        print("Region is:" + region)
        print("AA in wt sequence " + match)
        print("Old AA in mutation description " + old)
        print(" ")
        print("Reference: " + str(aligned_bovine))
        print(" ")
        print("NewSeq___: " + str(mutated))

    if match != old :
         print("ERROR: The mutation label does not match the sequence")
         quit()
    else: 
         if v == 'yes':
             print("Changing: " + str(aligned_wt[mutsite+gaps:mutsite+gaps+1]))
             print("   To: " + new)
         mutated = str(mutated[:mutsite+gaps]) + new + str(mutated[mutsite+gaps+1:]) 

mutated = mutated.replace('-','')
newacc = accession + "_" + mutsuite
print(newacc, end='\t\t')
print(mutated)


