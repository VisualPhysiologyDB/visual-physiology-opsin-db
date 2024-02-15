import argparse
from posixpath import splitext
import re #regular expressions
from Bio import Entrez, SeqIO
from skbio import TabularMSA
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from Bio.Align import substitution_matrices

ap = argparse.ArgumentParser(description='Mutagenesis changes protein sequences with mutations named in the standard way')

ap.add_argument("-co","--chimera_opsins", required=True,
        help="Chimeric protein name in format  Acc#1_AARange|Acc#2_AARange")
# Naming convention - DX1234_1_306-BV1234_307_348

ap.add_argument("-ra", "--reference_accession", required=False, default = "NM_001014890", 
        help="accession number for reference sequence numbering. This should be the DNA accession and protein will be pulled from translation. Default is Bos taurus rh1")

ap.add_argument("-v", "--verbose", required=False, default="no", 
        help="-v yes will print more information")

ap.add_argument("-em", "--email", required=True,
        help = "enter school email address")

args = vars(ap.parse_args())
raccession = args["reference_accession"]
v = args["verbose"]
cops = args["chimera_opsins"]
bovine_seq = "NM_001014890"
Entrez.email = args["email"]

chunks = cops.split('-')
l = len(chunks)
    
i = 0
    
accessions = []
    
for acc in range(len(chunks)):
    
    y = chunks[i].split('_')

    accessions.extend(y[0:1])
    
    if accessions[i] == "NM":

        accessions[i] = input("Accession Starting With 'NM' Detected! \n Please Re-enter Accession Here: ")
       
    i += 1

print(accessions)

i = 0
    
ranges = []
    
for chunk in range(len(chunks)):
    
    y = chunks[i].split('_')
    
    n = y[1]

    k = 0

    z = 0

    if len(n) > 3 :

        k = y[2]

        k = int(k) - 1

        z = y[3]

        z = int(z) 

    else:

        k = y[1]

        k = int(k)-1

        z = y[2]

        z = int(z) 

    ranges.append(k)

    ranges.append(z)

    i += 1

print(ranges)

def getAcc(aa_seq):

    if aa_seq == "swsanc1":
        return("MSKMSEEEDFYLFGNISSVSPFEGPQYHLAPKWAFYLQAAFMGFVFFVGTPLNAIVLFVTVKYKKLRQPLNYILVNISLGGFLFCIFSVSTVFFSSLRGYFVFGHTVCALEAFLGSVAGLVTGWSLAVLAFERYIVICKPFGNFKFGSKHALMAVVLTWIIGIGCSTPPFFGWSRYIPEGLQCSCGPDWYTVNTEYNSESYTWFLFIFCFIIPLSIITFSYSQLLGALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCVCYVPYAAMALYMVNNRNHGLDLRLVTIPAFFSKSSCVYNPIIYAFMNKQFRACIMETVCGKPMSDDSDVSSSSQKTEVSSVSSSQVSPS")
    if aa_seq == "swsanc2":
        return("MSKMSEEEDFYLFKNISSVGPWDGPQYHIAPKWAFYLQAAFMGFVFFVGTPLNAIVLIVTVKYKKLRQPLNYILVNISLGGFLFCIFSVFTVFVSSSQGYFVFGRTVCALEAFLGSVAGLVTGWSLAFLAFERYIVICKPFGNFRFSSKHALMVVVATWIIGIGVSIPPFFGWSRYIPEGLQCSCGPDWYTVGTKYKSEYYTWFLFIFCFIIPLSLICFSYSQLLGALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCLCYVPYAAMAMYMVNNRNHGLDLRLVTIPAFFSKSSCVYNPIIYSFMNKQFRACIMETVCGKPMTDDSDVSSSSQKTEVSSVSSSQVSPS")
    if aa_seq == "swsanc3":
        return("MSKMSEEEDFYLFKNISNVGPWDGPQYHIAPKWAFYLQAAFMGFVFFVGTPLNAIVLIVTVKYKKLRQPLNYILVNISLGGFLFCIFSVFTVFVSSSQGYFVFGRTVCALEAFLGSVAGLVTGWSLAFLAFERYIVICKPMGNFRFSSKHALMVVVATWIIGIGVSIPPFFGWSRYIPEGLQCSCGPDWYTVGTKYKSEYYTWFLFIFCFIIPLSLICFSYSQLLGALRAVAAQQQESATTQKAEREVSRMVIVMVGSFCLCYVPYAAMAMYMVNNRNHGLDLRLVTIPAFFSKSSCVYNPIIYSFMNKQFRACIMETVCGKPMSDDSSVSSSSQKTEVSSVSSSQVSPS")
    if aa_seq == "swsanc4":
        return("MSKMSEEEDFYLFKNISSVGPWDGPQYHIAPMWAFYLQAAFMGFVFFVGTPLNAIVLIVTVKYKKLRQPLNYILVNISLGGFLFCIFSVFTVFVASSQGYFVFGRHVCALEAFLGSVAGLVTGWSLAFLAFERYIVICKPFGNFRFSSKHALMVVVATWIIGIGVSIPPFFGWSRYIPEGLQCSCGPDWYTVGTKYKSEYYTWFLFIFCFIVPLSLICFSYSQLLGALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCLCYVPYAALAMYMVNNRNHGLDLRLVTIPAFFSKSSCVYNPIIYCFMNKQFRACIMETVCGKPMTDDSDVSSSAQKTEVSSVSSSQVSPS")
    if aa_seq == "swsanc5":
        return("MSKMSEEEDFYLFKNISSVGPWDGPQYHIAPMWAFYLQTAFMGFVFFVGTPLNAIVLIVTVKYKKLRQPLNYILVNISFAGFLFCIFSVFTVFVASSQGYFVFGRHVCALEAFLGSVAGLVTGWSLAFLAFERYIVICKPFGNFRFNSKHALLVVVATWIIGIGVSIPPFFGWSRYIPEGLQCSCGPDWYTVGTKYKSEYYTWFLFIFCFIVPLSLIIFSYSQLLGALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCLCYVPYAALAMYMVNNRDHGLDLRLVTIPAFFSKSSCVYNPIIYCFMNKQFRACIMETVCGKPMTDDSDVSSSAQKTEVSSVSSSQVSPS")
    if aa_seq == "swsanc6":
        return("MSKMSEEEDFYLFKNISSVGPWDGPQYHIAPMWAFYLQTAFMGFVFVVGTPLNAIVLIVTVKYKKLRQPLNYILVNISFSGFISCIFSVFTVFVASSQGYFVFGKHVCALEAFVGATGGLVTGWSLAFLAFERYIVICKPFGNFRFNSKHALLVVVATWIIGVGVAIPPFFGWSRYIPEGLQCSCGPDWYTVGTKYKSEYYTWFLFIFCFIVPLSLIIFSYSQLLSALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCLCYVPYAALAMYMVNNRDHGLDLRLVTIPAFFSKSSCVYNPIIYCFMNKQFRACIMETVCGKPMTDDSDVSSSAQRTEVSSVSSSQVSPS")
    if aa_seq == "swsanc7":
        return("MSKMSEEEDFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGFVFFVGTPLNAIVLVATLRYKKLRQPLNYILVNVSLGGFLFCIFSVFTVFIASCHGYFVFGRHVCALEAFLGSVAGLVTGWSLAFLAFERYIVICKPFGNFRFSSKHALMVVVATWIIGIGVSIPPFFGWSRFIPEGLQCSCGPDWYTVGTKYRSEYYTWFLFIFCFIVPLSLICFSYSQLLRALRAVAAQQQESATTQKAEREVSHMVVVMVGSFCLCYVPYAALAMYMVNNRNHGLDLRLVTIPAFFSKSSCVYNPIIYCFMNKQFRACIMEMVCGKPMTDESDVSSSAQKTEVSTVSSSQVGPN")
    if aa_seq == "manual":
        manual = input("Manual Sequence Request Detected! \nEnter Sequence Here: ")
        return(manual)
    
    handle = Entrez.efetch(db="nucleotide", id=aa_seq, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    for i,feature in enumerate(record.features):
        if feature.type=='CDS':
            aa = feature.qualifiers['translation'][0]

    return(aa)

#Fetch sequences to manipulate and align
if raccession == "manual":
    refacc = Protein(input("Manual Reference Sequence Request Detected! \nEnter Sequence Here: "))
else: 
    refacc = Protein(getAcc(raccession))

bovine = Protein(getAcc(bovine_seq))

substitution_matrix = substitution_matrices.load("BLOSUM45")

##Simple example for testing
#wt=Protein("ABCDEF")
#bovine=Protein("ABCDEF")

j = 0

k = 0

chimera_opsin = ""


for entries in range(len(accessions)):

    mut_check = accessions[j]

    if "[" in mut_check :

        pl_hld = mut_check.split('[')

        accessions[j] = pl_hld[0:1]

    wt = Protein(getAcc(accessions[j]))

    alignment, score, start_end_positions = global_pairwise_align_protein(refacc, wt, substitution_matrix=substitution_matrix)
    dic = alignment.to_dict()
    aligned_wt = dic[1]
    
    if ranges[k+1] == len(aligned_wt):

        temp = aligned_wt[ranges[k]:]
        chimera_opsin += str(temp)

    else:
        temp = aligned_wt[ranges[k]: ranges[k+1]]

        chimera_opsin += str(temp)

    k += 2

    j += 1



if "[" in cops :

    i = 0

    chimera_opsin = chimera_opsin.replace('-','')
    chimera_protein = Protein(chimera_opsin)
    alignment, score, start_end_positions = global_pairwise_align_protein(refacc, chimera_protein, substitution_matrix=substitution_matrix)
    dic = alignment.to_dict()
    aligned_bovine = str(dic[0])
    chimera_opsin = str(dic[1])

    pl_hld = cops.split('[')

    temp = pl_hld[1]

    pl_hld_two = temp.split(']')

    pl_hld_three = str(pl_hld_two[0:1])

    print(pl_hld_three)

    pl_hld_four = pl_hld_three.split(',')

    print(pl_hld_four)

    for mutations in pl_hld_four:

        mut = str(pl_hld_four[i]) 

        mut = mut.replace("'","")
        mut = mut.replace("[","")
        mut = mut.replace("]","")

        print(mut)

        if len(mut) == 4 :

            a = int(mut[1:3]) - 1
            b = mut[4]
            gaps = aligned_bovine[:a].count('-')
            chimera_opsin = chimera_opsin[0:a+gaps] + b + chimera_opsin[a+gaps+1:]
        

        else:

            a = int(mut[1:4]) - 1
            b = mut[4]
            gaps = aligned_bovine[:a].count('-')
            chimera_opsin = chimera_opsin[0:a+gaps] + b + chimera_opsin[a+gaps+1:]

        i+=1
    

new_chimera = chimera_opsin.replace('-','')

print(new_chimera)






