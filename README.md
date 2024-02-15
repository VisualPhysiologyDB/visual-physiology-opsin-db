**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) 

# Visual Physiology Opsin Database (VPOD)
_VPOD_, a database opsin data and machine-learning models to predict phenotype.

## Summary
Here, we report a newly compiled database of all heterologously expressed opsin genes with λmax phenotypes called the Visual Physiology Opsin Database (VPOD). VPOD_1.0 contains 864 unique opsin genotypes and corresponding λmax phenotypes collected across all animals from 73 separate publications. 

We use VPOD data and _[deepBreaks](https://github.com/omicsEye/deepbreaks)_(an ML tool designed for exploring genotype-phenotype associations) to show regression-based machine learning (ML) models often reliably predict λmax, account for non-additive effects of mutations on function, and identify functionally critical amino acid sites. 

---

# VPOD user guide #

## Table of Contents

1. [Instructions](#Instructions)
2. [License](#license)
3. [Citation](#citation)
4. [Contact](#contact)

## Instructions

### Data and Code Structure
Instructions for navigating data VPOD files, including raw and curated data used to train models, and for navigating scripts/notebooks used to enter sequence data, perform tests, or train models.
  ### Data
  * Navigate to _vpod_data/VPOD_1.0_
  * Select _formatted_data_subsets_ folder for subsets of the data which can be directly used to train models without the need for using mySQL or sequence alignment
      - _vpod_2023-10-16_12-13-11_ contains the all data subsets
      - Those files marked _VPOD_xxx_1.0_ (ex. _VPOD_vert_het_1.0.fasta_) are the fully aligned data subsets 
      - Those files marked _xxx_meta_ (ex. _wds_meta.tsv_) are the corresponding metadata files for each subset
  * Select _raw_database_files_ folder for subsets of the data which will need to be loaded into mySQL and formatted using steps 0-2 of the _vpod_main_wf.ipynb_ jupyter notebook

  ### Scripts & Notebooks
  Instructions for... 
   ### Scripts

   ### Notebooks
   * Navigate to _scripts_n_notebooks_ and select 
   * dsfsfsd
   * sdfsd

---

### Training New ML Models With _deepBreaks_
Instructions for using VPOD and training ML models with _deepBreaks_.
  1. Follow directions and install _[deepBreaks](https://github.com/omicsEye/deepbreaks)_ using the guide provided on their GitHub.
  2. Refer to the _requirements.txt_ provided above and ensure all neccessary pakcages are downloaded on a dedicated enivronment (using Conda is recommended)
  3. Navigate to _scripts_n_notebooks/vpod_ml_workflows_ and open _vpod_main_wf.ipynb_
     - IF you want to train models using _formatted_data_subsets_ Scroll down to **Step 3: deepBreaks** and follow the markdown notes from there
     - IF you want to train models using _raw_database_files_ start from top of document and follow markdown notes
     

---

## License
All data and code is covered under a GNU General Public License (GPL)(Version 3), in accordance with Open Source Initiative (OSI)-policies

## Citation
IF citing this GitHub and its contents use the following DOI provided by Zenodo...

    place-holder
    
IF citing the paper _"Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD)"_ use the following citation...

    Seth A. Frazer, Mahdi Baghbanzadeh, Ali Rahnavard, Keith A. Crandall, & Todd H Oakley. (2024). Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD). bioRxiv, 2024.02.12.579993. https://doi.org/10.1101/2024.02.12.579993 [pre-print]
    
## Contact
Contact information for author questions or feedback.

**Todd H. Oakley**  
  * oakley@ucsb.edu
  * https://orcid.org/0000-0002-4478-915X

**Seth A. Frazer** 
  * sethfrazer@ucsb.edu
  * https://orcid.org/0000-0002-3800-212X


