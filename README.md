**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **DOI**: [![DOI](https://zenodo.org/badge/733656093.svg)](https://zenodo.org/doi/10.5281/zenodo.10667839)


# Visual Physiology Opsin Database (VPOD)
**_VPOD_, a database opsin data and machine-learning models to predict phenotype.**

![](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db/blob/main/scripts_n_notebooks/figure_making/figures/sep_histograms_with_scaled_kde_and_colorbar_1_23_24.png?raw=True) <!-- Alt text: Histogram distributions of Vertebrate and Invertebrate Opsin Light Sensitivity Data - λmax - from VPOD_het_1.0 with a scaled Kernel Density Estimate (KDE) curves overlaid to better visualize the general shape and characteristics of our λmax distributions. -->

  _Histogram distributions of Vertebrate and Invertebrate Opsin Light Sensitivity Data - λmax - from VPOD_het_1.0 with a scaled Kernel Density Estimate (KDE) curves overlaid to better visualize the general shape and characteristics of our λmax distributions_

--- 
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
Instructions for navigating VPOD data files, including raw and curated data used for training models, and for accessing scripts/notebooks used for sequence data manipulation, testing, and model training.

  ### Data
  * Navigate to the folder _vpod_data/VPOD_1.0_ 
  * Select the _formatted_data_subsets_ folder to access subsets of the database suitable for direct model training without requiring mySQL or sequence alignment.
      - The folder _vpod_2023-10-16_12-13-11_ contains all data subsets.
      - Files marked _VPOD_xxx_1.0_ (ex. _VPOD_vert_het_1.0.fasta_) are the fully aligned data subsets 
      - Files marked _xxx_meta_ (ex. _wds_meta.tsv_) are the corresponding metadata files for each subset
  * Select the _raw_database_files_ folder for the raw components of the database taht will need to be loaded into mySQL and formatted using steps 0-2 of the _vpod_main_wf.ipynb_ Jupyter notebook

  ### Scripts & Notebooks
  Instructions for accessing scripts and Jupyter notebooks used to create the database and train ML models.

   ### Scripts
   NOTE - It's recommended that you open these scripts in a compiler for a more detailed explination of how to properly use them* 
   * Navigate to the folder _scripts_n_notebooks/sequence_manipulation_
       - Select the script _mutagenesis.py_ to generate mutants from a sequence accessible via NCBI by accession number or manually enter the raw sequence.
       - Select the script _chimeras.py_ to generate chimeric sequences from sequences accessible via NCBI by accession or manually enter the raw sequence.
           - Opsin _chimeras_ are mutants where one or more transmembrane domains are copied from a different opsin to replace the original. 
           - Contact us for assistance with using this script.
   ### Notebooks
   * Navigate to the folder _scripts_n_notebooks_
      - Select the folder _vpod_ML_workflows_ to access notebooks used for training ML models.
          - The _substests_ folder contains notebooks used for subtests outlined in "_Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD)_".
      - Select the folder _figure_making_ to access the Jupyter notebook _figuremaking.ipynb_ used to generate some of the figures in the publication.
      - Select the folder _phylogenetic_imputation_ to access the R-notebook _Phylogenetic_Imputation.Rmd_ used to compare ML predictions to phylogenetic imputation.
---

### Training New ML Models With _deepBreaks_
Instructions for using VPOD and training ML models with _deepBreaks_.

1. Follow the directions and install _[deepBreaks](https://github.com/omicsEye/deepbreaks)_ using the provided guide.
2. Refer to the _requirements.txt_ provided above and ensure all necessary packages are installed in a dedicated environment (using Conda is recommended).
3. Navigate to _scripts_n_notebooks/vpod_ml_workflows_ and open _vpod_main_wf.ipynb_.
   - To train models using _formatted_data_subsets_, scroll down to **Step 3: deepBreaks** and follow the instructions.
   - To train models using _raw_database_files_, start from the top of the document and follow the instructions.
     
---

## License
All data and code is covered under a GNU General Public License (GPL)(Version 3), in accordance with Open Source Initiative (OSI)-policies

## Citation
IF citing this GitHub and its contents use the following DOI provided by Zenodo...

    10.5281/zenodo.10667840
    
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


