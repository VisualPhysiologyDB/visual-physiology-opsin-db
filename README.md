**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **DOI**: [![DOI](https://zenodo.org/badge/733656093.svg)](https://zenodo.org/doi/10.5281/zenodo.10667839)


# Visual Physiology Opsin Database (VPOD)
**_VPOD_: A database of opsin data and machine-learning models to predict λmax phenotypes.**

![](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db/blob/main/scripts_n_notebooks/figure_making/figures/VPOD_1.1/opsin_histogram_with_scaled_kde_and_colorbar_5_24_24.svg?raw=True) <!-- Alt text: Histogram distributions of Vertebrate and Invertebrate Opsin Light Sensitivity Data - λmax - from VPOD_het_1.1 with a scaled Kernel Density Estimate (KDE) curves overlaid to better visualize the general shape and characteristics of our λmax distributions. -->

  _Histogram distributions of Vertebrate and Invertebrate Opsin Light Sensitivity Data - λmax - from VPOD_het_1.1 with a scaled Kernel Density Estimate (KDE) curves overlaid to better visualize the general shape and characteristics of our λmax distributions_

--- 
## Summary
Here, we report a newly compiled database of all heterologously expressed opsin genes with λmax phenotypes called the Visual Physiology Opsin Database (VPOD). VPOD_1.1 contains 1123 unique opsin genotypes and corresponding λmax phenotypes collected across all animals from 90 separate publications. 

We use VPOD data and _[deepBreaks](https://github.com/omicsEye/deepbreaks)_ (an ML tool designed for exploring genotype-phenotype associations) to show regression-based machine learning (ML) models often reliably predict λmax, account for non-additive effects of mutations on function, and identify functionally critical amino acid sites. 


---

# VPOD user guide #

## Table of Contents

1. [Instructions](#Instructions)
2. [License](#license)
3. [Citation](#citation)
4. [Contact](#contact)
5. [Additional Resources](#additional_resources)

## Instructions

### Data, Results and Code Structure
Instructions for navigating VPOD data files, including raw and curated data used for training models, and for accessing scripts/notebooks used for sequence data manipulation, testing, and model training.

  ### Data
  * Navigate to the folder _vpod_data/VPOD_1.X_ (i.e. _vpod_data/VPOD_1.0_ or _vpod_data/VPOD_1.1_)
    * Select the _formatted_data_subsets_ folder to access subsets of the database suitable for direct model training without requiring mySQL or sequence alignment.
      - The folder _vpod_1.1_data_subsets_2024-05-02_ contains all data subsets for VPOD_1.1.
      - Files marked _VPOD_xxx_het_1.1_ (ex. _VPOD_vert_het_1.1.fasta_) are the fully aligned data subsets 
      - Files marked _xxx_meta_ (ex. _wds_meta.tsv_) are the corresponding metadata files for each subset
    * Select the _raw_database_files_ folder for the raw components of the database that you can load into a mySQL database and creat your own formatted dataset using steps 0-2 of the _vpod_main_wf.ipynb_ Jupyter notebook
  
  ### Results
  * Navigate to the folder _results_files_ - All subfolders contain results specific to different subtests outlined in our publication, "_Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD)_" (_Frazer et al. 2024_)
      - _epistasis_pred_test_
      - _full_iter_sample_results_
      - _imputation_tests_
      - _main_model_results_
      - _msp_tests_
      - _perf_v_tds_
      - _sws_ops_prediction_comparison_test_
      - _wt_model_predicting_all_mutants_test_
    
  ### Scripts & Notebooks
  Instructions for accessing scripts and Jupyter notebooks used to create the database and train ML models.

   ### Scripts
   NOTE - It's recommended that you open these scripts in a compiler for a more detailed explination of how to properly use them* 
   * Navigate to the folder _scripts_n_notebooks/sequence_manipulation_
       - _mutagenesis.py_ - Used to generate mutants from a sequence accessible via NCBI by accession number or manually enter the raw sequence.
       - _chimeras.py_ - Used to generate chimeric sequences from sequences accessible via NCBI by accession or manually enter the raw sequence.
           - Opsin _chimeras_ are mutants where one or more transmembrane domains are copied from a different opsin to replace the original. 
           - Contact us for assistance with using this script.
   ### Notebooks
   * Navigate to the folder _scripts_n_notebooks_
      - Select the folder _vpod_ML_workflows_ to access notebooks used for training ML models.
          - _vpod_main_wf.ipynb_ - Our primary notebook, with a full pipeline for everything from creating a local instance of VPOD using mySQL to formatting datasets and training ML models for λmax predictions. 
          - _substests_ folder - Contains notebooks used for subtests outlined in _Frazer et al. 2024_.
            - _vpod_main_wf_msp_iterate.ipynb_ - Used to iteratively train and test models from each dataset on thirty unseen wild-type invertebrate opsins from a separately curated MSP dataset.
            - _vpod_wf_imp_sample_test_iterate.ipynb_ - Used to iteratively train and test models from each dataset on randomly sampled subset of data; for comparison to predictions made by phylogentic imputation.
            - _vpod_wf_iterate_all_sample_t1_ops.ipynb_ 
            - _vpod_wf_iterate_epistasic_muts.ipynb_
            - _vpod_wf_iterate_model_perf_vs_tds.ipynb_
            - _vpod_wf_iterate_subsample.ipynb_
            - _vpod_wf_iterate_sws_comp.ipynb_
            - _vpod_wf_iterate_train_all_subsets.ipynb_
            - _vpod_wf_wt_mut_test.ipynb_
      - Select the folder _figure_making_ to access the Jupyter notebook _figuremaking.ipynb_ used to generate some of the figures used in _Frazer et al. 2024_.
          - _figures_ contains a collection of completed figures and drafts use in _Frazer et al. 2024_ - seperated by version of the database used to generate the figures (i.e VPOD_1.0 or VPOD_1.1)
      - Select the folder _phylogenetic_imputation_ to access all files used to predict λmax via phylogenetic imputation and compare with predictions made by ML, as detailed in  _Frazer et al. 2024_.
          - _Phylogenetic_Imputation.Rmd_ - Used to load tree files and make λmax predictions via phylogenetic imputation [Requires RStudio].
          - _trees_ - Contains all the tree files and λmax meta-data neccessary for predictions via phylogentic imputation.
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

## Additional_Resources 
  **[Here](https://tinyurl.com/mwmjf72f)** is a link to a bibliography of the publications used to build VPOD_1.1
  
  If you know of publications for training opsin ML models not included in the VPOD_1.1 database, please send them to us through **[this form](https://tinyurl.com/29afaxyr)**
