**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **VPOD_1.1 DOI**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12213246.svg)](https://doi.org/10.5281/zenodo.12213246)


# Visual Physiology Opsin Database (VPOD) - Version 1.3
**_VPOD_: A database of opsins and machine-learning models to predict λmax phenotypes.**

![](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db/blob/main/scripts_n_notebooks/figure_making/figures/VPOD_1.1/opsin_histogram_with_scaled_kde_and_colorbar_5_24_24.svg?raw=True) <!-- Alt text: Histogram distributions of Vertebrate and Invertebrate Opsin Light Sensitivity Data - λmax - from VPOD_het_1.1 with a scaled Kernel Density Estimate (KDE) curves overlaid to better visualize the general shape and characteristics of our λmax distributions. -->

  _Histogram distributions of Vertebrate and Invertebrate Opsin Light Sensitivity Data - λmax - from VPOD_het_1.1 with a scaled Kernel Density Estimate (KDE) curves overlaid to better visualize the general shape and characteristics of our λmax distributions_

--- 
## Summary
* We introduce the Visual Physiology Opsin Database, a newly compiled database for all heterologously expressed, and a partial collection of physiologically inferred, [opsin genes](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2005-6-3-213) with λmax phenotypes (wavelength of maximal absorbance; peak-senstivity). 

* **VPOD_1.3** contains 1714 unique opsin genotypes and corresponding λmax phenotypes collected across all animals from 120+ separate publications. 

* We use VPOD data and _[deepBreaks](https://github.com/omicsEye/deepbreaks)_ (an ML tool designed for exploring genotype-phenotype associations) to show regression-based machine learning (ML) models often reliably predict λmax, account for non-additive effects of mutations on function, and identify functionally critical amino acid sites. 

* We provide an approach that lays the groundwork for future robust exploration of molecular-evolutionary patterns governing phenotype, with potential broader applications to any family of genes with quantifiable and comparable phenotypes.

---

# VPOD user guide #

## Developer Notes
* **UPDATE!!!- Welcome to VPOD_1.3!** The latest release of the Visual Physiology Opsin Database. We have quite a bit of new data, results and tools to check-out.
* **Want an user-friendly way to predict opsin λmax with our models?**
  * Visit our GitHub for _[OPTICS](https://github.com/VisualPhysiologyDB/optics)_, an open-source tool, avaiable online or via command line, that predicts opsin phenotypes (λmax)
  * Read our most recent publication, detailing the creation of VPOD_1.3 and OPTICS - _[Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map](https://doi.org/10.1101/2025.08.22.671864)_ [pre-print]

## Table of Contents

1. [Instructions](#Instructions)
   * [Data](#data)
   * [Results](#results)
   * [Scripts & Notebooks](#scripts--notebooks)
   * [Training New ML Models With _deepBreaks_](#training-new-ml-models-with-deepbreaks)
3. [License](#license)
4. [Citation](#citation)
5. [Contact](#contact)
6. [Additional Resources](#additional_resources)

## Instructions

### Data, Results and Code Structure
Instructions for navigating VPOD data files, including raw and curated data used for training models, and for accessing scripts/notebooks used for sequence data manipulation, testing, and model training.

  ### Data
  * Navigate to the folder **_vpod_data/VPOD_1.X_** (i.e. _vpod_data/VPOD_1.1_ or _vpod_data/VPOD_1.3_)
    * Select the **_formatted_data_subsets_** folder to access subsets of the database suitable for direct model training without requiring mySQL or sequence alignment.
      - The folder _vpod_1.3_data_splits_2025-02-28_15-51-04_ contains all data subsets for VPOD_1.3.
      - Files marked _xxx.txt_ (ex. _vert.txt_) are the unaligned data subsets.
      - Files marked _xxx_aligned.txt_ (ex. _vert.txt_) are the aligned data subsets (not-formatted; will become obsolete in later versions of VPOD).
      - Files marked _VPOD_xxx_het_1.1_ (ex. _VPOD_vert_het_1.1.fasta_) are the fully aligned and formatted data subsets.
      - Files marked _xxx_meta_ (ex. _wds_meta.tsv_) are the corresponding metadata files for each subset (includes species, accession, gene names, λmax, etc).
    * Select the **_raw_database_files_** folder for the raw components of the database that you can load into a SQLite database and create your own formatted dataset using steps 0-2 outlined in the _vpod_main_wf.ipynb_ Jupyter notebook (for more information on the specific meaning of data columns, see _Frazer et al. 2024_).
      - _litsearch.csv_ - All literature search information relevant to the current version of VPOD.  
      - _references.csv_ - All publication references relevant to the current version of VPOD.
      - _opsins.csv_ - All opsin sequence data and taxonomic meta-data.
      - _heterologous.csv_ - All opsin phenotype data (λmax) and experiment related meta-data. 

---
  
  ### Results
  * Navigate to the folder _results_files_ - All subfolders contain results specific to different subtests outlined in _Frazer et al. 2024_ (for VPOD_1.0 and VPOD_1.1) or _Frazer and Oakley_ (for VPOD_1.3).
      - **_aa_prop_combos_** - Contains results for the combinatorial test to determine optimal, dataset specific, set of amino-acid properties to encode opsin amino acid sequences by [VPOD_1.3 only].
      - **_aa_prop_vs_one_hot_** - Contains comparisons of predictive power between ML models trained using either one-hot or amino-acid property encoding [VPOD_1.3 only].
      - **_epistasis_pred_test_** - Contains predictions by and comparisons between our WT and WDS models on 111 'epistatic mutations' (non-additive) to more generally the capabilities of our ML models to predict epistatic interactions between mutations.
        - _epistasis_analysis.ipynb_ - Contains code for running and results of _Wilcoxn Signed-Rank Tests_ - testing for statistically significant differences in the distributions of squared error between the WDS-epi model to WT model, WDS-epi model to the expected additive mutation λmax values (EAV), and WT model to EAV, respectively.
      - **_full_iter_sample_results_** - Contains results of our 'sample iterate test' on each dataset - where x number of datapoints are removed from the training data prior to training and used as test data. This is repeated until all points have been sampled once.
      - **_imputation_tests_** - Contains results comparing the predictions of phylogentic imputation and our ML models for the same data points (varies by dataset).
      - **_main_model_results_** - Contains model training results for each dataset, seperated by alignment method used (MAFFT, MUSCLE, and GBLOCKS following MUSCLE alignment) and then by database version (i.e. VPOD_1.0, VPOD_1.1, and VPOD_1.3).
      - **_msp_tests_** - Contains results for model predictions from each dataset on thirty unseen wild-type invertebrate opsins from a separately curated MSP dataset.
      - **_perf_v_tds_** - Contains results tracking the correlation between dataset size and model R^2, to better understand how training data relate to model performance.
      - **_phylo_weighted_cv_** - Contains results and all subsequent analysis of our systematic **_Phylogenetic-Weighted Cross-Validation_** method (including comparative and methodology analysis) [VPOD_1.3 only].
      - **_sws_ops_prediction_comparison_test_** - Contains results comparing predictive capabilities of models trained on different data subsets by randomly selecting and removing the same 25 wild-type ultraviolet or short-wave sensitive opsins from the training data of the WDS, Vertebrate, WT, and UVS/SWS.
      - **_wt_model_predicting_all_mutants_test_** - Contains results of the Wild-Type model (which lacks data from artificially mutated sequences) predictions on all experimentally mutated opsins in VPOD.
        - This folder contains results for several versions of VPOD (i.e. VPOD_1.0 and VPOD_1.1).
        - _mutant_pred_analysis.ipynb_ - Contains code to run and results for _Wilcoxn Signed-Rank Test_ - testing for statistically significant difference in the distributions of squared prediction errors for the whole-dataset model and wild-type model on all mutant data . Results are displayed directly in the notebook.

  ---

  ### Scripts & Notebooks
  Instructions for accessing scripts and Jupyter notebooks used to create the database and train ML models.

   ### Sequence Manipulation
   NOTE - It's recommended that you open these scripts in a compiler for a more detailed explination of how to properly use them* 
   * Navigate to the folder _scripts_n_notebooks/sequence_manipulation_
       - ```chimeras.py``` - Used to generate chimeric sequences from sequences accessible via NCBI by accession or manually enter the raw sequence.
           - Opsin _chimeras_ are mutants where one or more transmembrane domains are copied from a different opsin to replace the original. 
           - Contact us for assistance with using this script.
       - ```in_silico_dms.py``` - Used to generate mutants via an _in-silico_ 'site-saturated mutagenesis' // 'deep-mutational-scanning' (DMS).
       - ```mutagenesis.py``` - Used to generate site and amino-acid-specific mutants from a sequence accessible via NCBI by accession number or manually enter the raw sequence.
       - ```reciprocal_mutagenesis.py``` - Used to automatically generate site and amino-acid-specific mutants based on the sites that differ between two sequences the user wishes to compare (generates all single and multi-mutant combinations).
       - ```mutagenesis_ex_nb.ipynb``` - A notebook meant to demonstrate the functionalities for all the scripts detailed above ^^^ (We recommend trying this out first if you are familiar with Jupyter-Notebookss).
    
   ### Machine Learning Tests and Workflows
   * Navigate to the folder **_scripts_n_notebooks_**
      - Select the folder **_vpod_ML_workflows_** to access notebooks used for training ML models.
          - ```vpod_main_wf.ipynb``` - Primary notebook for users, with a full pipeline for everything from creating a local instance of VPOD using SQLite to formatting datasets and training ML models for λmax predictions.
          - **_boot_strap_model_gen_** folder - Contains scripts for generating ensembles of 100 bootstraped models for all VPOD dataset using one-hot or amino-acid property encoding methods [VPOD_1.3 only]
            - ```vpod_bootstrap_gen_aa_prop_wf.py``` - Generates 100 bootstraped models for each VPOD dataset using amino-acid property encoding.
            - ```vpod_bootstrap_gen_one_hot_wf.py``` - Generates 100 bootstraped models for each VPOD dataset using one-hot encoding.
          - **_mine_n_match_** folder - Contains scripts and data necessary for utilizing the _Mine-n-Match_ (MNM) workflow (as detailed in _Frazer & Oakley 2025_).
            - ```mine_n_match_workflow.ipynb``` - Primary notebook for running the MNM workflow. Note, to run this full workflow you should READ the notebook first and THEN follow directions to download [OPTICS](https://github.com/VisualPhysiologyDB/optics) from its corresponding _GitHub_ repository to make it available as a tool in the notebook (a totally seperate/dedicated Conda environment is recommended for using downloading OPTICS/using MNM). 
            - **_data_sources_** folder - Contains raw and formatted data sources used to compile sequence and _in_vivo_ λmax compendiums neccessary for running the main MNM workflow. For more information on these sources see _Frazer & Oakley 2025_ or read the markdown text in the main MNM workflow.
            - **_mnm_data_** folder - Contains the all the outputs from a single run of the MNM workflow.
              - Folder generated from running MNM are named and automatically timestamped by the user of the notebook in the format of 'mnm_on_XXX_date_time'.
                - Results pertaining the creation of VPOD_1.3 are in this folder, under the sub-folder - 'mnm_on_all_dbs_2025-02-24_16-29-54'
              - The final output of the MNM workflow follow the format of 'mnm_on_xxx_results_fully_filtered.csv'
              - For more information on the other outputs, kindly refer to the markdown text within the notebook first, then see _Frazer & Oakley 2025_ or contact us directly.
            - **_mnm_scripts_** folder - Contains the script necessary for running the main functions called in the main MNM workflow notebook (```mine_n_match_functions.py```).  This script is generally well commented for the sake of reproducability/reusability (although this script is more of a module holder, and not meant to be run via command-line). 
          - _substests_ folder - Contains notebooks used for subtests outlined in _Frazer et al. 2024_.
            - **_grid_search_opt_** folder - Contains scripts for running grid-search optimization ML models trained on all VPOD datasets (VPOD_1.3 only).
              - ```vpod_aa_prop_grid_search_iter.py``` - Iteratively trains models on all VPOD datasets utilizing amino-acid property encoding.
              - ```vpod_one_hot_grid_search_iter.py``` - Iteratively trains models on all VPOD datasets utilizing one-hot encoding.
            - **_phylo_weighted_cv_** folder - Contains scripts and trees needed to run _Phylogentically-Weighted Cross-Validation_ (PW-CV; see _Frazer & Oakley 2025_ for more detail)
              - ```run_phylo_weighted_cv.py``` - Main script for orchestrating the systematic creation phylogenetically weighted cross-validation folds and subsequent model training across several distance threshold.
              - ```phylo_weighted_cv_method_walkthrough.ipynb``` - Walkthrough of all the main steps from PW-CV for those who are interested.
              - **_trees_** folder - Contains all the tree files neccessary for calculating the phylogenetic distance matrix used in the main PW-CV script.
            - ```vpod_aa_prop_combos_wf.py``` - Iteratively trains ML models for all VPOD datasplits using all possible combinations of the available amino-acid properties to encode sequences by (output is a large .csv file with model metrics for each combination of amino acid properties).
              - **Warning** - this script can run for along time! We suggest running this script on a cluster where you can allot at least 50-100 hours to just this script.  
            - ```vpod_all_models_iter.py``` - Iteratively trains a single, unpreturbed, ML model for each of the VPOD datasets (VPOD_1.3 only); options for both one-hot and amino-acid property encoding. 
            - ```vpod_wf_iterate_train_all_subsets.ipynb``` - Iteratively train each dataset with no modifications to the training process; simply streamlines the training of all datasets in case of desire. 
            - ```vpod_wf_wt_mut_test.ipynb``` - Trains and tests the predictive capabilities of the Wild-Type model (which lacks data from artificially mutated sequences) on all experimentally mutated opsins in VPOD.
            - ```vpod_main_wf_msp_iterate.ipynb``` - Iteratively train and test models from each dataset on thirty unseen wild-type invertebrate opsins from a separately curated MSP dataset.
            - ```vpod_wf_imp_sample_test_iterate.ipynb``` - Iteratively train and test models from each dataset on randomly sampled subset of data; for comparison to predictions made by phylogentic imputation.
            - ```vpod_wf_iterate_all_sample_t1_ops.ipynb``` - Iteratively subsample T1 dataset, removing 'x' datapoints before training before training and using it as test-data until all datapoints are sampled once.
            - ```vpod_wf_iterate_epistasic_muts.ipynb``` - Iteratively subsample whole-dataset of mutations which demonstrate epistatic interactions between mutations, removing X datapoints before training before training and using it as test-data until all datapoints are sampled once (can also be modified to remove all epistatic mutants at once, as detailed in _Frazer et al. 2024_).
            - ```vpod_wf_iterate_model_perf_vs_tds.ipynb``` - Iteratively adds or substracts 'x' datapoints from dataset to track the correlation between dataset size and model performance.
            - ```vpod_wf_iterate_subsample.ipynb``` -  Iteratively subsample target dataset, removing 'x' datapoints before training before training and using it as test-data until all datapoints are sampled once.
            - ```vpod_wf_iterate_sws_comp.ipynb``` - Iteratively removes the same 25 randomly wild-type ultraviolet or short-wave sensitive opsins from the training data of the WDS, Vertebrate, WT, and UVS/SWS to compare predictive capabilities of models trained on different data subsets.
      - Select the folder **_figure_making_** to access the Jupyter notebook **_figuremaking.ipynb_** used to generate some of the figures used in _Frazer et al. 2024_.
          - _figures_ contains a collection of completed figures and drafts use in _Frazer et al. 2024_ - seperated by version of the database used to generate the figures (i.e VPOD_1.0 or VPOD_1.1)
            * Select the **_opsin_wt_tree_** folder for all files used to make the wild-type opsin gene-tree, supplementary figure 10 (S10), from _Frazer et al. 2024_.
            * We've also provided it as a _.svg_ file in this same folder or **[click here](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db/blob/main/scripts_n_notebooks/figure_making/figures/VPOD_1.1/opsin_wt_tree/wt_tree_w_ancestors.svg)** to download.
      - Select the folder **_phylogenetic_imputation_** to access all files used to predict λmax via phylogenetic imputation and compare with predictions made by ML, as detailed in  _Frazer et al. 2024_.
          - ```Phylogenetic_Imputation.Rmd``` - Used to load tree files and make λmax predictions via phylogenetic imputation **[Requires RStudio]**.
          - **_trees_** folder - Contains all the tree files and λmax meta-data neccessary for predictions via phylogentic imputation.
        
---

### Training New ML Models With _deepBreaks_
Instructions for using VPOD and training ML models with _deepBreaks_.

1. Follow the directions and install _[deepBreaks](https://github.com/omicsEye/deepbreaks)_ using the guide provided _deepBreaks_ on the GitHub.
2. Refer to the _requirements.txt_ provided in the base of this repostiory and ensure all necessary packages are installed in a dedicated environment (using Conda is recommended).
3. Navigate to _scripts_n_notebooks/vpod_ml_workflows_ and open _vpod_main_wf.ipynb_.
   - To train models using _raw_database_files_, start from the top of the document and follow the instructions provided in the notebook.
   - To train models using _formatted_data_subsets_, scroll down to **Step 3: deepBreaks** and follow the instructions provided in the notebook.
     
---

## License
All data and code is covered under a GNU General Public License (GPL)(Version 3), in accordance with Open Source Initiative (OSI)-policies

## Citation
As of 2025 we have two publications connected to VPOD. Choose which one to cite based on what information you use, or cite both!

IF citing this GitHub and its contents use the following DOI provided by Zenodo... [VPOD_1.2]

    10.5281/zenodo.10667840
    
IF elemments of this database connected to _"Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD)"_ use the following citation...

    Seth A. Frazer, Mahdi Baghbanzadeh, Ali Rahnavard, Keith A. Crandall, & Todd H Oakley. Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD). GigaScience, 2024.09.01. https://doi.org/10.1093/gigascience/giae073

IF citing element of this database (or OPTICS) connected to _"Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map"_ use the following citation...
  
    Seth A. Frazer, Todd H. Oakley. Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map. bioRxiv, 2025.08.22.671864. https://doi.org/10.1101/2025.08.22.671864 [pre-print]
    
## Contact
Contact information for author questions or feedback.

**Todd H. Oakley**  
  * oakley@ucsb.edu
  * https://orcid.org/0000-0002-4478-915X

**Seth A. Frazer** 
  * sethfrazer@ucsb.edu
  * https://orcid.org/0000-0002-3800-212X
---

## Additional_Resources 
  **[Here](https://tinyurl.com/u7hn9adm)** is a link to a bibliography of the publications used to build VPOD_1.1 (Full version not yet released)
  
  If you know of publications for training opsin ML models not included in the VPOD_1.2 database, please send them to us through **[this form](https://tinyurl.com/29afaxyr)**
