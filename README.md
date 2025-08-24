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
* **UPDATE!!!** - Welcome to VPOD_1.3! The latest release of the Visual Physiology Opsin Database. We have quite a bit of new data, results and tools to check-out.
* **Want an user-friendly way to predict opsin λmax with our models?**
  * Visit our GitHub for _[OPTICS](https://github.com/VisualPhysiologyDB/optics)_, an open-source tool, avaiable online or via command line, that predicts opsin phenotypes (λmax)

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
      - The folder _vpod_1.1_data_subsets_2024-05-02_ contains all data subsets for VPOD_1.1.
      - Files marked _xxx.txt_ (ex. _vert.txt_) are the unaligned data subsets.
      - Files marked _xxx_aligned.txt_ (ex. _vert.txt_) are the aligned data subsets (not-formatted; will become obsolete in later versions of VPOD).
      - Files marked _VPOD_xxx_het_1.1_ (ex. _VPOD_vert_het_1.1.fasta_) are the fully aligned and formatted data subsets.
      - Files marked _xxx_meta_ (ex. _wds_meta.tsv_) are the corresponding metadata files for each subset (includes species, accession, gene names, λmax, etc).
    * Select the **_raw_database_files_** folder for the raw components of the database that you can load into a SQLite database and create your own formatted dataset using steps 0-2 outlined in the _vpod_main_wf.ipynb_ Jupyter notebook. For more information on the specific meaning of data columns, see _Frazer et al. 2024_.
      - _litsearch.csv_ - All literature search information relevant to the current version of VPOD.  
      - _references.csv_ - All publication references relevant to the current version of VPOD.
      - _opsins.csv_ - All opsin sequence data and taxonomic meta-data.
      - _heterologous.csv_ - All opsin phenotype data (λmax) and experiment related meta-data. 
  
  ### Results
  * Navigate to the folder _results_files_ - All subfolders contain results specific to different subtests outlined in _Frazer et al. 2024_ (for VPOD_1.0 and VPOD_1.1) or _Frazer and Oakley_ (for VPOD_1.3).
      - **_aa_prop_combos_** - 
      - **_aa_prop_vs_one_hot_** - 
      - **_epistasis_pred_test_** - Contains predictions by and comparisons between our WT and WDS models on 111 'epistatic mutations' (non-additive) to more generally the capabilities of our ML models to predict epistatic interactions between mutations.
        - _epistasis_analysis.ipynb_ - Contains code for and results of _Wilcoxn Signed-Rank Tests_ - testing for statistically significant differences in the distributions of squared error between the WDS-epi model to WT model, WDS-epi model to the expected additive mutation λmax values (EAV), and WT model to EAV, respectively.
      - **_full_iter_sample_results_** - Contains results of our 'sample iterate test' on each dataset - where x number of datapoints are removed from the training data prior to training and used as test data. This is repeated until all points have been sampled once.
      - **_imputation_tests_** - Contains results comparing the predictions of phylogentic imputation and our ML models for the same data points (varies by dataset).
      - **_main_model_results_** - Contains model training results for each dataset, seperated by alignment method used (MAFFT, MUSCLE, and GBLOCKS following MUSCLE alignment) and then by database version (i.e. VPOD_1.0, VPOD_1.1, and VPOD_1.3).
      - **_msp_tests_** - Contains results for model predictions from each dataset on thirty unseen wild-type invertebrate opsins from a separately curated MSP dataset.
      - **_perf_v_tds_** - Contains results tracking the correlation between dataset size and model R^2, to better understand how training data relate to model performance.
      - **_phylo_weighted_cv_** - 
      - **_sws_ops_prediction_comparison_test_** - Contains results comparing predictive capabilities of models trained on different data subsets by randomly selecting and removing the same 25 wild-type ultraviolet or short-wave sensitive opsins from the training data of the WDS, Vertebrate, WT, and UVS/SWS.
      - **_wt_model_predicting_all_mutants_test_** - Contains results of the Wild-Type model (which lacks data from artificially mutated sequences) predictions on all experimentally mutated opsins in VPOD.
        - This folder contains results for several versions of VPOD (i.e. VPOD_1.0 and VPOD_1.1).
        - _mutant_pred_analysis.ipynb_ - Contains code to run and results for _Wilcoxn Signed-Rank Test_ - testing for statistically significant difference in the distributions of squared prediction errors for the whole-dataset model and wild-type model on all mutant data . Results are displayed directly in the notebook.

    
  ### Scripts & Notebooks
  Instructions for accessing scripts and Jupyter notebooks used to create the database and train ML models.

   ### Sequence Manipulation
   NOTE - It's recommended that you open these scripts in a compiler for a more detailed explination of how to properly use them* 
   * Navigate to the folder _scripts_n_notebooks/sequence_manipulation_
       - **_chimeras.py_** - Used to generate chimeric sequences from sequences accessible via NCBI by accession or manually enter the raw sequence.
           - Opsin _chimeras_ are mutants where one or more transmembrane domains are copied from a different opsin to replace the original. 
           - Contact us for assistance with using this script.
       - **_in_silico_dms.py_** - 
       - **_mutagenesis.py_** - Used to generate mutants from a sequence accessible via NCBI by accession number or manually enter the raw sequence.
       - **_reciprocal_mutagenesis.py_**
       - **_mutagenesis_ex_nb.ipynb_** - 
   ### Machine Learning Tests and Workflows
   * Navigate to the folder _scripts_n_notebooks_
      - Select the folder _vpod_ML_workflows_ to access notebooks used for training ML models.
          - **_vpod_main_wf.ipynb_** - Primary notebook for users, with a full pipeline for everything from creating a local instance of VPOD using SQLite to formatting datasets and training ML models for λmax predictions.
          - _boot_strap_model_gen_ folder - Contains...
            - **_vpod_bootstrap_gen_aa_prop_wf.py_** - 
            - **_vpod_bootstrap_gen_one_hot_wf.py_** -
          - _mine_n_match_ folder - Contains...
            - **_mine_n_match_workflow.ipynb_** - Primary notebook for MNM...
            - **_data_sources_** folder - Contains...
            - **_mnm_data_** folder - Contains...
            - **_mnm_scripts_** folder - Contains...
          - _substests_ folder - Contains notebooks used for subtests outlined in _Frazer et al. 2024_.
            - **_grid_search_opt_** folder - Contains...
              - **_vpod_aa_prop_grid_search_iter.py_** - Iteratively...
              - **_vpod_one_hot_grid_search_iter.py_** - Iteratively...
            - **_phylo_weighted_cv_** folder - Contains...
              - **_run_phylo_weighted_cv.py_** - Main script for orchestrating... 
              - **_phylo_weighted_cv_method_walkthrough.ipynb_** - Walkthrough of all the main steps from phylo CV...
            - **_vpod_aa_prop_combos_wf.py_** - Iteratively...
            - **_vpod_all_models_iter_** - Iteratively...
            - **_vpod_wf_iterate_train_all_subsets.ipynb_** - Iteratively train each dataset with no modifications to the training process; simply streamlines the training of all datasets in case of desire. 
            - **_vpod_wf_wt_mut_test.ipynb_** - Trains and tests the predictive capabilities of the Wild-Type model (which lacks data from artificially mutated sequences) on all experimentally mutated opsins in VPOD.
            - **_vpod_main_wf_msp_iterate.ipynb_** - Iteratively train and test models from each dataset on thirty unseen wild-type invertebrate opsins from a separately curated MSP dataset.
            - **_vpod_wf_imp_sample_test_iterate.ipynb_** - Iteratively train and test models from each dataset on randomly sampled subset of data; for comparison to predictions made by phylogentic imputation.
            - **_vpod_wf_iterate_all_sample_t1_ops.ipynb_** - Iteratively subsample T1 dataset, removing 'x' datapoints before training before training and using it as test-data until all datapoints are sampled once.
            - **_vpod_wf_iterate_epistasic_muts.ipynb_** - Iteratively subsample whole-dataset of mutations which demonstrate epistatic interactions between mutations, removing X datapoints before training before training and using it as test-data until all datapoints are sampled once (can also be modified to remove all epistatic mutants at once, as detailed in _Frazer et al. 2024_).
            - **_vpod_wf_iterate_model_perf_vs_tds.ipynb_** - Iteratively adds or substracts 'x' datapoints from dataset to track the correlation between dataset size and model performance.
            - **_vpod_wf_iterate_subsample.ipynb_** -  Iteratively subsample target dataset, removing 'x' datapoints before training before training and using it as test-data until all datapoints are sampled once.
            - **_vpod_wf_iterate_sws_comp.ipynb_** - Iteratively removes the same 25 randomly wild-type ultraviolet or short-wave sensitive opsins from the training data of the WDS, Vertebrate, WT, and UVS/SWS to compare predictive capabilities of models trained on different data subsets.
      - Select the folder **_figure_making_** to access the Jupyter notebook **_figuremaking.ipynb_** used to generate some of the figures used in _Frazer et al. 2024_.
          - _figures_ contains a collection of completed figures and drafts use in _Frazer et al. 2024_ - seperated by version of the database used to generate the figures (i.e VPOD_1.0 or VPOD_1.1)
            * Select the **_opsin_wt_tree_** folder for all files used to make the wild-type opsin gene-tree, supplementary figure 10 (S10), from _Frazer et al. 2024_.
            * We've also provided it as a _.svg_ file in this same folder or **[click here](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db/blob/main/scripts_n_notebooks/figure_making/figures/VPOD_1.1/opsin_wt_tree/wt_tree_w_ancestors.svg)** to download.
      - Select the folder _phylogenetic_imputation_ to access all files used to predict λmax via phylogenetic imputation and compare with predictions made by ML, as detailed in  _Frazer et al. 2024_.
          - **_Phylogenetic_Imputation.Rmd_** - Used to load tree files and make λmax predictions via phylogenetic imputation [Requires RStudio].
          - _trees_ - Contains all the tree files and λmax meta-data neccessary for predictions via phylogentic imputation.
---

### Training New ML Models With _deepBreaks_
Instructions for using VPOD and training ML models with _deepBreaks_.

1. Follow the directions and install _[deepBreaks](https://github.com/omicsEye/deepbreaks)_ using the guide provided _deepBreaks_ on the GitHub.
2. Refer to the _requirements.txt_ provided above and ensure all necessary packages are installed in a dedicated environment (using Conda is recommended).
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
  
    Seth A. Frazer & Todd H Oakley. Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map. bioRxiv, 2025.xx.xx. https://doi.org/xxx [Pre-Print]
    
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
