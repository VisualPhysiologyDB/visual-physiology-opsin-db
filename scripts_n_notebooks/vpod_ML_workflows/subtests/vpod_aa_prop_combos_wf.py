# %% [markdown]
# # <font color=green>deepBreaks Applications</font>
# ## Modeling spectral tuning sites of opsin proteins based on amino-acid sequence...  

# %%
# importing deepBreaks libraries 
from deepBreaks.utils_alt2 import get_models, get_scores, get_empty_params, make_pipeline
from deepBreaks.preprocessing import MisCare, ConstantCare, AminoAcidPropertyEncoder
from deepBreaks.preprocessing import FeatureSelection, CollinearCare
from deepBreaks.preprocessing import read_data
from deepBreaks.models import model_compare_cv, finalize_top, importance_from_pipe, mean_importance, summarize_results
from sklearn.utils import resample
from random import random
import numpy as np
import csv
import pandas as pd
import warnings
import datetime
import os
import shutil
import time
import itertools

warnings.filterwarnings("ignore")
warnings.simplefilter('ignore')


# %%
# defining user params, file pathes, analysis type

#assign your path to folder containing all the datasplits
path = './vpod_1.2_data_splits_2024-08-20_16-14-09'
meta_data_list = ['wds_meta.tsv','wt_meta.tsv','wt_vert_meta.tsv', 'inv_meta.tsv', 'vert_meta.tsv']
seq_data_list = ['wds_aligned_VPOD_1.2_het.fasta','wt_aligned_VPOD_1.2_het.fasta','wt_vert_aligned_VPOD_1.2_het.fasta', 'inv_only_aligned_VPOD_1.2_het.fasta', 'vert_aligned_VPOD_1.2_het.fasta']
ds_list = ['wds', 'wt', 'wt_vert', 'inv', 'vert']

# name of the phenotype
mt = 'Lambda_Max'

# type of the sequences
seq_type = 'aa'

# type of the analysis if it is a classification model, then we put cl instead of reg
ana_type = 'reg' 

gap_threshold = 0.5

#Specify which properties you want to keep for the amino-acid property encoding:
#We keep FIVE by deafult - 'H1, H3, P1, NCI, MASS' 
#But NINE total are avaliable -'H1, H2, H3, P1, P2, V, NCI, MASS, and SASA' 
#If you want to keep ALL aa props, just set props_to_keep = 'all'
# Or specify the properties in list format props_to_keep = ['H1', 'H3', 'P1', 'NCI', 'MASS']
encoding = 'aa_prop'

props_to_test = ['H1', 'H2', 'H3', 'P1', 'P2', 'V', 'NCI', 'MASS', 'SASA', 'SCT', 'PKA', 'PKB']
prop_combos = []
x = 1
# Get a list of all possible combinations of amino acid properties to train on
while x in range(len(props_to_test)+1):
    combinations = list(itertools.combinations(props_to_test, x))
    for combo in combinations:
        temp_list = []
        for item in combo:
            temp_list.append(item)
        prop_combos.append(temp_list)
    x+=1

timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

results_folder = f"aa_prop_combos_analysis_{timestamp}"
os.makedirs(results_folder, exist_ok=True)
# Create results file (choose CSV or Excel based on your preference)


# %%
for meta, seq, ds in zip(meta_data_list, seq_data_list, ds_list):
    # path to sequences of interest
    seqFileName = f'{path}/{seq}' 
    # path to corresponding metadata of interest
    metaDataFileName = f'{path}/{meta}' 
    # making a unique directory for saving the reports of the analysis
    #print('direcory preparation')
    
    #print('reading meta-data')
    # importing metadata
    meta_data = read_data(metaDataFileName, seq_type = None, is_main=False)
    metaFile = metaDataFileName.split('/')[1]
    # importing sequences data
    #print('reading fasta file')
    tr = read_data(seqFileName, seq_type = seq_type, is_main=True, gap_threshold=gap_threshold)
    #merging in lambda max values, simultaneously dropping all sequences without entries in metadata file
    tr = tr.merge(meta_data.loc[:, mt],  left_index=True, right_index=True)
    #tr.shape
    seqFile = seqFileName.split('/')[2]
    #print(seqFile)
    seqFile = seqFile.split('.')[0]+'.'+seqFile.split('.')[1]
    #write_fasta(dat = tr, fasta_file = f'{seqFile}_gap_dropped.fasta' , report_dir = report_dir)

    y = tr.loc[:, mt].values
    tr.drop(mt, axis=1, inplace=True)
    
    results_file = os.path.join(results_folder, f"{ds}_aa_prop_combos_results.csv")
    with open(results_file, "w") as f:
        f.write("Dataset,Props_Used,Model,R2,MAE,MAPE,MSE,RMSE\n") # Header

    for props_to_keep in prop_combos:
        props_used = ''
        for props in props_to_keep:
            props_used += props + '_'

        report_dir = str(f'./{results_folder}/{ds}_{props_used}')
        os.makedirs(report_dir)

        #settingthe paramaters for our ML pipeline
        prep_pipeline = make_pipeline(
            steps=[
                ('mc', MisCare(missing_threshold=0.05)),
                ('cc', ConstantCare()),
                ('aa_prop', AminoAcidPropertyEncoder(props_to_keep = props_to_keep)),
                ('feature_selection', FeatureSelection(model_type=ana_type, alpha=0.10, keep=False)),
                ('collinear_care', CollinearCare(dist_method='correlation', threshold=0.01, keep=False))
            ])

        #training models
        report, top = model_compare_cv(X=tr, y=y, preprocess_pipe=prep_pipeline,
                                    models_dict=get_models(ana_type=ana_type, dataset=ds, encoding=encoding),
                                    scoring=get_scores(ana_type=ana_type),
                                    report_dir=report_dir,
                                    cv=10, ana_type=ana_type, cache_dir=report_dir)
        
        # Record R^2 values for each model
        with open(results_file, "a") as f:
            for model_name, r2_value, mae_value, mape_value, mse_value, rmse_value in zip(report.index.to_list(), report["R2"],report["MAE"],report["MAPE"],report["MSE"],report["RMSE"]):
                f.write(f"{ds},{props_used},{model_name},{r2_value},{mae_value},{mape_value},{mse_value},{rmse_value}\n")
        
        
        try:
            shutil.rmtree(f'{report_dir}/joblib')
        except:
            pass
        try:
            os.remove(f'{report_dir}/joblib')
        except:
            pass
        try:
            shutil.rmtree(report_dir)
        except:
            pass
        try:
            os.remove(report_dir)
        except:
            pass

