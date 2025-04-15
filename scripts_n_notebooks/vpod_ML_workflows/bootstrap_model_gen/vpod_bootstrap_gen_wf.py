# %% [markdown]
# # <font color=green>deepBreaks Applications</font>
# ## Modeling spectral tuning sites of opsin proteins based on amino-acid sequence...  

# %%
# importing deepBreaks libraries 
from deepBreaks.utils_alt2 import get_models, get_scores, get_empty_params, make_pipeline
from deepBreaks.preprocessing import MisCare, ConstantCare, URareCare, CustomOneHotEncoder
from deepBreaks.preprocessing import FeatureSelection, CollinearCare
from deepBreaks.preprocessing import read_data
from deepBreaks.models import model_compare_cv, finalize_top, importance_from_pipe, mean_importance, summarize_results
from deepBreaks.visualization import plot_scatter, dp_plot, plot_imp_model, plot_imp_all
from deepBreaks.preprocessing import write_fasta
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

encoding = 'hot'

n_iterations = 100
rng = np.random.default_rng()  # Initialize NumPy's random number generator


# %%
for meta, seq, ds in zip(meta_data_list, seq_data_list, ds_list):
    # path to sequences of interest
    seqFileName = f'{path}/{seq}' 
    # path to corresponding metadata of interest
    metaDataFileName = f'{path}/{meta}' 
    # making a unique directory for saving the reports of the analysis
    #print('direcory preparation')
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    report_dir = str(f'{ds}_bootstrap_100_{dt_label}')
    os.makedirs(report_dir)

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

    full_tr = tr.copy()

    #setting the paramaters for our ML pipeline
    prep_pipeline = make_pipeline(
        steps=[
            ('mc', MisCare(missing_threshold=0.05)),
            ('cc', ConstantCare()),
            ('ur', URareCare(threshold=0.025)),
            ('cc2', ConstantCare()),
            ('one_hot', CustomOneHotEncoder()),
            ('feature_selection', FeatureSelection(model_type=ana_type, alpha=0.10, keep=False)),
            ('collinear_care', CollinearCare(dist_method='correlation', threshold=0.01, keep=False))
        ])

    for i in range(n_iterations):

        tr = full_tr
        random_state = rng.integers(0, 2**32 - 1)  # Generate a new random seed for perturbation
        X_res, y_res = resample(tr, y, random_state=random_state)
        
        #training models
        report, top = model_compare_cv(X=X_res, y=y_res, preprocess_pipe=prep_pipeline,
                                    models_dict=get_models(ana_type=ana_type, dataset=ds, encoding=encoding),
                                    scoring=get_scores(ana_type=ana_type),
                                    report_dir=report_dir,
                                    cv=10, ana_type=ana_type, cache_dir=report_dir)

                                    
        #setting parameters for tuning the top performing models
        prep_pipeline = make_pipeline(
            steps=[
                ('mc', MisCare(missing_threshold=0.05)),
                ('cc', ConstantCare()),
                ('ur', URareCare(threshold=0.025)),
                ('cc2', ConstantCare()),
                ('one_hot', CustomOneHotEncoder()),
                ('feature_selection', FeatureSelection(model_type=ana_type, alpha=0.10, keep=True)),
                ('collinear_care', CollinearCare(dist_method='correlation', threshold=0.01, keep=True))
            ])

        modified_top = []
        mtml = []
        for model in top:
            modified_top.append(make_pipeline(steps=[('prep', prep_pipeline), model.steps[-1]]))
            my_top_models = str(model[1:])
            #print(my_top_models)
            my_top_models = my_top_models.split("'")[3]
            mtml.append(my_top_models)

        #tuning the top 3 performing models 
        top = finalize_top(X=X_res, y=y_res, top_models=modified_top, grid_param=get_empty_params(),report_dir=report_dir, cv=10)
        #summarize the results by extracting feature importance and p-values and grouping correlated features.
        sr = summarize_results(top_models=top, report_dir=report_dir)
        mean_imp = mean_importance(top, report_dir=report_dir)
        
        try:
            original_file = f'{report_dir}/importance_report.csv'  
            new_file = f'{report_dir}/importance_report_iter_{str(i)}.csv' 
            shutil.copy(original_file, new_file)
            os.remove(original_file)
        except:
            raise Exception('Cannot copy or delete importance_report file either because it does not exist or the direcotry is incorrect')
        
        try:
            original_file = f'{report_dir}/model_report.csv'  
            new_file = f'{report_dir}/model_report_iter_{str(i)}.csv' 
            shutil.copy(original_file, new_file)
            os.remove(original_file)
        except:
            raise Exception('Cannot copy or delete model_report file either because it does not exist or the direcotry is incorrect')

        try:
            for model in mtml:
                original_file = f'{report_dir}/{model}.pkl'  
                new_file = f'{report_dir}/{model}_{str(i)}.pkl' 
                shutil.copy(original_file, new_file)
                os.remove(original_file)
        except:
            raise Exception('Cannot copy or delete model pkl file either because it does not exist or the direcotry is incorrect')

        try:
            shutil.rmtree(f'{report_dir}/joblib')
        except:
            pass
        try:
            os.remove(f'{report_dir}/joblib')
        except:
            pass


