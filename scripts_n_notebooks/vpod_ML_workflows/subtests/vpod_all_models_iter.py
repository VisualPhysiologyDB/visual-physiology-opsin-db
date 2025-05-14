# %% [markdown]
# # <font color=green>deepBreaks Applications</font>
# ## Modeling the phenotypes and spectral tuning sites of opsin proteins based on amino-acid sequence...  

# %% [markdown]
# # <font color=red>STEP 3: deepBreaks</font>
# ## THIS IS A LONG SECTION! 
# ### **Output** = folder containing all results from model training, including comparison of model performances, an amino-acid site importance report + figures, and the top 5 trained models in a .pkl file format.

# %%
# importing deepBreaks libraries 
from deepBreaks.utils_alt2 import get_models, get_scores, make_pipeline, get_best_aa_prop_combos
from deepBreaks.preprocessing import MisCare, ConstantCare
from deepBreaks.preprocessing import FeatureSelection, CollinearCare, AminoAcidPropertyEncoder, CustomOneHotEncoder, URareCare
from deepBreaks.preprocessing import read_data
from deepBreaks.models import model_compare_cv
import warnings
import datetime
import os
import shutil 
import numpy as np
import csv
import pandas as pd
import ast
import re
# %%
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# %%
# defining user params, file pathes, analysis type
print('Beginning script!\n')
#assign your path to folder containing all the datasplits
path = './vpod_1.2_data_splits_2025-02-28_15-51-04'
save_to = './vpod_1.3_all_models'
data_dict = {'wds_mnm':{'meta':'wds_mnm_meta.csv','seq':'wds_mnm_aligned_VPOD_1.2_het.fasta'},
             'wt_mnm':{'meta':'wt_mnm_meta.csv','seq':'wt_mnm_aligned_VPOD_1.2_het.fasta'},
             'wt_vert_mnm':{'meta':'wt_vert_mnm_meta.csv','seq':'wt_vert_mnm_aligned_VPOD_1.2_het.fasta'},
             'inv_mnm':{'meta':'inv_mnm_meta.csv','seq':'inv_mnm_aligned_VPOD_1.2_het.fasta'},
             'vert_mnm':{'meta':'vert_mnm_meta.csv','seq':'vert_mnm_aligned_VPOD_1.2_het.fasta'},
             'wds':{'meta':'wds_meta.tsv','seq':'wds_aligned_VPOD_1.2_het.fasta'},
             'wt':{'meta':'wt_meta.tsv','seq':'wt_aligned_VPOD_1.2_het.fasta'},
             'wt_vert':{'meta':'wt_vert_meta.tsv','seq':'wt_vert_aligned_VPOD_1.2_het.fasta'},
             'inv':{'meta':'inv_meta.tsv','seq':'inv_only_aligned_VPOD_1.2_het.fasta'},
             'vert':{'meta':'vert_meta.tsv','seq':'vert_aligned_VPOD_1.2_het.fasta'},
             't1':{'meta':'Karyasuyama_T1_ops_meta.tsv','seq':'Karyasuyama_T1_ops.fasta'}
            }


all_models_results_files = f'{save_to}/gs_model_report_all_datasets.csv'
try:
    all_model_results_df = pd.read_csv(filepath_or_buffer=all_models_results_files, index_col=0)
except:
    all_model_results_df = pd.DataFrame(data=data_dict.values(), index=data_dict.keys())
    
# name of the phenotype
mt = 'Lambda_Max'

# type of the sequences
seq_type = 'aa'

# type of the analysis if it is a classification model, then we put cl instead of reg
ana_type = 'reg' 

gap_threshold = 0.50
    
#Whether or not you want to drop the reference sequence from the training data- Usually 'Bovine' or 'Squid'
drop_ref = False

encoding_methods_list = ['one_hot','aa_prop']

if 'R2' not in all_model_results_df.columns:
    all_model_results_df['R2'] = None

for encoding in encoding_methods_list:
    for ds in all_model_results_df.index.to_list():
        try:
            #print(all_model_results_df.loc[ds])    
            if all_model_results_df.loc[ds]['R2'] is None or pd.isna(all_model_results_df.loc[ds]['R2']):
                print(f"Grid-Search Params Empty for Dataset {ds}\n")
                # path to sequences of interest
                seq = all_model_results_df.loc[ds]['seq']
                seqFileName = f'{path}/{seq}'
                # path to corresponding metadata of interest
                meta = all_model_results_df.loc[ds]['meta']
                metaDataFileName = f'{path}/{meta}' 


                # %%
                # making a unique directory for saving the reports of the analysis
                dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

                if encoding == 'aa_prop':
                    props_to_keep = get_best_aa_prop_combos(ds)
                    props_used = ''
                    for props in props_to_keep:
                        props_used += props + '_'
                    report_dir = str(save_to + '/' + ds + '_' + props_used + mt + '_gs_optimized_' + dt_label)
                else:
                    props_to_keep = None
                    report_dir = str(save_to + '/' + ds + '_gs_optimized_' + mt + '_' + dt_label)

                os.makedirs(report_dir)
                # %%
                #print('reading meta-data')
                # importing metadata
                meta_data = read_data(metaDataFileName, seq_type = None, is_main=False)
                # importing sequences data
                #print('reading fasta file')

                tr = read_data(seqFileName, seq_type = seq_type, is_main=True, gap_threshold=gap_threshold)

                #shutil.copy2(f'{seqFileName}',report_dir)
                #write_fasta(dat = tr, fasta_file = f'{seqFile}_gap_dropped.fasta' , report_dir = report_dir)

                # %%
                try:
                    reference_seq = tr.loc['Bovine'].copy()
                    ref_seq_name = 'bovine'
                    if drop_ref == True:
                        meta_data = meta_data.drop('Bovine')
                    #print(bovine)
                except:
                    try:
                        reference_seq = tr.loc['Squid'].copy()
                        ref_seq_name = 'squid'
                        #print(squid)
                    except:
                        pass
                
                try:
                    reference_seq.to_csv(path_or_buf= f'{report_dir}/ref_sequence.csv',index = True,mode="w")
                except:
                    pass

                # %%
                tr = tr.merge(meta_data.loc[:, mt],  left_index=True, right_index=True)
                tr.shape

                # %%
                y = tr.loc[:, mt].values
                tr.drop(mt, axis=1, inplace=True)
                #print('Shape of data is: ', tr.shape)
                
                
                # %%
                # Create prep_pipeline outside the loop to avoid re-creation
                if encoding == 'aa_prop':
                    prep_pipeline = make_pipeline(
                        steps=[
                            ('mc', MisCare(missing_threshold=0.05)),
                            ('cc', ConstantCare()),
                            ('aa_prop', AminoAcidPropertyEncoder(props_to_keep = props_to_keep)),
                            ('feature_selection', FeatureSelection(model_type=ana_type, alpha=0.10, keep=False)),
                            ('collinear_care', CollinearCare(dist_method='correlation', threshold=0.01, keep=False))
                        ])
                elif encoding == 'one_hot':
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
                else: 
                    raise Exception('You have provided an unsupported option for encoding. Please chhose either "one_hot" or "aa_prop"')


                report, top = model_compare_cv(X=tr, y=y, preprocess_pipe=prep_pipeline,
                                            models_dict=get_models(ana_type=ana_type, dataset=ds, encoding=encoding),
                                            scoring=get_scores(ana_type=ana_type),
                                            report_dir=report_dir,
                                            cv=10, ana_type=ana_type, cache_dir=report_dir)
                    
                # Record R^2 values for each model.
                for model_name, r2_value, mae_value, mape_value, mse_value, rmse_value in zip(report.index.to_list(), report["R2"],report["MAE"],report["MAPE"],report["MSE"],report["RMSE"]):
                    all_model_results_df.loc[ds]['encoding'] = encoding
                    all_model_results_df.loc[ds]['model'] = model_name
                    all_model_results_df.loc[ds]['R2'] = r2_value
                    all_model_results_df.loc[ds]['MAE'] = mae_value
                    all_model_results_df.loc[ds]['MAPE'] = mape_value
                    all_model_results_df.loc[ds]['MSE'] = mse_value
                    all_model_results_df.loc[ds]['RMSE'] = rmse_value
        except:
            try:
                all_model_results_df.to_csv(path_or_buf=all_models_results_files,index=True)

            except Exception as e:
                raise Exception(e)

all_model_results_df.to_csv(path_or_buf=all_models_results_files,index=True)

                    