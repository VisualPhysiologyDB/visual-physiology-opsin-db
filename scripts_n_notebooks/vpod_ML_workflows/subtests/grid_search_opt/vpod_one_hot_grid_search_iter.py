# %% [markdown]
# # <font color=green>deepBreaks Applications</font>
# ## Modeling the phenotypes and spectral tuning sites of opsin proteins based on amino-acid sequence...  

# %% [markdown]
# # <font color=red>STEP 3: deepBreaks</font>
# ## THIS IS A LONG SECTION! 
# ### **Output** = folder containing all results from model training, including comparison of model performances, an amino-acid site importance report + figures, and the top 5 trained models in a .pkl file format.

# %%
# importing deepBreaks libraries 
from deepBreaks.utils import get_models, get_scores, get_exp_params, make_pipeline
from deepBreaks.preprocessing import MisCare, ConstantCare, URareCare, CustomOneHotEncoder
from deepBreaks.preprocessing import FeatureSelection, CollinearCare
from deepBreaks.preprocessing import read_data
from deepBreaks.models import model_compare_cv, finalize_top,  mean_importance, summarize_results
from deepBreaks.visualization import plot_scatter, dp_plot
import warnings
import datetime
import os
import shutil 
import numpy as np
import pandas as pd
import ast
import re
# %%
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

def string_to_numpy_array(string_list):
    """
    Converts a string representation of a list of np.float64 values to a numpy array.

    Args:
        string_list (str): A string in the format '[np.float64(0.9673339346349513), np.float64(0.9621805105052479), np.float64(0.9688362369860075)]'.

    Returns:
        numpy.ndarray: A numpy array of float64 values.
    """
    try:
        # Remove 'np.float64(' and ')' using regular expressions
        cleaned_string = re.sub(r'np\.float64\(|\)', '', string_list)

        # Parse the cleaned string using ast.literal_eval
        parsed_list = ast.literal_eval(cleaned_string)

        # Convert the parsed list to a numpy array of float64
        numpy_array = np.array(parsed_list, dtype=np.float64)
        return numpy_array
    except (ValueError, SyntaxError) as e:
        print(f"Error parsing string: {string_list}. Error: {e}")
        return None

def string_to_list_of_dicts(string_dicts):
    """
    Converts a string representation of a list of dictionaries to a list of dictionaries.

    Args:
        string_dicts (str): A string in the format "[{'key1': value1, 'key2': value2}, {...}]".

    Returns:
        list: A list of dictionaries.
    """
    try:
        # Use ast.literal_eval to safely parse the string as a Python literal
        list_of_dicts = ast.literal_eval(string_dicts)
        return list_of_dicts
    except (ValueError, SyntaxError) as e:
        print(f"Error parsing string: {string_dicts}. Error: {e}")
        return None  # or handle the error in another appropriate way
    
# %%
# defining user params, file pathes, analysis type
#assign path to folder containing all the datasplits
path = './vpod_1.2_data_splits_2025-02-28_15-51-04'
save_to = 'mnm_one_hot_gs'
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

one_hot_gs_file = f'./{save_to}/vpod_1.2_one_hot_gs_all_datasets.csv'
try:
    one_hot_gs_df = pd.read_csv(filepath_or_buffer=one_hot_gs_file,index_col=0)
except:
    one_hot_gs_df = pd.DataFrame(data=data_dict.values(), index=data_dict.keys())
# name of the phenotype
mt = 'Lambda_Max'

# type of the sequences
seq_type = 'aa'

# type of the analysis if it is a classification model, then we put cl instead of reg
ana_type = 'reg' 

gap_threshold = 0.50
    
#Whether or not you want to drop the reference sequence from the training data- Usually 'Bovine' or 'Squid'
drop_ref = False


if 'gs_params' not in one_hot_gs_df.columns:
    one_hot_gs_df['gs_params'] = None

for ds in one_hot_gs_df.index.to_list():
    #print(one_hot_gs_df.loc[ds])    
    if one_hot_gs_df.loc[ds]['gs_params'] is None or pd.isna(one_hot_gs_df.loc[ds]['gs_params']):
        try:
            print(f"Grid-Search Params Empty for Dataset {ds}\n")
            # path to sequences of interest
            seq = one_hot_gs_df.loc[ds]['seq']
            seqFileName = f'{path}/{seq}'
            # path to corresponding metadata of interest
            meta = one_hot_gs_df.loc[ds]['meta']
            metaDataFileName = f'{path}/{meta}' 


            # %%
            # making a unique directory for saving the reports of the analysis
            #print('direcory preparation')
            dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
            seqFile = seqFileName.split('/')[2]
            #print(seqFile)
            seqFile = seqFile.split('.')[0]+'.'+seqFile.split('.')[1]
            #print(seqFile)
            report_dir = str(save_to + '/' + ds + '_' + mt + '_' + dt_label)
            os.makedirs(report_dir, exist_ok=True)

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
            #y_ev = 1239.8 / np.array(y)

            # %%
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

            # %%
            report, top = model_compare_cv(X=tr, y=y, preprocess_pipe=prep_pipeline,
                                        models_dict=get_models(ana_type=ana_type),
                                        scoring=get_scores(ana_type=ana_type),
                                        report_dir=report_dir,
                                        cv=10, ana_type=ana_type, cache_dir=report_dir, select_top=3)

            # %% [markdown]
            # MAE = Mean Absolute Error
            # 
            # MSE = Mean Squared Error
            # 
            # RMSE = Rooted Mean Square Error
            # 
            # MAPE = Mean Absolute % Error - the average magnitude of error produced by a model, or how far off predictions are on average. A MAPE value of 20% means that the average absolute percentage difference between the predictions and the actuals is 20%

            # %%
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

            # %%
            modified_top = []
            mtml = []
            for model in top:
                modified_top.append(make_pipeline(steps=[('prep', prep_pipeline), model.steps[-1]]))
                my_top_models = str(model[1:])
                #print(my_top_models)
                my_top_models = my_top_models.split("'")[3]
                mtml.append(my_top_models)
                #print(my_top_models)

            # %%
            top,params,r2_values = finalize_top(X=tr, y=y, top_models=modified_top, grid_param=get_exp_params(),report_dir=report_dir, cv=10, return_params = True)
            print(f'Here are the top params: {params}\n')
            print(f'Here is the R^2 values: {r2_values}\n')
            try:
                one_hot_gs_df.loc[ds,'gs_params'] = str(params)
                one_hot_gs_df.loc[ds,'gs_r2'] = str(r2_values)    
            except:   
                one_hot_gs_df.loc[ds]['gs_params'] = str(params)
                one_hot_gs_df.loc[ds]['gs_r2'] = str(r2_values)
            print(f'Here is the updated row: {one_hot_gs_df.loc[ds]}')
            one_hot_gs_df.to_csv(path_or_buf=one_hot_gs_file,index=True)
            print('Saving Updated GS Dataframe!')
            # %%
            sr = summarize_results(top_models=top, report_dir=report_dir)

            # %%
            scatter_plot = plot_scatter(summary_result=sr, report_dir=report_dir)

            # %%
            mean_imp = mean_importance(top, report_dir=report_dir)

            # %%
            dp_plot(importance=mean_imp,imp_col='mean', model_name='mean', report_dir=report_dir)

            try:
                shutil.rmtree(f'{report_dir}/joblib')
            except:
                pass
            try:
                os.remove(f'{report_dir}/joblib')
            except:
                pass
            
        except Exception as e:
            one_hot_gs_df.to_csv(path_or_buf=one_hot_gs_file,index=True)
            raise Exception(print(e))     
    
one_hot_gs_df['gs_r2'] = one_hot_gs_df['gs_r2'].apply(string_to_numpy_array)
one_hot_gs_df['gs_params'] = one_hot_gs_df['gs_params'].apply(string_to_list_of_dicts)

one_hot_gs_df['best_gs_r2'] = ''
one_hot_gs_df['best_gs_params'] = ''
one_hot_gs_df['best_gs_model'] = str('')
for i,row in enumerate(one_hot_gs_df['gs_r2']):
    best_r2 = max(row)
    best_r2_index = np.argmax(row)
    one_hot_gs_df.loc[i, ('best_gs_r2')] = best_r2
    
    param_list = one_hot_gs_df.loc[i, ('gs_params')]
    best_params = param_list[best_r2_index]
    one_hot_gs_df.loc[i, ('best_gs_params')] = str(best_params)
    
    keys = list(best_params.keys())
    best_model = keys[0].split('_')[0]
    one_hot_gs_df.loc[i, ('best_gs_model')] = best_model
    
one_hot_gs_df.to_csv(path_or_buf=one_hot_gs_file,index=True)
