# %%
# importing deepBreaks libraries 
from vpod_scripts.deepBreaks.utils_alt2 import get_models, get_best_aa_prop_combos, get_scores, get_empty_params, make_pipeline
from vpod_scripts.deepBreaks.preprocessing import MisCare, ConstantCare, URareCare, CustomOneHotEncoder, AminoAcidPropertyEncoder
from vpod_scripts.deepBreaks.preprocessing import FeatureSelection, CollinearCare
from vpod_scripts.deepBreaks.preprocessing import read_data
from vpod_scripts.deepBreaks.models import model_compare_cv, finalize_top, importance_from_pipe, aaprop_importance_from_pipe,  mean_importance, summarize_results
from vpod_scripts.deepBreaks.visualization import dp_plot, dp_aa_prop_plot, plot_imp_model, plot_imp_all
import numpy as np
import pandas as pd
import warnings
import datetime
import os
import shutil 
# %%
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# %%
# defining user params, file pathes, analysis type
print('Beginning script!\n')
#assign your path to folder containing all the datasplits
path = './vpod_1.3_data_splits_2025-10-06_16-50-06'
save_to = './vpod_1.3_all_models'
data_dict = {'wds_mnm':{'meta':'wds_mnm_meta.csv','seq':'wds_mnm_aligned_VPOD_1.3_het.fasta'},
             'wt_mnm':{'meta':'wt_mnm_meta.csv','seq':'wt_mnm_aligned_VPOD_1.3_het.fasta'},
             'wt_vert_mnm':{'meta':'wt_vert_mnm_meta.csv','seq':'wt_vert_mnm_aligned_VPOD_1.3_het.fasta'},
             'inv_mnm':{'meta':'inv_mnm_meta.csv','seq':'inv_mnm_aligned_VPOD_1.3_het.fasta'},
             'vert_mnm':{'meta':'vert_mnm_meta.csv','seq':'vert_mnm_aligned_VPOD_1.3_het.fasta'},
             'wds':{'meta':'wds_meta.tsv','seq':'wds_aligned_VPOD_1.3_het.fasta'},
             'wt':{'meta':'wt_meta.tsv','seq':'wt_aligned_VPOD_1.3_het.fasta'},
             'wt_vert':{'meta':'wt_vert_meta.tsv','seq':'wt_vert_aligned_VPOD_1.3_het.fasta'},
             'inv':{'meta':'inv_meta.tsv','seq':'inv_aligned_VPOD_1.3_het.fasta'},
             'vert':{'meta':'vert_meta.tsv','seq':'vert_aligned_VPOD_1.3_het.fasta'},
            # 't1':{'meta':'Karyasuyama_T1_ops_meta.tsv','seq':'Karyasuyama_T1_ops.fasta'}
            }


all_models_results_files = f'{save_to}/gs_model_report_all_datasets'
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

encoding_methods_list = ['hot','aa_prop']

if 'R2' not in all_model_results_df.columns:
    all_model_results_df['R2'] = None
    all_model_results_df['encoding'] = None
    all_model_results_df['model'] = None
    all_model_results_df['MAE'] = None
    all_model_results_df['MAPE'] = None
    all_model_results_df['MSE'] = None
    all_model_results_df['RMSE'] = None
    
for encoding_method in encoding_methods_list:
    all_model_results_df_copy=all_model_results_df.copy()
    for ds in all_model_results_df.index.to_list():
        #print(all_model_results_df.loc[ds])    
        if all_model_results_df.loc[ds]['R2'] is None or pd.isna(all_model_results_df.loc[ds]['R2']):
            # path to sequences of interest
            seq = all_model_results_df.loc[ds]['seq']
            seqFileName = f'{path}/{seq}'
            # path to corresponding metadata of interest
            meta = all_model_results_df.loc[ds]['meta']
            metaDataFileName = f'{path}/{meta}' 


            # %%
            # making a unique directory for saving the reports of the analysis
            dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

            if encoding_method == 'aa_prop':
                props_to_keep = get_best_aa_prop_combos(ds)
                props_used = ''
                for props in props_to_keep:
                    props_used += props + '_'
                report_dir = str(save_to + '/' + ds + '_' + props_used + 'gs_optimized_' + dt_label)
            else:
                props_to_keep = None
                report_dir = str(save_to + '/' + ds + '_gs_optimized_' + dt_label)

            os.makedirs(report_dir)
            # %%
            #print('reading meta-data')
            # importing metadata
            meta_data = read_data(metaDataFileName, seq_type = None, is_main=False)
            # importing sequences data
            #print('reading fasta file')

            tr = read_data(seqFileName, seq_type = seq_type, is_main=True, gap_threshold=gap_threshold)
            
            try:
                reference_seq = tr.loc['Bovine'].copy()
                ref_seq_name = 'bovine'
                if drop_ref == True:
                    meta_data = meta_data.drop('Bovine')
                #print(bovine)
            except:
                reference_seq = tr.loc['Squid'].copy()
                ref_seq_name = 'squid'
                #print(squid)
            reference_seq.to_csv(path_or_buf= f'{report_dir}/ref_sequence.csv',index = True,mode="w")
            # %%
            tr = tr.merge(meta_data.loc[:, mt],  left_index=True, right_index=True)
            tr.shape

            # %%
            y = tr.loc[:, mt].values
            tr.drop(mt, axis=1, inplace=True)
            #print('Shape of data is: ', tr.shape)
            
            
            # %%
            if encoding_method == 'aa_prop':
                prep_pipeline = make_pipeline(
                    steps=[
                        ('mc', MisCare(missing_threshold=0.05)),
                        ('cc', ConstantCare()),
                        ('aa_prop', AminoAcidPropertyEncoder(props_to_keep = props_to_keep)),
                        ('feature_selection', FeatureSelection(model_type=ana_type, alpha=0.10, keep=False)),
                        ('collinear_care', CollinearCare(dist_method='correlation', threshold=0.01, keep=False))
                    ])
            elif encoding_method == 'one_hot':
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
                                        models_dict=get_models(ana_type=ana_type, dataset=ds, encoding=encoding_method),
                                        scoring=get_scores(ana_type=ana_type),
                                        report_dir=report_dir,
                                        cv=10, ana_type=ana_type, cache_dir=report_dir)
            
            
            #print(report)    
            # Record R^2 values for each model.
            for model_name, r2_value, mae_value, mape_value, mse_value, rmse_value in zip(report.index.to_list(), report["R2"],report["MAE"],report["MAPE"],report["MSE"],report["RMSE"]):
                all_model_results_df_copy.loc[ds]['encoding'] = str(encoding_method)
                all_model_results_df_copy.loc[ds]['model'] = str(model_name)
                all_model_results_df_copy.loc[ds]['R2'] = float(r2_value)
                all_model_results_df_copy.loc[ds]['MAE'] = float(mae_value)
                all_model_results_df_copy.loc[ds]['MAPE'] = float(mape_value)
                all_model_results_df_copy.loc[ds]['MSE'] = float(mse_value)
                all_model_results_df_copy.loc[ds]['RMSE'] = float(rmse_value)
                
            
            if encoding_method == 'hot':
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
            elif encoding_method == 'aa_prop':
                prep_pipeline = make_pipeline(
                    steps=[
                        ('mc', MisCare(missing_threshold=0.05)),
                        ('cc', ConstantCare()),
                        ('aa_prop', AminoAcidPropertyEncoder(props_to_keep = props_to_keep)),
                        ('feature_selection', FeatureSelection(model_type=ana_type, alpha=0.10, keep=True)),
                        ('collinear_care', CollinearCare(dist_method='correlation', threshold=0.01, keep=True))
                    ])
            else:
                raise ValueError('Encoding method not recognized')

            modified_top = []
            mtml = []
            for model in top:
                modified_top.append(make_pipeline(steps=[('prep', prep_pipeline), model.steps[-1]]))
                my_top_models = str(model[1:])
                my_top_models = my_top_models.split("'")[3]
                mtml.append(my_top_models)
                #print(my_top_models)

            top = finalize_top(X=tr, y=y, top_models=modified_top, grid_param=get_empty_params(),report_dir=report_dir, cv=10)
            #summarize the results by extracting feature importance and p-values and grouping correlated features.
            sr = summarize_results(top_models=top, report_dir=report_dir)
            mean_imp = mean_importance(top, report_dir=report_dir)
            
            meta_var='λmax (nm)'            
            if encoding_method == 'hot':
                tr_copy = tr.copy()
                encoded_seqs = prep_pipeline[:4].fit_transform(tr)
                for model in top:
                    model_name = model.steps[-1][0]
                    dp_plot(importance=importance_from_pipe(model),
                            imp_col='standard_value',
                            model_name = model_name, report_dir=report_dir)
                    
                    plot_imp_model(importance=importance_from_pipe(model), 
                            X_train=encoded_seqs, y_train=y, model_name=model_name,
                                meta_var=meta_var, model_type=ana_type, report_dir=report_dir)

                pl = plot_imp_all(final_models=top,
                                X_train=tr, y_train=y,
                                model_type = ana_type,
                                report_dir=report_dir, max_plots=100,
                                figsize=(2.5, 3), meta_var=meta_var)
            elif encoding_method == 'aa_prop':
                for model in top:
                    encoded_seqs = model.named_steps['prep']['aa_prop'].aa_encoded_seqs_
                    model_name = model.steps[-1][0]
                    
                    dp_aa_prop_plot(importance=aaprop_importance_from_pipe(model),
                            imp_col='standard_value',
                            model_name = model_name, report_dir=report_dir, props_to_keep = props_to_keep)
                    
                    plot_imp_model(importance=aaprop_importance_from_pipe(model), 
                            X_train=encoded_seqs, y_train=y, model_name=model_name,
                                meta_var='λmax', model_type=ana_type, report_dir=report_dir)
                    
                pl = plot_imp_all(final_models=top,
                            X_train=encoded_seqs, y_train=y,
                            model_type = ana_type,
                            report_dir=report_dir, max_plots=100,
                            figsize=(2.5, 3), meta_var=meta_var)
                
            from vpod_scripts.translate_candidate_sts import translate_candidate_sts, translate_candidate_sts_aa_props

            if encoding_method == 'hot':
                trans_imp_report = translate_candidate_sts(report_dir, ref_seq_name, reference_seq)
            elif encoding_method == 'aa_prop':
                trans_imp_report = translate_candidate_sts_aa_props(report_dir, reference_seq, ref_seq_name)

                
    all_model_results_df_copy.to_csv(path_or_buf=f'{all_models_results_files}_{encoding_method}.csv',index=True)

            


                    