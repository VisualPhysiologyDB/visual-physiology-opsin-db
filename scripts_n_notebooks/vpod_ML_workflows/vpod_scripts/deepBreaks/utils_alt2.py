import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import pandas as pd
import seaborn as sns
from sklearn.pipeline import Pipeline
import joblib
from typing import List, Tuple

def get_best_aa_prop_combos(dataset):
    aa_props_dict = {
        'wds': ['H1','H2','H3','P2','V','SCT','PKA'],
        'wt': ['H2','H3','P1','NCI','PKA'],
        'vert': ['H2','H3','NCI','SCT','PKB'],
        'wt_vert': ['H2','P2','V','MASS'],
        'inv': ['H1','H3'],
        't1': ['H3','P1','PKB'],
        'wds_mnm': ['H2','H3','NCI','MASS'],
        'wt_mnm': ['H1','H2','H3','NCI','MASS','PKA'],
        'vert_mnm': ['H3','P2','SCT','PKA','PKB'],
        'wt_vert_mnm': ['H2','H3','PKB'],
        'inv_mnm': ['H1','P1','SCT']
    }
    try:
        aa_props = aa_props_dict[dataset]
    except:
        raise Exception(f'Given dataset name, {dataset}, not recognized.\n Here is your list of options: wds, wt, vert, wt_vert, inv, t1, wds_mnm, wt_mnm, vert_mnm, wt_vert_mnm, inv_mnm')

    return aa_props
        
        
def get_models(ana_type, dataset, encoding='hot'):
    if ana_type == 'reg':
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.ensemble import GradientBoostingRegressor
        from sklearn.linear_model import BayesianRidge
        from xgboost import XGBRegressor
        
        if dataset == 'wt':
            if encoding == 'aa_prop':
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.1, max_depth=3, max_features='log2', n_estimators=500, random_state=123),
                }
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=300, gamma=1.0, learning_rate=0.1, max_depth=3, reg_alpha=0, reg_lambda=0)
                }
        elif dataset == 'wt_vert':
            if encoding == 'aa_prop':
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.1, max_depth=3, max_features=None, n_estimators=300, random_state=123),
                }
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=200, gamma=1.0, learning_rate=0.1, max_depth=3, reg_alpha=0, reg_lambda=1.0)
                }
        elif dataset == 'inv':
            if encoding == 'aa_prop':
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.01, max_depth=3, max_features=None, n_estimators=500, random_state=123),
                }
            else:
                models = {
                    'BayesianRidge': BayesianRidge(compute_score=True, fit_intercept=True, alpha_1=0.01, alpha_2=1e-06, lambda_1=1e-06, lambda_2=0.01),
                }
        elif dataset == 'vert':
            if encoding == 'aa_prop':
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=200, gamma=0, learning_rate=0.1, max_depth=5, reg_alpha=0.1, reg_lambda=0.1, colsample_bytree=1.0, subsample=1.0)
                }
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=200, gamma=0, learning_rate=0.1, max_depth=3, reg_alpha=1.0, reg_lambda=1.0)
                }
        elif dataset == 'wds':
            if encoding == 'aa_prop':
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, colsample_bytree=0.8, n_estimators=100, gamma=0, learning_rate=0.2, max_depth=5, reg_alpha=0.1, reg_lambda=0.1, subsample=1.0)
                }
            else:
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.2, max_depth=3, max_features='sqrt', n_estimators=500, random_state=123),
                }
        elif dataset == 't1':
            if encoding == 'aa_prop':
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, colsample_bytree=0.8, n_estimators=300, gamma=0.1, learning_rate=0.1, max_depth=5, reg_alpha=1.0, reg_lambda=1.0, subsample = 0.8)
                }
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, colsample_bytree=0.8, n_estimators=200, gamma=0, learning_rate=0.2, max_depth=3, reg_alpha=1.0, reg_lambda=0, subsample = 0.8)
                }
        elif dataset == 'wds_mnm':
            if encoding == 'aa_prop':
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=300, gamma=1.0, learning_rate=0.1, max_depth=5, reg_alpha=0, reg_lambda=0, colsample_bytree=0.8, subsample=1.0)
                }
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=200, gamma=0, learning_rate=0.1, max_depth=5, reg_alpha=1.0, reg_lambda=0.1, colsample_bytree=0.8, subsample=1.0)
                }
        elif dataset == 'wt_mnm':
            if encoding == 'aa_prop':
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.1, max_depth=3, max_features='sqrt', n_estimators=800, random_state=123),
                }  
            else:
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.1, max_depth=5, max_features='sqrt', n_estimators=500, random_state=123),
                }  
        elif dataset == 'vert_mnm':
            if encoding == 'aa_prop':
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=300, gamma=1.0, learning_rate=0.1, max_depth=5, reg_alpha=0, reg_lambda=0.1, colsample_bytree=0.8, subsample=0.8)
                } 
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=150, gamma=1.0, learning_rate=0.2, max_depth=5, reg_alpha=0.1, reg_lambda=1.0, colsample_bytree=0.8, subsample=1.0)
                }  
        elif dataset == 'wt_vert_mnm':
            if encoding == 'aa_prop':
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=200, gamma=0, learning_rate=0.1, max_depth=5, reg_alpha=1.0, reg_lambda=1.0, colsample_bytree=0.8, subsample=0.8)
                } 
            else:
                models = {
                    'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators=100, gamma=0.1, learning_rate=0.2, max_depth=5, reg_alpha=0.1, reg_lambda=1.0, colsample_bytree=1.0, subsample=1.0)
                }  
        elif dataset == 'inv_mnm':
            if encoding == 'aa_prop':
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.01, max_depth=5, max_features='sqrt', n_estimators=500, random_state=123),
                }
            else:
                models = {
                    'gbr': GradientBoostingRegressor(learning_rate=0.01, max_depth=5, max_features='sqrt', n_estimators=500, random_state=123),
                }  
        else:
            models = {
                'rf': RandomForestRegressor(n_jobs=-1, random_state=123),
                'gbr': GradientBoostingRegressor(random_state=123),
                'BayesianRidge': BayesianRidge(compute_score=True),
                'xgb': XGBRegressor(n_jobs=-1, random_state=123, n_estimators = 300)
            }
    else:
        from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
        from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier
        from sklearn.tree import DecisionTreeClassifier
        from sklearn.linear_model import LogisticRegression
        from xgboost import XGBClassifier
        from lightgbm import LGBMClassifier

        models = {
            'rf': RandomForestClassifier(n_jobs=-1, random_state=123),
            'Adaboost': AdaBoostClassifier(random_state=123),
            'et': ExtraTreesClassifier(n_jobs=-1, random_state=123),
            'lg': LogisticRegression(n_jobs=-1, random_state=123, max_iter=2000),
            'gbc': GradientBoostingClassifier(random_state=123),
            'dt': DecisionTreeClassifier(random_state=123),
            'xgb': XGBClassifier(n_jobs=-1, random_state=123),
            'lgbm': LGBMClassifier(n_jobs=-1, random_state=123)
        }

    return models


def get_scores(ana_type):
    scores = {'cl': {'Accuracy': 'accuracy',
                     'AUC': 'roc_auc_ovr',
                     'F1': 'f1_macro',
                     'Recall': 'recall_macro',
                     'Precision': 'precision_macro'
                     },
              'reg': {
                  'R2': 'r2',
                  'MAE': 'neg_mean_absolute_error',
                  'MSE': 'neg_mean_squared_error',
                  'RMSE': 'neg_root_mean_squared_error',
                  'MAPE': 'neg_mean_absolute_percentage_error'
              }
              }
    return scores[ana_type]

def get_empty_params():
    params = {
        'dt': {'dt__max_depth': [4, 6, 8]},
    }
    return params

def get_simp_params():
    params = {
        # --- NEW: BayesianRidge ---
        'BayesianRidge': {
            'BayesianRidge__compute_score': [True, False], # Whether to compute the log marginal likelihood
            'BayesianRidge__fit_intercept': [True, False], # Whether to fit an intercept
            'BayesianRidge__n_iter': [100,300,500]
        },
        'rf': {
            'rf__max_features': ["sqrt", "log2", None],  # Consider not limiting features
            'rf__n_estimators': [100, 200, 500],  # Vary number of trees
            'rf__max_depth': [None, 10, 20, 30],  # Control tree depth
            'rf__bootstrap': [True, False],  # Use bootstrapping or not
        },
        'Adaboost': {'Adaboost__learning_rate': np.linspace(0.001, 0.1, num=2),
                     'Adaboost__n_estimators': [200, 400]},
        #'gbr': {'gbr__max_depth': range(3, 6), #original is (3, 6)
        #        'gbr__max_features': ['sqrt', 'log2', None],
        #        'gbr__n_estimators': [200, 300, 500], #200,500,800 original
        #        'gbr__learning_rate': [0.001,0.01,0.1]},
        'et': {'et__max_depth': [4, 6, 8, 10, 12, 20],
               'et__n_estimators': [500, 1000]},
        'dt': {'dt__max_depth': [4, 6, 8]},
        'Lasso': {'Lasso__alpha': np.linspace(0.01, 100, num=5)},
        'LassoLars': {'LassoLars__alpha': np.linspace(0.01, 100, num=5)},
        'LassoLars': {'LassoLars__alpha': np.linspace(0.01, 100, num=5)},
        #'xgb': {
        #    'xgb__max_depth': [3, 5, 7],
        #    'xgb__learning_rate': [0.01, 0.1, 0.2],
        #    'xgb__n_estimators': [100, 200, 300],
         # Minimum loss reduction for split
        #},
    }
    return params

def get_params():
    params = {
        'BayesianRidge': {
            'BayesianRidge__alpha_1': [1e-6, 1e-4, 1e-2],  # Prior for alpha
            'BayesianRidge__alpha_2': [1e-6, 1e-4, 1e-2],  # Prior for alpha
            'BayesianRidge__lambda_1': [1e-6, 1e-4, 1e-2], # Prior for lambda
            'BayesianRidge__lambda_2': [1e-6, 1e-4, 1e-2], # Prior for lambda
            'BayesianRidge__compute_score': [True, False], # Whether to compute the log marginal likelihood
            'BayesianRidge__fit_intercept': [True, False], # Whether to fit an intercept
            'BayesianRidge__n_iter': [100,300,500]
        },
        'rf': {
            'rf__max_features': ["sqrt", "log2", None],  # Consider not limiting features
            'rf__n_estimators': [100, 200, 500],  # Vary number of trees
            'rf__max_depth': [None, 10, 20, 30],  # Control tree depth
            'rf__min_samples_split': [2, 5, 10],  # Minimum samples for split
            'rf__min_samples_leaf': [1, 2, 4],  # Minimum samples per leaf
            'rf__bootstrap': [True, False],  # Use bootstrapping or not
        },
        'Adaboost': {'Adaboost__learning_rate': np.linspace(0.001, 0.1, num=2),
                     'Adaboost__n_estimators': [200, 400]},
        'gbr': {'gbc__max_depth': range(3, 12), #original is (3, 6)
                'gbc__max_features': ['sqrt', 'log2', None],
                'gbc__n_estimators': [200, 300, 500, 800], #200,500,800 original
                'gbc__learning_rate': [0.001,0.01,0.1]},
        'et': {'et__max_depth': [4, 6, 8, 10, 12, 20],
               'et__n_estimators': [500, 1000]},
        'dt': {'dt__max_depth': [4, 6, 8]},
        'Lasso': {'Lasso__alpha': np.linspace(0.01, 100, num=5)},
        'LassoLars': {'LassoLars__alpha': np.linspace(0.01, 100, num=5)},
        'lgbm': {
            'lgbm__num_leaves': [31, 63, 127],        # Controls tree complexity
            'lgbm__learning_rate': [0.01, 0.1, 0.2], 
            'lgbm__n_estimators': [100, 200, 300],
            'lgbm__max_depth': [-1, 5, 10],           # -1 means no limit
            'lgbm__reg_alpha': [0, 0.1, 1.0],        # L1 regularization
            'lgbm__reg_lambda': [0, 0.1, 1.0],       # L2 regularization
        },
        # --- NEW: XGBoost ---
        'xgb': {
            'xgb__max_depth': [3, 5, 7],
            'xgb__learning_rate': [0.01, 0.1, 0.2],
            'xgb__n_estimators': [100, 200, 300],
            'xgb__reg_alpha': [0, 0.1, 1.0],
            'xgb__reg_lambda': [0, 0.1, 1.0],
            'xgb__gamma': [0, 0.1, 1.0]              # Minimum loss reduction for split
        },
    }
    return params

def get_exp_params():
    params = {
        'rf': {
            'rf__max_features': ["sqrt", "log2", None],  # Consider not limiting features
            'rf__n_estimators': [100, 200, 500],  # Vary number of trees
            'rf__max_depth': [None, 10, 20, 30],  # Control tree depth
            'rf__min_samples_split': [2, 5, 10],  # Minimum samples for split
            'rf__min_samples_leaf': [1, 2, 4],  # Minimum samples per leaf
            'rf__bootstrap': [True, False],  # Use bootstrapping or not
        },
        'Adaboost': {
            'Adaboost__learning_rate': [0.01, 0.1, 1.0],  # Broader range
            'Adaboost__n_estimators': [50, 100, 200],  # More granular
        },
        'gbr': {
            'gbc__max_depth': [3, 5, 7, 9, 10, 20, None],  # Try unlimited depth
            'gbc__max_features': ["sqrt", "log2", "Auto", None],  # More flexibility
            'gbc__n_estimators': [200, 300, 500, 800],  # Reasonable range
            'gbc__learning_rate': [0.001, 0.01, 0.1, 0.2],  # Common values
            #'gbc__subsample': [0.8, 1.0],  # Consider subsampling 
        },
        'et': {
            'et__max_depth': [None, 10, 20, 30],  # More flexibility
            'et__n_estimators': [100, 200, 300],
            'et__min_samples_split': [2, 5, 10],
            'et__min_samples_leaf': [1, 2, 4],
        },
        'dt': {
            'dt__max_depth': [None, 5, 10, 15, 20], # More options for max_depth
            'dt__min_samples_split': [2, 5, 10],
            'dt__min_samples_leaf': [1, 2, 4],
        },
        'Lasso': {
            'Lasso__alpha': [0.001, 0.01, 0.1, 1, 10, 100]  # Wider logarithmic range
        },
        'LassoLars': {  # Might not be necessary to tune if similar to Lasso
            'LassoLars__alpha': [0.001, 0.01, 0.1, 1, 10, 100]  # Wider logarithmic range
        },
        'lgbm': {
            'lgbm__num_leaves': [31, 63, 127],        # Controls tree complexity
            'lgbm__learning_rate': [0.01, 0.1, 0.2], 
            'lgbm__n_estimators': [100, 200, 300],
            'lgbm__max_depth': [-1, 5, 10],           # -1 means no limit
            'lgbm__subsample': [0.8, 1.0],
            'lgbm__colsample_bytree': [0.8, 1.0],     # Feature subsampling
            'lgbm__reg_alpha': [0, 0.1, 1.0],        # L1 regularization
            'lgbm__reg_lambda': [0, 0.1, 1.0],       # L2 regularization
        },
        # --- NEW: XGBoost ---
        'xgb': {
            'xgb__max_depth': [3, 5, 7],
            'xgb__learning_rate': [0.01, 0.1, 0.2],
            'xgb__n_estimators': [100, 200, 300],
            'xgb__subsample': [0.8, 1.0],
            'xgb__colsample_bytree': [0.8, 1.0],
            'xgb__reg_alpha': [0, 0.1, 1.0],
            'xgb__reg_lambda': [0, 0.1, 1.0],
            'xgb__gamma': [0, 0.1, 1.0]              # Minimum loss reduction for split
        },
        # --- NEW: BayesianRidge ---
        'BayesianRidge': {
            'BayesianRidge__alpha_1': [1e-6, 1e-4, 1e-2],  # Prior for alpha
            'BayesianRidge__alpha_2': [1e-6, 1e-4, 1e-2],  # Prior for alpha
            'BayesianRidge__lambda_1': [1e-6, 1e-4, 1e-2], # Prior for lambda
            'BayesianRidge__lambda_2': [1e-6, 1e-4, 1e-2], # Prior for lambda
            'BayesianRidge__compute_score': [True, False], # Whether to compute the log marginal likelihood
            'BayesianRidge__fit_intercept': [True, False] # Whether to fit an intercept
        },
    }
    return params

def get_color_palette(char_list):
    """
    Generates a dictionary of colors for each character in a given character list.

    Parameters:
    -----------
    char_list : list
        List of characters to generate color palette for.

    Returns:
    --------
    color_dic : dict
        Dictionary of color codes for each character in char_list.
    """

    # Define key and color lists
    key_list = ['A', 'C', 'G', 'R',
                'T', 'N', 'D', 'E',
                'Q', 'H', 'I', 'L',
                'K', 'M', 'F', 'P',
                'S', 'W', 'Y', 'V', 'GAP']
    color_list = ['#0273b3', '#de8f07', '#029e73', '#d55e00',
                  '#cc78bc', '#ca9161', '#fbafe4', '#ece133',
                  '#56b4e9', '#bcbd21', '#aec7e8', '#ff7f0e',
                  '#ffbb78', '#98df8a', '#d62728', '#ff9896',
                  '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
                  '#dbdb8d']  # sns.color_palette('husl', 21)

    # Create dictionary of character-to-color mappings
    color_dic = {}
    for n, key in enumerate(key_list):
        color_dic[key] = color_list[n]
    color_dic['U'] = color_dic['T']

    # Add gray color for mixed combinations
    for let in set(char_list):
        if let not in color_dic and let is not np.nan:
            color_dic[let.upper()] = '#808080'  # hex code for gray color

    return color_dic


def kruskal_test(data: pd.DataFrame, group_col: str, response_var: str) -> float:
    """
    Performs the Kruskal-Wallis H test on a given dataset to determine if there are significant differences between groups
    in terms of a given response variable.

    Args:
    - data (pd.DataFrame): A pandas DataFrame containing the data to be tested.
    - group_col (str): The name of the column in the DataFrame containing the grouping variable.
    - response_var (str): The name of the column in the DataFrame containing the response variable.

    Returns:
    - p (float): The p-value resulting from the Kruskal-Wallis H test.
    """
    k, p = stats.kruskal(*[group[response_var].values for name, group in data.groupby(group_col)])
    return p


def chi2_test(cross_table=None, data=None, group_col=None, response_var=None):
    """Perform a chi-square test for independence of two categorical variables.

    Args:
        cross_table (pandas.DataFrame, optional): A contingency table. Defaults to None.
        data (pandas.DataFrame, optional): Data to create contingency table from. Defaults to None.
        group_col (str, optional): Column name for grouping variable. Defaults to None.
        response_var (str, optional): Column name for response variable. Defaults to None.

    Returns:
        float: p-value of chi-square test.
    """
    if cross_table is None:
        # create contingency table from data
        cross_tab = pd.crosstab(data[response_var], data[group_col])
    else:
        cross_tab = cross_table
    # perform chi-square test
    chi2, p, dof, expected = stats.chi2_contingency(cross_tab)
    return p


def box_plot(data, group_col, response_var, figsize=(3.2, 3.2), ax=None, p=None):
    """
    Create a box plot of response variable stratified by group_col,
    optionally with an associated Kruskal-Wallis test p-value.

    Args:
    - data (pandas DataFrame): The data to be plotted.
    - group_col (str): The name of the column in `data` containing group identifiers.
    - response_var (str): The name of the column in `data` containing the response variable to be plotted.
    - ax (matplotlib Axes object, optional): The Axes object to be plotted on.
    - p (float, optional): The p-value from a Kruskal-Wallis test.

    Returns:
    - ax (matplotlib Axes object): The plotted Axes object.
    """
    tmp = data.loc[:, [group_col, response_var]]
    tmp = tmp.sort_values(by=group_col)  # sort by grouping variable
    n_groups = len(set(data.loc[:, group_col]))
    # Compute p-value if not provided
    if p is None:
        p = kruskal_test(data=tmp, group_col=group_col, response_var=response_var)

    # Create plot
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=350)
    ax.grid(color='gray', linestyle='-', linewidth=0.2, axis='y')
    ax.set_axisbelow(True)

    try:
        sns.boxplot(ax=ax, x=group_col, y=response_var, data=tmp,
                    showfliers=False, dodge=False,
                    width=.6, linewidth=.4,
                    palette=get_color_palette(char_list=tmp.loc[:, group_col]))
        sns.despine(ax=ax)
        sns.stripplot(ax=ax, x=group_col, y=response_var, data=tmp,
                    size=5-np.log(n_groups-1), alpha=0.3, linewidth=.2, hue=group_col,
                    palette=get_color_palette(char_list=tmp.loc[:, group_col]))
        ax.get_legend().remove() 
        ax.set_xlabel('', fontsize=8)
        ax.set_ylabel(response_var, fontsize=8)
        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)

        ax.set_title(group_col + ', P-value of KW test: ' + str(round(p, 3)), fontsize=8)
        ax.set_xlabel('')
    
    except:
        sns.boxplot(ax=ax, x=group_col, y=response_var, data=tmp,
                showfliers=False, dodge=False,
                width=.6, linewidth=.4,)
        sns.despine(ax=ax)
        sns.stripplot(ax=ax, x=group_col, y=response_var, data=tmp,
                    size=5-np.log(n_groups-1), alpha=0.3, linewidth=.2, hue=group_col)
        ax.get_legend().remove() 
        ax.set_xlabel('', fontsize=8)
        ax.set_ylabel(response_var, fontsize=8)
        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)

        ax.set_title(group_col + ', P-value of KW test: ' + str(round(p, 3)), fontsize=8)
        ax.set_xlabel('')

    return ax


def stacked_barplot(cross_table=None, data=None, group_col=None, response_var=None, ax=None, figsize=(3.2, 3.2)):
    """
    Generate a stacked bar plot for categorical data.

    Args:
    - cross_table (pandas.DataFrame): a contingency table of counts for categorical data.
    - data (pandas.DataFrame): the input data.
    - group_col (str): the name of the column in the input data that contains the grouping variable.
    - response_var (str): the name of the column in the input data that contains the response variable.
    - ax (matplotlib.axes.Axes): the target axes object for plotting. If None, a new figure will be created.

    Returns:
    - ax (matplotlib.axes.Axes): the generated stacked bar plot axes object.
    """
    # if cross_table is not provided, generate one from the input data
    if cross_table is None:
        cross_tb = pd.crosstab(data[response_var], data[group_col])
    else:
        cross_tb = cross_table

    # calculate chi-square test p-value
    p = chi2_test(cross_table=cross_tb)

    # get color palette for the plot
    color_dic = get_color_palette(char_list=cross_tb.columns.tolist())

    # Create plot
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=350)

    # generate the stacked bar plot
    cross_tb.plot(kind="bar", stacked=True, rot=0, ax=ax,
                  color=color_dic, width=.3)

    # set plot title and axis labels
    ax.set_title(group_col + ', P-value of Chi-square test: ' + str(round(p, 3)), fontsize=6)
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=6, rotation=90)
    ax.set_ylabel('Counts', fontsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.legend(title=None, fontsize=6)

    return ax


def make_pipeline(steps: List[Tuple[str, object]], cache_dir: str = None) -> Pipeline:
    """
    Creates a scikit-learn pipeline with the given steps and memory cache.

    Parameters:
    -----------
    steps : list of tuples
        The pipeline steps as a list of tuples, where each tuple contains the step name and the corresponding estimator.
    cache_dir : str or None, optional (default=None)
        The directory to use as a memory cache for the pipeline.

    Returns:
    --------
    pipeline : Pipeline object
        The scikit-learn pipeline created with the given steps and memory cache.
    """
    # Create a scikit-learn pipeline with the given steps and memory cache
    pipeline = Pipeline(memory=cache_dir, steps=steps)
    # Return the pipeline
    return pipeline


def save_obj(obj: object, file_name: str) -> str:
    """
    Saves a Python object to a file in the pickle format.

    Parameters:
    -----------
    obj : object
        The Python object to be saved.
    file_name : str
        The name of the file to be created.

    Returns:
    --------
    str
        A string confirming that the object has been saved to the file.
    """
    # Extract the file extension from the file name
    extension = file_name.split('.')[-1]
    # Check if the file extension is '.pkl'
    assert extension == 'pkl', 'File name should be saved as a .pkl file. Please modify your file_name'
    # Save the object to the file in the pickle format using joblib.dump()
    joblib.dump(obj, file_name)
    # Return a confirmation message
    return 'Object saved'


def load_obj(file_name: str) -> object:
    """
    Loads a Python object from a file in the pickle format.

    Parameters:
    -----------
    file_name : str
        The name of the file to be loaded.

    Returns:
    --------
    object
        The Python object loaded from the file.
    """
    # Load the object from the file in the pickle format using joblib.load()
    obj = joblib.load(file_name)
    # Return the loaded object
    return obj
