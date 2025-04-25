import numpy as np
import pandas as pd
import random
from Bio import Phylo
import seaborn as sns
import matplotlib.pyplot as plt
import random

def get_dist_matrix_from_tree(tree_file):
    """
    Performs phylogenetic block k-fold cross-validation.

    Args:
        tree_file (str): Path to the phylogenetic tree file (Newick format).
        fold_list (list): List of the number of folds to test.

    Returns:
        list: List of dictionaries, each containing fold assignments and distance metrics.
    """

    # Load phylogenetic tree
    tree = Phylo.read(tree_file, "newick")
    # Get tip names
    tip_names = [terminal.name for terminal in tree.get_terminals()]
    #dist_matrix = tree.distance_matrix().values

    #Create distance matrix from the tree
    dist_matrix = np.zeros((len(tree.get_terminals()), len(tree.get_terminals())))
    for i, terminal1 in enumerate(tree.get_terminals()):
        for j, terminal2 in enumerate(tree.get_terminals()):
            if i < j:  # Only calculate upper triangle to avoid redundancy
                dist_matrix[i, j] = tree.distance(terminal1, terminal2)
                dist_matrix[j, i] = dist_matrix[i, j]  # Mirror for symmetry
    
    #tips_for_tr = []
    #for names in tip_names:
    #    tips_for_tr.append(f'd_{names}')
    #print(tips_for_tr)
    #dist_df =  pd.DataFrame(dist_matrix, index=tip_names, columns=tip_names)

            
    return dist_matrix, tip_names

def percentile_threshold(distance_matrix, percentile=5):
    """Calculates a distance threshold based on a given percentile of pairwise distances."""
    distances = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]  # Extract upper triangle
    return np.percentile(distances, percentile)

def farthest_points(distance_matrix, k=10):
    n = distance_matrix.shape[0]

    # Calculate average distances (initially to all points)
    
    avg_distances = np.mean(distance_matrix, axis=1)
    avg_distances_hold = np.mean(distance_matrix, axis=1)
    #print(f'These are the average distances for each leaf:\n{avg_distances}')
    
    # Initialize cluster (here, we start with the point with max avg distance)
    point_list = [np.argmax(avg_distances)]
    avg_distances = np.delete(avg_distances, np.argmax(avg_distances))

    while len(point_list) < k:
        # Find the point with the highest average distance to existing cluster members
        max_avg_dist = -1
        best_point = -1
        for i in range(n):
            if i in point_list:
                continue
            avg_dist_to_cluster = np.mean(distance_matrix[i, point_list])
            if avg_dist_to_cluster > max_avg_dist:
                max_avg_dist = avg_dist_to_cluster
                best_point = i
        
        # Add the best point to the cluster
        point_list.append(best_point)
        avg_distances = np.delete(avg_distances, best_point)

    #print(f'Length of avg. distances after intial point assignment is:{len(avg_distances)}')

    return point_list, avg_distances, avg_distances_hold

def phylo_weighted_cv(distance_matrix, tip_names, n_folds, distance_threshold, relation_mode='leave_out', shuffle_unassigned=False, reverse_sort=True):
    
    """
    Clusters terminal leaves of a phylogenetic tree into bins based on distance,
    prioritizing phylogenetic relationships.

    Args:
        distance_matrix (np.ndarray): Square matrix of pairwise distances between leaves.
        n_folds (int): The desired number of folds (bins).
        tip_names (list): List of names corresponding to the rows/columns of the distance_matrix.
        distance_threshold (float): Minimum distance required for assigning a leaf to a bin.
        relation_mode (str, optional): Strategy for assigning leaves that don't meet
                                       the threshold ('random', 'merge', 'max_mean', 'leave_out').
                                       Defaults to 'leave_out'.
        reverse_sort (bool, optional): If True and not shuffling, process unassigned leaves
                                       with larger avg_distances first. Defaults to True.
        shuffle_unassigned (bool, optional): If True, shuffle the order of processing
                                             unassigned leaves randomly, ignoring avg_distances sort.
                                             Defaults to False.
    """

    n_leaves = distance_matrix.shape[0]
    class_assignments = np.full(n_leaves, -1, dtype=int)  # Initialize as -1 (unassigned)

    # 1. Initialize Bins with Most Distant Points:
    initial_points, avg_distances, avg_distances_hold = farthest_points(distance_matrix=distance_matrix , k=n_folds)
    #print(f'These are the intial points: {initial_points}')
    #print(f'This is the length of the intial points: {len(initial_points)}')

    for i, idx in enumerate(initial_points):
        class_assignments[idx] = i

    # 2. Iteratively Add Points to Bins:
    unassigned_leaves = np.where(class_assignments == -1)[0]

    if shuffle_unassigned == True:
        random.shuffle(unassigned_leaves)
        sorted_avg_distances = avg_distances
    else:  
        # Zip the values and labels together
        zipped_pairs = zip(avg_distances, unassigned_leaves)

        # Sort the zipped pairs in descending order based on the values
        if reverse_sort == True:
            sorted_pairs = sorted(zipped_pairs, reverse=True)
        else:
            sorted_pairs = sorted(zipped_pairs, reverse=False)  

        # Unzip the sorted pairs back into separate lists
        sorted_avg_distances, unassigned_leaves = zip(*sorted_pairs)
    
    #print(f'These are the unassigned leaves: {unassigned_leaves}')
    #print(f'This is the length of unassigned leaves: {len(unassigned_leaves)}')

    for leaf_idx in unassigned_leaves:
        mean_distances = [np.mean(distance_matrix[leaf_idx, class_assignments == bin_num])
                          for bin_num in range(n_folds)]
        #print(f'Just checking that these are distance: {mean_distances}')
        # Find best bin (highest mean distance above threshold)
        best_bin = np.argmax(mean_distances)
        if mean_distances[best_bin] >= distance_threshold:
            # Check individual distances within the best bin
            if all(distance_matrix[leaf_idx, class_assignments == best_bin] >= distance_threshold):
                class_assignments[leaf_idx] = best_bin
            else:
                # Try next best bin until a suitable bin is found or none exist
                sorted_bins = np.argsort(mean_distances)[::-1]  # Descending order
                #print(f'This is the sorted bins by mean distance: {sorted_bins}')
                for bin_idx in sorted_bins[1:]:  # Skip the already checked best bin
                    if all(distance_matrix[leaf_idx, class_assignments == bin_idx] >= distance_threshold):
                        class_assignments[leaf_idx] = bin_idx
                        break
        if class_assignments[leaf_idx] == -1:
            if relation_mode == 'random':
                sorted_bins = np.argsort(mean_distances)[::-1]  # Descending order
                best_bin = random.choice(sorted_bins)
                class_assignments[leaf_idx] = best_bin
                #print(f'The bin chosen using the random method is {best_bin}')            
            elif relation_mode == 'merge':
                #best_bin = np.argmin(mean_distances)
                min_distances = [min(distance_matrix[leaf_idx, class_assignments == bin_num])
                          for bin_num in range(n_folds)]
                min_dist = min(min_distances)
                best_bin = min_distances.index(min_dist)  
                #print(f'The best bin using the merge method is {best_bin}')            
                class_assignments[leaf_idx] = best_bin
            elif relation_mode == 'max_mean':
                class_assignments[leaf_idx] = best_bin
            else:
                #if relation_mode is set to 'leave_out' we simply keep the fold assignment as -1 
                pass


    #tip_to_fold = {tip_name: {'Fold' : assignment, 'Mean_Distance' : dist} for tip_name, assignment, dist in zip(tip_names, class_assignments, avg_distances_hold)}
    tip_to_fold_dict = {'Tip_Names' : tip_names , 'Fold' : class_assignments, 'Mean_Distance' : avg_distances_hold}
    tip_to_fold_df = pd.DataFrame(tip_to_fold_dict)
    tip_to_fold_df = tip_to_fold_df.set_index('Tip_Names')
    #tip_to_dist = dict(zip(tip_names, ))
    return tip_to_fold_df

def plot_phylo_cv_line_graphs(report_dir, results_file, atts_of_intrst=['R2', 'MAE', 'MAPE', 'MSE', 'RMSE']):
    
    results_df = pd.read_csv(results_file)
    model_att = results_df[results_df['Model'] != 'lr']
    # Get unique relation handling methods
    relation_handling_methods = model_att['Relation_Handling'].unique()
    for attr in atts_of_intrst:
        for method in relation_handling_methods:
            # Filter for 'leave_out' method and drop NaN values
            method_df = model_att[model_att['Relation_Handling'] == method].dropna(subset=attr)
            # Create the plot using Seaborn
            plt.figure(figsize=(16, 8))  # Increase width for better clarity
            plt.rcParams['font.family'] = 'Century Gothic'
            sns.lineplot(data=method_df, x='Threshold', y=attr, hue='Model')

            # Add labels and title with larger font size
            plt.xlabel(f'Percentile Threshold\n(Using "{method}" Threshold Handling)', fontsize=20)
            plt.ylabel(attr, fontsize=20)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)

            # Place the legend in the bottom left corner
            plt.legend(loc='lower left', fontsize=16)

            # Show the plot
            plt.tight_layout()

            # Save the figure as an SVG file
            plt.savefig(f'./{report_dir}/{attr}_{method}_handling_performance_trend.svg', format='svg')  # You can change the filename if needed

            plt.show()

def plot_phylo_cv_indv_model_graphs(report_dir, results_file):

    results_df = pd.read_csv(results_file)
    results_df = results_df[results_df['Model'] != 'lr']
    model_list = results_df['Model'].unique().tolist()

    for model in model_list:
        # Filter for method and drop NaN values
        model_df = results_df[results_df['Model'] == model].dropna(subset=['R2'])
        model_df = model_df[model_df['R2'] >= 0]

        # Create the plot using Seaborn
        plt.figure(figsize=(16, 8))  # Increase width for better clarity
        sns.lineplot(data=model_df, x='Threshold', y='R2', hue='Relation_Handling')

        # Add labels and title with larger font size
        plt.xlabel(f'Percentile Threshold\n(All Methods of Threshold Handling for {model} Model)', fontsize=20)
        plt.ylabel('R^2', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        # Place the legend in the bottom left corner
        plt.legend(loc='lower left', fontsize=16)

        # Show the plot
        plt.tight_layout()

        # Save the figure as an SVG file
        plt.savefig(f'./{report_dir}/{model}_performance_vs_handling.svg', format='svg')  # You can change the filename if needed

        plt.show()

def plot_phylo_cv_bar_graphs(report_dir, results_file, atts_of_intrst=['R2', 'MAE', 'MAPE', 'MSE', 'RMSE']):

    results_df = pd.read_csv(results_file)
    model_att = results_df[results_df['Model'] != 'lr']
    
    for attr in atts_of_intrst:
        # Group by `Relation_Handling` and `Model`, calculate the mean of `R2`, and reset the index
        model_att = results_df.groupby(['Relation_Handling', 'Model'])[attr].mean().reset_index()
        model_att = model_att[model_att['Model'] != 'lr']
        # Get unique relation handling methods
        relation_handling_methods = model_att['Relation_Handling'].unique()

        # Create a separate bar plot for each relation handling method
        for method in relation_handling_methods:
            # Filter data for the current method
            method_data = model_att[model_att['Relation_Handling'] == method]
            method_data = method_data.sort_values([attr], ascending=[False])

            # Create the bar plot
            plt.figure(figsize=(10, 6))  # Adjust figure size as needed
            plt.rcParams['font.family'] = 'Century Gothic'
            sns.barplot(data=method_data, x='Model', y=attr)

            # Add labels and title
            plt.xlabel('Model')
            plt.ylabel(f'Mean {attr}')
            plt.title(f'Mean {attr} for Relation Handling Method: {method}')
            plt.xticks(rotation=45)  # Rotate x-axis labels if needed

            # Show the plot
            plt.tight_layout()

            # Save the figure as an SVG file
            plt.savefig(f'./{report_dir}/avg_{attr}_perf_{method}_handling_bar_graph.svg', format='svg')  # You can change the filename if needed

            # Display the plot
            plt.show()
            
def plt_fold_phylo_distributions(tip_to_fold, handeling_method, threshold=5, n_folds=10, model='wt'):
    tip_to_fold_copy = tip_to_fold.copy()
    if handeling_method == 'leave_out':
        rounds = 2
    else:
        rounds = 1
    
    drop_unused_fold = False
    for round in range(rounds):
        if drop_unused_fold == True:
            tip_to_fold_copy = tip_to_fold_copy[tip_to_fold_copy['Fold']!=-1]
        # Determine the number of unique labels for subplot layout
        num_labels = tip_to_fold_copy['Fold'].nunique()
        num_cols = 1  # You can adjust the number of columns as desired
        num_rows = int(np.ceil(num_labels / num_cols))  # Calculate rows based on columns
        # Generate Color Palette
        colors = sns.color_palette("hls", num_labels)
        # Plotting with Subplots
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(8, 12), sharex=True, sharey=True)  # Create subplot grid
        plt.rcParams['font.family'] = 'Century Gothic'
        fig.suptitle('Density Plot of Phylogenetic Distances by Fold Using Phylo-Weighted CV')

        # Flatten the axes array for easier iteration
        axes = axes.flatten()
        folds = tip_to_fold_copy['Fold'].unique().tolist()
        folds.sort()
        # Calculate the mean of the first class
        mean_class = tip_to_fold_copy['Mean_Distance'].max()
        # Plot each label's density on a separate subplot
        for i, label in enumerate(folds):
            ax = axes[i]  # Get the current subplot axis
            label_data = tip_to_fold_copy[tip_to_fold_copy['Fold'] == label]['Mean_Distance']
            sns.kdeplot(label_data, fill=True, alpha=0.5,color=colors[i], ax=ax)
            # Add label text centered on the x-axis at the mean of Mean_Distance
            ylim = ax.get_ylim()
            ax.text(mean_class, ylim[0] + (ylim[1] - ylim[0])/2, f'Fold {label}',
                    horizontalalignment='right', verticalalignment='center', fontsize=15)


        # Set shared labels only once (for the first subplot in each column and row)
        for ax in axes[:]:
            ax.set_ylabel('Density', fontsize=12)
        
        axes[-1].set_xlabel('Mean Distance', fontsize=12)
            
        # Turn off unused subplots
        for i in range(num_labels, len(axes)):
            axes[i].set_axis_off()

        plt.tight_layout()
        if drop_unused_fold == True:
            plt.savefig(f'./{model}_percentile{threshold}_{handeling_method}_cv_{n_folds}fold_phylo_distributions_dropneg.svg', format='svg')  # You can change the filename if needed
        else:
            plt.savefig(f'./{model}_percentile{threshold}_{handeling_method}_cv_{n_folds}fold_phylo_distributions.svg', format='svg')  # You can change the filename if needed
        
        plt.show()
        
        drop_unused_fold = True