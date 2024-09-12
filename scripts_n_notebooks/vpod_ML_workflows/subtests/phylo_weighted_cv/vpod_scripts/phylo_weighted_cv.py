import numpy as np
import pandas as pd
from Bio import Phylo
import numpy as np
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

def phylo_weighted_cv(distance_matrix, tip_names, n_folds, distance_threshold, relation_mode='leave_out'):
    
    """
    Clusters terminal leaves of a phylogenetic tree into bins based on distance,
    prioritizing phylogenetic relationships.

    Args:
        distance_matrix: A square numpy array representing pairwise distances.
        n_folds: The number of initial bins to create with the most distant points.
        distance_threshold: The minimum allowable distance for bin assignment.
        relation_mode: The method used to deal with nodes which fall beneath the distance threshold.
        
    Returns:
        A dataframe of fold assignments (bin numbers) and mean 'one-to-all' distances for each terminal node of a phylogenetic tree.
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

    # Zip the values and labels together
    zipped_pairs = zip(avg_distances, unassigned_leaves)

    # Sort the zipped pairs in descending order based on the values
    sorted_pairs = sorted(zipped_pairs, reverse=True)  

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




