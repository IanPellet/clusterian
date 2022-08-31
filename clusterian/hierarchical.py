import pandas as pd
#import numpy as np

#ward_g2v = cluster.ward_tree(g2v)
#ward_child = ward_g2v[0]
#n_samples = g2v.shape[0]
#genes = g2v.index

def ward_tree_2_label_mat(ward_child, genes):
    """Create a label matrix from ward tree children array.
    
    Create a matrix containing all possible clustering solutions from the Ward tree
    for all number of clusters `k` in `[1, n_samples]`.
    
    Parameters :
    ----------
    ward_child : int 2D array, shape=(n_samples-1, 2)
        Chilrend from Ward tree computation, each row gives the nodes being merged 
        into a new node at each iteration.
        
    genes : sting array, shape=(n_samples,)
        List of all genes in the dataset.
        
    Returns :
    -------
    lab_df : pandas DataFrame, shape=(n_samples, n_samples)
        Data frame containing all possible solutions for the different cuts of the 
        Ward tree in column. Column names gives the number of clusters in the 
        clustering solution stored as a label list (corresponding column content).
    """
    n_samples = len(genes) 
    # initialy 1 gene is 1 cluster in itself
    n_cl = n_samples
    lab = pd.Series([i for i in range(n_samples)],index=genes) 

    # number of clusters at each iteration
    n_cl = [n_samples-i for i in range(len(ward_child)+1)] 
    # matrix to fill with labels at each iteration
    lab_df = pd.DataFrame(0, index=genes, columns=n_cl) 
    lab_df.iloc[:,0] = lab 
    for i in range(len(ward_child)):
        a,b = ward_child[i] 
        # at iteration i, clusters a and b are merged into one cluster
        new_node = n_samples + i # label of the new cluster
        # the new label list is inialised as equal to the last one
        new_lab = lab_df.iloc[:,i].copy(deep=True) 
        
        # merged cluster lables are modified
        new_lab[lab_df.iloc[:,i]==a] = new_node
        new_lab[lab_df.iloc[:,i]==b] = new_node
        lab_df.iloc[:,i+1] = new_lab # new label list is added to the matrix

    return lab_df