import pandas as pd
import numpy as np
import os

def getAll_mm_def(mm_path = './ModuleMatrix/'):
    """Get informations about membership matrices in `mm_path`.
    
    Parameters
    ----------
    mm_path : string
        Path to the directory containing all the membership matrices.
        
    Returns
    -------
    all_mm : pandas.DataFrame, shape=[n_mm, 5]
        Data frame with each row cotaining information about one cluster solution
        found in `mm_path`. The fields are : 
            - `Alg` : the algorithm used,
            - `Dataset` : the data set used,
            - `Param` : parameters values for the algorithm,
            - `Prefix` : file prefix as `Alg_Dataset_Param`,
            - `File` : path to the membership matrix file.
    """
    mm_file_names = []
    #dir_in_mm = []
    for (dirpath, dirnames, filenames) in os.walk(mm_path):
        mm_file_names.extend(filenames)
        #dir_in_mm.extend(dirnames)
        break

    algo = []   
    dataset = []
    parm = []
    prefix = []
    for filename in mm_file_names:
        #path = os.path.join(mm_path,filename)
        attr = filename.split('_')
        algo.append(attr[0])
        dataset.append(attr[1])
        if len(attr)>2:
            parm.append('_'.join(attr[2:])[:-4])
        else:
            parm.append(None)
        prefix.append(filename[:-4])

    all_mm = pd.DataFrame([algo, dataset, parm, prefix, 
                           [os.path.join(mm_path,f) for f in mm_file_names]], 
                          index=['Alg', 'Dataset', 'Param', 'Prefix', 'File']).T
    return all_mm

def lab2mat(lab, genes):
    """Create membership matrix from label list.
    
    Parameters
    ----------
    lab : int array, shape=(n_samples,)
        List of labels for each sample in the dataset.
        
    genes : string array, shape=(n_samples,)
        List of gene names to use as membership matrix index.
        
    Returns
    -------
    mm : pandas.DataFrame, shape=(n_samples, n_clusters)
        Membership matrix with samples (genes) in row and clusters in column.
    """
    cl = np.unique(lab)
    mm = pd.DataFrame(0, columns=cl, index=genes)
    for i in range(len(lab)):
        mm.loc[genes[i],lab[i]] = 1
    return mm