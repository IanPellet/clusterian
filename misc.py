import pandas as pd
import numpy as np
import os
from statistics import mean, median

def get_mm_def(mm_file):
    """Get informations about membership matrix saved as `mm_file`.
    
    Parameters
    ----------
    mm_file : string
        Path to the file containing the membership matrix.
        
    Returns
    -------
    mm_def : pandas.DataFrame, shape=(1, 5)
        Data frame with fields : 
            - `Alg` : the algorithm used,
            - `Dataset` : the data set used,
            - `Param` : parameters values for the algorithm,
            - `Prefix` : file prefix as `Alg_Dataset_Param`,
            - `File` : path to the membership matrix file.
    """
    filename = mm_file.split('/')[-1]
    attr = filename.split('_')
    
    mm_def_list = [attr[0], attr[1]]

    if len(attr)>2:
        mm_def_list.append('_'.join(attr[2:])[:-4])
    else:
        mm_def_list.append(None)
    mm_def_list.append(filename[:-4])
    mm_def_list.append(mm_file)

    mm_def = pd.DataFrame(mm_def_list, index=['Alg', 'Dataset', 'Param', 'Prefix',
                                              'File']).T
    return mm_def
    
def getAll_mm_def(mm_path = './ModuleMatrix/'):
    """Get informations about membership matrices in `mm_path`.
    
    Parameters
    ----------
    mm_path : string
        Path to the directory containing all the membership matrices.
        
    Returns
    -------
    all_mm_def : pandas.DataFrame, shape=[n_mm, 5]
        Data frame with each row containing information about one cluster solution
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
    mm_files = [os.path.join(mm_path,f) for f in mm_file_names]
    
    all_mm_def_list = []
    for i in range(len(mm_files)):
        i_def = get_mm_def(mm_files[i])
        i_def.index = [i]
        all_mm_def_list.append(i_def)

    all_mm_def = pd.concat(all_mm_def_list)
    return all_mm_def

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

def summary(mm=None, file_name=None):
    """Get a summary of cluster number, size, composition.
    
    Parameters
    ----------
    mm : pandas DataFrame, shape=(n_samples, n_clusters), default=None
        Membership matrix of the clustering solution to analyse. If specified, 
        `file_name` is not needed.
    
    file_name : string, default=None
        Path to the file containing the membership matrix to analyse. Must be given
        only if `mm` is not specified.
    
    Returns
    -------
    mm_summary : pandas DataFrame, shape=(,14)
        Contains the values for :
            - `n_clusters` : number of clusters
            - `max_cl_size` : maximal cluster size (# of samples)
            - `min_cl_size` : minimal cluster size (# of samples)
            - `mean_cl_size` : mean cluster size (# of samples)
            - `med_cl_size` : median cluster size (# of samples)
            - `n_empty_cl` : number of empty clusters
            - `n_1_cl` : number of clusters containing only one sample
            - `n_genes` : total number of samples
            - `max_cl` : maximal number of cluster one sample is part of (# of
            clusters)
            - `min_cl` : minimal number of cluster one sample is part of (# of
            clusters)
            - `mean_cl` : mean number of cluster one sample is part of (# of
            clusters)
            - `med_cl` : median number of cluster one sample is part of (# of
            clusters)
            - `uncl_part` : part of unclustered samples (%)
            - `uncl_WGCNA` : part of unclustered samples for algo WGCNA (%)
            - `in1cl_part` : part of samples in only one cluster (%)
    """
        
    if mm==None and file_name==None:
        raise('At least one of `mm` or `file_name` must be specified')
        
    if file_name!=None:
        mm = load_mm(file_name)
    elif mm!=None and method==None:
        raise('If `file_name` is not specified, `method` must be specified')
            
    n_clusters = mm.shape[1]
    cl_sizes = np.sum(mm, axis=0)
    max_cl_size = int(max(cl_sizes))
    min_cl_size = int(min(cl_sizes))
    mean_cl_size = mean(cl_sizes)
    med_cl_size = median(cl_sizes)
    n_empty_cl = int(sum(cl_sizes==0))
    n_1_cl = int(sum(cl_sizes==1))
    
    n_genes = mm.shape[0]
    cl_per_genes = np.sum(mm, axis=1)
    max_cl = int(max(cl_per_genes))
    min_cl = int(min(cl_per_genes))
    mean_cl = mean(cl_per_genes)
    med_cl = median(cl_per_genes)
    uncl_part = sum(cl_per_genes==0)/n_genes
    cl0_part = cl_sizes.iat[0]/n_genes
    in1cl_part = sum(cl_per_genes==1)/n_genes
    
    
    ppt = [n_clusters, max_cl_size, min_cl_size, mean_cl_size, med_cl_size, 
           n_empty_cl, n_1_cl, n_genes, max_cl, min_cl, mean_cl, med_cl, uncl_part,
           cl0_part, in1cl_part]
    names = ['n_clusters', 'max_cl_size', 'min_cl_size', 'mean_cl_size', 
             'med_cl_size', 'n_empty_cl', 'n_1_cl', 'n_genes', 'max_cl', 'min_cl',
             'mean_cl', 'med_cl', 'uncl_part', 'uncl_WGCNA', 'in1cl_part']
    mm_summary = pd.DataFrame(ppt, index=names).T
    return mm_summary

def summary_all(mm_path='./MembMatrix'):
    """Get a summary of cluster number, size, composition for all mm in `mm_path`.
    
    Parameters
    ----------    
    mm_path : string, default='./MembMatrix'
        Path to the directory containing all the membership matrix to analyse.
    
    Returns
    -------
    all_summary : pandas DataFrame, shape=(n_mm,14)
        Each row represents one clustering solution.
        Contains the values for :
             - `Alg` : the algorithm used,
            - `Dataset` : the data set used,
            - `Param` : parameters values for the algorithm,
            - `Prefix` : file prefix as `Alg_Dataset_Param`,
            - `File` : path to the membership matrix file.
            - `n_clusters` : number of clusters
            - `max_cl_size` : maximal cluster size (# of samples)
            - `min_cl_size` : minimal cluster size (# of samples)
            - `mean_cl_size` : mean cluster size (# of samples)
            - `med_cl_size` : median cluster size (# of samples)
            - `n_empty_cl` : number of empty clusters
            - `n_1_cl` : number of clusters containing only one sample
            - `n_genes` : total number of samples
            - `max_cl` : maximal number of cluster one sample is part of (# of
            clusters)
            - `min_cl` : minimal number of cluster one sample is part of (# of
            clusters)
            - `mean_cl` : mean number of cluster one sample is part of (# of
            clusters)
            - `med_cl` : median number of cluster one sample is part of (# of
            clusters)
            - `uncl_part` : part of unclustered samples (%)
            - `in1cl_part` : part of samples in only one cluster (%)
    """
    all_def = getAll_mm_def(mm_path)
    all_sum_list = []
    fname_col = all_def.columns[-1]
    for i in all_def.index:
        i_sum = summary(file_name=all_def.loc[i,fname_col])
        i_sum.index = [i]
        all_sum_list.append(i_sum)

    all_sum_ppt = pd.concat(all_sum_list)
    all_summary = pd.concat([all_def,all_sum_ppt], axis=1)
    all_summary.loc[all_summary.iloc[:,0]=='WGCNA','uncl_part'] = all_summary.loc[all_summary.iloc[:,0]=='WGCNA','uncl_WGCNA']
    all_summary = all_summary.drop('uncl_WGCNA', axis=1)
    return all_summary

def load_mm(file_name):
    """Load membership matrix.
    
    Parameters
    ----------
    file_name : string
        Path to the file containing the membership matrix to load. 
    
    Returns
    -------
    mm : pandas DataFrame, shape=(n_samples, n_clusters)
        Loaded membership matrix
    """
    mm = pd.read_csv(file_name, header=None, index_col=0)
    if type(mm.index[0])!=str or len(mm.index[0])<2:
        mm = pd.read_csv(file_name, header=0, index_col=0)        
    return mm
    
#def cleanup(mm=None, file_name=None):
   