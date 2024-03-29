import pandas as pd
import numpy as np
import os
from statistics import mean, median
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import random
import seaborn as sns

def get_mm_def(mm_file, extension='.csv'):
    """Get informations about membership matrix saved as `mm_file`.
    
    Parameters
    ----------
    mm_file : string
        Path to the file containing the membership matrix.
        
    extension : string or None, default='.csv'
        File extension.
        
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
    
    if type(extension)==str:
        if extension[0]!='.':
            extension = '.' + extension
        filename = filename.removesuffix(extension)

    if len(attr)>2:
        mm_def_list.append('_'.join(attr[2:]))
    else:
        mm_def_list.append(None)
    mm_def_list.append(filename)
    mm_def_list.append(mm_file)

    mm_def = pd.DataFrame(mm_def_list, index=['Alg', 'Dataset', 'Param', 'Prefix',
                                              'File']).T
    return mm_def
    
def getAll_mm_def(mm_path = './MembMatrix'):
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
        
    if type(mm)==type(None) and file_name==None:
        raise Exception('At least one of `mm` or `file_name` must be specified') 
    if file_name!=None:
        mm = load_mm(file_name)
            
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
           cl0_part, in1cl_part, np.unique(cl_sizes)]
    names = ['n_clusters', 'max_cl_size', 'min_cl_size', 'mean_cl_size', 
             'med_cl_size', 'n_empty_cl', 'n_1_cl', 'n_genes', 'max_cl', 'min_cl',
             'mean_cl', 'med_cl', 'uncl_part', 'uncl_WGCNA', 'in1cl_part', 
             'cl_sizes']
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
        
    mm.columns = [i for i in range(mm.shape[1])]
    return mm
    
def cleanup(mm=None, file_name=None, save=False):
    """Remove the empty cluster from membership matrix.
    
    Parameters
    ----------
    mm : pandas DataFrame, shape=(n_samples, n_clusters), default=None
        Membership matrix of the clustering solution to analyse. If specified and 
        `save`=False, `file_name` is not needed.
    
    file_name : string, default=None
        Path to the file containing the membership matrix to analyse. Must be given
        if `mm` is not specified or if `save`=True.
    
    save : bool, default=False
        If True, the cleaned membership matrix will be saved to `file_name`.
    
    Returns
    -------
    mm : pandas DataFrame, shape=(n_samples, n_clusters)
        Cleaned membership matrix.
    """
    if type(mm)==type(None) and type(file_name)==type(None):
        raise Exception('At least one of `mm` or `file_name` must be specified')
    if type(file_name)!=type(None):
        mm = load_mm(file_name)
    
    cl_sizes = np.sum(mm, axis=0)
    empty_cl = np.where(cl_sizes==0)[0]
    
    if len(empty_cl)==0:
        print("No empty cluster to remove")
        return mm
    
    rm_col = mm.iloc[:,empty_cl].columns
    mm = mm.drop(rm_col, axis=1)
    
    if save:
        if type(file_name)==type(None):
            raise Exception('`file_name` must be passed to save cleaned matrix')
        mm.to_csv(file_name, header=None, index=True)
    return mm

def plotClust(mm, tSNE, title=None):
    # Filter out filled markers and marker settings that do nothing.
    markers = [m for m,func in Line2D.markers.items() if func != 'nothing']
    random.shuffle(markers)

    combColMark = comb(marker=markers, color=sns.color_palette("tab10"))
    random.shuffle(combColMark)
    fsizeinit = plt.rcParams['figure.figsize']
    plt.rcParams['figure.figsize'] = 15,15
    for ic in range(mm.shape[1]):
        inClust = mm.iloc[:,ic]==1
        genes = mm.iloc[list(inClust),ic].index
        coord = tSNE.loc[genes,:]
        plt.scatter(x=coord.loc[:,'x'],y=coord.loc[:,'y'], **combColMark[ic][1], s=30, alpha=0.7)
        
    if title != None:
        plt.title(title)
    plt.show()
    plt.rcParams['figure.figsize'] = fsizeinit
    
    
def comb(*args, **kwargs):
    Li = [[i for i in range(len(a))] for a in args]
    for arg in kwargs:
        Li.append([i for i in range(len(kwargs[arg]))])

    all_combinations = []
    past = []
    future = [len(l) for l in Li]

    for i in range(len(Li)):
        tmp = future.pop(0)
        rep = np.repeat(Li[i], np.prod(future), axis=0)
        nrep = np.prod(past, dtype=int)
        if nrep>1:
            til = np.tile(rep, nrep)
        else:
            til = rep

        all_combinations.append(til)
        past.append(tmp)
    
    comb = np.vstack(all_combinations).T

    comb2 = []
    for ic in range(comb.shape[0]):
        ia = comb[ic, 0:len(args)]
        ikwa = comb[ic, len(args):(len(args)+len(kwargs))]
        
        a = []
        for i in range(len(args)):
            a.append(args[i][ia[i]])
        a = tuple(a)
        kwa = {key:kwargs[key][ival] for key,ival in zip(kwargs, ikwa)}
        
        comb2.append([a,kwa])
    return comb2

def plotCl(mm, tSNE, cl, title=None, s=20, scl=50, alpha=0.5, save=False, fname=None):
    
    plt.scatter(x=tSNE.loc[:,'x'],y=tSNE.loc[:,'y'], c='grey', s=s, alpha=alpha)
    inClust = mm.iloc[:,cl]==1
    genes = mm.iloc[list(inClust),cl].index
    coord = tSNE.loc[genes,:]
    plt.scatter(x=coord.loc[:,'x'],y=coord.loc[:,'y'], c='red', s=scl)    
    if title != None:
        plt.title(title)
    if save:
        plt.savefig(fname, bbox_inches='tight')
    plt.show()

def plotCls(mm, tSNE, cl_, title=None, s=20, scl=50, alpha=0.5, save=False, fname=None):
    plt.scatter(x=tSNE.loc[:,'x'],y=tSNE.loc[:,'y'], c='grey', s=s, alpha=alpha) # plot t-SNE in grey
    for cl in cl_:
        inClust = mm.iloc[:,cl]==1
        genes = mm.iloc[list(inClust),cl].index
        coord = tSNE.loc[genes,:]
        plt.scatter(x=coord.loc[:,'x'],y=coord.loc[:,'y'], s=scl, label=str(cl))    
    plt.legend(title='Cluster N°', fontsize='xx-large', title_fontsize='xx-large')
    if title != None:
        plt.title(title)
    if save:
        plt.savefig(fname, bbox_inches='tight')
    plt.show()