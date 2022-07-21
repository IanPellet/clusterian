import gseapy as gp
import pandas as pd
import numpy as np
import os
from . import misc
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

def runEnrichr(modules, desc, gene_sets='GO_Biological_Process_2021', out_dir='~/Documents/Clustering/ModuleMatrix/Enrich/'):
    """Run enrichment analysis on all clusters, returns output directory.
    
    A new directory named as `desc` is created in the output durectory to contain the analysis of all the clusters.
    
    Arguments:
    modules -- membership matrix with genes as index
    desc -- prefix to use to name created directories
    Keyword arguments:
    gene_sets -- gene library to use for analysis (default 'GO_Biological_Process_2021')
    out_dir -- output directory (default '~/Documents/Clustering/ModuleMatrix/Enrich/')
    Return value:
    od -- output directory
    """
    od = os.path.join(out_dir, desc)
    os.makedirs(od, exist_ok=True)
    mm_enrich = []
    for i in range(modules.shape[1]):
        # list, dataframe, series inputs are supported
        genes = modules.loc[modules.iloc[:,i]==1,i].index.tolist()
        if len(genes)>1:
            #print(genes)
            try :
                enr = gp.enrichr(gene_list=genes,
                                 gene_sets=gene_sets,
                                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                                 description=desc,
                                 outdir=os.path.join(od, desc+'_'+str(i)),
                                 no_plot=True,
                                 cutoff = 1
                                )
                mm_enrich.append(enr)
            except Exception :
                print("Cluster", i)
            #else :
            #    print("Cluster", i, ":", genes)
            #minAdjPval.append(enr.results.loc[0,"Adjusted P-value"])
        else:
            mm_enrich.append(None)
            #minAdjPval.append(1)
        print("\r",round((i+1)/modules.shape[1]*100),"%",end="\r")
        
    #return minAdjPval
    return od

def getEnriched(enr_dir, score_th=300, genes=None, pval=None):
    """Create and return a DataFrame of all the enriched terms in each cluster.
    
    Parameters :
    ----------
    enr_dir : string
        Path to the directory where to find the subdirectories created by enrichr
        analysis for each cluster.
        
    score_th : int, default=300
        Combined Score threshold above which a term is considered as enriched.
    
    genes : array, shape = (n_genes,)
        List of all genes in the dataset.
        
    pval : float, default = None
        If `pval` if given, `score_th` is not used, instead the adjusted p-value
        threshold is used to evaluate enriched terms.
        
    Returns :
    -------
    enriched : pandas DataFrame, shape(n_enriched_term, 11)
        All the enriched terms in each cluster with : 
            - `Cluster` : correponding cluster
            - `Gene_set` : data base used to run enrichment analysis
            - `Term` : enriched term name
            - `Overlap`
            - `P-value`
            - `Adjusted P-value`
            - `Old P-value`
            - `Old Adjusted P-value`
            - `Odds Ratio`
            - `Combined Score`
            - `Genes` : set of genses in the cluster and in the enriched term
    """
    enriched_ = []
    n_term_cl = []
    for root, dirs, files in os.walk(enr_dir, topdown=False):
        path_dirs = root.split('/')
        if path_dirs[-1][0]=='.':
            files = []
        files_res = []
        for f in files:
            if "enrichr.reports" in f:
                files_res.append(f)

        for name in files_res:
            cl = int(root.split('_')[-1])
            path = os.path.join(root, name)
            full_tab = pd.read_csv(path, sep='\t')
            try :
                if type(pval)==type(None):
                    enr_tab = full_tab.loc[full_tab.loc[:,'Combined Score']>score_th,:]
                else:
                    enr_tab = full_tab.loc[full_tab.loc[:,'Adjusted P-value']<pval,:]
                new_col = pd.Series([cl]*enr_tab.shape[0], index=enr_tab.index,
                                dtype=int)
                enriched_.append(pd.concat([new_col,enr_tab], axis=1))
                n_term_cl.append([cl, len(new_col)])
            except KeyError :
                print('Invalid report for', path)
                #print(full_tab)
            
            
    if len(enriched_)==0:
        return None, None, None, None
    
    enriched = pd.concat(enriched_, axis=0).sort_values(0)
    col = enriched.columns.tolist()
    col[0] = 'Cluster'
    enriched.columns = col
    n_term_cl = pd.DataFrame(n_term_cl, columns=['Cluster',
                                                 'n_enrich']).sort_values('Cluster')
    enr_cl_part = sum(n_term_cl.loc[:,'n_enrich']>0)/n_term_cl.shape[0]
    genes_in_enr = np.unique(';'.join(enriched.loc[:,'Genes']).split(';'))
    if type(genes)!=type(None):
        genes_in_enr_part = len(genes_in_enr)/len(genes)
    return enriched, n_term_cl, enr_cl_part, genes_in_enr_part if type(genes)!=type(None) else genes_in_enr

def runEnrichr_directory(mm_path = './MembMatrix', mm_enrich_dir = 'Enrich'):
    """Run enrichment analysis for all membership matrices.
    
    Keyword arguments:
    mm_path -- Path of the directory containing the membership matrices (default './ModuleMatrix')
    mm_enrich_dir -- Name of the directory where to output results (default 'mm_enrich_dir')
    """
    mm_file_names = []
    dir_in_mm = []
    for (dirpath, dirnames, filenames) in os.walk(mm_path):
        mm_file_names.extend(filenames)
        dir_in_mm.extend(dirnames)
        break
    
    enr_path = os.path.join(mm_path,mm_enrich_dir)
    if mm_enrich_dir not in dir_in_mm:
        os.makedirs(enr_path, exist_ok=True)
    
    enr_dir_content = []
    for (dirpath, dirnames, filenames) in os.walk(enr_path):
        enr_dir_content.extend(dirnames)
        break

    for filename in mm_file_names:
        path = os.path.join(mm_path,filename)
        #attr = filename.split('_')
        #algo = attr[0]
        #dataset = attr[1]
        #if len(attr)>2:
        #    parm = '_'.join(attr[2:])[:-4]

        desc = filename.removesuffix('.csv')

        if desc not in enr_dir_content:
            print('Runing analysis on ', desc) 

            mm = pd.read_csv(path, header=None, index_col=0)
            if type(mm.index[0])!=str or len(mm.index[0])<2:
                mm = pd.read_csv(path, header=0, index_col=0)
            mm.columns = range(0,mm.shape[1])

            print('out_dir :', enr_path)
            where_ = runEnrichr(mm, desc, gene_sets='GO_Biological_Process_2021', 
                                       out_dir= enr_path)
            print('Enrichment results in :', where_)
            
def getAll_enr_def(enr_path = './MembMatrix/Enrich/'):
    """Get informations about enrichment directories in `enr_path`.
    
    Parameters
    ----------
    enr_path : string
        Path to the directory containing all the enrichment results.
        
    Returns
    -------
    all_def : pandas.DataFrame, shape=[n_mm, 5]
        Data frame with each row containing information about one cluster solution
        found in `mm_path`. The fields are : 
            - `Alg` : the algorithm used,
            - `Dataset` : the data set used,
            - `Param` : parameters values for the algorithm,
            - `Prefix` : file prefix as `Alg_Dataset_Param`,
            - `File` : path to the enrichment directory.
    """
    dir_in_path = []
    for (dirpath, dirnames, filenames) in os.walk(enr_path):
        dir_in_path.extend(dirnames)
        break
    enr_dirs = [os.path.join(enr_path,d) for d in dir_in_path if d[0]!='.']
    #return enr_dirs
    
    all_def_list = []
    for i in range(len(enr_dirs)):
        i_def = misc.get_mm_def(enr_dirs[i], extension=None)
        #print(enr_dirs[i], ':', i_def, '\n')
        i_def.index = [i]
        all_def_list.append(i_def)

    all_def = pd.concat(all_def_list)
    return all_def

def evalEnriched_all(datasets, enr_path='./MembMatrix/Enrich', score_th = 300, pval=None):
    """Get cluster evaluation according to enrichment results in `enr_path`.

    Parameters
    ----------
    datasets : dict {string : pandas.DataFrame}
        Dictionnary containing the different datasets used in the study referened by
        its name as defined in all enrichment file names.
    
    enr_path : string
        Path to the directory containing all the enrichment results.
        
    score_th : int, default=300
        Combined Score threshold above which a term is considered as enriched.
        
    pval : float, default = None
        If `pval` if given, `score_th` is not used, instead the adjusted p-value
        threshold is used to evaluate enriched terms.
        
    Returns
    -------
    evalEnr_all : pandas.DataFrame, shape=[n_dir, 8]
        Data frame with each row containing information about one cluster solution
        found in `enr_path` and the evaluation metrics computed. The fields are : 
            - `Alg` : the algorithm used,
            - `Dataset` : the data set used,
            - `Param` : parameters values for the algorithm,
            - `Prefix` : file prefix as `Alg_Dataset_Param`,
            - `File` : path to the enrichment directory,
            - `rC` : ratio of clusters with at least one enriched term,
            - `rG` : ratio of genes included in an enriched term,
            - `Zc` : z-score for `rC`,
            - `Zg` : z-score for `rG`,
            - `Z-score` : z-score for `Zc + Zg`,.
    """
    all_enr =  getAll_enr_def(enr_path)
    rC = []
    rG = []
    ri = []
    for i in all_enr.index:
        data = datasets[all_enr.at[i,'Dataset']]
        _, _, rCi, rGi = getEnriched(all_enr.at[i,'File'],score_th,genes=data.index, pval=pval)
        
#        if False:
#            try :
#                _, _, rCi, rGi = getEnriched(all_enr.at[i,'File'],
#                                                        score_th,genes=data.index)
#                if type(rCi)==type(None) or type(rGi)==type(None):
#                    raise Exception('Empty enrichment directory')
#            except KeyError or Exception :
#                print(all_enr.at[i,'File'])
                
        if type(rCi)==type(None) or type(rGi)==type(None):
            rC.append(0)
            rG.append(0)
        else:
            rC.append(rCi)
            rG.append(rGi)
        ri.append(i)
    
    
    Zc = pd.Series((rC-np.mean(rC))/np.std(rC))
    Zg = pd.Series((rG - np.mean(rG))/np.std(rG))
    evalEnr = pd.DataFrame([pd.Series(rC), pd.Series(rG),Zc, Zg, (Zc+Zg)/sqrt(2)], 
                           columns=ri, index=['rC', 'rG','Zc','Zg', 'Z-score']).T
    evalEnr_all = pd.concat([all_enr, evalEnr], axis=1)
    #evalEnr_sorted = evalEnr_all.sort_values('Z', ascending=False)
    #evalEnr_sorted.index = [i for i in range(evalEnr_sorted.shape[0])]
    #return evalEnr_sorted
    return evalEnr_all

def getBest_cluster_list(enr_dir, pval=False):
    """Get the list of most enriched clusters according to combined score.
    
    Parameters
    ----------
    enr_dir : string
        Path to the directory containing the enrichment analysis results for all 
        clusters.
    pval : Boolean, default = False
        If True, the enrichment terms are sorted according to the adjusted p-value
        instead of the combined score.
        
    Returns
    -------
    bestCl_list : int array, shape = (n_clusters,)
        List of cluster index, from most enriched to least.
    """
    enriched,_,_,_ = getEnriched(enr_dir)
    
    if pval:
        bestScore_sort = enriched.sort_values('Adjusted P-value', ascending=True)
    else:
        bestScore_sort = enriched.sort_values('Combined Score', ascending=False)
    items = bestScore_sort.iloc[:,0]
    bestCl_list = list(dict.fromkeys(items))
    return bestCl_list

def getEnrich_cl(enr_dir, cl, score_th=300, pval=None):
    """Get the enrichment analysis results for one cluster.
    
    Parameters
    ----------
    enr_dir : string
        Path to the directory containing the enrichment analysis results for all 
        clusters.
    
    cl : int
        Index of the cluster for which to get the enrichment results.
        
    score_th : int or None, default = 300
        Combined Score threshold above which a term is considered as enriched.
        If `None` all terms are returned.
        
    pval : float, default = None
        If `pval` if given, `score_th` is not used, instead the adjusted p-value
        threshold is used to evaluate enriched terms.
        
    Returns
    -------
    enr_cl : pandas DataFrame,  shape(n_enriched_term, 10)
        All the enriched terms in the cluster with : 
            - `Gene_set` : data base used to run enrichment analysis
            - `Term` : enriched term name
            - `Overlap`
            - `P-value`
            - `Adjusted P-value`
            - `Old P-value`
            - `Old Adjusted P-value`
            - `Odds Ratio`
            - `Combined Score`
            - `Genes` : set of genses in the cluster and in the enriched term
            
    enr_genes : string array, shape=(n_genes_enrich,)
        Only returned if `score_th != None`. List of genes in the contributing to 
        the enriched terms.
        
    cl_genes : string array, shape=(n_genes_cl,)
        List of genes in the cluster.
    """
    enr_cl_dir = enr_dir.split('/')[-1]+'_'+str(int(cl))
    enr_cl_path = os.path.join(enr_dir, enr_cl_dir)

    for root, dirs, files in os.walk(enr_cl_path, topdown=False):
        for f in files:
            if "enrichr.reports" in f:
                break

    enr_cl_file = os.path.join(enr_cl_path,f)
    enr_cl = pd.read_csv(enr_cl_file, sep='\t').sort_values('Combined Score',
                                                            ascending=False)
    
    mm_file_list = enr_dir.split('/')
    mm_file_list.pop(-2)
    mm_file = '/'.join(mm_file_list)+'.csv'
    mm = pd.read_csv(mm_file, index_col=0)
    cl_genes = np.sort(mm.loc[mm.iloc[:,cl]==1,:].index)
    
    if type(score_th)==type(None) and type(pval)==type(None):
        return enr_cl, cl_genes
    
    if type(pval)==type(None):
        enr_cl = enr_cl.loc[enr_cl.loc[:,'Combined Score']>score_th,:]
    else:
        enr_cl = enr_cl.loc[enr_cl.loc[:,'Adjusted P-value']<pval,:]
    enr_genes = np.unique(';'.join(enr_cl.iloc[:,-1]).split(';'))
    return enr_cl, enr_genes, cl_genes

def getBest_cluster_list_all(enr_dir='./MembMatrix/Enrich/', pval=None, cover=False):
    """Get the list of most enriched clusters for all clustering solutions.
    
    Parameters
    ----------
    enr_dir : string
        Path to the directory containing the enrichment analysis results for all 
        clustering solutions.
    pval : float, default = None
        If `pval` if given, `score_th` is not used, instead the adjusted p-value
        threshold is used to evaluate enriched terms.
        
    Returns
    -------
    bestCl_dict : dict {string:int array}, length = n_solutions
        List of cluster index, from most enriched to least for each clustering 
        solution.
    """
    all_def = getAll_enr_def(enr_dir)
    
    bestCl_dict = {}
    
    for i in all_def.index:
        Prefix = all_def.at[i,'Prefix']
        if cover:
            bestCl_dict[Prefix] = getBest_cluster_cover(all_def.at[i,'File'], 
                                                   pval=pval)
        else:
            bestCl_dict[Prefix] = getBest_cluster_list(all_def.at[i,'File'], 
                                                   pval=pval)
        #print(Prefix,':',bestCl_dict[Prefix])
        
    return bestCl_dict

def getBest_cl_enr(keys, bestCl_dict, n_cl=1, enr_dir='./MembMatrix/Enrich/',
                   score_th=1000, pval=None):
    """Get the enriched terms of the `n_cl` best clusters for solutions in `keys`.
    
    Parameters
    ----------
    keys : string array, shape = (n_solutions,)
        Prefixes of the files for the solutions to analyse.
    bestCl_dict : dict {string:int array}, length = n_solutions
        List of cluster index, from most enriched to least for each clustering 
        solution. 
    n_cl : int, default = 1
        Number of clusters to analyse per solution.
    enr_dir : string
        Path to the directory containing the enrichment analysis results for all 
        clustering solutions.    
    score_th : int or None, default = 300
        Combined Score threshold above which a term is considered as enriched.
        If `None` all terms are returned.
    pval : float, default = None
        If `pval` if given, `score_th` is not used, instead the adjusted p-value
        threshold is used to evaluate enriched terms.
       
    Returns
    -------
    enr : 2D pandas DataFrame array, shape = (n_solutions, n_cl)
        Dataframes giving the enrichment analysis results for the enriched terms of
        the `n_cl` best clusters of each solution in `keys`.
    cover : 2D float array, shape = (n_solutions, n_cl)
        Number of genes in enriched terms compared to total number of genes in the
        cluster for each of the `n_cl` best clusters of each solution in `keys`.
    GO : 2D string array, shape = (n_solutions, n_cl)
        GeneOntology ids of enriched terms for each of the `n_cl` best clusters of 
        each solution in `keys`.
    """
    all_enr_def = getAll_enr_def(enr_dir)
    
    enr = []
    cover = [] 
    GO = []

    prefix = all_enr_def.loc[:,'Prefix']
    for key in keys:
        try:
            dir_ = all_enr_def.loc[prefix==key,'File']

            enr_ = []
            cover_ = []
            GO_ = []
            for i in range(n_cl):
                if i>=len(bestCl_dict[key]):
                    break
                cl = bestCl_dict[key][i]
                if type(score_th)==type(None) and type(pval)==type(None):
                    enr_cli, _ = getEnrich_cl(dir_.values[0], cl, score_th, pval=pval)
                else:
                    enr_cli, enr_genes, cl_genes = getEnrich_cl(dir_.values[0], cl,
                                                                score_th, pval=pval)
                    cover_.append(len(enr_genes)/len(cl_genes))
                    if len(enr_genes[0])==0:
                        enr_genes = []
                
                GO_i = []
                for t in enr_cli.loc[:,'Term']:
                    GO_i.append(t.split(' ')[-1][1:-1])

                enr_.append(enr_cli)
                
                GO_.append(','.join(GO_i))
            enr.append(enr_)
            cover.append(cover_)
            GO.append(GO_)
        except KeyError:
            print(f'Best cluster not found for {key}.')
            enr.append([])
            cover.append([])
            GO.append([])
        
    return enr, cover, GO


def getBest_cluster_cover(enr_dir, pval=0.05):
    """Get the list of most enriched clusters according to coverage.
    
    Parameters
    ----------
    enr_dir : string
        Path to the directory containing the enrichment analysis results for all 
        clusters.
    pval : Boolean, default = False
        If True, the enrichment terms are sorted according to the adjusted p-value
        instead of the combined score.
        
    Returns
    -------
    bestCl_list : int array, shape = (n_clusters,)
        List of cluster index, from most enriched to least.
    """

    cover = []
    cl_ = []
    cl = 0
    cond = True
    while cond:
        #print(cl)
        try :
            enr_cl, enr_genes, cl_genes = getEnrich_cl(enr_dir, cl, pval=pval)
            cover.append(len(enr_genes)/len(cl_genes))
            cl_.append(cl)
            cl += 1
        except UnboundLocalError:
            cond = False
            
    df = pd.DataFrame(cover, columns=['cv'], index=cl_)
    bestScore_sort = df.sort_values('cv', ascending=False)
    items = bestScore_sort.index
    bestCl_list = [int(i) for i in items]
    #print(bestCl_list)
    return bestCl_list


def CLEAN(sol, mm_dir='./MembMatrix/', enr_dir='Enrich', pval=0.1,
          clustCLEAN=False):
    """Get gene-wise CLEAN score for clustering solution in sol.
    
    Parameters :
    ----------
    sol : string
        Prefixe of the clustering solution in `mm_dir` for which to compute CLEAN 
        score.
        
    mm_dir : string, default = './MembMatrix/'
        Path to the directory in which to find the membership matrics for the 
        solution to evaluate. The membership matrix files must be named as 
        `prefix.csv`, with prefix corresponding to `sol`.
        
    enr_dir : string, default = 'Enrich'
        Name of the subdirectory of `mm_dir` in which to find enrichment analysis 
        results of solution `sol`. 
        
    pval : float, default = 0.1
        Cutoff p-value above which the enrichment in one term is not considered as
        statisticaly significant.
        
    clustCLEAN : bool, default = False
        If `True`, the CLEAN score of clusters is also computed.
        
    
    Returns : 
    -------
    CLEAN : float array, len = n_genes
        List of gene-wise CLEAN score given in the same order as the genes in the
        membership matrix.
        
    cCLEAN : pandas DataFrame, shape = (n_genes, n_clusters)
        Dataframe of CLEAN score given for each gene in each cluster. Only returned 
        if `clustCLEAN = True`.
    """
    mm_file = os.path.join(mm_dir,sol+'.csv')
    enr_path = os.path.join(mm_dir,enr_dir,sol)
    mm = misc.load_mm(mm_file)
    mm.columns = [i for i in range(mm.shape[1])]
    all_genes = mm.index
    
    enr = getEnriched(enr_path, score_th=None, genes=all_genes, pval=pval)[0]
    
    if type(enr)==type(None):
        raise Exception(f'No enrichment analysis results found for {sol}.')
    
    func_cat = np.unique(enr.loc[:,'Term'])
    clusters = mm.columns
    p = pd.DataFrame(1, columns=func_cat, index=[int(c) for c in clusters])
    gene_cat = pd.DataFrame(0, columns=func_cat, index=all_genes)
    
    for i in range(enr.shape[0]):
        line = enr.iloc[i,:]
        cl = line['Cluster']
        term = line['Term']
        genes = line['Genes'].split(';')
        q = line['Adjusted P-value']
        p.at[cl,term]=q
        for g in genes:
            gene_cat.at[g,term]=1
            
    CLEAN = []
    for g in all_genes:
        inCl = mm.columns[mm.loc[g,:]==1]
        inCat = gene_cat.columns[gene_cat.loc[g,:]==1]
        pij = []
        for i in inCl:
            for j in inCat:
                pij.append(p.loc[i,j])
        if len(pij)==0:
            CLEAN.append(0)
        else:
            CLEAN.append(np.max(-np.log10(pij)))
            
    if not clustCLEAN:
        return CLEAN
    
    cCLEAN = mm.replace(0,None).multiply(CLEAN, axis=0)

    return CLEAN, cCLEAN

def multiCLEAN(sol_, mm_dir='./MembMatrix/', enr_dir='Enrich', pval=0.1,
               clustCLEAN=False, genes=None):
    """Get CLEAN score for all clustering solutions in sol.
    
    Parameters :
    ----------
    sol : string array, len = n_sol
        List of prefixes of the clustering solutions in `mm_dir` for which to
        compute CLEAN score.
        
    mm_dir : string, default = './MembMatrix/'
        Path to the directory in which to find the membership matrices for the 
        solutions to evaluate. The membership matrix files must be named as 
        `prefix.csv`, with prefix corresponding to what's given in list `sol`.
        
    enr_dir : string, default = 'Enrich'
        Name of the subdirectory of `mm_dir` in which to find enrichment analysis 
        results of solutions in `sol`. 
        
    pval : float, default = 0.1
        Cutoff p-value above which the enrichment in one term is not considered as
        statisticaly significant.
        
    clustCLEAN : bool, default = False
        If `True`, the CLEAN score of clusters is also computed.
        
    
    Returns : 
    -------
    CLEAN_df : pandas DataFrame, shape = (n_sol, n_genes)
        List of gene-wise CLEAN scores given in the same order as the genes in the
        membership matrix.
        
    cCLEAN : pandas DataFrame list, shape = (n_sol, n_genes, n_clusters)
        List of dataframes containing cluster-wise CLEAN scores for each solution.
        Only returned if `clustCLEAN = True`.
    """
    CLEAN_ = []
    cCLEAN_ = []
    for sol in sol_:
        if clustCLEAN:
            clean, cclean  = CLEAN(sol, mm_dir=mm_dir, enr_dir=enr_dir, pval=pval,
                                   clustCLEAN=clustCLEAN)
            
            cCLEAN_.append(cclean)
        else:
            clean  = CLEAN(sol, mm_dir=mm_dir, enr_dir=enr_dir, pval=pval,
                           clustCLEAN=clustCLEAN)
        CLEAN_.append(clean)
        
    if type(genes)!=type(None):
        CLEAN_df = pd.DataFrame(CLEAN_, index=sol_, 
                                columns=pd.Series(genes, name='Genes'))
    else:
        CLEAN_df = pd.DataFrame(CLEAN_, index=sol_)
    
    if not clustCLEAN:
        return CLEAN_df
    
    
    return CLEAN_df, cCLEAN_
        
    
def plot_multiCLEAN(CLEAN_df, n_pts=100, yscale='log', alpha=0.9):
    """Plot the number of genes with CLEAN score above x for each solution.
    
    
    Parameters :
    ----------
    CLEAN_df : pandas DataFrame, shape = (n_sol, n_genes)
        List of gene-wise CLEAN scores.
        
    n_pts : int, default = 100
        Number of points x on which to plot.
        
    yscale : {"linear", "log", "symlog", "logit", ...}, default = 'log'
        The axis scale type to apply to axis y.
        
    alpha : float in [0,1], default = 0.9
        Alpha parameter for the plot function. 
    """
    n = CLEAN_df.shape[0]
    colors = pl.cm.jet(np.linspace(0,1,n))


    Mclean = CLEAN_df.to_numpy().max()
    #print(Mclean)
    x = np.linspace(1,Mclean,n_pts)
    y = []
    for i in range(n):
        clean = np.array(CLEAN_df.iloc[i,:])
        y.append([sum(clean>=i) for i in x])
        
    for i in range(n):
        plt.plot(x,y[i], alpha=alpha, color=colors[i])
    plt.legend([i for i in CLEAN_df.index])
    plt.yscale(yscale)
    plt.xlabel('CLEAN score')
    plt.ylabel('# genes with score >= x')
    plt.show()
    

def plot_clustCLEAN(cCLEAN_df, n_pts=100, yscale='log', alpha=0.5, leg=False):
    """Plot the number of genes with CLEAN score above x for each cluster.
    
    Parameters :
    ----------
    cCLEAN_df : pandas DataFrame list, shape = (n_sol, n_genes, n_clusters)
        List of dataframes containing cluster-wise CLEAN scores for each solution.
        
    n_pts : int, default = 100
        Number of points x on which to plot.
        
    yscale : {"linear", "log", "symlog", "logit", ...}, default = 'log'
        The axis scale type to apply to axis y.
        
    alpha : float in [0,1], default = 0.9
        Alpha parameter for the plot function. 
        
    lef : bool, default = False
        Whether or not to display legend.
        
    Returns :
    -------
    x : float array, len = n_pts
        Plot points for CLEAN score (x axis).
        
    ydf : pandas DataFrame, shape = (n_clusters, n_pts)
        Number of genes with CLEAN score greater or equal to x.
    """
    n = cCLEAN_df.shape[1]
    colors = pl.cm.jet(np.linspace(0,1,n))
    
    Mclean = cCLEAN_df.fillna(0).to_numpy().max()
    x = np.linspace(1,Mclean,n_pts)
    y = []
    for i in range(n):
        cclean = cCLEAN_df.iloc[:,i]
        n_genes = sum(abs(np.isnan(cclean.to_list())-1))
        clean = np.array(cclean)
        y.append([sum(clean>i)/n_genes for i in x])

    ydf = pd.DataFrame(y, index=cCLEAN_df.columns)
    ydf = ydf.replace(0,None)
    for i in range(n):
        plt.plot(x,ydf.iloc[i,:], alpha=alpha, color=colors[i])
    if leg:
        plt.legend([i for i in cCLEAN_df.columns])
    plt.yscale(yscale)
    plt.xlabel('CLEAN score')
    plt.ylabel('# genes with score >= x')
    plt.show()
    
    return x,ydf