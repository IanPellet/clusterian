import gseapy as gp
import pandas as pd
import numpy as np
import os

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
            enr = gp.enrichr(gene_list=genes,
                             gene_sets=gene_sets,
                             organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                             description=desc,
                             outdir=os.path.join(od, desc+'_'+str(i)),
                             no_plot=True,
                             cutoff = 1
                            )
            mm_enrich.append(enr)
            #minAdjPval.append(enr.results.loc[0,"Adjusted P-value"])
        else:
            mm_enrich.append(None)
            #minAdjPval.append(1)
        print("\r",round((i+1)/modules.shape[1]*100),"%",end="\r")
        
    #return minAdjPval
    return od

def getEnriched(enr_dir, score_th = 300):
    """Create and return a DataFrame of all the enriched terms in each cluster.
    
    Arguments:
    enr_dir -- path to the directory where to find the subdirectories created by enrichr analysis for each cluster
    Keyword arguments:
    score_th -- Combined Score threshold above which a term is considered as enriched (default 300)
    Return value:
    enriched -- pandas DataFrame with all the enriched terms in each cluster, first column is the coresponding cluster number
    """
    enriched_ = []
    for root, dirs, files in os.walk(enr_dir, topdown=False):
        for name in files:
            cl = int(root.split('_')[-1][3:])
            path = os.path.join(root, name)
            full_tab = pd.read_csv(path, sep='\t').sort_values('Adjusted P-value', ascending=True)
            enr_tab = full_tab.loc[full_tab.loc[:,'Combined Score']>score_th,:]
            new_col = pd.Series([cl]*enr_tab.shape[0], index=enr_tab.index, dtype=int)
            enriched_.append(pd.concat([new_col,enr_tab], axis=1))

    enriched = pd.concat(enriched_, axis=0).sort_values(0)
    col = enriched.columns.tolist()
    col[0] = 'Cluster'
    enriched.columns = col
    return enriched

