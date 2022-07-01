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
                print("Cluster", i, ":", genes)
            #else :
            #    print("Cluster", i, ":", genes)
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

def runEnrichr_directory(mm_path = './ModuleMatrix', mm_enrich_dir = 'Enrich'):
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

        desc = filename[:-4]

        if desc not in enr_dir_content:
            print('Runing analysis on ', desc) 

            mm = pd.read_csv(path, header=None, index_col=0)
            if type(mm.index[0])!=str or len(mm.index[0])<2:
                mm = pd.read_csv(path, header=0, index_col=0)
            mm.columns = range(0,mm.shape[1])

            print('out_dir :', enr_path)
            where_ = enrich.runEnrichr(mm, desc, gene_sets='GO_Biological_Process_2021', 
                                       out_dir= enr_path)
            print('Enrichment results in :', where_)