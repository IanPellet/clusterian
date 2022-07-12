import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#from matplotlib.patches import Ellipse
from scipy import stats
from scipy.linalg import svd
#from sklearn.manifold import MDS
from sklearn.decomposition import FastICA
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.packages import importr
import rpy2.rinterface as ri

def randDataSVD(data, K=10):
    """Create K random data sets, compute SVD and return K singular value vectors.
    
    The K random data sets are drawn from a normal distribution with same means and standard deviation as the original data set.
    
    Arguments:
    data -- original data set
    Keyword arguments:
    K -- number of random data set SVD to perform (default 10)
    Return value:
    sRand -- list of singular value vectors
    """
    dataFlat = np.array(data).flatten()
    scale = np.std(data)
    loc = np.mean(data)
                
    sRand = np.empty_like(None,shape=(K,data.shape[1]))
    for k in range(K):
        X = np.random.normal(loc, scale, data.shape)
        _,sRand[k,:],_ = svd(X)
    return sRand

def plotSVDscreeplot(s, sRand, plot=True):
    """SVD screeplot with original and mean random data, returns number of components to use.
    
    The number of components to use is defined as the maximal number of components for which the explained variance of the original data is greater than the mean explained variance of the random data sets.
    
    Arguments:
    s -- singular value vector from the SVD on original data set
    sRand -- list of K singular value vectors from the SVD on K random data sets
    Keyword arguments:
    plot -- plot the screeplot (default True)
    Return value:
    nComp -- number of components to use
    """
    sMean = np.mean(sRand, axis=0)
    expVar = s/sum(s)*100
    expVarR = sMean/sum(sMean)*100
    expVarRcor = np.zeros(shape=expVarR.shape)
    for i in range(len(expVarRcor)):
        expVarRcor[i] = expVarR[i]-expVar[i]
    if plot:
        fig, ax1 = plt.subplots()

        color = 'tab:blue'
        ax1.set_ylabel("explained variance (%)", color=color)
        ax1.plot(expVar, color=color)
        ax1.plot(expVarRcor, color='tab:purple')
        ax1.set_yscale("log")
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis


        color = 'tab:orange'
        ax2.set_ylabel("cumulative explained variance (%)", color=color)
        ax2.plot(np.cumsum(expVar), color=color)
        ax2.tick_params(axis='y', labelcolor=color)

    try:
        nComp = np.where((expVarRcor>=expVar)==True)[0][1]
        if plot:
            ax1.vlines(nComp, ls='--', ymin=ax1.get_ylim()[0], ymax=ax1.get_ylim()[1], 
                       color = '0', lw=1)
    except IndexError:
        nComp = len(expVar+1)

    if plot:
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.title("SVD screeplot")
        plt.show()
        print("Optimal number of components :", nComp)
    
    return nComp

def delIndivSpe(S_, threshold=0.10, verbose=True):
    """Remove components characterized by a single gene explaining most of the variability.
    
    Arguments:
    S_ -- components matrix
    Keyword arguments:
    threshold -- max part of explained variance for one gene from which the corresponding component is removed (default 0.10)
    verbose -- print number of components removed and remaining (default True)
    Return value:
    S_delIndiv -- components matrix after removing individual-specific components
    delComp -- indexes of the removed individual-specific components
    """
    delComp = np.where((np.amax(np.abs(S_), axis=0)>threshold)==True)[0]
    S_delIndiv = np.delete(S_, delComp, 1)
    if verbose:
        print(len(delComp), "components are characterized by a single gene who explain more than",
          threshold*100,"% of the variability of the pattern in the population.")
        print("After romving those components,", S_delIndiv.shape[1], "components remain.")
    return S_delIndiv, delComp

def plotKurt(S_delIndiv, kurt, ktsmin = 3):
    """Plot components with min, max and given kurtosis measure.
    
    Arguments:
    S_delIndiv -- components matrix
    kurt -- kurtosis measure for all components
    Keyword arguments:
    ktsmin -- a given kurtosis measure for which the component with the closest kurtosis will be ploted (default 3)
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 3))
    fig.suptitle('Distribution of components with minimal and maximal kurtosis and kurtosis='+
                 str(ktsmin), y = 1.05)

    sns.histplot((S_delIndiv[:,np.argmin(kurt)]), stat="percent", ax=axes[0])
    minKurt = f"{np.min(kurt):.1f}"
    axes[0].set_title("Kurtosis = " + minKurt)

    sns.histplot((S_delIndiv[:,np.argmin(abs(kurt-ktsmin))]), stat="percent", ax=axes[1])
    Kurt3 = f"{kurt[np.argmin(abs(kurt-ktsmin))]:.1f}"
    axes[1].set_title("Kurtosis = " + Kurt3)

    sns.histplot((S_delIndiv[:,np.argmax(kurt)]), stat="percent", ax=axes[2])
    maxKurt = f"{np.max(kurt):.1f}"
    axes[2].set_title("Kurtosis = " + maxKurt)

    plt.show()
    
def delKurtosis(S_, ktsmin=3, verbose=True):
    """Remove components characterized by a kurtosis measure above threashold ktsmin.
    
    Arguments:
    S_ -- components matrix
    Keyword arguments:
    ktsmin -- kurtosis threashold under which the component will be removed (default 3)
    verbose -- print number of components removed and remaining (default True)
    Return value:
    S_delKurt -- components matrix after removing low-kurtosis components
    delKurt -- indexes of the removed components
    """
    kurt = stats.kurtosis(S_)
    delKurt = np.where((kurt>=ktsmin)==False)[0]
    S_delKurt = np.delete(S_, delKurt, 1)
    if verbose:
        print(len(delKurt), "components are characterized by a kurtosis <", 
              ktsmin, "and will be removed.")
        print("After romving those components,", S_delKurt.shape[1], "components remain.")
    return S_delKurt, delKurt

def runFDRanalysis(S_):
    """Compute False Discovery Rate on components matrix.
    
    Arguments:
    S_ -- components matrix
    Return value:
    S_fdr -- matrix of FDR
    """
    fdrtools = importr('fdrtool')

    #sort = robjects.r['sort']
    #fdr = robjects.r['fdrtool']

    S_fdr = np.zeros_like(S_)
    for i in range(S_.shape[1]):
        #res = robjects.FloatVector(S_[:,i])
        res = FloatVector(S_[:,i])
        try:
            #S_fdr[:,i] = fdr(res, plot=False, verbose=False)[0]
            S_fdr[:,i] = fdrtools.fdrtool(res, plot=False, verbose=False)[0]
        except BaseException:
            S_fdr[:,i] = abs(S_[:,i]-1)
    return S_fdr

def plotFDR(flatFDR, FDRmax = 1e-3):
    """Plot FDR distribution and FDRmax line.
    
    Arguments:
    flatFDR -- list of FDR (flattent FDR matrix)
    Keyword arguments:
    FDRmax -- threshold before which the genes will not be accepted into the module (default 1e-3)
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 3))
    fig.suptitle('Distribution of FDR', y = 1.05)

    sns.histplot(flatFDR, stat="percent", ax=axes[0])

    sns.histplot(flatFDR, stat="percent", ax=axes[1], binrange=(0,0.02), binwidth=FDRmax/4)
    axes[1].vlines(FDRmax, ls='--', ymin=axes[1].get_ylim()[0], ymax=axes[1].get_ylim()[1], 
               color = 'r', lw=1)
    axes[1].text(x=1.5*FDRmax, y=axes[1].get_ylim()[1]/2, s= str(FDRmax), color='r')

    plt.show()
    
def moduleSummary(modules):
    """Plot modules statistics and distributions.
    
    Plots the distribution of number of genes per module and number of modules per gene.
    Prints the mean, median, max and min number of genes per module and number of modules per gene.
    Prints the part of genes non-included in any module, included in only one module and included in more than one module.
    
    Arguments:
    modules -- membership matrix
    """
    GenesPerModule = np.sum(modules, axis=0)
    ModulesPerGene = np.sum(modules, axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(15, 3))

    sns.histplot(GenesPerModule, ax=axes[0])
    axes[0].set_title("Distribution of the number of genes in a module")
    axes[0].set_xlabel("Number of genes")

    sns.histplot(ModulesPerGene, stat="percent", ax=axes[1])
    axes[1].set_title("Distribution of the number of modules one gene is a member")
    axes[1].set_xlabel("Number of modules")

    plt.show()

    print("Number of genes in a module :",
         "\n- Mean :", f"{np.mean(GenesPerModule):.2f}",
         "\n- Median :", np.median(GenesPerModule),
         "\n- Max :", np.max(GenesPerModule),
         "\n- Min :", np.min(GenesPerModule),"\n")

    print("Number of modules per gene :",
         "\n- Mean :", f"{np.mean(ModulesPerGene):.2f}",
         "\n- Median :", np.median(ModulesPerGene),
         "\n- Max :", np.max(ModulesPerGene),
         "\n- Min :", np.min(ModulesPerGene),"\n")
    
    unclust = sum(ModulesPerGene==0)/len(ModulesPerGene)*100
    print(f"{unclust:.2f}", 
          "% of genes are not included in any module.")
    print(f"{sum(ModulesPerGene==1)/len(ModulesPerGene)*100:.2f}", 
          "% of genes are included in only one module.")
    print(f"{sum(ModulesPerGene>1)/len(ModulesPerGene)*100:.2f}", 
          "% of genes are included in more than one module.")

    
def runICAmodDetection(data, IndivSpe=0.10, ktsmin=3, FDRmax = 1e-3, verbose=True, nComp=None, plot=True):
    """Run the whole ICA module detection method on a data set, return modules membership matrix.
    
    Arguments:
    data -- data set on which to run the analysis
    Keyword arguments:
    IndivSpe -- max part of explained variance for one gene from which the corresponding component is removed, if False no component is removed (default 0.10)
    ktsmin -- kurtosis threashold under which the component will be removed, if False no component is removed (default 3)
    verbose -- print number of components removed and remaining (default True)
    nComp -- number of components to use for the ICA, if None SVD is run to determine optimal number of components (default None)
    plot -- plot SVD screeplot if nComp not given (default True)
    Return value:
    modules -- modules membership matrix
    """
    # outlier removal
    #distMDS = MDS(n_components=2, dissimilarity='precomputed').fit_transform(dist)
    #outliersIndex = plotOutliersMDS(distMDS)
    #data_noOut = np.delete(data, outliersIndex, 0)
    #data_noOut = data
    
    if nComp==None:
        # SVD
        _,s,_ = svd(data)
        sRand = randDataSVD(data, K=10)
        nComp = plotSVDscreeplot(s, sRand, plot=verbose)
        del s
    
    # ICA
    ICA = FastICA(n_components=nComp, fun='logcosh')
    S_ = ICA.fit_transform(data_noOut)  # Reconstruct signals
    del ICA, data_noOut
    
    if IndivSpe!=False:
        # Individual-specific components removal
        S_delIndiv,_ = delIndivSpe(S_, threshold=IndivSpe, verbose=verbose)
        del S_
    
    if ktsmin!=False:
        # Gaussian component removal
        S_delKurt,_ = delKurtosis(S_delIndiv, ktsmin=ktsmin, verbose=verbose)
        del S_delIndiv
    
    # FDR module definition
    S_fdr = runFDRanalysis(S_delKurt)
    del S_delKurt
    modules = pd.DataFrame((S_fdr<FDRmax).astype(int), index=data.index)
    if plot:
        moduleSummary(modules)
    
    return modules