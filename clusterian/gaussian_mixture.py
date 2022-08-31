import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
import pickle

def choose_cv_type(data, n_components_range=[150, 200, 250, 300], 
                   cv_types = ["spherical", "tied", "diag", "full"], 
                   save_gmms=False, file_name=None):
    """Evaluate the score for different parameters for Gaussian Mixture clustering
    
    Parameters
    ----------
    data : float matrix, shape = (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds to a single 
        data point.
        
    n_components_range : int array, default=[150, 200, 250, 300]
        Different values for `n_components` to test for Gaussian Mixture computation.
        
    cv_types : {‘full’, ‘tied’, ‘diag’, ‘spherical’} array, 
               default=["spherical", "tied", "diag", "full"]
        Strings describing the type of covariance parameters to test.
        Possible values :
        - 'full': each component has its own general covariance matrix.
        - 'tied': all components share the same general covariance matrix.
        - 'diag': each component has its own diagonal covariance matrix.
        - 'spherical': each component has its own single variance.
        
    save_gmms : bool, default=False
        If `True`, all computed gaussian mixtures, tested covariance types, tested
        number of components and scores are saved to `file_name`.
        
    file_name : string, default=None
        Path and file name to save the results.
        Must be specified if `save_gmms=True`.
        
    
    Returns
    -------
    scores_df : pandas DataFrame, shape=(n_tested_cv_types, n_tested_components)
        Likelyhood scores computed for each computed gaussian mixture stored i a
        DataFrame with the covarience type used as index and number of gaussian
        components as column.
        
    gmms : 2D array, shape=(n_tested_cv_types*n_tested_components, 3)
        Each row contains the covariance type, number of components and the 
        corresponding fitted gaussian mixture object.
    """
    n_components_range.reverse()
    gmms = []
    for cv_type in cv_types:
        print('\n',cv_type)
        for n_components in n_components_range:
            print(n_components, end = ' ')
            # Fit a Gaussian mixture with EM
            gmm = GaussianMixture(n_components=n_components, covariance_type=cv_type,
                                  init_params= "kmeans", n_init=10)
            gmm.fit(data)
            gmms.append([cv_type,n_components,gmm])
            print("ok", end = ' / ')
            
    scores = []
    for cv,n,gmm in gmms:
        s = gmm.score(data)
        scores.append([cv,n,s])
        
    if save_gmms:
        if file_name==None:
            raise("To save GMMs a file name needs to be specified")
        if file_name[-4:] != '.txt':
            file_name += '.txt'
        file = open(file_name, "wb")
        pickle.dump([gmms, cv_types, n_components_range, scores], file)
        file.close()
        
    scores_df = pd.DataFrame(0, index=cv_types, columns=n_components_range)
    for cv,n,s in scores:
        scores_df.loc[cv,n] = s
        
    return scores_df, gmms