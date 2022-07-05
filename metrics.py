import numpy as np
import scipy.sparse as sp

def expanded_contingency_matrix(A, B, sparse=False, full=False):
    """Build an expanded contingency matrix to assess the degree of agreement 
    between two clustering.
    
    Information about the number of object pairs found in the same number of 
    clusters in the two clustering.
    
    Parameters
    ----------
    A : bool matrix, shape = (n_samples,n_classes_A)
        First membership matrix to compare.
        
    B : bool matrix, shape = (n_samples,n_classes_B)
        Second membership matrix to compare.
        
    sparse : bool, default=False
        Must be set to `True` if matrices `A` and `B` are sparse matrices.
        
    full : bool, default=False
        If `True`, returns the full contingency matrix with shape 
        (n_classes_A, n_classes_B). If `False`, only the submatix containing all the
        non-zero coefficients of the full matrix is returned.
           
    Returns
    -------
    contingency : int array, shape=[n_classes_A, n_classes_B]
        Matrix `C` such that `C_{i, j}` is the number of samples pairs that is found
        is the same `i` classes in clustering A and in the same `j` classes in 
        clustering B. Will be of shape [j,k]<=[n_classes_A, n_classes_B] if 
        `full=False`.
        
    References
    ----------
    Collins L. M., Dent C. W. (1988). Omega: A General Formulation of the Rand Index 
    of Cluster Recovery Suitable for Non-disjoint Solutions. Multivariate Behavioral 
    Research, 23(2), 231–242. https://doi.org/10.1207/s15327906mbr2302_6
    """
    n = A.shape[0]
    
    AA = np.dot(A,A.T)
    BB = np.dot(B,B.T)
    if sparse:
        AA = AA.toarray()
        BB = BB.toarray()
    
    tr_i = np.triu_indices(n=n,k=1)
    row = AA[tr_i]
    col = BB[tr_i]
    data = [1]*len(row)
    
    if full:
        return sp.coo_matrix((data, (row, col)), 
                             shape=(A.shape[1],B.shape[1])).toarray()
        
    J = np.max(row, axis=None)
    K = np.max(col, axis=None)
    return sp.coo_matrix((data, (row, col)), shape=(J+1, K+1)).toarray()

def omega_index(A,B):
    """Compute the Omega index between clustering A and B.
    
    The Omega index is a generalisation of the Adjusted Rand Index in the case of
    overlapping cluster solutions.
    
    Parameters
    ----------
    A : bool matrix, shape = (n_samples,n_classes_A)
        First membership matrix to compare.
        
    B : bool matrix, shape = (n_samples,n_classes_B)
        Second membership matrix to compare.

    Returns
    -------
    omega : float
        Similarity score between `-1.0` and `1.0`. Random labelings have an `omega`
        close to `0.0`. `1.0` stands for perfect match.

    References
    ----------
    Collins L. M., Dent C. W. (1988). Omega: A General Formulation of the Rand Index 
    of Cluster Recovery Suitable for Non-disjoint Solutions. Multivariate Behavioral 
    Research, 23(2), 231–242. https://doi.org/10.1207/s15327906mbr2302_6
    """
    n = A.shape[0]
    N = n*(n-1)/2
    
    a = sp.csr_matrix(A)
    b = sp.csr_matrix(B)

    expanded_contingency_AB = expanded_contingency_matrix(a,b,sparse=True)
    N_i = np.sum(expanded_contingency_AB, axis=1)
    N_j = np.sum(expanded_contingency_AB, axis=0)
    J = min(len(N_i),len(N_j))

    obs_i = sum(np.diag(expanded_contingency_AB))/N
    exp_i = sum(N_i[:J]*N_j[:J])/N**2

    return (obs_i-exp_i)/(1-exp_i)

def soft_omega_index(A,B):
    """Compute the Omega index between clustering A and B.
    
    The Omega index is a generalisation of the Adjusted Rand Index in the case of
    overlapping cluster solutions.
    
    Parameters
    ----------
    A : bool matrix, shape = (n_samples,n_classes_A)
        First membership matrix to compare.
        
    B : bool matrix, shape = (n_samples,n_classes_B)
        Second membership matrix to compare.

    Returns
    -------
    soft_omega : float
        Similarity score between `-1.0` and `1.0`. Random labelings have a        
        `soft_omega` close to `0.0`. `1.0` stands for perfect match.
        
    References
    ----------
    Lutov A., Khayati M., Cudré-Mauroux P. (2019, April 1). Accuracy Evaluation of
    Overlapping and Multi-Resolution Clustering Algorithms on Large Datasets. 2019
    IEEE International Conference on Big Data and Smart Computing, BigComp 2019 - 
    Proceedings. https://doi.org/10.1109/BIGCOMP.2019.8679398
    """
    n = A.shape[0]
    N = n*(n-1)/2
    
    a = sp.csr_matrix(A)
    b = sp.csr_matrix(B)

    expanded_contingency_AB = expanded_contingency_matrix(a,b,sparse=True,full=True)
    del a,b
    N_i = np.sum(expanded_contingency_AB, axis=1)
    N_j = np.sum(expanded_contingency_AB, axis=0)
    
    szmax = max(expanded_contingency_AB.shape)
    tui = np.triu_indices(szmax)
    tli = np.tril_indices(szmax)
    nom = np.ones((szmax,szmax))
    nom[tui] = tui[0]
    nom[tli] = tli[1]
    nom[0,0] = 1
    nom = nom[:expanded_contingency_AB.shape[0],:expanded_contingency_AB.shape[1]]
    denom = np.ones((szmax,szmax))
    denom[tui] = tui[1]
    denom[tli] = tli[0]
    denom[0,0] = 1
    denom = denom[:expanded_contingency_AB.shape[0],
                  :expanded_contingency_AB.shape[1]]
    del tui,tli
    
    soft_fact = nom/denom
    del nom,denom
    nobs = np.sum(expanded_contingency_AB*soft_fact)/N
    del soft_fact
    
    J = min(len(N_i),len(N_j))
    exp_i = sum(N_i[:J]*N_j[:J])/N**2
    return (nobs - exp_i)/(1-exp_i)