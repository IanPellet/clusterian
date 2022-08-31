import numpy as np
import random
import pandas as pd

def getMetric(res_df, res_, m):
    m_score = pd.DataFrame(None, index=res_, columns=[i for i in range(res_df[0].shape[0])], 
                           dtype=str)
    for i,df in enumerate(res_df):
        sol = res_[i]
        m_score.loc[sol,:] = df.loc[:,m]
    return m_score

# Sorts a (portion of an) array, divides it into partitions, then sorts those
def quicksort(A, lo, hi, test):
    # Ensure indices are in correct order
    if lo >= hi or lo < 0 :
        return
    
    # Partition array and get the pivot index
    p = partition(A, lo, hi, test) 
    #print(A[p])  
    # Sort the two partitions
    #print(A)
    quicksort(A, lo, p - 1, test) # Left side of pivot
    #print(A)
    quicksort(A, p + 1, hi, test) # Right side of pivot
    #print(A)

# Divides array into two partitions
def partition(A, lo, hi, test):
    pivot = A[hi] # Choose the last element as the pivot
    #print(pivot)
    # Temporary pivot index
    i = lo - 1

    for j in range(lo, hi):
        # If the current element is less than or equal to the pivot
        a = A[j]
        iaa = pd.Series(test.loc[:,'a']==a, dtype=int)
        iab = pd.Series(test.loc[:,'b']==a, dtype=int)
        ipa = pd.Series(test.loc[:,'a']==pivot, dtype=int)
        ipb = pd.Series(test.loc[:,'b']==pivot, dtype=int)

        compAB = test.loc[(iaa+ipb)==2,['P(a>b)','P(rope)','P(b>a)']].values.flatten().tolist()
        compBA = test.loc[(iab+ipa)==2,['P(a>b)','P(rope)','P(b>a)']].values.flatten().tolist()

        # a<=pivot
        if len(compAB)>len(compBA):
            cond = random.choices([0,1,1], weights=compAB, k=1)[0]
        else:
            cond = random.choices([1,1,0], weights=compBA, k=1)[0]
            
        if cond:
            #print('swap')
            # Move the temporary pivot index forward
            i = i + 1
            # Swap the current element with the element at the temporary pivot index
            tmp = A[i]
            A[i] = A[j]
            A[j] = tmp

    # Move the pivot element to the correct pivot position (between the smaller and larger elements)
    i = i + 1
    tmp = A[i]
    A[i] = A[hi]
    A[hi] = tmp
    return i # the pivot index