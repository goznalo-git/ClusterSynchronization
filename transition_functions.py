"""
Author: Gonzalo Contreras Aso
Title: Functions characterizing the transition to synchronization
"""

import numpy as np
from itertools import product

def E_lambda(V, i):
    """Given the eigenvector matrix V (with them as columns),
    compute the matrix E_lambda_i of the differences between componets.
    """
    E = np.zeros(V.shape)
    for j, k in product(range(V.shape[0]), range(V.shape[1])):
        E[j,k] = (V[j,i] - V[k,i])**2
    return E
   
def Sn_matrices(V):
    """Given an eigenvector matrix V (with them as columns), compute 
    the Sn matrices whose entry i,j is 2 if those nodes for a cluster.
    """
    Si = {}
    Sn = np.zeros(eigenVectors.shape)
    for i in range(eigenVectors.shape[1])[::-1]:
        Sn += E_lambda(eigenVectors, i)
        Si[i] = np.copy(Sn)

    return Si 
