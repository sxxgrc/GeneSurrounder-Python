"""
An implementation of the Decay of Differential Expression (DDE) component of the 
GeneSurrounder algorithm.

Main driver function is at the bottom, describing more about the program.
"""

import networkx as nx
import pandas as pd
import random as rand
from scipy.stats import ttest_ind, kendalltau as kendall_corr

# Helper Functions

def ttest(expr):
    high_val = expr.loc['high'].values
    low_val = expr.loc['low'].values
    
    return ttest_ind(high_val, low_val)[1]

def computeG(expr, grade):
    exp = expr.transpose()
    exp['grade'] = grade
    exp = exp.set_index('grade')
    exp = exp.apply(ttest)
    
    return exp

"""
Computes the decay of differential expression value for the graph.
"""
def computeD(g_vals, dists):
    return kendall_corr(g_vals, dists)[0]

"""
Computes the distribution of DDE across the given dataset and network.
"""
def computeDistD(g_vals, dists):
    # Amount of permutations to do.
    perm_count = 10 ** 3
    g_copy = g_vals.copy()
    
    # Generate distribution across all permutations.
    distribution = []
    
    for k in range(perm_count):        
        for j in range(50):
            swap1 = rand.choice(g_copy.index)
            swap2 = rand.choice(g_copy.index)
            
            temp = g_copy[swap1]
            g_copy[swap1] = g_copy[swap2]
            g_copy[swap2] = temp
        
        tau = kendall_corr(g_copy, dists)[0]
        distribution.append(tau)

    return distribution


# Main Driver

"""
Main function for the Decay Differential Expression algorithm.

Parameters:
    - expr : DataFrame for the expression values over the genes in the data
    - grade : DataFrame for the grade values over the experiments in the data
    - dist : Distance matrix of all geodesic distances between assayed node pairs in G
    - diameter : Diameter of the network
"""
def decayDE(expr, grade, dist, diameter):
    # Computation of p value.    
    pvals = [0 for _ in range(diameter)]
    for thresh in range(diameter):
        dist_filter = list(zip(*[(val, dist[val]) for val in dist.keys() if 
                                 dist[val] < thresh + 1]))
        expr_radius = expr.loc[list(dist_filter[0])]
        g_vals = computeG(expr_radius, grade)

        trueval = computeD(g_vals, dist_filter[1])
        distribution = computeDistD(g_vals, dist_filter[1])
        val = float(len([x <= trueval for x in distribution])) / len(distribution)
        pvals[thresh] = val if val > 0 else (1.0 / ((10 ** 3) + 1))

    return pvals
