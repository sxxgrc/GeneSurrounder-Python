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
TODO: Shouldn't this use the putative genes?
"""
def computeD(thresh, expr, grade, dist):
    g_vals = computeG(expr, grade)
    discordance = pd.Series()
    genes = dist.keys()
    
    for gene in genes:
        nearby = [(val, dist[gene][val]) for val in genes if dist[gene][val] < thresh]
        d = []
        g = []
        
        for node in nearby:
            if node[0] in g_vals.index:
                d.append(node[1])
                g.append(g_vals[node[0]])

        tau = kendall_corr(g, d)[0]
        discordance[gene] = tau
    
    return discordance

"""
Computes the distribution of DDE across the given dataset and network.
"""
def computeDistD(thresh, expr, grade, dist):
    # Amount of permutations to do.
    perm_count = 10 ** 3
    
    # Generate distribution across all permutations.
    discordance = pd.DataFrame(columns=range(perm_count))
    grade_var = pd.DataFrame(index=grade.index)
    
    for k in range(perm_count):
        l = [rand.random() for i in range(len(grade))]
        t = sum(grade['x'] == 'high') / len(grade)
        perm_grade = []
        
        for i in l:
            if i < t:
                perm_grade.append('high')
            else:
                perm_grade.append('low')
        
        grade_var['x'] = perm_grade
        sample = computeD(thresh, expr, grade_var, dist)
        discordance.iloc[:, k] = sample

    return discordance


# Main Driver

"""
Main function for the Decay Differential Expression algorithm.

Parameters:
    - expr : DataFrame for the expression values over the genes in the data
    - grade : DataFrame for the grade values over the experiments in the data
    - thresh : Threshold value for the algorithm
    - dist : Distance matrix of all geodesic distances between node pairs in G
"""
def decayDE(expr, grade, thresh, dist):
    trueval = computeD(thresh, expr, grade, dist)
    distribution = computeDistD(thresh, expr, grade, dist)
    pvals = distribution.subtract(trueval, 0).where(distribution < 0).count(1)
    pvals /= distribution.shape[1]
    
    return pvals
