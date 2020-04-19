#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import networkx as nx
from import_cur_data import import_cur_data
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import kendalltau as kendall_corr
import random as rand

def ttest(expr):
    high_val = expr.loc['high'].values
    low_val = expr.loc['low'].values
    return ttest_ind(high_val, low_val)[1]

def compute_G(expr, grade):
    exp = expr.transpose()
    exp['grade'] = grade
    exp = exp.set_index('grade')
    exp = exp.apply(ttest)
    
    return exp

def dist_mat(net):
    dist = nx.floyd_warshall_numpy(net)
    dist_df = pd.DataFrame(dist, columns = net.nodes, index = net.nodes)
    return dist_df

# dist = dist_mat(G)

def computeD(thresh, expr, grade, dist):
    g_vals = compute_G(expr, grade)
    discordance = pd.Series()
    for gene in dist.index:
        nearby = dist[gene][dist[gene] < thresh]
        d = []
        g = []
        for node in nearby.index:
            if node in g_vals.index:
                d.append(nearby[node])
                g.append(g_vals[node])
            
        tau = kendall_corr(g, d)[0]
        discordance[gene] = tau
    return discordance
# dist = dist_mat(G)
# out = computeD(4, expr[0][0:1000], grade[0][0:1000], dist)
# out.head

def computeDdist(thresh, expr, grade, G):
    perm_count = 10**3
    # should be 10**3 in official
    
    discordance = pd.DataFrame(columns = range(perm_count))
    dist = dist_mat(G)
    grade_var = pd.DataFrame(index = grade.index)
    for k in range(perm_count):
#         l = [rand.randint(0, 1) for i in range(len(grade))]
#         perm_grade = []
#         for i in l:
#             if i == 1:
#                 perm_grade.append('high')
#             else:
#                 perm_grade.append('low')
                
        l = [rand.random() for i in range(len(grade))]
        perm_grade = []
        t = sum(grade['x'] == 'high')/len(grade)
        for i in l:
            if i < t:
                perm_grade.append('high')
            else:
                perm_grade.append('low')
                
        
        grade_var['x'] = perm_grade
        sample = computeD(thresh, expr, grade_var, dist)
        discordance.iloc[:, k] = sample
    return discordance
        
def decayDE(expr, grade, thresh, G):    
    dist = dist_mat(G)
    trueval = computeD(thresh, expr, grade, dist)
    distribution = computeDdist(thresh, expr, grade, G)
    pvals = distribution.subtract(trueval, 0).where(distribution < 0).count(1)
    pvals = pvals/distribution.shape[1]
    return pvals

### driver
[expr, grade] = import_cur_data('data_auto_import/CurOvGradeKEGGnets.RData')
datasets = [2, 5, 7]

def get_items_list(obj, ind):
    return [obj[x] for x in ind]

grade = get_items_list(grade, datasets)
expr = get_items_list(expr, datasets)

nodes = pd.read_csv('data/Sample3.txt')
edgelist = pd.read_csv('data/largestCompKEGGigraph.txt', sep=' ', header = None)
edgelist = edgelist.replace(nodes.index, nodes.values)
# edgelist = edgelist[0:1000]
G = nx.from_pandas_edgelist(edgelist, 0, 1)

decay = decayDE(expr[0], grade[0], 4, G)
    # expr, grade, thresh, G
decay.head()
