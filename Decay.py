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



def dist_vect(gene_ID, desired_nodes, net):
    dist = {}
    for node in desired_nodes:
        dist[node] = [nx.shortest_path_length(net, gene_ID, node)]
    dist_df = pd.DataFrame(dist.values(), index = dist.keys())
    return dist_df

# dist = dist_mat(G)

def computeD(thresh, expr, grade, dist):
    dist_filter = dist[dist[0] < thresh]
    expr_radius = expr.loc[dist_filter.index]
    g_vals = compute_G(expr_radius, grade)  
    tau = kendall_corr(g_vals, dist_filter[0])[0]
    return tau
# dist = dist_mat(G)
# out = computeD(4, expr[0][0:1000], grade[0][0:1000], dist)
# out.head

def computeDdist(thresh, expr, dist, grade, G):
    perm_count = 10**3
    # should be 10**3 in official
    
    dist_filter = dist[dist[0] < thresh]
    expr_radius = expr.loc[dist_filter.index]
    g_vals = compute_G(expr_radius, grade)
    distribution = []
    for k in range(perm_count):
        for j in range(50):
            swap1 = rand.choice(g_vals.index)
            swap2 = rand.choice(g_vals.index)
        
            temp = g_vals[swap1]
            g_vals[swap1] = g_vals[swap2]
            g_vals[swap2] = temp
            
        tau = kendall_corr(g_vals, dist_filter[0])[0]
        distribution.append(tau)

    return distribution
        
        
    


def decayDE(gene_ID, expr, grade, thresh, G):
    """
    expr is dataframe of expression of each gene (row) in each experiment (column)
    """
    
    overlap = set(expr.index).intersection(set(G.nodes))
    if gene_ID not in overlap:
        return "Error: Target gene not in overlap"


    dist = dist_vect(gene_ID, overlap, G)

    trueval = computeD(thresh, expr, grade, dist)
    distribution = computeDdist(thresh, expr, dist, grade, G)
    pval = sum([x < trueval for x in distribution])/len(distribution)
    return pval
   
#decay = decayDE('hsa:1365', expr[0], grade[0], 4, G)
    # expr, grade, thresh, G





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



# decay = decayDE('hsa:1365', expr[0], grade[0], 4, G)
#     # expr, grade, thresh, G


# # to run through all datasets, iterate over expr and grade as below
decay = []
for expr_sample, grade_sample in zip(expr, grade):
    decay.append({})
    for R in range(1, 10):
        decay[-1][R] = decayDE('hsa:1365', expr_sample, grade_sample, R, G)