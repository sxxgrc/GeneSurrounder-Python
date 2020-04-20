"""
An implementation of the Decay of Differential Expression (DDE) component of the
GeneSurrounder algorithm.

Main driver function is at the bottom, describing more about the program.
"""

import networkx as nx
import pandas as pd
import random as rand
from scipy.stats import ttest_ind, kendalltau as kendall_corr
from scipy.stats import combine_pvalues
import pandas as pd
import subprocess
import os


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
    pvals = {}
    for thresh in range(1, diameter + 1):
        print(thresh, 'of', diameter)
        nodes_in_rad = []
        for x, val in dist.items():
            if val < thresh:
                nodes_in_rad.append(x)

        if len(nodes_in_rad) == 0:
            pvals[thresh] = 0.5
        else:
            dist_filter = [dist[x] for x in nodes_in_rad]
            expr_radius = expr.loc[nodes_in_rad]
            g_vals = computeG(expr_radius, grade)

            trueval = computeD(g_vals, dist_filter)
            distribution = computeDistD(g_vals, dist_filter)
            pvals[thresh] = sum([x < trueval for x in distribution]) / len(distribution)

    return pvals






"""
Main function for the Decay Differential Expression algorithm.

Parameters:
    - expr : DataFrame for the expression values over the genes in the data
    - grade : DataFrame for the grade values over the experiments in the data
    - dist : Distance matrix of all geodesic distances between assayed node pairs in G
    - diameter : Diameter of the network


def decayDE(expr, grade, dist, diameter):

"""

"""
An implementation of the Sphere of Influence component of the GeneSurrounder
algorithm.

Main driver function is at the bottom, explaining more about the computation.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scipy.stats import spearmanr

# Helper Functions

"""
Builds the SI component which determines if a gene is more strongly correlated with its
neighbors.
"""


def SIObserve(distance, correlation, diameter, overlap_genes):
    SI = [0 for _ in range(diameter)]

    rho_sum = 0.0
    for i in range(diameter):
        index = [j for j in range(len(overlap_genes)) if distance[overlap_genes[j]] == i + 1]
        rho_sum += sum(abs(correlation[z]) for z in index)
        SI[i] = rho_sum

    return SI


"""
Generates resamples and calculates resulting SI score for the gene.
"""


def SIPermute(overlap_genes, distance, correlation, diameter, resample, gene_of_interest):
    SI_permute = np.zeros((resample, diameter))
    neighbors = []

    rho_sum = [0.0 for _ in range(resample)]
    for i in range(diameter):
        for k in range(resample):
            index = [j for j in range(len(overlap_genes)) if distance[overlap_genes[j]] == i + 1]
            a = [v for v in range(len(overlap_genes))]
            a.remove(overlap_genes.index(gene_of_interest))
            permute_index = np.random.choice(a, size=len(index), replace=True, p=None)
            rho_sum[k] += sum([abs(correlation[z]) for z in permute_index])
            SI_permute[k, i] = rho_sum[k]

    return SI_permute


# Main Driver

"""
Main function for the Sphere of Influence computation.

Parameters:
    - dist : Geodesic distance matrix for the assayed genes general network
    - expr : Expression data for genes in experiments
    - diameter : The diameter of the network
    - overlap_genes : Overlapping genes between network and expression data
    - resample : Number of resamples desired
    - gene_of_interest : The specific gene to run this computation on
"""


def sphereOfInf(dist, expr, diameter, overlap_genes, resample, gene_of_interest):
    # Compute correlation data frame.
    cor_index = list(expr.index)
    cor = pd.DataFrame(data=spearmanr(expr, axis=1)[0], index=cor_index, columns=cor_index)

    # Compute SI scores.
    correlation = [cor[gene][gene_of_interest] for gene in overlap_genes]
    SI_matrix = SIObserve(dist, correlation, diameter, overlap_genes)
    SI_permute = SIPermute(overlap_genes, dist, correlation, diameter, resample, gene_of_interest)

    p_SI = [0 for _ in range(len(SI_matrix))]
    for i in range(len(SI_matrix)):
        SI_permute_r = list(SI_permute[:, i])
        overlap = [j for j in range(len(SI_permute_r)) if SI_permute_r[j] >= SI_matrix[i]]
        p_SI[i] = (len(overlap) * 1.0) / (len(SI_permute_r) * 1.0)

        if p_SI[i] == 0:
            p_SI[i] = 1.0 / (resample + 1)

    return p_SI


def import_cur_data(filename = 'CurOvGradeKEGGnets.RData', folder = 'data'):
    path = os.getcwd()
    command = 'Rscript'
    path2script = path + '/read_cur_data.R'
    filename = folder + '/' + filename
    args = [filename]
    cmd = [command, path2script] + args
    subprocess.check_output(cmd, universal_newlines=True)

    files = []
    # r=root, d=directories, f = files
    exprs = {}
    grades = {}
    for r, d, f in os.walk(path):
        for file in f:
            if '.csv' in file:
                if 'expr' in file:
                    loc = os.path.join(r, file)
                    df = pd.read_csv(loc, index_col = 0)
                    exprs[file[4]] = df
                if 'grade' in file:
                    loc = os.path.join(r, file)
                    df = pd.read_csv(loc, index_col = 0)
                    grades[file[5]] = df

    expr = []
    grade = []
    for i in range(len(exprs)):
        expr.append(exprs[str(i)])
        grade.append(grades[str(i)])

    return expr, grade

[expr, grade] = import_cur_data('CurOvGradeKEGGnets.RData')
datasets = [2, 5, 7]
print('line 246')


def get_items_list(obj, ind):
    return [obj[x] for x in ind]

grade = get_items_list(grade, datasets)
expr = get_items_list(expr, datasets)

nodes = pd.read_csv('data_conversion/Sample3.txt')
edgelist = pd.read_csv('data_conversion/largestCompKEGGigraph.txt', sep=' ', header = None)
edgelist = edgelist.replace(nodes.index, nodes.values)
# edgelist = edgelist[0:1000]
G = nx.from_pandas_edgelist(edgelist, 0, 1)



# decay = decayDE('hsa:1365', expr[0], grade[0], 4, G)
#     # expr, grade, thresh, G

gene_ID = 'hsa:1365'
print('line 267')
# # to run through all datasets, iterate over expr and grade as below
DE = []
SI = []
comb_vals = []
def dist_vect(gene_ID, desired_nodes, net):
    dist = {}
    for node in desired_nodes:
        dist[node] = nx.shortest_path_length(net, gene_ID, node)

    return dist


count = 0
for expr_sample, grade_sample in zip(expr, grade):
    count += 1
    print('dataset', count)
    overlap_genes = set(expr_sample.index).intersection(set(G.nodes))
    dist = dist_vect(gene_ID, overlap_genes, G)
    diameter = max(dist.values())
    if gene_ID not in overlap_genes:
        print("Error, gene not in overlap")
    expr_filter = expr_sample.loc[overlap_genes]
    DE.append(decayDE(expr_filter, grade_sample, dist, diameter))
    SI.append(sphereOfInf(dist, expr_sample, diameter, list(overlap_genes), 1000, gene_ID))
    comb = {}
    for thresh in DE[-1].keys():
        print(thresh)
        comb[thresh] = combine_pvalues([DE[-1][thresh], SI[-1][thresh - 1]])
    comb_vals.append(comb)
print(DE, SI)
print(comb_vals)

