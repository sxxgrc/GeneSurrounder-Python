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
    
    rho_sum = 0
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

    rho_sum = [0 for _ in range(resample)]
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
        p_SI[i] = round(len(overlap) / len(SI_permute_r), 8)
        if p_SI[i] == 0:
            p_SI[i] = 1/(resample+1)
    
    return p_SI
