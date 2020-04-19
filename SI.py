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
def SIObserve(distance, correlation, diameter, overlap_genes, gene_of_interest):
    SI = np.zeros(1, diameter)
    neigbors = []
    
    rho_sum = 0
    for i in range(diameter):
        index = [j for j in range(len(distance)) if distance[j] == i + 1]
        rho_sum += sum(abs(correlation[z]) for z in index)
        SI[0, i] = rho_sum
        
    return SI

"""
Generates resamples and calculates resulting SI score for the gene.
"""
def SIPermute(overlap_genes, distance, correlation, diameter, resample, gene_of_interest):
    SI_permute = np.zeros(resample, diameter)
    neighbors = []

    rho_sum = np.zeros(1, resample)
    for i in range(diameter):
        for k in range(resample):
            index = [j for j in range(len(distance)) if distance[j] == i + 1]
            a = [v for v in range(len(overlap_genes))]
            a.remove(overlap_genes.index(gene_of_interest))
            permute_index = np.random.choice(a, size=len(index), replace=True, p=None)
            rho_sum[0, k] += sum(abs(correlation[z]) for z in permute_index)
            SI_permute[k, i] = rho_sum[0, k]
            
    return SI_permute


# Main Driver

"""
Main function for the Sphere of Influence computation.

Parameters:
    - dist : Geodesic distance matrix for the general network
    - expr : Expression data for genes in experiments
    - diameter : The diameter of the network
    - overlap_genes : Overlapping genes between network and expression data
    - resample : Number of resamples desired    
"""
def sphereOfInf(dist, expr, diameter, overlap_genes, resample):
    # Compute correlation data frame.
    cor_index = list(expr.index)
    cor = pd.DataFrame(data=spearmanr(expr, axis=1)[0], index=cor_index, columns=cor_index)
    
    # Compute SI scores. TODO: How to integrate gene of interest so as to get overall values?
    gene_of_interest = None
    distance = [dist[gene_of_interest][gene] for gene in overlap_genes]
    correlation = [cor[gene_of_interest][gene] for gene in overlap_genes]
    SI_matrix = SIObserve(distance, correlation, diameter, overlap_genes, gene_of_interest)
    SI_permute = SIPermute(overlap_genes, distance, correlation, diameter, resample, 
                           gene_of_interest)
    
    p_SI = [0 for _ in range(len(SI_matrix[0, :]))]
    for i in range(len(SI_matrix[0, :])):
        SI_permute_r = list(SI_permute[:, i])
        overlap = [j for j in range(len(SI_permute_r)) if SI_permute_r[j] >= SI_matrix[:, i]]
        p_SI[i] = len(overlap) / len(SI_permute_r)
    
    return p_SI
