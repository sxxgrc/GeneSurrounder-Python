import pandas as pd
edgelist=pd.read_csv('largestCompKEGGigraph.txt',sep=' ')
ovarian_expression1=pd.read_csv('Ovarian_expression_1.txt', sep='\t')
ovarian_expression2=pd.read_csv('Ovarian_expression_2.txt', sep='\t')
ovarian_expression3=pd.read_csv('Ovarian_expression_3.txt', sep='\t')
ovarian_expression4=pd.read_csv('Ovarian_expression_4.txt', sep='\t')
ovarian_expression5=pd.read_csv('Ovarian_expression_5.txt', sep='\t')
ovarian_expression6=pd.read_csv('Ovarian_expression_6.txt', sep='\t')
ovarian_expression7=pd.read_csv('Ovarian_expression_7.txt', sep='\t')
ovarian_expression8=pd.read_csv('Ovarian_expression_8.txt', sep='\t')
ovarian_expression9=pd.read_csv('Ovarian_expression_9.txt', sep='\t')
ovarian_expression10=pd.read_csv('Ovarian_expression_10.txt', sep='\t')

nodes=pd.read_csv('Sample3.txt')
Kegg_gene=list(nodes.iloc[:,0])

# Calculate the expression data from the 3 datasets
import scipy as sp
import scipy.stats
ova_genes=list(ovarian_expression1.index)

def calculate_spearman_corr(expression):
    rho, pval = sp.stats.spearmanr(expression,axis=1)
    return(rho)

#construct the matrix of spearman correlation of the genes from each experiment
# NOTE: authors used dataset 2,3 and 7 (GSE14764 GSE17260 GSE9891)
rho1=calculate_spearman_corr(ovarian_expression1)
rho2=calculate_spearman_corr(ovarian_expression2)
rho3=calculate_spearman_corr(ovarian_expression3)
rho4=calculate_spearman_corr(ovarian_expression4)
rho5=calculate_spearman_corr(ovarian_expression5)
rho6=calculate_spearman_corr(ovarian_expression6)
rho7=calculate_spearman_corr(ovarian_expression7)
rho8=calculate_spearman_corr(ovarian_expression8)
rho9=calculate_spearman_corr(ovarian_expression9)
rho10=calculate_spearman_corr(ovarian_expression10)

# Build a dataframe using one of the ovarian dataset GSE14764
rho2_matrix=pd.DataFrame(data=rho2,index=list(ovarian_expression2.index),columns=list(ovarian_expression2.index))

# Find overlapping gene between KEGG and the ovarian data
overlap_genes=list(set(ova_genes).intersection(set(Kegg_gene)))

# Rebuild the KEGG network
import networkx as nx
G=nx.Graph()
for i in range(0,len(edgelist)):
    index1=edgelist.iloc[i,0]
    index2=edgelist.iloc[i,1]
    G.add_edge(nodes.iloc[index1,0],nodes.iloc[index2,0])

import numpy as np
# Calculate the distance matrix (geodesic distance between each of the two nodes)
# In the distance matrix, row=source, column=target
length=dict(nx.all_pairs_shortest_path_length(G))
distance=np.ones((len(list(G.nodes())),len(list(G.nodes()))))
for i in range(0,len(list(G.nodes()))):
    for j in range(i,len(list(G.nodes()))):
        distance[i,j]=length[list(G.nodes())[i]][list(G.nodes())[j]]
        distance[j,i]=length[list(G.nodes())[j]][list(G.nodes())[i]]

# Build the distance matrix into a dataframe
distance_matrix=pd.DataFrame(data=distance, index=list(G.nodes()),columns=list(G.nodes()))

# Calculate diameter of Kegg network
diameter=nx.diameter(G)

#build the SI component
#the SI component determines if a gene is more strongly correlated with its neighbors
def SIobserve(distance_matrix,cor_matrix,diameter,genes_assayed_network,gene_of_interest):
    SI=np.zeros((1,diameter))
    neighbors=[]
    distance=distance_matrix[gene_of_interest][genes_assayed_network]
    correlation=cor_matrix[gene_of_interest][genes_assayed_network]
    Kegg_gene_name=list(distance_matrix.columns)
    #print(Kegg_gene_name)
    rho_sum=0 
    for i in range(0,diameter):
        index = [j for j in range(len(list(distance))) if list(distance)[j] == i+1] 
        rho_sum=rho_sum+sum(abs(list(correlation)[z]) for z in index)
        SI[0,i]=rho_sum                   
    return(SI)
        
# Generate resamples and calculate the SI score for a given gene of interest
def SIpermutation(genes_assayed_network,distance_matrix,cor_matrix,diameter,resample,gene_of_interest):
    SIpermute=np.zeros((resample,diameter))
    neighbors=[]
    distance=distance_matrix[gene_of_interest][genes_assayed_network]
    correlation=cor_matrix[gene_of_interest][genes_assayed_network]
    Kegg_gene_name=list(distance_matrix.columns)
    #print(Kegg_gene_name)
    rho_sum=np.zeros((1,resample))
    for i in range(0,diameter):
        for k in range(0,resample):
            index = [j for j in range(len(list(distance))) if list(distance)[j] == i+1] 
            a=list(range(0,len(genes_assayed_network)))
            a.remove(genes_assayed_network.index(gene_of_interest))
            permute_index=np.random.choice(a, size=len(index), replace=True, p=None)
            rho_sum[0,k]=rho_sum[0,k]+sum(abs(list(correlation)[z]) for z in permute_index)
            SIpermute[k,i]=rho_sum[0,k]  
    return(SIpermute)
            
# Call SIobserve and SIpermutation
# Input: distance matrix, correlation matrix, overlapping genes, gene of interest, number of resample needed
# Output: SI matrix (SI score at each distance) and p_SI (p-value of gene's SI score based on null distribution)
def SI(distance_matrix,cor_matrix,dimeter,genes_assayed_network,gene_of_interest,resample):
    SI_matrix=SIobserve(distance_matrix,cor_matrix,diameter,genes_assayed_network,gene_of_interest)
    SIpermute=SIpermutation(overlap_genes,distance_matrix,rho2_matrix,diameter,resample,gene_of_interest)
    p_SI=np.zeros((1,len(SI_matrix[0,:])))
    for i in range(0,len(SI_matrix[0,:])):
        SIpermute_r=list(SIpermute[:,i])
        overlap=[j for j in range(0,len(SIpermute_r)) if SIpermute_r[j]>=SI_matrix[:,i]]
        p_SI[0,i]=len(overlap)/len(SIpermute_r)
    return(p_SI,SI_matrix)


# # Below are tests, corresponding to supplement 12859_2019_2829_MOESM1_ESM
hsa2720_p_SI=SI(distance_matrix,rho2_matrix,diameter,overlap_genes,'hsa:2720',200)
print(hsa2720_p_SI)

hsa2562_p_SI=SI(distance_matrix,rho2_matrix,diameter,overlap_genes,'hsa:2562',200)
print(hsa2562_p_SI)

hsa2316_p_SI,hsa2316_Si_matrix=SI(distance_matrix,rho2_matrix,diameter,overlap_genes,'hsa:2316',200)

print(hsa2316_p_SI)
print(hsa2316_Si_matrix)

p_SI=np.zeros((1,len(SI_matrix[0,:])))
for i in range(0,len(SI_matrix[0,:])):
    SIpermute_r=list(SIpermute[:,i])
    overlap=[j for j in range(0,len(SIpermute_r)) if SIpermute_r[j]>=SI_matrix[:,i]]
    p_SI[0,i]=len(overlap)/len(SIpermute_r)
print(p_SI)

