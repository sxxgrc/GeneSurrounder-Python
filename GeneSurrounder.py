"""
Main program for the GeneSurrounder algorithm.
"""

import argparse
import networkx as nx
import pandas as pd
from scipy.stats import combine_pvalues

from data_conversion.import_cur_data import import_cur_data
from data_conversion.store_dict import save_obj, load_obj
from DecayDE import decayDE
from SI import sphereOfInf

# Take in command line arguments.
parser = argparse.ArgumentParser()
parser.add_argument("-v", action="store_true", help="Verbosity of the program.")
parser.add_argument("-netfile", nargs="?", default="data/KEGGNetwork.txt", 
                    help="Path to the network model of cellular interactions. Should " +
                    "be in the format 'gene1 gene2' for each line, where gene1 and gene2 " +
                    "make up an edge in the network. Defaults to " +
                    "'data/KEGGNetwork.txt' if none given.")
parser.add_argument("-geneRData", required=True, help="Path to the RData file " +
                    "containing all of the gene data.")
parser.add_argument("-items", nargs="*", type=int, help="Use to define which items from " +
                    "geneRData to use in the analysis. Note that this implies that " +
                    "the data must be ordered as desired.")
parser.add_argument("-s", action="store_true", help="Save initial geodesic matrix and parsed " +
                    "expression/grade data created for future uses.")
parser.add_argument("-l", action="store_true", help="Load initial geodesic matrix and parsed " +
                    " expression/grade data from previous use.")
args = parser.parse_args()

# Little logger.
def verboseprint(*print_args, **kwargs):
    print(*print_args, **kwargs) if args.v else lambda *a, **k: None

# Display what's going on.
print("Running GeneSurrounder with network " + args.netfile + " and data " + args.geneRData + ".")

if (len(args.items) > 0):
    out = "Looking at only experiments "
    for i in range(len(args.items) - 1):
        out += str(args.items[i]) + ", "
    out += "and " + str(args.items[-1]) + "."
    print(out)

if args.l:
    print("Loading geodesic matrix and expression/grade data from previous use.")
elif args.s:
    print("Saving geodesic matrix and expression/grade data for future uses.")

# Load in the network file.
network = nx.read_edgelist(args.netfile)
verboseprint("Finished loading network!")

# Calculate geodesic distances for the network.
if args.l:
    distances = load_obj("distances")
    verboseprint("Finished loading distances.")
else:
    distances = nx.floyd_warshall(network)
    verboseprint("Finished computing geodesic distance matrix!")
    
    if args.s:
        save_obj(distances, "distances")
        verboseprint("Finished saving geodesic distance matrix.")
    
print(distances)

# Load in the gene expression data sets along with the phenotype grade values.
if args.l:
    expr = load_obj("expr")
    grade = load_obj("grade")
    verboseprint("Finished loading expr and grade data.")
else:
    [expr, grade] = import_cur_data(args.geneRData)
    verboseprint("Finished loading gene expression and grade data.")
    
    if args.s:
        save_obj(expr, "expr")
        save_obj(grade, "grade")
        verboseprint("Finished saving gene expression and grade data.")

# Keep only the experiments desired.
if (len(args.items) > 0):
    expr = [expr[i] for i in args.items]
    grade = [grade[i] for i in args.items]

# Get the putative genes from the intersection of those studied and those in network.
overlap_genes = list(network.nodes.intersection(expr[0].index))

# Get the diameter of the graph.
diameter = nx.diameter(network)

# Run each part of the algorithm to get the respective p values.
p_decay = [0 for _ in range(len(expr))]

for i in range(len(expr)):
    p_decay[i] = decayDE(expr[i], grade[i], 4, distances)
verboseprint("Finished computing Decay of Differential expression value!")

p_SI = [0 for _ in range(len(expr))]

for i in range(len(expr)):
    p_SI[i] = sphereOfInf(distances, expr[i], diameter, overlap_genes, 200)
verboseprint("Finished computing Sphere of Influence value!")

# Compute the final p values.
for i in range(len(expr)):
    print(combine_pvalues([p_decay[i], p_SI[i]]))
