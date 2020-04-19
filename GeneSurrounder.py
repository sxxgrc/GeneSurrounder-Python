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
parser.add_argument("geneRData", help="Path to the RData file " +
                    "containing all of the gene data.")
parser.add_argument("gene", help="Target gene for the computation.")
parser.add_argument("-v", action="store_true", help="Verbosity of the program.")
parser.add_argument("-netfile", nargs="?", default="data/KEGGNetwork.txt", 
                    help="Path to the network model of cellular interactions. Should " +
                    "be in the format 'gene1 gene2' for each line, where gene1 and gene2 " +
                    "make up an edge in the network. Defaults to " +
                    "'data/KEGGNetwork.txt' if none given.")
parser.add_argument("-items", nargs="*", type=int, help="Use to define which items from " +
                    "geneRData to use in the analysis. Note that this implies that " +
                    "the data must be ordered as desired.")
parser.add_argument("-l", action="store_true", help="Load parsed expression/grade data " +
                    "from previous use.")
parser.add_argument("-s", action="store_true", help="Save parsed expression/grade data " + 
                    "created for future uses.")
parser.add_argument("-lg", action="store_true", help="Load distance matrix and diameter " +
                    "for current gene from previous use.")
parser.add_argument("-sg", action="store_true", help="Save distance matrix and diameter " +
                    "for current gene for future use.")
args = parser.parse_args()

# Little logger.
def verboseprint(*print_args, **kwargs):
    print(*print_args, **kwargs) if args.v else lambda *a, **k: None

# Display what's going on.
print("Running GeneSurrounder with network " + args.netfile + " and data " + args.geneRData + 
      " for the gene " + args.gene + ".")

if (args.items != None):
    out = "Looking at only experiments "
    for i in range(len(args.items) - 1):
        out += str(args.items[i]) + ", "
    out += "and " + str(args.items[-1]) + "."
    print(out)

if args.l:
    print("Loading expression/grade data from previous use.")
elif args.s:
    print("Saving expression/grade data for future uses.")
    
if args.lg:
    print("Loading distance matrix and network diameter for current gene.")
elif args.sg:
    print("Saving distance matrix and network diameter for current gene.")

# Load in the network file.
network = nx.read_edgelist(args.netfile)
verboseprint("Finished loading network!")

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
overlap_genes = set(network.nodes()).intersection(expr[0].index)

if args.gene not in overlap_genes:
    raise ValueError("Desired gene not present in both the network and the expression data!")

overlap_genes = list(overlap_genes)
verboseprint("Finished computing overlapping genes between network and data.")

# Calculate geodesic distances for the network.
if args.lg:
    distances = load_obj("distances")
    verboseprint("Finished loading distances.")
else:
    distances = {}
    
    for node in overlap_genes:
        distances[node] = nx.shortest_path_length(network, source=args.gene, target=node)
        
    verboseprint("Finished computing geodesic distance matrix over assayed genes!")

    if args.sg:
        save_obj(distances, "distances")
        verboseprint("Finished saving geodesic distance matrix over assayed genes.")

# Get the diameter of the graph.
if args.lg:
    diameter = load_obj("diameter")
    verboseprint("Finished loading diameter " + str(diameter) + ".")
else:
    diameter = max(distances.values())
    verboseprint("Finished computing network diameter!")
    
    if args.sg:
        save_obj(diameter, "diameter")
        verboseprint("Finished saving network diameter.")

# Run each part of the algorithm to get the respective p values.
p_decay = [0 for _ in range(len(expr))]

for i in range(len(expr)):
    p_decay[i] = decayDE(expr[i], grade[i], distances, diameter)
verboseprint("Finished computing Decay of Differential expression values!")
verboseprint("Values: " + str(p_decay))

p_SI = [0 for _ in range(len(expr))]

for i in range(len(expr)):
    p_SI[i] = sphereOfInf(distances, expr[i], diameter, overlap_genes, 200, args.gene)
verboseprint("Finished computing Sphere of Influence values!")
verboseprint("Values: " + str(p_SI))

# Compute the final p values.
print("Final p values:")
for i in range(len(expr)):
    combined_vals = [combine_pvalues([p_decay[i][j], p_SI[i][j]]) for j in range(diameter)]
    verboseprint(combined_vals)
    
    # Display output.
    item = i if args.items == None else args.items[i]
    print("p value for data set " + str(item) + " : " + str(min(combined_vals)))
