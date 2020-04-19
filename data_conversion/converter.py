import pandas as pd

nodes = pd.read_csv('data/Sample3.txt')
edgelist = pd.read_csv('data/largestCompKEGGigraph.txt', sep=' ', header = None)
edgelist = edgelist.replace(nodes.index, nodes.values)

with open("KEGGNetwork.txt", "w") as f:
    cur_row = None
    for i in range(edgelist.shape[0]):
        cur_row = edgelist.iloc[i]
        f.write(cur_row[0] + " " + cur_row[1] + "\n")
