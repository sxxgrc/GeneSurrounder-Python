import pandas as pd
import os
import networkx as nx
from import_cur_data import import_cur_data

[expr, csv] = import_cur_data('/data_auto_import/CurOvGradeKEGGnets.RData')



print(expr[0].head())
print(csv[0].head())
