#!/usr/bin/env python
# coding: utf-8

# In[1]:





# In[2]:





# In[3]:


def import_cur_data(filename = 'CurOvGradeKEGGnets.RData'):
    import pandas as pd
    import subprocess
    import os
    import csv
    path = os.getcwd()
    command = 'Rscript'
    path2script = path + '/read_cur_data.R'
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

import_cur_data()


# In[4]:


# print(len(expr[0]))
# expr[0].head()


# In[5]:


# grade[1].head()


# In[ ]:




