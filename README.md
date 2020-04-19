# GeneSurrounder-Python
An implementation of gene surrounder (https://github.com/sahildshah1/gene-surrounder) in Python for our COMP 572 final project.

## Runnning the Program

Can run using : python3 GeneSurrounder.py [RData file] [gene name] [-args]

The RData file is the path to the desired file to use for expression data and grade data, and the gene name
is the specific gene to run the algorithm on.

The optional command line arguments are as follows:

    -v : Verbose flag which enables printing of progress in the program.
    -netfile [path to network file] : Path to a network file to use for the program, defaults to the KEGGNetwork in the data folder.
    -items [space-separated integers] : The indices of experiments from the RData file to run the program on, defaults to all of them.
    -l : Loads the distance matrix and expression/grade data from a previous run. Takes priority over -s.
    -s : Stores the distance matrix and the expression/grade data into a file so they can be loaded in subsequent runs.

Example of running the program for the CurOvGradeKEGGnets.RData file in the data folder, for gene hsa:2316, and experiments 2, 3, and 7, also
with verbosity and saving the generated items (s can be switched to l to load instead):

`python3 GeneSurrounder.py data/CurOvGradeKEGGnets.RData hsa:2316 -items 2 3 7 -v -s`
