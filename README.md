=====================
  ENCORE
=====================

ENCORE is a Python package designed to calculate similarity measures between protein ensembles, as reported in: 

       Lindorff-Larsen K, Ferkinghoff-Borg J (2009) 
       Similarity Measures for Protein Ensembles. 
       PLoS ONE 4(1): e4203. doi:10.1371/journal.pone.0004203

The package includes facilities for handling ensembles and trajectories, performing clustering or dimensionality reduction of the ensemble space, estimating multivariate probability distributions from the input data, calculate conformational entropies, principal component analysis and more. The package was designed as a Python 2.6 library, but most library files include their own __main__ section which allows them to be run as user scripts that accept command line arguments. Usually, the included help strings (python encore/script.py -h) is pretty much self-explanatory. An example on how the similarity measures are calculated on a number of ensembles is also available.

ENCORE is able to use as input data files deriving from simulations and experimental data (e.g. from NMR structure resolution experiments). The software is able to handle the most popular trajectory formats (as xtc or dcd files), although the periodic boundaries condition must be solved before being used, as well as multimodel PDB files.

