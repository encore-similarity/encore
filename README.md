===================== ENCORE =====================

ENCORE is a Python package designed to  quantify the similarity
between protein ensembles, using three different methods detailed in:

   Lindorff-Larsen K, Ferkinghoff-Borg J (2009) 
   Similarity Measures for Protein Ensembles. 
   PLoS ONE 4(1): e4203. doi:10.1371/journal.pone.0004203

The package includes facilities for handling ensembles and
trajectories, performing clustering or dimensionality reduction of the
ensemble space, estimating multivariate probability distributions from
the input data, and more. ENCORE can be used to compare experimental
and simulation-derived ensembles, as well as estimate the convergence
of trajectories from time-dependent simulations. The package was
designed as a Python 2.6 (or any higher 2.X version) library. Some
library files can also be run as user scripts that accept command line
arguments. Usually, the help text included for each script (obtained
running “python encore/script.py -h”) is self-explanatory. Examples on
how the similarity measures can be calculated using ENCORE on a number
of ensembles are also available.

The similarity measures implemented in ENCORE are based on three
different methods, which all rely on the following assumption. Given
two or more ensembles of the same topology (i.e. structure), we view
the particular set of conformations from each ensemble as a sample
from an underlying unknown distribution. We use this sample to model
the probability density function of said distributions. Then we
compare the modeled distributions using standard measures of the
similarity between two probability densities, such as the
Jensen-Shannon divergence.
  
In the ENCORE package, we have  implemented three methods to estimate the density:
Harmonic ensembles similarity: we assume that each ensemble is derived
from a multivariate normal distribution. We, thus, estimate the
parameters for the distribution of each ensemble (mean and covariance
matrix) and compare them using a symmetrized version of the
Kullback-Leibler divergence. For each ensemble, the mean conformation
is estimated as the average over the ensemble, and the covariance
matrix is calculated by default using a shrinkage estimate method (or
by a maximum-likelihood method, optionally).

Clustering-based similarity: We use the affinity propagation method
for clustering to partition the whole space of conformations in
clusters of structures. After the structures are clustered we take the
population of each ensemble in each cluster as a probability
distribution of conformations. We then compare the obtained
probability distribution using the Jensen-Shannon divergence measure
between probability distributions.

Dimensionality reduction-based similarity: We use the gaussian
kernel-based density estimation to estimate the probability density,
and use that as probability function in order to compare different
ensembles. Before doing that, however, due to the limited size of the
sample, it is necessary to reduce the dimensionality of the input
space. This is performed by using the Stochastic Proximity Embedding
algorithm.

ENCORE is able to use, as input data, structural ensembles deriving
from simulations and experimental data (e.g. from NMR structure
resolution experiments, as PDB files). The software is able to handle
the most popular trajectory formats (files such as DCD, XTC, TRR, XYZ,
TRJ, MDCRD), although periodic boundaries conditions must be removed
before use. A topology file is also required.

Together with the software, we also provide three examples that showcase three typical cases of study:

* comparing simulation trajectories with other trajectories 
* estimating convergence of trajectories from MD simulations
* comparing  experimental  data from PDB

See the examples themselves for more information.

If you use ENCORE for your scientific work, please cite:
   Matteo Tiberti, Elena Papaleo, Tone Bengtsen, Wouter Boomsma and Kresten Lindorff-Larsen,	
   ENCORE: Software for quantitative ensemble comparison
   submitted to PLoS Computational Biology

