===================== ENCORE =====================

ENCORE is a Python package designed to calculate similarity measures between
protein ensembles, as reported in:

       Lindorff-Larsen K, Ferkinghoff-Borg J (2009) 
       Similarity Measures for Protein Ensembles. 
       PLoS ONE 4(1): e4203. doi:10.1371/journal.pone.0004203

The package includes facilities for handling ensembles and trajectories,
performing clustering or dimensionality reduction of the ensemble space,
estimating multivariate probability distributions from the input data, calculate
conformational entropies, principal component analysis and more. The package was
designed as a Python 2.6 library, but most library files include their own
__main__ section which allows them to be run as user scripts that accept command
line arguments. Usually, the included help strings (python encore/script.py -h)
is pretty much self-explanatory. An example on how the similarity measures are
calculated on a number of ensembles is also available.

The similarity measures implemented in ENCORE are based on three different
methods. Given two or more ensembles of the same topology (i.e. structure),
our strategy is first to estimate a probability density for each ensemble and
subsequently to compare these densities. We thus view a particular set of
conformations as a sample from an underlying distribution and aim to model
this distribution based on the sample at hand. We have thus devised three
methods to estimate the density:

* Harmonic ensembles similarity (hes):    we assume that each ensemble is
derived from a multivariate normal distribution. We thus estimate the
parameters for the distribution of each ensemble (mean and covariance matrix)
and comapre them using a symmetrized version of the Kullback-Leibler
divergence. For each ensemble, the mean conformation was estimated as the
average over the ensemble, and the covariance matrix is calculated by default
using a shrinkage estimate method (or by a maximum-likelihood method,
optionally0

* Clustering-based similarity: We use a clustering method (Affinity propagation)
in order to partition the whole space of conformations in clusters of
structures. After the structures are clustered we take the population of each
ensemble in each cluster as a probability distribution of conformations. We
then compare the obtained PDFs using the well-known Jensen-Shannon divergence
measure between probability distributions.

* Dimensionality reduction-based similarity: We use the gaussian kernel-based
density estimation (KDE) to estimate the probability density, and use that as
probability function in order to compare different ensemble. Before doing
that, however, due to the limited number of sample, it is necessary to reduce
the dimensionality of the input space. This is performed by using the
Stochasic Proximity Embedding (SPE) algorithm.

ENCORE is able to use as input data files deriving from simulations and
experimental data (e.g. from NMR structure resolution experiments). The
software is able to handle the most popular trajectory formats  (files such as
DCD, XTC, TRR, XYZ, TRJ, MDCRD), although periodic boundaries conditions must
be removed before use. A topology file is also required.

Together with the software, we also provide three examples that showcase three
typical usage cases:

* comparing simulation trajectories with other trajectories alone, 
* comparing simulation trajectories with other simulation trajectories and 
experimental PDB data, 
* estimating convergence of trajectories from MD simulations

See the examples themselves for more information.

If you use ENCORE for your scientific work, please cite:

	Matteo Tiberti, Elena Papaleo, Wouter Boomsma and Kresten Lindorff-Larsen,
	ENCORE: Software for quantitative ensemble comparison
	submitted to PLoS Computational Biology



1