Introduction
===============

ENCORE is a Python package designed to quantify the similarity between
conformational ensembles of proteins (or in principle other
macromolecules), using three different methods originally described
in:
::
   Kresten Lindorff-Larsen, Jesper Ferkinghoff-Borg (2009) 
   Similarity Measures for Protein Ensembles. 
   PLoS ONE 4(1): e4203. doi:10.1371/journal.pone.0004203

A description of ENCORE and a number of application can be found in:
::
   Matteo Tiberti, Elena Papaleo, Tone Bengtsen, Wouter Boomsma and 
   Kresten Lindorff-Larsen, 
   ENCORE: Software for quantitative ensemble comparison
   Submitted

The package includes facilities for handling ensembles and
trajectories, performing clustering or dimensionality reduction of the
ensemble space, estimating multivariate probability distributions from
the input data, and more. ENCORE can be used to compare experimental
and simulation-derived ensembles, as well as estimate the convergence
of trajectories from time-dependent simulations. The package was
designed as a Python 2.6 (or any higher 2.X version) library. The user
may also use some of the library files as scripts that accept command
line arguments. Usually, the help text included for each script
(obtained running “python encore/script.py -h”) is
self-explanatory. Examples are also available on how ENCORE may be
used to calculate the similarity measures on a number of ensembles.

The similarity measures implemented in ENCORE are based on three
different methods, which all rely on the following idea: Given two or
more conformational ensembles of the same topology (i.e. structure),
we view the particular set of conformations from each ensemble as a
sample from an underlying, but unknown, probability distribution. We
use this sample to model the probability density function of said
distribution. Then we compare the modeled distributions using standard
measures of the similarity between two probability densities, such as
the Jensen-Shannon divergence.
  
In the ENCORE package, we have implemented three methods to estimate
the density:

* Harmonic ensembles similarity (HES): we assume that each ensemble is
  derived from a multivariate normal distribution. We, thus, estimate
  the parameters for the distribution of each ensemble (mean and
  covariance matrix) and compare them using a symmetrized version of
  the Kullback-Leibler divergence. For each ensemble, the mean
  conformation is estimated as the average over the ensemble, and the
  covariance matrix is calculated by default using a shrinkage
  estimate method (or by a maximum-likelihood method, optionally).

* Clustering-based similarity (CES): We use the affinity propagation
  method for clustering to partition the whole space of conformations
  in to clusters of structures. After the structures are clustered we
  take the population of each ensemble in each cluster as a
  probability distribution of conformations. We then compare the
  obtained probability distribution using the Jensen-Shannon
  divergence measure between probability distributions.

* Dimensionality reduction-based similarity (DRES): We use a gaussian
  kernel-based density estimation method to estimate the probability
  density, and use that as probability function in order to compare
  different ensembles. Before doing that, however, due to the limited
  size of the sample, it is necessary to reduce the dimensionality of
  the input space. Thus, the method first projects the ensembles into
  lower dimensions by using the Stochastic Proximity Embedding
  algorithm.

ENCORE is able to use, as input data, structural ensembles deriving
both from molecular simulations (e.g. molecular dynamics or Monte
Carlo methods) or experimental structural ensembles (e.g. NMR
structures as PDB files). The software is able to handle the most
popular trajectory formats (files such as DCD, XTC, TRR, XYZ, TRJ,
MDCRD), although periodic boundaries conditions must be removed before
use. A topology file is also required.

Together with the software, we also provide three examples that
showcase three typical cases of study:

* comparing simulation trajectories with other trajectories 
* estimating convergence of trajectories from molecular dynamics simulations
* comparing  experimentally-derived ensembles from the PDB

See the examples themselves for more information.
If you use ENCORE for your scientific work, please cite:
::
   Matteo Tiberti, Elena Papaleo, Tone Bengtsen, Wouter Boomsma 
   and Kresten Lindorff-Larsen,
   ENCORE: Software for quantitative ensemble comparison
   Submitted

