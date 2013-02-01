
This is an extension to the MDAnalysis python package for
analysing ensembles for molecular simulation. So far, the
package contains support for calculating similarity scores, 
covariance matrices (including shrinkage estimation), PCA,
and entropy.

The .py files are libraries, but can also serve as 
executables. For instance: 

python ensemble/similarity.py --trajectory1 traj1.dcd \
                              --trajectory2 traj2.dcd \
                              --topology1 top1.pdb \
                              --topology2 top2.pdb 

When running them as executables, use the --help option to see 
the available command line options.

