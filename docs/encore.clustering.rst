encore.clustering package
=========================

Submodules
----------

encore.clustering.Cluster module
--------------------------------

.. automodule:: encore.clustering.Cluster
    :members:
    :undoc-members:
    :show-inheritance:

encore.clustering.affinityprop module
-------------------------------------

.. automodule:: encore.clustering.affinityprop
    :members:
    :undoc-members:
    :show-inheritance:

.. py:class:: class AffinityPropagation:
   :module: encore.clustering.affinityprop

   Affinity propagation clustering algorithm. This class is a Cython wrapper around the Affinity propagation algorithm, which is implement as a C library (see ap.c). The implemented algorithm is described in the paper:
        
        Clustering by Passing Messages Between Data Points.
        Brendan J. Frey and Delbert Dueck, University of Toronto
        Science 315, 972–976, February 2007 

       .. py:method:: run(self, s, preference, double lam, int max_iterations, int convergence, int noise=1)
          :module: encore.clustering.affinityprop

        Run the clustering algorithm. 

        **Arguments:**
        
        `s` : encore.utils.TriangularMatrix object
                Triangular matrix containing the similarity values for each pair of clustering elements. Notice that the current implementation does not allow for asymmetric values (i.e. similarity(a,b) is assumed to be equal to similarity(b,a))

        `preference` : numpy.array of floats or float
                Preference values, which the determine the number of clusters. If a single value is given, all the preference values are set to that. Otherwise, the list is used to set the preference values (one value per element, so the list must be of the same size as the number of elements)
        `lam` : float
                Floating point value that defines how much damping is applied to the solution at each iteration. Must be ]0,1]

        `max_iterations` : int 
                Maximum number of iterations

        `convergence` : int
                Number of iterations in which the cluster centers must remain the same in order to reach convergence

        `noise` : int
                Whether to apply noise to the input s matrix, such there are no equal values. 1 is for yes, 0 is for no. 

        **Returns:**
        
        `elements` : list of int or None
                List of cluster-assigned elements, which can be used by encore.utils.ClustersCollection to generate Cluster objects. See these classes for more details.

