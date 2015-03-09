encore.dimensionality_reduction package
=======================================

Module contents
---------------

.. automodule:: encore.dimensionality_reduction
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: encore.dimensionality_reduction.stochasticproxembed
    :members:
    :undoc-members:
    :show-inheritance:

.. py:class:: StochasticProximityEmbedding
   :module: encore.dimensionality_reduction.stochasticproxembed

    Stochastic proximity embedding dimensionality reduction algorithm. The algorithm implemented here is described in this paper:

        Dmitrii N. Rassokhin, Dimitris K. Agrafiotis
        A modified update rule for stochastic proximity embedding
        Journal of Molecular Graphics and Modelling 22 (2003) 133–140

    This class is a Cython wrapper for a C implementation (see spe.c)

    .. py:method:: run(self, s, double rco, int dim, double maxlam, double minlam, int ncycle, int nstep, int stressfreq)
       :module: encore.dimensionality_reduction.stochasticproxembed
    Run stochastic proximity embedding.

        **Arguments:**
        
        `s` : encore.utils.TriangularMatrix object
                Triangular matrix containing the distance values for each pair of elements in the original space.       

        `rco` : float
                neighborhood distance cut-off

        `dim` : int
                number of dimensions for the embedded space

        `minlam` : float
                final learning parameter

        `maxlam` : float
                starting learning parameter

        `ncycle` : int
                number of cycles. Each cycle is composed of nstep steps. At the end of each cycle, the lerning parameter lambda is updated.

        `nstep` : int
                number of coordinate update steps for each cycle
        
        **Returns:**

        `space` : (float, numpy.array)
                float is the final stress obtained; the array are the coordinates of the elements in the embedded space 
 
        `stressfreq` : int
                calculate and report stress value every stressfreq cycle

.. py:class:: kNNStochasticProximityEmbedding
   :module: encore.dimensionality_reduction.stochasticproxembed

    k-Nearest Neighbours Stochastic proximity embedding dimensionality reduction algorithm. 
    This is a variation of the SPE algorithm in which neighbourhood is not defined by a distance cut-off; instead, at each step, when a point is randomly chosen to perform coordinate updates, the coordinates of its k nearest neighbours are updated as well.
    This class is a Cython wrapper for a C implementation (see spe.c)

    .. py:method:: run(self, s, int kn, int dim, double maxlam, double minlam, int ncycle, int nstep, int stressfreq)
       :module: encore.dimensionality_reduction.stochasticproxembed

    Run kNN-SPE.

         **Arguments:**

        `s` : encore.utils.TriangularMatrix object
                Triangular matrix containing the distance values for each pair of elements in the original space.

        `kn` : int
                number of k points to be used as neighbours, in the original space

        `dim` : int
                number of dimensions for the embedded space

        `minlam` : float
                final learning parameter

        `maxlam` : float
                starting learning parameter

        `ncycle` : int
                number of cycles. Each cycle is composed of nstep steps. At the end of each cycle, the lerning parameter lambda is updated.

        `nstep` : int
                number of coordinate update steps for each cycle

        **Returns:**

        `space` : (float, numpy.array)
                float is the final stress obtained; the array are the coordinates of the elements in the embedded space

        `stressfreq` : int
                calculate and report stress value every stressfreq cycle
