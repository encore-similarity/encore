# similarity.py --- Simularity measures between protein ensembles
# Copyright (C) 2012 Wouter Boomsma
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Ensemble similarity calculations --- :mod:`MDAnalysis.analysis.ensemble.similarity`
=====================================================================

:Author: Wouter Boomsma
:Year: 2012
:Copyright: GNU Public License v3

The module contains implementations of similary measures between
protein ensembles described in:

     Similarity Measures for Protein Ensembles. Lindorff-Larsen, K.; 
     Ferkinghoff-Borg, J. PLoS ONE 2009, 4, e4203.

"""

import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.analysis.align
import numpy

from Ensemble import Ensemble
from covariance import covariance_matrix, EstimatorShrinkage, EstimatorML


def harmonic_ensemble_similarity(ensemble1,
                                 ensemble2,
                                 mass_weighted=True,
                                 covariance_estimator = EstimatorShrinkage()):
    ''' 
    Calculate the harmonic ensemble similarity measure
    as defined in 

    Similarity Measures for Protein Ensembles. Lindorff-Larsen, K.; 
    Ferkinghoff-Borg, J. PLoS ONE 2009, 4, e4203.

    '''


    # Extract coordinates from ensembles
    coordinates_system1 = ensemble1.get_coordinates()
    coordinates_system2 = ensemble2.get_coordinates()
    
    # Average coordinates in the two systems
    x1 = numpy.average(coordinates_system1, axis=0).flatten()
    x2 = numpy.average(coordinates_system2, axis=0).flatten()

    # Covariance matrices in the two systems
    sigma1 = covariance_matrix(ensemble1, 
                               mass_weighted=mass_weighted,
                               estimator = covariance_estimator)
    sigma2 = covariance_matrix(ensemble2, 
                               mass_weighted=mass_weighted,
                               estimator = covariance_estimator)

    # Inverse covariance matrices
    sigma1_inv = numpy.linalg.pinv(sigma1)
    sigma2_inv = numpy.linalg.pinv(sigma2)

    # Difference between average vectors
    d_avg = x1 - x2

    # Sigma
    sigma = sigma1_inv + sigma2_inv

    # Distance measure
    trace = numpy.trace(numpy.dot(sigma1_inv, sigma2) + 
                        numpy.dot(sigma1, sigma2_inv)
                        - 2*numpy.identity(sigma1.shape[0]))
    d_hes = 0.25*(numpy.dot(numpy.transpose(d_avg), 
                            numpy.dot(sigma1_inv + sigma2_inv,
                                      d_avg)) + trace)

    return d_hes
    
    


if __name__ == "__main__":

    import optparse
    import utils

    usage = "usage %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option("--atom-selection", dest="atom_selection_string", default="(name CA)",
                      help="CHARMM-style atom selection")
    parser.add_option("--topology1", dest="topology1_filename", default="",
                      help="Ensemble 1: Path to a topology file (supported formats: PDB,PSF,CRD,GRO)")
    parser.add_option("--trajectory1", dest="trajectory1_filename", action="callback", type="string", nargs=0, default=[],
                      callback=utils.vararg_callback,
                      metavar="TRAJECTORY_FILENAME(S)",
                      help="Ensemble 1: Add trajectory filenames")
    parser.add_option("--frame-interval1", dest="frame_interval1", default=1,
                      help="Ensemble1: The interval between frames used for the analysis")
    parser.add_option("--topology2", dest="topology2_filename", default="",
                      help="Ensemble 1: Path to a topology file (supported formats: PDB,PSF,CRD,GRO)")
    parser.add_option("--trajectory2", dest="trajectory2_filename", action="callback", type="string", nargs=0, default=[],
                      callback=utils.vararg_callback,
                      metavar="TRAJECTORY_FILENAME(S)",
                      help="Ensemble 2: Add trajectory filenames")
    parser.add_option("--frame-interval2", dest="frame_interval2", default=1,
                      help="Ensemble2: The interval between frames used for the analysis")
    parser.add_option("--no-align", dest="align", action="store_false", default=True,
                      help="Whether to align ensembles internally before calculating similarity")
    parser.add_option("--covariance-estimator", type="choice", dest="covariance_estimator", default="shrinkage",
                      choices=["ml","shrinkage"],
                      help="Type of estomator (maximum likelihood (ml) or shrinkage")
    (options, args) = parser.parse_args()

    if (not options.trajectory1_filename or not options.topology1_filename or 
        not options.trajectory2_filename or not options.topology2_filename):
        parser.error("--trajectory and --topology options must be specified for both ensembles")

    # Construct reference systems
    try:
        reference1 = MDAnalysis.Universe(options.topology1_filename, 
                                         options.topology1_filename)
    except:
        reference1 = None

    try:
        reference2 = MDAnalysis.Universe(options.topology2_filename, 
                                         options.topology2_filename)
    except:
        reference2 = None


    # Construct ensembles
    ensemble1 = Ensemble(topology=options.topology1_filename,
                         trajectory=options.trajectory1_filename,
                         atom_selection_string=options.atom_selection_string,
                         frame_interval=options.frame_interval1)
    ensemble2 = Ensemble(topology=options.topology2_filename,
                         trajectory=options.trajectory2_filename,
                         atom_selection_string=options.atom_selection_string,
                         frame_interval=options.frame_interval2)

    if options.align:
        ensemble1.align(reference1)
        ensemble2.align(reference2)
    
    # Calculate covariance matrix
    covariance_estimator = EstimatorML()
    if options.covariance_estimator == "shrinkage":
        covariance_estimator = EstimatorShrinkage()    

    # Calculate similarity
    print harmonic_ensemble_similarity(ensemble1, ensemble2,
                                       covariance_estimator=covariance_estimator)

