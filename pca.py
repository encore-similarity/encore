# pca.py --- Principal Component Analysis
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
Covariance calculation --- :mod:`MDAnalysis.analysis.ensemble.pca`
=====================================================================

:Author: Wouter Boomsma
:Year: 2012
:Copyright: GNU Public License v3

The module contains functions to do an principle component analysis of
an ensemble of structures.

"""

import numpy

from Ensemble import Ensemble
from covariance import covariance_matrix, EstimatorShrinkage, EstimatorML


def pca(ensemble,
        estimator = EstimatorShrinkage(),
        mass_weighted=True,
        reference = None):
    
    '''
    Calculates (optionally mass weighted) covariance matrix

    * ensemble - the structural ensemble
    * estimator - which estimator type to use (maximum likelihood, shrinkage)
    * mass_weighted - whether to do a mass-weighted analysis
    * reference - use the distances to a specific reference structure rather
                  than the distance to the mean.

    '''
    
    # Calculate covariance matrix
    sigma = covariance_matrix(ensemble, estimator, mass_weighted, reference)

    # Find eigen values and vectors
    eigen_val,eigen_vec = numpy.linalg.eig(sigma)

    # Sort in decreasing eigen value size
    sorted_indices = eigen_val.argsort()[::-1]
    eigen_val = eigen_val[sorted_indices]
    eigen_vec = eigen_vec[sorted_indices]
    
    # Extract coordinates from ensemble
    coordinates = ensemble.get_coordinates()

    # Flatten coordinate matrix into n_frame x n_coordinates
    coordinates = numpy.reshape(coordinates, (coordinates.shape[0], -1))

    # Project all coordinates onto all eigen vectors
    projections = numpy.dot(coordinates, eigen_vec)

    return projections




if __name__ == "__main__":

    import optparse
    import utils
    import sys

    usage = "usage %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option("--topology", dest="topology_filename", default="",
                      help="Path to a topology file (supported formats: PDB,PSF,CRD,GRO)")
    parser.add_option("--trajectory", dest="trajectory_filename", action="callback", type="string", nargs=0, default=[],
                      callback=utils.vararg_callback,
                      metavar="TRAJECTORY_FILENAME(S)",
                      help="Add trajectory filenames")
    parser.add_option("--atom-selection", dest="atom_selection_string", default="(name CA)",
                      help="CHARMM-style atom selection")
    parser.add_option("--mass-weighted-analysis", dest="mass_weighted_analysis", action="store_true", default=False,
                      help="Use mass-weighting in the calculation of the covariance matrix.")
    parser.add_option("--no-align", dest="no_align", action="store_true", default=False,
                      help="Do not superimpose structures before calculating covariance.")
    parser.add_option("--frame-interval", dest="frame_interval", default=1,
                      help="The interval between frames used for the analysis")
    parser.add_option("--use-distance-to-reference", dest="use_distance_to_reference", action="store_true", default=False,
                      help="Whether to use the distance to the reference structure rather than the distance to the average structure when calculating covariance matrix.")
    parser.add_option("--output", dest="output_filename", default="", 
                      help="Output file for projection data")
    parser.add_option("--covariance-estimator", type="choice", dest="covariance_estimator", default="shrinkage",
                      choices=["ml","shrinkage"],
                      help="Type of estomator (maximum likelihood (ml) or shrinkage")
    (options, args) = parser.parse_args()

    if not options.trajectory_filename or not options.topology_filename:
        parser.error("--trajectory and --topology options must be specified")


    # Construct reference if available
    try:
        reference = MDAnalysis.Universe(options.topology_filename, 
                                        options.topology_filename)
    except:
        reference = None

    # Construct ensemble
    ensemble = Ensemble(topology=options.topology_filename,
                        trajectory=options.trajectory_filename, 
                        atom_selection_string=options.atom_selection_string,
                        frame_interval=options.frame_interval)

    # Align ensemble to reference
    if not options.no_align:
        ensemble.align(reference)

    # Select covariance estimator
    estimator = EstimatorML()
    if options.covariance_estimator == "shrinkage":
        estimator = EstimatorShrinkage()

    # Disable reference unless use_distance_to_reference is set
    if not options.use_distance_to_reference:
        reference = None

    # Do principle component analysis
    projections = pca(ensemble,
                      estimator = estimator,
                      mass_weighted=options.mass_weighted_analysis,
                      reference=reference)

    # Output
    if options.output_filename:
        print >> sys.stderr, projections
        output_file = open(options.output_filename, 'w')
        numpy.savetxt(output_file, projections, delimiter=" ")
        output_file.close()
    else:
        numpy.savetxt(sys.stdout, projections, delimiter=" ")
    
