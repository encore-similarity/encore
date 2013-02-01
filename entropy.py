# entropy.py --- Ensemble entropy calculation using the Schlitter approximation
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
Calculate entropy --- :mod:`MDAnalysis.analysis.ensemble.entropy`
=====================================================================

:Author: Wouter Boomsma
:Year: 2013
:Copyright: GNU Public License v3

The module contains functions to calculate the structural entropy of
an ensemble using the Schlitter approximation, optionally as it develops 
over time.

"""

from Ensemble import Ensemble
from covariance import covariance_matrix, EstimatorShrinkage, EstimatorML

import numpy

def calc_entropy(ensemble,
                 estimator = EstimatorShrinkage(),
                 align=True,
                 mass_weighted=True,
                 reference = None,
                 start=0,
                 end=None,
                 temperature=300, # Kelvin
                 omit_last_eigenvalues = 0):

    '''
    Calculate entropy from a covariance matrix, using the technique described in

     Estimation of absolute and relative entropies of macromolecules using 
     the covariance matrix. Jurgen Schlitter, Chemical Physics Letters, 215, 1993.

    * ensemble - the structural ensemble
    * estimator - which estimator type to use (maximum likelihood, shrinkage)
    * align - whether to align the structures first
    * mass_weighted - whether to do a mass-weighted analysis
    * reference - use the distances to a specific reference structure rather
                  than the distance to the mean.
    * start - Start frame index
    * end - End frame index
    * temperature - Temperature (Kelvin)
    * omit_last_eigenvalues - How many of the smallest eigen values to omit (these are
                              close to zero when aligning the ensemble and can be omitted)

    '''

    if align:
        ensemble.align(reference, start=start, end=end)        

    # Calculate covariance matrix
    sigma = covariance_matrix(ensemble, estimator, mass_weighted, reference, start=start, end=end)       # AA^2 * amu

    # print "trace: ", numpy.trace(sigma)

    # Find eigen values and vectors
    eigen_values,eigen_vectors = numpy.linalg.eig(sigma)

    # Sort in decreasing eigen value size
    sorted_indices = eigen_values.argsort()[::-1]
    eigen_values = eigen_values[sorted_indices]
    eigen_vectors = eigen_vectors[sorted_indices]

    # Planck's constant (Js)
    h = 6.62606957E-34
    h_bar = h/(2*numpy.pi)

    # Temperatur (K)
    T = temperature  

    # Boltzmann's constant (J/K) 
    k = 1.3806488E-23

    # Euler's number
    e = numpy.e

    # Atomic mass unit
    amu = 1.6605402E-27   # kg/amu

    aangstroem = 1E-10    # meter/AA
    
    prefactor = k*T*e*e/(h_bar*h_bar) * amu * aangstroem*aangstroem  #  J/(Js)^2 * kg/amu * meter^2/AA^2 = 1/(Js^2) * kg/amu * meter^2/AA^2

    # print prefactor

    avogadro = 6.0221367E23 # 1/mol

    last_index = None
    if omit_last_eigenvalues > 0:
        last_index = -omit_last_eigenvalues

    determinant = 0
    for i, eigen_value in enumerate(eigen_values[:last_index]):
        # print i, eigen_value, determinant
        determinant += numpy.log(prefactor*eigen_value + 1)        # J/(Js)^2 * kg/amu * meter^2/AA^2 * (AA^2*amu) = 1/J * 1/s^2 * kg * meter^2 = 1

    # print k, avogadro, k*avogadro, determinant

    entropy = 0.5*k*avogadro*determinant # J/K * 1/mol = J/(K mol)

    return entropy
    # print "entropy: ", entropy

    # print h_bar*h_bar

    # print numpy.linalg.det(k*T*e*e/(h_bar*h_bar)*sigma + numpy.identity(sigma.shape[0]))
    
    # return 0.5*k*numpy.log(numpy.linalg.det(k*T*e*e/(h_bar*h_bar)*sigma + 1))
    

def calc_entropy_progression(ensemble,
                             estimator = EstimatorShrinkage(),
                             align=True,
                             mass_weighted=True,
                             reference = None,
                             temperature=300, # Kelvin
                             omit_last_eigenvalues = 0,
                             realign=False,
                             steps = 10):
    
    '''
    Calculate, for 0-t, with increasing t - the entropy from a covariance matrix, using the
    technique described in

     Estimation of absolute and relative entropies of macromolecules using 
     the covariance matrix. Jurgen Schlitter, Chemical Physics Letters, 215, 1993.

    * ensemble - the structural ensemble
    * estimator - which estimator type to use (maximum likelihood, shrinkage)
    * align - whether to align the structures first
    * mass_weighted - whether to do a mass-weighted analysis
    * reference - use the distances to a specific reference structure rather
                  than the distance to the mean.
    * temperature - Temperature (Kelvin)
    * omit_last_eigenvalues - How many of the smallest eigen values to omit (these are
                              close to zero when aligning the ensemble and can be omitted)
    * realign - Whether to realign in each step
    * steps - How many steps to make in the progression

    '''

    if align and not realign:
        # Align only once
        ensemble.align(reference)
        align = False

    step_size = ensemble.size()/steps
    result = []
    for end in xrange(step_size, ensemble.size(), step_size):
        entropy = calc_entropy(ensemble,
                               estimator=estimator,
                               align=align,
                               mass_weighted=mass_weighted,
                               reference=reference,
                               temperature=temperature,
                               omit_last_eigenvalues=omit_last_eigenvalues,
                               start=0,
                               end=end)
        result.append([end,entropy])
    return numpy.array(result)


if __name__ == "__main__":

    import optparse
    import utils

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
                      help="Output file for covariance matrix")
    parser.add_option("--covariance-estimator", type="choice", dest="covariance_estimator", default="shrinkage",
                      choices=["ml","shrinkage"],
                      help="Type of estomator (maximum likelihood (ml) or shrinkage")
    parser.add_option("--omit-last-eigenvalues", dest="omit_last_eigenvalues", default=6,
                      help="How many of the smallest eigenvalues to skip in the entropy calculation. When aligning the structures, 6 eigenvalues will become very small, and can be omitted from the entropy calculation")
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

    # When aligning the ensemble, 6 eigen values become
    # very small, and can be omitted from the calculation.
    # When not doing alignment, we include all values
    if options.no_align:
        options.omit_last_eigenvalues = 0

    # Select covariance estimator
    estimator = EstimatorML()
    if options.covariance_estimator == "shrinkage":
        estimator = EstimatorShrinkage()

    # Disable reference unless use_distance_to_reference is set
    if not options.use_distance_to_reference:
        reference = None
    
    # Calculate entropy 
    entropy = calc_entropy(ensemble,
                           estimator = estimator,
                           align = (not options.no_align),
                           mass_weighted = options.mass_weighted_analysis,
                           reference = reference,
                           omit_last_eigenvalues = options.omit_last_eigenvalues)
    
    print "Entropy for complete trajectory: ", entropy, "J/(K mol)"

    calc_entropy_progression(ensemble,
                             estimator = estimator,
                             align = (not options.no_align),
                             mass_weighted = options.mass_weighted_analysis,
                             reference = reference,
                             omit_last_eigenvalues = options.omit_last_eigenvalues)
