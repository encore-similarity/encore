# Ensemble.py --- Representation of a protein ensemble
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
Ensemble representation --- :mod:`MDAnalysis.analysis.ensemble.ensemble`
=====================================================================

:Author: Wouter Boomsma
:Year: 2012
:Copyright: GNU Public License v3

The module contains the Ensemble class allowing for easy reading in 
and alignment of the ensemble contained in a trajectory file.

"""

import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.analysis.align
import numpy

class Ensemble:

    def __init__(self, 
                 universe=None,
                 topology=None,
                 trajectory=None,
                 atom_selection_string='(name CA)',
                 frame_interval=1):

        if not universe:

            # Chained trajectories cannot use TimeSeries functionality
            # and the analysis is therefore slower - we therefore use a 
            # single trajectory value when possible
            if len(trajectory) == 1:
                trajectory = trajectory[0]            

            universe = MDAnalysis.Universe(topology, 
                                           trajectory)

        # Construct atom selection from string
        self.atom_selection_string = atom_selection_string
        self.atom_selection = universe.selectAtoms(atom_selection_string)

        # Try to extract coordinates using Timeseries object
        # This is significantly faster, but only implemented for certain 
        # trajectory file formats
        try:
            self.coordinates = universe.trajectory.timeseries(self.atom_selection, skip=frame_interval, format='fac')
        except:                     # if the Timeseries extraction fails, fall back to a slower approach
            n_coordinates = 0
            for i,time_step in enumerate(universe.trajectory):
                if (i % frame_interval) == 0:
                    n_coordinates += 1
            self.coordinates = numpy.zeros(tuple([n_coordinates]) + self.atom_selection.coordinates().shape)
            for i, time_step in enumerate(universe.trajectory):
                if (i % frame_interval) == 0:
                    self.coordinates[time_step.frame-1] = self.atom_selection.coordinates(time_step)

    def size(self):
        return self.coordinates.shape[0]
                
    def get_coordinates(self, start=0, end=None):
        return self.coordinates[start:end]

    def get_atom_selection_string(self):
        return self.atom_selection_string

    def get_atom_selection(self, start=0, end=None):
        return self.atom_selection

    def align(self, reference=None, start=0, end=None):

        coordinates = self.get_coordinates(start=start, end=end)
        masses = self.get_atom_selection(start=start, end=end).masses()
        
        # Find center of mass for all frames
        coordinates_center_of_mass = numpy.average(coordinates,
                                                   axis=1,
                                                   weights=masses)

        # Move all structures to their center of mass
        coordinates -= coordinates_center_of_mass[:,numpy.newaxis]

        # If reference is given, align all to that
        offset = start
        if reference:

            # Select the same atoms in reference structure
            reference_atom_selection = reference.selectAtoms(self.atom_selection_string)
            reference_coordinates = reference_atom_selection.atoms.coordinates()[start:end]
            reference_masses = self.atom_selection.masses()[start:end]

            # Reference center of mass
            reference_center_of_mass = numpy.average(reference_coordinates, axis=0,
                                                     weights=reference_masses)

            # Move reference structure to its center of mass 
            reference_coordinates -= reference_center_of_mass
        
        else:
            
            reference_coordinates = coordinates[0]
            
            # Skip first frame when rotating
            offset = 1

        # Apply optimal rotations for each frame
        for frame_coordinates in coordinates[offset:]:

            rotation_matrix = MDAnalysis.analysis.align.rotation_matrix(frame_coordinates, 
                                                                        reference_coordinates, 
                                                                        masses)[0]
            frame_coordinates[:] = numpy.transpose(numpy.dot(rotation_matrix, 
                                                             numpy.transpose(frame_coordinates)))


