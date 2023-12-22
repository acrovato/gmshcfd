# -*- coding: utf-8 -*-

# GmshCFD
# Copyright (C) 2023 Adrien Crovato
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

def compute_wing_cfg(spans, tapers, sweeps, dihedrals, root_chord):
    """Compute the leading edge coordinates and the chord length of the airfoil at the tip of each planform
    
    Parameters:
    spans: list of float
        lengths of each planform
    tapers: list of float
        taper ratios of each planform
    sweeps: list of float
        sweep angle of leading edge (in degrees)
    dihedrals: list of float
        dihedral angle of leading edge (in degrees)
    root_chord: float
        chord length of wing root

    Return:
    le_coords: list of list of float
        leading edge coordinates
    chords: list of float
        chord lengths
    """
    le_coords = [[0., 0., 0.]]
    chords = [root_chord]
    for i in range(len(spans)):
        lec = [0., 0., 0.]
        lec[0] = le_coords[i][0] + np.tan(sweeps[i] * np.pi / 180) * spans[i]
        lec[1] = le_coords[i][1] + spans[i]
        lec[2] = le_coords[i][2] + np.tan(dihedrals[i] * np.pi / 180) * spans[i]
        le_coords.append(lec)
        chords.append(chords[i] * tapers[i])
    return le_coords, chords

def compute_ff_mesh_size(surface_size, domain_length, factor):
    """Compute the mesh size at the farfield boundary
    
    Parameters:
    surface_size: float
        mesh size at surface boundaries
    domain_length: float
        distance between surface boundaries and farfield boundary
    factor: float
        progression factor

    Return:
    _: float
        mesh size at farfield boundary
    """
    n = np.log(1 - (1 - factor) * domain_length / surface_size) / np.log(factor)
    return surface_size * pow(factor, n - 1)

def sharpen_te(fpath, n_change=10, gui=False):
    """Convert a blunt trailing edge to a sharp trailing edge and write coordinates

    Parameters:
    fpath: str
        path to file containing airfoil coordinates
    n_change: int
        number of points to be altered
    gui: bool
        whether to display the airfoil
    """
    # Read header and coordinates
    file = open(fpath)
    header = file.readline().strip()
    file.close()
    coords = np.loadtxt(fpath, skiprows=1)
    new_coords = coords.copy()

    # Convert blunt TE to sharp TE
    z_mean = (coords[0:n_change, 1] + coords[-1:-n_change-1:-1, 1]) / 2
    tck = (coords[0:n_change, 1] - coords[-1:-n_change-1:-1, 1]) / 2
    new_tck = np.linspace(0, tck[-1], n_change)
    new_coords[0:n_change, 1] = z_mean + new_tck;
    new_coords[-1:-n_change-1:-1, 1] = z_mean - new_tck;

    # Write file
    import ntpath
    fname = ntpath.basename(fpath).split('.')[0] + '_ste.dat'
    np.savetxt(fname, new_coords, fmt='%1.5e', header=header + ' sharp TE')

    # Display
    if gui:
        import matplotlib.pyplot as plt
        plt.plot(coords[:, 0], coords[:, 1], marker='o')
        plt.plot(new_coords[:, 0], new_coords[:, 1], marker='o', mfc=None)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
