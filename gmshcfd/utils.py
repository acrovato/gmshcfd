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
from scipy.interpolate import PchipInterpolator
import os, os.path

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

def add_section(le_coords, chords, incidences, airf_path, y_sec=0.99):
    """Add a section to the wing planform

    Parameters:
    le_coords: list of list of float
        leading edge coordinates of each cross-section
    chords: list of float
        chord length of each cross-section
    incidences: list of float
        angles of incidence of each cross-section
    airf_path : list of string
        path to files containing airfoil coordinates of each cross-section
    y_sec: float
        normalized y-coordinate of cross-section to add (between 0 and 1)

    Return:
    le_coords: list of list of float
        leading edge coordinates with added cross-section
    chords: list of float
        chord lengths with added cross-section
    incidences: list of float
        angles of incidence with added cross-section
    airf_path : list of string
        path to files containing airfoil coordinates with added cross-section
    """
    # Check if input is in range
    if y_sec <= 0. or y_sec >= 1.:
        raise RuntimeError('add_section: Parameter "y_sec" must be between 0 and 1!\n')

    # Compute y-coordinate of section and find sections between which to interpolate
    y = y_sec * le_coords[-1][1]
    for isec in range(len(le_coords)):
        if y == le_coords[isec][1]:
            print('Section to add corresponds to one of the input sections, nothing done.')
            return le_coords, chords, airf_path
        if y - le_coords[isec][1] < 0.:
            break

    # Read airfoil coordinates, header and file name
    fnames = []
    header = []
    airfoil_coords = []
    for ipth in [isec - 1, isec]:
        fnames.append(os.path.basename(airf_path[ipth]).split('.')[0])
        file = open(airf_path[ipth])
        header.append(file.readline().strip())
        file.close()
        airfoil_coords.append(np.loadtxt(airf_path[ipth], skiprows=1))

    # Compute interpolation parameter
    a = (y - le_coords[isec - 1][1]) / (le_coords[isec][1] - le_coords[isec - 1][1])
    # Interpolate leading edge coordinates
    if le_coords[isec - 1][0] == le_coords[isec][0] and le_coords[isec - 1][2] == le_coords[isec][2]:
        new_le_coords = [le_coords[isec - 1][0], y, le_coords[isec - 1][2]]
    else:
        new_le_coords = list((1 - a) * np.array(le_coords[isec - 1]) + a * np.array(le_coords[isec]))
    # Interpolate chord length
    if chords[isec - 1] == chords[isec]:
        new_chord = chords[isec - 1]
    else:
        new_chord = (1 - a) * chords[isec - 1] + a * chords[isec]
    # Interpolate incidence angle
    if incidences[isec - 1] == incidences[isec]:
        new_incidence = incidences[isec - 1]
    else:
        new_incidence = (1 - a) * incidences[isec - 1] + a * incidences[isec]
     # Interpolate airfoil coordinates
    if fnames[0] == fnames[1]:
        new_fpath = airf_path[0]
    else:
        new_airfoil_coords = (1 - a) * interpolate_coords(airfoil_coords[0]) + a * interpolate_coords(airfoil_coords[1])
        new_fname = f'{fnames[0]}_{fnames[1]}.dat'
        np.savetxt(new_fname, new_airfoil_coords, fmt='%1.5e', header=f'{header[0]} - {header[1]}')
        new_fpath = os.path.join(os.getcwd(), new_fname)

    # Add new section to lists
    le_coords.insert(isec, new_le_coords)
    chords.insert(isec, new_chord)
    incidences.insert(isec, new_incidence)
    airf_path.insert(isec, new_fpath)
    return le_coords, chords, incidences, airf_path

def interpolate_coords(coords, n_pts=80, flip=False):
    """Interpolate coordinates

    Parameters:
    coords: ndarray
        airfoil coordinates
    n_pts: int
        number of interpolation points along x-axis
    flip: bool
        whether to reverse order of coordinates along x-axis

    Return:
    new_coords: ndarray
        interpolated coordinates
    """
    # Revert coordinates order if not Selig
    if flip:
        coords = np.flipud(coords)

    # Split suction and pressure sides
    le = np.argmin(coords[:, 0])
    x_u = coords[:le+1, 0]
    y_u = coords[:le+1, 1]
    x_l = coords[le:, 0]
    y_l = coords[le:, 1]

    # Interpolate using half-cosine spacing
    x_new = 0.5 * (1 - np.cos(np.linspace(0., np.pi, n_pts)))
    y_u_new = PchipInterpolator(np.flipud(x_u), np.flipud(y_u))(x_new)
    y_l_new = PchipInterpolator(x_l, y_l)(x_new)

    # Concatenate and remove duplicated LE
    new_coords = np.vstack((np.hstack((np.flipud(x_new[1:]), x_new)),
                            np.hstack((np.flipud(y_u_new[1:]), y_l_new))))
    return new_coords.T

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
    fname = os.path.basename(fpath).split('.')[0] + '_ste.dat'
    np.savetxt(fname, new_coords, fmt='%1.5e', header=header + ' sharp TE')

    # Display
    if gui:
        import matplotlib.pyplot as plt
        plt.plot(coords[:, 0], coords[:, 1], marker='o')
        plt.plot(new_coords[:, 0], new_coords[:, 1], marker='o', mfc=None)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

def compute_le_sweep(sweep, chord, span, taper, chord_frac=0.25)
    """Compute the leading edge sweep angle from the sweep angle at another chord fraction

    Parameters:
    sweep: float
        sweep angle at given chord fraction
    chord: float
        chord length
    span: float
        semi-span length
    taper: float
        taper ratio
    chord_frac: float
        fraction of the chord at which sweep is defined (between 0 and 1)
    """
    return np.arctan(np.tan(sweep) + 0.25 * chord / span * (1 - taper))
