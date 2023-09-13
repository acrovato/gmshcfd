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

# ONERA M6 wing
# Adrien Crovato

import gmshcfd

def build_cfg():
    # Build path
    import os.path
    apath = os.path.abspath(os.path.join(os.path.dirname(__file__), 'airfoils', 'onera_m6.dat'))
    airf_path =  [apath, apath]
    # Compute wing leading edge coordinates and chord lengths
    le_coords, chords = gmshcfd.utils.compute_wing_cfg([1.196], [0.56], [30], [0], 0.8059)
    # Compute mesh sizes
    root_size = chords[0] / 100
    tip_size = chords[1] / 100
    ff_size = gmshcfd.utils.compute_ff_mesh_size(root_size, 20 * chords[0], 1.2)
    # Build cfg
    cfg = {
        'wings': {
            'wing': {
                'offset': [0., 0.], # x and z offset at the leading edge root
                'le_offsets': le_coords, # x, y and z coordinates of each airfoil's leading edge, relative to offset
                'airfoils': airf_path, # name of (Selig formatted) files containing airfoil coordinates
                'chords': chords, # chord of each airfoil
                'incidences': [0., 0.] # incidence angle of each airfoil
            }
        },
        'domain': {
            'type': 'potential', # nature of equations to be sovled
            'length': 20 * chords[0] # extent of domain
        },
        'mesh': {
            'wing_sizes': {
                'wing': [[root_size, root_size], [tip_size, tip_size]] # mesh sizes at TE and LE of each airfoil
            },
            'domain_size': ff_size # mesh size at farfield
        }
    }
    return cfg

def main():
    # Generate wing and domain geometry
    cfd = gmshcfd.GmshCFD('m6', build_cfg())
    cfd.generate_geometry()

    # Generate mesh
    cfd.generate_mesh()

    # Write geometry and mesh
    cfd.write_geometry()
    cfd.write_mesh('msh2')

    # eof
    print('')

if __name__ == '__main__':
    main()
