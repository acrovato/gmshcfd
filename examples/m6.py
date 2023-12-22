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
    sizes = [c / 100 for c in chords]
    ff_size = gmshcfd.utils.compute_ff_mesh_size(sizes[0], 20 * chords[0], 1.2)
    # Build cfg
    cfg = {
        'wings': {
            'wing': {
                'offset': [0., 0.],
                'le_offsets': le_coords,
                'airfoils': airf_path,
                'chords': chords,
                'incidences': [0., 0.]
            }
        },
        'domain': {
            'type': 'potential',
            'length': 20 * chords[0]
        },
        'mesh': {
            'wing_sizes': {
                'wing': {
                    'te': [s for s in sizes],
                    'le': [s for s in sizes]
                }
            },
            'domain_size': ff_size
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
