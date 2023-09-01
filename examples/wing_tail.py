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

# Generic wing/tail
# Adrien Crovato

import gmshcfd

def build_cfg():
    # Build path
    import os.path
    wairf_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'airfoils', 'rae_2822.dat'))
    tairf_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'airfoils', 'naca_0012.dat'))
    airf_path =  [[wairf_path, wairf_path, wairf_path], [tairf_path, tairf_path]]
    # Compute wing and tail leading edge coordinates and chord lengths
    wle_coords, wchords = gmshcfd.utils.compute_wing_cfg([0.5, 3.0], [0.82, 0.35], [20., 20.], [1., 3.], 1.0)
    tle_coords, tchords = gmshcfd.utils.compute_wing_cfg([1.0], [0.3], [25.], [2.], 0.5)
    # Compute mesh sizes
    wsizes = [c / 100 for c in wchords]
    tsizes = [c / 100 for c in tchords]
    ff_size = gmshcfd.utils.compute_ff_mesh_size(wsizes[0], 50 * wchords[0], 1.2)
    # Build cfg
    cfg = {
        'wings': {
            'wing': {
                'le_offset': [0., 0.],
                'le_coords': wle_coords,
                'coordinates': airf_path[0],
                'chords': wchords,
                'incidences': [2., 1., -1.]
            },
            'tail': {
                'le_offset': [4., 0.5],
                'le_coords': tle_coords,
                'coordinates': airf_path[1],
                'chords': tchords,
                'incidences': [0., 0.]
            }
        },
        'domain': {
            'type': 'euler',
            'length': 50 * wchords[0]
        },
        'mesh': {
            'wing_sizes': {
                'wing': [[s, s] for s in wsizes],
                'tail': [[s, s] for s in tsizes]
            },
            'domain_size': ff_size
        }
    }
    return cfg

def main():
    # Generate wing and domain geometry
    cfd = gmshcfd.GmshCFD('wing_tail', build_cfg())
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
