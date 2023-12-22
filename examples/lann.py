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

# LANN wing
# Adrien Crovato

import gmshcfd

def build_cfg():
    # Build path
    import os.path
    airf_path = []
    for i in range(8):
        airf_path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'airfoils', f'lann_{i}.dat')))
    # Define wing leading edge coordinates and chord lengths
    le_coords = [[0.0, 0.0, 0.0],
                 [0.10408, 0.2, 0.0],
                 [0.16913, 0.325, 0.0],
                 [0.24720, 0.475, 0.0],
                 [0.33827, 0.65, 0.0],
                 [0.42934, 0.825, 0.0],
                 [0.49439, 0.95, 0.0],
                 [0.52041, 1.0, 0.0]]
    chords = [0.3606, 0.31765, 0.29071, 0.25806, 0.22029, 0.18235, 0.15534, 0.14445]
    # Compute mesh sizes
    sizes = [c / 100 for c in chords]
    ff_size = gmshcfd.utils.compute_ff_mesh_size(sizes[0], 50 * chords[0], 1.2)
    # Build cfg
    cfg = {
        'wings': {
            'wing': {
                'offset': [0., 0.],
                'le_offsets': le_coords,
                'airfoils': airf_path,
                'chords': chords,
                'incidences': [0. for _ in chords]
            }
        },
        'domain': {
            'type': 'rans',
            'length': 50 * chords[0]
        },
        'mesh': {
            'wing_sizes': {
                'wing': {
                    'te': [s for s in sizes],
                    'le': [s for s in sizes]
                }
                #'wing': {
                #    'num_cell_chord': 125,
                #    'num_cell_span': [50, 30, 40, 40, 40, 30, 10],
                #    'prog_chord': 0.25
                #}
            },
            'boundary_layer': {
                'wing': {
                    'num_layer': 30,
                    'growth_ratio': 1.2,
                    'hgt_first_layer': 2e-6
                },
                'write_tags': True
            },
            'domain_size': ff_size
        }
    }
    return cfg

def main():
    # Generate wing and domain geometry
    cfd = gmshcfd.GmshCFD('lann', build_cfg())
    cfd.generate_geometry()

    # Generate mesh
    cfd.generate_mesh()

    # Write mesh (write geometry not supported for extruded boundary layer)
    cfd.write_mesh('msh')

    # eof
    print('')

if __name__ == '__main__':
    main()
