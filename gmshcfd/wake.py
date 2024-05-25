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

from .errors import GmshCFDError
import gmsh

class Wake:
    """Wake attached to the trailing edge of a lifting surface

    Parameters:
    pte_ids: dict
        tag and coordinates of trailing edge points
    cte_ids: list
        tag of trailing edge curves
    name: string
        name of the wake
    domain_cfg: dict
        parameters to configure domain
    mesh_size: float
        mesh size at farfield boundary

    Attribute:
    tags: dict
        gmsh tags of remarkable curves/surfaces of the lifting surface
    """
    def __init__(self, pte_ids, cte_ids, name, domain_cfg, mesh_size):
        # Create Gmsh wake model
        wake_length = domain_cfg['length']
        merge = domain_cfg.get('merge_wake')
        if merge not in [None, 'all', 'last']:
            raise GmshCFDError('merge must be either None, "all" or "last"!\n', self)
        self.__create_model(pte_ids, cte_ids, name, wake_length, merge, mesh_size)

    def __create_model(self, pte_ids, cte_ids, name, wake_length, merge, mesh_size):
        """Create wake points, curves and surfaces in Gmsh model
        """
        # Retain points and curves needed to build wake
        if merge == 'all':
            pte_tags = [list(pte_ids.keys())[0], list(pte_ids.keys())[-1]]
            pte_coords = [list(pte_ids.values())[0], list(pte_ids.values())[-1]]
            cte_tags = [[-i for i in cte_ids]]
        elif merge == 'last':
            pte_tags = [*list(pte_ids.keys())[:-2], list(pte_ids.keys())[-1]]
            pte_coords = [*list(pte_ids.values())[:-2], list(pte_ids.values())[-1]]
            cte_tags = [*[[-i] for i in cte_ids[:-2]], [-cte_ids[-2], -cte_ids[-1]]]
        else:
            pte_tags = list(pte_ids.keys())
            pte_coords = list(pte_ids.values())
            cte_tags = [[-i] for i in cte_ids]

        # Create points
        ptags = []
        for te_coord in pte_coords:
            ptags.append(gmsh.model.geo.add_point(wake_length, te_coord[1], te_coord[2]))

        # Create curves
        shed_ctags = []
        for i, te_tag in enumerate(pte_tags):
            shed_ctags.append(gmsh.model.geo.add_line(te_tag, ptags[i]))
        trail_ctags = []
        for i in range(len(ptags) - 1):
            trail_ctags.append(gmsh.model.geo.add_line(ptags[i], ptags[i+1]))

        # Create surfaces
        stags = []
        for i in range(len(cte_tags)):
            cltag = gmsh.model.geo.add_curve_loop([trail_ctags[i], -shed_ctags[i+1], *cte_tags[i], shed_ctags[i]])
            stags.append(gmsh.model.geo.add_plane_surface([cltag]))

        # Add physical groups
        gmsh.model.geo.synchronize()
        gmsh.model.add_physical_group(1, [shed_ctags[-1]], name=name+'Tip')
        gmsh.model.add_physical_group(2, stags, name=name)

        # Add meshing constraints
        for i in range(len(ptags)):
            gmsh.model.mesh.set_size([(0, ptags[i])], mesh_size)

        # Create tag dictionary
        self.tags = {'symmetry_point': ptags[0],
                     'symmetry_curve': shed_ctags[0],
                     'trailing_curves': trail_ctags,
                     'wake_surfaces': stags}
