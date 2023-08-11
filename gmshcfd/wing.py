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

from .wake import Wake
import gmsh
import numpy as np

class Wing:
    """Lifting surface

    Parameters:
    name: string
        name of the lifting surface
    cfg: dict
        parameters to configure lifting surface
    domain_cfg: dict
        parameters to configure domain
    mesh_cfg: dict
        parameters to configure mesh

    Attributes:
    height: float
        z_coordinate of the airfoil trailing edge on the symmetry plane
    wake: gmshcfd.Wake object
        wake attached to the lifting surface
    tags: dict
        gmsh tags of remarkable curves/surfaces of the lifting surface
    """
    def __init__(self, name, cfg, domain_cfg, mesh_cfg):
        # Compute airfoil coordinates
        n_airf, coords, le_idx = self.__compute_coordinates(cfg)
        # Create Gmsh wing and wake models
        self.__create_model(n_airf, coords, le_idx, name, cfg, domain_cfg, mesh_cfg)

    def __compute_coordinates(self, cfg):
        """Read airfoil coordinates from files and transform them using wing planform configuration
        """
        # Get number of airfoils
        n_airf = len(cfg['coordinates'])

        # Read coordinates
        coords = []
        for fname in cfg['coordinates']:
            coords.append(np.loadtxt(fname, skiprows=1))
        # Find local index of leading edge point
        le_idx = []
        for c in coords:
            le_idx.append(np.argmin(c[:, 0]))

        # Transform airfoil coordinates using chord length, incidence and offset
        for i in range(n_airf):
            # Reverse order of coordinates if not in standard Selig format
            if coords[i][int(coords[i].shape[0] / 4), 1] < coords[i][-int(coords[i].shape[0] / 4), 1]:
                coords[i] = np.flipud(coords[i])
            # Check that first point is the trailing edge
            if coords[i][0, 0] != 1.0:
                if coords[i][-1, 0] == 1.0:
                    coords[i] = np.roll(coords[i], 1, axis=0) # if last point is the trailing edge instead, insert it as first
                else:
                    raise RuntimeError(str(self.__class__) + ' airfoil coordinates must be ordered using Selig format!\n')
            # Delete last point if same as first point
            if coords[i][0].all() == coords[i][-1].all():
                coords[i] = np.delete(coords[i], (-1), axis=0)
            # Scale according to chord length
            coords[i] *= cfg['chords'][i]
            # Rotate around quarter-chord according to incidence angle
            cos = np.cos(cfg['incidences'][i] * np.pi / 180)
            sin = np.sin(cfg['incidences'][i] * np.pi / 180)
            coords[i] = np.dot(coords[i], np.array([[cos, -sin],[sin, cos]]))
            # Add x and z offset
            coords[i][:, 0] += cfg['le_coords'][i][0] + cfg['le_offset'][0]
            coords[i][:, 1] += cfg['le_coords'][i][2] + cfg['le_offset'][1]
            # Insert y-coordinates 
            coords[i] = np.hstack((coords[i], cfg['le_coords'][i][1] * np.ones([coords[i].shape[0], 1])))
            coords[i][:, [1, 2]] = np.fliplr(coords[i][:, [1, 2]])

        # Get the trailing edge z_coordinate of the airfoil on the symmetry plane
        self.height = coords[0][0, 0]

        return n_airf, coords, le_idx

    def __create_model(self, n_airf, coords, le_idx, name, cfg, domain_cfg, mesh_cfg):
        """Create wing points, curves and surfaces in Gmsh model
        """
        # Add airfoil points
        ptags = []
        for c in coords:
            ptag = []
            for row in c:
                ptag.append(gmsh.model.geo.add_point(row[0], row[1], row[2]))
            ptags.append(ptag)

        # Add airfoils splines
        airf_ctags = []
        for i in range(n_airf):
            ctag_up = gmsh.model.geo.add_spline(ptags[i][0:le_idx[i]+1])
            ctag_lw = gmsh.model.geo.add_spline(np.append(ptags[i][le_idx[i]:], ptags[i][0]))
            airf_ctags.append([ctag_up, ctag_lw])
        # Add TE and LE lines
        plan_ctags = []
        for i in range(n_airf - 1):
            ctag_te = gmsh.model.geo.add_line(ptags[i][0], ptags[i+1][0])
            ctag_le = gmsh.model.geo.add_line(ptags[i][le_idx[i]], ptags[i+1][le_idx[i]])
            plan_ctags.append([ctag_te, ctag_le])

        # Add planforms
        stags = []
        for i in range(n_airf - 1):
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[i][0], plan_ctags[i][1], -airf_ctags[i+1][0], -plan_ctags[i][0]])
            stag_up = gmsh.model.geo.add_surface_filling([-cltag])
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[i][1], plan_ctags[i][0], -airf_ctags[i+1][1], -plan_ctags[i][1]])
            stag_lw = gmsh.model.geo.add_surface_filling([-cltag])
            stags.append(stag_up)
            stags.append(stag_lw)
        # Add wingtip
        cltag = gmsh.model.geo.add_curve_loop([airf_ctags[-1][0], airf_ctags[-1][1]])
        stags.append(gmsh.model.geo.add_plane_surface([-cltag]))

        # Add physical groups
        gmsh.model.geo.synchronize()
        if domain_cfg['type'] == 'potential':
            gmsh.model.add_physical_group(1, [tags[0] for tags in plan_ctags], name=name+'Te')
            gmsh.model.add_physical_group(2, stags[0::2], name=name)
            gmsh.model.add_physical_group(2, stags[1:-1:2], name=name+'_')
        else:
            gmsh.model.add_physical_group(2, stags, name=name)

        # Add meshing constraint on TE and LE
        for i in range(n_airf):
            gmsh.model.mesh.set_size([(0, ptags[i][0])], mesh_cfg['wing_sizes'][name][i][0])
            gmsh.model.mesh.set_size([(0, ptags[i][le_idx[i]])], mesh_cfg['wing_sizes'][name][i][1])

        # Generate tag dictionary
        self.tags = {'airfoil_curves': airf_ctags, 'wing_surfaces': stags}

        # Create wake if requested
        if domain_cfg['type'] == 'potential':
            pte_ids = {} # TE point
            for i in range(n_airf):
                pte_ids[ptags[i][0]] = coords[i][0]
            cte_ids = [] # TE curves
            for i in range(n_airf - 1):
                cte_ids.append(plan_ctags[i][0])
            self.wake = Wake(pte_ids, cte_ids, name+'Wake', domain_cfg['length'], mesh_cfg['domain_size'])
