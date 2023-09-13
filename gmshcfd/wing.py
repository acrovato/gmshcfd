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
from .errors import GmshCFDError
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
        is_closed, n_airf, coords, le_idx = self.__compute_coordinates(cfg, domain_cfg)
        # Create Gmsh wing and wake models
        if is_closed:
            self.__create_model_closed(n_airf, coords, le_idx, name, domain_cfg, mesh_cfg)
        else:
            self.__create_model_open(n_airf, coords, le_idx, name, domain_cfg, mesh_cfg)

    def __compute_coordinates(self, cfg, domain_cfg):
        """Read airfoil coordinates from files and transform them using wing planform configuration
        """
        # Get number of airfoils
        n_airf = len(cfg['airfoils'])

        # Read coordinates
        coords = []
        for fname in cfg['airfoils']:
            coords.append(np.loadtxt(fname, skiprows=1))
        # Find local index of leading edge point
        le_idx = []
        for c in coords:
            le_idx.append(np.argmin(c[:, 0]))

        # Transform airfoil coordinates using chord length, incidence and offset
        for i in range(n_airf):
            # Check if TE is duplicated and if x-coordinate is 1
            if coords[i][0][0] != 1.0 or coords[i][-1][0] != 1.0:
                raise GmshCFDError('airfoil coordinates must be ordered using Selig format: TE point must be first and last, duplicated, and its x-coordinate must be equal to 1.0!\n', self)
            # Reverse order of coordinates if not in standard Selig format
            if coords[i][1, 1] < coords[i][-2, 1]:
                coords[i] = np.flipud(coords[i])
            # Check if airfoil is closed or open
            if coords[i][0][1] == coords[i][-1][1]:
                is_closed = True
                if domain_cfg['type'] == 'rans':
                    raise GmshCFDError('sharp trailing edge airfoils cannot be used with extruded boundary layers (rans domain type)!\n', self)
                coords[i] = np.delete(coords[i], (-1), axis=0) # delete duplicated last point
            else:
                is_closed = False
                if domain_cfg['type'] == 'potential':
                    raise GmshCFDError('blunt trailing edge airfoils cannot be used with wakes (potential domain type)!\n', self)
            # Scale according to chord length
            coords[i] *= cfg['chords'][i]
            # Rotate around quarter-chord according to incidence angle
            cos = np.cos(cfg['incidences'][i] * np.pi / 180)
            sin = np.sin(cfg['incidences'][i] * np.pi / 180)
            coords[i] = np.dot(coords[i], np.array([[cos, -sin],[sin, cos]]))
            # Add x and z offset
            coords[i][:, 0] += cfg['le_offsets'][i][0] + cfg['offset'][0]
            coords[i][:, 1] += cfg['le_offsets'][i][2] + cfg['offset'][1]
            # Insert y-coordinates
            coords[i] = np.hstack((coords[i], cfg['le_offsets'][i][1] * np.ones([coords[i].shape[0], 1])))
            coords[i][:, [1, 2]] = np.fliplr(coords[i][:, [1, 2]])

        # Get the trailing edge z_coordinate of the airfoil on the symmetry plane
        self.height = coords[0][0, 0]

        return is_closed, n_airf, coords, le_idx

    def __create_model_closed(self, n_airf, coords, le_idx, name, domain_cfg, mesh_cfg):
        """Create points, curves and surfaces for wing with sharp trailing edge in Gmsh model
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
            stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[i][1], plan_ctags[i][0], -airf_ctags[i+1][1], -plan_ctags[i][1]])
            stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
        # Add cutoff wingtip
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
        self.tags = {'symmetry_curves': airf_ctags[0], 'wing_surfaces': stags}

        # Create wake if requested
        if domain_cfg['type'] == 'potential':
            pte_ids = {} # TE point
            for i in range(n_airf):
                pte_ids[ptags[i][0]] = coords[i][0]
            cte_ids = [] # TE curves
            for i in range(n_airf - 1):
                cte_ids.append(plan_ctags[i][0])
            self.wake = Wake(pte_ids, cte_ids, name+'Wake', domain_cfg['length'], mesh_cfg['domain_size'])

    def __create_model_open(self, n_airf, coords, le_idx, name, domain_cfg, mesh_cfg):
        """Create points, curves and surfaces for wing with blunt trailing edge in Gmsh model
        """
        # Add airfoil points
        ptags = []
        for c in coords:
            ptag = []
            for row in c:
                ptag.append(gmsh.model.geo.add_point(row[0], row[1], row[2]))
            ptags.append(ptag)
        # Add trailing edge and tip points and centers for blunt TE
        if domain_cfg['type'] == 'rans':
            # TE and TE center
            tec_ptags = []
            te_ptags = []
            for i in range(n_airf):
                c_coord = [0.5 * (coords[i][0, 0] + coords[i][-1, 0]), coords[i][0, 1], 0.5 * (coords[i][0, 2] + coords[i][-1, 2])]
                tec_ptags.append(gmsh.model.geo.add_point(c_coord[0], c_coord[1], c_coord[2]))
                te = gmsh.model.geo.copy([(0, ptags[i][0])])
                gmsh.model.geo.rotate(te, c_coord[0], c_coord[1], c_coord[2], 0, 1, 0, np.pi / 2)
                te_ptags.append(te[0][1])
            # Tip
            tip_ptags = []
            for i in range(1, 10):
                idx = le_idx[-1] * i // 10
                xc = 0.5 * (coords[-1][idx, 0] + coords[-1][-idx-1, 0])
                zc = 0.5 * (coords[-1][idx, 2] + coords[-1][-idx-1, 2])
                yc = coords[-1][idx, 1] + 0.5 * abs(coords[-1][idx, 2] - coords[-1][-idx-1, 2])
                tip_ptags.append(gmsh.model.geo.add_point(xc, yc, zc))

        # Add airfoils splines
        airf_ctags = []
        for i in range(n_airf):
            ctag_up = gmsh.model.geo.add_spline(ptags[i][0:le_idx[i]+1])
            ctag_lw = gmsh.model.geo.add_spline(ptags[i][le_idx[i]:])
            airf_ctags.append([ctag_up, ctag_lw])
        # Add TE and LE lines
        plan_ctags = []
        for i in range(n_airf - 1):
            ctag_teu = gmsh.model.geo.add_line(ptags[i][0], ptags[i+1][0])
            ctag_le = gmsh.model.geo.add_line(ptags[i][le_idx[i]], ptags[i+1][le_idx[i]])
            ctag_tel = gmsh.model.geo.add_line(ptags[i][-1], ptags[i+1][-1])
            plan_ctags.append([ctag_teu, ctag_le, ctag_tel])
        # Add trailing edge curves and tip spline
        if domain_cfg['type'] == 'rans':
            # TE airfoil circles
            airf_te_ctags = []
            for i in range(n_airf):
                ctag_teu = gmsh.model.geo.add_circle_arc(te_ptags[i], tec_ptags[i], ptags[i][0])
                ctag_tel = gmsh.model.geo.add_circle_arc(ptags[i][-1], tec_ptags[i], te_ptags[i])
                airf_te_ctags.append([ctag_teu, ctag_tel])
            # TE planform line
            te_ctags = []
            for i in range(n_airf-1):
                te_ctags.append(gmsh.model.geo.add_line(te_ptags[i], te_ptags[i+1]))
            # Tip
            tip_ctag = gmsh.model.geo.add_spline([te_ptags[-1]] + tip_ptags)
            tip_le_ctag = gmsh.model.geo.add_line(tip_ptags[-1], ptags[-1][le_idx[-1]])
        else:
            # TE line
            airf_te_ctags = []
            for i in range(n_airf):
                airf_te_ctags.append(gmsh.model.geo.add_line(ptags[i][-1], ptags[i][0]))

        # Add planforms
        stags = []
        for i in range(n_airf - 1):
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[i][0], plan_ctags[i][1], -airf_ctags[i+1][0], -plan_ctags[i][0]])
            stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[i][1], plan_ctags[i][2], -airf_ctags[i+1][1], -plan_ctags[i][1]])
            stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
        # Add TE and wingtip
        if domain_cfg['type'] == 'rans':
            # Rounded TE
            for i in range(n_airf - 1):
                cltag = gmsh.model.geo.add_curve_loop([airf_te_ctags[i][0], plan_ctags[i][0], -airf_te_ctags[i+1][0], -te_ctags[i]])
                stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
                cltag = gmsh.model.geo.add_curve_loop([airf_te_ctags[i][1], te_ctags[i], -airf_te_ctags[i+1][1], -plan_ctags[i][2]])
                stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
            # Rounded wingtip
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[-1][0], -tip_le_ctag, -tip_ctag, airf_te_ctags[-1][0]])
            stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
            cltag = gmsh.model.geo.add_curve_loop([airf_te_ctags[-1][1], tip_ctag, tip_le_ctag, airf_ctags[-1][1]])
            stags.append(gmsh.model.geo.add_surface_filling([-cltag]))
        else:
            # Cutoff TE
            for i in range(n_airf - 1):
                cltag = gmsh.model.geo.add_curve_loop([airf_te_ctags[i], plan_ctags[i][0], -airf_te_ctags[i+1], -plan_ctags[i][2]])
                stags.append(gmsh.model.geo.add_plane_surface([-cltag]))
            # Cutoff wingtip
            cltag = gmsh.model.geo.add_curve_loop([airf_ctags[-1][0], airf_ctags[-1][1], airf_te_ctags[-1]])
            stags.append(gmsh.model.geo.add_plane_surface([-cltag]))

        # Add boundary layer
        if domain_cfg['type'] == 'rans':
            # Extrude surfaces to create boundary layer
            N = mesh_cfg['boundary_layer']['wing']['num_layer'] # number of layers
            r = mesh_cfg['boundary_layer']['wing']['growth_ratio'] # ratio
            d = [mesh_cfg['boundary_layer']['wing']['hgt_first_layer']] # heigth of first layer
            for i in range(1, N):
                d.append(d[-1] + d[0] * r**i)
            bl = gmsh.model.geo.extrudeBoundaryLayer([(2, tag) for tag in stags], [1] * N, d, True)
            gmsh.model.geo.synchronize()
            # Get tags of top and symmetry surfaces, and volume
            bl_top_stags = []
            bl_sym_stags = []
            bl_vtags = []
            for i in range(len(bl)):
                if bl[i][0] == 3:
                    bl_top_stags.append(bl[i-1][1])
                    bl_vtags.append(bl[i][1])
                    # Get surfaces on symmetry plane (originate from 4 extruded surfaces, each giving 5 surfaces + 1 volume)
                    if i < 12:
                        bl_sym_stags.append(bl[i + 4])
                        bl_sym_stags.append(bl[2 * (n_airf - 1) * 6 + i + 4])
            # Get the bounding curves on symmetry plane, gives 2 closed countours (inner and outer), so remove inner
            bl_sym_ctags = [c[1] for c in gmsh.model.getBoundary(bl_sym_stags)]
            for i in range(2):
                bl_sym_ctags.remove(-airf_ctags[0][i])
                bl_sym_ctags.remove(-airf_te_ctags[0][i])

        # Add physical groups
        gmsh.model.geo.synchronize()
        gmsh.model.add_physical_group(2, stags, name=name)
        if domain_cfg['type'] == 'rans' and mesh_cfg['boundary_layer']['write_tags']:
            gmsh.model.add_physical_group(2, bl_top_stags, tag=9998, name=name+'BoundaryLayerSurface')
            gmsh.model.add_physical_group(3, bl_vtags, tag=9999, name=name+'BoundaryLayerVolume')

        # Add meshing constraint on TE and LE
        for i in range(n_airf):
            gmsh.model.mesh.set_size([(0, ptags[i][0]), (0, ptags[i][-1])], mesh_cfg['wing_sizes'][name][i][0])
            gmsh.model.mesh.set_size([(0, ptags[i][le_idx[i]])], mesh_cfg['wing_sizes'][name][i][1])
        if domain_cfg['type'] == 'rans':
            for i in range(n_airf):
                gmsh.model.mesh.set_size([(0, te_ptags[i])], mesh_cfg['wing_sizes'][name][i][0])
            gmsh.model.mesh.set_size([(0, tip_ptags[-1])], mesh_cfg['wing_sizes'][name][i][0])

        # Generate tag dictionary
        if domain_cfg['type'] == 'rans':
            self.tags = {'symmetry_curves': bl_sym_ctags,
                         'boundary_layer_top_surfaces': bl_top_stags,
                         'boundary_layer_symmetry_surfaces': [s[1] for s in bl_sym_stags],
                         'boundary_layer_volume': bl_vtags}
        else:
            self.tags = {'symmetry_curves': airf_ctags[0] + [airf_te_ctags[0]], 'wing_surfaces': stags}
