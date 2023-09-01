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

import gmsh

class Box:
    """Box-shaped domain

    Parameters:
    wings: list
        list of lifting surfaces
    cfg: dict
        parameters to configure domain
    mesh_cfg: dict
        parameters to configure mesh
    """
    def __init__(self, wings, cfg, mesh_cfg):
        # Create Gmsh box domain model
        wing_order = self.__order_wakes(wings)
        self.__create_model(wings, wing_order, cfg, mesh_cfg)

    def __order_wakes(self, wings):
        """Order the lifting surfaces according to the z-coordinate of their trailing edge on the symmetry plane
        """
        idx = []
        zc = []
        for i in range(len(wings)):
            idx.append(i)
            zc.append(wings[i].height)
        return [x for _, x in sorted(zip(zc, idx))]

    def __create_model(self, wings, wing_order, cfg, mesh_cfg):
        """Create domain points, curves, surfaces and volume in Gmsh model
        """
        # Add points
        ptags = []
        for i in range(2):
            ptag = []
            ptag.append(gmsh.model.geo.add_point(cfg['length'], i*cfg['length'], cfg['length']))
            ptag.append(gmsh.model.geo.add_point(-cfg['length'], i*cfg['length'], cfg['length']))
            ptag.append(gmsh.model.geo.add_point(-cfg['length'], i*cfg['length'], -cfg['length']))
            ptag.append(gmsh.model.geo.add_point(cfg['length'], i*cfg['length'], -cfg['length']))
            ptags.append(ptag)

        # Add curves
        oxz_ctags = []
        # oxz plane, symmetry side
        ctag = []
        for i in range(3):
            ctag.append(gmsh.model.geo.add_line(ptags[0][i], ptags[0][i+1]))
        ctag.append(gmsh.model.geo.add_line(ptags[0][3], wings[wing_order[0]].wake.tags['symmetry_point']))
        for i in range(len(wing_order) - 1):
            ctag.append(gmsh.model.geo.add_line(wings[wing_order[i]].wake.tags['symmetry_point'], wings[wing_order[i+1]].wake.tags['symmetry_point']))
        ctag.append(gmsh.model.geo.add_line(wings[wing_order[-1]].wake.tags['symmetry_point'], ptags[0][0]))
        oxz_ctags.append(ctag)
        # oxz plane, back side
        ctag = []
        for i in range(4):
            ctag.append(gmsh.model.geo.add_line(ptags[1][i], ptags[1][(i+1) % 4]))
        oxz_ctags.append(ctag)
        # oxy plane
        oyz_ctags = []
        for i in range(4):
            oyz_ctags.append(gmsh.model.geo.add_line(ptags[0][i], ptags[1][i]))

        # Add surfaces
        stags = []
        # symmetry
        cltag_wings = []
        for i in range(len(wings)):
            cltag_wings.append(gmsh.model.geo.add_curve_loop(wings[i].tags['symmetry_curves']))
        cltag = gmsh.model.geo.add_curve_loop(list(reversed(oxz_ctags[0])))
        stags.append(gmsh.model.geo.add_plane_surface([-cltag] + cltag_wings))
        # upstream
        cltag = gmsh.model.geo.add_curve_loop([oxz_ctags[0][1], oyz_ctags[2], -oxz_ctags[1][1], -oyz_ctags[1]])
        stags.append(gmsh.model.geo.add_plane_surface([cltag]))
        # downstream
        cltag = gmsh.model.geo.add_curve_loop([oyz_ctags[0], -oxz_ctags[1][3], -oyz_ctags[3]] + oxz_ctags[0][3:])
        stags.append(gmsh.model.geo.add_plane_surface([cltag]))
        # farfield
        cltag = gmsh.model.geo.add_curve_loop([oxz_ctags[0][0], oyz_ctags[1], -oxz_ctags[1][0], -oyz_ctags[0]])
        stags.append(gmsh.model.geo.add_plane_surface([cltag]))
        cltag = gmsh.model.geo.add_curve_loop([oyz_ctags[3], -oxz_ctags[1][2], -oyz_ctags[2], oxz_ctags[0][2]])
        stags.append(gmsh.model.geo.add_plane_surface([cltag]))
        cltag = gmsh.model.geo.add_curve_loop(oxz_ctags[1])
        stags.append(gmsh.model.geo.add_plane_surface([cltag]))

        # Add volume
        wing_tags = []
        for w in wings:
            wing_tags.extend(w.tags['wing_surfaces'])
        sltag = gmsh.model.geo.add_surface_loop(stags + wing_tags)
        vtag = gmsh.model.geo.add_volume([sltag])

        # Add embbeded entities
        symmetry_tags = []
        trailing_tags = []
        wake_tags = []
        for w in wings:
            symmetry_tags.append(w.wake.tags['symmetry_curve'])
            trailing_tags.extend(w.wake.tags['trailing_curves'])
            wake_tags.extend(w.wake.tags['wake_surfaces'])
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(1, symmetry_tags, 2, stags[0])
        gmsh.model.mesh.embed(1, trailing_tags, 2, stags[2])
        gmsh.model.mesh.embed(2, wake_tags, 3, vtag)

        # Add physical groups
        gmsh.model.geo.synchronize()
        gmsh.model.add_physical_group(2, [stags[0]], name='symmetry')
        gmsh.model.add_physical_group(2, [stags[1]], name='upstream')
        gmsh.model.add_physical_group(2, [stags[2]], name='downstream')
        gmsh.model.add_physical_group(2, stags[3:], name='farfield')
        gmsh.model.add_physical_group(3, [vtag], name='field')

        # Add meshing constraints
        for i in range(2):
            for j in range(len(ptags[i])):
                gmsh.model.mesh.set_size([(0, ptags[i][j])], mesh_cfg['domain_size'])

class Sphere:
    """Sphere-shaped domain

    Parameters:
    wings: list
        list of lifting surfaces
    cfg: dict
        parameters to configure domain
    mesh_cfg: dict
        parameters to configure mesh
    """
    def __init__(self, wings, cfg, mesh_cfg):
        # Create Gmsh sphere domain model
        self.__create_model(wings, cfg, mesh_cfg)

    def __create_model(self, wings, cfg, mesh_cfg):
        """Create domain points, curves, surfaces and volume in Gmsh model
        """
        # Add points
        c_ptag = gmsh.model.geo.add_point(0., 0., 0.) # center of the sphere
        ptags = []
        ptags.append(gmsh.model.geo.add_point(cfg['length'], 0., 0.))
        ptags.append(gmsh.model.geo.add_point(0., 0., cfg['length']))
        ptags.append(gmsh.model.geo.add_point(-cfg['length'], 0., 0.))
        ptags.append(gmsh.model.geo.add_point(0., 0., -cfg['length']))
        ptags.append(gmsh.model.geo.add_point(0., cfg['length'], 0.))

        # Add curves
        ctags = []
        # symmetry plane
        ctag = []
        for i in range(4):
            ctag.append(gmsh.model.geo.add_circle_arc(ptags[i], c_ptag, ptags[(i+1) % 4]))
        ctags.append(ctag)
        # outside symmetry plane
        ctag = []
        for i in range(4):
            ctag.append(gmsh.model.geo.add_circle_arc(ptags[i], c_ptag, ptags[-1]))
        ctags.append(ctag)

        # Add surfaces
        stags = []
        # symmetry
        cltag_wings = []
        for i in range(len(wings)):
            cltag_wings.append(gmsh.model.geo.add_curve_loop(wings[i].tags['symmetry_curves']))
        cltag = gmsh.model.geo.add_curve_loop(list(reversed(ctags[0])))
        stags.append(gmsh.model.geo.add_plane_surface([-cltag] + cltag_wings))
        # farfield
        for i in range(4):
            cltag = gmsh.model.geo.add_curve_loop([ctags[0][i], ctags[1][(i+1) % 4], -ctags[1][i]])
            stags.append(gmsh.model.geo.add_surface_filling([cltag]))

        # Add volume
        wing_tags = []
        tag_name = 'boundary_layer_top' if cfg['type'] == 'rans' else 'wing'
        for w in wings:
            wing_tags.extend(w.tags[tag_name + '_surfaces'])
        sltag = gmsh.model.geo.add_surface_loop(stags + wing_tags)
        vtag = gmsh.model.geo.add_volume([sltag])

        # Add physical groups
        gmsh.model.geo.synchronize()
        sym_stags = [stags[0]]
        if cfg['type'] == 'rans':
            for w in wings:
                sym_stags.extend(w.tags['boundary_layer_symmetry_surfaces'])
        vtags = [vtag]
        if cfg['type'] == 'rans':
            for w in wings:
                vtags.extend(w.tags['boundary_layer_volume'])
        gmsh.model.add_physical_group(2, sym_stags, name='symmetry')
        gmsh.model.add_physical_group(2, stags[1:], name='farfield')
        gmsh.model.add_physical_group(3, vtags, name='field')

        # Add meshing constraints
        for i in range(len(ptags)):
            gmsh.model.mesh.set_size([(0, ptags[i])], mesh_cfg['domain_size'])
