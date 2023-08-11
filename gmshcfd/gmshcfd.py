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

from .wing import Wing
from .domain import Box, Sphere
import gmsh

# TODOLIST:
# Check installed version + logger
# Add full support for blunt TE
# Add support for BL
# Refactoring: check/split method? call from driver instead of constructor?

class GmshCFD:
    """Main driver

    Parameters:
    name: string
        name of the model
    cfg: dict
        geometrical and mesh parameters

    Attributes:
    name: string
        name of the model
    wing_cfgs: dict
        geometrical parameters defining the lifting surfaces
    domain_cfg: dict
        geometrical parameters defining the domain
    mesh_cfg: dict
        parameters defining the mesh
    wings: list
        list of lifting surfaces
    domain: gmshcfd.Box or Sphere object
        domain
    """
    def __init__(self, name, cfg):
        # Initialize attributes
        self.__name = name
        self.__wing_cfgs = cfg['wings']
        self.__domain_cfg = cfg['domain']
        self.__mesh_cfg = cfg['mesh']
        self.__wings = []
        self.__domain = None
        # Start Gmsh logger
        gmsh.initialize()
        gmsh.logger.start()
        gmsh.model.add(name)

    def __del__(self):
        # Get log and stop Gmsh
        log_msgs = gmsh.logger.get()
        gmsh.logger.stop()
        file = open('log', 'a')
        for m in log_msgs:
            file.write(m + '\n')
        file.close()
        gmsh.finalize()

    def generate_geometry(self):
        """Generate the wings and the domain using the configurations
        """
        # Create wings
        for name, cfg in self.__wing_cfgs.items():
            self.__wings.append(Wing(name, cfg, self.__domain_cfg, self.__mesh_cfg))
        # Create domain
        if self.__domain_cfg['type'] == 'potential':
            self.__domain = Box(self.__wings, self.__domain_cfg, self.__mesh_cfg)
        else:
            self.__domain = Sphere(self.__wings, self.__domain_cfg, self.__mesh_cfg)
        # Synchronize model
        gmsh.model.geo.synchronize()

    def generate_mesh(self):
        """Generate mesh
        """
        gmsh.option.set_number('Mesh.Algorithm', 5) # Delaunay
        gmsh.option.set_number('Mesh.Algorithm3D', 1) # Delaunay
        gmsh.option.set_number('Mesh.Optimize', 1)
        gmsh.option.set_number('Mesh.Smoothing', 10)
        gmsh.option.set_number('Mesh.SmoothNormals', 1)
        gmsh.model.mesh.generate(3)

    def write_geometry(self):
        """Save geometry to disk and rename using .geo
        """
        import os
        gmsh.write(self.__name + '.geo_unrolled')
        nname = self.__name + '.geo'
        if os.path.isfile(nname):
            os.remove(nname)
        os.rename(self.__name + '.geo_unrolled', nname)

    def write_mesh(self, format):
        """Save mesh to disk
        """
        if format == 'msh2':
            gmsh.option.set_number('Mesh.MshFileVersion', 2.2)
            gmsh.write(self.__name + '.msh')
        else:
            gmsh.write(self.__name + '.' + format)
