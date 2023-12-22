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
import gmsh, os

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
        file = open(f'log_{self.__name}', 'w')
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

    def generate_mesh(self, algo_2d='delaunay', algo_3d='hxt'):
        """Generate mesh
        """
        algos_2d = {'delaunay': 5, 'frontal-delaunay': 6}
        algos_3d = {'delaunay': 1, 'hxt': 10}
        gmsh.option.set_number('Mesh.Algorithm', algos_2d[algo_2d])
        gmsh.option.set_number('Mesh.Algorithm3D', algos_3d[algo_3d])
        gmsh.option.set_number('Mesh.Optimize', 1)
        gmsh.option.set_number('Mesh.Smoothing', 10)
        gmsh.option.set_number('Mesh.SmoothNormals', 1)
        gmsh.option.set_number('General.NumThreads', os.cpu_count())
        try:
            gmsh.model.mesh.generate(3)
        except Exception as e:
            gmsh.write(self.__name + '.msh')
            raise Exception(e)

    def write_geometry(self):
        """Save geometry to disk and rename using .geo
        """
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
