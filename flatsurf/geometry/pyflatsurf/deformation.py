# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# ********************************************************************
from flatsurf.geometry.deformation import Deformation


class Deformation_to_pyflatsurf(Deformation):
    def __init__(self, domain, codomain, pyflatsurf_conversion):
        self._pyflatsurf_conversion = pyflatsurf_conversion
        super().__init__(domain, codomain)

    def _image_edge(self, label, edge):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: deformation = S.pyflatsurf()
            sage: deformation._image_edge((0, 0), 0)
            [((1, 2, 3), 0)]

        """
        half_edge = self._pyflatsurf_conversion((label, edge))
        face = tuple(self.codomain()._flat_triangulation.face(half_edge))
        label = type(self.codomain())._normalize_label(face)
        edge = label.index(half_edge)
        return [(label, edge)]

    def _image_saddle_connection(self, connection):
        from flatsurf.geometry.pyflatsurf.saddle_connection import SaddleConnection_pyflatsurf
        return SaddleConnection_pyflatsurf(self._pyflatsurf_conversion(connection))


class Deformation_from_pyflatsurf(Deformation):
    pass


class Deformation_pyflatsurf(Deformation):
    def __init__(self, domain, codomain, pyflatsurf_deformation):
        self._pyflatsurf_deformation = pyflatsurf_deformation
        super().__init__(domain, codomain)

    def _image_edge(self, label, edge):
        # TODO: We should probably mark this method as only being correct in homology.
        half_edge = label[edge]

        from pyflatsurf import flatsurf
        saddle_connection = flatsurf.SaddleConnection[type(self.domain()._flat_triangulation)](self.domain()._flat_triangulation, flatsurf.HalfEdge(half_edge))
        path = flatsurf.Path[type(self.domain()._flat_triangulation)](saddle_connection)

        path = self._pyflatsurf_deformation(path)

        if not path:
            raise NotImplementedError("cannot map edge through this deformation in pyflatsurf yet")

        path = path.value()

        image = []
        for step in path:
            chain = step.chain()
            for edge, coefficient in chain:
                from flatsurf.geometry.pyflatsurf_conversion import RingConversion
                coefficient = RingConversion.from_pyflatsurf_from_elements([coefficient]).section(coefficient)

                half_edge = edge.positive()
                if coefficient < 0:
                    half_edge = -half_edge
                    coefficient *= -1

                face = tuple(self.codomain()._flat_triangulation.face(half_edge))
                label = type(self.codomain())._normalize_label(face)
                edge = label.index(half_edge)

                for repetition in range(coefficient):
                    image.append((label, edge))

        return image
