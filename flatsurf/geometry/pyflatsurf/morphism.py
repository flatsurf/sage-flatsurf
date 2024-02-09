# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2024 Julian Rüth
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
from flatsurf.geometry.morphism import SurfaceMorphism


class Morphism_to_pyflatsurf(SurfaceMorphism):
    def __init__(self, parent, pyflatsurf_conversion):
        super().__init__(parent)
        self._pyflatsurf_conversion = pyflatsurf_conversion

    def _image_edge(self, label, edge):
        half_edge = self._pyflatsurf_conversion((label, edge))
        face = tuple(self.codomain()._flat_triangulation.face(half_edge))
        label = type(self.codomain())._normalize_label(face)
        edge = label.index(half_edge)
        return (label, edge)

    def _image_homology_edge(self, label, edge):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: deformation = S.pyflatsurf()
            sage: deformation._image_homology_edge((0, 0), 0)
            [(1, (1, 2, 3), 0)]

        """
        return [(1, *self._image_edge(label, edge))]

    def _image_saddle_connection(self, connection):
        from flatsurf.geometry.pyflatsurf.saddle_connection import SaddleConnection_pyflatsurf
        return SaddleConnection_pyflatsurf(self._pyflatsurf_conversion(connection), self.codomain())

    def section(self):
        return Morphism_from_pyflatsurf._create_morphism(self.codomain(), self.domain(), self._pyflatsurf_conversion)

    def _image_point(self, point):
        from flatsurf.geometry.pyflatsurf.surface_point import SurfacePoint_pyflatsurf
        # TODO: Category is unset.
        return SurfacePoint_pyflatsurf(self._pyflatsurf_conversion(point), self.codomain())

    # TODO: Implement __eq__


class Morphism_from_pyflatsurf(SurfaceMorphism):
    def __init__(self, parent, pyflatsurf_conversion):
        super().__init__(parent)
        self._pyflatsurf_conversion = pyflatsurf_conversion

    def _image_half_edge(self, half_edge):
        face = tuple(self.domain()._flat_triangulation.face(half_edge))
        label = type(self.domain())._normalize_label(face)
        edge = label.index(half_edge)
        return (label, edge)

    def _image_point(self, point):
        from flatsurf.geometry.pyflatsurf.surface_point import SurfacePoint_pyflatsurf
        assert isinstance(point, SurfacePoint_pyflatsurf)

        return self._pyflatsurf_conversion.section(point._point)

    def _image_saddle_connection(self, connection):
        return self._pyflatsurf_conversion._preimage_saddle_connection(connection._connection)

    def _image_homology_edge(self, label, edge):
        half_edge = label[edge]

        import pyflatsurf
        half_edge = pyflatsurf.flatsurf.HalfEdge(int(half_edge))

        return [(1, *self._pyflatsurf_conversion._preimage_half_edge(half_edge))]

    # TODO: Implement __eq__


class Morphism_from_Deformation(SurfaceMorphism):
    def __init__(self, parent, deformation):
        super().__init__(parent)
        self._deformation = deformation

    def _image_homology_edge(self, label, edge):
        # TODO: We should probably mark this method as only being correct in homology.
        half_edge = label[edge]

        from pyflatsurf import flatsurf
        saddle_connection = flatsurf.SaddleConnection[type(self.domain()._flat_triangulation)](self.domain()._flat_triangulation, flatsurf.HalfEdge(half_edge))
        path = flatsurf.Path[type(self.domain()._flat_triangulation)](saddle_connection)

        path = self._deformation(path)

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
                    image.append((1, label, edge))

        return image

    def _image_point(self, point):
        from flatsurf.geometry.pyflatsurf.surface_point import SurfacePoint_pyflatsurf
        assert isinstance(point, SurfacePoint_pyflatsurf)

        image = self._deformation(point._point)

        # TODO: Category is unset.
        return SurfacePoint_pyflatsurf(image, self.codomain())

    # TODO: Implement __eq__
