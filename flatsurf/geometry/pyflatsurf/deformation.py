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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: deformation = S.pyflatsurf()
            sage: deformation._image_edge(0, 0)
            [((1, 2, 3), 0)]

        """
        half_edge = self._pyflatsurf_conversion((label, edge))
        face = tuple(self.codomain()._flat_triangulation.face(half_edge))
        label = type(self.codomain())._normalize_label(face)
        edge = label.index(half_edge)
        return [(label, edge)]


class Deformation_from_pyflatsurf(Deformation):
    pass


class Deformation_pyflatsurf(Deformation):
    def __init__(self, domain, codomain, pyflatsurf_deformation):
        super().__init__(domain, codomain)
