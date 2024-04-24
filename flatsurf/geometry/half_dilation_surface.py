# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
#                          2023 Julian Rüth
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
# ****************************************************************************

from flatsurf.geometry.mappings import SurfaceMapping


class GL2RMapping(SurfaceMapping):
    r"""
    This class pushes a surface forward under a matrix.

    Note that for matrices of negative determinant we need to relabel edges (because
    edges must have a counterclockwise cyclic order). For each n-gon in the surface,
    we relabel edges according to the involution `e \mapsto n-1-e`.

    EXAMPLE::

        sage: from flatsurf import translation_surfaces
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.half_dilation_surface import GL2RMapping
        sage: mat=Matrix([[2,1],[1,1]])
        sage: m=GL2RMapping(s,mat)
        sage: TestSuite(m.codomain()).run()
    """

    def __init__(self, s, m, ring=None, category=None):
        r"""
        Hit the surface s with the 2x2 matrix m which should have positive determinant.
        """
        codomain = GL2RImageSurface(s, m, ring=ring, category=category or s.category())
        self._m = m
        self._im = ~m
        SurfaceMapping.__init__(self, s, codomain)

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        return self.codomain().tangent_vector(
            tangent_vector.polygon_label(),
            self._m * tangent_vector.point(),
            self._m * tangent_vector.vector(),
        )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        return self.domain().tangent_vector(
            tangent_vector.polygon_label(),
            self._im * tangent_vector.point(),
            self._im * tangent_vector.vector(),
        )
