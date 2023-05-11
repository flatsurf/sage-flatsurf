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

from flatsurf.geometry.surface import OrientedSimilaritySurface
from flatsurf.geometry.mappings import SurfaceMapping

from flatsurf.geometry.polygon import ConvexPolygons
from sage.misc.cachefunc import cached_method


class GL2RImageSurface(OrientedSimilaritySurface):
    r"""
    This is a lazy implementation of the SL(2,R) image of a translation surface.

    EXAMPLE::

        sage: import flatsurf
        sage: s=flatsurf.translation_surfaces.octagon_and_squares()
        sage: r=matrix(ZZ,[[0,1],[1,0]])
        sage: ss=r*s
        sage: TestSuite(ss).run()
        sage: s.canonicalize()==ss.canonicalize()
        True

    """

    def __init__(self, surface, m, ring=None, category=None):
        if surface.is_mutable():
            if surface.is_finite_type():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                self._s = MutableOrientedSimilaritySurface.from_surface(surface)
            else:
                raise ValueError("Can not apply matrix to mutable infinite surface.")
        else:
            self._s = surface

        det = m.determinant()

        if det > 0:
            self._det_sign = 1
        elif det < 0:
            self._det_sign = -1
        else:
            raise ValueError("Can not apply matrix with zero determinant to surface.")

        if m.is_mutable():
            from sage.all import matrix
            m = matrix(m)
            m.set_immutable()
        self._m = m

        if ring is None:
            if m.base_ring() == self._s.base_ring():
                base_ring = self._s.base_ring()
            else:
                from sage.structure.element import get_coercion_model

                cm = get_coercion_model()
                base_ring = cm.common_parent(m.base_ring(), self._s.base_ring())
        else:
            base_ring = ring

        self._P = ConvexPolygons(base_ring)

        if category is None:
            category = surface.category()

        super().__init__(base_ring, category=category)

    def roots(self):
        return self._s.roots()

    def is_compact(self):
        return self._s.is_compact()

    def is_mutable(self):
        return False

    def is_translation_surface(self, positive=True):
        return self._s.is_translation_surface(positive=positive)

    @cached_method
    def polygon(self, lab):
        if self._det_sign == 1:
            p = self._s.polygon(lab)
            edges = [self._m * p.edge(e) for e in range(p.num_edges())]
            return self._P(edges)
        else:
            p = self._s.polygon(lab)
            edges = [self._m * (-p.edge(e)) for e in range(p.num_edges() - 1, -1, -1)]
            return self._P(edges)

    def labels(self):
        return self._s.labels()

    def opposite_edge(self, p, e):
        if self._det_sign == 1:
            return self._s.opposite_edge(p, e)
        else:
            polygon = self._s.polygon(p)
            pp, ee = self._s.opposite_edge(p, polygon.num_edges() - 1 - e)
            polygon2 = self._s.polygon(pp)
            return pp, polygon2.num_edges() - 1 - ee

    def __repr__(self):
        if self.is_finite_type():
            from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
            return repr(MutableOrientedSimilaritySurface.from_surface(self))

        return f"GL2RImageSurface of {self._s}"

    def __hash__(self):
        return hash((self._s, self._m))

    def __eq__(self, other):
        r"""
        Return whether this image is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of inequality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: m = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = m * S
            sage: S == S
            True

        """
        if not isinstance(other, GL2RImageSurface):
            return False

        return self._s == other._s and self._m == other._m and self.base_ring() == other.base_ring()


class GL2RMapping(SurfaceMapping):
    r"""
    This class pushes a surface forward under a matrix.

    Note that for matrices of negative determinant we need to relabel edges (because
    edges must have a counterclockwise cyclic order). For each n-gon in the surface,
    we relabel edges according to the involution e mapsto n-1-e.

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
