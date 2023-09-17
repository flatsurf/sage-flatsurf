r"""
The category of polygons in the hyperbolic plane.

EXAMPLES::

    sage: from flatsurf import HyperbolicPlane
    sage: H = HyperbolicPlane()

    sage: P = H.polygon([
    ....:   H.vertical(1).left_half_space(),
    ....:   H.vertical(-1).right_half_space(),
    ....:   H.half_circle(0, 2).left_half_space(),
    ....:   H.half_circle(0, 4).right_half_space(),
    ....: ])

    sage: P.category()
    Category of facade convex simple hyperbolic polygons over Rational Field

    sage: from flatsurf.geometry.categories import HyperbolicPolygons
    sage: P in HyperbolicPolygons(QQ)
    True

"""
# ****************************************************************************
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
# ****************************************************************************
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

from flatsurf.geometry.categories.polygons import Polygons


class HyperbolicPolygons(Category_over_base_ring):
    r"""
    The category of polygons in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import HyperbolicPolygons
        sage: C = HyperbolicPolygons(QQ)

    TESTS::

        sage: TestSuite(C).run()

    """

    def super_categories(self):
        r"""
        Return the categories that a hyperbolic polygon is also a member of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import HyperbolicPolygons
            sage: C = HyperbolicPolygons(QQ)
            sage: C.super_categories()
            [Category of polygons over Rational Field]

        """
        return [Polygons(self.base_ring())]

    class ParentMethods:
        def is_equilateral(self):
            r"""
            Return whether all edges of this polygon have the same length.

            ALGORITHM:

            We construct isometries between the edges. If no such isometries
            exist, the sides are not of the same length.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: P = H.convex_hull(H(I), H(I + 1), H(2*I), H(I - 1))
                sage: P.is_equilateral()
                True

            ::

                sage: from flatsurf import HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: P = H.convex_hull(H(I), H(I + 1), H(3*I), H(I - 1))
                sage: P.is_equilateral()
                False

            """
            edges = list(self.edges())
            finite_edges = [edge for edge in edges if edge.is_finite()]

            if not finite_edges:
                return True

            if len(finite_edges) != len(edges):
                return False

            from itertools import combinations
            for e, f in combinations(edges, 2):
                try:
                    e.parent().isometry(e, f)
                except ValueError:
                    # TODO: Is this true? Can such an isometry always be constructed?
                    return False

            return True

        def is_equiangular(self):
            r"""
            Return whether this polygon is the convex hull of its vertices and
            the angle at each of its vertices is the same.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: P = H.convex_hull(H(I), H(2*I + 1), H(2*I - 1))
                sage: P.is_equiangular()
                True

            ::

                sage: P = H.convex_hull(H(I), H(I + 1), H(2*I), H(I - 1))
                sage: P.is_equiangular()
                False

            """
            if len(self.edges()) != len(self.vertices()):
                return False

            from itertools import pairwise

            for a, b in pairwise(self.angles()):
                if a != b:
                    return False

            return True

        def is_right_triangle(self):
            r"""
            Return whether this is a triangle with a π/2 angle.

            ALGORITHM:

            We construct a geodesic through each vertex perpendicular to an
            adjacent edge. If that geodesic is parallel to the other edge, the
            angle is π/2. See
            :meth:`hyperbolic.HyperbolicGeodesic.perpendicular`.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: P = H.intersection(
                ....:   H.vertical(0).left_half_space(),
                ....:   H.vertical(-1).right_half_space(),
                ....:   H.half_circle(0, 1).left_half_space())

                sage: P.is_right_triangle()
                True

            ::

                sage: P = H.intersection(
                ....:   H.vertical(1).left_half_space(),
                ....:   H.vertical(-1).right_half_space(),
                ....:   H.half_circle(0, 2).left_half_space())

                sage: P.is_right_triangle()
                False

            """
            if len(self.vertices()) != 3:
                return False

            if len(self.edges()) != 3:
                return False

            for edge, following_edge in self.edges().pairs():
                vertex = edge.end()

                if not vertex.is_finite():
                    continue

                if edge.geodesic().perpendicular(edge.end()) == following_edge.geodesic().unoriented():
                    return True

            return False

        def is_isosceles_triangle(self):
            r"""
            Return whether this is a triangle with two sides of equal length.

            ALGORITHM:

            We check whether two sides have equal length by searching for an
            isometry between them, see :meth:`HyperbolicPlane.isometry`.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: P = H.convex_hull(H(I), H(2*I + 1), H(2*I - 1))

                sage: P.is_isosceles_triangle()
                True

            ::

                sage: P = H.intersection(
                ....:   H.vertical(1).left_half_space(),
                ....:   H.vertical(-1).right_half_space(),
                ....:   H.half_circle(0, 2).left_half_space())

                sage: P.is_isosceles_triangle()
                False

            """
            if len(self.edges()) != 3:
                return False

            finite_edges = [segment for segment in self.edges() if segment.is_finite()]

            if len(finite_edges) <= 1:
                # When there are at least two infinite segments, then these are
                # of equal length and this triangle is isosceles.
                return True

            from itertools import combinations
            for (e, f) in combinations(finite_edges, 2):
                try:
                    e.parent().isometry(e, f)
                except ValueError:
                    # TODO: Is this true? Could such an isometry always be constructed?
                    continue
                else:
                    return True

            return False

        def get_point_position(self, point):
            point = self.parent()(point)

            from flatsurf.geometry.polygon import PolygonPosition
            for (i, v) in enumerate(self.vertices()):
                if point == v:
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=i)

            for (i, e) in enumerate(self.edges()):
                if point in e:
                    return PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)

            if point in self:
                return PolygonPosition(PolygonPosition.INTERIOR)

            return PolygonPosition(PolygonPosition.OUTSIDE)

    class Simple(CategoryWithAxiom_over_base_ring):
        class Convex(CategoryWithAxiom_over_base_ring):
            class ParentMethods:
                def area(self):
                    # Returned in multiples of π.
                    angles = self.angles()

                    # Gauss Bonnet, 1.4.2 in Katok
                    return len(angles) - sum(angles)
