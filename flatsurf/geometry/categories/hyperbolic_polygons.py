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
            # TODO
            return None

        def is_equiangular(self):
            # TODO
            return None

        def is_right(self):
            # TODO
            return False

        def is_isosceles(self):
            # TODO
            return False

        def get_point_position(self, point):
            point = self.parent()(point)

            from flatsurf.geometry.polygon import PolygonPosition
            for (i, v) in enumerate(self.vertices()):
                if point == v:
                    edges = [e for e, edge in enumerate(self.edges()) if edge.start() == v]
                    assert len(edges) == 1
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=i, edge=edges[0])

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
