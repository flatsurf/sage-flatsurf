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
