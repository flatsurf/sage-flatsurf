# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring, all_axioms

from sage.categories.all import Sets


class Polygons(Category_over_base_ring):
    def super_categories(self):
        return [Sets()]

    class ParentMethods:
        def _check(self):
            if self.area() < 0:
                raise ValueError("polygon has negative area; probably the vertices are not in counter-clockwise order?")

    class Convex(CategoryWithAxiom_over_base_ring):
        r"""
        The set of convex polygons with a fixed base field.

        EXAMPLES::

            sage: from flatsurf import ConvexPolygons
            sage: C = ConvexPolygons(QQ)
            sage: C(vertices=[(0,0), (2,0), (1,1)])
            polygon(vertices=[(0, 0), (2, 0), (1, 1)])
            sage: C(edges=[(1,0), (0,1), (-1,0), (0,-1)])
            polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        A set of polygons can also be created over non-fields::

            sage: ConvexPolygons(ZZ)
            Category of convex real projective polygons over Integer Ring

        TESTS::

            sage: ConvexPolygons(QQ) is ConvexPolygons(QQ)
            True
            sage: TestSuite(ConvexPolygons(QQ)).run()
            sage: TestSuite(ConvexPolygons(QQbar)).run()
            sage: TestSuite(ConvexPolygons(ZZ)).run()
        """

        class ParentMethods:
            def is_convex(self):
                return True

    class SubcategoryMethods:
        def Convex(self):
            return self._with_axiom("Convex")

        def Rational(self):
            return self._with_axiom("Rational")

        def field(self):
            r"""
            Return the field over which this polygon is defined.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: P = polygons(vertices=[(0,0),(1,0),(2,1),(-1,1)])
                sage: P.field()
                doctest:warning
                ...
                UserWarning: field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() instead
                Rational Field

            """
            import warnings
            warnings.warn("field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() or base_ring().fraction_field() instead")

            return self.base_ring().fraction_field()

    class Rational(CategoryWithAxiom_over_base_ring):
        pass


all_axioms += ("Convex", )

ConvexPolygons = Polygons.Convex
