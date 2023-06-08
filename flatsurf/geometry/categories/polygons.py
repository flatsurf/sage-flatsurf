r"""
The category of polyogons

This module provides shared functionality for all polygons in sage-flatsurf.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category of a polygon is automatically determined.

EXAMPLES::

    sage: from flatsurf.geometry.categories import Polygons
    sage: C = Polygons(QQ)

    sage: from flatsurf import polygons
    sage: polygons.square() in C
    True

"""
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
from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring, all_axioms

from sage.categories.all import Sets


class Polygons(Category_over_base_ring):
    r"""
    The category of polygons defined over a base ring.

    This comprises arbitrary base ring, e.g., this category contains Euclidean
    polygons and hyperbolic polygons.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import Polygons
        sage: Polygons(QQ)
        Category of polygons over Rational Field

    """
    def super_categories(self):
        r"""
        Return the categories polygons automatically belong to.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import Polygons
            sage: C = Polygons(QQ)
            sage: C.super_categories()
            [Category of sets]

        """
        return [Sets()]

    class ParentMethods:
        r"""
        Provides methods available to all polygons.

        If you want to add functionality to all polygons, independent of
        implementation, you probably want to put it here.
        """
        def _check(self):
            r"""
            Verify that this is a valid polygon.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: P = polygons.square()
                sage: P._check()

            """
            if self.area() < 0:
                raise ValueError("polygon has negative area; probably the vertices are not in counter-clockwise order?")

    class Convex(CategoryWithAxiom_over_base_ring):
        r"""
        The axiom satisfied by convex polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import Polygons
            sage: C = Polygons(QQ)
            sage: C.Convex()
            Category of convex polygons over Rational Field

        """
        class ParentMethods:
            r"""
            Provides methods available to all convex polygons.

            If you want to add functionality to all such polygons, you probably
            want to put it here.
            """

            def is_convex(self):
                r"""
                Return whether this is a convex polygon, which it is.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: P = polygons.square()
                    sage: P.is_convex()
                    True

                """
                return True

    class Rational(CategoryWithAxiom_over_base_ring):
        r"""
        The axiom satisfied by polygons whose inner angles are rational
        multiples of π.

        EXAMPLES::

        """
        pass

    class SubcategoryMethods:
        def Convex(self):
            r"""
            Return the subcategory of convex polygons.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import Polygons
                sage: Polygons(QQ).Convex()
                Category of convex polygons over Rational Field

            """
            return self._with_axiom("Convex")

        def Rational(self):
            r"""
            Return the subcategory of polygons with rational angles.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import Polygons
                sage: Polygons(QQ).Rational()
                Category of rational polygons over Rational Field

            """
            return self._with_axiom("Rational")

        def field(self):
            r"""
            Return the field over which these polygons are defined.

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: P = Polygon(vertices=[(0,0),(1,0),(2,1),(-1,1)])
                sage: P.category().field()
                doctest:warning
                ...
                UserWarning: field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() or base_ring().fraction_field() instead
                Rational Field

            """
            import warnings
            warnings.warn("field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() or base_ring().fraction_field() instead")

            return self.base_ring().fraction_field()

        @cached_method
        def base_ring(self):
            # Copied this trick from SageMath's Modules category.
            for C in self.super_categories():
                if hasattr(C, "base_ring"):
                    return C.base_ring()
            assert False

        def _test_base_ring(self, **options):
            tester = self._tester(**options)

            from sage.categories.all import Rings
            tester.assertTrue(self.base_ring() in Rings())


# Currently, there is no "Convex" axiom in SageMath so we make it known to the
# category framework. Note that "Rational" is already defined by topological
# surfaces so we don't need to declare it here again.
all_axioms += ("Convex", )
