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
from sage.categories.category_with_axiom import (
    CategoryWithAxiom_over_base_ring,
    all_axioms,
)
from sage.misc.abstract_method import abstract_method

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

    @staticmethod
    def _describe_polygon(num_edges, **kwargs):
        r"""
        Return a printable description of a polygon with ``num_edges`` edges
        and additional features encoded as keyword arguments.

        The returned strings form a triple (indeterminate article, singular, plural).

        Usually, you don't call this method directly, but
        :meth:`Polygons.ParentMethods.describe_polygon`.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import Polygons
            sage: Polygons._describe_polygon(3)
            ('a', 'triangle', 'triangles')
            sage: Polygons._describe_polygon(3, equiangular=True, equilateral=True)
            ('an', 'equilateral triangle', 'equilateral triangles')
            sage: Polygons._describe_polygon(9, equiangular=True, equilateral=True)
            ('a', 'regular nonagon', 'regular nonagons')
            sage: Polygons._describe_polygon(4, equiangular=False, equilateral=True)
            ('a', 'rhombus', 'rhombi')
            sage: Polygons._describe_polygon(64, equiangular=False, equilateral=False)
            ('a', '64-gon', '64-gons')

        """
        from sage.all import infinity

        # From https://en.wikipedia.org/wiki/Polygon#Naming
        ngon_names = {
            1: ("a", "monogon"),
            2: ("a", "digon"),
            3: ("a", "triangle"),
            4: ("a", "quadrilateral"),
            5: ("a", "pentagon"),
            6: ("a", "hexagon"),
            7: ("a", "heptagon"),
            8: ("an", "octagon"),
            9: ("a", "nonagon"),
            10: ("a", "decagon"),
            11: ("a", "hendecagon"),
            12: ("a", "dodecagon"),
            13: ("a", "tridecagon"),
            14: ("a", "tetradecagon"),
            15: ("a", "pentadecagon"),
            16: ("a", "hexadecagon"),
            17: ("a", "heptadecagon"),
            18: ("an", "octadecagon"),
            19: ("an", "enneadecagon"),
            20: ("an", "icosagon"),
            # Most people probably don't know the prefixes after that. We
            # keep a few easy/fun ones.
            100: ("a", "hectogon"),
            1000: ("a", "chiliagon"),
            10000: ("a", "myriagon"),
            1000000: ("a", "megagon"),
            infinity: ("an", "apeirogon"),
        }

        description = ngon_names.get(num_edges, ("a", f"{num_edges}-gon"))
        description = description + (description[1] + "s",)

        def augment(article, *attributes):
            nonlocal description
            description = (
                article,
                " ".join(attributes + (description[1],)),
                " ".join(attributes + (description[2],)),
            )

        def augment_if(article, attribute, *properties):
            if all(
                kwargs.get(property, False) for property in (properties or [attribute])
            ):
                augment(article, attribute)
                return True
            return False

        def augment_if_not(article, attribute, *properties):
            if all(
                kwargs.get(property, True) is False
                for property in (properties or [attribute])
            ):
                augment(article, attribute)
                return True
            return False

        if augment_if("a", "degenerate"):
            return description

        if num_edges == 3:
            if not augment_if("an", "equilateral", "equiangular"):
                if not augment_if("an", "isosceles"):
                    augment_if("a", "right")

            return description

        if num_edges == 4:
            if kwargs.get("equilateral", False) and kwargs.get("equiangular", False):
                return "a", "square", "squares"

            if kwargs.get("equiangular", False):
                return "a", "rectangle", "rectangles"

            if kwargs.get("equilateral", False):
                return "a", "rhombus", "rhombi"

        augment_if("a", "regular", "equilateral", "equiangular")

        augment_if_not("a", "non-convex", "convex")

        marked_vertices = kwargs.get("marked_vertices", 0)
        if marked_vertices:
            if marked_vertices == 1:
                suffix = "with a marked vertex"
            else:
                suffix = f"with {kwargs.get('marked_vertices')} marked vertices"
            description = (
                description[0],
                f"{description[1]} {suffix}",
                f"{description[2]} {suffix}",
            )

        return description

    class ParentMethods:
        r"""
        Provides methods available to all polygons.

        If you want to add functionality to all polygons, independent of
        implementation, you probably want to put it here.
        """

        @abstract_method
        def change_ring(self, ring):
            r"""
            Return a copy of this polygon which is defined over ``ring``.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2)**(1/2))
                sage: S.change_ring(K)
                Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

            """

        @abstract_method
        def is_convex(self, strict=False):
            r"""
            Return whether this is a convex polygon.

            INPUT:

            - ``strict`` -- whether to check for strict convexity, i.e., a
              polygon with a π angle is not considered convex.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.is_convex()
                True
                sage: S.is_convex(strict=True)
                True

            """

        @abstract_method
        def is_degenerate(self):
            r"""
            Return whether this polygon is considered degenerate.

            EXAMPLES:

            Polygons with zero area are considered degenerate::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (1, 0)], check=False)
                sage: p.is_degenerate()
                True

            Polygons with marked vertices are considered degenerate::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (4, 0), (2, 2)])
                sage: p.is_degenerate()
                True

            """

        @abstract_method
        def is_simple(self):
            r"""
            Return whether this polygon is not self-intersecting.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.is_simple()
                True

            """

        def base_ring(self):
            r"""
            Return the ring over which this polygon is defined.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.base_ring()
                Rational Field

            """
            return self.category().base_ring()

        def describe_polygon(self):
            r"""
            Return a textual description of this polygon for generating
            human-readable messages.

            The description is returned as a triple (indeterminate article,
            singular, plural).

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.describe_polygon()
                ('a', 'square', 'squares')

            """
            marked_vertices = set(self.vertices()).difference(
                self.vertices(marked_vertices=False)
            )

            if marked_vertices and self.area() != 0:
                self = self.erase_marked_vertices()

            properties = {
                "degenerate": self.is_degenerate(),
                "equilateral": self.is_equilateral(),
                "equiangular": self.is_equiangular(),
                "convex": self.is_convex(),
                "marked_vertices": len(marked_vertices),
            }

            if len(self.vertices()) == 3:
                slopes = self.slopes(relative=True)
                properties["right"] = any(slope[0] == 0 for slope in slopes)

                from flatsurf.geometry.euclidean import is_parallel

                properties["isosceles"] = (
                    is_parallel(slopes[0], slopes[1])
                    or is_parallel(slopes[0], slopes[2])
                    or is_parallel(slopes[1], slopes[2])
                )

            return Polygons._describe_polygon(len(self.vertices()), **properties)

        def _test_refined_category(self, **options):
            r"""
            Verify that the polygon is in the smallest category that can be
            easily determined.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: P = polygons.square()
                sage: P._test_refined_category()

            """
            tester = self._tester(**options)

            if self.is_convex():
                tester.assertTrue("Convex" in self.category().axioms())

            if self.is_simple():
                tester.assertTrue("Simple" in self.category().axioms())

            # We do not test for Rational and WithAngles since these are only
            # determined on demand (computing angles can be very costly).

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

            def is_convex(self, strict=False):
                r"""
                Return whether this is a convex polygon, which it is.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: P = polygons.square()
                    sage: P.is_convex()
                    True

                """
                if strict:
                    raise NotImplementedError(
                        "cannot decide strict convexity for this polygon yet"
                    )

                return True

    class Simple(CategoryWithAxiom_over_base_ring):
        r"""
        The axiom satisfied by polygons that are not self-intersecting.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import Polygons
            sage: C = Polygons(QQ)
            sage: C.Simple()
            Category of simple polygons over Rational Field

        """

        class ParentMethods:
            r"""
            Provides methods available to all simple polygons.

            If you want to add functionality to all such polygons, you probably
            want to put it here.
            """

            def is_simple(self):
                r"""
                Return whether this polygon is not self-intersecting, i.e.,
                return ``True``.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: s = polygons.square()
                    sage: s.is_simple()
                    True

                """
                return True

    class Rational(CategoryWithAxiom_over_base_ring):
        r"""
        The axiom satisfied by polygons whose inner angles are rational
        multiples of π.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import Polygons
            sage: C = Polygons(QQ)
            sage: C.Rational()
            Category of rational polygons over Rational Field

        """

        class ParentMethods:
            r"""
            Provides methods available to all rational polygons.

            If you want to add functionality to all such polygons, you probably
            want to put it here.
            """

            def is_rational(self):
                r"""
                Return whether all inner angles of this polygon are rational
                multiples of π, i.e., return ``True``.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: s = polygons.square()
                    sage: s.is_rational()
                    True

                """
                return True

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

        def Simple(self):
            r"""
            Return the subcategyr of simple polygons.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import Polygons
                sage: Polygons(QQ).Simple()
                Category of simple polygons over Rational Field

            """
            return self._with_axiom("Simple")

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

            warnings.warn(
                "field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() or base_ring().fraction_field() instead"
            )

            return self.base_ring().fraction_field()

        @cached_method
        def base_ring(self):
            r"""
            Return the ring over which the polygons in this category are
            defined.

                sage: from flatsurf.geometry.categories import Polygons
                sage: C = Polygons(QQ).Rational().Simple().Convex()

                sage: C.base_ring()
                Rational Field

            """
            # Copied this trick from SageMath's Modules category.
            for C in self.super_categories():
                if hasattr(C, "base_ring"):
                    return C.base_ring()
            assert False

        def _test_base_ring(self, **options):
            tester = self._tester(**options)

            from sage.categories.all import Rings

            tester.assertTrue(self.base_ring() in Rings())

        def change_ring(self, ring):
            r"""
            Return this category but defined over the ring ``ring``.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: C = s.category()
                sage: C
                Category of convex simple euclidean rectangles over Rational Field
                sage: C.change_ring(AA)
                Category of convex simple euclidean rectangles over Algebraic Real Field

            """
            from sage.categories.category import JoinCategory

            if isinstance(self, JoinCategory):
                from sage.categories.category import Category

                return Category.join(
                    [S.change_ring(ring) for S in self.super_categories()]
                )

            # This is a hack to make the change ring of EuclideanPolygonsWithAngles subcategories work
            if hasattr(self, "angles"):
                return type(self)(ring, self.angles())

            if isinstance(self, Category_over_base_ring):
                return type(self)(ring)

            raise NotImplementedError("cannot change_ring() of this category yet")


# Currently, there is no "Convex" and "Simple" axiom in SageMath so we make it
# known to the category framework. Note that "Rational" is already defined by
# topological surfaces so we don't need to declare it here again.
all_axioms += (
    "Convex",
    "Simple",
)
