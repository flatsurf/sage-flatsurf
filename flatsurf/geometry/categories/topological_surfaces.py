r"""
The category of topological surfaces.

This module provides a base category for all surfaces in sage-flatsurf.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: S = MutableOrientedSimilaritySurface(QQ)

    sage: from flatsurf.geometry.categories import TopologicalSurfaces
    sage: S in TopologicalSurfaces()
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

from sage.categories.category_with_axiom import all_axioms
from sage.categories.topological_spaces import TopologicalSpaces
from sage.misc.abstract_method import abstract_method
from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)


class TopologicalSurfaces(SurfaceCategory):
    r"""
    The category of topological surfaces, i.e., surfaces that are locally
    homeomorphic to the real plane or the closed upper half plane.

    This category does not provide much functionality but just a common base
    for all the other categories defined in sage-flatsurf.

    In particular, this does not require a topology since there is no general
    concept of open subsets of a surface in sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import TopologicalSurfaces
        sage: TopologicalSurfaces()
        Category of topological surfaces

    """

    def super_categories(self):
        r"""
        Return the categories a topological surface is also a member of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import TopologicalSurfaces
            sage: TopologicalSurfaces().super_categories()
            [Category of topological spaces]

        """
        return [TopologicalSpaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces in sage-flatsurf.

        If you want to add functionality for all surfaces you most likely want
        to put it here.
        """

        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in.

            The result of this method can be fed to ``_refine_category_`` to
            change the category of the surface (and enable functionality
            specific to the smaller classes of surfaces).

            Note that this method does have much effect for a general
            topological surface. Subcategories and implementations of surfaces
            should override this method to derive more features.

            EXAMPLES::

                sage: from flatsurf import MutableOrientedSimilaritySurface
                sage: S = MutableOrientedSimilaritySurface(QQ)

                sage: from flatsurf import polygons
                sage: S.add_polygon(polygons.square(), label=0)
                0
                sage: S.refined_category()
                Category of connected with boundary finite type translation surfaces

                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 3))
                sage: S.refined_category()
                Category of connected without boundary finite type translation surfaces

            """
            category = self.category()

            if self.is_orientable():
                category &= category.Orientable()

            if self.is_with_boundary():
                category &= category.WithBoundary()
            else:
                category &= category.WithoutBoundary()

            if self.is_compact():
                category &= category.Compact()

            if self.is_connected():
                category &= category.Connected()

            return category

        def _test_refined_category(self, **options):
            r"""
            Verify that all (immutable) surfaces are contained in their refined
            category automatically.

            To pass this test, surfaces should either set their ``category``
            explicitly or ensure to run `_refine_category_(refined_category())`
            at some point.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_refined_category()

            """
            tester = self._tester(**options)

            tester.assertTrue(self.category().is_subcategory(self.refined_category()))

        def _Hom_(self, Y, category=None):
            r"""
            Return the space of morphisms from this surface to ``Y``.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: End(S)
                Surface Endomorphisms of Translation Surface in H_1(0) built from a square

            """
            if Y in TopologicalSurfaces():
                from flatsurf.geometry.morphism import SurfaceMorphismSpace

                return SurfaceMorphismSpace(self, Y, category=category)

            return super()._Hom_(Y, category=category)

        @abstract_method
        def is_mutable(self):
            r"""
            Return whether this surface allows modifications.

            All surfaces in sage-flatsurf must implement this method.

            .. NOTE::

                We do not specify the interface of such mutations. Any mutable
                surface should come up with a good interface for its use case. The
                point of this method is to signal that is likely unsafe to use this
                surface in caches (since it might change later) and that the
                category of the surface might still change.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_refined_category()
            """

        @abstract_method
        def is_orientable(self):
            r"""
            Return whether this surface is orientable.

            All surfaces in sage-flatsurf must implement this method.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom
                :class:`.TopologicalSurfaces.Orientable`. Surfaces must
                override this method to perform specialized logic, see the note
                in :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_orientable()
                True

            """

        @abstract_method
        def is_with_boundary(self):
            r"""
            Return whether this a topological surface with boundary.

            All surfaces in sage-flatsurf must implement this method.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom :class:`.WithBoundary`
                or :class:`.TopologicalSurfaces.WithoutBoundary`. Surfaces must
                override this method to perform specialized logic, see the note
                in :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_with_boundary()
                False

            """

        @abstract_method
        def is_compact(self):
            r"""
            Return whether this surface is compact.

            All surfaces in sage-flatsurf must implement this method.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom of compactness.
                Surfaces can override this method to perform specialized logic,
                see the note in :mod:`flatsurf.geometry.categories` for
                performance considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_compact()
                True

            """

        @abstract_method
        def is_connected(self):
            r"""
            Return whether this surface is connected.

            All surfaces in sage-flatsurf must implement this method.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom of connectedness.
                Surfaces can override this method to perform specialized logic,
                see the note in :mod:`flatsurf.geometry.categories` for
                performance considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_connected()
                True

            """

        @abstract_method(optional=True)
        def genus(self):
            r"""
            Return the genus of this surface.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: translation_surfaces.octagon_and_squares().genus()
                3

            """

        def _an_element_(self):
            r"""
            Return a point on this surface.

            EXAMPLES::

                sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                sage: s = SimilaritySurfaceGenerators.example()
                sage: s.an_element()
                Point (4/3, -2/3) of polygon 0

            ::

                sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface

                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
                0
                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 3))

                sage: S.an_element()
                Point (1/2, 1/2) of polygon 0

            TESTS:

            Verify that this method works over non-fields (if 2 is
            invertible)::

              sage: from flatsurf import similarity_surfaces
              sage: from flatsurf import EuclideanPolygonsWithAngles
              sage: E = EuclideanPolygonsWithAngles((3, 3, 5))
              sage: from pyexactreal import ExactReals # optional: pyexactreal  # random output due to pkg_resources deprecation warnings in some contexts
              sage: R = ExactReals(E.base_ring()) # optional: pyexactreal
              sage: angles = (3, 3, 5)
              sage: slopes = EuclideanPolygonsWithAngles(*angles).slopes()
              sage: P = Polygon(angles=angles, edges=[R.random_element() * slopes[0]])  # optional: pyexactreal
              sage: S = similarity_surfaces.billiard(P) # optional: pyexactreal
              sage: S.an_element()  # optional: pyexactreal
              Point ((1/2 ~ 0.50000000)*ℝ(0.303644…), 0) of polygon 0

            """
            return next(iter(self.some_elements()))

    class Orientable(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that can be oriented.

        As of 2023, all surfaces in sage-flatsurf satisfy this axiom.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'Orientable' in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all orientable surfaces in
            sage-flatsurf.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def is_orientable(self):
                r"""
                Return whether this surface is orientable, i.e., return ``True``.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_orientable()
                    True

                """
                return True

    class WithBoundary(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that have a boundary, i.e., at some
        points this surface is homeomorphic to the closed upper half plane.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square(), label=0)
            0
            sage: S.set_immutable()
            sage: 'WithBoundary' in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all surfaces with boundary in
            sage-flatsurf.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def is_with_boundary(self):
                r"""
                Return whether this is a surface with boundary, i.e., return ``True``.

                EXAMPLES::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: S = MutableOrientedSimilaritySurface(QQ)

                    sage: from flatsurf import polygons
                    sage: S.add_polygon(polygons.square(), label=0)
                    0
                    sage: S.set_immutable()
                    sage: S.is_with_boundary()
                    True

                """
                return True

        class WithoutBoundary(SurfaceCategoryWithAxiom):
            r"""
            An impossible category, the surfaces with and without boundary.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import TopologicalSurfaces
                sage: C = TopologicalSurfaces()
                sage: C.WithBoundary().WithoutBoundary()
                Traceback (most recent call last):
                ...
                TypeError: a surface cannot be both with and without boundary
                sage: C.WithoutBoundary().WithBoundary()
                Traceback (most recent call last):
                ...
                TypeError: a surface cannot be both with and without boundary

            """

            def __init__(self, *args, **kwargs):
                raise TypeError("a surface cannot be both with and without boundary")

    class WithoutBoundary(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that have no boundary, i.e., the
        surface is everywhere homeomorphic to the real plane.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'WithoutBoundary' in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all surfaces without boundary in
            sage-flatsurf.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def is_with_boundary(self):
                r"""
                Return whether this is a surface with boundary, i.e., return ``False``.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_with_boundary()
                    False

                """
                return False

    class Connected(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that are topologically connected.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'Connected' in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all connected surfaces in
            sage-flatsurf.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def is_connected(self):
                r"""
                Return whether this surface is connected, i.e., return
                ``True``.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_connected()
                    True

                """
                return True

    class Compact(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that are compact as topological spaces.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'Compact' in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all compact surfaces in
            sage-flatsurf.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def is_compact(self):
                r"""
                Return whether this surface is compact, i.e., return ``True``.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_compact()
                    True

                """
                return True

    class SubcategoryMethods:
        def Orientable(self):
            r"""
            Return the subcategory of surfaces that can be oriented.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import TopologicalSurfaces
                sage: TopologicalSurfaces().Orientable()
                Category of orientable topological surfaces

            """
            return self._with_axiom("Orientable")

        def WithBoundary(self):
            r"""
            Return the subcategory of surfaces that have a boundary, i.e.,
            points at which they are homeomorphic to the closed upper half
            plane.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import TopologicalSurfaces
                sage: TopologicalSurfaces().WithBoundary()
                Category of with boundary topological surfaces

            """
            return self._with_axiom("WithBoundary")

        def WithoutBoundary(self):
            r"""
            Return the subcategory of surfaces that have no boundary, i.e.,
            they are everywhere isomorphic to the real plane.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import TopologicalSurfaces
                sage: TopologicalSurfaces().WithoutBoundary()
                Category of without boundary topological surfaces

            """
            return self._with_axiom("WithoutBoundary")


# Currently, there is no "Orientable", "WithBoundary", and "WithoutBoundary"
# axiom in SageMath so we make it known to the category framework.
all_axioms += ("Orientable", "WithBoundary", "WithoutBoundary")
