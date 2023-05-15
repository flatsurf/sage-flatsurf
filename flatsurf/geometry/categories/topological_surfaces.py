r"""
The category of topological surfaces.

This module provides a base category for all the surfaces in sage-flatsurf.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: C = MutableOrientedSimilaritySurface(QQ).category()

    sage: from flatsurf.geometry.categories import TopologicalSurfaces
    sage: C.is_subcategory(TopologicalSurfaces())
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

    In particular, this does not really require a topology since there is no
    general concept of open subsets of a surface in sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import TopologicalSurfaces
        sage: TopologicalSurfaces()
        Category of topological surfaces

    """

    def super_categories(self):
        return [TopologicalSpaces()]

    class ParentMethods:
        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in.

            The result of this method can be fed to ``_refine_category_`` to
            change the category of the surface (and enable functionality
            specific to the smaller classes of surfaces.)

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
            tester = self._tester(**options)

            tester.assertTrue(self.category().is_subcategory(self.refined_category()))

        @abstract_method
        def is_mutable(self):
            pass

        @abstract_method
        def is_orientable(self):
            r"""
            Return whether this surface is orientable.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom :class:`Orientable`.
                Surfaces can override this method to perform specialized logic,
                see the note in :mod:`flatsurf.geometry.categories` for
                performance considerations.

            EXAMPLES::

                sage: from flatsurf import polygon, similarity_surfaces
                sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_orientable()
                True

            """

        @abstract_method
        def is_with_boundary(self):
            r"""
            Return whether this a topological surface with boundary.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom :class:`WithBoundary`
                or :class:`WithoutBoundary`. Surfaces can override this method
                to perform specialized logic, see the note in
                :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import polygon, similarity_surfaces
                sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_with_boundary()
                False

            """

        @abstract_method
        def is_compact(self):
            r"""
            Return whether this surface is compact.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom of compactness.
                Surfaces can override this method to perform specialized logic,
                see the note in :mod:`flatsurf.geometry.categories` for
                performance considerations.

            EXAMPLES::

                sage: from flatsurf import polygon, similarity_surfaces
                sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_compact()
                True

            """

        @abstract_method
        def is_connected(self):
            r"""
            Return whether this surface is connected.

            .. NOTE::

                This method is used by :meth:`refined_category` to determine
                whether this surface satisfies the axiom of connectedness.
                Surfaces can override this method to perform specialized logic,
                see the note in :mod:`flatsurf.geometry.categories` for
                performance considerations.

            EXAMPLES::

                sage: from flatsurf import polygon, similarity_surfaces
                sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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

    class Orientable(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that can be oriented.

        As of 2023, all surfaces in sage-flatsurf satisfy this axiom.

        EXAMPLES::

            sage: from flatsurf import polygon, similarity_surfaces
            sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'Orientable' in S.category().axioms()
            True

        """

        class ParentMethods:
            def is_orientable(self):
                r"""
                Return whether this surface is orientable, i.e., return ``True``.

                EXAMPLES::

                    sage: from flatsurf import polygon, similarity_surfaces
                    sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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
            def __init__(self, *args, **kwargs):
                raise TypeError

    class WithoutBoundary(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that have no boundary, i.e., the
        surface is everywhere homeomorphic to the real plane.

        EXAMPLES::

            sage: from flatsurf import polygon, similarity_surfaces
            sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'WithoutBoundary' in S.category().axioms()
            True

        """

        class ParentMethods:
            def is_with_boundary(self):
                r"""
                Return whether this is a surface with boundary, i.e., return ``False``.

                EXAMPLES::

                    sage: from flatsurf import polygon, similarity_surfaces
                    sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_with_boundary()
                    False

                """
                return False

    class Connected(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def is_connected(self):
                return True

    class Compact(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def is_compact(self):
                r"""
                Return whether this surface is compact.

                EXAMPLES::

                    sage: from flatsurf import polygon, similarity_surfaces
                    sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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
