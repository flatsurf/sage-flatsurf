r"""
The category of topological surfaces.

This provides a base category for all the surfaces in sage-flatsurf.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict(QQ).category()

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

from sage.categories.category_with_axiom import CategoryWithAxiom, all_axioms
from sage.categories.topological_spaces import TopologicalSpaces
from sage.misc.abstract_method import abstract_method
from flatsurf.geometry.categories.surface_category import SurfaceCategory


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

                sage: from flatsurf import Surface_dict
                sage: S = Surface_dict(QQ)

                sage: from flatsurf import polygons
                sage: S.add_polygon(polygons.square(), label=0)
                0
                sage: S.refined_category()
                Category of compact connected with boundary finite type translation surfaces

                sage: S.set_edge_pairing(0, 0, 0, 2)
                sage: S.set_edge_pairing(0, 1, 0, 3)
                sage: S.refined_category()
                Category of compact connected without boundary finite type translation surfaces

            """
            category = self.category()

            # TODO: Force all surfaces to implement these methods and remove the try/except blocks.
            try:
                orientable = self.is_orientable()
            except NotImplementedError:
                pass
            else:
                if orientable:
                    category &= category.Orientable()
                else:
                    raise NotImplementedError("there is no axiom for non-orientable surfaces yet")

            try:
                with_boundary = self.is_with_boundary()
            except NotImplementedError:
                pass
            else:
                if with_boundary:
                    category &= category.WithBoundary()
                else:
                    category &= category.WithoutBoundary()

            try:
                compact = self.is_compact()
            except NotImplementedError:
                pass
            else:
                if compact:
                    category &= category.Compact()

            try:
                connected = self.is_connected()
            except NotImplementedError:
                pass
            else:
                if connected:
                    category &= category.Connected()

            return category

        @abstract_method
        def is_mutable(self):
            # TODO: Deprecate. We only really care whether a surface can be put into a cache here so we should test that instead. If it's mutable, it won't allow us to do that. (We should state that somewhere and test for it.)
            pass

        @abstract_method
        def is_orientable(self):
            r"""
            Return whether this surface is orientable.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_orientable()
                True

            """

        @abstract_method
        def is_with_boundary(self):
            r"""
            Return whether this a topological surface with boundary.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_with_boundary()
                False

            """

        @abstract_method
        def is_compact(self):
            r"""
            Return whether this surface is compact.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_compact()
                True

            """

        @abstract_method
        def is_connected(self):
            r"""
            Return whether this surface is connected.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_connected()
                True

            """

    class Orientable(CategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that can be oriented.

        As of 2023, all surfaces in sage-flatsurf satisfy this axiom.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'Orientable' in S.category().axioms()
            True

        """

        class ParentMethods:
            def is_orientable(self):
                r"""
                Return whether this surface is orientable, i.e., return ``True``.

                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_orientable()
                    True

                """
                return True

    class WithBoundary(CategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that have a boundary, i.e., at some
        points this surface is homeomorphic to the closed upper half plane.

        EXAMPLES::

            sage: from flatsurf import Surface_dict
            sage: S = Surface_dict(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square(), label=0)
            0
            sage: S.set_immutable()
            sage: 'WithBoundary' in S.category().axioms()
            True

        """

        def is_with_boundary(self):
            r"""
            Return whether this is a surface with boundary, i.e., return ``True``.

            EXAMPLES::

                sage: from flatsurf import Surface_dict
                sage: S = Surface_dict(QQ)

                sage: from flatsurf import polygons
                sage: S.add_polygon(polygons.square(), label=0)
                0
                sage: S.set_immutable()
                sage: S.is_with_boundary()
                True

            """
            return True

    # TODO: Can we somehow force that a surface can only be with XOR without boundary?
    class WithoutBoundary(CategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that have no boundary, i.e., the
        surface is everywhere homeomorphic to the real plane.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'WithoutBoundary' in S.category().axioms()
            True

        """

        def is_with_boundary(self):
            r"""
            Return whether this is a surface with boundary, i.e., return ``False``.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_with_boundary()
                False

            """
            return False

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


all_axioms += ("Orientable", "WithBoundary", "WithoutBoundary")
