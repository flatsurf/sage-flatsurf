r"""
The category of surfaces built from polygons.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from polygons (such as Euclidean polygons or hyperbolic polygons.)

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict(QQ).category()

    sage: from flatsurf.geometry.categories import PolygonalSurfaces
    sage: C.is_subcategory(PolygonalSurfaces())
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

from flatsurf.geometry.categories.surface_category import SurfaceCategory, SurfaceCategoryWithAxiom
from sage.categories.category_with_axiom import all_axioms
from sage.misc.abstract_method import abstract_method


class PolygonalSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing polygons defined in some space
    such as the projective plane (see
    :mod:`flatsurf.geometry.categories.real_projective_polygonal_surfaces`) or
    the hyperbolic plane (see
    :mod:`flatsurf.geometry.categories.hyperbolic_polygonal_surfaces`.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories import PolygonalSurfaces
        sage: PolygonalSurfaces()
        Category of polygonal surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.topological_surfaces import TopologicalSurfaces
        return [TopologicalSurfaces()]

    class ParentMethods:
        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in by
            consulting which edges are glued to each other.

            Note that this does not take into account how the edges are glued
            to each other exactly (i.e., by which similarity) since at this
            level (i.e., without knowing about the space in which the polygons
            live) the gluing is just described combinatorially.

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
            from flatsurf.geometry.categories.topological_surfaces import TopologicalSurfaces
            category = TopologicalSurfaces.ParentMethods.refined_category(self)

            try:
                finite_type = self.is_finite()
            except NotImplementedError:
                pass
            else:
                if finite_type:
                    category &= category.FiniteType()
                else:
                    category &= category.InfiniteType()

            return category

        def is_triangulated(self):
            if self.num_polygons() == 0:
                return True

            if self.polygon(self.base_label()).num_edges() != 3:
                return False

            raise NotImplementedError

        def walker(self):
            # TODO: Deprecate, use labels() instead
            from flatsurf.geometry.surface_legacy import LabelWalker
            return LabelWalker(self)

        def labels(self):
            r"""
            Return the labels used to enumerate the polygons that make up this surface.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.labels()
                (0,)

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.labels()
                (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

            """
            from flatsurf.geometry.surface import Labels
            return Labels(self)

        def num_polygons(self):
            # TODO: Deprecate
            if not self.is_finite():
                from sage.all import infinity
                return infinity
            return len(self.labels())

        def label_iterator(self, polygons=False):
            # TODO: Deprecate
            if polygons:
                for entry in self.label_polygon_iterator():
                    yield entry
            else:
                for entry in self.labels():
                    yield entry

        def edge_iterator(self, gluings=False):
            r"""
            Iterate over the edges of polygons, which are pairs (l,e) where l is a polygon label, 0 <= e < N and N is the number of edges of the polygon with label l.
            """
            # TODO: Deprecate
            if gluings:
                for entry in self.edge_gluing_iterator():
                    yield entry
                return
            for label, polygon in self.label_polygon_iterator():
                for edge in range(polygon.num_edges()):
                    yield label, edge

        def edge_gluing_iterator(self):
            r"""
            Iterate over the ordered pairs of edges being glued.
            """
            # TODO: Deprecate
            for label_edge_pair in self.edge_iterator():
                yield (
                    label_edge_pair,
                    self.opposite_edge(label_edge_pair[0], label_edge_pair[1]),
                )

        def label_polygon_iterator(self):
            r"""
            Iterate over pairs (label, polygon).

            Subclasses should consider overriding this method for increased
            performance.
            """
            # TODO: Deprecate
            for label in self.label_iterator():
                yield label, self.polygon(label)

        @abstract_method
        def base_label(self):
            pass

        @abstract_method
        def polygon(self, label):
            r"""
            Return the polygon with ``label``.

            INPUT:

            - ``label`` -- one of the labels included in :meth:`labels`

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.polygon(0)
                Polygon: (0, 0), (2, 0), (1, 4), (0, 5)

            """

        @abstract_method
        def opposite_edge(self, label, edge):
            r"""
            Return the polygon label and edge that is glued to the ``edge`` of
            the polygon with ``label``.

            INPUT:

            - ``label`` -- one of the labels included in :meth:`labels`

            - ``edge`` -- a non-negative integer to specify an edge (the edges
              of a polygon are numbered starting from zero.)

            OUTPUT:

            A tuple ``(label, edge)`` with the semantics as in the input.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.opposite_edge(0, 0)
                (0, 0)

            """

        # TODO: Deprecate is_finite everywhere since it clashes with the
        # notion of Sets(). Instead use is_finite_type().
        @abstract_method
        def is_finite(self):
            r"""
            Return whether this surface is constructed from finitely many polygons.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_finite()
                True

            """

        def is_with_boundary(self):
            r"""
            Return whether this surface has a boundary, i.e., unglued polygon edges.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_with_boundary()
                False

            """
            # Since this methods overrides the implementation from the axioms,
            # we reenable it.
            if 'WithBoundary' in self.category().axioms():
                return True
            if 'WithoutBoundary' in self.category().axioms():
                return False

            if not self.is_finite():
                raise NotImplementedError("cannot decide wether a surface has boundary for surfaces of infinite type")

            for label in self.label_iterator():
                for edge in range(self.polygon(label).num_edges()):
                    cross = self.opposite_edge(label, edge)
                    if cross is None:
                        return True

            return False

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
            if 'Compact' in self.category().axioms():
                return True

            if not self.is_finite():
                raise NotImplementedError("cannot decide whether this infinite type surface is compact")

            return True

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
            if 'Connected' in self.category().axioms():
                return True

            if not self.is_finite():
                raise NotImplementedError("cannot decide whether this infinite type surface is connected")

            # We use the customary union-find algorithm to identify the connected components.
            union_find = {label: label for label in self.label_iterator()}

            def find(label):
                if union_find[label] == label:
                    return label
                parent = find(union_find[label])
                union_find[label] = parent
                return parent

            for label in self.label_iterator():
                for edge in range(self.polygon(label).num_edges()):
                    cross = self.opposite_edge(label, edge)
                    if cross is None:
                        continue

                    cross_label, cross_edge = cross

                    x = find(label)
                    y = find(cross_label)
                    union_find[x] = y

            return len(union_find.values()) <= 1

        def _test_gluings(self, **options):
            # iterate over pairs with pair1 glued to pair2
            tester = self._tester(**options)

            if self.is_finite():
                it = self.label_iterator()
            else:
                from itertools import islice

                it = islice(self.label_iterator(), 30)

            for lab in it:
                p = self.polygon(lab)
                for k in range(p.num_edges()):
                    e = (lab, k)
                    f = self.opposite_edge(lab, k)
                    if f is None:
                        continue
                    g = self.opposite_edge(f[0], f[1])
                    tester.assertEqual(
                        e,
                        g,
                        "edge gluing is not a pairing:\n{} -> {} -> {}".format(e, f, g),
                    )

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces built from finitely many polygons.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'FiniteType' in S.category().axioms()
            True

        """

        class ParentMethods:
            def is_finite(self):
                return True

            def is_triangulated(self):
                for label in self.labels():
                    p = self.polygon(label)
                    if p.num_edges() != 3:
                        return False
                return True

    # TODO: Can we somehow force that a surface can only be finite XOR infinite type?
    class InfiniteType(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces built from infinitely many polygons.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: 'InfiniteType' in S.category().axioms()
            True

        """

        class ParentMethods:
            def is_finite(self):
                return False

    class Oriented(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by orientable surfaces with an orientation which
        is compatible with the orientation of the ambient space of the
        polygons (assuming that ambient space is orientable.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: 'Oriented' in S.category().axioms()
            True

        """

        def extra_super_categories(self):
            r"""
            Return the axioms that are automatically satisfied by a surfaces
            which is oriented, namely, that such a surface is orientable.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: 'Orientable' in S.category().axioms()
                True

            """
            from flatsurf.geometry.categories.topological_surfaces import TopologicalSurfaces
            return (TopologicalSurfaces().Orientable(),)

    class WithoutBoundary(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def _test_gluings_without_boundary(self, **options):
                # iterate over pairs with pair1 glued to pair2
                tester = self._tester(**options)

                if self.is_finite():
                    it = self.label_iterator()
                else:
                    from itertools import islice

                    it = islice(self.label_iterator(), 30)

                for lab in it:
                    p = self.polygon(lab)
                    for k in range(p.num_edges()):
                        f = self.opposite_edge(lab, k)
                        tester.assertFalse(
                            f is None, "edge ({}, {}) is not glued".format(lab, k)
                        )

    class SubcategoryMethods:
        def FiniteType(self):
            r"""
            Return the subcategory of surfaces built from finitely many polygons.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: PolygonalSurfaces().FiniteType()
                Category of finite type polygonal surfaces

            """
            return self._with_axiom("FiniteType")

        def InfiniteType(self):
            r"""
            Return the subcategory of surfaces built from infinitely many polygons.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: PolygonalSurfaces().InfiniteType()
                Category of infinite type polygonal surfaces

            """
            return self._with_axiom("InfiniteType")

        def Oriented(self):
            r"""
            Return the subcategory of surfaces with an orientation that is
            inherited from the polygons that it is built from.

            This assumes that the ambient space of the polygons is orientable.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: PolygonalSurfaces().Oriented()
                Category of oriented polygonal surfaces

            """
            return self._with_axiom("Oriented")


all_axioms += ("FiniteType", "InfiniteType", "Oriented")
