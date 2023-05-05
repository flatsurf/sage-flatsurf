r"""
The category of surfaces built from polygons.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from polygons (such as Euclidean polygons or hyperbolic polygons.)

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: C = MutableOrientedSimilaritySurface(QQ).category()

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

from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)
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
        from flatsurf.geometry.categories.topological_surfaces import (
            TopologicalSurfaces,
        )

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

                sage: from flatsurf import MutableOrientedSimilaritySurface, polygons
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
            from flatsurf.geometry.categories.topological_surfaces import (
                TopologicalSurfaces,
            )

            category = TopologicalSurfaces.ParentMethods.refined_category(self)

            try:
                finite_type = self.is_finite_type()
            except NotImplementedError:
                pass
            else:
                if finite_type:
                    category &= category.FiniteType()
                else:
                    category &= category.InfiniteType()

            return category

        def is_triangulated(self):
            if len(self.polygons()) == 0:
                return True

            if self.polygon(self.base_label()).num_edges() != 3:
                return False

            raise NotImplementedError

        def walker(self):
            r"""
            Return an iterable that walks the labels of the surface.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: walker = S.walker()
                doctest:warning
                ...
                UserWarning: walker() is deprecated and will be removed from a future version of sage-flatsurf; use labels() instead.
                sage: list(walker)
                [0]

            """
            import warnings

            warnings.warn(
                "walker() is deprecated and will be removed from a future version of sage-flatsurf; use labels() instead."
            )

            from flatsurf.geometry.surface_legacy import LabelWalker

            return LabelWalker(self, deprecation_warning=False)

        def labels(self):
            r"""
            Return the labels used to enumerate the polygons that make up this
            surface.

            The labels are returned in a breadth first search order starting at
            the :meth:`base_label`. This order is compatible with the order in
            which polygons are returned by :meth:`polygons`.

            .. NOTE::

                The generic implementation of this method returns a collection
                that is very slow at computing its length and deciding
                containment. To speed things up it is recommended to override
                this method.

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

        def polygons(self):
            r"""
            Return the polygons that make up this surface (in the same order as
            the labels are returned by :meth:`labels`)

            .. NOTE::

                Unlike with :meth:`labels`, this method should usually not be
                overriden. Things will be fast if :meth:`labels` is fast.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.polygons()
                (polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]),)

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.polygons()
                (polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), ...)

            """
            from flatsurf.geometry.surface import Polygons

            return Polygons(self)

        def _test_labels_polygons(self, **options):
            tester = self._tester(**options)

            labels = self.labels()
            polygons = self.polygons()

            if not self.is_finite_type():
                import itertools

                labels = itertools.islice(labels, 32)

            for label, polygon in zip(labels, polygons):
                tester.assertEqual(self.polygon(label), polygon)

        def num_polygons(self):
            import warnings

            warnings.warn(
                "num_polygons() is deprecated and will be removed in a future version of sage-flatsurf; use len(polygons()) instead."
            )

            # Note that using len(self.polygons()) on
            # MutableOrientedSimilaritySurface is only very slightly slower:
            # %timeit num_polygons()
            # 137ns
            # %timeit len(polygons())
            # 159ns

            # On other surfaces, the effect can be much more pronounced. The
            # overhead of calling through the category framework and creating a
            # Labels object can lead to runtimes of about 1μs.

            if not self.is_finite_type():
                from sage.all import infinity

                return infinity
            return len(self.polygons())

        def label_iterator(self, polygons=False):
            r"""
            TESTS::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.label_iterator())
                doctest:warning
                ...
                UserWarning: label_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use labels() instead
                [0]

            """
            import warnings

            if polygons:
                warnings.warn(
                    "label_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use zip(labels(), polygons()) instead"
                )
                for entry in zip(self.labels(), self.polygons()):
                    yield entry
            else:
                warnings.warn(
                    "label_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use labels() instead"
                )
                for entry in self.labels():
                    yield entry

        def edge_iterator(self, gluings=False):
            r"""
            Iterate over the edges of polygons, which are pairs (l,e) where l
            is a polygon label, 0 <= e < N and N is the number of edges of the
            polygon with label l.

            TESTS::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.edge_iterator())
                doctest:warning
                ...
                UserWarning: edge_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use edges() instead
                [(0, 0), (0, 1), (0, 2), (0, 3)]

            ::

                sage: from flatsurf import polygon, MutableOrientedSimilaritySurface
                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(polygon(edges=[(1,0),(0,1),(-1,-1)]))
                0
                sage: S.add_polygon(polygon(edges=[(-1,0),(0,-1),(1,1)]))
                1
                sage: S.glue((0, 0), (1, 0))
                sage: S.glue((0, 1), (1, 1))
                sage: S.glue((0, 2), (1, 2))
                sage: for edge in S.edge_iterator():
                ....:     print(edge)
                (0, 0)
                (0, 1)
                (0, 2)
                (1, 0)
                (1, 1)
                (1, 2)

            """
            import warnings

            if gluings:
                warnings.warn(
                    "edge_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use gluings() instead"
                )
                for entry in self.gluings():
                    yield entry
                return
            for label, polygon in zip(self.labels(), self.polygons()):
                warnings.warn(
                    "edge_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use edges() instead"
                )
                for edge in range(polygon.num_edges()):
                    yield label, edge

        def edges(self):
            from flatsurf.geometry.surface import Edges

            return Edges(self)

        def edge_gluing_iterator(self):
            r"""
            Iterate over the ordered pairs of edges being glued.

            TESTS::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.edge_gluing_iterator())
                doctest:warning
                ...
                UserWarning: edge_gluing_iterator() has been deprecated and will be removed in a future version of sage-flaturf; use gluings() instead
                [((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)), ((0, 3), (0, 3))]

            """
            import warnings

            warnings.warn(
                "edge_gluing_iterator() has been deprecated and will be removed in a future version of sage-flaturf; use gluings() instead"
            )

            for label_edge_pair in self.edges():
                yield (
                    label_edge_pair,
                    self.opposite_edge(label_edge_pair[0], label_edge_pair[1]),
                )

        def gluings(self):
            from flatsurf.geometry.surface import Gluings

            return Gluings(self)

        def label_polygon_iterator(self):
            r"""
            Iterate over pairs (label, polygon).

            Subclasses should consider overriding this method for increased
            performance.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.label_polygon_iterator())
                doctest:warning
                ...
                UserWarning: label_polygon_iterator() has been deprecated and will be removed from a future version of sage-flatsurf; use zip(labels(), polygons()) instead
                [(0, polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]))]

            """
            import warnings

            warnings.warn(
                "label_polygon_iterator() has been deprecated and will be removed from a future version of sage-flatsurf; use zip(labels(), polygons()) instead"
            )

            return zip(self.labels(), self.polygons())

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
                polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])

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

        def is_finite(self):
            r"""
            Return whether this surface is constructed from finitely many polygons.

            .. NOTE::

                The semantics of this function clash with the notion of finite
                sets inherited from the category of sets. Therefore
                :meth:`is_finite_type` should be used instead.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_finite_type()
                True

            """
            import warnings

            warnings.warn(
                "is_finite() has been deprecated and will be removed in a future version of sage-flatsurf; use is_finite_type() instead"
            )
            return self.is_finite_type()

        @abstract_method
        def is_finite_type(self):
            r"""
            Return whether this surface is constructed from finitely many polygons.

            .. NOTE::

                This method is used to determine whether this surface satisfies
                the :class:`FiniteType` axiom or the :class:`InfiniteType`
                axiom. Surfaces can override this method to perform specialized
                logic, see the note in :mod:`flatsurf.geometry.categories` for
                performance considerations.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_finite_type()
                True

            """

        def num_edges(self):
            r"""
            Return the total number of edges of all polygons used.
            """
            if self.is_finite_type():
                try:
                    return self._cache["num_edges"]
                except KeyError:
                    num_edges = sum(p.num_edges() for p in self.polygons())
                    self._cache["num_edges"] = num_edges
                    return num_edges
            else:
                from sage.rings.infinity import Infinity

                return Infinity

        def _test_gluings(self, **options):
            # iterate over pairs with pair1 glued to pair2
            tester = self._tester(**options)

            if self.is_finite_type():
                it = self.labels()
            else:
                from itertools import islice

                it = islice(self.labels(), 30)

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

        def extra_super_categories(self):
            from sage.categories.topological_spaces import TopologicalSpaces

            return (TopologicalSpaces().Compact(),)

        class InfiniteType(SurfaceCategoryWithAxiom):
            r"""
            TESTS::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: PolygonalSurfaces().FiniteType() & PolygonalSurfaces().InfiniteType()
                Traceback (most recent call last):
                ...
                TypeError: surface cannot be finite type and infinite type at the same time
                sage: PolygonalSurfaces().InfiniteType() & PolygonalSurfaces().FiniteType()
                Traceback (most recent call last):
                ...
                TypeError: surface cannot be finite type and infinite type at the same time

            """

            def __init__(self, *args, **kwargs):
                raise TypeError(
                    "surface cannot be finite type and infinite type at the same time"
                )

        class ParentMethods:
            def is_finite_type(self):
                return True

            def is_triangulated(self):
                for label in self.labels():
                    p = self.polygon(label)
                    if p.num_edges() != 3:
                        return False
                return True

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
                for label in self.labels():
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
                # We use the customary union-find algorithm to identify the connected components.
                union_find = {label: label for label in self.labels()}

                def find(label):
                    if union_find[label] == label:
                        return label
                    parent = find(union_find[label])
                    union_find[label] = parent
                    return parent

                for label in self.labels():
                    for edge in range(self.polygon(label).num_edges()):
                        cross = self.opposite_edge(label, edge)
                        if cross is None:
                            continue

                        cross_label, cross_edge = cross

                        x = find(label)
                        y = find(cross_label)
                        union_find[x] = y

                roots = set(find(label) for label in self.labels())

                return len(roots) <= 1

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
            def is_finite_type(self):
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
            from flatsurf.geometry.categories.topological_surfaces import (
                TopologicalSurfaces,
            )

            return (TopologicalSurfaces().Orientable(),)

        class WithoutBoundary(SurfaceCategoryWithAxiom):
            class Connected(SurfaceCategoryWithAxiom):
                class ParentMethods:
                    def genus(self):
                        r"""
                        Return the genus of this surface.

                        ALGORITHM:

                        We deduce the genus from the Euler characteristic.

                        EXAMPLES::

                            sage: from flatsurf import translation_surfaces
                            sage: translation_surfaces.octagon_and_squares().genus()
                            3

                        """
                        return 1 - self.euler_characteristic() / 2

        class FiniteType(SurfaceCategoryWithAxiom):
            class ParentMethods:
                def euler_characteristic(self):
                    r"""
                    Return the Euler characteristic of this surface.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.octagon_and_squares()
                        sage: S.euler_characteristic()
                        -4

                    .. [Massart2021] \D. Massart. A short introduction to translation
                    surfaces, Veech surfaces, and Teichműller dynamics.
                    https://hal.science/hal-03300179/document

                    """
                    # Count the vertices
                    union_find = {edge: edge for edge in self.edges()}

                    def find(node):
                        if union_find[node] == node:
                            return node
                        parent = find(union_find[node])
                        union_find[node] = parent
                        return parent

                    for label, edge in self.edges():
                        previous = (edge - 1) % self.polygon(label).num_edges()
                        cross = self.opposite_edge(label, previous)
                        if cross is None:
                            continue

                        union_find[find((label, edge))] = find(cross)

                    V = len(set(find((label, edge)) for (label, edge) in self.edges()))

                    # Count the edges
                    from sage.all import QQ, ZZ

                    E = QQ(0)
                    for label, edge in self.edges():
                        if self.opposite_edge(label, edge) is None:
                            E += 1
                        elif self.opposite_edge(label, edge) == (label, edge):
                            E += 1
                            V += 1
                        else:
                            E += 1 / 2
                    assert E in ZZ

                    # Count the faces
                    F = len(self.polygons())

                    return ZZ(V - E + F)

    class WithoutBoundary(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def _test_gluings_without_boundary(self, **options):
                # iterate over pairs with pair1 glued to pair2
                tester = self._tester(**options)

                if self.is_finite_type():
                    it = self.labels()
                else:
                    from itertools import islice

                    it = islice(self.labels(), 30)

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
