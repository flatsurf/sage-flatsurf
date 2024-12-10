r"""
The category of surfaces built from polygons.

This module provides shared functionality for all surfaces in sage-flatsurf
that are built from polygons (such as Euclidean polygons or hyperbolic
polygons).

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

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
#        Copyright (C) 2023-2024 Julian Rüth
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
    such as the real plane (see
    :mod:`~flatsurf.geometry.categories.euclidean_polygonal_surfaces`).

    EXAMPLES::

        sage: from flatsurf.geometry.categories import PolygonalSurfaces
        sage: PolygonalSurfaces()
        Category of polygonal surfaces

    """

    def super_categories(self):
        r"""
        Return the categories that a polygonal surface is always also a member
        of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import PolygonalSurfaces
            sage: PolygonalSurfaces().super_categories()
            [Category of topological surfaces]

        """
        from flatsurf.geometry.categories.topological_surfaces import (
            TopologicalSurfaces,
        )

        return [TopologicalSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces that are built from polygons.

        If you want to add functionality for such surfaces you most likely
        want to put it here.
        """

        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in by
            consulting which edges are glued to each other.

            Note that this does not take into account how the edges are glued
            to each other exactly (e.g., by which similarity) since at this
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

            if self.is_finite_type():
                category &= category.FiniteType()
            else:
                category &= category.InfiniteType()

            return category

        def is_triangulated(self, limit=None):
            r"""
            Return whether this surface is built from triangles.

            Surfaces of infinite type should override this method.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.is_triangulated()
                False

            """
            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

            roots = self.roots()

            if not roots:
                return True

            for root in roots:
                if len(self.polygon(root).vertices()) != 3:
                    return False

            raise NotImplementedError(
                "cannot decide whether this (potentially infinite type) surface is triangulated"
            )

        def walker(self):
            r"""
            Return an iterable that walks the labels of the surface.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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

            return Labels(self, finite=self.is_finite_type())

        def base_label(self):
            r"""
            Return the polygon label from which iteration by :meth:`labels`
            should start.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.base_label()
                doctest:warning
                ...
                UserWarning: base_label() has been deprecated and will be removed in a future version of sage-flatsurf; use root() instead for connected surfaces and roots() in general
                0

            """
            import warnings

            warnings.warn(
                "base_label() has been deprecated and will be removed in a future version of sage-flatsurf; use root() instead for connected surfaces and roots() in general"
            )

            return self.root()

        def root(self):
            r"""
            Return the polygon label from which iteration by :meth:`labels`
            should start on this connected surface.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.root()
                0

            When there are multiple connected components, :meth:`roots` must be
            used instead::

                sage: from flatsurf import MutableOrientedSimilaritySurface, polygons
                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(polygons.square())
                0
                sage: S.add_polygon(polygons.square())
                1

                sage: S.root()
                Traceback (most recent call last):
                ...
                Exception: surface has more than one root label, use roots() instead
                sage: S.roots()
                (0, 1)

            """
            roots = self.roots()

            if not roots:
                raise Exception(
                    "cannot return a root label for the connected component on an empty surface, use roots() instead"
                )

            if len(roots) > 1:
                raise Exception(
                    "surface has more than one root label, use roots() instead"
                )

            return next(iter(roots))

        def polygons(self):
            r"""
            Return the polygons that make up this surface (in the same order as
            the labels are returned by :meth:`labels`)

            .. NOTE::

                Unlike with :meth:`labels`, this method should usually not be
                overridden. Things will be fast if :meth:`labels` is fast.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.polygons()
                (Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]),)

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.polygons()
                (Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), ...)

            """
            from flatsurf.geometry.surface import Polygons

            return Polygons(self)

        def _test_polygons(self, **options):
            r"""
            Verify that polygons() has been implemented correctly.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_polygons()

            """
            tester = self._tester(**options)

            polygons = self.polygons()

            if not self.is_finite_type():
                import itertools

                polygons = itertools.islice(polygons, 32)

            for polygon in polygons:
                tester.assertEqual(polygon.base_ring(), self.base_ring())

        def _test_labels_polygons(self, **options):
            r"""
            Verify that :meth:`labels` and :meth:`polygons` are compatible.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_labels_polygons()

            """
            tester = self._tester(**options)

            labels = self.labels()
            polygons = self.polygons()

            if not self.is_finite_type():
                import itertools

                labels = itertools.islice(labels, 32)

            for label, polygon in zip(labels, polygons):
                tester.assertEqual(self.polygon(label), polygon)

        def num_polygons(self):
            r"""
            Return the number of polygons that make up this surface.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.num_polygons()
                doctest:warning
                ...
                UserWarning: num_polygons() is deprecated and will be removed in a future version of sage-flatsurf; use len(polygons()) instead (and is_finite_type() for potentially infinite surfaces).
                1
                sage: len(S.polygons())
                1

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.num_polygons()
                +Infinity
                sage: S.is_finite_type()
                False

            """
            import warnings

            warnings.warn(
                "num_polygons() is deprecated and will be removed in a future version of sage-flatsurf; use len(polygons()) instead (and is_finite_type() for potentially infinite surfaces)."
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

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.label_iterator())
                doctest:warning
                ...
                UserWarning: label_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use labels() instead
                [0]
                sage: S.labels()
                (0,)

            """
            import warnings

            if polygons:
                warnings.warn(
                    "label_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use zip(labels(), polygons()) instead"
                )
                yield from zip(self.labels(), self.polygons())
            else:
                warnings.warn(
                    "label_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use labels() instead"
                )
                yield from self.labels()

        def edge_iterator(self, gluings=False):
            r"""
            Iterate over the edges of polygons, which are pairs (l,e) where l
            is a polygon label, 0 <= e < N and N is the number of edges of the
            polygon with label l.

            TESTS::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.edge_iterator())
                doctest:warning
                ...
                UserWarning: edge_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use edges() instead
                [(0, 0), (0, 1), (0, 2), (0, 3)]
                sage: S.edges()
                ((0, 0), (0, 1), (0, 2), (0, 3))

            ::

                sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(Polygon(edges=[(1,0),(0,1),(-1,-1)]))
                0
                sage: S.add_polygon(Polygon(edges=[(-1,0),(0,-1),(1,1)]))
                1
                sage: S.glue((0, 0), (1, 0))
                sage: S.glue((0, 1), (1, 1))
                sage: S.glue((0, 2), (1, 2))
                sage: list(S.edge_iterator())
                [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
                sage: S.edges()
                ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2))

            """
            import warnings

            if gluings:
                warnings.warn(
                    "edge_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use gluings() instead"
                )
                yield from self.gluings()
                return
            for label, polygon in zip(self.labels(), self.polygons()):
                warnings.warn(
                    "edge_iterator() has been deprecated and will be removed in a future version of sage-flatsurf; use edges() instead"
                )
                for edge in range(len(polygon.vertices())):
                    yield label, edge

        def edges(self):
            r"""
            Return the edges of the polygons that make up this surface as pairs
            (polygon label, edge index).

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.edges()
                ((0, 0), (0, 1), (0, 2), (0, 3))

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.edges()
                ((0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (-1, 0), (-1, 1), (-1, 2), (-1, 3), (2, 0), (2, 1), (2, 2), (2, 3), …)

            """
            from flatsurf.geometry.surface import Edges

            return Edges(self, finite=self.is_finite_type())

        def edge_gluing_iterator(self):
            r"""
            Iterate over the ordered pairs of edges being glued.

            TESTS::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: list(S.edge_gluing_iterator())
                doctest:warning
                ...
                UserWarning: edge_gluing_iterator() has been deprecated and will be removed in a future version of sage-flaturf; use gluings() instead
                [((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)), ((0, 3), (0, 3))]
                sage: S.gluings()
                (((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)), ((0, 3), (0, 3)))

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
            r"""
            Return the pairs of edges being glued to each other.

            Each gluing is reported as a pair of pairs (polygon label, edge
            index) and (glued polygon label, glued edge index).

            Note that each gluing is reported twice (unless it is a
            self-gluing).

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.gluings()
                (((0, 0), (0, 2)), ((0, 1), (0, 3)), ((0, 2), (0, 0)), ((0, 3), (0, 1)))

            A surface with only self-gluings::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.gluings()
                (((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)), ((0, 3), (0, 3)))

            """
            from flatsurf.geometry.surface import Gluings

            return Gluings(self)

        def label_polygon_iterator(self):
            r"""
            Iterate over pairs (label, polygon).

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: print(list(S.label_polygon_iterator()))
                doctest:warning
                ...
                UserWarning: label_polygon_iterator() has been deprecated and will be removed from a future version of sage-flatsurf; use zip(labels(), polygons()) instead
                [(0, Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]))]
                sage: print(list(zip(S.labels(), S.polygons())))
                [(0, Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]))]

            """
            import warnings

            warnings.warn(
                "label_polygon_iterator() has been deprecated and will be removed from a future version of sage-flatsurf; use zip(labels(), polygons()) instead"
            )

            return zip(self.labels(), self.polygons())

        @abstract_method
        def polygon(self, label):
            r"""
            Return the polygon with ``label``.

            INPUT:

            - ``label`` -- one of the labels included in :meth:`labels`

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.polygon(0)
                Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])

            """

        @abstract_method
        def opposite_edge(self, label, edge):
            r"""
            Return the polygon label and edge that is glued to the ``edge`` of
            the polygon with ``label``.

            INPUT:

            - ``label`` -- one of the labels included in :meth:`labels`

            - ``edge`` -- a non-negative integer to specify an edge (the edges
              of a polygon are numbered starting from zero).

            OUTPUT:

            A tuple ``(label, edge)`` with the semantics as in the input.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
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
                the :class:`~.PolygonalSurfaces.FiniteType` axiom or the
                :class:`~.PolygonalSurfaces.InfiniteType` axiom. Surfaces can
                override this method to perform specialized logic, see the note
                in :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_finite_type()
                True

            """

        def num_edges(self):
            r"""
            Return the total number of edges of all polygons used.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.num_edges()
                doctest:warning
                ...
                UserWarning: num_edges() has been deprecated and will be removed from a future version of sage-flatsurf; use sum(len(p.vertices()) for p in polygons()) instead
                4

            """
            import warnings

            warnings.warn(
                "num_edges() has been deprecated and will be removed from a future version of sage-flatsurf; use sum(len(p.vertices()) for p in polygons()) instead"
            )

            if self.is_finite_type():
                return sum(len(p.vertices()) for p in self.polygons())

            from sage.rings.infinity import Infinity

            return Infinity

        def _test_gluings(self, **options):
            r"""
            Verify that the gluings of this surface are consistent.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_gluings()

            """
            tester = self._tester(**options)

            if self.is_finite_type():
                it = self.labels()
            else:
                from itertools import islice

                it = islice(self.labels(), 30)

            for lab in it:
                p = self.polygon(lab)
                for k in range(len(p.vertices())):
                    e = (lab, k)
                    f = self.opposite_edge(lab, k)
                    if f is None:
                        continue
                    g = self.opposite_edge(*f)
                    tester.assertEqual(
                        e,
                        g,
                        "edge gluing is not a pairing:\n{} -> {} -> {}".format(e, f, g),
                    )

        @abstract_method
        def roots(self):
            r"""
            Return root labels for the polygons forming the connected
            components of this surface.

            A root label is just any of the :meth:`labels` of the surface.
            However, the iteration of :meth:`labels` starts from those root
            labels so for some surfaces they might have been specifically
            chosen.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.roots()
                (0,)

            """

        def is_connected(self):
            r"""
            Return whether this surface is connected.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_connected()
                True

            """
            return len(self.roots()) <= 1

        def component(self, root):
            r"""
            Return the labels contained in the connected component containing
            the label ``root``.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.component(0)
                (0,)

            """
            from flatsurf.geometry.surface import ComponentLabels

            return ComponentLabels(self, root)

        def components(self):
            r"""
            Return the connected components that make up this surface.

            OUTPUT:

            A sequence of connected components where each component is in turn
            a sequence of the polygon labels contained in that component.

            EXAMPLES::

                sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(polygons.square())
                0
                sage: S.add_polygon(polygons.square())
                1
                sage: S.add_polygon(polygons.square())
                2
                sage: S.glue((0, 0), (1, 0))
                sage: S.components()
                ((0, 1), (2,))

            """
            return tuple(self.component(root) for root in self.roots())

        def _test_components(self, **options):
            r"""
            Verify that :meth:`components` is compatible with :meth:`roots`.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_components()

            """
            tester = self._tester(**options)

            tester.assertEqual(len(self.components()), len(self.roots()))

    class ElementMethods:
        r"""
        Provides methods for all points on surfaces built from polygons.

        If you want to add functionality for such surfaces, you most likely
        want to put it here.
        """

        def is_in_edge_interior(self):
            r"""
            Return whether this point is on an edge (but not at a vertex) of
            one of the polygons that make up this surface.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

            A vertex is not contained in the interior of an edge::

                sage: S(0, 0).is_in_edge_interior()
                False

            A point on the edge that is not a vertex::

                sage: S(0, (1, 0)).is_in_edge_interior()
                True

            An inner point of a polygon::

                sage: S(0, (1, 1)).is_in_edge_interior()
                False

            """
            label, coordinates = self.representative()
            return (
                self.parent()
                .polygon(label)(coordinates)
                .position()
                .is_in_edge_interior()
            )

        def is_in_polygon_interior(self):
            r"""
            Return whether this point is in the interior of one of the
            polygons that make up this surface and not on an edge or at a
            vertex.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

            A vertex is not contained in the interior of a polygon::

                sage: S(0, 0).is_in_polygon_interior()
                False

            A point on an edge::

                sage: S(0, (1, 0)).is_in_polygon_interior()
                False

            An inner point of a polygon::

                sage: S(0, (1, 1)).is_in_polygon_interior()
                True

            """
            label, coordinates = self.representative()
            return self.parent().polygon(label)(coordinates).position().is_in_interior()

        def edges(self):
            r"""
            Return the edges of the polygons that contain this point.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

            For an inner point, no edges are reported::

                sage: S(0, (1, 1)).edges()
                set()

            For a point on a self-glued edge, one edge is reported::

                sage: S(0, (1, 0)).edges()
                {(0, 0)}

            For a point on a non self-glued edge, two edges are reported, for
            the two sides of the edge::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S(0, (1/2, 0)).edges()
                {(0, 0), (0, 2)}

            For a point on an unglued edge, a single edge is reported::

                sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (0, 1)]))
                0
                sage: S(0, (1/2, 0)).edges()
                {(0, 0)}

            All edges are reported for the vertex of this square torus::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S(0, 0).edges()
                {(0, 0), (0, 1), (0, 2), (0, 3)}

            """
            raise NotImplementedError(
                "points on this surface cannot determine which edges they are contained in yet"
            )

        def _test_edges(self, **options):
            r"""
            Verify that :meth:`edges` has been implemented correctly.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S(0, 0)._test_edges()

            """
            tester = self._tester(**options)

            edges = self.edges()

            if not edges:
                tester.assertTrue(self.is_in_polygon_interior())

            if len(edges) == 1:
                opposite = self.parent().opposite_edge(*edges[0])
                tester.assertTrue(opposite is None or opposite == edges[0])

            if self.is_vertex():
                tester.assertGreaterEqual(len(edges), 2)

            for edge in edges:
                opposite = self.parent().opposite_edge(*edge)
                if opposite is None:
                    continue

                tester.assertTrue(opposite in edges)

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces built from finitely many polygons.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'FiniteType' in S.category().axioms()
            True

        """

        def extra_super_categories(self):
            r"""
            Return the categories that surfaces built from finitely many
            polygons are additionally contained in; namely such a surface is a
            compact space.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: PolygonalSurfaces().FiniteType().extra_super_categories()
                (Category of compact topological spaces,)

            """
            from sage.categories.topological_spaces import TopologicalSpaces

            return (TopologicalSpaces().Compact(),)

        class InfiniteType(SurfaceCategoryWithAxiom):
            r"""
            The axiom satisfied by surfaces that are built from finitely and
            infinitely many polygons at the same time.

            This axiom does not exist and it is an error to create it.

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

        class Connected(SurfaceCategoryWithAxiom):
            r"""
            The axiom satisfied by connected surfaces built from finitely many polygons.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: PolygonalSurfaces().FiniteType().Connected()
                Category of connected finite type polygonal surfaces

            """

            class ParentMethods:
                r"""
                Provides methods available to all connected surfaces built from
                finitely many polygons.

                If you want to add functionality for such surfaces you most likely want
                to put it here.
                """

                def _test_roots(self, **options):
                    r"""
                    Verify that :meth:`roots` only reports a single connected
                    component.

                    EXAMPLES::

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                        sage: S = similarity_surfaces.self_glued_polygon(P)
                        sage: S._test_roots()

                    """
                    tester = self._tester(**options)

                    roots = self.roots()

                    for root in roots:
                        label = [label for label in self.labels() if label == root]
                        tester.assertEqual(len(label), 1)
                        tester.assertEqual(type(label[0]), type(root))

                    if not roots:
                        tester.assertTrue(not any(True for label in self.labels()))
                    else:
                        tester.assertTrue(next(iter(self.labels())) in roots)

        class Oriented(SurfaceCategoryWithAxiom):
            r"""
            The axiom satisfied by orientable surfaces with an orientation
            which is compatible with the orientation of the ambient space of
            the finitely many polygons that define the surface (assuming that
            that ambient space is orientable).

            EXAMPLES::

                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.octagon_and_squares()
                sage: S.category().is_subcategory(PolygonalSurfaces().FiniteType().Oriented())
                True

            """

            class ParentMethods:
                r"""
                Provides methods available to all surfaces built from finitely
                many oriented polygons.

                If you want to add functionality for such surfaces you most likely want
                to put it here.
                """

                def euler_characteristic(self):
                    r"""
                    Return the Euler characteristic of this surface.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.octagon_and_squares()
                        sage: S.euler_characteristic()
                        -4

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
                        previous = (edge - 1) % len(self.polygon(label).vertices())
                        cross = self.opposite_edge(label, previous)
                        if cross is None:
                            continue

                        union_find[find((label, edge))] = find(cross)

                    V = len({find((label, edge)) for (label, edge) in self.edges()})

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

        class ParentMethods:
            r"""
            Provides methods available to all surfaces built from finitely many polygons.

            If you want to add functionality for such surfaces you most likely want
            to put it here.
            """

            def is_finite_type(self):
                r"""
                Return whether this surface is built from finitely many polygons.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_finite_type()
                    True

                """
                return True

            def is_triangulated(self, limit=None):
                r"""
                Return whether this surfaces is built from triangles.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_triangulated()
                    False

                """
                if limit is not None:
                    import warnings

                    warnings.warn(
                        "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                        "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                    )

                for p in self.polygons():
                    if len(p.vertices()) != 3:
                        return False

                return True

            def is_with_boundary(self):
                r"""
                Return whether this surface has a boundary, i.e., unglued polygon edges.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_with_boundary()
                    False

                """
                for label in self.labels():
                    for edge in range(len(self.polygon(label).vertices())):
                        cross = self.opposite_edge(label, edge)
                        if cross is None:
                            return True

                return False

            def vertices(self):
                r"""
                Return the equivalence classes of the vertices of the polygons
                that make up this surface.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.vertices()
                    {Vertex 0 of polygon 0}

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.regular_octagon()
                    sage: S.vertices()
                    {Vertex 0 of polygon 0}

                """
                return {
                    # pylint: disable-next=not-callable
                    self(label, vertex)
                    for (label, vertex) in self.edges()
                }

            def _test_labels(self, **options):
                r"""
                Verify that :meth:`labels` has been implemented correctly by
                this surface.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S._test_labels()

                """
                tester = self._tester(**options)

                tester.assertEqual(len(list(self.labels())), len(self.labels()))

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
            r"""
            Provides methods available to all surfaces built from infinitely
            many polygons.

            If you want to add functionality for such surfaces you most likely want
            to put it here.
            """

            def is_finite_type(self):
                r"""
                Return whether this surfaces has been built from finitely many
                polygons which it has not.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.infinite_staircase()
                    sage: S.is_finite_type()
                    False

                """
                return False

    class Oriented(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by orientable surfaces with an orientation which is
        compatible with the orientation of the ambient space of the polygons
        (assuming that that ambient space is orientable).

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
            r"""
            The axiom satisfied by oriented surfaces built from polygons that
            have no unglued polygon edges.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: from flatsurf.geometry.categories import PolygonalSurfaces
                sage: S.category().is_subcategory(PolygonalSurfaces().Oriented().WithoutBoundary())
                True

            """

            class Connected(SurfaceCategoryWithAxiom):
                r"""
                The axiom satisfied by oriented connected surfaces built from
                polygons that have no unglued polygon edges.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.infinite_staircase()
                    sage: from flatsurf.geometry.categories import PolygonalSurfaces
                    sage: S.category().is_subcategory(PolygonalSurfaces().Oriented().WithoutBoundary().Connected())
                    True

                """

                class ParentMethods:
                    r"""
                    Provides methods available to all oriented connected
                    surfaces built from polygons without unglued edges.

                    If you want to add functionality for such surfaces you most
                    likely want to put it here.
                    """

                    def genus(self):
                        r"""
                        Return the genus of this surface.

                        ALGORITHM:

                        We deduce the genus from the Euler characteristic.

                        EXAMPLES::

                            sage: from flatsurf import translation_surfaces
                            sage: translation_surfaces.octagon_and_squares().genus()
                            3

                        This method might not be functional if the Euler
                        characteristic has not been implemented for the surface::

                            sage: S = translation_surfaces.infinite_staircase()
                            sage: S.genus()
                            Traceback (most recent call last):
                            ...
                            AttributeError: ... has no attribute 'euler_characteristic'...

                        """
                        return 1 - self.euler_characteristic() / 2

    class WithoutBoundary(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces built from polygons without any unglued
        edges.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: from flatsurf.geometry.categories import PolygonalSurfaces
            sage: S.category().is_subcategory(PolygonalSurfaces().WithoutBoundary())
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all surfaces that are built from
            polygons without unglued edges.

            If you want to add functionality for such surfaces you most likely want
            to put it here.
            """

            def _test_gluings_without_boundary(self, **options):
                r"""
                Verify that this surface has no unglued edges.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.infinite_staircase()
                    sage: S._test_gluings_without_boundary()

                """
                tester = self._tester(**options)

                if self.is_finite_type():
                    it = self.labels()
                else:
                    from itertools import islice

                    it = islice(self.labels(), 30)

                for lab in it:
                    p = self.polygon(lab)
                    for k in range(len(p.vertices())):
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


# Currently, there is no "FiniteType", "InfiniteType", and "Oriented"
# axiom in SageMath so we make them known to the category framework.
all_axioms += ("FiniteType", "InfiniteType", "Oriented")
