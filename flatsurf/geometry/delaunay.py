r"""
Triangulations, Delaunay triangulations, and Delaunay decompositions of
infinite surfaces.

EXAMPLES:

Typically, you don't need to create these surfaces directly, they are created
by invoking methods on the underlying surfaces::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.infinite_staircase()
    sage: S.triangulate()
    Triangulation of The infinite staircase

    sage: S.delaunay_triangulation()
    Delaunay triangulation of The infinite staircase

    sage: S.delaunay_decomposition()
    Delaunay cell decomposition of The infinite staircase

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
#                          2023 Julian Rüth
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
# ********************************************************************

from sage.misc.cachefunc import cached_method

from flatsurf.geometry.surface import (
    MutableOrientedSimilaritySurface_base,
    OrientedSimilaritySurface,
    Labels,
)


class LazyTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    Surface class used to triangulate an infinite surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()
        sage: S = S.triangulate()

    TESTS::

        sage: from flatsurf.geometry.delaunay import LazyTriangulatedSurface
        sage: isinstance(S, LazyTriangulatedSurface)
        True
        sage: TestSuite(S).run()  # long time (1s)

    """

    def __init__(self, similarity_surface, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyTriangulatedSurface; use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyTriangulatedSurface()"
                )

        if similarity_surface.is_mutable():
            raise ValueError("Surface must be immutable.")

        self._reference = similarity_surface

        OrientedSimilaritySurface.__init__(
            self,
            similarity_surface.base_ring(),
            category=category or self._reference.category(),
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.is_mutable()
            False

        """
        return False

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.is_compact()
            False

        """
        return self._reference.is_compact()

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.roots()
            ((0, (0, 1, 2)),)

        """
        return tuple(
            (reference_label, self._triangulation(reference_label)[(0, 1)])
            for reference_label in self._reference.roots()
        )

    def _triangulation(self, reference_label):
        r"""
        Return a triangulated of the ``reference_label`` in the underlying
        (typically non-triangulated) reference surface.

        INPUT:

        - ``reference_label`` -- a polygon label in the reference surface that
          we are triangulating.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S._triangulation(0)
            {(0, 1): (0, 1, 2),
             (0, 2): (0, 2, 3),
             (1, 2): (0, 1, 2),
             (2, 0): (0, 1, 2),
             (2, 3): (0, 2, 3),
             (3, 0): (0, 2, 3)}

        """
        reference_polygon = self._reference.polygon(reference_label)

        outer_edges = [
            (vertex, (vertex + 1) % len(reference_polygon.vertices()))
            for vertex in range(len(reference_polygon.vertices()))
        ]
        inner_edges = reference_polygon.triangulation()
        inner_edges.extend([(w, v) for (v, w) in inner_edges])

        edges = outer_edges + inner_edges

        def triangle(edge):
            v, w = edge
            next_edges = [edge for edge in edges if edge[0] == w]
            previous_edges = [edge for edge in edges if edge[1] == v]

            next_vertices = [edge[1] for edge in next_edges]
            previous_vertices = [edge[0] for edge in previous_edges]

            other_vertex = set(next_vertices).intersection(set(previous_vertices))

            assert len(other_vertex) == 1

            other_vertex = next(iter(other_vertex))

            vertices = v, w, other_vertex
            while vertices[0] != min(vertices):
                vertices = vertices[1:] + vertices[:1]

            return vertices

        return {edge: triangle(edge) for edge in edges}

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.polygon((0, (0, 1, 2)))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])

        """
        reference_label, vertices = label
        reference_polygon = self._reference.polygon(reference_label)

        from flatsurf import Polygon

        return Polygon(
            vertices=[reference_polygon.vertex(v) for v in vertices],
            category=reference_polygon.category(),
        )

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.opposite_edge((0, (0, 1, 2)), 0)
            ((1, (0, 2, 3)), 1)

        """
        reference_label, vertices = label
        reference_polygon = self._reference.polygon(reference_label)

        if vertices[(edge + 1) % 3] == (vertices[edge] + 1) % len(
            reference_polygon.vertices()
        ):
            # This is an edge of the reference surface
            cross_reference_label, cross_reference_edge = self._reference.opposite_edge(
                reference_label, vertices[edge]
            )
            cross_reference_polygon = self._reference.polygon(cross_reference_label)
            cross_vertices = self._triangulation(cross_reference_label)[
                (
                    cross_reference_edge,
                    (cross_reference_edge + 1)
                    % len(cross_reference_polygon.vertices()),
                )
            ]

            cross_edge = cross_vertices.index(cross_reference_edge)

            return (cross_reference_label, cross_vertices), cross_edge

        # This is an edge that was added by the triangulation
        edge = (vertices[edge], vertices[(edge + 1) % 3])
        cross_edge = (edge[1], edge[0])
        cross_vertices = self._triangulation(reference_label)[cross_edge]
        return (reference_label, cross_vertices), cross_vertices.index(cross_edge[0])

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.triangulate()) == hash(S.triangulate())
            True

        """
        return hash(self._reference)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.triangulate() == S.triangulate()
            True

        """
        if not isinstance(other, LazyTriangulatedSurface):
            return False

        return self._reference == other._reference

    def labels(self):
        r"""
        Return the labels of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.labels()
            ((0, (0, 1, 2)), (1, (0, 2, 3)), (-1, (0, 2, 3)), (0, (0, 2, 3)), (1, (0, 1, 2)), (2, (0, 1, 2)), (-1, (0, 1, 2)), (-2, (0, 1, 2)), (2, (0, 2, 3)), (3, (0, 2, 3)),
             (-2, (0, 2, 3)), (-3, (0, 2, 3)), (3, (0, 1, 2)), (4, (0, 1, 2)), (-3, (0, 1, 2)), (-4, (0, 1, 2)), …)

        """

        class LazyLabels(Labels):
            def __contains__(self, label):
                reference_label, vertices = label
                if reference_label not in self._surface._reference.labels():
                    return False

                return (
                    vertices in self._surface._triangulation(reference_label).values()
                )

        return LazyLabels(self, finite=self._reference.is_finite_type())

    def __repr__(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S
            Triangulation of The infinite staircase

        """
        return f"Triangulation of {self._reference!r}"


class LazyMutableOrientedSimilaritySurface(MutableOrientedSimilaritySurface_base):
    r"""
    A helper surface for :class:`LazyDelaunayTriangulatedSurface`.

    A mutable wrapper of an (infinite) reference surface. When a polygon is not
    present in this wrapper yet, it is taken from the reference surface and can
    then be modified.

    .. NOTE::

        This surface does not implement the entire surface interface correctly.
        It just supports the operations in the way that they are necessary to
        make :class:`LazyDelaunayTriangulatedSurface` work.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()

        sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
        sage: T = LazyMutableOrientedSimilaritySurface(S)
        sage: p = T.polygon(0)
        sage: p
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
        sage: q = p * 2

        sage: S.replace_polygon(0, q)
        Traceback (most recent call last):
        ...
        AttributeError: '_InfiniteStaircase_with_category' object has no attribute 'replace_polygon'

        sage: T.replace_polygon(0, q)
        sage: T.polygon(0)
        Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)])

    """

    def __init__(self, surface, category=None):
        from flatsurf.geometry.categories import SimilaritySurfaces

        if surface not in SimilaritySurfaces().Oriented().WithoutBoundary():
            raise NotImplementedError("cannot handle surfaces with boundary yet")

        self._reference = surface

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        self._surface = MutableOrientedSimilaritySurface(surface.base_ring())

        super().__init__(surface.base_ring(), category=category or surface.category())

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        .. NOTE::

            This assumes that :meth:`glue` is never called to glue components.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.roots()
            (0,)

        """
        return self._reference.roots()

    def labels(self):
        r"""
        Return the labels of this surface which are just the labels of the
        underlying reference surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.labels()
            (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

        """
        return self._reference.labels()

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``True``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.is_mutable()
            True

        """
        return True

    def replace_polygon(self, label, polygon):
        r"""
        Swap out the polygon with the label ``label`` with ``polygon``.

        The polygons must have the same number of sides since gluings are kept.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.replace_polygon(0, T.polygon(0))

        """
        self._ensure_polygon(label)
        return self._surface.replace_polygon(label, polygon)

    def glue(self, x, y):
        r"""
        Glue the (label, edge) pair ``x`` with the pair ``y`` in this surface.

        This unglues any existing gluings of these edges.

        .. NOTE::

            After a sequence of such glue operations, no edges must be unglued.
            Otherwise, gluings get copied over from the underlying surface with
            confusing side effects.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.gluings()
            (((0, 0), (1, 2)), ((0, 1), (-1, 3)), ((0, 2), (1, 0)), ((0, 3), (-1, 1)), ((1, 0), (0, 2)), ((1, 1), (2, 3)), ((1, 2), (0, 0)), ((1, 3), (2, 1)), ((-1, 0), (-2, 2)),
             ((-1, 1), (0, 3)), ((-1, 2), (-2, 0)), ((-1, 3), (0, 1)), ((2, 0), (3, 2)), ((2, 1), (1, 3)), ((2, 2), (3, 0)), ((2, 3), (1, 1)), …)
            sage: T.glue((0, 0), (1, 0))
            sage: T.glue((1, 2), (0, 2))
            sage: T.gluings()
            (((0, 0), (1, 0)), ((0, 1), (-1, 3)), ((0, 2), (1, 2)), ((0, 3), (-1, 1)), ((1, 0), (0, 0)), ((1, 1), (2, 3)), ((1, 2), (0, 2)), ((1, 3), (2, 1)), ((-1, 0), (-2, 2)),
             ((-1, 1), (0, 3)), ((-1, 2), (-2, 0)), ((-1, 3), (0, 1)), ((2, 0), (3, 2)), ((2, 1), (1, 3)), ((2, 2), (3, 0)), ((2, 3), (1, 1)), …)

        """
        return self._surface.glue(x, y)

    def _ensure_gluings(self, label):
        r"""
        Make sure that the surface used to internally represent this surface
        has copied over all the gluings for ``label`` from the underlying
        surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T._ensure_polygon(0)
            sage: T._ensure_gluings(0)

        """
        self._ensure_polygon(label)
        for edge in range(len(self._surface.polygon(label).vertices())):
            cross = self._surface.opposite_edge(label, edge)
            if cross is None:
                cross_label, cross_edge = self._reference.opposite_edge(label, edge)
                self._ensure_polygon(cross_label)

                assert (
                    self._surface.opposite_edge(cross_label, cross_edge) is None
                ), "surface must not have a boundary"

                # Note that we cannot detect whether something has been
                # explicitly unglued. So we just reestablish any gluings of
                # this edge.
                self._surface.glue((label, edge), (cross_label, cross_edge))

    def _ensure_polygon(self, label):
        r"""
        Make sure that the surface used to internally represent this surface
        has copied over the polygon ``label`` from the underlying surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T._ensure_polygon(0)

        """
        if label not in self._surface.labels():
            self._surface.add_polygon(self._reference.polygon(label), label=label)

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.polygon(0)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        self._ensure_polygon(label)
        return self._surface.polygon(label)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.delaunay import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.opposite_edge(0, 0)
            (1, 2)

        """
        self._ensure_polygon(label)
        self._ensure_gluings(label)
        cross_label, cross_edge = self._surface.opposite_edge(label, edge)
        self._ensure_polygon(cross_label)
        self._ensure_gluings(cross_label)
        return cross_label, cross_edge


class LazyDelaunayTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    Delaunay triangulation of an (infinite type) surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
        sage: len(S.polygon(S.root()).vertices())
        3
        sage: TestSuite(S).run()  # long time (.8s)
        sage: S.is_delaunay_triangulated(limit=10)
        True

        sage: from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface
        sage: isinstance(S, LazyDelaunayTriangulatedSurface)
        True

    ::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: S = chamanara_surface(QQ(1/2))
        sage: m = matrix([[2,1],[1,1]])**4
        sage: S = (m*S).delaunay_triangulation()
        sage: TestSuite(S).run()  # long time (1s)
        sage: S.is_delaunay_triangulated(limit=10)
        True
        sage: TestSuite(S).run()  # long time (.5s)

        sage: from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface
        sage: isinstance(S, LazyDelaunayTriangulatedSurface)
        True

    """

    def __init__(self, similarity_surface, direction=None, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyDelaunayTriangulatedSurface; use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyDelaunayTriangulatedSurface()"
                )

        if similarity_surface.is_mutable():
            raise ValueError("surface must be immutable")

        if not similarity_surface.is_connected():
            raise NotImplementedError("surface must be connected")

        self._reference = similarity_surface

        # This surface will converge to the Delaunay Triangulation
        self._surface = LazyMutableOrientedSimilaritySurface(
            LazyTriangulatedSurface(similarity_surface)
        )

        self._direction = (self._surface.base_ring() ** 2)(direction or (0, 1))
        self._direction.set_immutable()

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels = set()

        # Triangulate the base polygon
        root = self._surface.root()

        # Certify the base polygon (or apply flips...)
        while not self._certify_or_improve(root):
            pass

        OrientedSimilaritySurface.__init__(
            self,
            self._surface.base_ring(),
            category=category or self._surface.category(),
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.is_mutable()
            False

        """
        return False

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.is_compact()
            False

        """
        return self._reference.is_compact()

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.roots()
            ((0, (0, 1, 2)),)

        """
        return self._surface.roots()

    @cached_method
    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.polygon((0, (0, 1, 2)))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])

        """
        if label not in self.labels():
            raise ValueError("no polygon with this label")

        if label not in self._certified_labels:
            for certified_label in self._walk():
                if label == certified_label:
                    assert label in self._certified_labels
                    break

        return self._surface.polygon(label)

    def _walk(self):
        visited = set()
        from collections import deque

        next = deque(
            [(self.root(), 0), (self.root(), 1), (self.root(), 2)],
        )
        while next:
            label, edge = next.popleft()
            if label in visited:
                continue

            yield label

            visited.add(label)

            for edge in range(3):
                next.append(self.opposite_edge(label, edge))

    @cached_method
    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.opposite_edge((0, (0, 1, 2)), 0)
            ((1, (0, 2, 3)), 1)

        """
        self.polygon(label)
        while True:
            cross_label, cross_edge = self._surface.opposite_edge(label, edge)
            if self._certify_or_improve(cross_label):
                break

        return self._surface.opposite_edge(label, edge)

    def _certify_or_improve(self, label):
        r"""
        This method attempts to develop the circumscribing disk about the polygon
        with label ``label`` into the surface.

        The method returns True if this is successful. In this case the label
        is added to the set _certified_labels. It returns False if it failed to
        develop the disk into the surface. (In this case the original polygon was
        not a Delaunay triangle.

        The algorithm divides any non-certified polygon in self._s it encounters
        into triangles. If it encounters a pair of triangles which need a diagonal
        flip then it does the flip.
        """
        if label in self._certified_labels:
            # Already certified.
            return True
        p = self._surface.polygon(label)
        assert len(p.vertices()) == 3

        c = p.circumscribing_circle()

        # Develop through each of the 3 edges:
        for e in range(3):
            edge_certified = False
            # This keeps track of a chain of polygons the disk develops through:
            edge_stack = []

            # We repeat this until we can verify that the portion of the circle
            # that passes through the edge e developes into the surface.
            while not edge_certified:
                if len(edge_stack) == 0:
                    # Start at the beginning with label l and edge e.
                    # The 3rd coordinate in the tuple represents what edge to develop
                    # through in the triangle opposite this edge.
                    edge_stack = [(label, e, 1, c)]
                ll, ee, step, cc = edge_stack[len(edge_stack) - 1]

                lll, eee = self._surface.opposite_edge(ll, ee)

                if lll not in self._certified_labels:
                    ppp = self._surface.polygon(lll)
                    assert len(ppp.vertices()) == 3

                    if self._surface._delaunay_edge_needs_flip(ll, ee):
                        # Perform the flip
                        self._surface.triangle_flip(
                            ll, ee, in_place=True, direction=self._direction
                        )

                        # If we touch the original polygon, then we return False.
                        if label == ll or label == lll:
                            return False
                        # We might have flipped a polygon from earlier in the chain
                        # In this case we need to trim the stack down so that we recheck
                        # that polygon.
                        for index, tup in enumerate(edge_stack):
                            if tup[0] == ll or tup[0] == lll:
                                edge_stack = edge_stack[:index]
                                break
                        # The following if statement makes sure that we check both subsequent edges of the
                        # polygon opposite the last edge listed in the stack.
                        if len(edge_stack) > 0:
                            ll, ee, step, cc = edge_stack.pop()
                            edge_stack.append((ll, ee, 1, cc))
                        continue

                    # If we reach here then we know that no flip was needed.
                    ccc = self._surface.edge_transformation(ll, ee) * cc

                    # Check if the disk passes through the next edge in the chain.
                    lp = ccc.line_segment_position(
                        ppp.vertex((eee + step) % 3), ppp.vertex((eee + step + 1) % 3)
                    )
                    if lp == 1:
                        # disk passes through edge and opposite polygon is not certified.
                        edge_stack.append((lll, (eee + step) % 3, 1, ccc))
                        continue

                    # We reach this point if the disk doesn't pass through the edge eee+step of polygon lll.

                # Either lll is already certified or the disk didn't pass
                # through edge (lll,eee+step)

                # Trim off unnecessary edges off the stack.
                # prune_count=1
                ll, ee, step, cc = edge_stack.pop()
                if step == 1:
                    # if we have just done step 1 (one edge), move on to checking
                    # the next edge.
                    edge_stack.append((ll, ee, 2, cc))
                # if we have pruned an edge, continue to look at pruning in the same way.
                while step == 2 and len(edge_stack) > 0:
                    ll, ee, step, cc = edge_stack.pop()
                    # prune_count= prune_count+1
                    if step == 1:
                        edge_stack.append((ll, ee, 2, cc))
                if len(edge_stack) == 0:
                    # We're done with this edge
                    edge_certified = True
        self._certified_labels.add(label)
        return True

    def labels(self):
        r"""
        Return the labels of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.labels()
            ((0, (0, 1, 2)), (1, (0, 2, 3)), (-1, (0, 2, 3)), (0, (0, 2, 3)), (1, (0, 1, 2)), (2, (0, 1, 2)), (-1, (0, 1, 2)), (-2, (0, 1, 2)), (2, (0, 2, 3)), (3, (0, 2, 3)),
             (-2, (0, 2, 3)), (-3, (0, 2, 3)), (3, (0, 1, 2)), (4, (0, 1, 2)), (-3, (0, 1, 2)), (-4, (0, 1, 2)), …)

        """
        return self._surface.labels()

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.delaunay_triangulation()) == hash(S.delaunay_triangulation())
            True

        """
        return hash((self._reference, self._direction))

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.delaunay_triangulation() == S.delaunay_triangulation()
            True

        """
        if not isinstance(other, LazyDelaunayTriangulatedSurface):
            return False

        return (
            self._reference == other._reference and self._direction == other._direction
        )

    def __repr__(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()

        """
        return f"Delaunay triangulation of {self._reference!r}"


class LazyDelaunaySurface(OrientedSimilaritySurface):
    r"""
    Delaunay cell decomposition of an (infinite type) surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()
        sage: m = matrix([[2, 1], [1, 1]])
        sage: S = (m * S).delaunay_decomposition()

        sage: S.polygon(S.root())
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: S.is_delaunay_decomposed(limit=10)  # long time (.7s)
        True

        sage: TestSuite(S).run()  # long time (2s)

        sage: from flatsurf.geometry.delaunay import LazyDelaunaySurface
        sage: isinstance(S, LazyDelaunaySurface)
        True

    ::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: S = chamanara_surface(QQ(1/2))
        sage: m = matrix([[3, 4], [-4, 3]]) * matrix([[4, 0],[0, 1/4]])
        sage: S = (m * S).delaunay_decomposition()
        sage: S.is_delaunay_decomposed(limit=10)  # long time (.5s)
        True

        sage: TestSuite(S).run()  # long time (1.5s)

        sage: from flatsurf.geometry.delaunay import LazyDelaunaySurface
        sage: isinstance(S, LazyDelaunaySurface)
        True

    """

    def __init__(self, similarity_surface, direction=None, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyDelaunaySurface; use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyDelaunaySurface()"
                )

        if similarity_surface.is_mutable():
            raise ValueError("surface must be immutable.")

        self._reference = similarity_surface

        self._delaunay_triangulation = LazyDelaunayTriangulatedSurface(
            self._reference, direction=direction, relabel=relabel
        )

        super().__init__(
            similarity_surface.base_ring(),
            category=category or similarity_surface.category(),
        )

    @cached_method
    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.polygon((0, (0, 1, 2)))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        if label not in self._delaunay_triangulation.labels():
            raise ValueError("no polygon with this label")

        cell, edges = self._cell(label)

        if label != self._label(cell):
            raise ValueError("no polygon with this label")

        edges = [
            self._delaunay_triangulation.polygon(edge[0]).edge(edge[1])
            for edge in edges
        ]

        from flatsurf import Polygon

        return Polygon(edges=edges)

    @cached_method
    def _label(self, cell):
        r"""
        Return a canonical label for the Delaunay cell that is made up by the
        Delaunay triangles ``cell``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S._label(frozenset({
            ....:     ((0, (0, 1, 2))),
            ....:     ((0, (0, 2, 3)))}))
            (0, (0, 1, 2))

        """
        for label in self._delaunay_triangulation.labels():
            if label in cell:
                return label

    @cached_method
    def _normalize_label(self, label):
        r"""
        Return a canonical label for the Delaunay cell that contains the
        Delaunay triangle ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S._normalize_label((0, (0, 1, 2)))
            (0, (0, 1, 2))
            sage: S._normalize_label((0, (0, 2, 3)))
            (0, (0, 1, 2))

        """
        cell, _ = self._cell(label)
        return self._label(cell)

    @cached_method
    def _cell(self, label):
        r"""
        Return the labels of the Delaunay triangles that contain the Delaunay
        triangle ``label`` together with the interior edges in that cell.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()

        This cell (a square) is formed by two triangles that form a cylinder,
        i.e., the two triangles are glued at two of their edges::

            sage: S._cell((0, (0, 1, 2)))
            (frozenset({(0, (0, 1, 2)), (0, (0, 2, 3))}),
             [((0, (0, 1, 2)), 0),
              ((0, (0, 1, 2)), 1),
              ((0, (0, 2, 3)), 1),
              ((0, (0, 2, 3)), 2)])

        """
        edges = []
        cell = set()
        explore = [(label, 2), (label, 1), (label, 0)]

        while explore:
            triangle, edge = explore.pop()

            cell.add(triangle)

            delaunay = self._delaunay_triangulation._delaunay_edge_needs_join(
                triangle, edge
            )

            if not delaunay:
                edges.append((triangle, edge))
                continue

            cross_triangle, cross_edge = self._delaunay_triangulation.opposite_edge(
                triangle, edge
            )

            for shift in [2, 1]:
                next_triangle, next_edge = cross_triangle, (cross_edge + shift) % 3

                if (next_triangle, next_edge) in edges:
                    raise NotImplementedError
                if (next_triangle, next_edge) in explore:
                    raise NotImplementedError

                explore.append((next_triangle, next_edge))

        cell = frozenset(cell)
        normalized_label = self._label(cell)
        if normalized_label != label:
            return self._cell(normalized_label)

        return cell, edges

    @cached_method
    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.opposite_edge((0, (0, 1, 2)), 0)
            ((1, (0, 2, 3)), 2)

        """
        if label not in self._delaunay_triangulation.labels():
            raise ValueError

        cell, edges = self._cell(label)

        if label != self._label(cell):
            raise ValueError

        edge = edges[edge]

        cross_triangle, cross_edge = self._delaunay_triangulation.opposite_edge(*edge)

        cross_cell, cross_edges = self._cell(cross_triangle)
        cross_label = self._label(cross_cell)
        return cross_label, cross_edges.index((cross_triangle, cross_edge))

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.roots()
            ((0, (0, 1, 2)),)

        """
        return self._delaunay_triangulation.roots()

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.is_compact()
            False

        """
        return self._reference.is_compact()

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.is_mutable()
            False

        """
        return False

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.delaunay_decomposition()) == hash(S.delaunay_decomposition())
            True

        """
        return hash(self._delaunay_triangulation)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.delaunay import LazyDelaunaySurface
            sage: S = translation_surfaces.infinite_staircase()
            sage: m = matrix([[2, 1], [1, 1]])
            sage: S = LazyDelaunaySurface(m*S)
            sage: S == S
            True

        """
        if not isinstance(other, LazyDelaunaySurface):
            return False

        return self._delaunay_triangulation == other._delaunay_triangulation

    def __repr__(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.delaunay_decomposition()
            Delaunay cell decomposition of The infinite staircase

        """
        return f"Delaunay cell decomposition of {self._reference!r}"
