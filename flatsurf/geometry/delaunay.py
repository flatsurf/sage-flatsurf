r"""
This file contains classes implementing Surface which are used for
triangulating, Delaunay triangulating, and Delaunay decomposing infinite
surfaces.
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

from flatsurf.geometry.surface import MutableOrientedSimilaritySurface_base, OrientedSimilaritySurface, Labels


class LazyTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    Surface class used to triangulate an infinite surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.delaunay import LazyTriangulatedSurface
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=LazyTriangulatedSurface(s)
        sage: ss.polygon(ss.root()).num_edges()
        3
        sage: TestSuite(ss).run()
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
        return False

    def is_compact(self):
        return self._reference.is_compact()

    def roots(self):
        return tuple((reference_label, self._triangulation(reference_label)[(0, 1)]) for reference_label in self._reference.roots())

    def _triangulation(self, reference_label):
        reference_polygon = self._reference.polygon(reference_label)

        outer_edges = [
            (vertex, (vertex + 1) % reference_polygon.num_edges())
            for vertex in range(reference_polygon.num_edges())
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
        reference_label, vertices = label
        reference_polygon = self._reference.polygon(reference_label)

        from flatsurf import polygon
        return polygon(
            vertices=[reference_polygon.vertex(v) for v in vertices],
            category=reference_polygon.category()
        )

    def opposite_edge(self, label, edge):
        reference_label, vertices = label
        reference_polygon = self._reference.polygon(reference_label)

        if (
            vertices[(edge + 1) % 3]
            == (vertices[edge] + 1) % reference_polygon.num_edges()
        ):
            # This is an edge of the reference surface
            cross_reference_label, cross_reference_edge = self._reference.opposite_edge(
                reference_label, vertices[edge]
            )
            cross_reference_polygon = self._reference.polygon(cross_reference_label)
            cross_vertices = self._triangulation(cross_reference_label)[
                (
                    cross_reference_edge,
                    (cross_reference_edge + 1) % cross_reference_polygon.num_edges(),
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
        return hash(self._reference)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of inequality.

        EXAMPLES::

            sage: from flatsurf.geometry.delaunay import LazyTriangulatedSurface
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S = LazyTriangulatedSurface(S)
            sage: S == S
            True

        """
        if not isinstance(other, LazyTriangulatedSurface):
            return False

        return self._reference == other._reference

    def labels(self):
        return LazyTriangulatedSurface.LazyLabels(self)

    class LazyLabels(Labels):
        def __contains__(self, label):
            reference_label, vertices = label
            if reference_label not in self._surface._reference.labels():
                return False

            return vertices in self._surface._triangulation(reference_label).values()


class LazyMutableOrientedSimilaritySurface(MutableOrientedSimilaritySurface_base):
    def __init__(self, surface, category=None):
        self._reference = surface

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        self._surface = MutableOrientedSimilaritySurface(surface.base_ring())

        super().__init__(surface.base_ring(), category=category or surface.category())

    def roots(self):
        return self._reference.roots()

    def labels(self):
        return self._reference.labels()

    def is_mutable(self):
        return True

    def replace_polygon(self, label, polygon):
        return self._surface.replace_polygon(label, polygon)

    def glue(self, x, y):
        return self._surface.glue(x, y)

    def _ensure_gluings(self, label):
        assert label in self._surface.labels()
        for edge in range(self._surface.polygon(label).num_edges()):
            cross = self._surface.opposite_edge(label, edge)
            if cross is None:
                cross_label, cross_edge = self._reference.opposite_edge(label, edge)
                self._ensure_polygon(cross_label)
                self._surface.glue((label, edge), (cross_label, cross_edge))

    def _ensure_polygon(self, label):
        if label not in self._surface.labels():
            self._surface.add_polygon(self._reference.polygon(label), label=label)

    def polygon(self, label):
        self._ensure_polygon(label)
        return self._surface.polygon(label)

    def opposite_edge(self, label, edge):
        self._ensure_polygon(label)
        self._ensure_gluings(label)
        cross_label, cross_edge = self._surface.opposite_edge(label, edge)
        self._ensure_polygon(cross_label)
        self._ensure_gluings(cross_label)
        return cross_label, cross_edge


class LazyDelaunayTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    Surface class used to find a Delaunay triangulation of an infinite surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=LazyDelaunayTriangulatedSurface(s)
        sage: ss.polygon(ss.root()).num_edges()
        3
        sage: TestSuite(ss).run()
        sage: ss.is_delaunay_triangulated(limit=100)
        True

    Chamanara example::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: s=chamanara_surface(QQ(1/2))
        sage: m=matrix([[2,1],[1,1]])**4
        sage: ss=(m*s).delaunay_triangulation()
        sage: TestSuite(ss).run()
        sage: ss.is_delaunay_triangulated(limit=100)
        True
        sage: TestSuite(ss).run()
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
        self._surface = LazyMutableOrientedSimilaritySurface(LazyTriangulatedSurface(similarity_surface))

        self._direction = (self._surface.base_ring()**2)(direction or (0, 1))
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
        return False

    def is_compact(self):
        return self._reference.is_compact()

    def roots(self):
        return self._surface.roots()

    @cached_method
    def polygon(self, label):
        if label not in self.labels():
            raise ValueError
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
        assert p.num_edges() == 3

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
                    assert ppp.num_edges() == 3

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
        return self._surface.labels()

    def __hash__(self):
        return hash((self._reference, self._direction))

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of inequality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface
            sage: S = translation_surfaces.infinite_staircase()
            sage: S = LazyDelaunayTriangulatedSurface(S)
            sage: S == S
            True

        """
        if not isinstance(other, LazyDelaunayTriangulatedSurface):
            return False

        return self._reference == other._reference and self._direction == other._direction


class LazyDelaunaySurface(OrientedSimilaritySurface):
    r"""
    This is an implementation of Surface. It takes a surface (typically
    infinite) from the constructor. This class represents the Delaunay
    decomposition of this surface. We compute this decomposition lazily so that
    it works for infinite surfaces.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.delaunay import LazyDelaunaySurface
        sage: s=translation_surfaces.infinite_staircase()
        sage: m=matrix([[2,1],[1,1]])
        sage: ss=LazyDelaunaySurface(m*s)
        sage: ss.polygon(ss.root())
        polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
        sage: ss.is_delaunay_decomposed(limit=100)
        True
        sage: TestSuite(ss).run()

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: s=chamanara_surface(QQ(1/2))
        sage: m=matrix([[3,4],[-4,3]])*matrix([[4,0],[0,1/4]])
        sage: ss=LazyDelaunaySurface(m*s)
        sage: ss.is_delaunay_decomposed(limit=100)
        True
        sage: TestSuite(ss).run()
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
            raise ValueError("Surface must be immutable.")

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
        if label not in self._delaunay_triangulation.labels():
            raise ValueError

        cell, edges = self._cell(label)

        if label != self._label(cell):
            raise ValueError

        edges = [
            self._delaunay_triangulation.polygon(edge[0]).edge(edge[1])
            for edge in edges
        ]

        from flatsurf import polygon
        return polygon(edges=edges, category=self._delaunay_triangulation.polygon(label).parent())

    @cached_method
    def _label(self, cell):
        for label in self._delaunay_triangulation.labels():
            if label in cell:
                return label

    @cached_method
    def _normalize_label(self, label):
        cell, _ = self._cell(label)
        return self._label(cell)

    @cached_method
    def _cell(self, label):
        edges = []
        cell = set()
        explore = [(label, 2), (label, 1), (label, 0)]

        while explore:
            triangle, edge = explore.pop()

            cell.add(triangle)

            delaunay = self._delaunay_triangulation._delaunay_edge_needs_join(triangle, edge)

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
        return self._delaunay_triangulation.roots()

    def is_compact(self):
        return self._reference.is_compact()

    def is_mutable(self):
        return False

    def __hash__(self):
        return hash(self._delaunay_triangulation)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of inequality.

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
