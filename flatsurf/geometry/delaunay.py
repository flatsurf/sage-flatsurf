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

from flatsurf.geometry.surface import OrientedSimilaritySurface, Labels


class LazyTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    Surface class used to triangulate an infinite surface.

    EXAMPLES:

    Example with relabel=False::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=LazyTriangulatedSurface(s,relabel=False)
        sage: ss.polygon(ss.base_label()).num_edges()
        3
        sage: TestSuite(ss).run()

    Example with relabel=True::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=LazyTriangulatedSurface(s,relabel=True)
        sage: ss.polygon(ss.base_label()).num_edges()
        3
        sage: TestSuite(ss).run()
    """

    def __init__(self, similarity_surface, relabel=True, category=None):
        if similarity_surface.is_mutable():
            raise ValueError("Surface must be immutable.")

        self._reference = similarity_surface

        OrientedSimilaritySurface.__init__(
            self,
            similarity_surface.base_ring(),
            category=category or self._reference.category())

    def is_mutable(self):
        return False

    def base_label(self):
        reference_label = self._reference.base_label()
        return (reference_label, self._triangulation(reference_label)[(0, 1)])

    def _triangulation(self, reference_label):
        reference_polygon = self._reference.polygon(reference_label)

        outer_edges = [(vertex, (vertex + 1) % reference_polygon.num_edges()) for vertex in range(reference_polygon.num_edges())]
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
        return reference_polygon.parent()(vertices=[reference_polygon.vertex(v) for v in vertices])

    def opposite_edge(self, label, edge):
        reference_label, vertices = label
        reference_polygon = self._reference.polygon(reference_label)

        if vertices[(edge + 1) % 3] == (vertices[edge] + 1) % reference_polygon.num_edges():
            # This is an edge of the reference surface
            cross_reference_label, cross_reference_edge = self._reference.opposite_edge(reference_label, vertices[edge])
            cross_reference_polygon = self._reference.polygon(cross_reference_label)
            cross_vertices = self._triangulation(cross_reference_label)[(cross_reference_edge, (cross_reference_edge + 1) % cross_reference_polygon.num_edges())]

            cross_edge = cross_vertices.index(cross_reference_edge)

            return (cross_reference_label, cross_vertices), cross_edge

        # This is an edge that was added by the triangulation
        edge = (vertices[edge], vertices[(edge + 1) % 3])
        cross_edge = (edge[1], edge[0])
        cross_vertices = self._triangulation(reference_label)[cross_edge]
        return (reference_label, cross_vertices), cross_vertices.index(cross_edge[0])

    def _cache_key(self):
        return (type(self), self._reference)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.delaunay import LazyTriangulatedSurface
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S = LazyTriangulatedSurface(S)
            sage: S == S
            True

        """
        if isinstance(other, LazyTriangulatedSurface):
            if self._reference == other._reference:
                return True

        return super().__eq__(other)

    def labels(self):
        return LazyTriangulatedSurface.LazyLabels(self)

    class LazyLabels(Labels):
        def __contains__(self, label):
            reference_label, vertices = label
            if reference_label not in self._surface._reference.labels():
                return False

            return vertices in self._surface._triangulation(reference_label).values()


class LazyMutableSurface(OrientedSimilaritySurface):
    def __init__(self, surface, category=None):
        self._reference = surface

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
        self._surface = MutableOrientedSimilaritySurface(surface.base_ring())

        super().__init__(surface.base_ring(), category=category or surface.category())

    def base_label(self):
        return self._reference.base_label()

    def labels(self):
        return self._reference.labels()

    def is_mutable(self):
        return True

    def change_polygon(self, label, polygon):
        return self._surface.change_polygon(label, polygon)

    def change_edge_gluing(self, label0, edge0, label1, edge1):
        return self._surface.change_edge_gluing(label0, edge0, label1, edge1)

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

    EXAMPLES:

    Example with relabel=False::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=LazyDelaunayTriangulatedSurface(s,relabel=False)
        sage: ss.polygon(ss.base_label()).num_edges()
        3
        sage: TestSuite(ss).run()
        sage: ss.is_delaunay_triangulated(limit=100)
        True

    Example with relabel=True::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=LazyDelaunayTriangulatedSurface(s,relabel=True)
        sage: ss.polygon(ss.base_label()).num_edges()
        3
        sage: TestSuite(ss).run()
        sage: ss.is_delaunay_triangulated(limit=100)
        True

    Chamanara example::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.chamanara import *
        sage: s=chamanara_surface(QQ(1/2))
        sage: m=matrix([[2,1],[1,1]])**4
        sage: ss=(m*s).delaunay_triangulation()
        sage: TestSuite(ss).run()
        sage: ss.is_delaunay_triangulated(limit=100)
        True
        sage: TestSuite(ss).run()
    """

    def __init__(self, similarity_surface, direction=None, relabel=True, category=None):
        if similarity_surface.is_mutable():
            raise ValueError("Surface must be immutable.")

        self._reference = similarity_surface

        # This surface will converge to the Delaunay Triangulation
        self._surface = LazyMutableSurface(LazyTriangulatedSurface(similarity_surface))

        self._direction = self._surface.vector_space()(direction or (0, 1))

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels = set()

        # Triangulate the base polygon
        base_label = self._surface.base_label()

        # Certify the base polygon (or apply flips...)
        while not self._certify_or_improve(base_label):
            pass

        OrientedSimilaritySurface.__init__(
            self,
            self._surface.base_ring(),
            category=category or self._surface.category()
        )

    def is_mutable(self):
        return False

    def base_label(self):
        return self._surface.base_label()

    def polygon(self, label):
        if label not in self._surface.labels():
            raise ValueError
        if label not in self._certified_labels:
            for certified_label in self.labels():
                if label == certified_label:
                    assert label in self._certified_labels
                    break

        return self._surface.polygon(label)

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

                    if self._surface._edge_needs_flip(ll, ee):
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

    def __hash__(self):
        return super().__hash__()

    def _cache_key(self):
        return (type(self), self._reference, self._direction)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface
            sage: S = translation_surfaces.infinite_staircase()
            sage: S = LazyDelaunayTriangulatedSurface(S)
            sage: S == S
            True

        """
        if isinstance(other, LazyDelaunayTriangulatedSurface):
            if self._reference == other._reference and self._direction == other._direction:
                return True

        return super().__eq__(other)


class LazyDelaunaySurface(OrientedSimilaritySurface):
    r"""
    This is an implementation of Surface. It takes a surface (typically
    infinite) from the constructor. This class represents the Delaunay
    decomposition of this surface. We compute this decomposition lazily so that
    it works for infinite surfaces.

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: m=matrix([[2,1],[1,1]])
        sage: ss=LazyDelaunaySurface(m*s,relabel=False)
        sage: ss.polygon(ss.base_label())
        Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
        sage: ss.is_delaunay_decomposed(limit=100)
        True
        sage: TestSuite(ss).run()

        sage: from flatsurf import *
        sage: from flatsurf.geometry.chamanara import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=chamanara_surface(QQ(1/2))
        sage: m=matrix([[3,4],[-4,3]])*matrix([[4,0],[0,1/4]])
        sage: ss=LazyDelaunaySurface(m*s)
        sage: ss.is_delaunay_decomposed(limit=100)
        True
        sage: TestSuite(ss).run()
    """

    def __init__(self, similarity_surface, direction=None, relabel=True, category=None):
        r"""
        Construct a lazy Delaunay triangulation of the provided similarity_surface.

        """
        if similarity_surface.underlying_surface().is_mutable():
            raise ValueError("Surface must be immutable.")

        self._reference = similarity_surface

        # This surface will converge to the Delaunay Decomposition
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, mutable=True)

        self._direction = self._s.vector_space()(direction or (0, 1))

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels = set()
        self._decomposition_certified_labels = set()

        base_label = self._s.base_label()

        # We will now try to get the base_polygon.
        # Certify the base polygon (or apply flips...)
        while not self._certify_or_improve(base_label):
            pass
        self._certify_decomposition(base_label)

        OrientedSimilaritySurface.__init__(
            self,
            self._s.base_ring(),
            category=category or similarity_surface.category()
        )

    def _certify_decomposition(self, label):
        if label in self._decomposition_certified_labels:
            return
        assert label in self._certified_labels
        changed = True
        while changed:
            changed = False
            p = self._s.polygon(label)
            for e in range(p.num_edges()):
                ll, ee = self._s.opposite_edge(label, e)
                while not self._certify_or_improve(ll):
                    ll, ee = self._s.opposite_edge(label, e)
                if self._s._edge_needs_join(label, e):
                    # ll should not have already been certified!
                    assert ll not in self._decomposition_certified_labels
                    self._s.join_polygons(label, e, in_place=True)
                    changed = True
                    break
        self._decomposition_certified_labels.add(label)

    def polygon(self, label):
        if label in self._decomposition_certified_labels:
            return self._s.polygon(label)
        else:
            raise ValueError(
                "Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface."
            )

    def opposite_edge(self, label, edge):
        if label in self._decomposition_certified_labels:
            ll, ee = self._s.opposite_edge(label, edge)
            if ll in self._decomposition_certified_labels:
                return ll, ee
            self._certify_decomposition(ll)
            return self._s.opposite_edge(label, edge)
        else:
            raise ValueError(
                "Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface."
            )

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
        p = self._s.polygon(label)
        if p.num_edges() > 3:
            # not triangulated!
            self._s.triangulate(in_place=True, label=label)
            p = self._s.polygon(label)
            # Made major changes to the polygon with label l.
            return False
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

                lll, eee = self._s.opposite_edge(ll, ee)

                if lll not in self._certified_labels:
                    ppp = self._s.polygon(lll)
                    if ppp.num_edges() > 3:
                        # not triangulated!
                        self._s.triangulate(in_place=True, label=lll)
                        lll, eee = self._s.opposite_edge(ll, ee)
                        ppp = self._s.polygon(lll)
                    # now ppp is a triangle

                    if self._s._edge_needs_flip(ll, ee):
                        # Perform the flip
                        self._s.triangle_flip(
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
                    ccc = self._s.edge_transformation(ll, ee) * cc

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

    def base_label(self):
        return self._s.base_label()

    def __hash__(self):
        return super().__hash__()

    def _cache_key(self):
        return (type(self), self._reference, self._direction)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.delaunay import LazyDelaunaySurface
            sage: S = translation_surfaces.infinite_staircase()
            sage: m = matrix([[2, 1], [1, 1]])
            sage: S = LazyDelaunaySurface(m*S)
            sage: S == S
            True

        """
        if isinstance(other, LazyDelaunaySurface):
            if self._reference == other._reference and self._direction == other._direction:
                return True

        return super().__eq__(other)

    def is_mutable(self):
        return False
