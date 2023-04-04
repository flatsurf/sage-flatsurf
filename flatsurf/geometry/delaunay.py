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

from flatsurf.geometry.surface import Surface


class LazyTriangulatedSurface(Surface):
    r"""
    Surface class used to triangulate an infinite surface.

    EXAMPLES:

    Example with relabel=False::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyTriangulatedSurface(s,relabel=False))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run()

    Example with relabel=True::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyTriangulatedSurface(s,relabel=True))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run()
    """

    def __init__(self, similarity_surface, relabel=True):
        if similarity_surface.is_mutable():
            raise ValueError("Surface must be immutable.")

        # This surface will converge to the Delaunay Triangulation
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, mutable=True)

        Surface.__init__(
            self,
            self._s.base_ring(),
            self._s.base_label(),
            mutable=False,
            finite=self._s.is_finite(),
        )

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        polygon = self._s.polygon(lab)
        if polygon.num_edges() > 3:
            self._s.triangulate(in_place=True, label=lab)
            return self._s.polygon(lab)
        else:
            return polygon

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        pp, ee = self._s.opposite_edge(p, e)
        polygon = self._s.polygon(pp)
        if polygon.num_edges() > 3:
            self.polygon(pp)
            return self._s.opposite_edge(p, e)
        else:
            return (pp, ee)

    def __eq__(self, other):
        if isinstance(other, LazyTriangulatedSurface):
            if self._s == other._s:
                return True

        return super().__eq__(other)


class LazyDelaunayTriangulatedSurface(Surface):
    r"""
    Surface class used to find a Delaunay triangulation of an infinite surface.

    EXAMPLES:

    Example with relabel=False::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyDelaunayTriangulatedSurface(s,relabel=False))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run()
        sage: ss.is_delaunay_triangulated(limit=100)
        True

    Example with relabel=True::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyDelaunayTriangulatedSurface(s,relabel=True))
        sage: ss.polygon(0).num_edges()
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

    def _setup_direction(self, direction):
        # Our Delaunay will respect the provided direction.
        if direction is None:
            self._direction = self._s.vector_space()(
                (self._s.base_ring().zero(), self._s.base_ring().one())
            )
        else:
            self._direction = self._ss.vector_space()(direction)

    def __init__(self, similarity_surface, direction=None, relabel=True):
        r"""
        Construct a lazy Delaunay triangulation of the provided similarity_surface.
        """
        if similarity_surface.underlying_surface().is_mutable():
            raise ValueError("Surface must be immutable.")

        # This surface will converge to the Delaunay Triangulation
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, mutable=True)

        self._setup_direction(direction)

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels = set()

        # Triangulate the base polygon
        base_label = self._s.base_label()
        self._s.triangulate(in_place=True, label=base_label)

        # Certify the base polygon (or apply flips...)
        while not self._certify_or_improve(base_label):
            pass

        Surface.__init__(
            self,
            self._s.base_ring(),
            base_label,
            finite=self._s.is_finite(),
            mutable=False,
        )

    def polygon(self, label):
        if label in self._certified_labels:
            return self._s.polygon(label)
        else:
            raise ValueError(
                "Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface."
            )

    def opposite_edge(self, label, edge):
        if label in self._certified_labels:
            ll, ee = self._s.opposite_edge(label, edge)
            if ll in self._certified_labels:
                return ll, ee
            while not self._certify_or_improve(ll):
                ll, ee = self._s.opposite_edge(label, edge)
            return self._s.opposite_edge(label, edge)
        else:
            raise ValueError(
                "Asked for an edge of a polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface."
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

    def __eq__(self, other):
        if isinstance(other, LazyDelaunayTriangulatedSurface):
            if self._s == other._s:
                return True

        return super().__eq__(other)


class LazyDelaunaySurface(LazyDelaunayTriangulatedSurface):
    # We just inherit to use some methods.

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
        sage: ss=TranslationSurface(LazyDelaunaySurface(m*s,relabel=False))
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
        sage: ss=TranslationSurface(LazyDelaunaySurface(m*s))
        sage: ss.is_delaunay_decomposed(limit=100)
        True
        sage: TestSuite(ss).run()
    """

    def __init__(self, similarity_surface, direction=None, relabel=True):
        r"""
        Construct a lazy Delaunay triangulation of the provided similarity_surface.

        """
        if similarity_surface.underlying_surface().is_mutable():
            raise ValueError("Surface must be immutable.")

        # This surface will converge to the Delaunay Decomposition
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, mutable=True)

        self._setup_direction(direction)

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels = set()
        self._decomposition_certified_labels = set()

        base_label = self._s.base_label()

        # We will now try to get the base_polygon.
        # Certify the base polygon (or apply flips...)
        while not self._certify_or_improve(base_label):
            pass
        self._certify_decomposition(base_label)

        Surface.__init__(
            self,
            self._s.base_ring(),
            base_label,
            finite=self._s.is_finite(),
            mutable=False,
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

    def __eq__(self, other):
        if isinstance(other, LazyDelaunaySurface):
            if self._s == other._s:
                return True

        return super().__eq__(other)
