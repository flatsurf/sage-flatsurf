r"""
EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
    sage: s = translation_surfaces.veech_double_n_gon(5)
    sage: IsoDelaunayTessellation(s)
    IsoDelaunay Tessellation of TranslationSurface built from 2 polygons

TODO: Write something smart here.
TODO: Implement this surface in sage-flatsurf.

::

    sage: def prym3(w, h, t, e):
    ....:     from flatsurf import Surface_dict, TranslationSurface, polygons
    ....:     k.<rtD> = QuadraticField(e**2 + 8 * w * h)
    ....:     lmbd = (e + rtD)/2
    ....:
    ....:     square = polygons.square(lmbd, field=k)
    ....:     parallelogram_top = polygons((w - lmbd, 0), (lmbd, 0), (t, h), (-lmbd, 0), (lmbd - w, 0), (-t, -h), field=k)
    ....:     parallelogram_bottom = polygons((lmbd, 0), (w - lmbd, 0), (t, h), (lmbd - w, 0), (-lmbd, 0), (-t, -h), field=k)
    ....:
    ....:     surface = Surface_dict(base_ring=k)
    ....:
    ....:     surface.add_polygon(parallelogram_bottom, label=0)
    ....:     surface.add_polygon(square, label=1)
    ....:     surface.add_polygon(parallelogram_top, label=2)
    ....:
    ....:     surface.change_base_label(1)
    ....:
    ....:     surface.change_edge_gluing(0, 0, 2, 3)
    ....:     surface.change_edge_gluing(0, 1, 0, 3)
    ....:     surface.change_edge_gluing(0, 2, 0, 5)
    ....:     surface.change_edge_gluing(0, 4, 1, 0)
    ....:     surface.change_edge_gluing(1, 1, 1, 3)
    ....:     surface.change_edge_gluing(1, 2, 2, 1)
    ....:     surface.change_edge_gluing(2, 0, 2, 4)
    ....:     surface.change_edge_gluing(2, 2, 2, 5)
    ....:
    ....:     return TranslationSurface(surface)
    ....:
    sage: s = prym3(5, 1, 0, 0)
    sage: t = s.delaunay_triangulation(in_place=False)

    sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
    sage: idt = IsoDelaunayTessellation(t)
    sage: idt.explore()
    sage: idt.plot()
    Graphics object consisting of 110 graphics primitives

REFERENCES:

TODO: Write a text citing some references.

.. [JB2009] \J. Bowman, "Flat Structures and Complex Structures in Teichmüller
            Theory", PhD Thesis, https://hdl.handle.net/1813/13979

"""
# *********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Sam Freedman
#                      2022 Julian Rüth
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
# *********************************************************************

from sage.structure.parent import Parent
from flatsurf.geometry.hyperbolic import HyperbolicPlane
from sage.misc.cachefunc import cached_method


class IsoDelaunayTessellation(Parent):
    def __init__(self, surface):
        from sage.all import Graph
        self._surface_original = surface

        self._hyperbolic_plane = HyperbolicPlane(surface.base_ring())

        self._surface = self._nondegenerate_delaunay_triangulation(surface)
        self._surface.set_immutable()

        self._dual_graph = Graph(multiedges=True, loops=True)

        self._ensure_dual_graph_vertex(self.root(), self._surface, self.root().edges()[0])

    def _repr_(self):
        return f"IsoDelaunay Tessellation of {self._surface_original}"

    def explore(self, limit=None, tessellation_face=None):
        r"""
        Explore the dual graph of the IsoDelaunay tessellation up to the combinatorial ``limit`` where you start from ``vertex`` and then first cross ``edge``.
        When ``vertex`` is ``None``, start exploring from the vertex bounded by ``edge``.
        When ``edge`` is ``None``, explore all edges of the IsoDelaunay regions bounding ``vertex``.
        When both ``vertex`` and ``edge`` are ``None``, start exploring from a vertex that contains the point ``i``.


        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()

        """
        from sage.all import oo

        limit = oo if limit is None else limit
        if limit <= 0:
            return

        if tessellation_face is None:
            tessellation_face = self.root()

        if tessellation_face not in self._dual_graph:
            raise ValueError("tessellation_face must be a polygon of the explored tessellation")

        queued_tessellation_faces = set()

        from collections import deque
        queue = deque([(0, tessellation_face)])
        queued_tessellation_faces.add(tessellation_face)

        while queue:
            distance, tessellation_face = queue.popleft()

            for tessellation_edge in tessellation_face.edges():
                other_tessellation_face = self._explore(tessellation_face, tessellation_edge)
                if distance + 1 < limit and other_tessellation_face not in queued_tessellation_faces:
                    queued_tessellation_faces.add(other_tessellation_face)
                    queue.append((distance + 1, other_tessellation_face))

    def _explore(self, tessellation_face, tessellation_edge):
        assert tessellation_face.dimension() == 2
        assert tessellation_edge.dimension() == 1

        cross_tessellation_face = self._cross(tessellation_face, tessellation_edge)

        # nothing to do if we have already explored across this polygon edge
        if cross_tessellation_face is not None:
            return cross_tessellation_face

        _, source_triangulation = self._dual_graph.get_vertex(tessellation_face)
        target_triangulation = source_triangulation.copy()

        while True:
            for triangulation_edge in target_triangulation.edge_iterator():
                half_plane = self._half_plane(target_triangulation, triangulation_edge)
                if half_plane is None:
                    continue
                if half_plane.boundary() == tessellation_edge.geodesic():
                    target_triangulation = target_triangulation.triangle_flip(
                        *triangulation_edge)
                    break

                # glue triangles across triangulation_edge that flips to get mock Delaunay cells
                # solve for self-isomorphism on level of mock cells?
            else:
                break

        cross_tessellation_face = self._hyperbolic_plane.polygon(
            self._iso_delaunay_region(target_triangulation))

        assert -tessellation_edge in cross_tessellation_face.edges(), f"edge {-tessellation_edge} is not in the polygon {cross_tessellation_face} after crossing {tessellation_edge} from {tessellation_face}"

        target_triangulation.set_immutable()
        cross_tessellation_face, cross_tessellation_edge, is_new = self._ensure_dual_graph_vertex(cross_tessellation_face, target_triangulation, -tessellation_edge)

        self._dual_graph.add_edge(tessellation_face, cross_tessellation_face, label={tessellation_edge, cross_tessellation_edge})

        return cross_tessellation_face

    def _cross(self, tessellation_face, tessellation_edge):
        mod, _ = self._dual_graph.get_vertex(tessellation_face)

        if mod is not None:
            edges = list(tessellation_face.edges())
            edge = set(edges[edges.index(tessellation_edge) % mod::mod])
        else:
            edge = {tessellation_edge}

        for v, w, edges in self._dual_graph.edges(tessellation_face, labels=True):
            if any(e in edges for e in edge):
                if v == tessellation_face:
                    return w
                if w == tessellation_face:
                    return v
                assert False

        return None

    def insert_orbifold_points(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.mcmullen_genus2_prototype(1, 1, 0, -1)
            sage: t = s.delaunay_triangulation(in_place=False)
            sage: idt = IsoDelaunayTessellation(t)
            sage: z = idt._hyperbolic_plane(i)
            sage: idt.explore()

            sage: idt._dual_graph.vertices()
            [{2*l*(x^2 + y^2) + (-8*l + 4)*x ≥ 0} ∩ {(8*l - 4)*x - 4*l + 2 ≤ 0} ∩ {x ≥ 0}]
            sage: idt.insert_orbifold_points()
            sage: idt._dual_graph.vertices()
            [{2*l*(x^2 + y^2) + (-8*l + 4)*x ≥ 0} ∩ {(8*l - 4)*x - 4*l + 2 ≤ 0} ∩ {x ≥ 0} ∪ {I}]

        """
        # TODO: Should this mutate the tessellation or create a copy instead?
        for tessellation_face in self._dual_graph.vertices():
            for source_tessellation_face, target_tessellation_face, tessellation_edges in self._dual_graph.edges(tessellation_face, labels=True):
                # crossing edge of the polygon cycles back to the very edge in the
                # polygon, so there is an orbifold point on that edge.
                # We patch the polygon by inserting a marked point.
                if len(tessellation_edges) == 1:
                    assert source_tessellation_face == target_tessellation_face
                    tessellation_face, _ = self._insert_orbifold_point(tessellation_face, next(iter(tessellation_edges)).midpoint())

        # TODO: Insert orbifold points in the interior of a polygon, i.e., the ones detected with isomorphism()

    def _insert_orbifold_point(self, tessellation_face, point):
        r"""
        Insert ``point`` as a marked point on an edge of the polygon ``vertex``
        (and update the graph representing the explored part of the fundamental
        domain.)

        Return the new polygon and the edges adjacent to point.
        """
        for tessellation_edge in tessellation_face.edges():
            if point in tessellation_edge:
                break
        else:
            assert False

        tessellation_face_with_marked_vertices = tessellation_face.parent().polygon(
            tessellation_face.half_spaces(),
            marked_vertices=tuple(tessellation_face.vertices()) + (point,))

        assert tessellation_face_with_marked_vertices != tessellation_face

        self._dual_graph.add_vertex(tessellation_face_with_marked_vertices)
        self._dual_graph.set_vertex(tessellation_face_with_marked_vertices, self._dual_graph.get_vertex(tessellation_face))
        for source_tessellation_face, target_tessellation_face, tessellation_edges in list(self._dual_graph.edges(tessellation_face, labels=True)):
            self._dual_graph.delete_edge(source_tessellation_face, target_tessellation_face, tessellation_edges)
            if source_tessellation_face == tessellation_face:
                source_tessellation_face = tessellation_face_with_marked_vertices
            if target_tessellation_face == tessellation_face:
                target_tessellation_face = tessellation_face_with_marked_vertices

            self._dual_graph.add_edge(source_tessellation_face, target_tessellation_face, label=tessellation_edges)

        self._dual_graph.delete_vertex(tessellation_face)

        for e in tessellation_face_with_marked_vertices.edges():
            if e.start() == tessellation_edge.start():
                edge_before_marked_point = e
                break
        else:
            assert False

        for e in tessellation_face_with_marked_vertices.edges():
            if e.end() == edge_before_marked_point.end():
                edge_after_marked_point = e
                break
        else:
            assert False

        return tessellation_face_with_marked_vertices, (edge_before_marked_point, edge_after_marked_point)

    @cached_method
    def _to_pyflatsurf(self, triangulation):
        # TODO: Clear this cache when explore() is finished.
        from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
        return to_pyflatsurf(triangulation)

    def _ensure_dual_graph_vertex(self, tessellation_face, surface, tessellation_edge):
        r"""
        Return vertex and edge of hyperbolic polygon TODO
        """
        for tessellation_face_ in self._dual_graph:
            if tessellation_face_ == tessellation_face:
                return tessellation_face_, tessellation_edge, False

        surface_ = self._to_pyflatsurf(surface)

        def filter_matrix(a, b, c, d):
            # TODO fix interface for isomorphisms
            from pyeantic import RealEmbeddedNumberField
            k = RealEmbeddedNumberField(a.parent())
            a, b, c, d = k(a), k(b), k(c), k(d)
            k = k.number_field
            a, b, c, d = k(a), k(b), k(c), k(d)

            if (a, b, c, d) in isomorphisms or a * d - b * c == -1:
                return False
            isomorphisms[-1] = (a, b, c, d)
            return True

        isomorphisms = [(1, 0, 0, 1), (-1, 0, 0, -1)]

        while True:
            isomorphisms.append(())
            if not surface_.isomorphism(surface_, filter_matrix=filter_matrix).has_value():
                isomorphisms.pop()
                break

        if set(isomorphisms) != {(1, 0, 0, 1), (-1, 0, 0, -1)}:
            assert len(isomorphisms) % 2 == 0
            order = len(isomorphisms) // 2
            assert len(tessellation_face.edges()) % order == 0
            mod = len(tessellation_face.edges()) // order

            # add new vertex
            self._dual_graph.add_vertex(tessellation_face)
            assert surface is not None
            self._dual_graph.set_vertex(tessellation_face, (mod, surface))
            return tessellation_face, tessellation_edge, True

        for tessellation_face_ in self._dual_graph:
            isomorphism = None

            def capture_matrix(a, b, c, d):
                if a * d - b * c != 1:
                    return False

                from pyeantic import RealEmbeddedNumberField
                k = RealEmbeddedNumberField(a.parent())
                a, b, c, d = k(a), k(b), k(c), k(d)
                k = k.number_field
                a, b, c, d = k(a), k(b), k(c), k(d)

                nonlocal isomorphism
                isomorphism = (a, b, c, d)
                return True

            if len(tessellation_face_.edges()) != len(tessellation_face.edges()):
                # TODO: This optimization is not correct after insert_orbifold_points() has been called.
                continue

            surface__ = self._to_pyflatsurf(self._dual_graph.get_vertex(tessellation_face_)[1])
            if surface_.isomorphism(surface__, filter_matrix=capture_matrix).has_value():
                assert isomorphism is not None
                a, b, c, d = isomorphism
                from sage.all import matrix
                mob = matrix(2, [a, -b, -c, d])
                image_edge = tessellation_edge.apply_isometry(mob, model='half_plane')

                assert image_edge in tessellation_face_.edges()
                return tessellation_face_, image_edge, False

        # add new vertex
        self._dual_graph.add_vertex(tessellation_face)
        assert surface is not None
        self._dual_graph.set_vertex(tessellation_face, (None, surface))
        return tessellation_face, tessellation_edge, True

    def is_vertex(self, translation_surface):
        r"""
        Return whether this is a vertex.
        """
        raise NotImplementedError

    def root(self):
        r"""
        Return the tessellation face from which we started to build the IsoDelaunay tessellation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.root()
            {(6*a + 8)*(x^2 + y^2) + (-20*a - 28)*x + 14*a + 20 ≥ 0} ∩ {(2*a + 3)*(x^2 + y^2) + (-4*a - 6)*x - 2*a - 3 ≤ 0} ∩ {(4*a + 6)*(x^2 + y^2) - 4*a - 6 ≥ 0}
        """
        from sage.all import I
        return self.face(I)

    def point(self, surface):
        r"""
        Return the point in the :class:`HyperbolicPlane` corresponding to ``surface``.

        INPUT:

        - ``surface`` -- A surface in the SL(2, R)-orbit of the defining surface
        """
        raise NotImplementedError

    def edge(self, translation_surface):
        r"""
        Return the unoriented tessellation edge this translation surface is on.
        """
        raise NotImplementedError

    @classmethod
    def _face(cls, surface, point):
        A = cls._point_to_matrix(point)
        A_T = surface.apply_matrix(
            A, in_place=False).delaunay_triangulation(in_place=False)
        T = A_T.apply_matrix(~A, in_place=False)
        half_planes = cls._iso_delaunay_region(T)
        iso_delaunay_region = point.parent().polygon(half_planes)

        if iso_delaunay_region.dimension() < 2:
            return None

        return iso_delaunay_region

    def face(self, point):
        r"""
        Return a tessellation face containing the hyperbolic ``point``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.face(i)
            {(6*a + 8)*(x^2 + y^2) + (-20*a - 28)*x + 14*a + 20 ≥ 0} ∩ {(2*a + 3)*(x^2 + y^2) + (-4*a - 6)*x - 2*a - 3 ≤ 0} ∩ {(4*a + 6)*(x^2 + y^2) - 4*a - 6 ≥ 0}

        """
        point = self._hyperbolic_plane(point)
        x, y = point.coordinates()

        from sage.all import QQ
        shift = QQ(1)/2

        while True:
            shifted = point.parent().point(x + shift, y, model="half_plane")
            face = self._face(self._surface, shifted)
            if point in face:
                return face
            shift /= 2

    def surface(self, point):
        r"""
        Return a translation surface corresponding to this ``point`` in the
        hyperbolic plane.
        """
        # This should be shared with the IDR code.
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return self._surface.apply_matrix(matrix([[1, x], [0, y]]), in_place=False)

    def fundamental_domain(self):
        r"""
        Return the fundamental domain as a polygon with edge pairings.
        """
        # The result would be a fundamental domain of the group generated by the
        # symmetries discovered so far.
        raise NotImplementedError

    def plot(self):
        # TODO: Why is plotting so slow?
        return sum(idr.plot() for idr in self._dual_graph)

    def polygon(self, vertex_or_edge):
        r"""
        Return the polygon obtained as the union of the
        triangles bounded by this edge and its reverse /
        by the edges adjacent to this vertex.
        """
        raise NotImplementedError

    def geodesics(self, vertex):
        r"""
        Return the geodesics through ``vertex``.

        These are the geodesics underlying the :meth:`segments`, with end points removed.
        """

        raise NotImplementedError

    @classmethod
    def _point_to_matrix(cls, point):
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return matrix(2, [1, x, 0, y])

    @classmethod
    def _iso_delaunay_region(cls, triangulation):
        return [half_plane for edge in triangulation.edge_iterator()
                if (half_plane := cls._half_plane(triangulation, edge)) is not None]

    @classmethod
    def _half_plane(cls, surface, edge):
        r"""
        Build the halfplane associated to ``edge``.
        If the hinge is not convex, return ``None``.
        """
        # check if hinge is convex
        if not surface.triangle_flip(*edge, test=True):
            return None

        index_triangle, index_edge = edge
        index_opposite_triangle, index_opposite_edge = surface.opposite_edge(
            edge)
        v0 = surface.polygon(index_opposite_triangle).edge(
            (index_opposite_edge + 1) % 3)
        v1 = surface.polygon(index_triangle).edge(index_edge)
        v2 = -surface.polygon(index_triangle).edge((index_edge - 1) % 3)

        x0, y0 = v0
        x1, y1 = v1
        x2, y2 = v2

        a = x0 * y1 * y2 * (y2 - y1) + x1 * y0 * y2 * \
            (y0 - y2) + x2 * y0 * y1 * (y1 - y0)
        b = x0 * y0 * (x1 * y2 - x2 * y1) + x1 * y1 * \
            (x2 * y0 - x0 * y2) + x2 * y2 * (x0 * y1 - x1 * y0)
        c = x0 * x1 * y2 * (x0 - x1) + x1 * x2 * y0 * \
            (x1 - x2) + x0 * x2 * y1 * (x2 - x0)

        H = HyperbolicPlane(surface.base_ring())
        return H.geodesic(a, 2 * b, c, model="half_plane").left_half_space()

    @classmethod
    def _nondegenerate_delaunay_triangulation(cls, surface):
        r"""
        Return a Delaunay triangulation for `_surface` whose IDR is 2-dimensional
        """
        # Let M_i be starting surf with arbitrary DT producing potentially degen IDR
        # Consider T_e * M_i until it has a nondegenerate IDR
        # let M_{i+e} be the Delaunay triangulated surface with a nondegen IDR
        # then T_{-e} * M_{i + e} is a new DT of M_i with a nondegen IDR

        # If M_i has a degenerate IDR, do flips until not the case
        # e.g. only 8 / 250 of the triangulations of the regular octagon lead to a nondegen IDR

        surface = surface.delaunay_triangulation(in_place=False)

        H = HyperbolicPlane(surface.base_ring())

        shift = 0
        while True:
            shifted = H.point(shift, 1, model="half_plane")
            face = cls._face(surface, shifted)

            from sage.all import I
            if face is not None and I in face:
                break

            if shift == 0:
                from sage.all import QQ
                shift = QQ(1)

            shift /= 2

        perturbation = cls._point_to_matrix(shifted)
        return surface.apply_matrix(perturbation, in_place=False).delaunay_triangulation(in_place=False).apply_matrix(~perturbation, in_place=False)

    def vertices(self):
        r"""
        Return the vertices of the completed fundamental domain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()
            sage: len(idt.vertices())
            3

        """
        # TODO: Check that the fundamental domain has been computed.

        half_edges = set((tessellation_face, tessellation_edge) for tessellation_face in self._dual_graph.vertices() for tessellation_edge in tessellation_face.edges())

        vertices = []

        while half_edges:
            vertex = [next(iter(half_edges))]

            while True:
                (tessellation_face, tessellation_edge) = vertex[-1]

                assert tessellation_edge in tessellation_face.edges()

                previous_tessellation_edge = tessellation_face.edges()[tessellation_face.edges().index(tessellation_edge) - 1]

                # TODO: Merge this with the code in _cross maybe.
                for source_tessellation_face, target_tessellation_face, tessellation_edges in self._dual_graph.edges(tessellation_face, labels=True):

                    if previous_tessellation_edge in tessellation_edges:
                        if len(tessellation_edges) == 1:
                            cross_tessellation_edge = previous_tessellation_edge
                        else:
                            tessellation_edges = list(tessellation_edges)
                            if previous_tessellation_edge == tessellation_edges[0]:
                                cross_tessellation_edge = tessellation_edges[1]
                            else:
                                cross_tessellation_edge = tessellation_edges[0]           
                    else:
                        continue

                    cross_tessellation_edge = (target_tessellation_face if source_tessellation_face == tessellation_face else source_tessellation_face, cross_tessellation_edge)

                    vertex.append(cross_tessellation_edge)
                    half_edges.remove(cross_tessellation_edge)
                    break
                else:
                    assert False, f"{previous_tessellation_edge} of {tessellation_face} not glued to another edge"

                if vertex[0] == vertex[-1]:
                    vertex.pop()
                    vertices.append(vertex)
                    break

        return vertices

    def topological_euler_characteristic(self):
        # return V - E + F of the fundamental domain
        v = len(self.vertices())
        e = len(self._dual_graph.edges())
        f = len(self._dual_graph.vertices())

        return v - e + f - len(self.cusps())

    def orbifold_euler_characteristic(self):
        r"""
        chi_top = V - E + F in ZZ
        chi_orb = chi_top + \sum_{orbifold points} (-1 + 1/n) in QQ

        eg for H^2/SL(2, Z)
        chi_top = 2 - 2 + 1 = euler char of a punctured sphere

        chi_orb = 1 + (-1 + 1/2) + (-1 + 1/3) = -1 + 5/6 = -1/6
        """
        raise NotImplementedError

    def genus(self):
        r"""
        chi_top = 2 - n - 2g ==> g = (2 - n - chi_top)/2
        """
        pass

    def cusps(self):
        r"""
        Return the punctures of the quotient H^2/Gamma
        """
        return [vertex_equivalence_class for vertex_equivalence_class in self.vertices()
                       if vertex_equivalence_class[0][1].start().is_ideal()]

    def orbifold_points(self, order=None):
        r"""
        Return the set of orbifold points, i.e., the fixed points of finite-order rotations in Gamma.

        When ``order = k``, return only the orbifold points with total angle ``2pi/k``.
        """
        # TODO: implement by developing in a circle around vertex
        pass
