r"""
EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
    sage: s = translation_surfaces.veech_double_n_gon(5)
    sage: IsoDelaunayTessellation(s)
    IsoDelaunay Tessellation of TranslationSurface built from 2 polygons

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


class IsoDelaunayTessellation(Parent):
    r"""
    A HyperbolicPlaneTessellation + more data
    """

    def __init__(self, surface):
        from sage.all import Graph
        self._surface_original = surface
        self._hyperbolic_plane = HyperbolicPlane(surface.base_ring())
        self._surface = surface.delaunay_triangulation()
        self._surface = self._delaunay_triangulation(surface)
        self._faces = Graph(multiedges=True, loops=True)
        self._ensure_vertex(self.root(), self._surface, self.root().edges()[0])

    def _repr_(self):
        return f"IsoDelaunay Tessellation of {self._surface_original}"

    def explore(self, limit=None, vertex=None, edge=None):
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
        # TODO: distinguish between triangulation edges, hyperbolic polygon edges, graph edges, ... in variable naming.
        # We could try to adopt bowman notation: he calls our _faces "tiles"
        # and uses "tessellation edge" and "tessellation vertex" for the
        # subsets of H.

        from sage.all import oo

        limit = oo if limit is None else limit
        if limit <= 0:
            return

        if vertex is None:
            if edge is None:
                vertex = self.root()
            else:
                for vertex in self._faces.vertices():
                    if edge in vertex.edges():
                        break
                else:
                    raise ValueError("edge must be an edge of a polygon in the explored tesselation")

        if vertex not in self._faces:
            raise ValueError("vertex must be a polygon of the explored tesselation")

        if edge is None:
            for edge in vertex.edges():
                self.explore(limit=limit, vertex=vertex, edge=edge)
            return

        if edge not in vertex.edges():
            raise ValueError("edge must be an edge of the vertex polygon")

        # nothing to do if we have already explored across this polygon edge
        for _, _, edges in self._faces.edges(vertex, labels=True):
            if edge in edges:
                return

        source_triangulation = self._faces.get_vertex(vertex)
        target_triangulation = source_triangulation.copy()

        while True:
            for triangulation_edge in target_triangulation.edge_iterator():
                half_plane = self._half_plane(
                    target_triangulation, triangulation_edge)
                if half_plane is None:
                    continue
                if half_plane.boundary() == edge.geodesic():
                    target_triangulation = target_triangulation.triangle_flip(
                        *triangulation_edge)
                    break

                # glue triangles across triangulation_edge that flips to get mock Delaunay cells
                # solve for self-isomorphism on level of mock cells?
            else:
                break

        target = self._hyperbolic_plane.polygon(
            self._iso_delaunay_region(target_triangulation))

        assert -edge in target.edges(), f"edge {-edge} is not in the polygon {target} after crossing {edge} from {vertex}"

        target, polygon_edge, is_new = self._ensure_vertex(target, target_triangulation, -edge)

        self._faces.add_edge(vertex, target, label={edge, polygon_edge})

        if is_new:
            for edge in target.edges():
                self.explore(limit=limit-1, vertex=target, edge=edge)

    def insert_orbifold_points(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.mcmullen_genus2_prototype(1, 1, 0, -1)
            sage: t = s.delaunay_triangulation()
            sage: idt = IsoDelaunayTessellation(t)
            sage: z = idt._hyperbolic_plane(i)
            sage: idt.explore()

            sage: idt._faces.vertices()
            [{2*l*(x^2 + y^2) + (-8*l + 4)*x ≥ 0} ∩ {(8*l - 4)*x - 4*l + 2 ≤ 0} ∩ {x ≥ 0}]
            sage: idt.insert_orbifold_points()
            sage: idt._faces.vertices()
            [{2*l*(x^2 + y^2) + (-8*l + 4)*x ≥ 0} ∩ {(8*l - 4)*x - 4*l + 2 ≤ 0} ∩ {x ≥ 0} ∪ {I}]

        """
        # TODO: Should this mutate the tesselation or create a copy instead?
        for vertex in self._faces.vertices():
            for source, target, edges in self._faces.edges(vertex, labels=True):
                # crossing edge of the polygon cycles back to the very edge in the
                # polygon, so there is an orbifold point on that edge.
                # We patch the polygon by inserting a marked point.
                if len(edges) == 1:
                    assert source == target
                    vertex, _ = self._insert_orbifold_point(vertex, next(iter(edges)).midpoint())

        # TODO: Insert orbifold points in the interior of a polygon, i.e., the ones detected with isomorphism()

    def _insert_orbifold_point(self, vertex, point):
        r"""
        Insert ``point`` as a marked point on an edge of the polygon ``vertex``
        (and update the graph representing the explored part of the fundamental
        domain.)

        Return the new polygon and the edges adjacent to point.
        """
        for edge in vertex.edges():
            if point in edge:
                break
        else:
            assert False

        polygon = vertex.parent().polygon(
            vertex.half_spaces(),
            marked_vertices=tuple(vertex.vertices()) + (point,))

        assert polygon != vertex

        self._faces.add_vertex(polygon)
        self._faces.set_vertex(polygon, self._faces.get_vertex(vertex))
        for source, target, label in list(self._faces.edges(vertex, labels=True)):
            self._faces.delete_edge(source, target, label)
            if source == vertex:
                source = polygon
            if target == vertex:
                target = polygon

            self._faces.add_edge(source, target, label=label)

        self._faces.delete_vertex(vertex)

        for e in polygon.edges():
            if e.start() == edge.start():
                edge = e
                break
        else:
            assert False

        for e in polygon.edges():
            if e.end() == edge.end():
                polygon_edge = e
                break
        else:
            assert False

        return polygon, (edge, polygon_edge)

    def _ensure_vertex(self, target_polygon, target_triangulation, target_polygon_edge):
        r"""
        Return vertex and edge of hyperbolic polygon TODO
        """
        for vertex in self._faces:
            if vertex == target_polygon:
                return vertex, target_polygon_edge, False

        from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
        flat_target_triangulation = to_pyflatsurf(target_triangulation)

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

        isomorphisms = []

        while True:
            isomorphisms.append(())
            if not flat_target_triangulation.isomorphism(flat_target_triangulation, filter_matrix=filter_matrix).has_value():
                isomorphisms.pop()
                break

        # TODO get "primitive rotation"
        if set(isomorphisms) != {(1, 0, 0, 1), (-1, 0, 0, -1)}:
            raise NotImplementedError("subdivide IDR not implemented")

        # TODO vertex is an awful name
        for vertex in self._faces:
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

            triangulation = to_pyflatsurf(self._faces.get_vertex(vertex))
            if flat_target_triangulation.isomorphism(triangulation, filter_matrix=capture_matrix).has_value():
                assert isomorphism is not None
                a, b, c, d = isomorphism
                from sage.all import matrix
                mob = matrix(2, [a, -b, -c, d])
                image_edge = target_polygon_edge.apply_isometry(mob, model='half_plane')

                assert image_edge in vertex.edges()
                return vertex, image_edge, False

        # add new vertex
        self._faces.add_vertex(target_polygon)
        self._faces.set_vertex(target_polygon, target_triangulation)
        return target_polygon, target_polygon_edge, True

    def is_vertex(self, translation_surface):
        r"""
        Return whether this is a vertex.
        """
        raise NotImplementedError

    def root(self):
        r"""
        Return the vertex from which we started to build the IsoDelaunay tessellation.

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

    def face(self, point_or_edge):
        r"""
        Return a tessellation face this translation surface is in.

        If ``edge`` is oriented, return the bounding face.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: i = idt._hyperbolic_plane(i) # TODO Automatically coerce
            sage: idt.face(i)
            {(6*a + 8)*(x^2 + y^2) + (-20*a - 28)*x + 14*a + 20 ≥ 0} ∩ {(2*a + 3)*(x^2 + y^2) + (-4*a - 6)*x - 2*a - 3 ≤ 0} ∩ {(4*a + 6)*(x^2 + y^2) - 4*a - 6 ≥ 0}

        """

        if point_or_edge in self._hyperbolic_plane:
            point_or_edge = self._hyperbolic_plane(point_or_edge)
            if point_or_edge.dimension() == 1:
                z = point_or_edge.an_element()
            else:
                z = point_or_edge

        else:  # point_or_edge is a surface
            raise NotImplementedError

        A = self._point_to_matrix(z)
        A_T = self._surface.apply_matrix(
            A, in_place=False).delaunay_triangulation(in_place=False)
        T = A_T.apply_matrix(~A, in_place=False)
        half_planes = self._iso_delaunay_region(T)
        iso_delaunay_region = self._hyperbolic_plane.polygon(half_planes)

        # Let M_i be starting surf with arbitrary DT producing potentially degen IDR
        # Consider T_e * M_i until it has a nondegenerate IDR
        # let M_{i+e} be the Delaunay triangulated surface with a nondegen IDR
        # then T_{-e} * M_{i + e} is a new DT of M_i with a nondegen IDR

        # If M_i has a degenerate IDR, do flips until not the case
        # e.g. only 8 / 250 of the triangulations of the regular octagon lead to a nondegen IDR
        if iso_delaunay_region.dimension() < 2:
            from sage.all import QQ
            epsilon = QQ(1)/2
            while True:
                x, y = z.coordinates()
                z1 = self._hyperbolic_plane.point(
                    x + epsilon, y, model="half_plane")
                face1 = self.face(z1)
                if z in face1:
                    return face1
                epsilon /= 2

        return iso_delaunay_region

    def surface(self, point):
        r"""
        Return a translation surface corresponding to this ``point`` in the
        hyperbolic plane.
        """
        # This should be shared with the IDR code.
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return self._surface.apply_matrix(matrix([[1, x], [0, y]], in_place=False))

    def fundamental_domain(self):
        r"""
        Return the fundamental domain as a polygon with edge pairings.
        """
        # The result would be a fundamental domain of the group generated by the
        # symmetries discovered so far.
        raise NotImplementedError

    def plot(self):
        return sum(idr.plot() for idr in self._faces)

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

    def _surface_to_point(self, surface):
        t0 = self._surface.polygon(0)
        t1 = surface.p
        v0, v1, _ = t0.edges()
        w0, w1, _ = t1.edges()

        assert False

    def _point_to_matrix(self, point):
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return matrix(2, [1, x, 0, y])

    def _point_to_surface(self, point):
        m = self._point_to_matrix(point)
        return self._surface.apply_matrix(m)

    def _iso_delaunay_region(self, triangulation):
        return [half_plane for edge in triangulation.edge_iterator()
                if (half_plane := self._half_plane(triangulation, edge)) is not None]

    def _half_plane(self, triangulation, edge):
        r"""
        Build the halfplane associated to ``edge``.
        If the hinge is not convex, return ``None``.
        """
        # check if hinge is convex
        if not triangulation.triangle_flip(*edge, test=True):
            return None

        index_triangle, index_edge = edge
        index_opposite_triangle, index_opposite_edge = triangulation.opposite_edge(
            edge)
        v0 = triangulation.polygon(index_opposite_triangle).edge(
            (index_opposite_edge + 1) % 3)
        v1 = triangulation.polygon(index_triangle).edge(index_edge)
        v2 = -triangulation.polygon(index_triangle).edge((index_edge - 1) % 3)

        x0, y0 = v0
        x1, y1 = v1
        x2, y2 = v2

        a = x0 * y1 * y2 * (y2 - y1) + x1 * y0 * y2 * \
            (y0 - y2) + x2 * y0 * y1 * (y1 - y0)
        b = x0 * y0 * (x1 * y2 - x2 * y1) + x1 * y1 * \
            (x2 * y0 - x0 * y2) + x2 * y2 * (x0 * y1 - x1 * y0)
        c = x0 * x1 * y2 * (x0 - x1) + x1 * x2 * y0 * \
            (x1 - x2) + x0 * x2 * y1 * (x2 - x0)

        return self._hyperbolic_plane.geodesic(a, 2 * b, c, model="half_plane").left_half_space()

    def _delaunay_triangulation(self, _surface):
        r"""
        Return a Delaunay triangulation for `_surface` whose IDR is 2-dimensional
        """
        from sage.all import QQ, I

        epsilon = QQ(1/2)
        while epsilon > QQ(1/32):
            I1 = self._hyperbolic_plane.point(epsilon, 1, model="half_plane")
            face1 = self.face(I1)
            # if z1 in face1.interior() and I in face1: TODO add interior()
            if I in face1:
                break
            epsilon /= 2

        nondegenerate_triangulation = self._point_to_surface(
            I1).delaunay_triangulation()
        perturbation = self._point_to_matrix(I1)
        return nondegenerate_triangulation.apply_matrix(perturbation.inverse())
