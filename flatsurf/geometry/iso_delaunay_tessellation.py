r"""
EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
    sage: s = translation_surfaces.veech_double_n_gon(5)
    sage: IsoDelaunayTessellation(s)
    IsoDelaunay Tessellation of TranslationSurface built from 2 polygons

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
        self._faces = DiGraph([[self.root()], []], multiedges=True, loops=True)
        self._faces.set_vertex(self.root(), self._surface)

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
        limit = 3 if limit is None else limit  # TODO: explore until finding FD
        if limit <= 0:
            return

        if edge is None:
            vertex = vertex or self.root()
            for edge in vertex.edges():
                self.explore(limit=limit, vertex=vertex, edge=edge)
            return

        source = vertex or self.face(edge)
        source_triangulation = self._faces.get_vertex(source)
        target_triangulation = source_triangulation.copy()
        edges_to_flip = []
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
            else:
                break

        target = self._hyperbolic_plane.polygon(
            self._iso_delaunay_region(target_triangulation))

        assert -edge in target.edges()

        kind, vertices = self._classify_idr(target, target_triangulation)

        if kind == "OLD":
            assert len(vertices) == 1
            self._faces.add_edge(source, vertices[0], label=(edge, -edge))
            self._faces.add_edge(vertices[0], source, label=(-edge, edge))
        elif kind == "ISOMORPHIC":
            assert len(vertices) == 1
            self._faces.add_edge(source, vertices[0], label=(edge, -edge))
            self._faces.add_edge(vertices[0], source, label=(-edge, edge))
        elif kind == "SYMMETRIC":
            assert len(vertices) >= 3
            raise NotImplementedError
        elif kind == "NEW":
            self._faces.add_vertex(target)
            self._faces.set_vertex(target, target_triangulation)
            self._faces.add_edge(source, target, label=(edge, -edge))
            self._faces.add_edge(target, source, label=(-edge, edge))

            for edge in target.edges():
                self.explore(limit=limit-1, vertex=target, edge=edge)

        else:
            raise ValueError(f"kind {kind} was unexpected")

    def _classify_idr(self, target_polygon, target_triangulation):
        r"""
        Classify whether ``target`` is identical to a vertex in ``self._faces``, or is isomorphic to a vertex in ``self._faces``, or has a self-isomorphism, or is a new IDR.
        """
        for vertex in self._faces:
            if vertex == target_polygon:
                return "OLD", [vertex]

        from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
        target_triangulation = to_pyflatsurf(target_triangulation)

        def filter_matrix(a, b, c, d):
            from pyeantic import RealEmbeddedNumberField
            k = RealEmbeddedNumberField(a.parent())
            a, b, c, d = k(a), k(b), k(c), k(d)
            l = k.number_field
            a, b, c, d = l(a), l(b), l(c), l(d)

            if (a, b, c, d) in isomorphisms or a * d - b * c == -1:
                return False
            isomorphisms[-1] = (a, b, c, d)
            return True

        isomorphisms = []

        while True:
            isomorphisms.append(())
            if not target_triangulation.isomorphism(target_triangulation, filter_matrix=filter_matrix).has_value():
                isomorphisms.pop()
                break
        print(isomorphisms)

        # TODO get "primitive rotation"

        for vertex in self._faces:
            triangulation = to_pyflatsurf(self._faces.get_vertex(vertex))
            if triangulation.isomorphism(target_triangulation, filter_matrix=lambda a, b, c, d: a*d - b*c == 1).has_value():
                return "ISOMORPHIC", [vertex]

        return "NEW", None

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
                if (half_plane :=self._half_plane(triangulation, edge)) is not None]

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
