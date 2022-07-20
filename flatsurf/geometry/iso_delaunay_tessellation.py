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
        self._surface = self._surface_original.delaunay_triangulation()
        self._hyperbolic_plane = HyperbolicPlane(surface.base_ring())
        self._faces = Graph([[self.root()], []])
        self._faces.set_vertex(self.root(), self._surface)

    def _repr_(self):
        return f"IsoDelaunay Tessellation of {self._surface_original}"

    def explore(self, limit=None, vertex=None, edge=None):
        r"""
        Explore the dual graph of the IsoDelaunay tessellation up to the combinatorial ``limit`` where you first cross ``edge``.
        When ``edge`` is ``None``, explore all edges of the IsoDelaunay regions containing `i`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()

        """
        limit = limit or 3  # TODO: explore until finding FD
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
        for source_triangulation_edge in source_triangulation.edge_iterator():
            half_plane = self._half_plane(source_triangulation, source_triangulation_edge)
            if half_plane is None:
                continue
            if half_plane.boundary() == edge:
                edges_to_flip.append(source_triangulation_edge)
        for edge_to_flip in edges_to_flip:
            target_triangulation.flip_edge(edge_to_flip)
        
        target = self._hyperbolic_plane.polygon(self._iso_delaunay_region(target_triangulation))

        assert -edge in target.edges()

        is_new_idr = target not in self._faces
        if is_new_idr:
            self._faces.add_vertex(target)

        if not self._faces.has_edge(source, target):
            self._faces.add_edge(source, target)

        for edge in target.edges():
            self.explore(limit=limit-1, edge=edge)

    def is_vertex(self, translation_surface):
        r"""
        Return whether this is a vertex.
        """
        raise NotImplementedError

    def root(self):
        r"""
        Return the vertex from which we started to build the IsoDelaunay tessellation.
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
        if iso_delaunay_region.dimension() < 2:
            epsilon = 1/2
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
        raise NotImplementedError

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
        t1 = surface.polygon(0)
        print(t1)
        v0, v1, _ = t0.edges()
        w0, w1, _ = t1.edges()

        assert False

    def _point_to_matrix(self, point):
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return matrix(2, [1, x, 0, y])

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
