r"""
EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.spine_tessellation import SpineTessellation
    sage: s = translation_surfaces.veech_double_n_gon(5)
    sage: SpineTessellation(s)
    Spine Tessellation of TranslationSurface built from 2 polygons

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


class SpineTessellation(Parent):
    r"""
    A HyperbolicPlaneTessellation + more data
    """

    def __init__(self, surface):
        self._surface_original = surface
        self._surface = self._surface_original.delaunay_triangulation()
        self._hyperbolic_plane = HyperbolicPlane(surface.base_ring())

        shortest_directions = self.shortest_directions(self._surface)
        if len(shortest_directions) == 1:
            raise NotImplementedError("deformation retract onto spine")
        if len(shortest_directions) == 2:
            raise NotImplementedError("move along edge to a vertex")
        assert len(shortest_directions) >= 3

    def _repr_(self):
        return f"Spine Tessellation of {self._surface_original}"

    def explore(self, limit=None, surface=None):
        r"""
        Explore the spine tree up to the combinatorial ``limit`` starting from
        ``translation_surface`` (moving it as in :meth:`__init__`).
        """
        raise NotImplementedError

    def is_vertex(self, translation_surface):
        r"""
        Return whether this is a vertex.
        """
        raise NotImplementedError

    def root(self):
        r"""
        Return the vertex from which we started to build the spine tree.
        """
        return self._hyperbolic_plane.point(0, 1, model="half_plane")

    def point(self, surface):
        r"""
        Return the point in the :class:`HyperbolicPlane` corresponding to ``surface``.

        INPUT:

        - ``surface`` -- A surface in the SL(2, R)-orbit of the defining surface
        """
        raise NotImplementedError

    def edge(self, translation_surface):
        r"""
        Return the spine edge this translation surface is on.

        The edge is oriented such that it points from the :meth:`root`.
        """
        raise NotImplementedError

    def face(self, translation_surface_or_edge):
        r"""
        Return the spine face this translation surface is in.

        Return the face this oriented edge bounds.
        """
        raise NotImplementedError

    def shortest_edges(self, point_or_surface, check=True):
        r"""
        Return the shortest periods for the vertex ``translation_surface``, sorted by slope. There is only one representative of each unoriented edge.

        INPUT:
        - ``check`` -- check whether ``point_or_surface`` is Delaunay triangulated

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.spine_tessellation import SpineTessellation
            sage: s = translation_surfaces.veech_double_n_gon(5)
            sage: t = s.delaunay_triangulation()
            sage: SpineTessellation(s).shortest_edges(t)
            [(0, 2), (5, 2), (2, 1), (1, 1), (2, 2)]


        TODO:: Should computing shortest period be a method for translation surfaces?


        """
        def slope(v):
            from sage.all import oo
            # TODO: Requires ring has division (implement for ExactReal?)
            return v[1]/v[0] if v[0] else oo


        if point_or_surface in self._hyperbolic_plane:
            point_or_surface = self.surface(point_or_surface)

        if check and not point_or_surface.is_delaunay_triangulated():
            raise ValueError("surface must be Delaunay triangulated")

        vectors = {(i, j): point_or_surface.polygon(i).edge(j)
                   for (i, j), _ in point_or_surface.edge_iterator(gluings=True)}
        vectors = {key: value for key, value in vectors.items()
                   if value[0] > 0 or value[0] == 0 and value[1] > 0}

        minimum_length = min(v[0]**2 + v[1]**2 for v in vectors.values())

        shortest_edges = [edge for edge, v in vectors.items()
                            if v[0]**2 + v[1]**2 == minimum_length]

        return sorted(shortest_edges,
                      key=lambda edge: slope(point_or_surface.polygon(edge[0]).edge(edge[1])))

    def shortest_directions(self, point_or_surface, check=True):
        from sage.all import vector

        if point_or_surface in self._hyperbolic_plane:
            point_or_surface = self.surface(point_or_surface)
        shortest_edges = self.shortest_edges(point_or_surface, check=check)
        shortest_directions = {vector(point_or_surface.polygon(t).edge(e), immutable=True): (t, e) for t, e in shortest_edges}
        return list(shortest_directions.keys())

    def standard_form_surface(self, edge):
        r"""
        Return the standard form surface determined by the ``edge``
        as a point on the underlying geodesic.
        """
        raise NotImplementedError

    def surface(self, point):
        r"""
        Return a translation surface corresponding to this ``point`` in the
        hyperbolic plane.
        """
        # This should be shared with the IDR code.
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return self._surface.apply_matrix(matrix([[1, x], [0, y]]))

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

        a, b are saddle connection periods on a starting surface M_0
        {x + iy : |A_z a|^2 = |A_z b|^2, where A_z := [1 x | 0 y]}
        = gamma((b_x - a_x)/(a_y - b_y), -(b_x + a_x)/(a_y + b_y))

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.spine_tessellation import SpineTessellation
            sage: s = translation_surfaces.veech_double_n_gon(5)
            sage: T = SpineTessellation(s)
            sage: T.geodesics(T.root())
            [{-a*(x^2 + y^2) + (2*a^2 - 6)*x + a = 0},
             {(-a^3 + 3*a)*(x^2 + y^2) + (-2*a^2 + 4)*x + a^3 - 3*a = 0},
             {(a^3 - 3*a)*(x^2 + y^2) + (-2*a^2 + 4)*x - a^3 + 3*a = 0},
             {a*(x^2 + y^2) + (2*a^2 - 6)*x - a = 0},
             {4*x = 0}]
        """
        from sage.all import oo, matrix, vector

        # vertex = A * M_0, gamma(t) = sigma(A)^-1 * Rot * (e^{-2t} i)
        # t -> oo: sigma(A)^-1 * Rot * 0
        # t -> -oo: sigma(A)^-1 * Rot * oo
        # sigma([a b | c d])^{-1} = [d b | c a]
        # v, w are the shortest periods on A * M_0
        # [Rot90Clockwise(v + w) | v + w]^T

        shortest_directions = self.shortest_directions(vertex)
        x, y = vertex.coordinates(model="half_plane")
        geodesics = []
        for v, w in zip(shortest_directions, shortest_directions[1:] + [-shortest_directions[0]]):
            def rotation90clockwise(v): return vector([v[1], -v[0]])
            rotation = matrix([rotation90clockwise(v + w), v + w]).transpose()
            sigmaAinv = matrix([[y, x], [0, 1]])
            [[a, b], [c, d]] = sigmaAinv * rotation
            start = a/c if c != 0 else oo
            end = b/d if d != 0 else oo
            geodesic = self._hyperbolic_plane.geodesic(start, end)
            assert vertex in geodesic
            geodesics.append(geodesic)

        return geodesics

    def segment(self, vertex, geodesic):
        r"""
        Return the segment starting at ``vertex`` in direction of the oriented ``geodesic``.

        The segment consists of surfaces where exactly two periods of ``vertex`` are the shortest.
        """
        raise NotImplementedError

    def segments(self, vertex):
        r"""
        Return the segments starting at ``vertex``.

        Each segment consists of surfaces where exactly two periods of ``vertex`` are the shortest.

        """
        r"""
        (g_t*R) * (A * M_0)
        (g_t*R)v0 = (g_t*R)v

        g_t * i = e^(2t) i , t->oo get oo, t->-oo get 0, 
        """
        return [self.segment(vertex, geodesic) for geodesic in self.geodesics(vertex)]
