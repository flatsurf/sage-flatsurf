r"""
EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.spine_tessellation import SpineTessellation
    sage: s = translation_surfaces.mcmullen_L(1, 1, 1, 1)
    sage: SpineTessellation(s)
    Spine Tessellation of TranslationSurface built from 3 polygons

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

        shortest_periods = self.shortest_periods(self._surface)
        if len(shortest_periods) == 1:
            raise NotImplementedError("deformation retract onto spine")
        if len(shortest_periods) == 2:
            raise NotImplementedError("move along edge to a vertex")
        assert len(shortest_periods) >= 3

        self._hyperbolic_plane = HyperbolicPlane(surface.base_ring())

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
        return self._hyperbolic_plane.point(0, 1)

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

    def shortest_periods(self, point_or_surface, check=True):
        r"""
        Return the shortest periods for the vertex ``translation_surface``, sorted by slope.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.spine_tessellation import SpineTessellation
            sage: s = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: t = s.delaunay_triangulation()
            sage: SpineTessellation(s).shortest_periods(t)
            [(3, 1), (4, 1), (5, 1), (3, 2), (4, 2), (5, 2)]

        TODO:: Should computing shortest period be a method for translation surfaces?


        """
        def slope(v):
            from sage.all import oo
            # TODO: Requires ring has division (do we want to implement for ExactReal?)
            return v[1]/v[0] if v[0] else oo

        if check and not point_or_surface.is_delaunay_triangulated():
            raise ValueError("surface must be Delaunay triangulated")

        vectors = {(i, j): point_or_surface.polygon(i).edge(j)
                   for (i, j), _ in point_or_surface.edge_iterator(gluings=True)}
        vectors = {key: value for key, value in vectors.items()
                   if value[0] > 0 or value[0] == 0 and value[1] > 0}

        minimum_length = min(v[0]**2 + v[1]**2 for v in vectors.values())

        shortest_periods = [edge for edge, v in vectors.items()
                            if v[0]**2 + v[1]**2 == minimum_length]

        return sorted(shortest_periods,
                      key=lambda edge: slope(point_or_surface.polygon(edge[0]).edge(edge[1])))

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
        raise NotImplementedError

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
            sage: s = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = SpineTessellation(s)
            sage: T.geodesics(T.root())
        """
        from sage.all import oo

        # TODO: assert that the vertex is a point
        shortest_periods = self.shortest_periods(vertex)
        geodesics = []
        for (ax, ay), (bx, by) in zip(shortest_periods, shortest_periods[1:] + shortest_periods[:1]):
            # TODO: determine correct orientation
            endpoint0 = (bx - ax)/(ay - by) if ay != by else oo
            endpoint1 = -(bx + ax)/(ay + by) if ay != -by else oo
            geodesic = self._hyperbolic_plane.geodesic(endpoint0, endpoint1)
            assert vertex in geodesic
            geodesics.append(geodesic)
        
        return geodesics

    def segments(self, vertex):
        r"""
        Return the segments starting at ``vertex``.

        Each segment consists of surfaces where two periods of ``vertex`` are the shortest.

        """
        raise NotImplementedError
