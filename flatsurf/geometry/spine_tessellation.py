#*********************************************************************
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
#*********************************************************************
class SpineTessellation:
    r"""
    A HyperbolicPlaneTessellation + more data
    """
    def __init__(self, translation_surface):
        r"""
        
        """
        # Move to a vertex canonically.
        # We'll explore the spine tree from there.
        
    def explore(self, limit=None, translation_surface=None):
        r"""
        Explore the spine tree up to the combinatorial ``limit`` starting from
        ``translation_surface`` (moving it as in :meth:`__init__`).
        """
        
    def is_vertex(self, translation_surface):
        r"""
        Return whether this is a vertex.
        """
        
    def root(self):
        r"""
        Return the vertex from which we started to build the spine tree.
        """

    def edge(self, translation_surface):
        r"""
        Return the spine edge this translation surface is on.
        
        The edge is oriented such that it points from the :meth:`root`.
        """
        
    def face(self, translation_surface_or_edge):
        r"""
        Return the spine face this translation surface is in.
        
        Return the face this oriented edge bounds.
        """

    def shortest_periods(self, translation_surface):
        r"""
        Return the shortest periods for the vertex ``translation_surface``.
        """
        
    def standard_form_surface(self, edge):
        r"""
        Return the standard form surface determined by the ``edge``
        as a point on the underlying geodesic.
        """
        
    def surface(self, point):
        r"""
        Return a translation surface corresponding to this ``point``.
        """
        # This should be shared with the IDR code.
        
    def fundamental_domain(self):
        r"""
        Return the fundamental domain as a polygon with edge pairings
        
        The result would be a fundamental domain of the group generated by the
        symmetries discovered so far.
        """
        
    def plot(self):
        r"""
        Return a plot.
        """

    def polygon(self, vertex_or_edge):
        r"""
        Return the polygon obtained as the union of the
        triangles bounded by this edge and its reverse /
        by the edges adjacent to this vertex.
        """
        
    def other_endpoint(self, vertex, edge_initial):
        r"""
        Return the new vertex of the enriched spine that is connected to ``vertex``
        in the direction of ``edge_initial``.
        """
        # Compute half-planes for when hinges in ``vertex`` will flip, and for when two vectors in ``vertex`` become same length
        # Compute the geodesic that begins flowing along the edge
        # If a half-plane flips first, change the triangulation (and compute new halfplanes?)
        # If an edge ties first, change the set of systoles