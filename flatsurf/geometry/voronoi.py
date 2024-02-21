r"""
Voronoi diagrams on cone surfaces

A Voronoi diagram on a cone surface partitions a surface into Voronoi cells
(and their boundary which is shared between adjacent cells.) Each vertex of the
surface is surrounded by a Voronoi cell. Each point of the surface belongs to
the Voronoi cells whose center it is "closest" to.

Clasically, the distance between points is measured as the length of the
shortest path between two points (as induced by the Euclidean distance in the
polygons that make up the surface.) With this notion of distance the
boundaries of the Voronoi cells can be described by line segments.

Points on the boundary of a Voronoi cell are the points for which two different
shortest paths of equal length to a vertex exist. Note that these two paths may
end at the same vertex. For example on the square torus, i.e., the surface
built from a square with opposite sides glued, the boundary of the single
Voronoi cell is given by a + shaped set centered in the square.

We also support a more generalized notion of distance which allows a positive
weight to be assigned to each vertex. The distance of a point to a vertex is
then given by the usual distance divided by that weight. The Voronoi cells in
this setup are delimited by line segments and segments of circular arcs.

EXAMPLES::

    sage: from flatsurf import translation_surfaces, VoronoiTessellation  # random output due to deprecation warnings
    sage: S = translation_surfaces.square_torus()
    sage: V = VoronoiTessellation(S)
    sage: V
    Voronoi Tessellation of Translation Surface in H_1(0) built from a square
    sage: cell = next(iter(V.cells()))
    sage: cell.boundary()
    `+`

"""
# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023-2024 Julian Rüth
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
# ####################################################################


from flatsurf.geometry.surface_objects import SurfacePoint


class VoronoiTessellation:
    r"""
    A Voronoi diagram on a ``surface`` with vertex ``weights``.

    ALGORITHM:

    Internally, a Voronoi diagram is composed of :class:`VoronoiCell`s. Each
    cell is the union of :class:`VoronoiPolygonCell`s which represents the
    intersection of a Voronoi cell with an Euclidean polygon making up the
    surface. To compute these cells on the level of polygons, we determine the
    curves that bound the Voronoi cell.

    When there are no weights, these curves are line segments situated half way
    between the center vertex of the cell and the center vertex of an adjacent
    cell. The segments are orthogonal to the line connecting the two vertices.

    When there are weights, the cells are bounded by line segments of the same
    kind, when two adjacent vertices have equal weights. Otherwise, the
    boundary is given by a segment of a circular arc. Indeed, the distance from
    a point to a vertex is scaled by the inverse of that weight. So the points
    that are at equal distance to two vertices must have a classical distance
    whose ratio is the same as the ratio of the corresponding weights. Hence
    these points are on a circle of Apollonius with the vertices as foci. (Note
    that none of the vertices is at the center of the circle.)

    INPUT:

    - ``surface`` -- an immutable cone surface

    - ``weights`` -- a dict mapping vertices of the ``surface`` to their
      weights or ``None`` (default: ``None``); if ``None``, all vertices have
      the same weight. Weights can be elements of the base ring, floating
      point numbers, or Euclidean distances.

    EXAMPLES::

        sage: from flatsurf import VoronoiTessellation, translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: center = S(0, S.polygon(0).centroid())
        sage: V = VoronoiTessellation(S)
        sage: V.plot()
        sage: V = VoronoiTessellation(S, {vertex: VoronoiTessellation.radius_of_convergence(vertex) for vertex in S.vertices})
        sage: V.plot()

    The same Voronoi diagram but starting from a more complicated description
    of the octagon::

        sage: from flatsurf import Polygon, translation_surfaces, polygons, VoronoiTessellation
        sage: S = translation_surfaces.regular_octagon()
        sage: S = S.subdivide().codomain()
        sage: V = VoronoiTessellation(S)
        sage: V.plot()

        sage: V = VoronoiTessellation(S, {v: VoronoiTessellation.radius_of_convergence(v) for v in S.vertices()})
        sage: V.plot()

    .. SEEALSO::

        :meth:`radius_of_convergence` to create ``weights``

    """

    def __init__(self, surface, weights=None):
        if surface.is_mutable():
            raise ValueError("surface must be immutable")

        if weights is None:
            weights = {vertex: 1 for vertex in surface.vertices()}

        self._surface = surface
        self._weights = weights

        self._cells = {vertex: VoronoiCell(self, vertex) for vertex in surface.vertices()}

    def surface(self):
        r"""
        Return the surface which this diagram is breaking up into cells.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: V.surface() is S
            True

        """
        return self._surface

    def __repr__(self):
        return f"Voronoi Tessellation of {self._surface}"

    def plot(self, graphical_surface=None):
        r"""
        Return a graphical representation of this Voronoi diagram.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiTessellation(S)
            sage: V.plot()

        The underlying surface is not plotted automatically when it is provided
        as a keyword argument::

            sage: V.plot(graphical_surface=S.graphical_surface())

        """
        plot_surface = graphical_surface is None

        if graphical_surface is None:
            graphical_surface = self._surface.graphical_surface(edge_labels=False, polygon_labels=False)

        plot = []
        if plot_surface:
            plot.append(graphical_surface.plot())

        for cell in self._cells.values():
            plot.append(cell.plot(graphical_surface))

        return sum(plot)

    def cells(self, point=None):
        r"""
        Return the Voronoi cells that contain ``point``.

        If no ``point`` is given, return all cells.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: S = S.insert_marked_points(S(0, S.polygon(0).centroid())).codomain()
            sage: V = VoronoiDiagram(S)
            sage: cells = V.cells(S(0, (1/2, 0)))
            sage: list(cells)
            [Voronoi cell at Vertex 0 of polygon 0]

        """
        for cell in self._cells.values():
            if point is None or cell.contains_point(point):
                yield cell

    def polygon_cells(self, label, coordinates):
        r"""
        Return the Voronoi cells that contain the point ``coordinates`` of the
        polygon with ``label``.

        This is similar to :meth:`cells` but it returns the
        :class:`VoronoiPolygonCell`s instead of the full :class:`VoronoiCell`s.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: V.polygon_cells(0, (1/2, 0))
            [Voronoi cell in polygon 0 at (0, 0), Voronoi cell in polygon 0 at (1, 0)]

        """
        raise NotImplementedError


class VoronoiCell:
    def __init__(self, tessellation: VoronoiTessellation, center: SurfacePoint):
        self._tessellation = tessellation
        self._center = center

    def plot(self, graphical_surface=None):
        r"""
        Return a graphical representation of this cell.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiTessellation(S)
            sage: cell = V.cell(S(0, 0))
            sage: cell.plot()

        """
        plot_surface = graphical_surface is None

        if graphical_surface is None:
            graphical_surface = self.surface().graphical_surface()

        plot = []
        if plot_surface:
            plot.append(graphical_surface.plot())

        for polygon_cell in self.polygon_cells():
            plot.append(polygon_cell.plot(graphical_polygon=graphical_surface.graphical_polygon(polygon_cell.label())))

        return sum(plot)

    def surface(self):
        r"""
        Return the surface which this cell is a subset of.

    Voronoi Tessellation of Translation Surface in H_1(0) built from a square
        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiTessellation(S)
            sage: cell = V.cell(S(0, 0))
            sage: cell.surface() is S
            True

        """
        return self._tessellation.surface()

    def boundary(self):
        raise NotImplementedError

    def contains_point(self, point):
        r"""
        Return whether this cell contains the ``point`` of the surface.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.cell(S(0, 0))
            sage: point = S(0, (0, 0))
            sage: cell.contains_point(point)
            True

        """
        label, coordinates = point.representative()
        return any(cell.contains_point(coordinates) for cell in self.polygon_cells(label=label))

    def polygon_cells(self, label=None):
        r"""
        Return this cell broken into pieces that live entirely inside a
        polygon.

        If ``label`` is specified, only the bits that are inside the polygon
        with ``label`` are returned.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.cell(S(0, 0))
            sage: cell.polygon_cells()
            [Voronoi cell in polygon 0 at (0, a + 1),
             Voronoi cell in polygon 0 at (1, a + 1),
             Voronoi cell in polygon 0 at (0, 0),
             Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a + 1),
             Voronoi cell in polygon 0 at (1, 0),
             Voronoi cell in polygon 0 at (-1/2*a, 1/2*a),
             Voronoi cell in polygon 0 at (-1/2*a, 1/2*a + 1),
             Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a)]

        """
        cells = []

        for (lbl, coordinates) in self._center.representatives():
            if label is not None and label != lbl:
                continue

            for c in self._tessellation.polygon_cells(lbl, coordinates):
                cells.append(c)

        assert cells

        return cells


