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

Here, we implement a more generalized notion of distance which assigns to each
vertex a weight. The distance of a point to a vertex is then given by the usual
distance divided by that weight. The Voronoi cells in this setup are delimited
by line segments and segments of circular arcs.

EXAMPLES::

    sage: from flatsurf import translation_surfaces, VoronoiTessellation
    sage: S = translation_surfaces.square_torus()
    sage: V = VoronoiTessellation(S)
    sage: V
    Voronoi Tessellation of Translation Surface in H_1(0) built from a square
    sage: cell = next(iter(V.cells()))
    sage: cell.boundary()
    `+`

Using a weight at the sole vertex does not change the picture here::

    sage: V = VoronoiTessellation(S, weights="radius_of_convergence")
    sage: cell = next(iter(V.cells()))
    sage: cell.boundary()

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
from sage.misc.cachefunc import cached_method


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

    - ``weights`` -- ``"radius_of_convergence"`` or ``None`` (default:
      ``None``); the weights to use at each vertex if any.

    EXAMPLES::

        sage: from flatsurf import VoronoiTessellation, translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: center = S(0, S.polygon(0).centroid())
        sage: S = S.insert_marked_points(center).codomain()
        sage: V = VoronoiTessellation(S)
        sage: V.plot()
        sage: V = VoronoiTessellation(S, weights="radius_of_convergence")
        sage: V.plot()

    The same Voronoi diagram but starting from a more complicated description
    of the octagon::

        sage: from flatsurf import Polygon, translation_surfaces, polygons, VoronoiTessellation
        sage: S = translation_surfaces.regular_octagon()
        sage: S = S.subdivide().codomain()
        sage: V = VoronoiTessellation(S)
        sage: V.plot()

        sage: V = VoronoiTessellation(S, weights="radius_of_convergence")
        sage: V.plot()

    """

    def __init__(self, surface, weights=None):
        if surface.is_mutable():
            raise ValueError("surface must be immutable")
        # TODO: Surface must also be a Cone Surface and connected

        if weights is None:
            pass
        elif weights == "radius_of_convergence":
            pass
        else:
            raise ValueError("unsupported weights for Voronoi tessellation")


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

    def cell(self, point):
        r"""
        Return a Voronoi cell that contains ``point``.
m
        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: cell = V.cell(center)
            sage: cell.contains_point(center)
            True

        """
        return next(iter(self.cells(point)))

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

    @cached_method
    def radius_of_convergence(self, center):
        r"""
        Return the radius of convergence when developing a power series
        at ``center``.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S)

            sage: V.radius_of_convergence2(S(0, S.polygon(0).centroid()))
            1/2*a + 1
            sage: V.radius_of_convergence2(next(iter(S.vertices())))
            1
            sage: V.radius_of_convergence2(S(0, (1/2, 1/2)))
            1/2
            sage: V.radius_of_convergence2(S(0, (1/2, 0)))
            1/4
            sage: V.radius_of_convergence2(S(0, (1/4, 0)))
            1/16

        """
        norm = self.surface().euclidean_plane().norm()

        if all(vertex.angle() == 1 for vertex in self._surface.vertices()):
            return norm.infinite()

        erase_marked_points = self._surface.erase_marked_points()
        center = erase_marked_points(center)

        if not center.is_vertex():
            insert_marked_points = center.surface().insert_marked_points(center)
            center = insert_marked_points(center)

        surface = center.surface()

        for connection in surface.saddle_connections():
            start = surface(*connection.start())
            end = surface(*connection.end())
            if start == center and end.angle() != 1:
                x, y = connection.holonomy()
                return norm.from_norm_squared(x**2 + y**2)

        assert False, "unreachable on a surface without boundary"


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
        r"""
        Return the line segments and circular arc segments bounding this
        Voronoi cell.

        The segments are returned in order of a counterclockwise walk (as seen
        from the center of the cell) along the boundary.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiTessellation(S)
            sage: cell = V.cell(S(0, 0))
            sage: cell.boundary()

        """
        return [self.surface()(polygon_cell.label(), boundary) for polygon_cell in self.polygon_cells() for boundary in polygon_cell.boundary()]

    @cached_method
    def _boundary_saddle_connections(self):
        r"""
        Return saddle connections to the vertices that have an influence on the
        shape of this cell.

        ALGORITHM:

        For good choices of weights, such as equal weights or weights
        corresponding to radii of convergence, only vertices that are connected
        to the center by a saddle connection (possibly crossing over a marked
        point) influence the Voronoi cell.

        Again, for good choices of weights, we only have to consider saddle
        connections up to a certain length if we assume that all cells are
        contained in the their disk of convergence. In that case we can stop at
        twice the maximum radius of convergence.

        TODO: Put more of a proof for all this here. See my handwritten notes.

        TODO: Currently, we do not include vertices hidden by a marked point here.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiTessellation(S)
            sage: cell = V.cell(S(0, 0))
            sage: cell._boundary_saddle_connections()

        """
        R = max(self._tessellation.radius_of_convergence(v) for v in self.surface().vertices())

        if not R.is_finite():
            # When there are no singularities, we don't know a priori where to
            # stop the saddle connection search here. So we cheat a bit and ask
            # the Delaunay triangulation (which is dual to the Voronoi diagram
            # in this case) for the maximum edge length which bounds the
            # diameter of any cell.
            T = self.surface().delaunay_triangulation()

            norm = R.parent()
            R = max(norm(edge) for label in T.labels() for edge in T.polygon(label).edges())

        boundary_saddle_connections = list(self.surface().saddle_connections(length_bound=2*R, initial_vertex=self._center))
        assert all(self.surface()(*c.start()) == self._center for c in boundary_saddle_connections)

        return boundary_saddle_connections

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

        Otherwise, all parts of the cell are returned in a counterclockwise
        order around the center of the cell.

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

        initial_label, initial_vertex = self._center.representative()
        initial_vertex = self.surface().polygon(label).get_point_position(initial_vertex).get_vertex()
        label, vertex = initial_label, initial_vertex

        # Walk counterclockwise around the center vertex.
        while True:
            yield self.polygon_cell(label, vertex)

            polygon = self.surface().polygon(label)

            previous_edge = (vertex - 1) % len(polygon.vertices())
            label, vertex = self.surface().opposite_edge(label, previous_edge)

            if label == initial_label and vertex == initial_vertex:
                break

    @cached_method
    def polygon_cell(self, label, vertex):
        r"""
        Return the :class:`VoronoiPolygonCell` centered at the ``vertex`` of
        the polygon with ``label``.

        ALGORITHM:

        # TODO

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, VoronoiTessellation
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: S = S.insert_marked_points(center).codomain()
            sage: V = VoronoiTessellation(S)
            sage: V.cell(S((0, 0), 0)).polygon_cell(0, 0)

        """
        polygon = self.surface().polygon(label)
        start_edge = polygon.edge(vertex)
        end_edge = polygon.edge(vertex + 1)

        precone = self.surface().cone(label, polygon.vertex(vertex), (start_edge[1], -start_edge[0]), start_edge)
        cone = self.surface().cone(label, polygon.vertex(vertex), start_edge, end_edge)
        postcone = self.surface().cone(label, polygon.vertex(vertex), end_edge, (-end_edge[1], end_edge[0]))

        preconnections = [c for c in self._boundary_saddle_connections() if precone.contains_tangent_vector(c.start_tangent_vector())]
        connections = [c for c in self._boundary_saddle_connections() if cone.contains_tangent_vector(c.start_tangent_vector())]
        postconnections = [c for c in self._boundary_saddle_connections() if postcone.contains_tangent_vector(c.start_tangent_vector())]

        preconnections = sorted(preconnections, key=lambda connection: precone.tangent_vector_sort_key(connection.start_tangent_vector()))
        connections = sorted(connections, key=lambda connection: cone.tangent_vector_sort_key(connection.start_tangent_vector()))
        postconnections = sorted(postconnections, key=lambda connection: postcone.tangent_vector_sort_key(connection.start_tangent_vector()))

        # TODO: Map holonomies into polygon if not a translation surface.
        return VoronoiPolygonCell(self, label, [self._polygon_cell_boundary(connection) for connection in preconnections + connections + postconnections if self._is_saddle_connection_defining_cell(connection)])

    def _polygon_cell_boundary(self, connection):
        # TODO: If the visibility from the connection.end to the center is
        # limited by singularities, we would have to crop the half space
        # accordingly. (In practice it probably makes no difference.)
        raise NotImplementedError

    def _is_saddle_connection_defining_cell(self, connection):
        source = self.surface()(*connection.start())
        target = self.surface()(*connection.end())

        R = lambda vertex: self._tessellation.radius_of_convergence(vertex)

        # TODO: Add a proof that this is the case iff the half way point is in
        # both cells (under which conditions?)
        midpoint = connection.bisect(R(self._center), R(self.surface()(*connection.end())))

        distances = midpoint.distances()
        distances = {vertex: distance / R(vertex) for (vertex, distance) in distances.items()}

        min_distance = min(distances.values())

        if distances[source] != min_distance:
            return False

        if distances[target] != min_distance:
            return False

        return min_distance * (R(source) + R(target)) == connection.length()


class VoronoiPolygonCell:
    def __init__(self, cell, boundary):
        pass
