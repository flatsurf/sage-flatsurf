# TODO: Change the inheritance structure so it's clear which implementations
# use the fact that the cells are convex, which use the fact that the cells are
# given additionally by line segments (so given by their corners.)


from sage.misc.cachefunc import cached_method


class CellDecomposition:
    Cell = None

    def __init__(self, surface):
        self._surface = surface

    def surface(self):
        return self._surface

    @cached_method
    def centers(self):
        r"""
        Return the center points of this cell decomposition.
        """
        return frozenset(self._surface.vertices())

    def cell(self, point):
        r"""
        Return a cell containing ``point``.
        """
        return next(iter(self.cells(point)))

    @cached_method
    def cell_at_center(self, center):
        r"""
        Return the cell centered at ``center``.
        """
        if self.Cell is None:
            raise NotImplementedError("this decomposition does not implement cells yet")

        return self.Cell(self, center)

    def cells(self, point=None):
        r"""
        Return the cells containing ``point``.
        """
        if point is None:
            yield from self
            return


        for cell in self:
            if cell.contains_point(point):
                yield cell

    @cached_method
    def polygon_cells(self, label):
        return [polygon_cell for cell in self for polygon_cell in cell.polygon_cells() if polygon_cell.label() == label]

    def split_segment_at_cells(self, segment):
        r"""
        Return the ``segment`` split into shorter segments that each lie in a
        single cell.

        Returns a sequence of pairs, consisting of subsegment and a cell that
        contains the segment.
        """
        raise NotImplementedError("this decomposition does not implement splitting of segments yet")

    def split_segment_at_polygon_cells(self, segment):
        r"""
        Return the ``segment`` split into shorter segments that each lie in a
        single polygon cell.

        Returns a sequence of pairs, consisting of subsegment and a cell that
        contains the segment.
        """
        from flatsurf.geometry.euclidean import OrientedSegment
        return [(OrientedSegment(*s), polygon_cell) for (label, subsegment, _) in segment.split() for (s, polygon_cell) in self._split_segment_at_polygon_cells(label, subsegment)]

    def _split_segment_at_polygon_cells(self, label, segment):
        start, end = segment

        if start == end:
            return []

        for start_cell in self.polygon_cells(label=label):
            if not start_cell.contains_point(start):
                continue

            polygon = start_cell.polygon()

            try:
                exit = polygon.flow_to_exit(start, end - start)
            except ValueError:
                continue

            return [((start, exit), start_cell)] + self._split_segment_at_polygon_cells(label, (exit, end))

        assert False

    def __iter__(self):
        r"""
        Return an iterator over the cells of this decomposition.
        """
        for center in self.centers():
            yield self.cell_at_center(center)

    def _neg_boundary_segment(self, boundary_segment):
        r"""
        Return the boundary segment that has the same segment as
        ``boundary_segment`` but with the opposite orientation.
        """
        search = -boundary_segment.segment()

        for cell in self:
            for boundary in cell.boundary():
                if boundary.segment() == search:
                    return boundary

        assert False, "boundary segment has no negative in this decomposition"

    def __repr__(self):
        raise NotImplementedError("this cell decomposition cannot print itself yet")

    def plot(self, graphical_surface=None):
        r"""
        Return a graphical representation of this cell decomposition.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: S = S.insert_marked_points(center).codomain()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: V.plot()
            Graphics object consisting of 73 graphics primitives

        The underlying surface is not plotted automatically when it is provided
        as a keyword argument::

            sage: V.plot(graphical_surface=S.graphical_surface())
            Graphics object consisting of 48 graphics primitives

        """
        plot_surface = graphical_surface is None

        if graphical_surface is None:
            graphical_surface = self._surface.graphical_surface(edge_labels=False, polygon_labels=False)

        plot = []
        if plot_surface:
            plot.append(graphical_surface.plot())

        for cell in self:
            plot.append(cell.plot(graphical_surface))

        return sum(plot)

    def _test_boundary_consistency(self):
        # TODO: Use SageMath test framework.
        segments = [boundary.segment() for cell in self for boundary in cell.boundary()]

        for cell in self:
            for boundary in cell.boundary():
                assert -boundary.segment() in segments, f"{boundary} has no negative {-boundary.segment()}"


class Cell:
    BoundarySegment = None
    PolygonCell = None

    def __init__(self, decomposition, center):
        self._decomposition = decomposition
        self._center = center

    def surface(self):
        return self._decomposition.surface()

    def decomposition(self):
        return self._decomposition

    def center(self):
        r"""
        Return the point of the surface that should be considered as the center
        of this cell.
        """
        return self._center

    def polygon_cells(self, label=None):
        r"""
        Return restrictions of this cell to all polygons of the surface.

        If ``label`` is given, return only the restrictions to that polygon.
        """
        if label is None:
            return sum([self.polygon_cells(label=label) for label in self.surface().labels()], start=())

        return self._polygon_cells(label=label)

    def _polygon_cells(self, label):
        r"""
        Return the restrictions of this cell to its connected components in the
        polygon with ``label``.

        Note that a cell can be separated by a boundary from itself so a single
        cell can induce several polygon cells in a polygon even if these
        polygon cells touch at their boundaries.
        """
        raise NotImplementedError("this cell decomposition cannot restrict cells to polygons yet")

    @cached_method
    def boundary(self):
        r"""
        Return the boundary of this cell in a counterclockwise walk.

        This method is only available if the boundary is connected, i.e. the
        cell is simply connected.
        """
        raise NotImplementedError("this cell decomposition cannot compute boundaries of cells yet")

    def radius(self):
        r"""
        Return the distance from the center to the furthest point in this cell.
        """
        return max(boundary.radius() for boundary in self.boundary())

    def point_from_complex(self, z):
        from math import pi

        arg = z.arg()
        if arg < 0:
            arg += 2*pi

        branch = int(arg * self._center.angle() // (2*pi))

        from sage.all import vector
        xy = vector(list(z ** self._center.angle()))

        # TODO: This is all not too robust since we are feeding
        # floating point numbers into exact polygons here. As a
        # result roots could show up twice or never.
        for polygon_cell in self.polygon_cells():
            polygon = self.surface().polygon(polygon_cell.label())
            center = polygon_cell.center()
            vertex = polygon.get_point_position(center).get_vertex()

            from flatsurf.geometry.euclidean import ccw
            if ccw(polygon.edge(vertex), xy) >= 0 and ccw(-polygon.edge(vertex - 1), xy) < 0:
                from flatsurf.geometry.euclidean import OrientedSegment
                p = center + xy
                p = p.change_ring(self.surface().base_ring())
                if polygon_cell.root_branch(OrientedSegment(p, p + xy)) == branch:
                    if polygon.get_point_position(p).is_outside():
                        return None

                    if not polygon_cell.contains_point(p):
                        return None

                    print(f"Root at distance {xy.norm():.5} of {self.center()} (angle {2*self.center().angle()}π)")

                    return self.surface()(polygon_cell.label(), p)

        assert False

    def __eq__(self, other):
        if not isinstance(other, Cell):
            return False
        return self._decomposition == other._decomposition and self._center == other._center

    def __hash__(self):
        return hash(self._center)

    def __repr__(self):
        return f"Cell at {self.center()}"

    def plot(self, graphical_surface=None):
        r"""
        Return a graphical representation of this cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.cell(S(0, 0))
            sage: cell.plot()
            Graphics object consisting of 34 graphics primitives

        """
        plot_surface = graphical_surface is None

        if graphical_surface is None:
            graphical_surface = self.surface().graphical_surface(edge_labels=False, polygon_labels=False)

        plot = []
        if plot_surface:
            plot.append(graphical_surface.plot())

        for boundary in self.boundary():
            plot.append(boundary.plot(graphical_surface=graphical_surface))

        return sum(plot)

    def contains_point(self, point):
        for label, coordinates in point.representatives():
            for polygon_cell in self.polygon_cells(label=label):
                if polygon_cell.contains_point(coordinates):
                    return True

        return False


class LineSegmentCell(Cell):
    r"""
    A cell whose boundary is given by line segments.
    """

    @cached_method
    def corners(self):
        return [segment.segment().start() for segment in self.boundary()]

    @cached_method
    def _polygon_cells(self, label):
        surface = self.surface()
        polygon = surface.polygon(label)

        cell_boundary = self.boundary()

        if any(cell_boundary[b].segment().start() != cell_boundary[b - 1].segment().end() for b in range(len(cell_boundary))):
            raise NotImplementedError("boundary of cell must be connected")

        from collections import namedtuple
        PolygonCellBoundary = namedtuple("PolygonCellBoundary", ("cell_boundary", "label", "segment", "center"))

        # Filter out the bits of the boundary that lie in the polygon with label.
        polygon_cells_boundary = []

        for boundary in cell_boundary:
            for lbl, subsegment, start_segment in boundary.segment().split(label=label):
                assert lbl == label
                center = subsegment[0] - start_segment - boundary._center_to_start.holonomy()
                center.set_immutable()
                polygon_cells_boundary.append(PolygonCellBoundary(boundary, label, subsegment, center))

        polygon_cells = []

        if not polygon_cells_boundary:
            # TODO: Here we use the assumption that entire polygons cannot be contained in a cell.
            return ()

        unused_polygon_cells_boundaries = set(polygon_cells_boundary)

        from collections import defaultdict
        polygon_cells_boundary_from = defaultdict(lambda: [])
        for boundary in polygon_cells_boundary:
            start, end = boundary.segment
            polygon_cells_boundary_from[start].append(boundary)

        while unused_polygon_cells_boundaries:
            # Build a new polygon cell starting from a random boundary segment.
            polygon_cell = [unused_polygon_cells_boundaries.pop()]
            center = polygon_cell[0].center

            # Walk the boundary of the polygon cell until it closed up.
            while True:
                # Find the first segment at the end point of the previous
                # segment that is at a clockwise turn but still within the polygon.
                start, end = polygon_cell[-1].segment

                strictly_clockwise_from = start - end

                end_position = polygon.get_point_position(end)
                if end_position.is_vertex():
                    counterclockwise_from = polygon.edge(end_position.get_vertex())
                elif end_position.is_in_edge_interior():
                    counterclockwise_from = polygon.edge(end_position.get_edge())
                else:
                    counterclockwise_from = start - end

                # Find candidate next boundaries that are starting at the end
                # point of the previous boundary and in the correct sector.
                from flatsurf.geometry.euclidean import is_between, is_parallel
                polygon_cells_boundary_from_in_sector = [boundary for boundary in polygon_cells_boundary_from[end] if 
                    (not is_parallel(counterclockwise_from, strictly_clockwise_from) and is_parallel(boundary.segment[1] - boundary.segment[0], counterclockwise_from)) or 
                    is_between(counterclockwise_from, strictly_clockwise_from, boundary.segment[1] - boundary.segment[0])]

                # Pick the first such boundary turning clockwise from the
                # previous boundary.
                def angle_key(boundary):
                    class AngleKey:
                        def __init__(self, vector, sector):
                            self._vector = vector
                            self._sector = sector
                            assert is_parallel(self._sector[0], self._vector) or is_between(*self._sector, self._vector)

                        def __gt__(self, other):
                            from flatsurf.geometry.euclidean import is_parallel
                            assert not is_parallel(self._vector, other._vector), "cell must not have equally oriented parallel boundaries"

                            if is_parallel(self._sector[0], self._vector):
                                return False
                            if is_parallel(self._sector[0], other._vector):
                                return True
                            return is_between(self._sector[0], self._vector, other._vector)

                    return AngleKey(boundary.segment[1] - boundary.segment[0], (counterclockwise_from, strictly_clockwise_from))

                next_polygon_cell_boundary = max(polygon_cells_boundary_from_in_sector, key=angle_key, default=None)

                if next_polygon_cell_boundary is None:
                    assert len(polygon_cells_boundary_from_in_sector) == 0

                    # This boundary touches the polygon edges but no explicit
                    # boundary starts at the point where it touches.
                    if end_position.is_vertex():
                        edge = end_position.get_vertex()

                    else:
                        # Part of a polygon edge needs to be added to the boundary
                        # of this cell.
                        # Walk the polygon edges until we find a starting point of
                        # a boundary segment (or hit a vertex.)
                        assert end_position.is_in_edge_interior(), "boundary segments ends in interior of polygon but no other segment continues here"
                        edge = end_position.get_edge()

                    from flatsurf.geometry.euclidean import time_on_ray
                    polygon_cells_boundary_on_edge = [boundary for boundary in polygon_cells_boundary if
                        polygon.get_point_position(boundary.segment[0]).is_in_edge_interior() and
                        polygon.get_point_position(boundary.segment[0]).get_edge() == edge and
                        time_on_ray(end, polygon.edge(edge), boundary.segment[0])[0] > 0
                    ]

                    next_end = min(((time_on_ray(end, polygon.edge(edge), boundary.segment[0])[0], boundary.segment[0]) for boundary in polygon_cells_boundary_on_edge), default=(None, None))[1]
                    if next_end is None:
                        # No segment starts on this edge. Go to the vertex.
                        next_end = polygon.vertex(edge + 1)

                    next_polygon_cell_boundary = PolygonCellBoundary(None, label, (end, next_end), center)
                else:
                    if next_polygon_cell_boundary == polygon_cell[0]:
                        break
                    assert next_polygon_cell_boundary in unused_polygon_cells_boundaries, f"boundary segment present in multiple polygon cells"
                    unused_polygon_cells_boundaries.remove(next_polygon_cell_boundary)

                assert next_polygon_cell_boundary not in polygon_cell, "boundary segment must not repeat in polygon cell boundary"

                assert next_polygon_cell_boundary.center == center, f"segments in cell boundary refer to inconsistent cell centers, ({next_polygon_cell_boundary.center} != {center})"

                polygon_cell.append(next_polygon_cell_boundary)
                
            # Build an actual polygon cell from the boundary segments
            assert center is not None
            polygon_cells.append(self.PolygonCell(cell=self, label=label, center=center, boundary=tuple([boundary.segment for boundary in polygon_cell])))

        return tuple(polygon_cells)


class BoundarySegment:
    def __init__(self, cell, segment):
        self._cell = cell
        self._segment = segment

    def surface(self):
        return self._cell.surface()

    def cell(self):
        return self._cell

    def __bool__(self):
        return bool(self._holonomy())

    def __neg__(self):
        return self.cell().decomposition()._neg_boundary_segment(self)

    def segment(self):
        return self._segment

    def polygon_cell_boundaries(self):
        r"""
        Return this segment split into subsegments that each live entirely
        within a polygon.

        Returns the subsegments as triples (polygon cell, subsegment in polygon
        cell, opposite polygon cell).
        """
        raise NotImplementedError("this cell decomposition cannot restrict its boundaries to polygons yet")

    def __eq__(self, other):
        if not isinstance(other, BoundarySegment):
            return False

        return self.cell() == other.cell() and self.segment() == other.segment()

    def __hash__(self):
        return hash(self.segment())

    def __repr__(self):
        return f"Boundary at {self.segment()}"


class LineBoundarySegment(BoundarySegment):
    r"""
    A segment on the boundary of a cell that is a line segment.

    This segment and the segments connecting the end points of the segment to
    the center of the cell form a triangle in the plane.

    INPUT:

    - ``center_to_start`` -- a segment from the center of the Voronoi cell to
      the start of ``segment``

    """
    def __init__(self, cell, segment, center_to_start):
        super().__init__(cell, segment)

        
        self._center_to_start = center_to_start

    @cached_method
    def polygon_cell_boundaries(self):
        boundaries = []
        for (label, subsegment, start_holonomy) in self.segment().split():
            polygon = self.surface().polygon(label)
            for polygon_cell in self._cell.polygon_cells(label):
                if subsegment in polygon_cell.boundary():
                    midpoint = (subsegment[0] + subsegment[1]) / 2
                    midpoint_position = polygon.get_point_position(midpoint)
                    assert midpoint_position.is_inside()
                    assert not midpoint_position.is_vertex()

                    opposite_label = None
                    if not midpoint_position.is_in_edge_interior():
                        opposite_label = label
                    else:
                        edge = midpoint_position.get_edge()
                        opposite_label, opposite_edge = self.surface().opposite_edge(label, edge)
                        opposite_polygon = self.surface().polygon(opposite_label)
                        midpoint += (opposite_polygon.vertex(opposite_edge + 1) - polygon.vertex(edge))

                    opposite_polygon_cell = [opposite_polygon_cell for opposite_polygon_cell in (-self)._cell.polygon_cells(opposite_label) if opposite_polygon_cell.contains_point(midpoint) and opposite_polygon_cell != polygon_cell]

                    assert opposite_polygon_cell, "no polygon cell on the other side of this boundary"
                    assert len(opposite_polygon_cell) == 1, "more than one polygon cell on the other side of this boundary"
                    opposite_polygon_cell = opposite_polygon_cell[0]

                    boundaries.append((polygon_cell, subsegment, opposite_polygon_cell))
                    break
            else:
                assert False, "subsegment of boundary must also be a boundary on the level of polygons"
    
        return boundaries

    def plot(self, graphical_surface=None):
        return self.segment().plot(graphical_surface=graphical_surface)

    def radius(self):
        norm = self.surface().euclidean_plane().norm()
        return max([
            norm.from_vector(self._center_to_start.holonomy()),
            norm.from_vector(self._center_to_start.holonomy() + self._segment.holonomy()),
        ])


class PolygonCell:
    r"""
    A :class:`Cell` restricted to a polygon of the surface.

    INPUT:

    - ``cell`` -- the :class:`Cell` which this polygon cell is a restriction of

    - ``label`` -- the polygon to which this was restricted

    """
    def __init__(self, cell, label):
        self._cell = cell
        self._label = label

    def label(self):
        return self._label

    def polygon(self):
        return self.cell().surface().polygon(self._label)

    def surface(self):
        return self.cell().surface()

    def cell(self):
        return self._cell

    def contains_point(self, point):
        raise NotImplementedError

    def contains_segment(self, segment):
        raise NotImplementedError

    def boundary(self):
        raise NotImplementedError("this cell decomposition does not know how to compute the boundary segments of a polygon cell yet")

    def __eq__(self, other):
        raise NotImplementedError

    def __hash__(self):
        raise NotImplementedError

    def corners(self):
        raise NotImplementedError


class PolygonCellWithCenter(PolygonCell):
    r"""
    A :class:`Cell` restricted to a polygon of the surface.

    INPUT:

    - ``cell`` -- the :class:`Cell` which this polygon cell is a restriction of

    - ``label`` -- the polygon to which this was restricted

    - ``center`` -- the position of the center of this cell; if outside of the
      polygon, then a segment from any interior point of the polygon to this point is
      assumed to correspond to an actual line segment in the surface connecting
      that point to the center.

    """
    def __init__(self, cell, label, center):
        super().__init__(cell, label)

        if center is None:
            raise ValueError

        self._center = center

    def center(self):
        r"""
        Return the point (not necessarily of the polygon) that should be
        considered the center of this cell.
        """
        return self._center

    def split_segment_with_constant_root_branches(self, segment):
        r"""
        Return the ``segment`` split into smaller segments such that these
        segments do not cross the horizontal line left of the center of the
        cell (if that center is an actual singularity and not just a marked
        point.)

        On such a shorter segment, we can then develop an n-th root
        consistently where n-1 is the order of the singularity.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (1, 0))

            sage: from flatsurf.geometry.euclidean import OrientedSegment
            sage: cell.split_segment_uniform_root_branch(OrientedSegment((0, -1), (0, 1)))
            [OrientedSegment((0, -1), (0, 0)), OrientedSegment((0, 0), (0, 1))]

        """
        from flatsurf.geometry.euclidean import OrientedSegment

        d = self.cell().center().angle()

        assert d >= 1
        if d == 1:
            return [segment]

        if segment.contains_point(self.center()) and segment.start() != self.center() and segment.end() != self.center():
            return [OrientedSegment(segment.start(), self.center()), OrientedSegment(self.center(), segment.end())]

        from flatsurf.geometry.euclidean import Ray
        ray = Ray(self.center(), (-1, 0))

        if ray.contains_point(segment.start()) or ray.contains_point(segment.end()):
            return [segment]

        intersection = ray.intersection(segment)

        if intersection is None:
            return [segment]

        return [OrientedSegment(segment.start(), intersection), OrientedSegment(intersection, segment.end())]

    def root_branch(self, segment):
        r"""
        Return which branch can be taken consistently along the ``segment``
        when developing an n-th root at the center of this Voronoi cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))

            sage: from flatsurf.geometry.euclidean import OrientedSegment
            sage: cell.root_branch(OrientedSegment((0, -1), (0, 1)))
            Traceback (most recent call last):
            ...
            ValueError: segment does not permit a consistent choice of root

            sage: cell.root_branch(OrientedSegment((0, 0), (0, 1/2)))
            0

            sage: cell = V.polygon_cell(0, (1, 0))
            sage: cell.root_branch(OrientedSegment((0, 0), (0, 1/2)))
            1

            sage: a = S.base_ring().gen()
            sage: cell = V.polygon_cell(0, (1 + a/2, a/2))
            sage: cell.root_branch(OrientedSegment((1, 1/2 + a/2), (1 + a/2, 1/2 + a/2)))
            2

        """
        if self.split_segment_with_constant_root_branches(segment) != [segment]:
            raise ValueError("segment does not permit a consistent choice of root")

        S = self.surface()

        center = self.cell().center()

        angle = center.angle()

        assert angle >= 1
        if angle == 1:
            return 0

        # Choose a horizontal ray to the right, that defines where the
        # principal root is being used. We use the "smallest" vertex in the
        # "smallest" polygon containing such a ray.
        from flatsurf.geometry.euclidean import ccw
        primitive_label, primitive_vertex = min((label, vertex) for (label, _) in center.representatives() for vertex in range(len(S.polygon(label).vertices()))
            if S(label, vertex) == center and
               ccw((1, 0), S.polygon(label).edge(vertex)) <= 0 and
               ccw((1, 0), -S.polygon(label).edge(vertex - 1)) >= 0)

        # Walk around the vertex to determine the branch of the root for the
        # (midpoint of) the segment.
        point = segment.midpoint()

        branch = 0
        label = primitive_label
        vertex = primitive_vertex

        while True:
            polygon = S.polygon(label)
            if label == self.label() and polygon.vertex(vertex) == self.center():
                low = ccw((-1, 0), polygon.edge(vertex)) <= 0 and ccw((-1, 0), point - polygon.vertex(vertex)) > 0
                if low:
                    return (branch + 1) % angle
                return branch

            if ccw((-1, 0), polygon.edge(vertex)) <= 0 and ccw((-1, 0), -polygon.edge(vertex - 1)) > 0:
                branch += 1
                branch %= angle

            label, vertex = S.opposite_edge(label, (vertex - 1) % len(polygon.vertices()))

    def __repr__(self):
        return f"{self.cell()} restricted to polygon {self.label()} at {self._center}"


class LineSegmentPolygonCell(PolygonCellWithCenter):
    r"""
    A :class:`Cell` restricted to a polygon of the surface.

    INPUT:

    - ``cell`` -- the :class:`Cell` which this polygon cell is a restriction of

    - ``label`` -- the polygon to which this was restricted

    - ``center`` -- the position of the center of this cell; if outside of the
      polygon, then a segment from any interior point of the polygon to this point is
      assumed to correspond to an actual line segment in the surface connecting
      that point to the center.

    - ``boundary`` -- segments that are restrictions of the boundary of the
      cell to this polygon, in counterclockwise order around the ``center``.
    """
    def __init__(self, cell, label, center, boundary):
        super().__init__(cell, label, center)

        self._boundary = boundary

    def boundary(self):
        return self._boundary

    def contains_point(self, point):
        return self.polygon().contains_point(point)

    def polygon(self):
        r"""
        Return a polygon (subset of :meth:`polygon`) that describes this cell.
        """
        from flatsurf import Polygon
        return Polygon(vertices=[segment[0] for segment in self._boundary])

    def __eq__(self, other):
        if not isinstance(other, PolygonCell):
            return False

        return self.cell() == other.cell() and self.label() == other.label() and self.center() == other.center()

    def __hash__(self):
        return hash((self.label(), self.center()))

    def plot(self, graphical_polygon=None):
        r"""
        Return a graphical representation of this cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.plot()
            Graphics object consisting of 3 graphics primitives

        """
        plot_polygon = graphical_polygon is None

        if graphical_polygon is None:
            graphical_polygon = self.surface().graphical_surface().graphical_polygon(self.label())

        shift = graphical_polygon.transformed_vertex(0) - self.polygon().vertex(0)

        from flatsurf.geometry.euclidean import OrientedSegment
        plot = sum(OrientedSegment(*segment).translate(shift).plot(point=False) for segment in self.boundary())

        if plot_polygon:
            pass # plot = graphical_polygon.plot_polygon() + plot

        return plot

    def corners(self):
        return tuple(segment[0] for segment in self._boundary)


class ConvexLineSegmentPolygonCell(LineSegmentPolygonCell):
    def contains_segment(self, segment):
        return self.contains_point(segment.start()) and self.contains_point(segment.end())


class VoronoiCellBoundarySegment(LineBoundarySegment):
    pass

 
class VoronoiPolygonCell(ConvexLineSegmentPolygonCell):
    pass

class VoronoiCell(LineSegmentCell):
    BoundarySegment = VoronoiCellBoundarySegment
    PolygonCell = VoronoiPolygonCell

    @cached_method
    def boundary(self):
        surface = self.surface()

        (label, coordinates) = self.center().representative()
        vertex = surface.polygon(label).get_point_position(coordinates).get_vertex()

        boundary = []

        initial_label, initial_vertex = label, vertex
        while True:
            polygon = surface.polygon(label)

            center = polygon.circumscribing_circle().center().vector()

            next_label, next_vertex = surface.opposite_edge(label, (vertex - 1) % len(polygon.vertices()))
            next_polygon = surface.polygon(next_label)

            next_center = next_polygon.circumscribing_circle().center().vector()
            
            # Bring next_center into the coordinate system of polygon
            next_center += (polygon.vertex(vertex) - next_polygon.vertex(next_vertex))

            holonomy = next_center - center

            if not holonomy:
                # Ambiguous Delaunay triangulation. The two circumscribing circles coincide.
                pass
            else:
                # Construct the start point of the boundary segment and a segment from the center to that start point.
                start_label, start_edge = (label, vertex)
                start_holonomy = center - polygon.vertex(start_edge)

                from flatsurf.geometry.euclidean import ccw
                if ccw(-polygon.edge(start_edge - 1), start_holonomy) > 0:
                    start_label, start_edge = surface.opposite_edge(start_label, (start_edge - 1) % len(polygon.vertices()))
                    polygon = surface.polygon(start_label)
                elif ccw(polygon.edge(start_edge), start_holonomy) <= 0:
                    start_label, start_edge = surface.opposite_edge(start_label, start_edge)
                    polygon = surface.polygon(start_label)
                    start_edge = (start_edge + 1) % len(polygon.vertices())

                assert ccw(-polygon.edge(start_edge), start_holonomy) <= 0 and ccw(polygon.edge(start_edge), start_holonomy) > 0, "center of Voronoi cell must be within a neighboring triangle"

                start_segment = SurfaceLineSegment(surface, start_label, polygon.vertex(start_edge), start_holonomy)

                assert not start_segment.end().is_vertex(), "boundary of a Voronoi cell cannot go through a vertex"

                segment = SurfaceLineSegment(surface, *start_segment.end().representative(), holonomy)

                boundary.append(VoronoiCellBoundarySegment(self, segment, start_segment))

            label, vertex = next_label, next_vertex

            if (label, vertex) == (initial_label, initial_vertex):
                break

            label, vertex = next_label, next_vertex

        return tuple(boundary)


class VoronoiCellDecomposition(CellDecomposition):
    r"""
    EXAMPLES::

        sage: from flatsurf import translation_surfaces, VoronoiCellDecomposition

        sage: S = translation_surfaces.regular_octagon().subdivide().codomain().delaunay_triangulation()
        sage: V = VoronoiCellDecomposition(S)

        sage: V.plot()

    ::

        sage: from flatsurf import Polygon, similarity_surfaces, VoronoiCellDecomposition

        sage: S = similarity_surfaces.billiard(Polygon(angles=[3, 4, 13])).minimal_cover("translation").erase_marked_points().codomain().delaunay_triangulation()
        sage: V = VoronoiCellDecomposition(S)
        
        sage: V.plot()

    """
    Cell = VoronoiCell

    def __init__(self, surface):
        if not surface.is_delaunay_triangulated():
            raise NotImplementedError("surface must be Delaunay triangulated")
        if not surface.is_translation_surface():
            raise NotImplementedError("surface must be a translation surface")

        super().__init__(surface)

    def __repr__(self):
        return f"Voronoi cell decomposition of {self.surface()}"

    def change(self, surface=None):
        if surface is not None:
            self = VoronoiCellDecomposition(surface)

        return self

class ApproximateWeightedCellBoundarySegment(LineBoundarySegment):
    pass


class ApproximateWeightedVoronoiPolygonCell(ConvexLineSegmentPolygonCell):
    pass


class ApproximateWeightedVoronoiCell(LineSegmentCell):
    BoundarySegment = ApproximateWeightedCellBoundarySegment
    PolygonCell = ApproximateWeightedVoronoiPolygonCell

    @cached_method
    def boundary(self):
        surface = self.surface()

        (label, coordinates) = self.center().representative()
        vertex = surface.polygon(label).get_point_position(coordinates).get_vertex()

        boundary = []

        initial_label, initial_vertex = label, vertex

        while True:
            polygon = surface.polygon(label)

            corners = self.decomposition()._split_polygon_at_vertex(label, vertex)

            assert len(corners) > 1, "cannot create segments in this polygon from a single corner"

            next_label, next_vertex = surface.opposite_edge(label, (vertex - 1) % len(polygon.vertices()))

            for corner, next_corner in zip(corners, corners[1:]):
                segment = SurfaceLineSegment(surface, label, corner, next_corner - corner)

                start_segment_holonomy = corner - polygon.vertex(vertex)
                from flatsurf.geometry.euclidean import is_anti_parallel
                if is_anti_parallel(polygon.edge(vertex - 1), start_segment_holonomy):
                    start_segment = SurfaceLineSegment(surface, next_label, surface.polygon(next_label).vertex(next_vertex), start_segment_holonomy)
                else:
                    start_segment = SurfaceLineSegment(surface, label, polygon.vertex(vertex), start_segment_holonomy)
                boundary.append(ApproximateWeightedCellBoundarySegment(self, segment, start_segment))

                # TODO: Why is this needed?
                if len(boundary) > 1:
                    if boundary[-1] == boundary[-2]:
                        boundary.pop()

            ## Now taken care of by split_polygon_at_vertex()
            ## corners_next_polygon = self.decomposition()._split_polygon_at_vertex(next_label, next_vertex)
            ## if surface(label, corners[-1]) != surface(next_label, corners_next_polygon[0]):
            ##     # The split of the polygons along the edge is not compatible,
            ##     # add a short segment along that edge to connect the last
            ##     # corner of this polygon to the first corner of the next
            ##     # polygon.
            ##     next_polygon = surface.polygon(next_label)
            ##     corner_next_polygon = corners_next_polygon[0] + (polygon.vertex(vertex) - next_polygon.vertex(next_vertex))

            ##     segment = SurfaceLineSegment(surface, label, corners[-1], corner_next_polygon - corners[-1])
            ##     start_segment = SurfaceLineSegment(surface, next_label, next_polygon.vertex(next_vertex), corners[-1] - polygon.vertex(vertex))

            ##     boundary.append(ApproximateWeightedCellBoundarySegment(self, segment, start_segment))

            label, vertex = next_label, next_vertex

            if (label, vertex) == (initial_label, initial_vertex):
                break

        # TODO: Why is this needed?
        if len(boundary) > 1:
            if boundary[-1] == boundary[0]:
                boundary.pop()

        return tuple(boundary)


class ApproximateWeightedVoronoiCellDecomposition(CellDecomposition):
    r"""
    EXAMPLES::

        sage: from flatsurf import translation_surfaces, ApproximateWeightedVoronoiCellDecomposition

        sage: S = translation_surfaces.regular_octagon().subdivide().codomain().delaunay_triangulation()
        sage: V = ApproximateWeightedVoronoiCellDecomposition(S)

        sage: V.plot()

    ::

        sage: from flatsurf import Polygon, similarity_surfaces, ApproximateWeightedVoronoiCellDecomposition

        sage: S = similarity_surfaces.billiard(Polygon(angles=[3, 4, 13])).minimal_cover("translation").erase_marked_points().codomain().delaunay_triangulation()
        sage: V = ApproximateWeightedVoronoiCellDecomposition(S)
        
        sage: V.plot()

    """
    Cell = ApproximateWeightedVoronoiCell

    def __init__(self, surface):
        if not surface.is_delaunay_triangulated():
            raise NotImplementedError("surface must be triangulated")
        if not surface.is_translation_surface():
            raise NotImplementedError("surface must be a translation surface")

        super().__init__(surface)

    def __repr__(self):
        return f"Approximate weighted Voronoi cell decomposition of {self.surface()}"

    def change(self, surface=None):
        if surface is not None:
            self = ApproximateWeightedVoronoiCellDecomposition(surface)

        return self

    def _exactify(self, x):
        R = self.surface().base_ring()
        return R(round(float(x), 3))

    @cached_method
    def _split_polygon_at_vertex(self, label, vertex):
        # Return a refined version of _split_polygon(label)[vertex] that adds points for the split polygon that happened across any edges (so each segment has a negative.)
        surface = self.surface()
        polygon = surface.polygon(label)

        split = [polygon.vertex(vertex)]

        split.extend(self._split_polygon(label)[vertex])
        
        # previous_label, previous_vertex = surface.opposite_edge(label, vertex)
        # previous_polygon = surface.polygon(previous_label)
        # previous_vertex = (previous_vertex + 1) % len(previous_polygon.vertices())
        # previous_corner = self._split_polygon(previous_label)[previous_vertex][-1]
        # previous_corner += (polygon.vertex(vertex) - previous_polygon.vertex(previous_vertex))

        # assert polygon.get_point_position(previous_corner).is_in_edge_interior() and polygon.get_point_position(previous_corner).get_edge() == vertex

        # if split[1] != previous_corner:
        #     split = split[:1] + [previous_corner] + split[1:]

        next_label, next_vertex = surface.opposite_edge(label, (vertex - 1) % len(polygon.vertices()))
        next_polygon = surface.polygon(next_label)
        next_corner = self._split_polygon(next_label)[next_vertex][0]
        next_corner += (polygon.vertex(vertex) - next_polygon.vertex(next_vertex))

        assert polygon.get_point_position(next_corner).is_in_edge_interior() and polygon.get_point_position(next_corner).get_edge() == (vertex - 1) % len(polygon.vertices())

        if split[-1] != next_corner:
            split.append(next_corner)

        refined_split = []
        for corner, next_corner in zip(split, split[1:] + split[:1]):
            assert corner != next_corner
            midpoint = (corner + next_corner) / 2
            midpoint_position = polygon.get_point_position(midpoint)
            if midpoint_position.is_in_edge_interior():
                edge = midpoint_position.get_edge()

                opposite_label, opposite_edge = surface.opposite_edge(label, edge)
                opposite_polygon = surface.polygon(opposite_label)

                edge_points = self._split_polygon_edge(opposite_label, opposite_edge)
                edge_points = [p + (polygon.vertex(edge) - opposite_polygon.vertex(opposite_edge + 1)) for p in edge_points]

                assert all(polygon.get_point_position(p).get_edge() == edge for p in edge_points)

                from flatsurf.geometry.euclidean import time_on_segment
                edge_points = [p for p in edge_points if time_on_segment((corner, next_corner), p) not in [0, 1, None]]
                for p in sorted(edge_points, key=lambda p: time_on_segment((corner, next_corner), p)):
                    assert p not in refined_split
                    assert p != next_corner
                    refined_split.append(p)

            assert next_corner not in refined_split
            refined_split.append(next_corner)

        return tuple(refined_split[:-1])

    @cached_method
    def _split_polygon_edge(self, label, edge):
        polygon = self.surface().polygon(label)

        points = []

        for split in self._split_polygon(label):
            for point in split:
                position = polygon.get_point_position(point)
                if not position.is_in_edge_interior():
                    continue
                if position.get_edge() != edge:
                    continue
                if point not in points:
                    points.append(point)

        return tuple(points)

    @cached_method
    def _split_polygon(self, label):
        surface = self.surface()
        polygon = surface.polygon(label)
        nvertices = len(polygon.vertices())
        assert nvertices == 3

        weights = [float(surface(label, v).radius_of_convergence()) for v in range(nvertices)]

        splits = [self._split_segment((polygon.vertex(v), polygon.vertex(v + 1)), polygon.vertex(v + 2), weights[v:] + weights[:v]) for v in range(nvertices)]

        lens = [len(split) for split in splits]

        if lens == [1, 1, 1]:
            # The edges are far away from the opposite corners.
            # A---+---C
            # |  /|  /
            # | / | /
            # |/C |/
            # +---+
            # |  /
            # | /
            # |/
            # B
            # We give each small triangle in this picture to its vertex and
            # split the central triangle C between the three vertices. (So each
            # cell will be a quadrilateral.)
            center = sum(splits[v][0] * self._exactify((weights[v] + weights[(v+1) % nvertices]) / (2 * sum(weights))) for v in range(nvertices))

            return [[splits[0][0], center, splits[2][0]], [splits[1][0], center, splits[0][0]], [splits[2][0], center, splits[1][0]]]
        if lens == [2, 1, 1]:
            return [[splits[0][0], splits[2][0]], [splits[1][0], splits[0][1]], [splits[2][0], splits[0][0], splits[0][1], splits[1][0]]]
        if lens == [1, 1, 2]:
            return [[splits[0][0], splits[2][1]], [splits[1][0], splits[2][0], splits[2][1], splits[0][0]], [splits[2][0], splits[1][0]]]
        if lens == [1, 2, 1]:
            return [[splits[0][0], splits[1][0], splits[1][1], splits[2][0]], [splits[1][0], splits[0][0]], [splits[2][0], splits[1][1]]]

        raise NotImplementedError(f"non-trivial split for {label} with {lens}")

    def _split_segment(self, AB, C, weights):
        A, B = AB
        wA, wB, wC = weights

        def circle_of_apollonius(P, Q, w, v):
            if w == v:
                raise NotImplementedError("circle of Apollonius is a line")

            C = (w * Q + v * P) / (w + v)
            D = (v * P - w * Q) / (v - w)

            center = (C + D) / 2

            radius = (D - C).norm() / 2

            return center, radius

        def circle_segment_intersection(center, radius, segment):
            # Shift the center to the origin.
            segment = segment[0] - center, segment[1] - center

            # Let P(t) be the point on the segment at time t in [0, 1].
            # We solve the quadratic equation for |P(t)| = radius^2.

            a = (segment[1][0] - segment[0][0])**2 + (segment[1][1] - segment[0][1])**2
            b = 2*(segment[0][0]*segment[1][0] + segment[0][1]*segment[1][1] - segment[0][0]**2 - segment[0][1]**2)
            c = segment[0][0]**2 + segment[0][1]**2 - radius**2

            from math import sqrt
            for t in [
                (-b + sqrt(float(b**2 - 4*a*c))) / (2*a),
                (-b - sqrt(float(b**2 - 4*a*c))) / (2*a),
            ]:
                if 0 <= t <= 1:
                    return self._exactify(t)

            assert False, f"intersection point must be on segment but time {t} is not in the unit interval, i.e., segment {segment} does not intersect circle at origin with radius {radius}"

        A_equal_B = A + (B - A) * self._exactify(wA / (wA + wB))
        if ((A - A_equal_B) / wA).norm() > ((C - A_equal_B) / wC).norm():
            # C is closer to A on the segment AB than B.
            # Construct the circle of Apollonius with foci A and C.
            if wA == wC:
                # Circle of Apollonius is a line.
                AC = C - A
                midpoint = (C + A) / 2
                from sage.all import vector
                orthogonal_ray = (midpoint, vector((AC[1], -AC[0])))
                
                from flatsurf.geometry.euclidean import ray_segment_intersection
                A_equal_C = ray_segment_intersection(*orthogonal_ray, (A, B))
                assert A_equal_C is not None
            else:
                center, radius = circle_of_apollonius(A, C, wA, wC)
                t = circle_segment_intersection(center, radius, (A, B))
                A_equal_C = (1 - t) * A + t * B
            if wB == wC:
                # Circle of Apollonius is a line.
                BC = C - B
                midpoint = (C + B) / 2
                from sage.all import vector
                orthogonal_ray = (midpoint, vector((-BC[1], BC[0])))
                
                from flatsurf.geometry.euclidean import ray_segment_intersection
                B_equal_C = ray_segment_intersection(*orthogonal_ray, (B, A))
                assert B_equal_C is not None
            else:
                center, radius = circle_of_apollonius(B, C, wB, wC)
                t = circle_segment_intersection(center, radius, (B, A))
                B_equal_C = (1 - t) * B + t * A

            return [A_equal_C, B_equal_C]

        return [A_equal_B]

# TODO: Move to surface objects.
class SurfaceLineSegment:
    def __init__(self, surface, label, start, holonomy):
        if not holonomy:
            raise ValueError

        polygon = surface.polygon(label)

        position = polygon.get_point_position(start)
        if position.is_outside():
            raise ValueError(f"start point of segment must be in the polygon but {start} is not in {polygon}")
        if position.is_in_edge_interior():
            edge = position.get_edge()
            from flatsurf.geometry.euclidean import ccw, is_anti_parallel
            if ccw(polygon.edge(edge), holonomy) < 0 or (ccw(polygon.edge(edge), holonomy) == 0 and is_anti_parallel(polygon.edge(edge), holonomy)):
                start -= polygon.vertex(edge)
                label, edge = surface.opposite_edge(label, edge)
                polygon = surface.polygon(label)
                start += polygon.vertex(edge + 1)
        if position.is_vertex():
            vertex = position.get_vertex()
            from flatsurf.geometry.euclidean import ccw
            if ccw(polygon.edge(vertex), holonomy) >= 0 and ccw(-polygon.edge(vertex - 1), holonomy) < 0:
                # holonomy points into the polygon
                pass
            else:
                if surface(label, vertex).angle() != 1:
                    raise ValueError("holonomy must point into the polygon at a singular point")
                raise NotImplementedError("cannot handle holonomy vector points out of the polygon at a marked point yet")

        self._surface = surface
        self._label = label
        self._start = start
        self._holonomy = holonomy

        # TODO: Instead, copy if mutable.
        self._start.set_immutable()
        self._holonomy.set_immutable()

    def surface(self):
        return self._surface

    def start(self):
        return self._surface(*self.start_representative())

    def start_representative(self):
        return self._label, self._start

    def an_inner_point(self):
        segment = self.split()[0]
        return self.surface()(segment[0], (segment[1][0] + segment[1][1]) / 2)

    def end(self):
        return self._surface(*self.end_representative())

    def end_representative(self):
        end = self.split()[-1]
        return end[0], end[1][1]

    def _flow_to_exit(self):
        r"""
        Split this segment into a segment that lies entirely in its source
        polygon, and a rest segment.

        If this segment lies entirely in its source segment, returns ``None``
        as the rest.

        The initial part is returned as a segment label and a Euclidean segment
        in the corresponding polygon.
        """
        # TDOO: This generic flowing foo should not probably be implemented in a better place.
        surface = self.surface()
        polygon = surface.polygon(self._label)

        end = self._start + self._holonomy
        end.set_immutable()

        if polygon.get_point_position(end).is_outside():
            exit = polygon.flow_to_exit(self._start, self._holonomy)
            exit.set_immutable()

            Δ = exit - self._start

            rest = SurfaceLineSegment(surface, self._label, exit, self._holonomy - Δ)

            return (self._label, (self._start, exit)), rest

        return (self._label, (self._start, end)), None
        
    def __bool__(self):
        return bool(self._holonomy)

    def holonomy(self, start=None):
        if start is None:
            start = self._label, self._start

        if start == (self._label, self._start):
            return self._holonomy

        # TODO: Why should anything be different?
        return self._holonomy
        raise NotImplementedError("cannot determine holonomy at vertex yet")

    def plot(self):
        # TODO: Better plotting
        return self.surface().plot() + self.start().plot(size=100) + self.end().plot(size=50)

    def __repr__(self):
        return f"{self.start()}→{self.end()}"

    def __neg__(self):
        return SurfaceLineSegment(self._surface, *self.end_representative(), -self._holonomy)

    def split(self, label=None):
        r"""
        Return this segment split into subsegments that live entirely in
        polygons of the surface.

        The subsegments are returned as triples (polygon label, segment in the
        euclidean plane, holonomy from the start point of the segment to the
        start point of the subsegment).

        If ``label`` is set, only return the subsegments that are in the
        polygon with ``label``.
        """
        surface = self.surface()

        segments = []

        from sage.all import vector
        holonomy = vector((0, 0))

        while True:
            if self is None:
                break

            (lbl, segment), self = self._flow_to_exit()

            if label is None or label == lbl:
                holonomy.set_immutable()
                segments.append((lbl, segment, holonomy))

            holonomy += (segment[1] - segment[0])

        return segments

    def __eq__(self, other):
        if self.start() != other.start():
            return False

        start = self.start_representative()

        return self.holonomy(start) == other.holonomy(start)

    def __hash__(self):
        return hash((self.start(), self._holonomy))

    def plot(self, graphical_surface=None):
        plot_surface = graphical_surface is None

        if graphical_surface is None:
            graphical_surface = self.surface().graphical_surface(edge_labels=False, polygon_labels=False)

        plot = []
        if plot_surface:
            plot.append(graphical_surface.plot())

        subsegments = self.split()
        for s, (label, subsegment, _) in enumerate(subsegments):
            polygon = self.surface().polygon(label)

            graphical_polygon = graphical_surface.graphical_polygon(label)

            # TODO: This should be implemented in graphical polygon probably.
            # TODO: This assumes that the surface is a translation surface.
            vertex = graphical_polygon.transformed_vertex(0)
            graphical_subsegment = [vertex + (subsegment[0] - polygon.vertex(0)), vertex + (subsegment[1] - polygon.vertex(0))]

            if s == len(subsegments) - 1:
                from sage.all import arrow2d
                plot.append(arrow2d(*graphical_subsegment, arrowsize=1.5, width=.4, color="green"))
            else:
                from sage.all import line2d
                plot.append(line2d(graphical_subsegment))

        return sum(plot)
