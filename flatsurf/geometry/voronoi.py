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
        return [(OrientedSegment(*s), polygon_cell) for (label, subsegment) in segment.split() for (s, polygon_cell) in self._split_segment_at_polygon_cells(label, subsegment)]

    def _split_segment_at_polygon_cells(self, label, segment):
        start, end = segment

        if start == end:
            return []

        for start_cell in self.polygon_cells(label=label):
            if not start_cell.contains_point(start):
                continue

            polygon = start_cell._polygon()

            try:
                exit = polygon.flow_to_exit(start, end - start)
            except ValueError:
                continue

            return [((start, exit), start_cell)] + self._split_segment_at_polygon_cells(label, (exit, end))

        assert False

        ## # TODO: Just use flow_to_exit() in the polygons instead.

        ## polygon = self.surface().polygon(label)

        ## polygon_cells = [polygon_cell for polygon_cell in self.polygon_cells(label=label) if polygon_cell.contains_point(start)]

        ## boundaries = {boundary: polygon_cell for polygon_cell in polygon_cells for boundary in polygon_cell.boundaries()}

        ## events = {}

        ## from collections import namedtuple
        ## Event = namedtuple("Event", ("enter", "boundary", "intersection"))

        ## ray = (start, end - start)

        ## # TODO: Instead we should intersect with the boundary shape (i.e., the entire polygon and not roll our own intersection algorithm yet again.)
        ## for boundary in boundaries:
        ##     # TODO: This makes assumptions about boundaries being segments.
        ##     from flatsurf.geometry.euclidean import ray_segment_intersection, time_on_ray
        ##     intersection = ray_segment_intersection(*ray, boundary)
        ##     if intersection is None:
        ##         continue

        ##     def record_event(intersection, boundary, enter: bool):
        ##         if polygon.get_point_position(intersection).is_outside():
        ##             return

        ##         intersection_time, length = time_on_ray(*ray, intersection)

        ##         if intersection_time > length:
        ##             return

        ##         if intersection_time == 0 and not enter:
        ##             return

        ##         if intersection_time not in events:
        ##             events[intersection_time] = []
        ##         events[intersection_time].append(Event(enter, boundary, intersection))

        ##     if intersection[0].parent() is start.parent():
        ##         # Intersection in a segment (and not only a point.)
        ##         record_event(intersection[0], boundary, enter=True)
        ##         record_event(intersection[0], boundary, enter=False)
        ##     else:
        ##         from flatsurf.geometry.euclidean import ccw
        ##         enter = ccw(boundary[1] - boundary[0], end - start)
        ##         assert enter != 0
        ##         record_event(intersection, boundary, enter=enter > 0)

        ## if not events:
        ##     assert len(polygon_cells) == 1
        ##     return segment, next(iter(polygon_cells))

        ## first_events = events[min(events)]

        ## if len(events) == 1:
        ##     if min(events) == 0:
        ##         assert all(event.enter for event in first_events)

        ##         if len(first_events) != 1:
        ##             raise NotImplementedError

        ##         event = next(iter(first_events))
        ##         polygon_cell = boundaries[event.boundary]
        ##         assert polygon_cell.contains_point(end)
        ##         return [(segment, polygon_cell)]

        ##     assert all(not event.enter for event in first_events)

        ##     if len(first_events) != 1:
        ##         raise NotImplementedError

        ##     event = first_events[0]
        ##     intersection = event.intersection
        ##     polygon_cell = boundaries[event.boundary]
        ##     return [((start, intersection), polygon_cell)] + self._split_segment_at_polygon_cells(label, (intersection, end))

        ## raise NotImplementedError

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


class Cell:
    BoundarySegment = None

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

    def radius(self):
        r"""
        Return the distance to the point in the cell furthest from the
        :meth:`center` of the cell.
        """
        raise NotImplementedError("this cell decomposition does not implement radius()")

    # TODO: Maybe add a label parameter to get only the ones that live in label?
    def polygon_cells(self):
        r"""
        Return the restrictions to the polygons that contain an inner point of
        the cell.
        """
        raise NotImplementedError("this cell decomposition cannot restrict cells to polygons yet")

    @cached_method
    def boundary(self):
        r"""
        Return the boundary of this cell in a counterclockwise walk around the
        center.
        """
        boundary = []

        corners = self._corners()

        for c, corner in enumerate(corners):
            next_corner = corners[(c + 1) % len(corners)]
            segment = self.BoundarySegment(self, corner, next_corner)
            boundary.append(segment)

        return tuple(boundary)

    def _radius_corners(self):
        r"""
        Return the distance from the center to the furthest corner in this
        cell.
        """
        norm = self.surface().euclidean_plane().norm()
        return max(norm.from_vector(self.surface().polygon(label).vertex(vertex) - corner) for (label, vertex, corner) in self._corners())

    def _corners(self):
        raise NotImplementedError("this cell decomposition cannot compute segments to each cells corners yet")

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

        for polygon_cell in self.polygon_cells():
            plot.append(polygon_cell.plot(graphical_polygon=graphical_surface.graphical_polygon(polygon_cell.label())))

        return sum(plot)


class BoundarySegment:
    def __init__(self, cell, corner, next_corner):
        self._cell = cell
        self._corners = (corner, next_corner)

    def surface(self):
        return self._cell.surface()

    def cell(self):
        return self._cell

    def __bool__(self):
        return bool(self.segment())

    def __neg__(self):
        return self.cell().decomposition()._neg_boundary_segment(self)

    def segment(self):
        raise NotImplementedError("this cell decomposition cannot make its boundary segments explicit yet")

    def _polygon_cell_boundaries(self):
        boundaries = {}

        segments = self.segment().split()

        for polygon_cell in self.cell().polygon_cells():
            for (label, segment) in segments:
                if label == polygon_cell.label() and segment in polygon_cell.boundaries():
                    boundaries[(polygon_cell.label(), segment)] = polygon_cell


        assert len(segments) == len(boundaries)

        return boundaries

    def polygon_cell_boundaries(self):
        r"""
        Return this boundary segment as a sequence of subsegments that live
        entirely in a :class:`PolygonCell`.

        The subsegments are returned as triples (polygon cell, segment in the
        Euclidean plane, opposite polygon cell).
        """
        boundaries = self._polygon_cell_boundaries()
        opposite_boundaries = (-self)._polygon_cell_boundaries()

        polygon_cell_boundaries = []

        for (label, (start, end)), polygon_cell in boundaries.items():
            polygon_cell_boundaries.append((
                polygon_cell,
                (start, end),
                opposite_boundaries[(label, (end, start))]))

        return polygon_cell_boundaries

    def __eq__(self, other):
        if not isinstance(other, BoundarySegment):
            return False

        return self.cell() == other.cell() and self.segment() == other.segment()

    def __hash__(self):
        return hash(self.segment())

    def __repr__(self):
        return f"Boundary at {self.segment()}"


class PolygonCell:
    def __init__(self, cell, label, center, boundaries):
        self._cell = cell
        self._label = label
        self._center = center

        for rotation in range(len(boundaries)):
            boundaries = boundaries[1:] + boundaries[:1]
            if all(boundary[1] == next_boundary[0] for (boundary, next_boundary) in zip(boundaries, boundaries[1:])):
                break
        else:
            raise NotImplementedError("boundaries must be connected")

        self._boundaries = boundaries

    def label(self):
        return self._label

    def polygon(self):
        return self.cell().surface().polygon(self._label)

    def surface(self):
        return self.cell().surface()

    def cell(self):
        return self._cell

    def center(self):
        r"""
        Return the point (not necessarily of the polygon) that should be
        considered the center of this cell.
        """
        return self._center

    # TODO: Unify naming of boundary() and boundaries()
    def boundaries(self):
        return self._boundaries

    def contains_point(self, point):
        return self._polygon().contains_point(point)

    def _polygon(self):
        r"""
        Return a polygon (subset of :meth:`polygon`) that describes this cell.
        """
        from flatsurf import Polygon
        polygon = self.polygon()

        for point in [self._boundaries[0][0], self._boundaries[-1][1]]:
            position = polygon.get_point_position(point)
            if position.is_vertex():
                continue

            assert position.is_in_edge_interior()

            edge = position.get_edge()
            vertices = list(polygon.vertices())
            vertices = vertices[:edge + 1] + [point] + vertices[edge + 1:]

            polygon = Polygon(vertices=vertices)

        start = polygon.get_point_position(self._boundaries[0][0]).get_vertex()
        end = polygon.get_point_position(self._boundaries[-1][1]).get_vertex()

        vertices = [boundary[0] for boundary in self._boundaries]

        while True:
            vertices.append(polygon.vertex(end))
            end = (end + 1) % len(polygon.vertices())
            if end == start:
                break

        return Polygon(vertices=vertices)

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

        d = self.surface()(self.label(), self.center()).angle()

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

    def __eq__(self, other):
        if not isinstance(other, PolygonCell):
            return False

        return self.cell() == other.cell() and self.label() == other.label() and self.center() == other.center()

    def __hash__(self):
        return hash((self.label(), self.center()))

    def __repr__(self):
        return f"{self.cell()} restricted to polygon {self.label()} at {self._center}"

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
        plot = sum(OrientedSegment(*segment).translate(shift).plot(point=False) for segment in self.boundaries())

        if plot_polygon:
            plot = graphical_polygon.plot_polygon() + plot

        return plot


class VoronoiCellBoundarySegment(BoundarySegment):
    def segment(self):
        return SurfaceLineSegment(self.surface(), *self._start(), self._holonomy())

    def _start(self):
        ((label, center, corner), _) = self._corners

        surface = self.surface()
        polygon = surface.polygon(label)

        if polygon.get_point_position(corner).is_outside():
            corner -= polygon.vertex(center)
            assert corner, "corner cannot be at a vertex"

            from flatsurf.geometry.euclidean import ccw
            if ccw(polygon.edge(center), corner) > 0:
                label, center = surface.opposite_edge(label, center)
                polygon = surface.polygon(label)
                center = (center + 1) % len(polygon.vertices())
            else:
                assert ccw(-polygon.edge(center - 1), corner) > 0, "center of circumcircle cannot be opposite of vertex edge"
                label, center = surface.opposite_edge(label, (center - 1) % len(polygon.vertices()))
                polygon = surface.polygon(label)

            corner += polygon.vertex(center)

        assert not polygon.get_point_position(corner).is_outside(), "center of circumcircle must be in this or a neighboring triangle"

        return (label, corner)

    def _holonomy(self):
        ((label, center, corner), (next_label, next_center, next_corner)) = self._corners

        surface = self.surface()

        polygon = surface.polygon(next_label)
        next_corner -= polygon.vertex(next_center)

        while next_label != label:
            next_label, next_center = surface.opposite_edge(next_label, next_center)
            polygon = surface.polygon(next_label)
            next_center = (next_center + 1) % len(polygon.vertices())

        next_corner += polygon.vertex(next_center)

        return next_corner - corner


class VoronoiCell(Cell):
    BoundarySegment = VoronoiCellBoundarySegment

    def radius(self):
        return self._radius_corners()

    @cached_method
    def _corners(self):
        # Return (label, vertex, corner) where label, vertex single out the
        # polygon and its vertex from which the corner can be reached, and
        # corner is a vector from that vertex to the corner.

        surface = self.surface()

        (label, coordinates) = self.center().representative()
        vertex = surface.polygon(label).get_point_position(coordinates).get_vertex()

        corners = []

        # TODO: This is a very common operation (walking around a vertex) and
        # should probably be implemented generically.
        initial_label, initial_vertex = label, vertex
        while True:
            polygon = surface.polygon(label)

            corners.append((label, vertex, polygon.circumscribing_circle().center().vector()))

            if len(corners) > 1 and not self.BoundarySegment(self, corners[-2], corners[-1]):
                corners.pop()

            label, vertex = surface.opposite_edge(label, (vertex - 1) % len(polygon.vertices()))
            if (label, vertex) == (initial_label, initial_vertex):
                break

        if not self.BoundarySegment(self, corners[-1], corners[0]):
            corners.pop()

        if len(corners) < 2:
            raise NotImplementedError("cannot create cell from less than 2 corners")

        return tuple(corners)

    @cached_method
    def polygon_cells(self):
        surface = self.surface()

        cells = {}
        for boundary in self.boundary():
            (label, center, corner) = boundary._corners[0]
            for (lbl, segment) in boundary.segment().split():
                assert label == lbl

                cell = (label, center)
                if cell not in cells:
                    cells[cell] = []
                cells[cell].append(segment)

                label, center = surface.opposite_edge(label, (center - 1) % len(surface.polygon(label).vertices()))

        return tuple(VoronoiPolygonCell(self, label, center, tuple(boundaries)) for ((label, center), boundaries) in cells.items())


class VoronoiPolygonCell(PolygonCell):
    def __init__(self, cell, label, vertex, boundaries):
        super().__init__(cell, label, cell.surface().polygon(label).vertex(vertex), boundaries)

    def contains_segment(self, segment):
        return self.contains_point(segment.start()) and self.contains_point(segment.end())


class VoronoiCellDecomposition(CellDecomposition):
    Cell = VoronoiCell

    def __init__(self, surface):
        if not surface.is_delaunay_triangulated():
            raise NotImplementedError("surface must be Delaunay triangulated")
        if not surface.is_translation_surface():
            raise NotImplementedError("surface must be a translation surface")

        super().__init__(surface)

    def __repr__(self):
        return f"Voronoi cell decomposition of {self.surface()}"


# TODO: Move to surface objects.
class SurfaceLineSegment:
    def __init__(self, surface, label, start, holonomy):
        if not holonomy:
            raise ValueError

        polygon = surface.polygon(label)

        position = polygon.get_point_position(start)
        if position.is_outside():
            raise ValueError("start point of segment must be in the polygon")
        if position.is_in_edge_interior():
            edge = position.get_edge()
            from flatsurf.geometry.euclidean import ccw
            if ccw(polygon.edge(edge), holonomy) < 0:
                start -= polygon.vertex(edge + 1)
                label, edge = surface.opposite_edge(label, edge)
                polygon = surface.polygon(label)
                start += polygon.vertex(edge)
        if position.is_vertex():
            vertex = position.get_vertex()
            from flatsurf.geometry.euclidean import ccw
            if ccw(polygon.edge(vertex), holonomy) >= 0 and ccw(-polygon.edge(vertex - 1), holonomy) < 0:
                # holonomy points into the polygon
                pass
            else:
                raise NotImplementedError

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

        raise NotImplementedError("cannot determine holonomy at vertex yet")

    def __repr__(self):
        return f"{self.start()}→{self.end()}"

    def __neg__(self):
        return SurfaceLineSegment(self._surface, *self.end_representative(), -self._holonomy)

    def split(self):
        r"""
        Return this segment split into subsegments that live entirely in
        polygons of the surface.

        The subsegments are returned as pairs (polygon label, segment in the
        euclidean plane).
        """
        segments = []

        while True:
            (label, segment), self = self._flow_to_exit()
            segments.append((label, segment))
            if self is None:
                break

        return segments

    def __eq__(self, other):
        if self.start() != other.start():
            return False

        start = self.start_representative()

        return self.holonomy(start) == other.holonomy(start)

    def __hash__(self):
        return hash((self.start(), self._holonomy))


### # TODO: Drop all use of coordinates. Use EuclideanPolygonPoint instead of (label, coordinates)
### # TODO: SageMath already has a VoronoiDiagram so maybe not use that name. Explain that we cannot use the builtin version since it has no weights.
### # TODO: The classical Voronoi cells are not correct currently.
### 
### from sage.misc.cachefunc import cached_method
### 
### 
### # TODO: Make these unique for each surface (but they cannot be unique
### # representation because the same surface might have several representatiions)
### class VoronoiDiagram:
###     r"""
###     ALGORITHM:
### 
###     Internally, a Voronoi diagram is composed by :class:`VoronoiCell`s. Each
###     cell is the union of :class:`VoronoiPolygonCell` which represents the
###     intersection of a Voronoi cell with an Euclidean polygon making up the
###     surface. To compute these cells on the level of polygons, we determine the
###     segments that bound the Voronoi cell where each such segment is half way (up
###     to optional weights) between the centers of two Voronoi cells.
### 
###     TODO: Reducing the computation to the level of polygons is not correct in
###     general. When centers are close to edges, the polygon on the other side of
###     the edge does not see that center. Can this be fixed? Or should we just say
###     that this is an approximation of the Voronoi diagram? And is it correct
###     without weighing and with centers at the vertices of a Delaunay cell
###     decomposition?
### 
###     TODO: Explain how these cells are not minimal. The same point that exists
###     twice in a polygon produces a boundary of the Voronoi cell with itself. Add
###     a method that gets rid of these boundaries (and throw a NotImplementedError
###     since it is unclear how to represent cells that contain an entire polygon?)
### 
###     EXAMPLES::
### 
###         sage: from flatsurf.geometry.voronoi import VoronoiDiagram  # random output due to deprecation warnings
###         sage: from flatsurf import translation_surfaces
###         sage: S = translation_surfaces.regular_octagon()
###         sage: center = S(0, S.polygon(0).centroid())
###         sage: S = S.insert_marked_points(center).codomain()
###         sage: V = VoronoiDiagram(S, S.vertices())
###         sage: V.plot()
###         Graphics object consisting of 73 graphics primitives
###         sage: V = VoronoiDiagram(S, S.vertices(), weight="radius_of_convergence")
###         sage: V.plot()
###         Graphics object consisting of 73 graphics primitives
### 
###     The same Voronoi diagram but starting from a more complicated description
###     of the octagon::
### 
###         sage: from flatsurf import Polygon, translation_surfaces, polygons
###         sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###         sage: S = translation_surfaces.regular_octagon()
###         sage: S = S.subdivide().codomain().delaunay_decomposition().codomain()
###         sage: V = VoronoiDiagram(S, S.vertices())
###         sage: V.plot()
###         Graphics object consisting of 73 graphics primitives
### 
###         sage: V = VoronoiDiagram(S, S.vertices(), weight="radius_of_convergence")
###         sage: V.plot()
###         Graphics object consisting of 73 graphics primitives
### 
###     """
### 
###     def __init__(self, surface, points, weight=None):
###         if surface.is_mutable():
###             raise ValueError("surface must be immutable")
### 
###         if weight is None:
###             weight = "classical"
### 
###         self._surface = surface
###         self._centers = frozenset(points)
###         self._weight = weight
### 
###         if not surface.vertices().issubset(self._centers):
###             # This probably essentially works. However, we cannot represent
###             # cells that contain an entire polygon yet, so this is not
###             # implemented in general.
###             raise NotImplementedError("can only compute Voronoi diagrams when all vertices are centers")
### 
###         if set(surface.vertices()) != set(self._centers):
###             raise NotImplementedError("non-vertex centers are not supported anymore")
### 
###         self._cells = {center: VoronoiCell(self, center) for center in self._centers}
### 
###     def surface(self):
###         r"""
###         Return the surface which this diagram is breaking up into cells.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.surface() is S
###             True
### 
###         """
###         return self._surface
### 
###     def plot(self, graphical_surface=None):
###         r"""
###         Return a graphical representation of this Voronoi diagram.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.plot()
###             Graphics object consisting of 73 graphics primitives
### 
###         The underlying surface is not plotted automatically when it is provided
###         as a keyword argument::
### 
###             sage: V.plot(graphical_surface=S.graphical_surface())
###             Graphics object consisting of 48 graphics primitives
### 
###         """
###         plot_surface = graphical_surface is None
### 
###         if graphical_surface is None:
###             graphical_surface = self._surface.graphical_surface(edge_labels=False, polygon_labels=False)
### 
###         plot = []
###         if plot_surface:
###             plot.append(graphical_surface.plot())
### 
###         for cell in self._cells.values():
###             plot.append(cell.plot(graphical_surface))
### 
###         return sum(plot)
### 
###     def cell(self, point):
###         r"""
###         Return a Voronoi cell that contains ``point``.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: insert = S.insert_marked_points(center)
###             sage: S = insert.codomain()
###             sage: center = insert(center)
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.cell(center)
###             sage: cell.contains_point(center)
###             True
### 
###         """
###         return next(iter(self.cells(point)))
### 
###     def cells(self, point=None):
###         r"""
###         Return the Voronoi cells that contain ``point``.
### 
###         If no ``point`` is given, return all cells.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cells = V.cells(S((0, 0), (1/2, 0)))
###             sage: list(cells)
###             [Voronoi cell at Vertex 0 of polygon (0, 0)]
### 
###         """
###         for cell in self._cells.values():
###             if point is None or cell.contains_point(point):
###                 yield cell
### 
###     def corners(self):
###         corners = []
###         for cell in self.cells():
###             corners.extend(cell.corners())
###         return set(corners)
### 
###     def polygon_cell(self, label, coordinates):
###         r"""
###         Return a Voronoi cell that contains the point ``coordinates`` of the
###         polygon with ``label``.
### 
###         This is to :meth:`polygon_cells` what :meth:`cell` is to :meth:`cells`.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: insert = S.insert_marked_points(center)
###             sage: S = insert.codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.polygon_cell((0, 0), (1/2, 1))
###             Voronoi cell in polygon (0, 0) at (1/2, 1/2*a + 1/2)
### 
###         """
###         return next(iter(self.polygon_cells(label, coordinates)))
### 
###     def polygon_cells(self, label, coordinates):
###         r"""
###         Return the Voronoi cells that contain the point ``coordinates`` of the
###         polygon with ``label``.
### 
###         This is similar to :meth:`cells` but it returns the
###         :class:`VoronoiPolygonCell`s instead of the full :class:`VoronoiCell`s.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.polygon_cells((0, 0), (1/2, 0))
###             [Voronoi cell in polygon (0, 0) at (0, 0), Voronoi cell in polygon (0, 0) at (1, 0)]
### 
###         """
###         if self.surface().polygon(label).get_point_position(coordinates).is_outside():
###             raise ValueError(f"coordinates must be inside polygon but {coordinates} are not in polygon with label {label}")
###         return self._diagram_polygon(label).polygon_cells(coordinates)
### 
###     @cached_method
###     def _diagram_polygon(self, label):
###         return VoronoiDiagram_Polygon(self, label, self._weight)
### 
###     def split_segment(self, label, segment):
###         r"""
###         Return the ``segment`` split into shorter segments that each lie in a
###         single Voronoi cell.
### 
###         The segments are returned as a dict mapping the Voronoi cell to the
###         shorter segments.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf.geometry.euclidean import OrientedSegment
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.split_segment((0, 0), OrientedSegment((0, 0), (1/2, 1)))
###             {Voronoi cell in polygon (0, 0) at (0, 0): OrientedSegment((0, 0), (-1/2*a + 1, -a + 2)),
###              Voronoi cell in polygon (0, 0) at (1/2, 1/2*a + 1/2): OrientedSegment((-1/2*a + 1, -a + 2), (1/2, 1))}
### 
###         When there are multiple ways to split the segment, the choice is made
###         randomly. Here, the first segment, is on the boundary of two cells::
### 
###             sage: V.split_segment((0, 0), OrientedSegment((1/2, 0), (1/2, 1)))
###             {Voronoi cell in polygon (0, 0) at (...): OrientedSegment((1/2, 0), (1/2, 1/2)),
###              Voronoi cell in polygon (0, 0) at (1/2, 1/2*a + 1/2): OrientedSegment((1/2, 1/2), (1/2, 1))}
### 
###         TESTS::
### 
###             sage: from flatsurf import similarity_surfaces, Polygon
###             sage: S = similarity_surfaces.billiard(Polygon(angles=[3/8, 1/2, 1/8], lengths=[1/2])).minimal_cover('translation')
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: centers = [S(label, S.polygon(label).centroid()) for label in S.labels()]
###             sage: S = S.insert_marked_points(*centers).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices(), weight="radius_of_convergence")
### 
###             sage: from flatsurf.geometry.euclidean import OrientedSegment
###             sage: c = S.base_ring().gen()
### 
###             # sage: segment = OrientedSegment((1/2, -c/2 - 1/2), (1/2, 0))
###             # sage: V.split_segment((1, 1, 0), segment)
###             # {Voronoi cell in polygon (1, 1, 0) at (1/2, -1/2*c0 - 1/2): OrientedSegment((1/2, -1/2*c0 - 1/2), (1/2, ...)),
###             #  Voronoi cell in polygon (1, 1, 0) at (1/2, 0): OrientedSegment((1/2, ...), (1/2, 0)),
###             #  Voronoi cell in polygon (1, 1, 0) at (1/3, -1/6*c0 - 1/6): OrientedSegment((1/2, ...), (1/2, ...))}
### 
###             # sage: segment = OrientedSegment((c/4, c/4), (-1/2, c/2 + 1/2))
###             # sage: V.split_segment((0, c/2, c/2), segment)
###             # {Voronoi cell in polygon (0, 1/2*c0, 1/2*c0) at (-1/2, 1/2*c0 + 1/2): OrientedSegment((..., ...), (-1/2, 1/2*c0 + 1/2)),
###             #  Voronoi cell in polygon (0, 1/2*c0, 1/2*c0) at (1/12*c0 - 1/6, 1/4*c0 + 1/6): OrientedSegment((1/4*c0 - ..., 1/4*c0 + ...), (..., ...)),
###             #  Voronoi cell in polygon (0, 1/2*c0, 1/2*c0) at (1/4*c0, 1/4*c0): OrientedSegment((1/4*c0, 1/4*c0), (1/4*c0 - ..., 1/4*c0 + ...))}
### 
###         """
###         segments = {}
### 
###         a = segment.start()
###         b = segment.end()
### 
###         while a != b:
###             a_cells = set(self.polygon_cells(label, a))
###             b_cells = set(self.polygon_cells(label, b))
### 
###             shared_cells = set(a_cells).intersection(b_cells)
### 
###             if shared_cells:
###                 # Both endpoints are in the same Voronoi cell.
###                 cell = next(iter(shared_cells))
###                 assert cell not in segments
### 
###                 from flatsurf.geometry.euclidean import OrientedSegment
###                 segments[cell] = OrientedSegment(a, b)
### 
###                 break
### 
###             # Bisect the segment [a, b] until we find a point c such that [a,
###             # c] crosses the boundary between two neighboring Voronoi cells.
###             c = b
###             while True:
###                 split = self._split_segment_at_boundary_point(label, a, c)
### 
###                 if split is None:
###                     c = (a + c) / 2
###                     continue
### 
###                 segment, _ = split
### 
###                 cell = next(iter(cell for cell in a_cells if cell.contains_segment(segment)))
### 
###                 assert cell not in segments
###                 segments[cell] = segment
###                 a = segment.end()
###                 break
### 
###         return segments
### 
###     def _split_segment_at_boundary_point(self, label, a, b):
###         r"""
###         Return the point that lies on the segment [a, b] and on the boundary
###         between the Voronoi cells containing a and b respectively.
### 
###         TODO: This does not return a point.
### 
###         Return ``None`` if no such point exists.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V._split_segment_at_boundary_point((0, 0), (0, 0), (1, 0))
###             (OrientedSegment((0, 0), (1/2, 0)), OrientedSegment((1/2, 0), (1, 0)))
###             sage: V._split_segment_at_boundary_point((0, 0), (0, 0), (1/2, 1))
###             (OrientedSegment((0, 0), (-1/2*a + 1, -a + 2)), OrientedSegment((-1/2*a + 1, -a + 2), (1/2, 1)))
### 
###         ::
### 
###             # TODO: We need an example with two cells in a triangle that do not touch.
###             # sage: V._split_segment_at_boundary_point((0, 0), (0, 0), (5/4, 5/4))
### 
###         """
###         a_cells = set(self.polygon_cells(label, a))
###         a_cell_centers = set(cell.center() for cell in a_cells)
###         b_cells = set(self.polygon_cells(label, b))
### 
###         from flatsurf.geometry.euclidean import OrientedSegment
###         ab = OrientedSegment(a, b)
### 
###         for b_cell in b_cells:
###             for opposite_center, boundary in b_cell.boundary().items():
###                 if opposite_center in a_cell_centers:
###                     intersection = ab.intersection(boundary)
###                     if intersection is None:
###                         continue
###                     if isinstance(intersection, OrientedSegment):
###                         raise NotImplementedError  # would need to extract the endpoint closest to c here.
### 
###                     return OrientedSegment(a, intersection), OrientedSegment(intersection, b)
### 
###         return None
### 
###     def boundaries(self):
###         r"""
###         Return the boundaries between Voronoi cells.
### 
###         The bondary segments are indexed by the two :class:`VoronoiPolygonCell`
###         defining the segment.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: center = S(0, S.polygon(0).centroid())
###             sage: S = S.insert_marked_points(center).codomain()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: boundaries = V.boundaries()
###             sage: len(boundaries)
###             24
### 
###             sage: key = frozenset([V.polygon_cell((0, 0), (0, 0)), V.polygon_cell((0, 0), (1, 0))])
###             sage: boundaries[key]
###             OrientedSegment((1/2, 1/2), (1/2, 0))
### 
###         """
###         boundaries = {}
### 
###         for cell in self.cells():
###             for polygon_cell in cell.polygon_cells():
###                 for opposite_center, boundary_segment in polygon_cell.boundary().items():
###                     if self.surface().polygon(polygon_cell.label()).get_point_position(opposite_center).is_outside():
###                         print("missing a segment of the boundary")
###                         continue
###                     boundaries[frozenset([polygon_cell, self.polygon_cell(polygon_cell.label(), opposite_center)])] = boundary_segment
### 
###         return boundaries
### 
###     @cached_method
###     def radius_of_convergence2(self, center):
###         r"""
###         Return the radius of convergence squared when developing a power series
###         at ``center``.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: V.radius_of_convergence2(S(0, S.polygon(0).centroid()))
###             1/2*a + 1
###             sage: V.radius_of_convergence2(next(iter(S.vertices())))
###             1
###             sage: V.radius_of_convergence2(S(0, (1/2, 1/2)))
###             1/2
###             sage: V.radius_of_convergence2(S(0, (1/2, 0)))
###             1/4
###             sage: V.radius_of_convergence2(S(0, (1/4, 0)))
###             1/16
### 
###         """
###         # TODO: This is useful more generally and should not be limited to Voronoi cells.
###         # TODO: Return an Euclidean distance.
### 
###         if all(vertex.angle() == 1 for vertex in self._surface.vertices()):
###             from sage.all import oo
###             return oo
### 
###         erase_marked_points = self._surface.erase_marked_points()
###         center = erase_marked_points(center)
### 
###         if not center.is_vertex():
###             insert_marked_points = center.surface().insert_marked_points(center)
###             center = insert_marked_points(center)
### 
###         surface = center.surface()
### 
###         for connection in surface.saddle_connections():
###             start = surface(*connection.start())
###             end = surface(*connection.end())
###             if start == center and end.angle() != 1:
###                 x, y = connection.holonomy()
###                 return x**2 + y**2
### 
###         assert False
### 
###     # TODO: Make comparable and hashable.
### 
### 
### class VoronoiCell:
###     r"""
###     A cell of a :class:`VoronoiDiagram`.
### 
###     EXAMPLES::
### W
###         sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###         sage: from flatsurf import translation_surfaces
###         sage: S = translation_surfaces.regular_octagon()
###         sage: V = VoronoiDiagram(S, S.vertices())
###         sage: cells = V.cells()
###         sage: list(cells)
###         [Voronoi cell at Vertex 0 of polygon 0]
### 
###     """
### 
###     def __init__(self, diagram, center):
###         self._parent = diagram
###         self._center = center
### 
###     def plot(self, graphical_surface=None):
###         r"""
###         Return a graphical representation of this cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.cell(S(0, 0))
###             sage: cell.plot()
###             Graphics object consisting of 34 graphics primitives
### 
###         """
###         plot_surface = graphical_surface is None
### 
###         if graphical_surface is None:
###             graphical_surface = self.surface().graphical_surface()
### 
###         plot = []
###         if plot_surface:
###             plot.append(graphical_surface.plot())
### 
###         for polygon_cell in self.polygon_cells():
###             plot.append(polygon_cell.plot(graphical_polygon=graphical_surface.graphical_polygon(polygon_cell.label())))
### 
###         return sum(plot)
### 
###     def surface(self):
###         r"""
###         Return this surface which this cell is a subset of.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.cell(S(0, 0))
###             sage: cell.surface() is S
###             True
### 
###         """
###         return self._parent.surface()
### 
###     def polygon_cells(self, label=None):
###         r"""
###         Return this cell broken into pieces that live entirely inside a
###         polygon.
### 
###         If ``label`` is specified, only the bits that are inside the polygon
###         with ``label`` are returned.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.cell(S(0, 0))
###             sage: cell.polygon_cells()
###             [Voronoi cell in polygon 0 at (0, a + 1),
###              Voronoi cell in polygon 0 at (1, a + 1),
###              Voronoi cell in polygon 0 at (0, 0),
###              Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a + 1),
###              Voronoi cell in polygon 0 at (1, 0),
###              Voronoi cell in polygon 0 at (-1/2*a, 1/2*a),
###              Voronoi cell in polygon 0 at (-1/2*a, 1/2*a + 1),
###              Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a)]
### 
###         """
###         cells = []
### 
###         for (lbl, coordinates) in self._center.representatives():
###             if label is not None and label != lbl:
###                 continue
### 
###             for c in self._parent.polygon_cells(lbl, coordinates):
###                 cells.append(c)
### 
###         # TODO: Why should we assert this here?
###         # assert cells
### 
###         return cells
### 
###     def contains_point(self, point):
###         r"""
###         Return whether this cell contains the ``point`` of the surface.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.cell(S(0, 0))
###             sage: point = S(0, (0, 0))
###             sage: cell.contains_point(point)
###             True
### 
###         """
###         label, coordinates = point.representative()
###         return any(cell.contains_point(coordinates) for cell in self.polygon_cells(label=label))
### 
###     def radius_of_convergence(self):
###         return self._parent.radius_of_convergence2(self._center)
### 
###     @cached_method
###     def radius(self):
###         # Returns the distance the point furthest from the center.
###         return max(cell.radius() for cell in self.polygon_cells())
### 
###     @cached_method
###     def furthest_point(self):
###         cell = max(self.polygon_cells(), key=lambda cell: cell.radius())
###         return self.surface()(cell.label(), cell.furthest_point())
### 
###     def __repr__(self):
###         return f"Voronoi cell at {self._center}"
### 
###     def corners(self):
###         corners = set()
###         for cell in self.polygon_cells():
###             polygon_corners = {}
###             for segment in cell.boundary().values():
###                 for p in [segment.start(), segment.end()]:
###                     p = self.surface()(cell.label(), p)
###                     if p not in polygon_corners:
###                         polygon_corners[p] = 1
###                     else:
###                         corners.add(p)
### 
###         return corners
### 
###     # TODO: Make comparable and hashable.
### 
### 
### class VoronoiDiagram_Polygon:
###     r"""
###     The part of a :class:`VoronoiDiagram" inside the polygon with ``label``.
### 
###     EXAMPLES::
### 
###         sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###         sage: from flatsurf import translation_surfaces
###         sage: S = translation_surfaces.regular_octagon()
###         sage: V = VoronoiDiagram(S, S.vertices())
### 
###         sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###         sage: VoronoiDiagram_Polygon(V, 0)
###         Voronoi diagram in polygon 0
### 
###     """
### 
###     def __init__(self, parent, label, weight=None):
###         self._parent = parent
###         self._label = label
###         self._weight = weight or "classical"
### 
###     @cached_method
###     def _cells(self):
###         r"""
###         Return the cells that make up this Voronoi diagram indexed by their
###         centers.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VD = VoronoiDiagram_Polygon(V, 0)
###             sage: VD._cells()
###             {(-1/2*a, 1/2*a): Voronoi cell in polygon 0 at (-1/2*a, 1/2*a),
###              (-1/2*a, 1/2*a + 1): Voronoi cell in polygon 0 at (-1/2*a, 1/2*a + 1),
###              (0, 0): Voronoi cell in polygon 0 at (0, 0),
###              (0, a + 1): Voronoi cell in polygon 0 at (0, a + 1),
###              (1, 0): Voronoi cell in polygon 0 at (1, 0),
###              (1, a + 1): Voronoi cell in polygon 0 at (1, a + 1),
###              (1/2*a + 1, 1/2*a): Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a),
###              (1/2*a + 1, 1/2*a + 1): Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a + 1)}
### 
###         """
###         return {center: VoronoiPolygonCell(self, center) for center in self.centers()}
### 
###     def label(self):
###         r"""
###         Return the label of the polygon that this diagram breaks into Voronoi
###         cells.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VD = VoronoiDiagram_Polygon(V, 0)
###             sage: VD.label()
###             0
### 
###         """
###         return self._label
### 
###     @cached_method
###     def centers(self):
###         r"""
###         Return the coordinates of the centers of Voronoi cells in this polygon.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VD = VoronoiDiagram_Polygon(V, 0)
###             sage: VD.centers()
###             [(0, a + 1),
###              (1, a + 1),
###              (0, 0),
###              (1/2*a + 1, 1/2*a + 1),
###              (1, 0),
###              (-1/2*a, 1/2*a),
###              (-1/2*a, 1/2*a + 1),
###              (1/2*a + 1, 1/2*a)]
### 
###         """
###         return [center_coordinates for center in self._parent._centers for (label, center_coordinates) in center.representatives() if label == self._label]
### 
###     def surface(self):
###         r"""
###         Return the surface containing this polygon.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VD = VoronoiDiagram_Polygon(V, 0)
###             sage: VD.surface() is S
###             True
### 
###         """
###         return self._parent.surface()
### 
###     def polygon(self):
###         r"""
###         Return the polygon that this diagram breaks up into cells.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VD = VoronoiDiagram_Polygon(V, 0)
###             sage: VD.polygon() == S.polygon(0)
###             True
### 
###         """
###         return self.surface().polygon(self._label)
### 
###     def polygon_cells(self, coordinates):
###         r"""
###         Return the :class:`VoronoiPolygonCell`s that contain the point at
###         ``coordinates``.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VP = VoronoiDiagram_Polygon(V, 0)
### 
###             sage: VP.polygon_cells((0, 0))
###             [Voronoi cell in polygon 0 at (0, 0)]
###             sage: VP.polygon_cells((1, 1))
###             [Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a)]
### 
###         """
###         cells = [cell for cell in self._cells().values() if cell.contains_point(coordinates)]
###         assert cells
###         return cells
### 
###     def _half_space(self, center, opposite_center):
###         r"""
###         Return the half space containing ``center`` but not containing
###         ``opposite_center`` that sits half-way (up to weighting) between these
###         two points.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VP = VoronoiDiagram_Polygon(V, 0)
###             sage: VP._half_space(vector((0, 0)), vector((1, 0)))
###             {-x ≥ -1/2}
### 
###         """
###         if self._weight == "classical":
###             return self._half_space_weighted(center, opposite_center, 1)
###         elif self._weight == "radius_of_convergence":
###             return self._half_space_radius_of_convergence(center, opposite_center)
###         else:
###             raise NotImplementedError("unsupported weight for Voronoi cells")
### 
###     def _half_space_weighted(self, center, opposite_center, weight=1):
###         r"""
###         Return the half space containing ``center`` but not containing
###         ``opposite_center``.
### 
###         This produces a weighted version of that half space, i.e., the boundary
###         point on the segment between the two centers is shifted according to
###         the weight towards the ``opposite_center``.
### 
###         TODO: Explain how just using half spaces might leave some empty space.
### 
###         Note that this is not a natural generalization of Voronoi cells.
###         Normally, one would all points on the boundary to have a distance that
###         is weighted in that way. However, this is a bit more complicated to
###         implement as you get a more complicated curve and then also harder to
###         integrate along later. So we opted for not implementing that.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VP = VoronoiDiagram_Polygon(V, 0)
###             sage: VP._half_space_weighted(vector((0, 0)), vector((1, 0)), 1)
###             {-x ≥ -1/2}
###             sage: VP._half_space_weighted(vector((0, 0)), vector((1, 0)), 2)
###             {-x ≥ -2/3}
###             sage: VP._half_space_weighted(vector((0, 0)), vector((1, 0)), 1/2)
###             {-x ≥ -1/3}
### 
###         """
###         if weight <= 0:
###             raise ValueError("weight must be positive")
### 
###         a, b = center - opposite_center
###         midpoint = (weight * opposite_center + center) / (weight + 1)
###         c = a * midpoint[0] + b * midpoint[1]
### 
###         from flatsurf.geometry.euclidean import HalfSpace
###         return HalfSpace(-c, a, b)
### 
###     def _half_space_radius_of_convergence(self, center, opposite_center):
###         r"""
###         Return the half space containing ``center`` but not containing
###         ``opposite_center``.
### 
###         The point on the boundary of that half space and the segment connecting
###         the two centers is shifted relative to the respective radius of
###         convergence at these centers.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VP = VoronoiDiagram_Polygon(V, 0)
###             sage: VP._half_space_radius_of_convergence(vector((0, 0)), vector((1, 0)))
###             {-x ≥ -1/2}
###             sage: VP._half_space_radius_of_convergence(vector((0, 0)), S.polygon(0).vertices()[4])
###             {-x - (-a - 1) * y ≥ -a - 2}
### 
###         """
###         return self._half_space_weighted(center, opposite_center, self._half_space_radius_of_convergence_weight(center, opposite_center))
### 
###     def _half_space_radius_of_convergence_weight(self, center, opposite_center):
###         r"""
###         Return the quotient of the radii of convergence at ``center`` and
###         ``opposite_center`` (when restricted to their polygon.)
###         """
###         surface = self._parent.surface()
###         weight = self._parent.radius_of_convergence2(surface(self._label, center))
###         opposite_weight = self._parent.radius_of_convergence2(surface(self._label, opposite_center))
### 
###         from sage.all import oo
###         if weight == oo:
###             assert opposite_weight == oo
###             return 1
### 
###         relative_weight = weight / opposite_weight
###         try:
###             return relative_weight.parent()(relative_weight.sqrt())
###         except Exception:
###             # TODO: This blows up coefficients too much. We added some rounding but that's also a hack.
###             # When the weight does not exist in the base ring we take an
###             # approximation (with possibly huge coefficients.)
###             if relative_weight > 1:
###                 # Make rounding errors symmetric so that two neighboring half
###                 # space are actually the negative of each other.
###                 return 1 / self._half_space_radius_of_convergence_weight(opposite_center, center)
### 
###             from math import sqrt
###             return relative_weight.parent()(round(sqrt(float(relative_weight)), 4))
### 
###     def half_spaces(self, center):
###         r"""
###         Return the half spaces that define the Voronoi cell centered at
###         ``center`` in this polygon, indexed by the other center that is
###         defining the half space.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VP = VoronoiDiagram_Polygon(V, 0)
###             sage: VP.half_spaces((0, 0))
###             {(-1/2*a, 1/2*a): {(1/4*a + 1/2) * x - (-1/4*a - 1/2) * y ≥ -1/4*a - 1/4},
###              (1, 0): {(-1/2*a - 1/2) * x ≥ -1/4*a - 1/4}}
### 
###         """
###         return {opposite_center: segment.left_half_space() for (opposite_center, segment) in self.boundary(center).items()}
### 
###     @cached_method
###     def boundary(self, center):
###         r"""
###         Return the boundary segment that define the Voronoi cell centered at
###         ``center`` in this polygon, indexed by the other center that is
###         defining the segment.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VP = VoronoiDiagram_Polygon(V, 0)
###             sage: VP.boundary((0, 0))
###             {(-1/2*a, 1/2*a): OrientedSegment((1/2, 1/2*a + 1/2), (-1/4*a, 1/4*a)),
###              (1, 0): OrientedSegment((1/2, 0), (1/2, 1/2*a + 1/2))}
### 
###         """
###         from sage.all import vector
###         center = vector(center)
### 
###         if center not in self.centers():
###             raise ValueError("center must be a center of a Voronoi cell")
### 
###         voronoi_half_spaces = {opposite_center: self._half_space(center, opposite_center) for opposite_center in self.centers() if opposite_center != center}
### 
###         if self._weight == "classical":
###             # TODO: Hack in the next polygon if it affects this cell.
###             vertex = self.polygon().get_point_position(center).get_vertex()
###             opposite_label, opposite_edge = self.surface().opposite_edge(self.label(), vertex)
###             opposite_center = center + self.surface().polygon(opposite_label).edge(opposite_edge + 1)
###             opposite_center.set_immutable()
###             voronoi_half_spaces[opposite_center] = self._half_space(center, opposite_center)
### 
###             # TODO: Hack in the previous polygon if it affects this cell.
###             opposite_label, opposite_edge = self.surface().opposite_edge(self.label(), vertex - 1)
###             opposite_center = center - self.surface().polygon(opposite_label).edge(opposite_edge - 1)
###             opposite_center.set_immutable()
###             voronoi_half_spaces[opposite_center] = self._half_space(center, opposite_center)
### 
###         # The half spaces whose intersection is the entire polygon.
###         from flatsurf.geometry.euclidean import HalfSpace
###         polygon = self.polygon()
###         polygon_half_spaces = [HalfSpace(vertex[0] * edge[1] - vertex[1] * edge[0], -edge[1], edge[0]) for (vertex, edge) in zip(polygon.vertices(), polygon.edges())]
### 
###         # Each segment defining this Voronoi cell is on the boundary of one of
###         # the half spaces defining this Voronoi cell. Namely, if the half space
###         # is not trivial in the intersection, then it contributes a segment to
###         # the boundary of the cell.
### 
###         segments = HalfSpace.compact_intersection(*voronoi_half_spaces.values(), *polygon_half_spaces)
### 
###         assert segments, "Voronoi cell is empty"
### 
###         # Filter out segments that are mearly edges of the polygon; they
###         # are an artifact of how we computed the segments here.
###         # TODO: This is very expensive in comparison to a simple:
###         #   segments = [segment for segment in segments if center not in segment]
###         # But that misses some degenerate cases. But do these cases actually
###         # make any sense?
### 
###         from flatsurf.geometry.euclidean import OrientedSegment
###         polygon_edges = [OrientedSegment(polygon.vertex(i), polygon.vertex(i + 1)) for i in range(len(polygon.vertices()))]
###         voronoi_segments = [segment for segment in segments if not any(segment.is_subset(edge) for edge in polygon_edges)]
### 
###         assert voronoi_segments
### 
###         # Associate with each segment which other center produced it.
###         boundary = {}
###         for segment in voronoi_segments:
###             opposite_center = [c for (c, half_space) in voronoi_half_spaces.items() if half_space == segment.left_half_space()]
###             assert len(opposite_center) == 1, "segment must be induced by exactly one other center of a Voronoi cell"
###             opposite_center = next(iter(opposite_center))
###             assert opposite_center not in boundary
### 
###             boundary[opposite_center] = segment
### 
###         return boundary
### 
###     def __repr__(self):
###         return f"Voronoi diagram in polygon {self.label()}"
### 
###     def __eq__(self, other):
###         r"""
###         Return whether this Voronoi diagram is indistinguishable from ``other``.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
###             sage: VoronoiDiagram_Polygon(V, 0) == VoronoiDiagram_Polygon(V, 0)
###             True
### 
###         """
###         if not isinstance(other, VoronoiDiagram_Polygon):
###             return False
### 
###         return self._parent == other._parent and self._label == other._label and self._weight == other._weight
### 
###     def __ne__(self, other):
###         return not (self == other)
### 
###     def __hash__(self):
###         return hash((self._label, self._weight))
### 
### 
### class VoronoiPolygonCell:
###     r"""
###     The part of a Voronoi cell that lives entirely within a single polygon that
###     makes up a surface.
### 
###     EXAMPLES::
### 
###         sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###         sage: from flatsurf import translation_surfaces
###         sage: S = translation_surfaces.regular_octagon()
###         sage: V = VoronoiDiagram(S, S.vertices())
###         sage: V.polygon_cell(0, (0, 0))
###         Voronoi cell in polygon 0 at (0, 0)
### 
###     """
### 
###     def __init__(self, parent, center):
###         self._parent = parent
###         self._center = center
### 
###     @cached_method
###     def half_spaces(self):
###         r"""
###         Return the half spaces that delimit this cell inside its polygon,
###         indexed by the two centers defining the half space.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.half_spaces()
###             {(-1/2*a, 1/2*a): {(1/4*a + 1/2) * x - (-1/4*a - 1/2) * y ≥ -1/4*a - 1/4},
###              (1, 0): {(-1/2*a - 1/2) * x ≥ -1/4*a - 1/4}}
### 
###         """
###         return self._parent.half_spaces(self._center)
### 
###     def boundary(self):
###         r"""
###         Return the segments that delimit this cell inside its polygon indexed
###         by the other centers defining the half space.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.boundary()
###             {(-1/2*a, 1/2*a): OrientedSegment((1/2, 1/2*a + 1/2), (-1/4*a, 1/4*a)),
###              (1, 0): OrientedSegment((1/2, 0), (1/2, 1/2*a + 1/2))}
### 
###         """
###         return self._parent.boundary(self._center)
### 
###     def polygon(self):
###         r"""
###         Return the polygon which contains this cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.polygon()
###             Polygon(vertices=[(0, 0), (1, 0), (1/2*a + 1, 1/2*a), (1/2*a + 1, 1/2*a + 1), (1, a + 1), (0, a + 1), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a)])
### 
###         """
###         return self._parent.polygon()
### 
###     def label(self):
###         r"""
###         Return the label of the polygon this cell is a part of.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.label()
###             0
### 
###         """
###         return self._parent.label()
### 
###     def center(self):
###         r"""
###         Return the coordinates of the center of this Voronoi cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.center()
###             (0, 0)
### 
###         """
###         return self._center
### 
###     def surface(self):
###         r"""
###         Return the surface containing the :meth:`polygon`.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.surface() is S
###             True
### 
###         """
###         return self._parent.surface()
### 
###     def contains_point(self, point):
###         r"""
###         Return whether this cell contains the ``point``.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
### 
###             sage: cell.contains_point((0, 0))
###             True
###             sage: cell.contains_point((1/2, 1/2))
###             True
###             sage: cell.contains_point((1, 1))
###             False
### 
###         """
###         if self.polygon().get_point_position(point).is_outside():
###             return False
### 
###         return all(half_space.contains_point(point) for half_space in self.half_spaces().values())
### 
###     def contains_segment(self, segment):
###         r"""
###         Return whether the ``segment`` is a subset of this cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
### 
###             sage: from flatsurf.geometry.euclidean import OrientedSegment
###             sage: cell.contains_segment(OrientedSegment((0, 0), (1, 0)))
###             False
### 
###         """
###         if isinstance(segment, tuple):
###             from flatsurf.geometry.euclidean import OrientedSegment
###             segment = OrientedSegment(*segment)
### 
###         return self.contains_point(segment.start()) and self.contains_point(segment.end())
### 
###     def plot(self, graphical_polygon=None):
###         r"""
###         Return a graphical representation of this cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
###             sage: cell.plot()
###             Graphics object consisting of 3 graphics primitives
### 
###         """
###         plot_polygon = graphical_polygon is None
### 
###         if graphical_polygon is None:
###             graphical_polygon = self.surface().graphical_surface().graphical_polygon(self.label())
### 
###         shift = graphical_polygon.transformed_vertex(0) - self.polygon().vertex(0)
### 
###         plot = sum((segment.translate(shift)).plot(point=False) for segment in self.boundary().values())
### 
###         if plot_polygon:
###             plot = graphical_polygon.plot_polygon() + plot
### 
###         return plot
### 
###     @cached_method
###     def radius(self):
###         Δ = self._center - self.furthest_point()
###         return Δ.dot_product(Δ)
### 
###     @cached_method
###     def furthest_point(self):
###         vertices = []
###         for segment in self.boundary().values():
###             vertices.append(segment.start())
###             vertices.append(segment.end())
### 
###         def norm2(x):
###             return x.dot_product(x)
### 
###         return max(vertices, key=lambda v: norm2(v - self._center))
### 
###     def __repr__(self):
###         r"""
###         Return a printable representation of this Voronoi cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.polygon_cell(0, (0, 0))
###             Voronoi cell in polygon 0 at (0, 0)
### 
###         """
###         return f"Voronoi cell in polygon {self.label()} at {self._center}"
### 
###     def __eq__(self, other):
###         r"""
###         Return whether this cell is indistinguishable from ``other``.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: V.polygon_cell(0, (0, 0)) == V.polygon_cell(0, (0, 1/3))
###             True
### 
###         """
###         if not isinstance(other, VoronoiPolygonCell):
###             return False
###         return self._center == other._center and self._parent == other._parent
### 
###     def __ne__(self, other):
###         return not (self == other)
### 
###     def __hash__(self):
###         return hash((self._parent, self._center))
### 
###     def split_segment_uniform_root_branch(self, segment):
###         r"""
###         Return the ``segment`` split into smaller segments such that these
###         segments do not cross the horizontal line left of the center of the
###         cell (if that center is an actual singularity and not just a marked
###         point.)
### 
###         On such a shorter segment, we can then develop an n-th root
###         consistently where n-1 is the order of the singularity.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (1, 0))
### 
###             sage: from flatsurf.geometry.euclidean import OrientedSegment
###             sage: cell.split_segment_uniform_root_branch(OrientedSegment((0, -1), (0, 1)))
###             [OrientedSegment((0, -1), (0, 0)), OrientedSegment((0, 0), (0, 1))]
### 
###         """
###         d = self.surface()(self.label(), self.center()).angle()
### 
###         assert d >= 1
###         if d == 1:
###             return [segment]
### 
###         from flatsurf.geometry.euclidean import OrientedSegment
### 
###         if segment.contains_point(self.center()) and segment.start() != self.center() and segment.end() != self.center():
###             return [OrientedSegment(segment.start(), self.center()), OrientedSegment(self.center(), segment.end())]
### 
###         from flatsurf.geometry.euclidean import Ray
###         ray = Ray(self.center(), (-1, 0))
### 
###         if ray.contains_point(segment.start()) or ray.contains_point(segment.end()):
###             return [segment]
### 
###         intersection = ray.intersection(segment)
### 
###         if intersection is None:
###             return [segment]
### 
###         return [OrientedSegment(segment.start(), intersection), OrientedSegment(intersection, segment.end())]
### 
###     def root_branch(self, segment):
###         r"""
###         Return which branch can be taken consistently along the ``segment``
###         when developing an n-th root at the center of this Voronoi cell.
### 
###         EXAMPLES::
### 
###             sage: from flatsurf.geometry.voronoi import VoronoiDiagram
###             sage: from flatsurf import translation_surfaces
###             sage: S = translation_surfaces.regular_octagon()
###             sage: V = VoronoiDiagram(S, S.vertices())
###             sage: cell = V.polygon_cell(0, (0, 0))
### 
###             sage: from flatsurf.geometry.euclidean import OrientedSegment
###             sage: cell.root_branch(OrientedSegment((0, -1), (0, 1)))
###             Traceback (most recent call last):
###             ...
###             ValueError: segment does not permit a consistent choice of root
### 
###             sage: cell.root_branch(OrientedSegment((0, 0), (0, 1/2)))
###             0
### 
###             sage: cell = V.polygon_cell(0, (1, 0))
###             sage: cell.root_branch(OrientedSegment((0, 0), (0, 1/2)))
###             1
### 
###             sage: a = S.base_ring().gen()
###             sage: cell = V.polygon_cell(0, (1 + a/2, a/2))
###             sage: cell.root_branch(OrientedSegment((1, 1/2 + a/2), (1 + a/2, 1/2 + a/2)))
###             2
### 
###         """
###         if self.split_segment_uniform_root_branch(segment) != [segment]:
###             raise ValueError("segment does not permit a consistent choice of root")
### 
###         S = self.surface()
### 
###         center = S(self.label(), self.center())
### 
###         angle = center.angle()
### 
###         assert angle >= 1
###         if angle == 1:
###             return 0
### 
###         # Choose a horizontal ray to the right, that defines where the
###         # principal root is being used. We use the "smallest" vertex in the
###         # "smallest" polygon containing such a ray.
###         from flatsurf.geometry.euclidean import ccw
###         primitive_label, primitive_vertex = min((label, vertex) for (label, _) in center.representatives() for vertex in range(len(S.polygon(label).vertices()))
###             if S(label, vertex) == center and
###                ccw((1, 0), S.polygon(label).edge(vertex)) <= 0 and
###                ccw((1, 0), -S.polygon(label).edge(vertex - 1)) >= 0)
### 
###         # Walk around the vertex to determine the branch of the root for the
###         # (midpoint of) the segment.
###         point = segment.midpoint()
### 
###         branch = 0
###         label = primitive_label
###         vertex = primitive_vertex
### 
###         while True:
###             polygon = S.polygon(label)
###             if label == self.label() and polygon.vertex(vertex) == self.center():
###                 low = ccw((-1, 0), polygon.edge(vertex)) <= 0 and ccw((-1, 0), point - polygon.vertex(vertex)) > 0
###                 if low:
###                     return (branch + 1) % angle
###                 return branch
### 
###             if ccw((-1, 0), polygon.edge(vertex)) <= 0 and ccw((-1, 0), -polygon.edge(vertex - 1)) > 0:
###                 branch += 1
###                 branch %= angle
### 
###             label, vertex = S.opposite_edge(label, (vertex - 1) % len(polygon.vertices()))
