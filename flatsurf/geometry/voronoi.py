# TODO: Drop all use of coordinates. Use EuclideanPolygonPoint instead of (label, coordinates)
# TODO: SageMath already has a VoronoiDiagram so maybe not use that name. Explain that we cannot use the builtin version since it has no weights.

from sage.misc.cachefunc import cached_method


class VoronoiDiagram:
    r"""
    ALGORITHM:

    Internally, a Voronoi diagram is composed by :class:`VoronoiCell`s. Each
    cell is the union of :class:`VoronoiPolygonCell` which represents the
    intersection of a Voronoi cell with an Euclidean polygon making up the
    surface. To compute these cells on the level of polygons, we determine the
    segments that bound the Voronoi cell where each such segment is half way (up
    to optional weights) between the centers of two Voronoi cells.

    TODO: Reducing the computation to the level of polygons is not correct in
    general. When centers are close to edges, the polygon on the other side of
    the edge does not see that center. Can this be fixed? Or should we just say
    that this is an approximation of the Voronoi diagram? And is it correct
    without weighing and with centers at the vertices of a Delaunay cell
    decomposition?

    TODO: Explain how these cells are not minimal. The same point that exists
    twice in a polygon produces a boundary of the Voronoi cell with itself. Add
    a method that gets rid of these boundaries (and throw a NotImplementedError
    since it is unclear how to represent cells that contain an entire polygon?)

    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram  # random output due to deprecation warnings
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: center = S(0, S.polygon(0).centroid())
        sage: V = VoronoiDiagram(S, S.vertices().union([center]))
        sage: V.plot()
        Graphics object consisting of 41 graphics primitives
        sage: V = VoronoiDiagram(S, S.vertices().union([center]), weight="radius_of_convergence")
        sage: V.plot()
        Graphics object consisting of 41 graphics primitives

    The same Voronoi diagram but starting from a more complicated description
    of the octagon::

        sage: from flatsurf import Polygon, translation_surfaces, polygons
        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: S = translation_surfaces.regular_octagon()
        sage: S = S.subdivide().codomain()
        sage: V = VoronoiDiagram(S, S.vertices())
        sage: V.plot()
        Graphics object consisting of 73 graphics primitives

        sage: V = VoronoiDiagram(S, S.vertices(), weight="radius_of_convergence")
        sage: V.plot()
        Graphics object consisting of 73 graphics primitives

    """

    def __init__(self, surface, points, weight=None):
        if surface.is_mutable():
            raise ValueError("surface must be immutable")

        if weight is None:
            weight = "classical"

        self._surface = surface
        self._centers = set(points)
        self._weight = weight

        if not surface.vertices().issubset(self._centers):
            # This probably essentially works. However, we cannot represent
            # cells that contain an entire polygon yet, so this is not
            # implemented in general.
            raise NotImplementedError("can only compute Voronoi diagrams when all vertices are centers")

        self._cells = {center: VoronoiCell(self, center) for center in self._centers}

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

    def plot(self, graphical_surface=None):
        r"""
        Return a graphical representation of this Voronoi diagram.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: V.plot()
            Graphics object consisting of 41 graphics primitives

        The underlying surface is not plotted automatically when it is provided
        as a keyword argument::

            sage: V.plot(graphical_surface=S.graphical_surface())
            Graphics object consisting of 32 graphics primitives

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

    def cells(self, point=None):
        r"""
        Return the Voronoi cells that contain ``point``.

        If no ``point`` is given, return all cells.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: cells = V.cells(S(0, (1/2, 0)))
            sage: list(cells)
            [Voronoi cell at Vertex 0 of polygon 0]

        """
        for cell in self._cells.values():
            if point is None or cell.contains_point(point):
                yield cell

    def polygon_cell(self, label, coordinates):
        r"""
        Return a Voronoi cell that contains the point ``coordinates`` of the
        polygon with ``label``.

        This is to :meth:`polygon_cells` what :meth:`cell` is to :meth:`cells`.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: V.polygon_cell(0, S.polygon(0).centroid())
            Voronoi cell in polygon 0 at (1/2, 1/2*a + 1/2)

        """
        return next(iter(self.polygon_cells(label, coordinates)))

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
        diagram = VoronoiDiagram_Polygon(self, label, self._weight)
        return diagram.polygon_cells(coordinates)

    def split_segment(self, label, a, b):
        r"""
        Return the segment [a, b] split into shorter segments that each lie in
        a single Voronoi cell.

        The segments are returned as a dict mapping the Voronoi cell to the
        shorter segments.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: V.split_segment(0, (0, 0), (1, 1))
            {Voronoi cell in polygon 0 at (0, 0): OrientedSegment((0, 0), (1/2, 1/2)),
             Voronoi cell in polygon 0 at (1/2, 1/2*a + 1/2): OrientedSegment((1/2, 1/2), (1, 1))}

        When there are multiple ways to split the segment, the choice is made
        randomly. Here, the first segment, is on the boundary of two cells::

            sage: V.split_segment(0, (1/2, 0), (1/2, 1))
            {Voronoi cell in polygon 0 at (...): OrientedSegment((1/2, 0), (1/2, 1/2)),
             Voronoi cell in polygon 0 at (1/2, 1/2*a + 1/2): OrientedSegment((1/2, 1/2), (1/2, 1))}

        """
        segments = {}

        while a != b:
            a_cells = set(self.polygon_cells(label, a))
            b_cells = set(self.polygon_cells(label, b))

            shared_cells = set(a_cells).intersection(b_cells)

            if shared_cells:
                # Both endpoints are in the same Voronoi cell.
                cell = next(iter(shared_cells))
                assert cell not in segments

                from flatsurf.geometry.euclidean import OrientedSegment
                segments[cell] = OrientedSegment(a, b)

                break

            # Bisect the segment [a, b] until we find a point c such that [a,
            # c] crosses the boundary between two neighboring Voronoi cells.
            c = b
            while True:
                split = self._split_segment_at_boundary_point(label, a, c)

                if split is None:
                    c = (a + c) / 2
                    continue

                segment, _ = split

                cell = next(iter(cell for cell in a_cells if cell.contains_segment(segment)))

                assert cell not in segments
                segments[cell] = segment
                a = segment.end()
                break

        return segments

    def _split_segment_at_boundary_point(self, label, a, b):
        r"""
        Return the point that lies on the segment [a, b] and on the boundary
        between the Voronoi cells containing a and b respectively.

        Return ``None`` if no such point exists.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: center = S(0, S.polygon(0).centroid())
            sage: V = VoronoiDiagram(S, S.vertices().union([center]))
            sage: V._split_segment_at_boundary_point(0, (0, 0), (1, 0))
            (OrientedSegment((0, 0), (1/2, 0)), OrientedSegment((1/2, 0), (1, 0)))
            sage: V._split_segment_at_boundary_point(0, (0, 0), (1, 1))
            (OrientedSegment((0, 0), (1/2, 1/2)), OrientedSegment((1/2, 1/2), (1, 1)))
            sage: V._split_segment_at_boundary_point(0, (0, 0), (5/4, 5/4))

        """
        a_cells = set(self.polygon_cells(label, a))
        a_cell_centers = set(cell.center() for cell in a_cells)
        b_cells = set(self.polygon_cells(label, b))

        from flatsurf.geometry.euclidean import OrientedSegment
        ab = OrientedSegment(a, b)

        for b_cell in b_cells:
            for opposite_center, boundary in b_cell.boundary().items():
                if opposite_center in a_cell_centers:
                    intersection = ab.intersection(boundary)
                    if intersection is None:
                        continue
                    if isinstance(intersection, OrientedSegment):
                        raise NotImplementedError  # would need to extract the endpoint closest to c here.

                    return OrientedSegment(a, intersection), OrientedSegment(intersection, b)

        return None


class VoronoiCell:
    r"""
    A cell of a :class:`VoronoiDiagram`.

    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: V = VoronoiDiagram(S, S.vertices())
        sage: cells = V.cells()
        sage: list(cells)
        [Voronoi cell at Vertex 0 of polygon 0]

    """

    def __init__(self, diagram, center):
        self._parent = diagram
        self._center = center

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
            graphical_surface = self.surface().graphical_surface()

        plot = []
        if plot_surface:
            plot.append(graphical_surface.plot())

        for polygon_cell in self.polygon_cells():
            plot.append(polygon_cell.plot(graphical_polygon=graphical_surface.graphical_polygon(polygon_cell.label())))

        return sum(plot)

    def surface(self):
        r"""
        Return this surface which this cell is a subset of.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.cell(S(0, 0))
            sage: cell.surface() is S
            True

        """
        return self._parent.surface()

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

            for c in self._parent.polygon_cells(lbl, coordinates):
                cells.append(c)

        return cells

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

    def __repr__(self):
        return f"Voronoi cell at {self._center}"


class VoronoiDiagram_Polygon:
    r"""
    The part of a :class:`VoronoiDiagram" inside the polygon with ``label``.

    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: V = VoronoiDiagram(S, S.vertices())

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
        sage: VoronoiDiagram_Polygon(V, 0)
        Voronoi diagram in polygon 0

    """

    def __init__(self, parent, label, weight=None):
        self._parent = parent
        self._label = label
        self._weight = weight or "classical"

    @cached_method
    def _cells(self):
        r"""
        Return the cells that make up this Voronoi diagram indexed by their
        centers.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VD = VoronoiDiagram_Polygon(V, 0)
            sage: VD._cells()
            {(-1/2*a, 1/2*a): Voronoi cell in polygon 0 at (-1/2*a, 1/2*a),
             (-1/2*a, 1/2*a + 1): Voronoi cell in polygon 0 at (-1/2*a, 1/2*a + 1),
             (0, 0): Voronoi cell in polygon 0 at (0, 0),
             (0, a + 1): Voronoi cell in polygon 0 at (0, a + 1),
             (1, 0): Voronoi cell in polygon 0 at (1, 0),
             (1, a + 1): Voronoi cell in polygon 0 at (1, a + 1),
             (1/2*a + 1, 1/2*a): Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a),
             (1/2*a + 1, 1/2*a + 1): Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a + 1)}

        """
        return {center: VoronoiPolygonCell(self, center) for center in self.centers()}

    def label(self):
        r"""
        Return the label of the polygon that this diagram breaks into Voronoi
        cells.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VD = VoronoiDiagram_Polygon(V, 0)
            sage: VD.label()
            0

        """
        return self._label

    @cached_method
    def centers(self):
        r"""
        Return the coordinates of the centers of Voronoi cells in this polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VD = VoronoiDiagram_Polygon(V, 0)
            sage: VD.centers()
            [(0, a + 1),
             (1, a + 1),
             (0, 0),
             (1/2*a + 1, 1/2*a + 1),
             (1, 0),
             (-1/2*a, 1/2*a),
             (-1/2*a, 1/2*a + 1),
             (1/2*a + 1, 1/2*a)]

        """
        return [center_coordinates for center in self._parent._centers for (label, center_coordinates) in center.representatives() if label == self._label]

    def surface(self):
        r"""
        Return the surface containing this polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VD = VoronoiDiagram_Polygon(V, 0)
            sage: VD.surface() is S
            True

        """
        return self._parent.surface()

    def polygon(self):
        r"""
        Return the polygon that this diagram breaks up into cells.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VD = VoronoiDiagram_Polygon(V, 0)
            sage: VD.polygon() == S.polygon(0)
            True

        """
        return self.surface().polygon(self._label)

    def polygon_cells(self, coordinates):
        r"""
        Return the :class:`VoronoiPolygonCell`s that contain the point at
        ``coordinates``.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VP = VoronoiDiagram_Polygon(V, 0)

            sage: VP.polygon_cells((0, 0))
            [Voronoi cell in polygon 0 at (0, 0)]
            sage: VP.polygon_cells((1, 1))
            [Voronoi cell in polygon 0 at (1/2*a + 1, 1/2*a)]

        """
        return [cell for cell in self._cells().values() if cell.contains_point(coordinates)]

    def _half_space(self, center, opposite_center):
        r"""
        Return the half space containing ``center`` but not containing
        ``opposite_center`` that sits half-way (up to weighting) between these
        two points.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VP = VoronoiDiagram_Polygon(V, 0)
            sage: VP._half_space(vector((0, 0)), vector((1, 0)))
            {-x ≥ -1/2}

        """
        if self._weight == "classical":
            return self._half_space_weighted(center, opposite_center, 1)
        elif self._weight == "radius_of_convergence":
            return self._half_space_radius_of_convergence(center, opposite_center)
        else:
            raise NotImplementedError("unsupported weight for Voronoi cells")

    def _half_space_weighted(self, center, opposite_center, weight=1):
        r"""
        Return the half space containing ``center`` but not containing
        ``opposite_center``.

        This produces a weighted version of that half space, i.e., the boundary
        point on the segment between the two centers is shifted according to
        the weight towards the ``opposite_center``.

        Note that this is not a natural generalization of Voronoi cells.
        Normally, one would all points on the boundary to have a distance that
        is weighted in that way. However, this is a bit more complicated to
        implement as you get a more complicated curve and then also harder to
        integrate along later. So we opted for not implementing that.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VP = VoronoiDiagram_Polygon(V, 0)
            sage: VP._half_space_weighted(vector((0, 0)), vector((1, 0)), 1)
            {-x ≥ -1/2}
            sage: VP._half_space_weighted(vector((0, 0)), vector((1, 0)), 2)
            {-x ≥ -2/3}
            sage: VP._half_space_weighted(vector((0, 0)), vector((1, 0)), 1/2)
            {-x ≥ -1/3}

        """
        if weight <= 0:
            raise ValueError("weight must be positive")

        a, b = center - opposite_center
        midpoint = (weight * opposite_center + center) / (weight + 1)
        c = a * midpoint[0] + b * midpoint[1]

        from flatsurf.geometry.euclidean import HalfSpace
        return HalfSpace(-c, a, b)

    def _half_space_radius_of_convergence(self, center, opposite_center):
        r"""
        Return the half space containing ``center`` but not containing
        ``opposite_center``.

        The point on the boundary of that half space and the segment connecting
        the two centers is shifted relative to the respective radius of
        convergence at these centers.

        Note that we do not really compute the radius of convergence but just
        pretend that it is the minimum distance to any other vertex in the
        polygon (even if it is not a singularity of the surface.) Also, there
        might be closer singularities in other nearby polygons that we ignore
        for this computation.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VP = VoronoiDiagram_Polygon(V, 0)
            sage: VP._half_space_radius_of_convergence(vector((0, 0)), vector((1, 0)))
            {-x ≥ -1/2}
            sage: VP._half_space_radius_of_convergence(vector((0, 0)), S.polygon(0).vertices()[4])
            {-x - (-a - 1) * y ≥ -a - 2}

        """
        weight = min((v[0] - center[0])**2 + (v[1] - center[1])**2 for v in self.polygon().vertices() if v != center)
        opposite_weight = min((v[0] - opposite_center[0])**2 + (v[1] - opposite_center[1])**2 for v in self.polygon().vertices() if v != opposite_center)

        relative_weight = weight / opposite_weight
        try:
            relative_weight = relative_weight.parent()(relative_weight.sqrt())
        except Exception:
            # When the weight does not exist in the base ring we take an
            # approximation (with possibly huge coefficients.)
            relative_weight = relative_weight.parent()(float(relative_weight.sqrt()))

        return self._half_space_weighted(center, opposite_center, relative_weight)

    def half_spaces(self, center):
        r"""
        Return the half spaces that define the Voronoi cell centered at
        ``center`` in this polygon, indexed by the other center that is
        defining the half space.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VP = VoronoiDiagram_Polygon(V, 0)
            sage: VP.half_spaces((0, 0))
            {(-1/2*a, 1/2*a): {(1/4*a + 1/2) * x - (-1/4*a - 1/2) * y ≥ -1/4*a - 1/4},
             (1, 0): {(-1/2*a - 1/2) * x ≥ -1/4*a - 1/4}}

        """
        return {opposite_center: segment.left_half_space() for (opposite_center, segment) in self.boundary(center).items()}

    @cached_method
    def boundary(self, center):
        r"""
        Return the boundary segment that define the Voronoi cell centered at
        ``center`` in this polygon, indexed by the other center that is
        defining the segment.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VP = VoronoiDiagram_Polygon(V, 0)
            sage: VP.boundary((0, 0))
            {(-1/2*a, 1/2*a): OrientedSegment((1/2, 1/2*a + 1/2), (-1/4*a, 1/4*a)),
             (1, 0): OrientedSegment((1/2, 0), (1/2, 1/2*a + 1/2))}

        """
        from sage.all import vector
        center = vector(center)

        if center not in self.centers():
            raise ValueError("center must be a center of a Voronoi cell")

        voronoi_half_spaces = {opposite_center: self._half_space(center, opposite_center) for opposite_center in self.centers() if opposite_center != center}

        # The half spaces whose intersection is the entire polygon.
        from flatsurf.geometry.euclidean import HalfSpace
        polygon = self.polygon()
        polygon_half_spaces = [HalfSpace(vertex[0] * edge[1] - vertex[1] * edge[0], -edge[1], edge[0]) for (vertex, edge) in zip(polygon.vertices(), polygon.edges())]

        # Each segment defining this Voronoi cell is on the boundary of one of
        # the half spaces defining this Voronoi cell. Namely, if the half space
        # is not trivial in the intersection, then it contributes a segment to
        # the boundary of the cell.

        segments = HalfSpace.compact_intersection(*voronoi_half_spaces.values(), *polygon_half_spaces)

        assert segments, "Voronoi cell is empty"

        # Filter out segments that are mearly edges of the polygon; they
        # are an artifact of how we computed the segments here.
        # TODO: This is very expensive in comparison to a simple:
        #   segments = [segment for segment in segments if center not in segment]
        # But that misses some degenerate cases. But do these cases actually
        # make any sense?

        from flatsurf.geometry.euclidean import OrientedSegment
        polygon_edges = [OrientedSegment(polygon.vertex(i), polygon.vertex(i + 1)) for i in range(len(polygon.vertices()))]
        voronoi_segments = [segment for segment in segments if not any(segment.is_subset(edge) for edge in polygon_edges)]

        assert voronoi_segments

        # Associate with each segment which other center produced it.
        boundary = {}
        for segment in voronoi_segments:
            opposite_center = [c for (c, half_space) in voronoi_half_spaces.items() if half_space == segment.left_half_space()]
            assert len(opposite_center) == 1, "segment must be induced by exactly one other center of a Voronoi cell"
            opposite_center = next(iter(opposite_center))
            assert opposite_center not in boundary

            boundary[opposite_center] = segment

        return boundary

    def __repr__(self):
        return f"Voronoi diagram in polygon {self.label()}"

    def __eq__(self, other):
        r"""
        Return whether this Voronoi diagram is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon
            sage: VoronoiDiagram_Polygon(V, 0) == VoronoiDiagram_Polygon(V, 0)
            True

        """
        if not isinstance(other, VoronoiDiagram_Polygon):
            return False

        return self._parent == other._parent and self._label == other._label and self._weight == other._weight

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self._parent, self._label, self._weight))


class VoronoiPolygonCell:
    r"""
    The part of a Voronoi cell that lives entirely within a single polygon that
    makes up a surface.

    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: V = VoronoiDiagram(S, S.vertices())
        sage: V.polygon_cell(0, (0, 0))
        Voronoi cell in polygon 0 at (0, 0)

    """

    def __init__(self, parent, center):
        self._parent = parent
        self._center = center

        assert self.contains_point(center)

    def half_spaces(self):
        r"""
        Return the half spaces that delimit this cell inside its polygon,
        indexed by the two centers defining the half space.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.half_spaces()
            {(-1/2*a, 1/2*a): {(1/4*a + 1/2) * x - (-1/4*a - 1/2) * y ≥ -1/4*a - 1/4},
             (1, 0): {(-1/2*a - 1/2) * x ≥ -1/4*a - 1/4}}

        """
        return self._parent.half_spaces(self._center)

    def boundary(self):
        r"""
        Return the segments that delimit this cell inside its polygon indexed
        by the other centers defining the half space.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.boundary()
            {(-1/2*a, 1/2*a): OrientedSegment((1/2, 1/2*a + 1/2), (-1/4*a, 1/4*a)),
             (1, 0): OrientedSegment((1/2, 0), (1/2, 1/2*a + 1/2))}

        """
        return self._parent.boundary(self._center)

    def polygon(self):
        r"""
        Return the polygon which contains this cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.polygon()
            Polygon(vertices=[(0, 0), (1, 0), (1/2*a + 1, 1/2*a), (1/2*a + 1, 1/2*a + 1), (1, a + 1), (0, a + 1), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a)])

        """
        return self._parent.polygon()

    def label(self):
        r"""
        Return the label of the polygon this cell is a part of.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.label()
            0

        """
        return self._parent.label()

    def center(self):
        r"""
        Return the coordinates of the center of this Voronoi cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.center()
            (0, 0)

        """
        return self._center

    def surface(self):
        r"""
        Return the surface containing the :meth:`polygon`.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))
            sage: cell.surface() is S
            True

        """
        return self._parent.surface()

    def contains_point(self, point):
        r"""
        Return whether this cell contains the ``point``.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))

            sage: cell.contains_point((0, 0))
            True
            sage: cell.contains_point((1/2, 1/2))
            True
            sage: cell.contains_point((1, 1))
            False

        """
        if self.polygon().get_point_position(point).is_outside():
            return False

        return all(half_space.contains_point(point) for half_space in self.half_spaces().values())

    def contains_segment(self, segment):
        r"""
        Return whether the ``segment`` is a subset of this cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: cell = V.polygon_cell(0, (0, 0))

            sage: from flatsurf.geometry.euclidean import OrientedSegment
            sage: cell.contains_segment(OrientedSegment((0, 0), (1, 0)))
            False

        """
        if isinstance(segment, tuple):
            from flatsurf.geometry.euclidean import OrientedSegment
            segment = OrientedSegment(*segment)

        return self.contains_point(segment.start()) and self.contains_point(segment.end())

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

        plot = sum((segment.translate(shift)).plot(point=False) for segment in self.boundary().values())

        if plot_polygon:
            plot = graphical_polygon.plot_polygon() + plot

        return plot

    def __repr__(self):
        r"""
        Return a printable representation of this Voronoi cell.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: V.polygon_cell(0, (0, 0))
            Voronoi cell in polygon 0 at (0, 0)

        """
        return f"Voronoi cell in polygon {self.label()} at {self._center}"

    def __eq__(self, other):
        r"""
        Return whether this cell is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.voronoi import VoronoiDiagram
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: V = VoronoiDiagram(S, S.vertices())
            sage: V.polygon_cell(0, (0, 0)) == V.polygon_cell(0, (0, 1/3))
            True

        """
        if not isinstance(other, VoronoiPolygonCell):
            return False
        return self._center == other._center and self._parent == other._parent

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self._parent, self._center))
