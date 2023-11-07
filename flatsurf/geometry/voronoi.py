# TODO: Documentation
# TODO: Drop all use of coordinates. Use EuclideanPolygonPoint instead of (label, coordinates)
# TODO: SageMath already has a VoronoiDiagram so maybe not use that name. Explain that we cannot use the builtin version since it has no weights.

from sage.misc.cachefunc import cached_method

class VoronoiDiagram:
    r"""
    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram  # random output due to deprecation warnings
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: center = S(0, S.polygon(0).centroid())
        sage: V = VoronoiDiagram(S, S.vertices().union([center]))
        sage: V.plot()
        Graphics object consisting of 43 graphics primitives

        sage: def weight(label, center):
        ....:     if center == S.polygon(0).centroid():
        ....:         return QQ(center.norm().n())
        ....:     if center in S.polygon(0).vertices():
        ....:         return 1
        ....:     raise NotImplementedError
        sage: from flatsurf.geometry.voronoi import FixedVoronoiWeight
        sage: V = VoronoiDiagram(S, S.vertices().union([center]), weight=FixedVoronoiWeight(weight))
        sage: V.plot()
        Graphics object consisting of 43 graphics primitives

    The same Voronoi diagram but starting from a more complicated description
    of the octagon::

        sage: from flatsurf import Polygon, translation_surfaces, polygons
        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: S = translation_surfaces.regular_octagon()
        sage: S = S.subdivide().codomain()
        sage: V = VoronoiDiagram(S, S.vertices())
        sage: V.plot()
        Graphics object consisting of 75 graphics primitives

        sage: def weight(label, center):
        ....:     P = S.polygon(label)
        ....:     vertex = list(P.vertices()).index(center)
        ....:     angle = P.angle(vertex)
        ....:     if angle == 1/8:
        ....:         return polygons.regular_ngon(8).centroid().norm().n()
        ....:     assert angle == 3/16, angle
        ....:     return 1

        sage: from flatsurf.geometry.voronoi import FixedVoronoiWeight
        sage: V = VoronoiDiagram(S, S.vertices(), weight=FixedVoronoiWeight(weight))
        sage: V.plot()
        Graphics object consisting of 75 graphics primitives

    We can also choose strange such that there is a saddle connection between a
    vertex and itself without leaving the cell::

        sage: def weight(label, center):
        ....:     P = S.polygon(label)
        ....:     vertex = list(P.vertices()).index(center)
        ....:     angle = P.angle(vertex)
        ....:     if angle == 1/8:
        ....:         return polygons.regular_ngon(8).centroid().norm().n()
        ....:     assert angle == 3/16, angle
        ....:     return 1/6

        sage: from flatsurf.geometry.voronoi import FixedVoronoiWeight
        sage: V = VoronoiDiagram(S, S.vertices(), weight=FixedVoronoiWeight(weight))
        sage: V.plot()
        Graphics object consisting of 59 graphics primitives
    """

    def __init__(self, surface, points, weight=None):
        if weight is None:
            weight = FixedVoronoiWeight()

        points = set(points)

        if not surface.vertices().issubset(points):
            raise NotImplementedError("can only compute Voronoi diagrams when all vertices are centers")

        self._surface = surface
        self._points = points

        # TODO: These mapping structures are too complex.
        self._polygon_diagrams = {label: VoronoiDiagram_Polygon(surface.polygon(label), [coordinates for point in points for (lbl, coordinates) in point.representatives() if lbl == label], create_half_space=lambda *args: weight.create_half_space(label, *args)) for label in surface.labels()}
        self._segments = {point: [
            VoronoiCellBoundarySegment(point, surface(label, other_center), label, center, other_center, segment) for (label, center) in point.representatives() for (other_center, segment) in self._polygon_diagrams[label].boundary_segments(center).items()]
            for point in points}

    def plot(self, gs=None):
        # TODO: Generalize
        colors = ["red", "green", "orange", "purple", "cyan", "blue"]
        if gs is None:
            gs = self._surface.graphical_surface(edge_labels=False, polygon_labels=False)
        return gs.plot() + \
            sum(center.plot(gs, color=colors[i]) for i, center in enumerate(self._segments)) + \
            sum(segment.plot(gs=gs, color=colors[i]) for i, point in enumerate(self._segments) for segment in self._segments[point])

    def cell(self, center):
        return self._segments[center]

    def relativize(self, point):
        r"""
        Return which Voronoi cell ``point`` lies in.

        Returns the center of the cell as (label, coordinates) and a vector Δ
        from that center to the ``point``.
        """
        for label, coordinates in point.representatives():
            center = self._polygon_diagrams[label].relativize(coordinates)
            return (label, center), coordinates - center


class FixedVoronoiWeight:
    def __init__(self, weight=None):
        if weight is None:
            def weight(_, __): return 1

        self._weight = weight

    def create_half_space(self, label, center, point):
        a, b = center - point
        # TODO: Generalize this to get half spaces that are weighted by the radii of convergence.
        weights = self._weight(label, point), self._weight(label, center)
        midpoint = (weights[1] * point + weights[0] * center) / sum(weights)
        c = a * midpoint[0] + b * midpoint[1]

        from sage.all import Polyhedron
        return Polyhedron(ieqs=[(-c, a, b)])


class VoronoiDiagram_Polygon:
    r"""
    The part of the Voronoi diagram inside ``polygon`` with cells centered at
    ``points``.

    TESTS::

        sage: from flatsurf import Polygon
        sage: K.<a> = QuadraticField(2)
        sage: P = Polygon(vertices=([(1, 0), (1/2*a + 1, 1/2*a), (1/2, 1/2*a + 1/2)]))
        sage: points = list(P.vertices())
        sage: points = points[1:] + points[:1]

        sage: def weight(_, center):
        ....:     assert center in points
        ....:     if center == points[2]:
        ....:         return 1/8
        ....:     return 1

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram_Polygon, FixedVoronoiWeight
        sage: _ = VoronoiDiagram_Polygon(P, points, lambda *args: FixedVoronoiWeight(weight).create_half_space(None, *args))

    """

    def __init__(self, polygon, centers, create_half_space):
        self._polygon = polygon
        self._centers = centers
        self._create_half_space = create_half_space

        # Determine the segments that bound each Voronoi cell within this
        # polygon, indexed by the coordinates of the center of the Voronoi
        # cell.
        self._cell_boundary = {}
        for center in centers:
            self._cell_boundary[center] = {}

            # Each segment is on the boundary of one of the half spaces
            # defining this Voronoi cell. Namely, if the half space is not
            # trivial in the intersection, then it contributes a segment to the
            # boundary of the cell.
            half_spaces = list(self._boundary_half_spaces(center).values()) + self._polygon_half_spaces()

            from itertools import chain
            from sage.all import Polyhedron
            intersection = Polyhedron(ieqs=chain.from_iterable(half_space.inequalities() for half_space in half_spaces))

            assert not intersection.is_empty()
            assert intersection.dimension() == 2

            segments = intersection.faces(1)
            assert all(segment.is_compact() for segment in segments)

            segments = [Polyhedron(vertices=segment.vertices()) for segment in segments]

            # Filter out segments that are nearly edges of the polygon; they
            # are an artifact of how we computed the segments here.
            # TODO: This is very expensive in comparison to a simple:
            #   segments = [segment for segment in segments if center not in segment]
            # But is it actually correct? Does it make any sense to have such
            # Voronoi cells? (I mean, they are not really Voronoi cells since
            # the boundaries should not be straight line segments but curved
            # anyway…
            segments = [segment for segment in segments if all(segment.intersection(edge_half_space.faces(1)[0].as_polyhedron()) != segment for edge_half_space in self._polygon_half_spaces())]

            assert segments

            # Associate with each segment which pair of vertices produced it,
            # i.e., apart from center which other point.
            for segment in segments:
                opposite_centers = [opposite_center for opposite_center, half_space in self._boundary_half_spaces(center).items() if segment.intersection(half_space.faces(1)[0].as_polyhedron()) == segment]
                assert len(opposite_centers) == 1, "segment must be induced by exactly one other center of a Voronoi cell"

                opposite_center = opposite_centers[0]
                assert opposite_center not in self._cell_boundary[center]

                self._cell_boundary[center][opposite_center] = segment

    @cached_method
    def _boundary_half_spaces(self, center):
        r"""
        Return half spaces whose intersection is the Voronoi cell at
        ``center`` one for each other center of a Voronoi cell.
        """
        return {opposite_center: self._create_half_space(center, opposite_center) for opposite_center in self._centers if opposite_center != center}

    @cached_method
    def _polygon_half_spaces(self):
        r"""
        Return the half spaces whose intersection is the underlying polygon.
        """
        # TODO: Implement this in polygon?
        from sage.all import Polyhedron
        return [Polyhedron(ieqs=[(vertex[0] * edge[1] - vertex[1] * edge[0], -edge[1], edge[0])]) for (vertex, edge) in zip(self._polygon.vertices(), self._polygon.edges())]

    def relativize(self, coordinates):
        r"""
        Return the center of the Voronoi cell which ``coordinates`` lies in.
        """
        for center in self._half_spaces:
            if all([half_space.contains(coordinates) for half_space in self._half_spaces[center].values()]):
                return center

        assert False

    def boundary_segments(self, center):
        r"""
        Return the segments bounding the Voronoi cell centered at the
        coordinates ``center`` in this polygon.
        """
        return self._cell_boundary[center]


class VoronoiCellBoundarySegment:
    def __init__(self, center, other, label, coordinates, other_coordinates, segment):
        self._center = center
        self._other = other
        self._label = label
        self._coordinates = coordinates
        self._other_coordinates = other_coordinates
        self._segment = segment

    def plot(self, gs=None, **kwargs):
        shift = None
        if gs is not None:
            graphical_polygon = gs.graphical_polygon(self._label)
            polygon = graphical_polygon.base_polygon()
            shift = graphical_polygon.transformed_vertex(0) - polygon.vertex(0)
        return (self._segment + shift).plot(point=False, **kwargs)

    def center(self):
        return self._center, self._label, self._coordinates

    def opposite_center(self):
        return self._other, self._label, self._other_coordinates

    def segment(self):
        return self._label, self._segment

    def segments_with_uniform_roots(self):
        segments = [self]

        import itertools

        if self._center.is_vertex():
            segments = itertools.chain.from_iterable(segment.split_vertically(self._coordinates[1]) for segment in segments)
        if self._other.is_vertex():
            segments = itertools.chain.from_iterable(segment.split_vertically(self._other_coordinates[1]) for segment in segments)

        return list(segments)

    def split_vertically(self, y):
        return [VoronoiCellBoundarySegment(self._center, self._other, self._label, self._coordinates, self._other_coordinates, segment) for segment in self._segment.split_vertically(y)]

    def __ne__(self, other):
        return not (self == other)

    def __eq__(self, other):
        if not isinstance(other, VoronoiCellBoundarySegment):
            return False

        return self._center == other._center and self._other == other._other and self._label == other._label and self._coordinates == other._coordinates and self._other_coordinates == other._other_coordinates and self._segment == other._segment

    def __neg__(self):
        return VoronoiCellBoundarySegment(self._other, self._center, self._label, self._other_coordinates, self._coordinates, -self._segment)


# TODO: Deleteme
class HalfSpace:
    r"""
    The half space ax + by + c ≥ 0.
    """

    def __init__(self, a, b, c):
        self._a = a
        self._b = b
        self._c = c

    def plot(self, polygon):
        endpoints = self.polygon_boundary_intersection(polygon)

        from sage.all import line
        return line(endpoints)

    def __repr__(self):
        return f"({self._a})*x + ({self._b})*y + {self._c} ≥ 0"

    def half_space_boundary_intersections(self, other):
        if not isinstance(other, HalfSpace):
            raise NotImplementedError

        from sage.all import matrix, vector
        A = matrix([
            [self._a, self._b],
            [other._a, other._b],
        ])
        b = vector((-self._c, -other._c))
        try:
            intersection = A.solve_right(b)
        except ValueError:
            return []

        intersection.set_immutable()

        return [intersection]

    def polygon_boundary_intersection(self, polygon):
        intersections = set()
        for (vertex, edge) in zip(polygon.vertices(), polygon.edges()):
            a, b = -edge[1], edge[0]
            c = - a*vertex[0] - b*vertex[1]

            from sage.all import matrix, vector
            A = matrix([
                [self._a, self._b],
                [a, b],
            ])
            b = vector((-self._c, -c))

            try:
                intersection = A.solve_right(b)
            except ValueError:
                continue

            if not polygon.contains_point(intersection):
                continue

            intersection.set_immutable()

            intersections.add(intersection)

        assert len(intersections) <= 2 or not polygon.is_convex()
        return intersections

    def contains(self, point):
        return self._a * point[0] + self._b * point[1] + self._c >= 0

    def __eq__(self, other):
        if not isinstance(other, HalfSpace):
            return False

        from sage.all import sgn

        if sgn(self._a) != sgn(other._a):
            return False

        if sgn(self._b) != sgn(other._b):
            return False

        if sgn(self._c) != sgn(other._c):
            return False

        return self._a * other._b == self._b * other._a and self._a * other._c == self._c * other._a and self._b * other._c == self._c * other._b

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        return HalfSpace(-self._a, -self._b, -self._c)


class Segment:
    def __init__(self, half_space, endpoints):
        if len(endpoints) != 2:
            raise NotImplementedError

        # Sort endpoints so that they are in counterclockwise order when seen
        # from the center.
        if not isinstance(half_space, HalfSpace):
            raise NotImplementedError
        endpoints.sort(key=lambda xy: xy[0] * half_space._b - xy[1] * half_space._a)

        self._half_space = half_space
        self._endpoints = endpoints

    def plot(self, shift=None, **kwargs):
        from sage.all import arrow
        endpoints = [e + shift for e in self._endpoints]
        return arrow(*endpoints, arrowsize=3, **kwargs)

    def split_vertically(self, y):
        if (self._endpoints[0][1] < y and self._endpoints[1][1] > y) or (self._endpoints[0][1] > y and self._endpoints[1][1] < y):
            split = self._half_space.half_space_boundary_intersections(HalfSpace(0, 1, -y))[0]
            return [Segment(self._half_space, [self._endpoints[0], split]), Segment(self._half_space, [split, self._endpoints[1]])]

        return [self]

    def endpoints(self):
        return self._endpoints

    def vector(self):
        if len(self._endpoints) != 2:
            raise NotImplementedError

        from sage.all import vector
        return vector((self._endpoints[1][0] - self._endpoints[0][0], self._endpoints[1][1] - self._endpoints[0][1]))

    def __eq__(self, other):
        if not isinstance(other, Segment):
            return False

        return self._half_space == other._half_space and self._endpoints == other._endpoints

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        return Segment(-self._half_space, self._endpoints[::-1])
