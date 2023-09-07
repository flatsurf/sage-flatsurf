# TODO: Documentation
# TODO: Drop all use of coordinates. Use EuclideanPolygonPoint instead of (label, coordinates)

class VoronoiDiagram:
    r"""
    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: center = S(0, S.polygon(0).centroid())
        sage: V = VoronoiDiagram(S, S.vertices().union([center]))
        sage: V.plot()
        Graphics object consisting of 43 graphics primitives

        sage: def weight(center):
        ....:     if center == S.polygon(0).centroid():
        ....:         return QQ(center.norm().n())
        ....:     if center in S.polygon(0).vertices():
        ....:         return 1
        ....:     raise NotImplementedError
        sage: from flatsurf.geometry.voronoi import FixedVoronoiWeight
        sage: V = VoronoiDiagram(S, S.vertices().union([center]), weight=FixedVoronoiWeight(weight))
        sage: V.plot()
        Graphics object consisting of 43 graphics primitives

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
        polygon_diagrams = {label: VoronoiDiagram_Polygon(surface.polygon(label), [coordinates for point in points for (lbl, coordinates) in point.representatives() if lbl == label], create_half_space=lambda *args: weight.create_half_space(label, *args)) for label in surface.labels()}
        self._segments = {point: [
            VoronoiCellBoundarySegment(point, surface(label, other_coordinates), label, coordinates, other_coordinates, segment) for (label, coordinates) in point.representatives() for (other_coordinates, segment) in polygon_diagrams[label].segments(coordinates).items()]
            for point in points}

    def plot(self):
        colors = ["red", "green"]
        return self._surface.plot(edge_labels=False, polygon_labels=False) + \
            sum(center.plot(color=colors[i]) for i, center in enumerate(self._segments)) + \
            sum(segment.plot(color=colors[i]) for i, point in enumerate(self._segments) for segment in self._segments[point])

    def cell(self, center):
        return self._segments[center]


class FixedVoronoiWeight:
    def __init__(self, weight=None):
        if weight is None:
            def weight(_): return 1

        self._weight = weight

    def create_half_space(self, label, center, point):
        a, b = center - point
        # TODO: Generalize this to get half spaces that are weighted by the radii of convergence.
        weights = self._weight(point), self._weight(center)
        midpoint = (weights[1] * point + weights[0] * center) / sum(weights)
        c = - a * midpoint[0] - b * midpoint[1]
        return HalfSpace(a, b, c)


class VoronoiDiagram_Polygon:
    def __init__(self, polygon, points, create_half_space):
        self._polygon = polygon

        half_spaces = {center: {point: create_half_space(center, point) for point in points if point != center} for center in points}
        self._segments = {center: self._segments(half_spaces[center]) for center in points}

    def segments(self, coordinates):
        return self._segments[coordinates]

    def _segments(self, half_spaces):
        segments = {}
        for point, half_space in half_spaces.items():
            from itertools import chain
            endpoints = set(
                list(chain.from_iterable(half_space.half_space_boundary_intersections(other) for other in half_spaces.values() if other is not half_space)) +
                list(half_space.polygon_boundary_intersection(self._polygon))
            )
            endpoints = [endpoint for endpoint in endpoints if self._polygon.contains_point(endpoint) and all(half_space.contains(endpoint) for half_space in half_spaces.values())]
            if endpoints:
                segments[point] = Segment(half_space, endpoints)

        return segments


class VoronoiCellBoundarySegment:
    def __init__(self, center, other, label, coordinates, other_coordinates, segment):
        self._center = center
        self._other = other
        self._label = label
        self._coordinates = coordinates
        self._other_coordinates = other_coordinates
        self._segment = segment

    def plot(self, **kwargs):
        return self._segment.plot(**kwargs)

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

    def plot(self, **kwargs):
        from sage.all import arrow
        return arrow(*self._endpoints, arrowsize=3, **kwargs)

    def split_vertically(self, y):
        if (self._endpoints[0][1] < y and self._endpoints[1][1] > y) or (self._endpoints[0][1] > y and self._endpoints[1][1] < y):
            split = self._half_space.half_space_boundary_intersections(HalfSpace(0, 1, -y))[0]
            return [Segment(self._half_space, [self._endpoints[0], split]), Segment(self._half_space, [split, self._endpoints[1]])]

        return [self]

    def endpoints(self):
        return self._endpoints
