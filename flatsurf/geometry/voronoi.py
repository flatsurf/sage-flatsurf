class VoronoiDiagram:
    r"""
    EXAMPLES::

        sage: from flatsurf.geometry.voronoi import VoronoiDiagram
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: center = S(0, S.polygon(0).centroid())
        sage: V = VoronoiDiagram(S, S.vertices().union([center]))

    """

    def __init__(self, surface, points, weight=None):
        points = set(points)

        if not surface.vertices().issubset(points):
            raise NotImplementedError("can only compute Voronoi diagrams when all vertices are centers")

        self._surface = surface
        self._points = points

        polygon_diagrams = {label: VoronoiDiagram_Polygon(surface.polygon(label), [coordinates for point in points for (lbl, coordinates) in point.representatives() if lbl == label], weight=weight) for label in surface.labels()}
        self._segments = {point: [(label, segment) for (label, coordinates) in point.representatives() for segment in polygon_diagrams(label).segments(coordinates)] for point in points}


class VoronoiDiagram_Polygon:
    def __init__(self, polygon, points, weight=None):
        self._polygon = polygon
        self._points = points
        self._weight = weight

        half_spaces = {center: self._half_spaces(center) for center in points}
        self._segments = {center: self._segments(self, center, half_spaces[center]) for center in points}

    def segments(self, coordinates):
        raise NotImplementedError()

    def _half_spaces(self, center):
        raise NotImplementedError()

    def _segments(self, center, half_spaces):
        raise NotImplementedError()
