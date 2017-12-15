from __future__ import absolute_import, print_function
from flatsurf.graphical.surface import GraphicalSurface
# The real vector space:
from flatsurf.geometry.surface_objects import SurfacePoint
from sage.plot.point import point2d

class GraphicalSurfacePoint:
    def __init__(self, graphical_surface, surface_point):
        r"""
        Create a graphical segment from a graphical surface and a SegmentInPolygon.
        """
        assert surface_point.surface() == graphical_surface.get_surface()
        self._gs = graphical_surface
        self._sp = surface_point

    def surface_point(self):
        r"""
        Return the underlying SurfacePoint.
        """
        return self._sp

    def points(self):
        r"""Return the list of points as RDF vectors."""
        point_list = []
        for label in self._sp.labels():
            if self._gs.is_visible(label):
                for coord in self._sp.coordinates(label):
                    point_list.append( self._gs.graphical_polygon(label).transform(coord) )
        return point_list

    def plot(self, **options):
        r"""
        Plot the point (which might involve drawing several dots.
        
        The options are passed to point2d.
        """
        return point2d(points=self.points(), **options)

