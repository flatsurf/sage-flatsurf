from flatsurf.graphical.surface import GraphicalSurface
# The real vector space:
from flatsurf.graphical.polygon import V

from flatsurf.geometry.straight_line_trajectory import SegmentInPolygon, StraightLineTrajectory

class GraphicalSegmentInPolygon:
    def __init__(self, graphical_surface, segment, color=None):
        r"""
        Create a graphical segment from a graphical surface and a SegmentInPolygon.
        """
        self._gs = graphical_surface
        self._seg = segment
        self._color = color
        label=self.polygon_label()
        self._start=self._gs.graphical_polygon(label).transform(segment.start().point())
        if self._seg.is_edge():
            self._end=self._gs.graphical_polygon(label).transform(self._seg.start().polygon().vertex(self._seg.edge()+1))
        else:
            self._end=self._gs.graphical_polygon(label).transform(segment.end().point())

    def polygon_label(self):
        return self._seg.polygon_label()

    def start(self):
        r"""Return the start point as a RDF Vector."""
        return self._start

    def end(self):
        r"""Return the end point as a RDF Vector."""
        return self._end

    def plot(self, color=None):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: v = s.tangent_vector(0, (1,-0.5), (3,-1))
            sage: from flatsurf.geometry.straight_line_trajectory import SegmentInPolygon
            sage: seg = SegmentInPolygon(v)
            sage: from flatsurf.graphical.straight_line_trajectory import *
            sage: gseg = GraphicalSegmentInPolygon(s.graphical_surface(), seg)
            sage: gseg.plot()                # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 1 graphics primitive
            sage: gseg.plot(color='red')     # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 1 graphics primitive
        """
        if self._gs.is_visible(self.polygon_label()):
            if color is None:
                if self._color is None:
                    color = "black"
                else:
                    color = self._color
            from sage.plot.line import line2d
            return line2d([self.start(), self.end()],color=color)
        else:
            from sage.plot.graphics import Graphics
            return Graphics()

class GraphicalStraightLineTrajectory:
    r"""
    Allows for the rendering of a straight-line trajectory through a graphical surface.
    """
    def __init__(self, graphical_surface, trajectory, color=None):
        self._gs = graphical_surface
        self._traj = trajectory
        self._segments = [GraphicalSegmentInPolygon(self._gs, s) for s in self._traj.segments()]
        self._color = color

    def plot(self, color=None):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface()
            sage: K.<sqrt2>=NumberField(x^2-2,embedding=1)
            sage: v = s.tangent_vector(0, (1,-1), (sqrt2,-1))
            sage: traj = v.straight_line_trajectory()
            sage: traj.flow(100)
            sage: traj.flow(-5)
            sage: gtraj = traj.graphical_trajectory(gs)
            sage: gs.plot() + gtraj.plot()      # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 119 graphics primitives
        """
        from sage.plot.graphics import Graphics
        p = Graphics()
        for seg in self._segments:
            p += seg.plot(color=(self._color if color is None else color))
        return p

