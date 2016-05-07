from flatsurf.graphical.surface import GraphicalSurface
# The real vector space:
from flatsurf.graphical.polygon import V

from flatsurf.geometry.straight_line_trajectory import SegmentInPolygon, StraightLineTrajectory

class GraphicalSegmentInPolygon:
    def __init__(self, graphical_surface, segment):
        r"""
        Create a graphical segment from a graphical surface and a SegmentInPolygon.
        """
        self._gs = graphical_surface
        self._seg=segment
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
        
            from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            s=SimilaritySurfaceGenerators.example()
            from flatsurf.graphical.surface import GraphicalSurface
            gs=GraphicalSurface(s)
            gs.make_visible(1)
            from flatsurf.geometry.tangent_bundle import *
            tb = SimilaritySurfaceTangentBundle(s)
            V=tb.surface().vector_space()
            v=SimilaritySurfaceTangentVector(tb, 0, V((1,-0.5)), V((3,-1)))
            from flatsurf.geometry.straight_line_trajectory import *
            seg = SegmentInPolygon(v)
            from flatsurf.graphical.straight_line_trajectory import *
            gseg = GraphicalSegmentInPolygon(gs, seg)
            show(gs.plot()+gseg.plot())
        """
        if self._gs.is_visible(self.polygon_label()):
            if color==None:
                color="black"
            from sage.plot.line import line2d
            return line2d([self.start(), self.end()],color=color)
        else:
            from sage.plot.graphics import Graphics
            return Graphics()

class GraphicalStraightLineTrajectory:
    r"""Allows for the rendering of a straight-line trajectory through a graphical surface."""
    
    def __init__(self, graphical_surface, trajectory):
        self._gs = graphical_surface
        self._traj = trajectory
        self._segments = [GraphicalSegmentInPolygon(self._gs, s) for s in self._traj.segments()]

    def plot(self, color=None):
        r"""
        EXAMPLES::
        
            from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            s=SimilaritySurfaceGenerators.example()
            from flatsurf.graphical.surface import GraphicalSurface
            gs=GraphicalSurface(s)
            gs.make_visible(1)
            from flatsurf.geometry.tangent_bundle import *
            K.<sqrt2>=NumberField(x^2-2,embedding=1)
            tb = SimilaritySurfaceTangentBundle(s)
            from flatsurf.sage.modules.free_module_element import vector
            v=SimilaritySurfaceTangentVector(tb, 0, vector((1,-1)), vector((sqrt2,-1)))
            from flatsurf.geometry.straight_line_trajectory import *
            traj = StraightLineTrajectory(v)
            traj.flow(100)
            traj.flow(-5)
            from flatsurf.graphical.straight_line_trajectory import *
            gtraj = GraphicalStraightLineTrajectory(gs, traj)
            show(gs.plot()+gtraj.plot())
        """        
        from sage.plot.graphics import Graphics
        p = Graphics()
        for seg in self._segments:
            p += seg.plot(color=color)
        return p

