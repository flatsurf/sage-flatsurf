from collections import deque

from flatsurf.geometry.tangent_bundle import *
from flatsurf.geometry.polygon import is_same_direction

class SegmentInPolygon:
    r""" 
    Stores a maximal segment in a polygon of a translation surface.
    """
    def __init__(self, tangent_vector, end_vector=None):
        r""" 
        Construct a segment associated to a vector which is 
        either inside or pointed into a polygon."""
        if not end_vector is None:
            self._start=tangent_vector
            self._end=end_vector
        else:
            self._end=tangent_vector.forward_to_polygon_boundary()
            if tangent_vector.is_in_boundary_of_polygon():
                self._start=tangent_vector
            else:
                self._start=self._end.forward_to_polygon_boundary()

    def __repr__(self):
        return "Segment in polygon "+repr(self.polygon_label())+" starting at "+\
            repr(self.start_point())+" and ending at "+repr(self.end_point())

    def start(self):
        r"""
        Return a TangentVector associated to the start of a trajectory pointed forward.
        """
        return self._start

    def start_point(self):
        return self._start.point()
        
    def start_direction(self):
        return self._start.vector()

    def start_is_singular(self):
        return self._start.is_based_at_singularity()

    def end(self):
        r"""
        Return a TangentVector associated to the end of a trajectory, pointed backward.
        """
        return self._end

    def end_point(self):
        return self._end.point()
        
    def end_direction(self):
        return self._end.vector()

    def end_is_singular(self):
        return self._end.is_based_at_singularity()

    def polygon_label(self):
        return self._start.polygon_label()

    def invert(self):
        return SegmentInPolygon(self._end, end_vector=self._start)

    def next(self):
        r"""
        Return the next segment obtained by continuing straight through the end point.
        
        EXAMPLES::
        
            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: print("Polygon 0 is "+str(s.polygon(0)))
            Polygon 0 is Polygon: (0, 0), (2, -2), (2, 0)
            sage: print("Polygon 1 is "+str(s.polygon(1)))
            Polygon 1 is Polygon: (0, 0), (2, 0), (1, 3)
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector
            sage: V = tb.surface().vector_space()
            sage: v = SimilaritySurfaceTangentVector(tb, 0, V((0,0)), V((3,-1)))
            sage: from flatsurf.geometry.straight_line_trajectory import *
            sage: seg = SegmentInPolygon(v)
            sage: print(seg)
            Segment in polygon 0 starting at (0, 0) and ending at (2, -2/3)
            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: print("Polygon 0 is "+str(s.polygon(0)))
            Polygon 0 is Polygon: (0, 0), (2, -2), (2, 0)
            sage: print("Polygon 1 is "+str(s.polygon(1)))
            Polygon 1 is Polygon: (0, 0), (2, 0), (1, 3)
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector
            sage: V = tb.surface().vector_space()
            sage: v = SimilaritySurfaceTangentVector(tb, 0, V((0,0)), V((3,-1)))
            sage: from flatsurf.geometry.straight_line_trajectory import *
            sage: seg = SegmentInPolygon(v)
            sage: print(seg)
            Segment in polygon 0 starting at (0, 0) and ending at (2, -2/3)
            sage: print(seg.next())
            Segment in polygon 1 starting at (2/3, 2) and ending at (14/9, 4/3)
        """        
        if self.end_is_singular():
            raise ValueError("Cannot continue from singularity")
        return SegmentInPolygon(self._end.invert())

    def previous(self):
        if self.end_is_singular():
            raise ValueError("Cannot continue from singularity")
        return SegmentInPolygon(self._start.invert()).invert()


class StraightLineTrajectory:
    r""" Abstract class for a straight-line trajectory. """

    def __init__(self, tangent_vector):
        self._segments=deque()
        seg = SegmentInPolygon(tangent_vector)
        self._segments.append(seg)
        self._setup_forward()
        self._setup_backward()

    def segments(self):
        return self._segments

    def combinatorial_length(self):
        return len(self.segments())

    def _setup_forward(self):
        v=self.terminal_tangent_vector()
        if v.is_based_at_singularity():
            self._forward=None
        else:
            self._forward=v.invert()

    def _setup_backward(self):
        v=self.initial_tangent_vector()
        if v.is_based_at_singularity():
            self._backward=None
        else:
            self._backward=v.invert()

    def initial_segment(self):
        return self._segments[0]

    def terminal_segment(self):
        return self._segments[-1]

    def initial_tangent_vector(self):
        return self.initial_segment().start()

    def terminal_tangent_vector(self):
        return self.terminal_segment().end()

    def is_forward_separatrix(self):
        return self._forward is None

    def is_backward_separatrix(self):
        return self._backward is None

    def is_saddle_connection(self):
        return (self._forward is None) and (self._backward is None)
    
    def is_closed(self):
        return (not self.is_forward_separatrix()) and \
            self._forward.differs_by_scaling(self.initial_tangent_vector())

    def __str__(self):
        return "StraightLineTrajectory"+str(self._segments)

    def flow(self, steps):
        r"""
        Append or preprend segments to the trajectory.
        If steps is positive, attempt to append this many segments.
        If steps is negative, attempt to prepent this many segments.
        Will fail gracefully the trajectory hits a singularity or closes up.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: print("Polygon 0 is "+str(s.polygon(0)))
            Polygon 0 is Polygon: (0, 0), (2, -2), (2, 0)
            sage: print("Polygon 1 is "+str(s.polygon(1)))
            Polygon 1 is Polygon: (0, 0), (2, 0), (1, 3)
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector
            sage: V = tb.surface().vector_space()
            sage: v = SimilaritySurfaceTangentVector(tb, 0, V((1,-0.5)), V((3,-1)))
            sage: from flatsurf.geometry.straight_line_trajectory import *
            sage: traj = StraightLineTrajectory(v)
            sage: print(traj)
            StraightLineTrajectorydeque([Segment in polygon 0 starting at (1/4, -1/4) and ending at (2, -5/6)])
            sage: traj.flow(1)
            sage: print(traj)
            StraightLineTrajectorydeque([Segment in polygon 0 starting at (1/4, -1/4) and ending at (2, -5/6), Segment in polygon 1 starting at (7/12, 7/4) and ending at (61/36, 11/12)])
            sage: traj.flow(-1)
            sage: print(traj)
            StraightLineTrajectorydeque([Segment in polygon 1 starting at (15/16, 45/16) and ending at (9/8, 21/8), Segment in polygon 0 starting at (1/4, -1/4) and ending at (2, -5/6), Segment in polygon 1 starting at (7/12, 7/4) and ending at (61/36, 11/12)])
        """
        while steps>0 and \
            (not self.is_forward_separatrix()) and \
            (not self.is_closed()):
                self._segments.append(SegmentInPolygon(self._forward))
                self._setup_forward()
                steps -= 1
        while steps<0 and \
            (not self.is_backward_separatrix()) and \
            (not self.is_closed()):
                self._segments.appendleft(SegmentInPolygon(self._backward).invert())
                self._setup_backward()
                steps += 1
    
    def graphical_trajectory(self, graphical_surface):
        r""" Returns a GraphicalStraightLineTrajectory corresponding to this trajectory in the provided GraphicalSurface. """
        from flatsurf.graphical.straight_line_trajectory import GraphicalStraightLineTrajectory
        return GraphicalStraightLineTrajectory(graphical_surface, self)
