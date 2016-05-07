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
        either inside or pointed into a polygon.
        """
        if not end_vector is None:
            self._start = tangent_vector
            self._end = end_vector
        else:
            self._end = tangent_vector.forward_to_polygon_boundary()
            if tangent_vector.is_in_boundary_of_polygon():
                self._start = tangent_vector
            else:
                self._start = self._end.forward_to_polygon_boundary()

    def __repr__(self):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.straight_line_trajectory import SegmentInPolygon
            sage: s = similarity_surfaces.example()
            sage: v = s.tangent_vector(0, (0,0), (3,-1))
            sage: SegmentInPolygon(v)
            Segment in polygon 0 starting at (0, 0) and ending at (2, -2/3)
        """
        return "Segment in polygon {} starting at {} and ending at {}".format(
                self.polygon_label(), self.start_point(), self.end_point())

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

    def is_edge(self):
        if not self.start_is_singular() or not self.end_is_singular():
            return False
        vv=self.start().vector()
        vertex=self.start().singularity()
        ww=self.start().polygon().edge(vertex)
        from flatsurf.geometry.polygon import is_same_direction
        return is_same_direction(vv,ww)

    def edge(self):
        if not self.is_edge():
            raise ValueError("Segment asked for edge when not an edge")
        return self.start().singularity()

    def polygon_label(self):
        return self._start.polygon_label()

    def invert(self):
        return SegmentInPolygon(self._end, end_vector=self._start)

    def next(self):
        r"""
        Return the next segment obtained by continuing straight through the end point.
        
        EXAMPLES::
        
            sage: from flatsurf import *
            sage: from flatsurf.geometry.straight_line_trajectory import SegmentInPolygon

            sage: s = similarity_surfaces.example()
            sage: s.polygon(0)
            Polygon: (0, 0), (2, -2), (2, 0)
            sage: s.polygon(1)
            Polygon: (0, 0), (2, 0), (1, 3)
            sage: v = s.tangent_vector(0, (0,0), (3,-1))
            sage: seg = SegmentInPolygon(v)
            sage: seg
            Segment in polygon 0 starting at (0, 0) and ending at (2, -2/3)
            sage: seg.next()
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
    r"""
    Straight-line trajectory in a translation surface.
    """
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

    def __repr__(self):
        start = self._segments[0]
        end = self._segments[-1]
        return "Straight line trajectory made of {} segments from {} in polygon {} to {} in polygon {}".format(
                len(self._segments),
                start.start_point(), start.start().polygon_label(),
                end.end_point(), end.end().polygon_label())

    def flow(self, steps):
        r"""
        Append or preprend segments to the trajectory.
        If steps is positive, attempt to append this many segments.
        If steps is negative, attempt to prepend this many segments.
        Will fail gracefully the trajectory hits a singularity or closes up.

        EXAMPLES::

            sage: from flatsurf import *

            sage: s = similarity_surfaces.example()
            sage: v = s.tangent_vector(0, (1,-1/2), (3,-1))
            sage: traj = v.straight_line_trajectory()
            sage: traj
            Straight line trajectory made of 1 segments from (1/4, -1/4) in polygon 0 to (2, -5/6) in polygon 0
            sage: traj.flow(1)
            sage: traj
            Straight line trajectory made of 2 segments from (1/4, -1/4) in polygon 0 to (61/36, 11/12) in polygon 1
            sage: traj.flow(-1)
            sage: traj
            Straight line trajectory made of 3 segments from (15/16, 45/16) in polygon 1 to (61/36, 11/12) in polygon 1
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
        r"""
        Returns a GraphicalStraightLineTrajectory corresponding to this trajectory in the provided GraphicalSurface.
        """
        from flatsurf.graphical.straight_line_trajectory import GraphicalStraightLineTrajectory
        return GraphicalStraightLineTrajectory(graphical_surface, self)
