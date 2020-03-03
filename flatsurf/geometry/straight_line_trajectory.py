from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

from collections import deque, defaultdict

from .polygon import is_same_direction, line_intersection
from .surface_objects import SaddleConnection

# Vincent question:
# using deque has the disadvantage of losing the initial points
# ideally doig
#  my_line[i]
# we should always access to the same element

# I wanted to be able to flow backward thus inserting at the beginning of a list.
# Perhaps it would be better to model this on a deque-like class that is indexed by
# all integers rather than just the non-negative ones? Do you know of such
# a class? Alternately, we could store an offset.

def get_linearity_coeff(u, v):
    r"""
    Given the two 2-dimensional vectors ``u`` and ``v``, return ``a`` so that
    ``v = a*u``

    If the vectors are not colinear, a ``ValueError`` is raised.

    EXAMPLES::

        sage: from flatsurf.geometry.straight_line_trajectory import get_linearity_coeff

        sage: V = VectorSpace(QQ,2)
        sage: get_linearity_coeff(V((1,0)), V((2,0)))
        2
        sage: get_linearity_coeff(V((2,0)), V((1,0)))
        1/2
        sage: get_linearity_coeff(V((0,1)), V((0,2)))
        2
        sage: get_linearity_coeff(V((0,2)), V((0,1)))
        1/2
        sage: get_linearity_coeff(V((1,2)), V((-2,-4)))
        -2

        sage: get_linearity_coeff(V((1,1)), V((-1,1)))
        Traceback (most recent call last):
        ...
        ValueError: non colinear
    """
    if u[0]:
        a = v[0]/u[0]
        if v[1] != a*u[1]:
            raise ValueError("non colinear")
        return a
    elif v[0]:
        raise ValueError("non colinear")
    elif u[1]:
        return v[1]/u[1]
    else:
        raise ValueError("zero vector")

class SegmentInPolygon:
    r"""
    Maximal segment in a polygon of a similarity surface

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.straight_line_trajectory import SegmentInPolygon
        sage: s = similarity_surfaces.example()
        sage: v = s.tangent_vector(0, (1/3,-1/4), (0,1))
        sage: SegmentInPolygon(v)
        Segment in polygon 0 starting at (1/3, -1/3) and ending at (1/3, 0)
    """
    def __init__(self, start, end=None):
        if not end is None:
            # WARNING: here we assume that both start and end are on the
            # boundary
            self._start = start
            self._end = end
        else:
            self._end = start.forward_to_polygon_boundary()
            self._start = self._end.forward_to_polygon_boundary()

    def __eq__(self, other):
        return type(self) is type(other) and \
               self._start == other._start and \
               self._end == other._end

    def __ne__(self, other):
        return type(self) is not type(other) or \
               self._start != other._start or \
               self._end != other._end

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
                self.polygon_label(), self.start().point(), self.end().point())

    def start(self):
        r"""
        Return the tangent vector associated to the start of a trajectory pointed forward.
        """
        return self._start

    def start_is_singular(self):
        return self._start.is_based_at_singularity()

    def end(self):
        r"""
        Return a TangentVector associated to the end of a trajectory, pointed backward.
        """
        return self._end

    def end_is_singular(self):
        return self._end.is_based_at_singularity()

    def is_edge(self):
        if not self.start_is_singular() or not self.end_is_singular():
            return False
        vv=self.start().vector()
        vertex=self.start().vertex()
        ww=self.start().polygon().edge(vertex)
        from flatsurf.geometry.polygon import is_same_direction
        return is_same_direction(vv,ww)

    def edge(self):
        if not self.is_edge():
            raise ValueError("Segment asked for edge when not an edge")
        return self.start().vertex()

    def polygon_label(self):
        return self._start.polygon_label()

    def invert(self):
        return SegmentInPolygon(self._end, self._start)

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


    # DEPRECATED STUFF THAT WILL BE REMOVED

    def start_point(self):
        from sage.misc.superseded import deprecation
        deprecation(1, "do not use start_point but start().point()")
        return self._start.point()

    def start_direction(self):
        from sage.misc.superseded import deprecation
        deprecation(1, "do not use start_direction but start().vector()")
        return self._start.vector()

    def end_point(self):
        from sage.misc.superseded import deprecation
        deprecation(1, "do not use end_point but end().point()")
        return self._end.point()

    def end_direction(self):
        from sage.misc.superseded import deprecation
        deprecation(1, "do not use end_direction but end().vector()")
        return self._end.vector()

class AbstractStraightLineTrajectory:
    r"""
    You need to implement:

    - ``def segment(self, i)``
    - ``def segments(self)``
    """
    def surface(self):
        raise NotImplementedError

    def __repr__(self):
        start = self.segment(0).start()
        end = self.segment(-1).end()
        return "Straight line trajectory made of {} segments from {} in polygon {} to {} in polygon {}".format(
                self.combinatorial_length(),
                start.point(), start.polygon_label(),
                end.point(), end.polygon_label())

    def plot(self, *args, **options):
        r"""
        Plot this trajectory by converting to a graphical trajectory.

        If any arguments are provided in `*args` it must be only one argument containing a GraphicalSurface.
        The keyword arguments in `**options` are passed on to :func:`GraphicalStraightLineTrajectory.plot`.

        EXAMPLES::

            sage: from flatsurf import *
            sage: T = translation_surfaces.square_torus()
            sage: v = T.tangent_vector(0, (0,0), (5,7))
            sage: L = v.straight_line_trajectory()
            sage: L.plot()               # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 1 graphics primitive
            sage: L.plot(color='red')    # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 1 graphics primitive
        """
        if len(args) > 1:
            raise ValueError("SimilaritySurface.plot() can take at most one non-keyword argument.")
        if len(args)==1:
            from flatsurf.graphical.surface import GraphicalSurface
            if not isinstance(args[0], GraphicalSurface):
                raise ValueError("If an argument is provided, it must be a GraphicalSurface.")
            return self.graphical_trajectory(graphical_surface = args[0]).plot(**options)
        return self.graphical_trajectory().plot(**options)

    def graphical_trajectory(self, graphical_surface=None, **options):
        r"""
        Returns a ``GraphicalStraightLineTrajectory`` corresponding to this
        trajectory in the provided  ``GraphicalSurface``.
        """
        from flatsurf.graphical.straight_line_trajectory import GraphicalStraightLineTrajectory
        if graphical_surface is None:
            graphical_surface = self.surface().graphical_surface()
        return GraphicalStraightLineTrajectory(self, graphical_surface, **options)

    def cylinder(self):
        r"""
        If this is a closed orbit, return the associated maximal cylinder.
        Raises a ValueError if this trajectory is not closed.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.regular_octagon()
            sage: v = s.tangent_vector(0,(1/2,0),(sqrt(2),1))
            sage: traj = v.straight_line_trajectory()
            sage: traj.flow(4)
            sage: traj.is_closed()
            True
            sage: cyl = traj.cylinder()
            sage: cyl.area()                # a = sqrt(2)
            a + 1
            sage: cyl.holonomy()
            (3*a + 4, 2*a + 3)
            sage: cyl.edges()
            (2, 3, 3, 2, 4)
        """
        # Note may not be defined.
        if not self.is_closed():
            raise ValueError("Cylinder is only defined for closed straight-line trajectories.")
        from .surface_objects import Cylinder
        coding = self.coding()
        label = coding[0][0]
        edges = [ e for l,e in coding[1:] ]
        edges.append(self.surface().opposite_edge(coding[0][0],coding[0][1])[1])
        return Cylinder(self.surface(), label, edges)

    def coding(self, alphabet=None):
        r"""
        Return the coding of this trajectory with respect to the sides of the
        polygons

        INPUT:

        - ``alphabet`` -- an optional dictionary ``(lab,nb) -> letter``. If some
          labels are avoided then these crossings are ignored.

        EXAMPLES::

            sage: from flatsurf import *
            sage: t = translation_surfaces.square_torus()

            sage: v = t.tangent_vector(0, (1/2,0), (5,6))
            sage: l = v.straight_line_trajectory()
            sage: alphabet = {(0,0): 'a', (0,1): 'b', (0,2):'a', (0,3): 'b'}
            sage: l.coding()
            [(0, 0), (0, 1)]
            sage: l.coding(alphabet)
            ['a', 'b']
            sage: l.flow(10); l.flow(-10)
            sage: l.coding()
            [(0, 2), (0, 1), (0, 2), (0, 1), (0, 2), (0, 1), (0, 2), (0, 1), (0, 2)]
            sage: print(''.join(l.coding(alphabet)))
            ababababa

            sage: v = t.tangent_vector(0, (1/2,0), (7,13))
            sage: l = v.straight_line_trajectory()
            sage: l.flow(10); l.flow(-10)
            sage: print(''.join(l.coding(alphabet)))
            aabaabaababaabaabaab

        For a closed trajectory, the last label (corresponding also to the
        starting point) is not considered::

            sage: v = t.tangent_vector(0, (1/5,1/7), (1,1))
            sage: l = v.straight_line_trajectory()
            sage: l.flow(10)
            sage: l.is_closed()
            True
            sage: l.coding(alphabet)
            ['a', 'b']

        Check that the saddle connections that are obtained in the torus get the
        expected coding::

            sage: for _ in range(10):
            ....:     x = ZZ.random_element(1,30)
            ....:     y = ZZ.random_element(1,30)
            ....:     x,y = x/gcd(x,y), y/gcd(x,y)
            ....:     v = t.tangent_vector(0, (0,0), (x,y))
            ....:     l = v.straight_line_trajectory()
            ....:     l.flow(200); l.flow(-200)
            ....:     w = ''.join(l.coding(alphabet))
            ....:     assert Word(w+'ab'+w).is_balanced()
            ....:     assert Word(w+'ba'+w).is_balanced()
            ....:     assert w.count('a') == y-1
            ....:     assert w.count('b') == x-1
        """
        ans = []

        segments = self.segments()

        s = segments[0]
        start = s.start()
        if start._position._position_type == start._position.EDGE_INTERIOR:
            p = s.polygon_label()
            e = start._position.get_edge()
            lab = (p,e) if alphabet is None else alphabet.get((p,e))
            if lab is not None:
                ans.append(lab)

        for i in range(len(segments)-1):
            s = segments[i]
            end = s.end()
            p = s.polygon_label()
            e = end._position.get_edge()
            lab = (p,e) if alphabet is None else alphabet.get((p,e))
            if lab is not None:
                ans.append(lab)

        s = segments[-1]
        end = s.end()
        if end._position._position_type == end._position.EDGE_INTERIOR and \
           end.invert() != start:
            p = s.polygon_label()
            e = end._position.get_edge()
            lab = (p,e) if alphabet is None else alphabet.get((p,e))
            if lab is not None:
                ans.append(lab)

        return ans

    def initial_tangent_vector(self):
        return self.segment(0).start()

    def terminal_tangent_vector(self):
        return self.segment(-1).end()

    def intersects(self, traj, count_singularities = False):
        r"""
        Return true if this trajectory intersects the other trajectory.
        """
        try:
            next(self.intersections(traj, count_singularities = count_singularities))
        except StopIteration:
            return False
        return True

    def intersections(self, traj, count_singularities = False, include_segments = False):
        r"""
        Return the set of SurfacePoints representing the intersections
        of this trajectory with the provided trajectory or SaddleConnection.

        Singularities will be included only if count_singularities is
        set to True.

        If include_segments is True, it iterates over triples consisting of the SurfacePoint,
        and two sets. The first set consists of segments of this trajectory that contain the point
        and the second set consists of segments of traj that contain the point.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.square_torus()
            sage: traj1 = s.tangent_vector(0,(1/2,0),(1,1)).straight_line_trajectory()
            sage: traj1.flow(3)
            sage: traj1.is_closed()
            True
            sage: traj2 = s.tangent_vector(0,(1/2,0),(-1,1)).straight_line_trajectory()
            sage: traj2.flow(3)
            sage: traj2.is_closed()
            True
            sage: sum(1 for _ in traj1.intersections(traj2))
            2
        """
        # Partition the segments making up the trajectories by label.
        if isinstance(traj,SaddleConnection):
            traj = traj.trajectory()
        lab_to_seg1 = {}
        for seg1 in self.segments():
            label = seg1.polygon_label()
            if label in lab_to_seg1:
                lab_to_seg1[label].append(seg1)
            else:
                lab_to_seg1[label] = [seg1]
        lab_to_seg2 = {}
        for seg2 in traj.segments():
            label = seg2.polygon_label()
            if label in lab_to_seg2:
                lab_to_seg2[label].append(seg2)
            else:
                lab_to_seg2[label] = [seg2]
        intersection_points = set()
        if include_segments:
            segments={}
        for label,seg_list_1 in iteritems(lab_to_seg1):
            if label in lab_to_seg2:
                seg_list_2 = lab_to_seg2[label]
                for seg1 in seg_list_1:
                    for seg2 in seg_list_2:
                        x = line_intersection(seg1.start().point(),
                                              seg1.start().point()+seg1.start().vector(),
                                              seg2.start().point(),
                                              seg2.start().point()+seg2.start().vector())
                        if x is not None:
                            pos = self._s.polygon(seg1.polygon_label()).get_point_position(x)
                            if pos.is_inside() and (count_singularities or not pos.is_vertex()):
                                new_point = self._s.surface_point(seg1.polygon_label(),x)
                                if new_point not in intersection_points:
                                    intersection_points.add(new_point)
                                    if include_segments:
                                        segments[new_point]=({seg1},{seg2})
                                    else:
                                        yield new_point
                                elif include_segments:
                                    segments[new_point][0].append(seg1)
                                    segments[new_point][1].append(seg2)
        if include_segments:
            for x in iteritems(segments):
                yield x



class StraightLineTrajectory(AbstractStraightLineTrajectory):
    r"""
    Straight-line trajectory in a similarity surface.

    EXAMPLES::

        # Demonstrate the handling of edges
        sage: from flatsurf import *
        sage: from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectory
        sage: p = SymmetricGroup(2)('(1,2)')
        sage: s = translation_surfaces.origami(p,p)
        sage: traj = StraightLineTrajectory(s.tangent_vector(1,(0,0),(1,0)))
        sage: traj
        Straight line trajectory made of 1 segments from (0, 0) in polygon 1 to (1, 1) in polygon 2
        sage: traj.is_saddle_connection()
        True
        sage: traj2 = StraightLineTrajectory(s.tangent_vector(1,(0,0),(0,1)))
        sage: traj2
        Straight line trajectory made of 1 segments from (1, 0) in polygon 2 to (0, 1) in polygon 1
        sage: traj2.is_saddle_connection()
        True
    """
    def __init__(self, tangent_vector):
        self._segments = deque()
        seg = SegmentInPolygon(tangent_vector)
        self._segments.append(seg)
        self._setup_forward()
        self._setup_backward()
        self._s=tangent_vector.surface()

    def surface(self):
        return self._s

    def segment(self, i):
        r"""
        EXAMPLES::

            sage: from flatsurf import *

            sage: O = translation_surfaces.regular_octagon()
            sage: v = O.tangent_vector(0, (1,1), (33,45))
            sage: L = v.straight_line_trajectory()
            sage: L.segment(0)
            Segment in polygon 0 starting at (4/15, 0) and ending at (11/26*a +
            1, 15/26*a + 1)
            sage: L.flow(-1)
            sage: L.segment(0)
            Segment in polygon 0 starting at (-1/2*a, 7/22*a + 7/11) and ending
            at (4/15, a + 1)
            sage: L.flow(1)
            sage: L.segment(2)
            Segment in polygon 0 starting at (-1/13*a, 1/13*a) and ending at
            (9/26*a + 11/13, 17/26*a + 15/13)
        """
        return self.segments()[i]

    def combinatorial_length(self):
        return len(self.segments())

    def segments(self):
        return self._segments

    def _setup_forward(self):
        v = self.terminal_tangent_vector()
        if v.is_based_at_singularity():
            self._forward = None
        else:
            self._forward = v.invert()

    def _setup_backward(self):
        v = self.initial_tangent_vector()
        if v.is_based_at_singularity():
            self._backward = None
        else:
            self._backward = v.invert()

    def is_forward_separatrix(self):
        return self._forward is None

    def is_backward_separatrix(self):
        return self._backward is None

    def is_saddle_connection(self):
        return (self._forward is None) and (self._backward is None)

    def is_closed(self):
        r"""
        Test whether this is a closed trajectory.

        By convention, by a closed trajectory we mean a trajectory without any
        singularities.

        .. SEEALSO::

            :meth:`is_saddle_connection`

        EXAMPLES:

        An example in a cone surface covered by the torus::

            sage: from flatsurf import *
            sage: p = polygons.square()
            sage: s = Surface_list(base_ring=p.base_ring())
            sage: s.add_polygon(p,[(0,3),(0,2),(0,1),(0,0)])
            0
            sage: s.set_immutable()
            sage: t = RationalConeSurface(s)

            sage: v = t.tangent_vector(0, (1/2,0), (1/3,7/5))
            sage: l = v.straight_line_trajectory()
            sage: l.is_closed()
            False
            sage: l.flow(100)
            sage: l.is_closed()
            True

            sage: v = t.tangent_vector(0, (1/2,0), (1/3,2/5))
            sage: l = v.straight_line_trajectory()
            sage: l.flow(100)
            sage: l.is_closed()
            False
            sage: l.is_saddle_connection()
            False
            sage: l.flow(-100)
            sage: l.is_saddle_connection()
            True
        """
        return (not self.is_forward_separatrix()) and \
            self._forward.differs_by_scaling(self.initial_tangent_vector())

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


class StraightLineTrajectoryTranslation(AbstractStraightLineTrajectory):
    r"""
    Straight line trajectory in a translation surface.

    This is similar to :class:`StraightLineTrajectory` but implemented using
    interval exchange maps. It should be faster than the implementation via
    segments and flowing in polygons.

    This class only stores a list of triples ``(p, e, x)`` where:

    - ``p`` is a label of a polygon

    - ``e`` is the number of some edge in ``p``

    - ``x`` is the position of the point in ``e`` (be careful that it is not
      necessarily a number between 0 and 1. It is given relatively to the length
      of the induced interval in the iet)

    (see the methods :meth:`_prev` and :meth:`_next`)
    """
    def __init__(self, tangent_vector):
        t = tangent_vector.polygon_label()
        self._vector = tangent_vector.vector()
        self._s = tangent_vector.surface()

        seg = SegmentInPolygon(tangent_vector)
        if seg.is_edge():
            self._points = None
            self._edge = seg
            return

        start = seg.start()
        pos = start._position
        if pos._position_type == pos.EDGE_INTERIOR:
            i = pos.get_edge()
        elif pos._position_type == pos.VERTEX:
            i = pos.get_vertex()
        else:
            raise RuntimeError("PROBLEM!")

        p = start.polygon_label()
        poly = self._s.polygon(p)

        T = self._get_iet(p)
        x = get_linearity_coeff(poly.vertex(i+1) - poly.vertex(i),
                                start.point() - poly.vertex(i))
        x *= T.length_bot(i)

        self._points = deque() # we store triples (lab, edge, rel_pos)
        self._points.append((p, i, x))

    def _next(self, p, e, x):
        r"""
        Return the image of ``(p, e, x)``

        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectoryTranslation
            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: v = o.tangent_vector(1, (1/3,1/7), (5,13))
            sage: L = StraightLineTrajectoryTranslation(v)
            sage: t0 = (1,0,1/3)
            sage: t1 = L._next(*t0)
            sage: t2 = L._next(*t1)
            sage: t0,t1,t2
            ((1, 0, 1/3), (3, 0, 16/3), (1, 0, 31/3))
            sage: assert L._previous(*t2) == t1
            sage: assert L._previous(*t1) == t0
        """
        e, x = self._get_iet(p).forward_image(e, x)
        p, e = self._s.opposite_edge(p, e)
        return (p, e, x)

    def _previous(self, p, e, x):
        r"""
        Return the preimage of ``(p, e, x)``
        """
        p, e = self._s.opposite_edge(p, e)
        e, x = self._get_iet(p).backward_image(e, x)
        return (p, e, x)

    def combinatorial_length(self):
        if self._points is None:
            return 1
        return len(self._points)

    def _get_iet(self, label):
        polygon = self._s.polygon(label)
        try:
            return self._iets[polygon]
        except AttributeError:
            self._iets = {polygon: polygon.flow_map(self._vector)}
        except KeyError:
            self._iets[polygon] = polygon.flow_map(self._vector)
        return self._iets[polygon]

    def segment(self, i):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectoryTranslation

            sage: O = translation_surfaces.regular_octagon()
            sage: v = O.tangent_vector(0, (1,1), (33,45))
            sage: L = StraightLineTrajectoryTranslation(v)
            sage: L.segment(0)
            Segment in polygon 0 starting at (4/15, 0) and ending at (11/26*a +
            1, 15/26*a + 1)
            sage: L.flow(-1)
            sage: L.segment(0)
            Segment in polygon 0 starting at (-1/2*a, 7/22*a + 7/11) and ending
            at (4/15, a + 1)
            sage: L.flow(1)
            sage: L.segment(2)
            Segment in polygon 0 starting at (-1/13*a, 1/13*a) and ending at
            (9/26*a + 11/13, 17/26*a + 15/13)
        """
        if self._points is None:
            return self._edge
        lab, e0, x0 = self._points[i]
        iet = self._get_iet(lab)
        e1, x1 = iet.forward_image(e0, x0)
        poly = self._s.polygon(lab)

        l0 = iet.length_bot(e0)
        l1 = iet.length_top(e1)

        point0 = poly.vertex(e0) + poly.edge(e0) * x0/l0
        point1 = poly.vertex(e1) + poly.edge(e1) * (l1-x1)/l1
        v0 = self._s.tangent_vector(lab, point0, self._vector, ring=self._vector.base_ring())
        v1 = self._s.tangent_vector(lab, point1, -self._vector, ring=self._vector.base_ring())
        return SegmentInPolygon(v0,v1)

    def segments(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectoryTranslation

            sage: s = translation_surfaces.square_torus()
            sage: v = s.tangent_vector(0, (0,0), (1,1+AA(5).sqrt()), ring=AA)
            sage: L = StraightLineTrajectoryTranslation(v)
            sage: L.flow(2)
            sage: L.segments()
            [Segment in polygon 0 starting at (0, 0) and ending at (0.3090169943749474?, 1),
             Segment in polygon 0 starting at (0.3090169943749474?, 0) and ending at (0.618033988749895?, 1),
             Segment in polygon 0 starting at (0.618033988749895?, 0) and ending at (0.9270509831248423?, 1)]
        """
        return [self.segment(i) for i in range(self.combinatorial_length())]

    def is_closed(self):
        if self._points is None:
            raise NotImplementedError
        return self._points[0] == self._next(*self._points[-1])

    def is_forward_separatrix(self):
        if self._points is None:
            return True
        p1,e1,x1 = self._next(*self._points[-1])
        return x1.is_zero()

    def is_backward_separatrix(self):
        return self._points is None or self._points[0][2].is_zero()

    def is_saddle_connection(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectoryTranslation

            sage: torus = translation_surfaces.square_torus()
            sage: v = torus.tangent_vector(0, (1/2,1/2), (1,1))
            sage: S = StraightLineTrajectoryTranslation(v)
            sage: S.is_saddle_connection()
            True

            sage: v = torus.tangent_vector(0, (1/3,2/3), (1,2))
            sage: S = StraightLineTrajectoryTranslation(v)
            sage: S.is_saddle_connection()
            False
            sage: S.flow(1)
            sage: S.is_saddle_connection()
            True
        """
        return self._points is None or (self.is_forward_separatrix() and self.is_backward_separatrix())

    def flow(self, steps):
        if self._points is None:
            return
        if steps > 0:
            t = self._points[-1]
            for i in range(steps):
                t = self._next(*t)
                if t == self._points[0] or t[2].is_zero():
                    break
                self._points.append(t)
        elif steps < 0:
            t = self._points[0]
            for i in range(-steps):
                if t[2].is_zero():
                    break
                t = self._previous(*t)
                if t == self._points[-1]:
                    # closed curve or backward separatrix
                    break
                self._points.appendleft(t)
