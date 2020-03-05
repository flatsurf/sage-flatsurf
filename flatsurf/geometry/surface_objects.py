r"""
Geometric objects on surfaces.

This includes singularities, saddle connections and cylinders.
"""

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

from sage.misc.cachefunc import cached_method
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.plot.graphics import Graphics
from sage.plot.polygon import polygon2d
from sage.rings.qqbar import AA
from sage.structure.sage_object import SageObject

from .polygon import dot_product, ConvexPolygons, wedge_product
from .similarity import SimilarityGroup

class Singularity(SageObject):
    r"""
    Represents a combinatorial singularity on a surface.

    Such a combinatorial singularity is an equivalence class of vertices of the polygons
    making up the surface. This is the coarsest equivalence relation where two vertices
    are equivalent if they are glued along an edge.

    EXAMPLES::

        sage: from flatsurf.geometry.similarity_surface_generators import TranslationSurfaceGenerators
        sage: s=TranslationSurfaceGenerators.veech_2n_gon(5)
        sage: from flatsurf.geometry.surface_objects import Singularity
        sage: sing=Singularity(s,0,1)
        sage: print(sing)
        singularity with vertex equivalence class frozenset(...)
        sage: TestSuite(sing).run()
    """

    def __init__(self, similarity_surface, l, v, limit=None):
        r"""
        Represents the singularity associated to the v-th vertex of the polygon with
        label l.

        If the surface is infinite, the limit needs to be set. In this case the construction
        of the singularity is successful if the sequence of vertices hit by passing through
        edges closes up in limit or less steps.
        """
        from .similarity_surface import SimilaritySurface
        self._ss=similarity_surface
        self._s=set()
        if not self._ss.is_finite() and limit is None:
            raise ValueError("need a limit when working with an infinite surface")
        start=(l,v)
        self._s.add(start)
        edge=self._ss.opposite_edge(l,v)
        next = (edge[0], (edge[1]+1)%self._ss.polygon(edge[0]).num_edges() )
        while start!=next:
            self._s.add(next)
            if not limit is None and len(self._s)>limit:
                raise ValueError("Number of vertices in singularities exceeds limit.")
            edge=self._ss.opposite_edge(next)
            next = (edge[0], (edge[1]+1)%self._ss.polygon(edge[0]).num_edges() )
        self._s=frozenset(self._s)

    def surface(self):
        r"""
        Return the SimilaritySurface where the singularity appears.
        """
        return self._ss

    def one_vertex(self):
        r"""
        Return a pair (l,v) from the equivalence class of this singularity.
        """
        return next(iter(self._s))

    def vertex_set(self):
        r"""
        Return the set of pairs (l,v) in the equivalence class of this singularity.
        """
        return self._s

    def contains_vertex(self, l, v=None):
        r"""
        Checks if the pair (l,v) is in the equivalence class returning true or false.

        If v is None, then check if the pair l is in the equivalence class.
        """
        if v is None:
            return l in self._s
        else:
            return (l,v) in self._s

    def _repr_(self):
        return "singularity with vertex equivalence class "+repr(self._s)

    def __eq__(self,other):
        if self is other:
            return True
        if not isinstance(other, Singularity):
            raise TypeError
        if not self._ss==other._ss:
            return False
        return self._s == other._s

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        # Hash only using the set of vertices (rather than including the surface)
        return hash(self._s)

class SurfacePoint(SageObject):
    r"""
    Represents a point on a SimilaritySurface.

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.surface_objects import SurfacePoint
        sage: p = SymmetricGroup(2)('(1,2)')
        sage: s = translation_surfaces.origami(p,p)
        sage: SurfacePoint(s, 1, (1/2,1/2))
        Surface point located at (1/2, 1/2) in polygon 1
        sage: sp1 = SurfacePoint(s, 1, (1/2,0))
        sage: sp1
        Surface point with 2 coordinate representations
        sage: sp2 = SurfacePoint(s, 2, (1/2,1))
        sage: sp1 == sp2
        True
        sage: hash(sp1) == hash(sp2)
        True
        sage: sp1.coordinates(2)
        frozenset({(1/2, 1)})
        sage: sp = SurfacePoint(s, 1, (0,0))
        sage: sp
        Surface point with 4 coordinate representations
    """
    def __init__(self, surface, label, point, ring = None, limit=None):
        self._s = surface
        if ring is None:
            self._ring = surface.base_ring()
        else:
            self._ring = ring
        p = surface.polygon(label)
        point = VectorSpace(self._ring,2)(point)
        point.set_immutable()
        pos = p.get_point_position(point)
        assert pos.is_inside(), \
            "Point must be positioned within the polygon with the given label."
        # This is the correct thing if point lies in the interior of the polygon with the given label.
        self._coordinate_dict = {label: {point}}
        if pos.is_in_edge_interior():
            label2,e2 = surface.opposite_edge(label, pos.get_edge())
            point2 = surface.edge_transformation(label, pos.get_edge())(point)
            point2.set_immutable()
            if label2 in self._coordinate_dict:
                self._coordinate_dict[label2].add(point2)
            else:
                self._coordinate_dict[label2]={point2}
        if pos.is_vertex():
            self._coordinate_dict = {}
            sing = surface.singularity(label, pos.get_vertex(), limit=limit)
            for l,v in sing.vertex_set():
                new_point = surface.polygon(l).vertex(v)
                new_point.set_immutable()
                if l in self._coordinate_dict:
                    self._coordinate_dict[l].add(new_point)
                else:
                    self._coordinate_dict[l] = {new_point}
        # Freeze the sets.
        for label,point_set in iteritems(self._coordinate_dict):
            self._coordinate_dict[label] = frozenset(point_set)

    def surface(self):
        r"""
        Return the surface containing this point.
        """
        return self._s

    def num_coordinates(self):
        r"""
        Return the number of different coordinate representations of the point.
        """
        try:
            return self._num_coordinates
        except AttributeError:
            count=0
            for label,point_set in iteritems(self._coordinate_dict):
                count += len(point_set)
            self._num_coordinates = count
            return count
        else:
            raise ValueError("Unable to return num_coordinates()")

    def labels(self):
        r"""
        Return the list of labels of polygons containing the point in its closure.

        This will be a list of one label if the point is the the interior of a polygon,
        at most two labels if it is on the interior of an edge, and can be lots of labels
        if the point is a singularity.
        """
        return list(self._coordinate_dict.keys())

    def coordinates(self, label):
        r"""
        Return a frozenset of coordinates for the closure of the point in the polygon
        with the provided label.

        The set will consist of one point if the point lies in the interior of a polygon,
        will be two points if the point lies in the interior of two edges of the polygon
        simultaneously and can be lots of points if the point is a singularity.
        """
        return self._coordinate_dict[label]

    def graphical_surface_point(self, graphical_surface = None):
        r"""
        Return the GraphicalSurfacePoint built from this SurfacePoint.
        """
        from flatsurf.graphical.surface_point import GraphicalSurfacePoint
        return GraphicalSurfacePoint(self, graphical_surface = graphical_surface)

    def plot(self, *args, **options):
        r"""
        Plot this point. There may be one argument which provides a graphical surface.
        All options are passed two the ploting method of GraphicalSurfacePoint.
        """
        if len(args) > 1:
            raise ValueError("SurfacePoint.plot() can take at most one argument.")
        if len(args) == 1:
            return self.graphical_surface_point(graphical_surface=args[0]).plot(**options)
        else:
            return self.graphical_surface_point().plot(**options)

    def __repr__(self):
        r"""
        TESTS::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1,1,1,1], [1, 1/2, 1/3, 1/4])
            sage: S.point(0, (1/2,1/2))
            Surface point located at (1/2, 1/2) in polygon 0
        """
        if self.num_coordinates() == 1:
            return "Surface point located at {} in polygon {}".format(
                next(iter(self.coordinates(self.labels()[0]))),self.labels()[0])
        else:
            return "Surface point with {} coordinate representations".format(
                self.num_coordinates())

    def __eq__(self,other):
        if self is other:
            return True
        if not isinstance(other, SurfacePoint):
            raise TypeError("Comparing SurfacePoint to an object of different type")
        if not self._s == other._s:
            raise ValueError("Comparing SurfacePoints on different surfaces")
        return self._coordinate_dict == other._coordinate_dict

    def __hash__(self):
        h=0
        for label,point_set in iteritems(self._coordinate_dict):
            h += 677*hash(label)+hash(point_set)
        return h

    def __ne__(self, other):
        return not self == other


class SaddleConnection(SageObject):
    r"""
    Represents a saddle connection on a SimilaritySurface.
    """

    def __init__(self, surface, start_data, direction,
            end_data=None, end_direction=None,
            holonomy=None, end_holonomy=None,
            check=True, limit=1000):
        r"""
        Construct a saddle connecton on a SimilaritySurface.

        The only necessary parameters are the surface, start_data, and direction
        (to start). If there is missing data that can not be inferred from the surface
        type, then a straight-line trajectory will be computed to confirm that this is
        indeed a saddle connection. The trajectory will pass through at most limit
        polygons before we give up.

        Details of the parameters are provided below.

        Parameters
        ----------
        surface : a SimilaritySurface
            which will contain the saddle connection being constructed.
        start_data : a pair
            consisting of the label of the polygon where the saddle connection starts
            and the starting vertex.
        direction : 2-dimensional vector with entries in the base_ring of the surface
            representing the direction the saddle connection is moving in (in the
            coordinates of the initial polygon).
        end_data : a pair
            consisting of the label of the polygon where the saddle connection terminates
            and the terminating vertex.
        end_direction : 2-dimensional vector with entries in the base_ring of the surface
            representing the direction to move backward from the end point (in the
            coordinates of the terminal polygon). If the surface is a DilationSurface
            or better this will be the negation of the direction vector. If the surface
            is a HalfDilation surface or better, then this will be either the direction
            vector or its negation. In either case the value can be inferred from the
            end_data.
        holonomy : 2-dimensional vector with entries in the base_ring of the surface
            the holonomy of the saddle connection measured from the start. To compute this
            you develop the saddle connection into the plane starting from the starting
            polygon.
        end_holonomy : 2-dimensional vector with entries in the base_ring of the surface
            the holonomy of the saddle connection measured from the end (with the opposite
            orientation). To compute this you develop the saddle connection into the plane
            starting from the terminating polygon. For a translation surface, this will be
            the negation of holonomy, and for a HalfTranslation surface it will be either
            equal to holonomy or equal to its negation. In both these cases the end_holonomy
            can be inferred and does not need to be passed to the constructor.
        check : boolean
            If all data above is provided or can be inferred, then when check=False this
            geometric data is not verified. With check=True the data is always verified
            by straight-line flow. Erroroneous data will result in a ValueError being thrown.
            Defaults to true.
        limit :
            The combinatorial limit (in terms of number of polygons crossed) to flow forward
            to check the saddle connection geometry.
        """
        from .similarity_surface import SimilaritySurface
        assert isinstance(surface,SimilaritySurface)
        self._s=surface

        # Sanitize the direction vector:
        V=self._s.vector_space()
        self._direction=V(direction)
        if self._direction==V.zero():
            raise ValueError("Direction must be nonzero.")
        # To canonicalize the direction vector we ensure its endpoint lies in the boundary of the unit square.
        xabs=self._direction[0].abs()
        yabs=self._direction[1].abs()
        if xabs>yabs:
            self._direction=self._direction/xabs
        else:
            self._direction=self._direction/yabs

        # Fix end_direction if not standard.
        if end_direction is not None:
            xabs=end_direction[0].abs()
            yabs=end_direction[1].abs()
            if xabs>yabs:
                end_direction=end_direction/xabs
            else:
                end_direction=end_direction/yabs

        self._start_data=tuple(start_data)

        if end_direction is None:
            from .half_dilation_surface import HalfDilationSurface
            from .dilation_surface import DilationSurface
            # Attempt to infer the end_direction.
            if isinstance(self._s,DilationSurface):
                end_direction=-self._direction
            elif isinstance(self._s,HalfDilationSurface) and end_data is not None:
                p=self._s.polygon(end_data[0])
                if wedge_product(p.edge(end_data[1]), self._direction)>=0 and \
                   wedge_product(p.edge( (p.num_edges()+end_data[1]-1)%p.num_edges() ), self._direction)>0:
                    end_direction=self._direction
                else:
                    end_direction=-self._direction

        if end_holonomy is None and holonomy is not None:
            # Attempt to infer the end_holonomy:
            from .half_translation_surface import HalfTranslationSurface
            from .translation_surface import TranslationSurface
            if isinstance(self._s,TranslationSurface):
                end_holonomy=-holonomy
            if isinstance(self._s,HalfTranslationSurface):
                if direction==end_direction:
                    end_holonomy=holonomy
                else:
                    end_holonomy=-holonomy

        if  end_data is None or end_direction is None or holonomy is None or end_holonomy is None or check:
            v=self.start_tangent_vector()
            traj=v.straight_line_trajectory()
            traj.flow(limit)
            if not traj.is_saddle_connection():
                raise ValueError("Did not obtain saddle connection by flowing forward. Limit="+str(limit))
            tv=traj.terminal_tangent_vector()
            self._end_data=(tv.polygon_label(), tv.vertex())
            if end_data is not None:
                if end_data!=self._end_data:
                    raise ValueError("Provided or inferred end_data="+str(end_data)+" does not match actual end_data="+str(self._end_data))
            self._end_direction=tv.vector()
            # Canonicalize again.
            xabs=self._end_direction[0].abs()
            yabs=self._end_direction[1].abs()
            if xabs>yabs:
                self._end_direction = self._end_direction / xabs
            else:
                self._end_direction = self._end_direction / yabs
            if end_direction is not None:
                if end_direction!=self._end_direction:
                    raise ValueError("Provided or inferred end_direction="+str(end_direction)+" does not match actual end_direction="+str(self._end_direction))

            if traj.segments()[0].is_edge():
                # Special case (The below method causes error if the trajectory is just an edge.)
                self._holonomy = self._s.polygon(start_data[0]).edge(start_data[1])
                self._end_holonomy = self._s.polygon(self._end_data[0]).edge(self._end_data[1])
            else:
                from .similarity import SimilarityGroup
                sim=SimilarityGroup(self._s.base_ring()).one()
                itersegs = iter(traj.segments())
                next(itersegs)
                for seg in itersegs:
                    sim = sim * self._s.edge_transformation(seg.start().polygon_label(),
                                                            seg.start().position().get_edge())
                self._holonomy = sim(traj.segments()[-1].end().point())- \
                    traj.initial_tangent_vector().point()
                self._end_holonomy = -( (~sim.derivative())*self._holonomy )

            if holonomy is not None:
                if holonomy!=self._holonomy:
                    print("Combinatorial length: "+str(traj.combinatorial_length()))
                    print("Start: "+str(traj.initial_tangent_vector().point()))
                    print("End: "+str(traj.terminal_tangent_vector().point()))
                    print("Start data:"+str(start_data))
                    print("End data:"+str(end_data))
                    raise ValueError("Provided holonomy "+str(holonomy)+
                                     " does not match computed holonomy of "+str(self._holonomy))
            if end_holonomy is not None:
                if end_holonomy!=self._end_holonomy:
                    raise ValueError("Provided or inferred end_holonomy "+str(end_holonomy)+
                                     " does not match computed end_holonomy of "+str(self._end_holonomy))
        else:
            self._end_data=tuple(end_data)
            self._end_direction=end_direction
            self._holonomy=holonomy
            self._end_holonomy=end_holonomy

        # Make vectors immutable
        self._direction.set_immutable()
        self._end_direction.set_immutable()
        self._holonomy.set_immutable()
        self._end_holonomy.set_immutable()

    def surface(self):
        return self._s

    def direction(self):
        r"""
        Returns a vector parallel to the saddle connection pointing from the start point.

        The will be normalized so that its l_\infty norm is 1.
        """
        return self._direction

    def end_direction(self):
        r"""
        Returns a vector parallel to the saddle connection pointing from the end point.

        The will be normalized so that its l_\infty norm is 1.
        """
        return self._end_direction

    def start_data(self):
        r"""
        Return the pair (l,v) representing the label and vertex of the corresponding polygon
        where the saddle connection originates.
        """
        return self._start_data

    def end_data(self):
        r"""
        Return the pair (l,v) representing the label and vertex of the corresponding polygon
        where the saddle connection terminates.
        """
        return self._end_data

    def holonomy(self):
        r"""
        Return the holonomy vector of the saddle connection (measured from the start).

        In a SimilaritySurface this notion corresponds to developing the saddle connection into the plane
        using the initial chart coming from the initial polygon.
        """
        return self._holonomy

    def length(self):
        r"""
        In a cone surface, return the length of this saddle connection. Since
        this may not lie in the field of definition of the surface, it is
        returned as an element of the Algebraic Real Field.
        """
        from .cone_surface import ConeSurface
        assert isinstance(self._s,ConeSurface), \
            "Length of a saddle connection only makes sense for cone surfaces."
        return vector(AA,self._holonomy).norm()

    def end_holonomy(self):
        r"""
        Return the holonomy vector of the saddle connection (measured from the end).

        In a SimilaritySurface this notion corresponds to developing the saddle connection into the plane
        using the initial chart coming from the initial polygon.
        """
        return self._end_holonomy


    def start_tangent_vector(self):
        r"""
        Return a tangent vector to the saddle connection based at its start.
        """
        return self._s.tangent_vector(self._start_data[0],
                                      self._s.polygon(self._start_data[0]).vertex(self._start_data[1]),
                                      self._direction)

    def trajectory(self, limit = 1000, cache = True):
        r"""
        Return a straight line trajectory representing this saddle connection.
        Fails if the trajectory passes through more than limit polygons.
        """
        try:
            return self._traj
        except AttributeError:
            pass
        v=self.start_tangent_vector()
        traj=v.straight_line_trajectory()
        traj.flow(limit)
        if not traj.is_saddle_connection():
            raise ValueError("Did not obtain saddle connection by flowing forward. Limit="+str(limit))
        if cache:
            self._traj = traj
        return traj

    def plot(self, *args, **options):
        r""" Equivalant to .trajectory().plot(*args, **options)
        """
        return self.trajectory().plot(*args, **options)

    def end_tangent_vector(self):
        r"""
        Return a tangent vector to the saddle connection based at its start.
        """
        return self._s.tangent_vector(self._end_data[0],
                                      self._s.polygon(self._end_data[0]).vertex(self._end_data[1]),
                                      self._end_direction)

    def invert(self):
        r"""
        Return this saddle connection but with opposite orientation.
        """
        return SaddleConnection(self._s,self._end_data, self._end_direction,
           self._start_data, self._direction,
           self._end_holonomy, self._holonomy,
           check=False)

    def intersections(self, traj, count_singularities = False, include_segments = False):
        r"""
        See documentation of :func:`~straight_line_trajectory.AbstractStraightLineTrajectory.intersections`
        """
        return self.trajectory().intersections(traj, count_singularities, include_segments)

    def intersects(self, traj, count_singularities = False):
        r"""
        See documentation of :func:`~straight_line_trajectory.AbstractStraightLineTrajectory.intersects`
        """
        return self.trajectory().intersects(traj, count_singularities=count_singularities)

    def __eq__(self,other):
        if self is other:
            return True
        if not isinstance(other, SaddleConnection):
            raise TypeError
        if not self._s==other._s:
            return False
        if not self._direction == other._direction:
            return False
        if not self._start_data == other._start_data:
            return False
        # Initial data should determine the saddle connection:
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return 41*hash(self._direction)-97*hash(self._start_data)

    def _test_geometry(self, **options):
        # Test that this saddle connection actually exists on the surface.
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        sc=SaddleConnection(self._s,self._start_data, self._direction,
                           self._end_data, self._end_direction,
                           self._holonomy, self._end_holonomy,
                           check=True)

    def __repr__(self):
        return "Saddle connection in direction {} with start data {} and end data {}".format(
            self._direction, self._start_data, self._end_data)

    def _test_inverse(self, **options):
        # Test that inverting works properly.
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        SaddleConnection(self._s,self._end_data, self._end_direction,
           self._start_data, self._direction,
           self._end_holonomy, self._holonomy,
           check=True)

class Cylinder(SageObject):
    r"""
    Represents a cylinder in a SimilaritySurface. A cylinder for these purposes is a
    topological annulus in a surface bounded by a finite collection of saddle connections
    meeting at 180 degree angles.

    To Do
    -----
    * Support cylinders whose monodromy is a dilation.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s = translation_surfaces.octagon_and_squares()
        sage: from flatsurf.geometry.surface_objects import Cylinder
        sage: cyl = Cylinder(s, 0, [2, 3, 3, 3, 2, 0, 1, 3, 2, 0])
        sage: cyl.initial_label()
        0
        sage: cyl.edges()
        (2, 3, 3, 3, 2, 0, 1, 3, 2, 0)
        sage: # a = sqrt(2) below.
        sage: cyl.area()
        2*a + 4
        sage: cyl.circumference().minpoly()
        x^4 - 680*x^2 + 400
        sage: cyl.holonomy()
        (8*a + 12, 4*a + 6)
    """
    def __init__(self, s, label0, edges):
        r"""
        Construct a cylinder on the surface `s` from an initial label and a
        sequence of edges crossed.

        Parameters
        ----------
        s: A SimilaritySurface
            the surface conaining the cylinder
        label0: An initial label
            representing a polygon the cylinder passes through.
        edges: a list
            giving the sequence of edges the cylinder crosses until it closes.
        """
        self._s = s
        self._label0 = label0
        self._edges = tuple(edges)
        ss = s.minimal_cover(cover_type="planar")
        SG = SimilarityGroup(s.base_ring())
        labels = [(label0, SG.one())] # labels of polygons on the cover ss.
        for e in edges:
            labels.append(ss.opposite_edge(labels[-1],e)[0])
        if labels[0][0] != labels[-1][0]:
            raise ValueError("Combinatorial path does not close.")
        trans = labels[-1][1]
        if not trans.is_translation():
            raise NotImplemented("Only cylinders with translational monodromy are currently supported")
        m = trans.matrix()
        v = vector(s.base_ring(),(m[0][2],m[1][2])) # translation vector
        from flatsurf.geometry.polygon import wedge_product

        p = ss.polygon(labels[0])
        e = edges[0]
        min_y = wedge_product(v, p.vertex(e))
        max_y = wedge_product(v, p.vertex((e+1)%p.num_edges()))
        if min_y >= max_y:
            raise ValueError("Combinatorial data does not represent a cylinder")

        # Stores the vertices where saddle connections starts:
        min_list = [0]
        max_list = [0]

        for i in range(1, len(edges)):
            e = edges[i]
            p = ss.polygon(labels[i])
            y = wedge_product(v, p.vertex(e))
            if y == min_y:
                min_list.append(i)
            elif y > min_y:
                min_list = [i]
                min_y = y
                if min_y >= max_y:
                    raise ValueError("Combinatorial data does not represent a cylinder")
            y = wedge_product(v, p.vertex((e+1)%p.num_edges()))
            if y == max_y:
                max_list.append(i)
            elif y < max_y:
                max_list = [i]
                max_y = y
                if min_y >= max_y:
                    raise ValueError("Combinatorial data does not represent a cylinder")

        # Extract the saddle connections on the right side:
        from flatsurf.geometry.surface_objects import SaddleConnection
        sc_set_right = set()
        vertices = []
        for i in min_list:
            l = labels[i]
            p = ss.polygon(l)
            vertices.append((i,p.vertex(edges[i])))
        i,vert_i = vertices[-1]
        vert_i = vert_i - v
        j,vert_j = vertices[0]
        if vert_i != vert_j:
            li = labels[i]
            li = (li[0], SG(-v)*li[1])
            lio = ss.opposite_edge(li,edges[i])
            lj = labels[j]
            sc = SaddleConnection(s,
                                  (lio[0][0], (lio[1]+1) % ss.polygon(lio[0]).num_edges()),
                                  (~lio[0][1])(vert_j)-(~lio[0][1])(vert_i))
            sc_set_right.add(sc)
        i = j
        vert_i = vert_j
        for j,vert_j in vertices[1:]:
            if vert_i != vert_j:
                li = labels[i]
                li = (li[0], SG(-v)*li[1])
                lio = ss.opposite_edge(li,edges[i])
                lj = labels[j]
                sc = SaddleConnection(s,
                                      (lio[0][0], (lio[1]+1) % ss.polygon(lio[0]).num_edges()),
                                      (~lio[0][1])(vert_j)-(~lio[0][1])(vert_i),
                                      limit = j-i)
                sc_set_right.add(sc)
            i = j
            vert_i =vert_j

        # Extract the saddle connections on the left side:
        sc_set_left = set()
        vertices = []
        for i in max_list:
            l = labels[i]
            p = ss.polygon(l)
            vertices.append((i,p.vertex((edges[i]+1)%p.num_edges())))
        i,vert_i = vertices[-1]
        vert_i = vert_i - v
        j,vert_j = vertices[0]
        if vert_i != vert_j:
            li = labels[i]
            li = (li[0], SG(-v)*li[1])
            lio = ss.opposite_edge(li,edges[i])
            lj = labels[j]
            sc = SaddleConnection(s,
                                  (lj[0], (edges[j]+1) % ss.polygon(lj).num_edges()),
                                  (~lj[1])(vert_i)-(~lj[1])(vert_j))
            sc_set_left.add(sc)
        i = j
        vert_i =vert_j
        for j,vert_j in vertices[1:]:
            if vert_i != vert_j:
                li = labels[i]
                lio = ss.opposite_edge(li,edges[i])
                lj = labels[j]
                sc = SaddleConnection(s,
                                      (lj[0], (edges[j]+1) % ss.polygon(lj).num_edges()),
                                      (~lj[1])(vert_i)-(~lj[1])(vert_j))
                sc_set_left.add(sc)
            i = j
            vert_i =vert_j
        self._boundary1 = frozenset(sc_set_right)
        self._boundary2 = frozenset(sc_set_left)
        self._boundary = frozenset(self._boundary1.union(self._boundary2))

        edge_intersections = []
        i = min_list[0]
        l = labels[i]
        p = ss.polygon(l)
        right_point = p.vertex(edges[i]) # point on the right boundary
        i = max_list[0]
        l = labels[i]
        p = ss.polygon(l)
        left_point = p.vertex((edges[i]+1)%p.num_edges())
        from flatsurf.geometry.polygon import solve
        for i in range(len(edges)):
            l = labels[i]
            p = ss.polygon(l)
            e = edges[i]
            v1 = p.vertex(e)
            v2 = p.vertex((e+1)%p.num_edges())
            a,b = solve(left_point, v, v1, v2-v1)
            w1 = (~(l[1]))(v1 + b*(v2-v1))
            a,b = solve(right_point, v, v1, v2-v1)
            w2 = (~(l[1]))(v1 + b*(v2-v1))
            edge_intersections.append((w1,w2))

        polygons = []
        P = ConvexPolygons(s.base_ring())
        pair1 = edge_intersections[-1]
        l1 = labels[-2][0]
        e1 = edges[-1]
        for i in range(len(edges)):
            l2 = labels[i][0]
            pair2 = edge_intersections[i]
            e2 = edges[i]
            trans = s.edge_transformation(l1,e1)
            pair1p = (trans(pair1[0]), trans(pair1[1]))
            polygon_verts = [pair1p[0], pair1p[1]]
            if pair2[1] != pair1p[1]:
                polygon_verts.append(pair2[1])
            if pair2[0] != pair1p[0]:
                polygon_verts.append(pair2[0])
            polygons.append((l2,P(vertices=polygon_verts)))
            l1 = l2
            pair1 = pair2
            e1 = e2
        self._polygons = tuple(polygons)

    def surface(self):
        return self._s

    def initial_label(self):
        r"""
        Return one label on the surface that the cylinder passes through.
        """
        return self._label0

    def edges(self):
        r"""
        Return a tuple of edge numbers representing the edges crossed
        when the cylinder leaves the polygon with `initial_label` until
        it returns by closing.
        """
        return self._edges

    def boundary(self):
        r"""
        Return the set of saddle connections in the boundary, oriented so that
        the surface is on the left.
        """
        return self._boundary

    def polygons(self):
        r"""
        Return a list of pairs each consisting of a label and a polygon.
        Each polygon represents a sub-polygon of the polygon on the surface
        with the given label. The union of these sub-polygons form the
        cylinder. The subpolygons are listed in cyclic order.
        """
        return self._polygons

    @cached_method
    def area(self):
        r"""
        Return the area of this cylinder if it is contained in a ConeSurface.
        """
        from .cone_surface import ConeSurface
        assert isinstance(self._s,ConeSurface), \
            "Area only makes sense for cone surfaces."
        area = 0
        for l,p in self.polygons():
            area += p.area()
        return area

    def plot(self, **options):
        r"""
        Plot this cylinder in coordinates used by a graphical surface. This
        plots this cylinder as a union of subpolygons. Only the intersections
        with polygons visible in the graphical surface are shown.

        Parameters other than `graphical_surface` are passed to `polygon2d`
        which is called to render the polygons.

        Parameters
        ----------
        graphical_surface : a GraphicalSurface
            If not provided or `None`, the plot method uses the default graphical
            surface for the surface.
        """
        if "graphical_surface" in options and options["graphical_surface"] is not None:
            gs = options["graphical_surface"]
            assert gs.get_surface() == self._s, "Graphical surface for the wrong surface."
            del options["graphical_surface"]
        else:
            gs = self._s.graphical_surface()
        plt = Graphics()
        for l,p in self.polygons():
            if gs.is_visible(l):
                gp = gs.graphical_polygon(l)
                t = gp.transformation()
                pp = t(p)
                poly = polygon2d(pp.vertices(), **options)
                plt += poly.plot()
        return plt

    @cached_method
    def labels(self):
        r"""
        Return the set of labels that this cylinder passes through.
        """
        polygons = self.polygons()
        return frozenset([l for l,p in polygons])

    def boundary_components(self):
        r"""
        Return a set of two elements: the set of saddle connections on
        the right and left sides. Saddle connections are oriented so that
        the surface is on the left.
        """
        return frozenset([self._boundary1,self._boundary2])

    def next(self, sc):
        r"""
        Return the next saddle connection as you move around the cylinder boundary
        moving from sc in the direction of its orientation.
        """
        assert sc in self._boundary
        v=sc.end_tangent_vector()
        v=v.clockwise_to(-v.vector())
        from flatsurf.geometry.polygon import is_same_direction
        for sc2 in self._boundary:
            if sc2.start_data()==(v.polygon_label(),v.vertex()) and \
                    is_same_direction(sc2.direction(), v.vector()):
                return sc2
        raise ValuError("Failed to find next saddle connection in boundary set.")

    def previous(self,sc):
        r"""
        Return the previous saddle connection as you move around the cylinder boundary
        moving from sc in the direction opposite its orientation.
        """
        assert sc in self._boundary
        v=sc.start_tangent_vector()
        v=v.counterclockwise_to(-v.vector())
        from flatsurf.geometry.polygon import is_same_direction
        for sc2 in self._boundary:
            if sc2.end_data()==(v.polygon_label(),v.vertex()) and \
                    is_same_direction(sc2.end_direction(), v.vector()):
                return sc2
        raise ValuError("Failed to find previous saddle connection in boundary set.")

    @cached_method
    def holonomy(self):
        r"""
        In a translation surface, return one of the two holonomy vectors of the cylinder,
        which differ by a sign.
        """
        from .translation_surface import TranslationSurface
        assert isinstance(self._s,TranslationSurface), \
            "Holonomy currently only computable for translation surfaces."
        V=self._s.vector_space()
        total=V.zero()
        for sc in self._boundary1:
            total += sc.holonomy()

        # Debugging:
        total2=V.zero()
        for sc in self._boundary2:
            total2 += sc.holonomy()
        assert total+total2==V.zero(), "Holonomy of the two boundary components should sum to zero."

        return total

    @cached_method
    def circumference(self):
        r"""
        In a cone surface, return the circumference, i.e., the length
        of a geodesic loop running around the cylinder. Since this may
        not lie in the field of definition of the surface, it is returned
        as an element of the Algebraic Real Field.
        """
        from .cone_surface import ConeSurface
        assert isinstance(self._s,ConeSurface), \
            "Circumference only makes sense for cone surfaces."
        total = 0
        for sc in self._boundary1:
            total += sc.length()
        return total

    #def width_vector(self):
    #    r"""
    #    In a translation surface, return a vector orthogonal to the holonomy vector which cuts
    #    across the cylinder.
    #    """
    #    from flatsurf.geometry.translation_surface import TranslationSurface
    #    assert isinstance(self._s,TranslationSurface), \
    #        "width_vector currently only computable for translation surfaces."
    #    w=self._across.holonomy()
    #    h=iter(self._boundary1).next().holonomy()
    #    from flatsurf.geometry.polygon import dot_product
    #    return w-(dot_product(w,h)/dot_product(h,h))*h
