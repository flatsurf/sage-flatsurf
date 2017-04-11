r"""
Geometric objects on surfaces.

This includes singularities, saddle connections and cylinders.
"""

from __future__ import absolute_import

from sage.structure.sage_object import SageObject

from flatsurf import * 
from .polygon import wedge_product, dot_product



class Singularity(SageObject):
    r"""
    Represents a combinatorial singularity on a surface.

    Such a combinatorial singularity is an equivalence class of vertices of the polygons 
    making up the surface. This is the coarsest equivalence relation where two vertices 
    are equivalent if they are glued along an edge.

    EXAMPLES::

        sage: from flatsurf.geometry.similarity_surface_generators import TranslationSurfaceGenerators
        sage: s=TranslationSurfaceGenerators.veech_2n_gon(5)
        sage: from flatsurf.geometry.singularity import Singularity
        sage: sing=Singularity(s,0,1)
        sage: print(sing)
        singularity with vertex equivalence class frozenset([(0, 1), (0, 9), (0, 3), (0, 5), (0, 7)])
        sage: TestSuite(sing).run(verbose=True)
        running ._test_category() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
    """
    
    def __init__(self, similarity_surface, l, v, limit=None):
        r"""
        Represents the singularity associated to the v-th vertex of the polygon with 
        label l.
        
        If the surface is infinite, the limit needs to be set. In this case the construction
        of the singularity is successful if the sequence of vertices hit by passing through
        edges closes up in limit or less steps.
        """
        self._ss=similarity_surface
        self._s=set()
        if not self._ss.is_finite() and limit is None:
            raise ValueError("Need a limit when working with an infinite surface.")
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
    
        self._start_data=tuple(start_data)
        
        if end_direction is None:
            # Attempt to infer the end_direction.
            if isinstance(s,DilationSurface):
                end_direction=-self._direction
            elif isinstance(s,HalfDilationSurface) and end_data is not None:
                p=self._s.polygon(end_data[0])
                if wedge_product(p.edge(end_data[1]), self._direction)>=0 and \
                   wedge_product(p.edge( (p.num_edges()+end_data[1]-1)%p.num_edges() ), self._direction)>0:
                    end_direction=self._direction
                else:
                    end_direction=-self._direction

        if end_holonomy is None and holonomy is not None:
            # Attempt to infer the end_holonomy:
            if isinstance(s,TranslationSurface):
                end_holonomy=-holonomy
            if isinstance(s,HalfTranslationSurface):
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
            self._end_data=(tv.polygon_label(), tv.singularity())
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

            from flatsurf.geometry.similarity import SimilarityGroup
            sim=SimilarityGroup(self._s.base_ring()).one()
            itersegs = iter(traj.segments())
            next(itersegs)
            for seg in itersegs:
                sim = sim * self._s.edge_transformation(seg.start().polygon_label(),
                                                        seg.start().position().get_edge())
            self._holonomy = sim(seg.end().point())-traj.segments()[0].start().point()
            if holonomy is not None:
                if holonomy!=self._holonomy:
                    raise ValueError("Provided holonomy "+str(holonomy)+
                                     " does not match computed holonomy of "+str(self._holonomy))
            self._end_holonomy = -( (~sim.derivative())*self._holonomy )
            if end_holonomy is not None:
                if end_holonomy!=self._end_holonomy:
                    raise ValueError("Provided or inferred end_holonomy "+str(end_holonomy)+
                                     " does not match computed end_holonomy of "+str(self._end_holonomy))
        else:
            self._end_data=tuple(end_data)
            self._end_direction=end_direction
            self._holonomy=holonomy
            self._end_holonomy=end_holonomy
        
    def surface(self):
        return self._s
            
    def direction(self):
        r"""
        Returns a vector paralell to the saddle connection. 
        
        The will be normalized so that its l_\infty norm is 1.
        """
        return self._direction

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
        return 41*hash(self._direction)-97*hash(start_data)
    
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
