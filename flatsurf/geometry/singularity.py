r"""
Combinatorial singularities of surfaces.
"""

from __future__ import absolute_import

from sage.structure.sage_object import SageObject


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

