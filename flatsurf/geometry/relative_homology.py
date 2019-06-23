r"""
This module contains a lazy implementation of a relative homology, 
$H_1(S,\Sigma; R)$, where $S$ is a similarity surface, $\Sigma$ is the singularities
or vertices, and $R$ is a ring. 

This implementation works for finite or infinite surfaces. For infinite surfaces,
we define relative homology formally. It is simply $R^E$ where $E$ is the edge set
modulo equivalences of two types:
1) If $e$ is an edge, and $e'$ is its opposite edge oriented counterclockwise from the polygon they bound then $e+e'=0$ in homology.
2) The sum of edges around a polygon is zero.
"""

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

from sage.structure.element import ModuleElement
from sage.modules.module import Module
from sage.rings.integer_ring import ZZ

from .similarity_surface import SimilaritySurface


class RelativeHomologyClass(ModuleElement):

    # Implementation notes:
    # self._d will be a dictionary mapping mapping pairs (label, edge) to the base_ring
    # By convention a pair will be in the dictionary only if its image is non-zero.
    def __init__(self, parent, d):
        r"""Do not call directly!"""
        # This should be a dict 
        if not isinstance(d,dict):
            raise ValueError("RelativeHomologyClass.__init__ must be passed a dictionary.")
        self._d = d
        ModuleElement.__init__(self, parent=parent)

    def _rmul_(self, c):
        if c==self.parent().base_ring().zero():
            return self.parent().zero()
        d=dict()
        r=self.parent().base_ring()
        for k,v in iteritems(self._d):
            d[k]=r(c*v)
        return self.parent()._element_from_dict(d)

    def _add_(self, other):
        d=dict()
        r=self.parent().base_ring()
        for k,v in iteritems(self._d):
            if k in other._d:
                total = v + other._d[k]
                if total != self.parent().base_ring().zero():
                    d[k] = r(total)
            else:
                d[k]=r(v)
        for k,v in iteritems(other._d):
            if k not in self._d:
                d[k]=r(v)
        return self.parent()._element_from_dict(d)
    
    def _neg_(self):
        return self._rmul_(-self.parent().base_ring().one())
        
    def __cmp__(self, other):
        # Construct a set of keys
        s=set()
        for k,v in iteritems(self._d):
            s.add(k)
        for k,v in iteritems(other._d):
            s.add(k)
        zero = self.parent().base_ring().zero()
        for k in s:
            c=cmp(self._d.get(k,zero), other._d.get(k,zero))
            if c!=0:
                return c
        return 0

    def __hash__(self):
        return hash(self._d)

    def _repr_(self):
        return repr(self._d)

    def weight(self, label, e):
        r"""Return the weight of the indexed edge."""
        return self._d.get( (label,e), self.parent().base_ring().zero() )
        
    def weighted_edges(self):
        r"""Return the set of pairs (label,e) representing edges with non-zero weights."""
        return list(self._d.keys())
    
    def edges_with_weights(self):
        r"""
        Returns a list of items of the form ((label,e),w) where (label,e) 
        represents and edge and w represents the non-zero weight assigned."""
        return self._d.items()

class RelativeHomology(Module):
    Element = RelativeHomologyClass
    def __init__(self, surface, base_ring=ZZ):
        self._base_ring=base_ring
        if not isinstance(surface,SimilaritySurface):
            raise ValueError("RelativeHomology only defined for SimilaritySurfaces (and better).")
        self._s=surface
        self._cached_edges=dict()
        Module.__init__(self, base_ring)
        
    def base_ring(self):
        return self._base_ring

    def _element_from_dict(self,d):
        return self.element_class(self, d)

    def _element_constructor_(self, x):
        if instanceof(x, RelativeHomologyClass):
            d=dict()
            for k,v in iteritems(x._d):
                v2=self._base_ring(v)
                if v2!=self._base_ring.zero():
                    d[k]=v2
            return self.element_class(self, d)

    def zero(self):
        return self.element_class(self, dict())

    def __cmp__(self, other):
        if not isinstance(other, RelativeHomology): 
            return cmp(type(other),RelativeHomology)
        c = cmp(self.base_ring(),other.base_ring())
        if c!=0:
            return c
        return cmp(self._s, other._s)

    def edge(self,label,e):
        r"""Return the homology class of the edge with the provided polygon label
        and edge index oriented counter-clockwise around the polygon."""
        try:
            # If already cached, return the cached copy.
            return self._cached_edges[(label,e)]
        except KeyError:
            # not cached!
            num_edges = self._s.polygon(label).num_edges()
            # Check to see if all other edges of the polygon are cached.
            has_all_others = True
            for i in range(1,num_edges):
                e2=(e+i)%num_edges
                if (label,e2) not in self._cached_edges:
                    has_all_others=False
                    break
            if has_all_others:
                # If all the other edges in the polygon are cached then we
                # know this edge's homology class is the negation of their sum.
                e2=(e+1)%num_edges
                total = -self._cached_edges[(label,e2)]
                for i in range(2,num_edges):
                    e2=(e+i)%num_edges
                    total -=  self._cached_edges[(label,e2)]
                # Cache this edge's value and the opposite edge's value.
                self._cached_edges[(label,e)] = total
                label2,e2 = self._s.opposite_edge(label,e)
                self._cached_edges[(label2,e2)] = -total
                return total
            else:
                # At least one other edge is not cached, so we can think of
                # the current edge as a generator.
                d = dict()
                d[(label,e)] = self._base_ring.one()
                v = self._element_from_dict(d)
                # Cache this edge's value and the opposite edge's value.
                self._cached_edges[(label,e)] = v
                label2,e2 = self._s.opposite_edge(label,e)
                self._cached_edges[(label2,e2)] = -v
                return v
