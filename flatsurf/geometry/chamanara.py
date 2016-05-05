r""" 
Construction of Chamanara's surfaces which depend on a parameter alpha less than one.
See the paper "Affine automorphism groups of surfaces of infinite type" in which the surface 
is called $X_\alpha$.

EXAMPLE::
    sage: from geometry.chamanara import GraphicalChamanaraSurface
    sage: s=GraphicalChamanaraSurface(QQ(1)/2,8)
    sage: s.plot()
    Launched png viewer for Graphics object consisting of 65 graphics primitives
"""

def ChamanaraPolygon(alpha):
    from sage.categories.fields import Fields
    field=alpha.parent()
    if not field in Fields:
        ValueError("The value of alpha must lie in a field.")
    if alpha<=0 or alpha>=1:
        ValueError("The value of alpha must be between zero and one.")
    # The value of x is $\sum_{n=0}^\infty \alpha^n$.
    x=1/(1-alpha)
    from geometry.polygon import PolygonCreator
    pc=PolygonCreator(field=field)
    pc.add_vertex((0,0))
    pc.add_vertex((1,0))
    pc.add_vertex((1-x,x))
    pc.add_vertex((1-x,x-1))
    return pc.get_polygon()

from geometry.similarity_surface import SimilaritySurface_generic
from sage.rings.integer_ring import ZZ

class ChamanaraSurface(SimilaritySurface_generic):
    r"""The ChamanaraSurface $X_{\alpha}$."""
    
    def __init__(self, alpha):
        self._p=ChamanaraPolygon(alpha)
        self._field=alpha.parent()
    
    def base_ring(self):
        return self._field
    
    def polygon_labels(self):
        return ZZ
        
    def polygon(self, lab):
        return self._p
    
    def opposite_edge(self, p, e):
        if e==0 or e==2:
            return 1-p,e
        elif e==1:
            if p<0:
                return p+1,3
            elif p>1:
                return p-1,3
            else:
                # p==0 or p==1
                return 1-p,1
        else:
            # e==3
            if p<=0:
                return p-1,1
            else:
                # p>=1
                return p+1,1

    def base_label(self):
        return ZZ(0)

def GraphicalChamanaraSurface(alpha,n):
    r"""Return a standard Graphical version of the ChamanaraSurface $X_{\alpha}$
    with $2n$ polygons shown."""
    s = ChamanaraSurface(alpha)
    from graphical.surface import GraphicalSurface
    gs = GraphicalSurface(s)
    # Make polygon 1 visible
    gs.make_adjacent_and_visible(0,1)
    for i in range(n-1):
        # Make polygon -i-1 visible
        gs.make_adjacent_and_visible(-i,3)
        # Make polygon i+2 visible
        gs.make_adjacent_and_visible(i+1,3)
    return gs

