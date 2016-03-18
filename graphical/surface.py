from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.real_double import RDF
from sage.rings.rational_field import QQ

from geometry.similarity_surface import SimilaritySurface_generic

from graphical.polygon import *


from graphical.edge_gluings import *


class GraphicalSurface:

    def __init__(self, similarity_surface, name=None):
        r"""
        Construct a GraphicalSurface from a similarity surface.
        """
        if name is None:
            name="Graphical Surface based on "+repr(similarity_surface)
        assert isinstance(similarity_surface, SimilaritySurface_generic)
        self._ss = similarity_surface
        self._field = self._ss.base_ring()
        # Polygons
        self._polygons = {}
        # Set of visible polygons
        self._visible = set()
        label=self._ss.base_label()
        self._visible.add(label)
        self._polygons[label]=GraphicalPolygon(self._ss.polygon(label))
        
    def __repr__(self):
        s = "Graphical version of Similarity Surface "+repr(self._ss)
        return s
        
    def is_visible(self,label):
        r"""
        Return true if the polygon with the given label is marked as visible.
        """
        return label in self._visible

    def make_visible(self,label):
        r"""
        Mark the polygon with the given label as visible.
        """
        self._visible.add(label)

    def get_surface(self):
        r"""
        Return the underlying similarity surface.
        """
        return self._ss
        
    def graphical_polygon(self, label):
        r"""
        Return the graphical_polygon with the given label.
        """
        if self._polygons.has_key(label):
            return self._polygons[label]
        else:
            p=GraphicalPolygon(self._ss.polygon(label))
            self._polygons[label]=p
            return p
    
    def make_adjacent(self, p, e):
        r"""Move the polygon across the prescribed edge so that is adjacent.
        
        EXAMPLES::
        
            sage: from geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s=SimilaritySurfaceGenerators.example()
            sage: from graphical.surface import GraphicalSurface
            sage: gs=GraphicalSurface(s)
            sage: print("Polygon 0: "+str(gs.graphical_polygon(0).vertices()))
            Polygon 0: [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
            sage: print("Polygon 1: "+str(gs.graphical_polygon(1).vertices()))
            Polygon 1: [(0.0, 0.0), (2.0, 0.0), (1.0, 3.0)]
            sage: print("Polygon 0, edge 0 is opposite "+str(gs.opposite_edge(0,0)))
            Polygon 0, edge 0 is opposite (1, 1)
            sage: gs.make_adjacent(0,0)
            sage: print("Polygon 0: "+str(gs.graphical_polygon(0).vertices()))
            Polygon 0: [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
            sage: print("Polygon 1: "+str(gs.graphical_polygon(1).vertices()))
            Polygon 1: [(0.4, -2.8), (2.0, -2.0), (0.0, 0.0)]
        """
        pp,ee=self._ss.opposite_edge(p,e)
        poly=self.graphical_polygon(pp)
        g=self._ss.edge_transformation(pp,ee)
        h=self.graphical_polygon(p).transformation()
        poly.set_transformation(h*g)

    def make_adjacent_and_visible(self, p, e):
        r"""Move the polygon across the prescribed edge so that is adjacent,
        and make the moved polygon visible."""
        self.make_adjacent(p, e)
        self.make_visible(self._ss.opposite_edge(p,e)[0])
        
    def is_adjacent(self,p,e):
        r"""
        Returns the truth value of the statement 
        'The polygon opposite edge (p,e) is adjacent to that edge.'
        """
        pp,ee=self.opposite_edge(p,e)
        return self.graphical_polygon(p).transformed_vertex(e)==self.graphical_polygon(pp).transformed_vertex(ee+1) and \
            self.graphical_polygon(p).transformed_vertex(e+1)==self.graphical_polygon(pp).transformed_vertex(ee)

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        return self._ss.opposite_edge(p,e)
    
    def plot(self):
        r"""Returns a plot of the GraphicalSurface
        
        EXAMPLES::
        
            sage: from geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s=SimilaritySurfaceGenerators.example()
            sage: from graphical.surface import GraphicalSurface
            sage: gs=GraphicalSurface(s)
            sage: gs.make_visible(1)
            sage: show(gs.plot())
            Launched png viewer for Graphics object consisting of 7 graphics primitives
        """
        i=iter(self._visible)
        label=i.next()
        polygon=self.graphical_polygon(label)
        p=polygon.plot()
        for e in range(polygon.num_edges()):
            if not self.is_adjacent(label,e):
                p += polygon.plot_edge(e,color="blue")
            else:
                pp,ee=self.opposite_edge(label,e)
                if label>pp or (label==pp and e>ee):
                    p += polygon.plot_edge(e,color="blue",dotted=True)
        while True:
            try:
                label=i.next()
                polygon=self.graphical_polygon(label)
                p+=polygon.plot()
                for e in range(polygon.num_edges()):
                    if not self.is_adjacent(label,e):
                       p += polygon.plot_edge(e,color="blue")
                    else:
                        pp,ee=self.opposite_edge(label,e)
                        if label>pp or (label==pp and e>ee):
                            p += polygon.plot_edge(e,color="blue",dotted=True)
            except StopIteration as e:
                break
        return p


