from geometry.polygon import Polygon
from geometry.similarity import SimilarityGroup
from geometry.translation import TranslationGroup
from sage.rings.real_double import RDF
from sage.modules.free_module import VectorSpace

V = VectorSpace(RDF, 2)

class GraphicalPolygon:
    r"""
    Stores data necessary to draw one of the polygons from a surface.
    """
    
    def __init__(self, polygon, transformation=None, outline_color=None, fill_color="#ccc"):
        self._p=polygon
        self.set_transformation(transformation)
        # Store colors
        self.set_outline_color(outline_color)
        self.set_fill_color(fill_color)
        
    def base_ring(self):
        return self._p.base_ring()

    field = base_ring

    def _repr_(self):
        r"""
        String representation.
        """
        return "Graphical Polygon based on "+repr(self._p)

    def base_polygon(self):
        return self._p

    def transformed_vertex(self, e):
        return self._transformation(self._p.vertex(e))

    def set_transformation(self,transformation):
        r"""Set the transformation to be applied to the polygon."""
        if transformation is None:
            self._transformation=TranslationGroup(self._p.base_ring()).one()
        else:
            self._transformation=transformation
        # Cache the location of vertices:
        self._v = [V(self._transformation(v)) for v in self._p.vertices()]

    def set_fill_color(self,fill_color):
        r"""
        Set the fill color.
        """
        self._fill_color=fill_color

    def set_outline_color(self,outline_color):
        r""" 
        Set the outline color.
        """
        self._outline_color=outline_color

    def num_edges(self):
        return self._p.num_edges()

    def vertices(self):
        r"""Return the vertices of the polygon as viewed through the transformation and converted to a
        list of points in VectorSpace(RDF, 2)."""
        return self._v
        
    def plot(self):
        r"""Returns a plot of the GraphicalPolygon.
        
        EXAMPLES::
        
            sage: from geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s=SimilaritySurfaceGenerators.example()
            sage: from graphical.surface import GraphicalSurface
            sage: gs=GraphicalSurface(s)
            sage: gs.graphical_polygon(0).set_fill_color("red")
            sage: show(gs.graphical_polygon(0).plot())
            Launched png viewer for Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.point import point2d
        from sage.plot.polygon import polygon2d
        v = self.vertices()
        from sage.plot.graphics import Graphics
        p = Graphics()
        if not self._fill_color is None:
            if self._outline_color is None:
                p = polygon2d(self.vertices(), color=self._fill_color,axes=False)
            else:
                p = polygon2d(self.vertices(), 
                    color=self._fill_color,edgecolor=self._outline_color,axes=False)
        elif not self._outline_color is None:
            p = polygon2d(self.vertices(), 
                    color=self._outline_color,axes=False,fill=False)
        return p

    def plot_edge(self, e, color=None, dotted=False):
        ne=self.num_edges()
        if color==None:
            color=self._outline_color
        if color==None:
            from sage.plot.graphics import Graphics
            return Graphics()
        from sage.plot.line import line2d
        v=self.vertices()[e]
        w=self.vertices()[(e+1)%ne]
        if dotted:
            return line2d([(v[0],v[1]), (w[0],w[1])],color=color,linestyle=":")
        else:
            return line2d([(v[0],v[1]), (w[0],w[1])],color=color)

