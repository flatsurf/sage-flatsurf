from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.real_double import RDF
from sage.rings.rational_field import QQ

from flatsurf.geometry.similarity_surface import SimilaritySurface_generic
from flatsurf.graphical.polygon import *
from flatsurf.graphical.edge_gluings import *


class GraphicalSurface:
    r"""
    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.graphical.surface import GraphicalSurface

        sage: s = similarity_surfaces.example()
        sage: gs = GraphicalSurface(s)
        sage: gs.graphical_polygon(0).set_fill_color("red")
        sage: gs.graphical_polygon(0).plot()
        Graphics object consisting of 2 graphics primitives
    """
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
        label = self._ss.base_label()
        self._visible.add(label)
        self._polygons[label] = GraphicalPolygon(self._ss.polygon(label), label=label)

    def __repr__(self):
        return "Graphical version of Similarity Surface {!r}".format(self._ss)

    def visible(self):
        r"""
        Return the set of visible labels.
        """
        return self._visible

    def is_visible(self,label):
        r"""
        Return whether the polygon with the given label is marked as visible.
        """
        return label in self._visible

    def make_visible(self,label):
        r"""
        Mark the polygon with the given label as visible.
        """
        self._visible.add(label)

    def make_all_visible(self):
        r"""Attempt to show all invisible polygons by walking over the surface."""
        assert self._ss.is_finite()
        for l in self._ss.polygon_labels():
            poly = self._ss.polygon(l)
            for e in range(poly.base_polygon().num_edges()):
                l2,e2 = self._ss.opposite_edge(l,e)
                if not self.is_visible(l2):
                    self.make_adjacent_and_visible(l,e)

    def get_surface(self):
        r"""
        Return the underlying similarity surface.
        """
        return self._ss

    def minx(self):
        r"""
        Return the minimal x-coordinate of a vertex of a visible graphical polygon.

        .. TODO::

            this should be xmin
        """
        return min([self.graphical_polygon(label).minx() for label in self.visible()])

    def miny(self):
        r"""
        Return the minimal y-coordinate of a vertex of a visible graphical polygon.

        .. TODO::

            this should be ymin
        """
        return min([self.graphical_polygon(label).miny() for label in self.visible()])

    def maxx(self):
        r"""
        Return the maximal x-coordinate of a vertex of a visible graphical polygon.

        .. TODO::

            this should be xmax
        """
        return max([self.graphical_polygon(label).maxx() for label in self.visible()])

    def maxy(self):
        r"""
        Return the minimal y-coordinate of a vertex of a visible graphical polygon.

        .. TODO::

            this should be ymax
        """
        return max([self.graphical_polygon(label).maxy() for label in self.visible()])

    def bounding_box(self):
        r"""
        Return the quadruple (x1,y1,x2,y2) where x1 and y1 are the minimal
        x- and y-coordinates of a visible graphical polygon and x2 and y2 are the
        maximal x-and y- cordinates  of a visible graphical polygon.
        """
        return self.minx(), self.miny(), self.maxx(), self.maxy()


    def graphical_polygon(self, label):
        r"""
        Return the graphical_polygon with the given label.
        """
        if self._polygons.has_key(label):
            return self._polygons[label]
        else:
            p = GraphicalPolygon(self._ss.polygon(label), label=label)
            self._polygons[label] = p
            return p

    def make_adjacent(self, p, e):
        r"""
        Move the polygon across the prescribed edge so that is adjacent.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.graphical.surface import GraphicalSurface
            sage: gs = GraphicalSurface(s)
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
        pp,ee = self._ss.opposite_edge(p,e)
        poly = self.graphical_polygon(pp)
        g = self._ss.edge_transformation(pp,ee)
        h = self.graphical_polygon(p).transformation()
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
        pp,ee = self.opposite_edge(p,e)
        return self.graphical_polygon(p).transformed_vertex(e)==self.graphical_polygon(pp).transformed_vertex(ee+1) and \
            self.graphical_polygon(p).transformed_vertex(e+1)==self.graphical_polygon(pp).transformed_vertex(ee)

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        return self._ss.opposite_edge(p,e)

    def plot(self):
        r"""
        Returns a plot of the GraphicalSurface

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.graphical.surface import GraphicalSurface
            sage: gs = GraphicalSurface(s)
            sage: gs.make_visible(1)
            sage: gs.plot()
            Graphics object consisting of 9 graphics primitives
        """
        i = iter(self._visible)
        label = i.next()
        polygon = self.graphical_polygon(label)
        p = polygon.plot()
        for e in range(polygon.base_polygon().num_edges()):
            if not self.is_adjacent(label,e):
                p += polygon.plot_edge(e,color="blue")
            else:
                pp,ee = self.opposite_edge(label,e)
                if label>pp or (label == pp and e > ee):
                    p += polygon.plot_edge(e,color="blue",dotted=True)
        for label in i:
            polygon = self.graphical_polygon(label)
            p += polygon.plot()
            for e in range(polygon.base_polygon().num_edges()):
                if not self.is_adjacent(label,e):
                   p += polygon.plot_edge(e,color="blue")
                else:
                    pp,ee = self.opposite_edge(label,e)
                    if label>pp or (label == pp and e > ee):
                        p += polygon.plot_edge(e,color="blue",dotted=True)
        return p


