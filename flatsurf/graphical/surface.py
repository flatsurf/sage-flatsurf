from __future__ import absolute_import

from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.rational_field import QQ

from sage.plot.graphics import Graphics

from flatsurf.geometry.similarity_surface import SimilaritySurface
from .polygon import *
from .edge_gluings import *

class GraphicalSurface:
    r"""
    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.graphical.surface import GraphicalSurface

        sage: s = similarity_surfaces.example()
        sage: gs = GraphicalSurface(s)
        sage: gs.graphical_polygon(0).set_fill_color("red")
        sage: gs.graphical_polygon(0).plot()
        Graphics object consisting of 5 graphics primitives
    """
    def __init__(self, similarity_surface, adjacencies=None, polygon_labels=True, edge_labels=True):
        r"""
        Construct a GraphicalSurface from a similarity surface.

        INPUT:

        - ``similarity_surface`` -- a similarity surface

        - ``polygon_labels`` -- a boolean (default ``True``) whether the label
          of polygons are displayed

        - ``edge_labels`` -- option to control the display of edge labels. It
          can be one of

            - ``False`` or ``None`` for no labels

            - ``'gluings'`` -- to put on each side of each non-adjacent edge, the
              name of the polygon to which it is glued

            - ``'number'`` -- to put on each side of each edge the number of the
              edge

            - ``'gluings and numbers'`` -- full information

        - ``adjacencies`` -- a list of pairs ``(p,e)`` to be used to set
          adjacencies of polygons. 
          
        If adjacencies is not defined and the surface is finite, make_all_visible()
        is called to make all polygons visible.
        """
        self._polygons_labels = polygon_labels

        self._edge_labels = 'gluings' if edge_labels is True else edge_labels

        assert isinstance(similarity_surface, SimilaritySurface)
        self._ss = similarity_surface

        self._polygons = {}
        self._visible = set([self._ss.base_label()])

        if adjacencies is None:
            if self._ss.is_finite():
                self.make_all_visible()
        else:
            for p,e in adjacencies:
                self.make_adjacent_and_visible(p,e)

    def process_options(self, adjacencies=None, polygon_labels=None, edge_labels=None):
        r"""
        Process the options listed as if the graphical_surface was first
        created.
        """
        if not adjacencies is None:
            for p,e in adjacencies:
                self.make_adjacent_and_visible(p,e)
        if not polygon_labels is None:
            self._polygons_labels = polygon_labels
        if not edge_labels is None:
            self._edge_labels = 'gluings' if edge_labels is True else edge_labels

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

    def make_visible(self, label):
        r"""
        Mark the polygon with the given label as visible.
        """
        self._visible.add(label)

    def make_all_visible(self, adjacent=True, limit=None):
        r"""
        Attempt to show all invisible polygons by walking over the surface.

        INPUT:

        - ``adjacent`` -- (default ``True``) whether the newly added polygon are
          set to be adjacent or not
        - ``limit`` -- (default ``None``) maximal number of additional polygons to make visible
        EXAMPLES::

            sage: from flatsurf import *

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface()
            sage: g.make_all_visible()
            sage: g.plot()
            Graphics object consisting of 13 graphics primitives

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(cached=False, adjacencies=[])
            sage: g.make_all_visible(adjacent=False)
            sage: g.plot()
            Graphics object consisting of 16 graphics primitives
        """
        if limit is None:
            assert self._ss.is_finite()
            if adjacent:
                for l,poly in self._ss.label_polygon_iterator():
                    for e in range(poly.num_edges()):
                        l2,e2 = self._ss.opposite_edge(l,e)
                        if not self.is_visible(l2):
                            self.make_adjacent_and_visible(l,e)
            else:
                from flatsurf.geometry.translation import TranslationGroup
                T = TranslationGroup(self._ss.base_ring())
                for l in self._ss.label_iterator():
                    if not self.is_visible(l):
                        poly = self._ss.polygon(l)
                        sxmax = self.xmax()
                        g = self.graphical_polygon(l)
                        pxmin = g.xmin()
                        t = T((self.xmax() - g.xmin() + 1,
                            -(g.ymin()+g.ymax())/2))
                        g.set_transformation(t)
                        self.make_visible(l)
        else:
            assert limit>0
            if adjacent:
                i = 0
                for l,poly in self._ss.label_polygon_iterator():
                    for e in range(poly.num_edges()):
                        l2,e2 = self._ss.opposite_edge(l,e)
                        if not self.is_visible(l2):
                            self.make_adjacent_and_visible(l,e)
                            i=i+1
                            if i>=limit:
                                return
            else:
                from flatsurf.geometry.translation import TranslationGroup
                T = TranslationGroup(self._ss.base_ring())
                i = 0
                for l in self._ss.label_iterator():
                    if not self.is_visible(l):
                        poly = self._ss.polygon(l)
                        sxmax = self.xmax()
                        g = self.graphical_polygon(l)
                        pxmin = g.xmin()
                        t = T((self.xmax() - g.xmin() + 1,
                            -(g.ymin()+g.ymax())/2))
                        g.set_transformation(t)
                        self.make_visible(l)
                        i=i+1
                        if i>=limit:
                            return

    def get_surface(self):
        r"""
        Return the underlying similarity surface.
        """
        return self._ss

    def xmin(self):
        r"""
        Return the minimal x-coordinate of a vertex of a visible graphical polygon.

        .. TODO::

            this should be xmin
        """
        return min([self.graphical_polygon(label).xmin() for label in self.visible()])

    def ymin(self):
        r"""
        Return the minimal y-coordinate of a vertex of a visible graphical polygon.
        """
        return min([self.graphical_polygon(label).ymin() for label in self.visible()])

    def xmax(self):
        r"""
        Return the maximal x-coordinate of a vertex of a visible graphical polygon.
        """
        return max([self.graphical_polygon(label).xmax() for label in self.visible()])

    def ymax(self):
        r"""
        Return the minimal y-coordinate of a vertex of a visible graphical polygon.
        """
        return max([self.graphical_polygon(label).ymax() for label in self.visible()])

    def bounding_box(self):
        r"""
        Return the quadruple (x1,y1,x2,y2) where x1 and y1 are the minimal
        x- and y-coordinates of a visible graphical polygon and x2 and y2 are the
        maximal x-and y- cordinates  of a visible graphical polygon.
        """
        return self.xmin(), self.ymin(), self.xmax(), self.ymax()


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

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface(adjacencies=[])
            sage: print gs.graphical_polygon(0)
            Polygon 0: [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
            sage: print gs.graphical_polygon(1)
            Polygon 1: [(0.0, 0.0), (2.0, 0.0), (1.0, 3.0)]
            sage: print("Polygon 0, edge 0 is opposite "+str(gs.opposite_edge(0,0)))
            Polygon 0, edge 0 is opposite (1, 1)
            sage: gs.make_adjacent(0,0)
            sage: print gs.graphical_polygon(0)
            Polygon 0: [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
            sage: print gs.graphical_polygon(1)
            Polygon 1: [(0.4, -2.8), (2.0, -2.0), (0.0, 0.0)]
        """
        pp,ee = self._ss.opposite_edge(p,e)
        poly = self.graphical_polygon(pp)
        g = self._ss.edge_transformation(pp,ee)
        h = self.graphical_polygon(p).transformation()
        poly.set_transformation(h*g)

    def make_adjacent_and_visible(self, p, e):
        r"""
        Move the polygon across the prescribed edge so that is adjacent,
        and make the moved polygon visible.
        """
        self.make_adjacent(p, e)
        self.make_visible(self._ss.opposite_edge(p,e)[0])

    def is_adjacent(self,p,e):
        r"""
        Returns the truth value of the statement
        'The polygon opposite edge (p,e) is adjacent to that edge.'

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(adjacencies=[])
            sage: g.is_adjacent(0,0)
            False
            sage: g.is_adjacent(0,1)
            False
            sage: g.make_all_visible(adjacent=True)
            sage: g.is_adjacent(0,0)
            True
            sage: g.is_adjacent(0,1)
            False
        """
        pp,ee = self.opposite_edge(p,e)
        if not self.is_visible(pp):
            return False
        g = self.graphical_polygon(p)
        gg = self.graphical_polygon(pp)
        return g.transformed_vertex(e) == gg.transformed_vertex(ee+1) and \
               g.transformed_vertex(e+1) == gg.transformed_vertex(ee)

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        return self._ss.opposite_edge(p,e)

    def edge_labels(self, lab):
        r"""
        Return the edge label to be used for the polygon with label ``lab``

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(cached=False, adjacencies=[])
            sage: g.edge_labels(0)
            ['1', '1', '1']
            sage: g.make_all_visible(adjacent=True)
            sage: g.edge_labels(0)
            [None, '1', '1']

            sage: g.make_adjacent(0,0)
            sage: g.edge_labels(0)
            [None, '1', '1']
            sage: g.edge_labels(1)
            ['0', None, '0']

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(cached=False, adjacencies=[], edge_labels='number')
            sage: g.edge_labels(0)
            ['0', '1', '2']

            sage: g = s.graphical_surface(cached=False, adjacencies=[], edge_labels='gluings and number')
            sage: g.edge_labels(0)
            ['0 -> (1, 1)', '1 -> (1, 2)', '2 -> (1, 0)']
            sage: g.make_all_visible(adjacent=True)
            sage: g.edge_labels(0)
            ['0', '1 -> (1, 2)', '2 -> (1, 0)']
        """
        if not self._edge_labels:
            return None

        s = self._ss
        g = self.graphical_polygon(lab)
        p = g.base_polygon()

        if self._edge_labels == 'gluings':
            ans = []
            for e in range(p.num_edges()):
                if self.is_adjacent(lab, e):
                    ans.append(None)
                else: 
                    llab,ee = s.opposite_edge(lab,e)
                    ans.append(str(llab))
        elif self._edge_labels == 'number':
            ans = map(str, range(p.num_edges()))
        elif self._edge_labels == 'gluings and number':
            ans = []
            for e in range(p.num_edges()):
                if self.is_adjacent(lab, e):
                    ans.append(str(e))
                else:
                    ans.append("{} -> {}".format(e, s.opposite_edge(lab,e)))
        else:
            raise RuntimeError("invalid option")

        return ans


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
            Graphics object consisting of 13 graphics primitives
        """
        p = Graphics()
        for label in self._visible:
            polygon = self.graphical_polygon(label)
            polygon.set_edge_labels(self.edge_labels(label))
            p += polygon.plot()
            for e in range(polygon.base_polygon().num_edges()):
                if not self.is_adjacent(label,e):
                    p += polygon.plot_edge(e,color="blue")
                else:
                    # we want to plot the edges only once!
                    pp,ee = self.opposite_edge(label,e)
                    if label>pp or (label == pp and e > ee):
                        p += polygon.plot_edge(e,color="blue",dotted=True)
        return p


