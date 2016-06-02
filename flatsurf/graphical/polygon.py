from __future__ import absolute_import

from sage.rings.real_double import RDF
from sage.modules.free_module import VectorSpace

from flatsurf.geometry.similarity import SimilarityGroup
from flatsurf.geometry.translation import TranslationGroup

V = VectorSpace(RDF, 2)

class GraphicalPolygon:
    r"""
    Stores data necessary to draw one of the polygons from a surface.
    """

    def __init__(self, polygon, transformation=None, outline_color=None,
            fill_color="#ccc", label=None, edge_labels=True):
        r"""
        INPUT:

        - ``polygon`` -- the actual polygon

        - ``transformation`` -- a transformation to be applied to the polygon

        - ``outline_color`` -- a color

        - ``fill_color`` -- another color

        - ``label`` -- an optional label for the polygon

        - ``edge_labels`` -- one of ``False``, ``True`` or a list of labels
        """
        self._p = polygon
        self._edge_labels = edge_labels

        # the following stores _transformation and _v
        self.set_transformation(transformation)

        # store colors
        self.set_outline_color(outline_color)
        self.set_fill_color(fill_color)
        self.set_label(label)

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface()
            sage: gs.graphical_polygon(0)
            Polygon 0: [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
        """
        return "Polygon {}: {}".format(self._label, self._v)

    def base_polygon(self):
        return self._p

    def transformed_vertex(self, e):
        return self._transformation(self._p.vertex(e))

    def xmin(self):
        r"""
        Return the minimal x-coordinate of a vertex.

        .. TODO::

            to fit with Sage conventions this should be xmin
        """
        return min([v[0] for v in self._v])

    def ymin(self):
        r"""
        Return the minimal y-coordinate of a vertex.

        .. TODO::

            to fit with Sage conventions this should be ymin
        """
        return min([v[1] for v in self._v])

    def xmax(self):
        r"""
        Return the maximal x-coordinate of a vertex.

        .. TODO::

            to fit with Sage conventions this should be xmax
        """
        return max([v[0] for v in self._v])

    def ymax(self):
        r"""
        Return the minimal y-coordinate of a vertex

        .. TODO::

            To fit with Sage conventions this should be ymax
        """
        return max([v[1] for v in self._v])

    def bounding_box(self):
        r"""
        Return the quadruple (x1,y1,x2,y2) where x1 and y1 are the minimal
        x- and y-coordinates and x2 and y2 are the maximal x-and y- cordinates.
        """
        return self.xmin(), self.ymin(), self.xmax(), self.ymax()

    def transform(self, point, field=None):
        r"""
        Return the transformation of point into graphical coordinates.
        """
        if self._transformation is None:
            return V(point)
        else:
            return V(self._transformation(point))

    def transformation(self):
        r"""
        Return the transformation (similarity) which converts from
        mathematical to graphical coordinates.
        """
        return self._transformation

    def set_transformation(self, transformation):
        r"""Set the transformation to be applied to the polygon."""
        if transformation is None:
            self._transformation = TranslationGroup(self._p.base_ring()).one()
        else:
            self._transformation = transformation
        # recompute the location of vertices:
        self._v = [V(self._transformation(v)) for v in self._p.vertices()]

    def set_fill_color(self, fill_color):
        r"""
        Set the fill color.
        """
        self._fill_color=fill_color

    def set_outline_color(self, outline_color):
        r"""
        Set the outline color.
        """
        self._outline_color=outline_color

    def set_label(self, label):
        self._label = label

    def set_edge_labels(self, edge_labels):
        self._edge_labels = edge_labels

    def polygon_options(self):
        d = {'axes': False}
        if self._fill_color is not None:
            d['color'] = self._fill_color
            if self._outline_color is not None:
                d['edgecolor'] = self._outline_color
        elif self._outline_color is not None:
            d['color'] = self._outline_color,
            d['fill'] = False

        return d

    def polygon_label_options(self):
        return {'color': 'black'}

    def edge_label_options(self):
        if self._outline_color is not None:
            return {'color': self._outline_color}
        return {}

    def plot(self, draw_polygon_label=True):
        r"""
        Returns a plot of the GraphicalPolygon.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: from flatsurf.graphical.surface import GraphicalSurface
            sage: gs = GraphicalSurface(s)
            sage: gs.graphical_polygon(0).set_fill_color("red")
            sage: gs.graphical_polygon(0).plot()
            Graphics object consisting of 5 graphics primitives
        """
        from sage.plot.point import point2d
        from sage.plot.polygon import polygon2d
        from sage.plot.graphics import Graphics
        from sage.plot.text import text

        p = polygon2d(self._v, **self.polygon_options())

        if self._label is not None and draw_polygon_label:
            p += text(str(self._label), sum(self._v) / len(self._v),
                    **self.polygon_label_options())

        if self._edge_labels:
            opt = self.edge_label_options()
            n = self.base_polygon().num_edges()
            for i in range(n):
                lab = str(i) if self._edge_labels is True else self._edge_labels[i]
                if lab:
                    e = self._v[(i+1)%n] - self._v[i]
                    no = V((-e[1], e[0]))
                    p += text(lab, self._v[i] + 0.3 * e + 0.05  * no, **self.edge_label_options())

        return p

    def plot_edge(self, e, color=None, dotted=False):
        ne = self.base_polygon().num_edges()
        if color is None:
            color = self._outline_color
        if color is None:
            from sage.plot.graphics import Graphics
            return Graphics()

        from sage.plot.line import line2d
        v = self._v[e]
        w = self._v[(e+1)%ne]
        if dotted:
            return line2d([v, w], color=color, linestyle=":")
        else:
            return line2d([v, w], color=color)

    # DEPRECATED METHODS THAT WILL BE REMOVED

    def vertices(self):
        r"""
        Return the vertices of the polygon as a list of floating point vectors.
        """
        from sage.misc.superseded import deprecation
        deprecation(1, "do not use vertices")
        return self._v

    def num_edges(self):
        from sage.misc.superseded import deprecation
        deprecation(1,"do not use num_edges but .base_polygon().num_edges()")
        return self._p.num_edges()

    def base_ring(self):
        from sage.misc.superseded import deprecation
        deprecation(1,"do not use .base_ring() but .base_polygon().base_ring()")
        return self._p.base_ring()

    field = base_ring


