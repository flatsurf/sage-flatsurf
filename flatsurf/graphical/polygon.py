#*****************************************************************************
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from sage.rings.real_double import RDF
from sage.modules.free_module import VectorSpace
from sage.plot.polygon import polygon2d
from sage.plot.graphics import Graphics
from sage.plot.text import text
from sage.plot.line import line2d
from sage.plot.point import point2d

from flatsurf.geometry.similarity import SimilarityGroup

V = VectorSpace(RDF, 2)

class GraphicalPolygon:
    r"""
    Stores data necessary to draw one of the polygons from a surface.

    Note that this involves converting between geometric coordinates, defined for the SimilaritySurface,
    and graphical coordinates. We do this with a similarity (called transformation below).
    """

    def __init__(self, polygon, transformation=None):
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

        # the following stores _transformation and _v
        self.set_transformation(transformation)

    def copy(self):
        r"""
        Return a copy of this GraphicalPolygon.
        """
        return GraphicalPolygon(self._p, self.transformation())

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface()
            sage: gs.graphical_polygon(0)
            GraphicalPolygon with vertices [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
        """
        return "GraphicalPolygon with vertices {}".format(self._v)

    def base_polygon(self):
        r"""
        Return the polygon of the surface in geometric coordinates.
        """
        return self._p

    def transformed_vertex(self, e):
        r"""
        Return the graphical coordinates of the vertex in double precision.
        """
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

    def transform(self, point, double_precision = True):
        r"""
        Return the transformation of point into graphical coordinates.

        By default returned point is in double precision. This can be changed
        to an exact representation by setting `double_precision` to False.
        """
        if self._transformation is None:
            if double_precision:
                return V(point)
            else:
                return point
        else:
            if double_precision:
                return V(self._transformation(point))
            else:
                return self._transformation(point)

    def transform_back(self, point):
        r"""
        Return the transformation of point from graphical coordinates to the geometric coordinates
        of the underlying SimilaritySurface.
        """
        if self._transformation is None:
            return point
        else:
            return (~self._transformation)(point)

    def contains(self, point):
        r"""
        Return the transformation of point from graphical coordinates to the geometric coordinates
        of the underlying SimilaritySurface.
        """
        return self._p.contains_point(self.transform_back(point))

    def transformation(self):
        r"""
        Return the transformation (similarity) which converts from
        mathematical to graphical coordinates.
        """
        return self._transformation

    def set_transformation(self, transformation = None):
        r"""Set the transformation to be applied to the polygon."""
        if transformation is None:
            self._transformation = SimilarityGroup(self._p.base_ring()).one()
        else:
            self._transformation = transformation
        # recompute the location of vertices:
        self._v = [V(self._transformation(v)) for v in self._p.vertices()]

    def plot_polygon(self, **options):
        r"""
        Returns only the filled polygon.

        Options are processed as in sage.plot.polygon.polygon2d except
        that by default axes=False.
        """
        if not "axes" in options:
            options["axes"] = False
        return polygon2d(self._v, **options)

    def plot_label(self, label, **options):
        r"""
        Write the label of the polygon as text.

        Set ``position`` to a pair (x,y) to determine where
        the label is drawn (in graphical coordinates). If this parameter
        is not provided, the label is positioned in the baricenter
        of the polygon.

        Other options are processed as in sage.plot.text.text.
        """
        if "position" in options:
            return text(str(label), options.pop("position"), **options)
        else:
            return text(str(label), sum(self._v) / len(self._v), **options)

    def plot_edge(self, e, **options):
        r"""
        Plot the edge e, with e a number 0,...,n-1 with n being the number
        of edges of the polygon.

        Options are processed as in sage.plot.line.line2d.
        """
        return line2d([self._v[e], self._v[(e+1)%self.base_polygon().num_edges()]], \
            **options)

    def plot_edge_label(self, i, label, **options):
        r"""
        Write label on the i-th edge.

        A parameter ``t`` in the interval [0,1] can be provided to position the
        label along the edge. A value of t=0 will position it at the starting
        vertex and t=1 will position it at the terminating vertex. Defaults to
        0.3.

        If the parameter ``position`` can take the values "outside", "inside"
        or "edge" to indicate if the label should be drawn outside the polygon,
        inside the polygon or on the edge. Defaults to "inside".

        A ``push_off`` perturbation parameter controls how far off the edge the label is pushed.

        Other options are processed as in sage.plot.text.text.
        """
        e = self._v[(i+1)%self.base_polygon().num_edges()] - self._v[i]

        if "position" in options:
            if options["position"] not in ["inside", "outside", "edge"]:
                raise ValueError("The 'position' parameter must take the value 'inside', 'outside', or 'edge'.")
            pos = options.pop("position")
        else:
            pos = "inside"

        if pos == "outside":
            # position outside polygon.
            if "horizontal_alignment" in options:
                pass
            elif e[1] > 0:
                options["horizontal_alignment"]="left"
            elif e[1] < 0:
                options["horizontal_alignment"]="right"
            else:
                options["horizontal_alignment"]="center"

            if "vertical_alignment" in options:
                pass
            elif e[0] > 0:
                options["vertical_alignment"]="top"
            elif e[0] < 0:
                options["vertical_alignment"]="bottom"
            else:
                options["vertical_alignment"]="center"

        elif pos == "inside":
            # position inside polygon.
            if "horizontal_alignment" in options:
                pass
            elif e[1] < 0:
                options["horizontal_alignment"]="left"
            elif e[1] > 0:
                options["horizontal_alignment"]="right"
            else:
                options["horizontal_alignment"]="center"

            if "vertical_alignment" in options:
                pass
            elif e[0] < 0:
                options["vertical_alignment"]="top"
            elif e[0] > 0:
                options["vertical_alignment"]="bottom"
            else:
                options["vertical_alignment"]="center"

        else:
            # centered on edge.
            if "horizontal_alignment" in options:
                pass
            else:
                options["horizontal_alignment"]="center"
            if "vertical_alignment" in options:
                pass
            else:
                options["vertical_alignment"]="center"

        if "t" in options:
            t = RDF(options.pop("t"))
        else:
            t = 0.3

        if "push_off" in options:
            push_off = RDF(options.pop("push_off"))
        else:
            push_off = 0.03
        if pos == "outside":
            push_off = -push_off
        # Now push_off stores the amount it should be pushed into the polygon

        no = V((-e[1], e[0]))
        return text(label, self._v[i] + t * e + push_off * no, **options)

    def plot_zero_flag(self, **options):
        r"""
        Draw a line segment from the zero vertex toward the baricenter.

        A real parameter ``t`` can be provided. If t=1, then the segment will
        go all the way to the baricenter.  The value of ``t`` is linear in the
        length of the segment. Defaults to t=0.5.

        Other options are processed as in sage.plot.line.line2d.
        """
        if "t" in options:
            t = RDF(options.pop("t"))
        else:
            t = 0.5

        return line2d([self._v[0], self._v[0] + t*(sum(self._v) / len(self._v) - self._v[0])],
            **options)

    def plot_points(self, points, **options):
        r"""
        Plot the points in the given collection of points.

        The options are passed to point2d.

        If no "zorder" option is provided then we set "zorder" to 50.

        By default coordinates are taken in the underlying surface. Call with coordinates="graphical"
        to use graphical coordinates instead.
        """
        if "zorder" not in options:
            options["zorder"]=50
        if "coordinates" not in options:
            points2 = [self.transform(point) for point in points]
        elif options["coordinates"]=="graphical":
            points2=[V(point) for point in points]
            del options["coordinates"]
        else:
            raise ValueError("Invalid value of 'coordinates' option")
        return point2d(points=points2, **options)
   

