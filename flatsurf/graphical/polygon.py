# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2016-2019 Vincent Delecroix
#                     2016-2018 W. Patrick Hooper
#                     2022-2023 Julian Rüth
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with sage-flatsurf. If not, see <http://www.gnu.org/licenses/>.
# ****************************************************************************
from sage.rings.real_double import RDF
from sage.modules.free_module import VectorSpace
from sage.plot.polygon import polygon2d
from sage.plot.text import text
from sage.plot.line import line2d
from sage.plot.point import point2d

from sage.misc.cachefunc import cached_method

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
        """
        self._polygon = polygon
        self._transformation = (
            transformation or SimilarityGroup(self._polygon.base_ring()).one()
        )

    def copy(self):
        r"""
        Return a copy of this GraphicalPolygon.
        """
        import warnings

        warnings.warn(
            "copy() has been deprecated as a method of GraphicalPolygon and will be removed in a future version of sage-flatsurf; create a copy manually instead"
        )

        return GraphicalPolygon(self._polygon, self.transformation())

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface()
            sage: gs.graphical_polygon(0)
            GraphicalPolygon(vertices=[(0, 0), (2, -2), (2, 0)])

        """
        return "Graphical" + repr(self.polygon())

    @cached_method
    def polygon(self):
        return self._transformation * self._polygon

    def base_polygon(self):
        r"""
        Return the polygon of the surface in geometric coordinates.
        """
        return self._polygon

    def transformed_vertex(self, v):
        r"""
        Return the graphical coordinates of the vertex in double precision.
        """
        return self._transformation(self._polygon.vertex(v))

    def transformed_side(self, e):
        return self._transformation * self._polygon.side(e)

    def xmin(self):
        r"""
        Return the minimal x-coordinate of a vertex.
        """
        return min(v[0] for v in self.polygon().vertices())

    def ymin(self):
        r"""
        Return the minimal y-coordinate of a vertex.
        """
        return min(v[1] for v in self.polygon().vertices())

    def xmax(self):
        r"""
        Return the maximal x-coordinate of a vertex.
        """
        return max(v[0] for v in self.polygon().vertices())

    def ymax(self):
        r"""
        Return the minimal y-coordinate of a vertex
        """
        return max(v[1] for v in self.polygon().vertices())

    def bounding_box(self):
        r"""
        Return the quadruple (x1,y1,x2,y2) where x1 and y1 are the minimal
        x and y coordinates and x2 and y2 are the maximal x and y coordinates.
        """
        return self.xmin(), self.ymin(), self.xmax(), self.ymax()

    def transform(self, point, double_precision=True):
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
        return self._polygon.contains_point(self.transform_back(point))

    def transformation(self):
        r"""
        Return the transformation (similarity) which converts from
        mathematical to graphical coordinates.
        """
        return self._transformation

    def set_transformation(self, transformation=None):
        r"""Set the transformation to be applied to the polygon."""
        import warnings

        warnings.warn(
            "set_transformation() has been deprecated and will be removed in a future version of sage-flatsurf; set the transformation when creating the graphical polygon instead"
        )

        if transformation is None:
            transformation = SimilarityGroup(self._polygon.base_ring()).one()

        self._transformation = transformation

    def plot_polygon(self, **options):
        r"""
        Returns only the filled polygon.

        Options are processed as in sage.plot.polygon.polygon2d except
        that by default axes=False.
        """
        if "axes" not in options:
            options["axes"] = False

        return self.polygon().plot(
            polygon_options=options, edge_options=None, vertex_options=None
        )

    def plot_label(self, label, **options):
        r"""
        Write the label of the polygon as text.

        Set ``position`` to a pair (x,y) to determine where
        the label is drawn (in graphical coordinates). If this parameter
        is not provided, the label is positioned in the baricenter
        of the polygon.

        Other options are processed as in sage.plot.text.text.
        """
        from sage.all import Graphics

        g = Graphics()

        def position(xlim, ylim):
            from flatsurf.graphical.hyperbolic import CartesianPathPlot

            path = CartesianPathPlot(self.polygon()._plot_commands())._create_path(
                xlim, ylim, fill=True
            )

            from matplotlib.transforms import Bbox

            path = path.clip_to_bbox(Bbox(((xlim[0], ylim[0]), (xlim[1], ylim[1]))))
            if path.codes is None:
                return None

            from sage.all import vector

            vertices = []
            for j, (vertex, code) in enumerate(zip(path.vertices, path.codes)):
                if code == 2:
                    if j > 0 and path.codes[j - 1] == 1:
                        vertices.append(path.vertices[j - 1])
                        if len(vertices) >= 2 and vector(vertices[-1]) == vector(
                            vertices[-2]
                        ):
                            vertices.pop()
                    vertices.append(vertex)
                    if len(vertices) >= 2 and vector(vertices[-1]) == vector(
                        vertices[-2]
                    ):
                        vertices.pop()

            vertices = list(map(vector, vertices))
            if vertices[0] == vertices[-1]:
                vertices.pop()

            if len(vertices) <= 1:
                return None

            from flatsurf import Polygon, EuclideanPlane
            from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
            from sage.all import RR

            # TODO: The check here is pretty evil. There is one vertex too many in closed polygons.
            centroid = (
                EuclideanPlane(RR, geometry=EuclideanEpsilonGeometry(RR, 1e-3))
                .polygon(vertices=vertices, check=False)
                .centroid()
            )
            return (float(centroid[0]), float(centroid[1]))

        position = options.pop("position", position)

        if "horizontal_alignment" not in options:
            options["horizontal_alignment"] = "center"
        if "vertical_alignment" not in options:
            options["vertical_alignment"] = "center"

        options = text(label, (0, 0), **options)[0].options()

        from flatsurf.graphical.hyperbolic import DynamicLabel

        g.add_primitive(DynamicLabel(str(label), position, options))
        return g

    def plot_edge(self, e, **options):
        r"""
        Plot the edge e, with e a number 0,...,n-1 with n being the number
        of edges of the polygon.

        Options are processed as in sage.plot.line.line2d.
        """
        return self.polygon().side(e).plot(**options)

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
        side = self.polygon().side(i)
        direction = side.direction()

        from sage.all import sgn

        direction_sgn = (sgn(direction[0]), sgn(direction[1]))

        if "position" in options:
            if options["position"] not in ["inside", "outside", "edge"]:
                raise ValueError(
                    "The 'position' parameter must take the value 'inside', 'outside', or 'edge'."
                )
            pos = options.pop("position")
        else:
            pos = "inside"

        if pos == "outside":
            # position outside polygon.
            if "horizontal_alignment" in options:
                pass
            elif direction[1] > 0:
                options["horizontal_alignment"] = "left"
            elif direction[1] < 0:
                options["horizontal_alignment"] = "right"
            else:
                options["horizontal_alignment"] = "center"

            if "vertical_alignment" in options:
                pass
            elif direction[0] > 0:
                options["vertical_alignment"] = "top"
            elif direction[0] < 0:
                options["vertical_alignment"] = "bottom"
            else:
                options["vertical_alignment"] = "center"

        elif pos == "inside":
            # position inside polygon.

            if "horizontal_alignment" not in options:
                options["horizontal_alignment"] = {
                    (-1, -1): "left",
                    (-1, 0): "right",
                    (-1, 1): "right",
                    (0, -1): "left",
                    (0, 1): "right",
                    (1, -1): "left",
                    (1, 0): "left",
                    (1, 1): "right",
                }[direction_sgn]

            if "vertical_alignment" not in options:
                options["vertical_alignment"] = {
                    (-1, -1): "top",
                    (-1, 0): "top",
                    (-1, 1): "top",
                    (0, -1): "top",
                    (0, 1): "bottom",
                    (1, -1): "bottom",
                    (1, 0): "bottom",
                    (1, 1): "bottom",
                }[direction_sgn]

        else:
            # centered on edge.
            if "horizontal_alignment" in options:
                pass
            else:
                options["horizontal_alignment"] = "center"
            if "vertical_alignment" in options:
                pass
            else:
                options["vertical_alignment"] = "center"

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

        normal = V((-direction[1], direction[0]))
        normal /= normal.norm()

        from sage.all import Graphics

        g = Graphics()

        def position(xlim, ylim):
            from flatsurf.graphical.hyperbolic import CartesianPathPlot

            path = CartesianPathPlot(
                self.polygon().side(i)._plot_commands()
            )._create_path(xlim, ylim, fill=False)

            from sage.all import vector

            vertices = []
            for j, (vertex, code) in enumerate(zip(path.vertices, path.codes)):
                if code == 2:
                    if j > 0 and path.codes[j - 1] == 1:
                        vertices.append(path.vertices[j - 1])
                        if len(vertices) >= 2 and vector(vertices[-1]) == vector(
                            vertices[-2]
                        ):
                            vertices.pop()
                    vertices.append(vertex)
                    if len(vertices) >= 2 and vector(vertices[-1]) == vector(
                        vertices[-2]
                    ):
                        vertices.pop()

            vertices = list(map(vector, vertices))
            if vertices[0] == vertices[-1]:
                vertices.pop()

            if len(vertices) < 2:
                raise NotImplementedError("edge is not visible")

            assert (
                len(vertices) == 2
            ), "path should render as exactly two points if it is visible"

            vertex = vertices[0]
            vector = vertices[1] - vertices[0]
            return vertex + t * vector + push_off * normal

        position = options.pop("position", position)

        options = text(label, (0, 0), **options)[0].options()

        from flatsurf.graphical.hyperbolic import DynamicLabel

        g.add_primitive(DynamicLabel(label, position, options))
        return g

    def plot_zero_flag(self, **options):
        r"""
        Draw a line segment from the zero vertex toward the baricenter.

        A real parameter ``t`` can be provided. If t=1, then the segment will
        go all the way to the baricenter.  The value of ``t`` is linear in the
        length of the segment. Defaults to t=0.5.

        Other options are processed as in sage.plot.line.line2d.
        """
        raise NotImplementedError
        # if "t" in options:
        #     t = RDF(options.pop("t"))
        # else:
        #     t = 0.5

        # return line2d(
        #     [self._v[0], self._v[0] + t * (sum(self._v) / len(self._v) - self._v[0])],
        #     **options,
        # )

    def plot_points(self, points, **options):
        r"""
        Plot the points in the given collection of points.

        The options are passed to point2d.

        If no "zorder" option is provided then we set "zorder" to 50.

        By default coordinates are taken in the underlying surface. Call with coordinates="graphical"
        to use graphical coordinates instead.
        """
        raise NotImplementedError
        if "zorder" not in options:
            options["zorder"] = 50
        if "coordinates" not in options:
            points2 = [self.transform(point) for point in points]
        elif options["coordinates"] == "graphical":
            points2 = [V(point) for point in points]
            del options["coordinates"]
        else:
            raise ValueError("Invalid value of 'coordinates' option")
        return point2d(points=points2, **options)
