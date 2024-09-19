r"""
The category of Euclidean polygons defined in the real plane.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category of a polygon is automatically determined.

EXAMPLES::

    sage: from flatsurf.geometry.categories import EuclideanPolygons
    sage: C = EuclideanPolygons(QQ)

    sage: from flatsurf import polygons
    sage: polygons.square() in C
    True

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.cachefunc import cached_method
from sage.all import FreeModule
from sage.misc.abstract_method import abstract_method
from sage.structure.element import get_coercion_model

from flatsurf.geometry.categories.polygons import Polygons
from flatsurf.geometry.euclidean import ccw

cm = get_coercion_model()


class EuclideanPolygons(Category_over_base_ring):
    r"""
    The category of Euclidean polygons defined in the real plane
    over a fixed base ring.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import EuclideanPolygons
        sage: EuclideanPolygons(QQ)
        Category of euclidean polygons over Rational Field

    """

    def super_categories(self):
        r"""
        Return the categories Euclidean polygons are also contained in, namely
        the polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: EuclideanPolygons(QQ).super_categories()
            [Category of polygons over Rational Field]

        """
        return [Polygons(self.base_ring())]

    class ParentMethods:
        r"""
        Provides methods available to all Euclidean polygons in the real plane.

        If you want to add functionality to such polygons, you probably want to
        put it here.
        """

        def vector_space(self):
            r"""
            Return the vector space of dimension 2 in which this polygon embeds.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.vector_space()
                doctest:warning
                ...
                UserWarning: vector_space() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring().fraction_field()**2 instead
                Vector space of dimension 2 over Rational Field

            """
            import warnings

            warnings.warn(
                "vector_space() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring().fraction_field()**2 instead"
            )

            return self.base_ring().fraction_field() ** 2

        def module(self):
            r"""
            Return the free module of rank 2 in which this polygon embeds.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.module()
                doctest:warning
                ...
                UserWarning: module() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring()**2 instead
                Vector space of dimension 2 over Rational Field

            """
            import warnings

            warnings.warn(
                "module() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring()**2 instead"
            )

            return self.base_ring() ** 2

        def field(self):
            r"""
            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.field()
                doctest:warning
                ...
                UserWarning: field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() or base_ring().fraction_field() instead
                Rational Field

            """
            import warnings

            warnings.warn(
                "field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() or base_ring().fraction_field() instead"
            )

            return self.base_ring().fraction_field()

        def _mul_(self, g, switch_sides=None):
            r"""
            Apply the 2x2 matrix `g` to this polygon.

            The matrix must have non-zero determinant. If the determinant is
            negative, then the vertices and edges are relabeled according to the
            involutions `v \mapsto (n-v)%n` and  `e \mapsto n-1-e` respectively.

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices = [(1,0),(0,1),(-1,-1)])
                sage: p
                Polygon(vertices=[(1, 0), (0, 1), (-1, -1)])

                sage: matrix(ZZ,[[0, 1], [1, 0]]) * p
                Polygon(vertices=[(0, 1), (-1, -1), (1, 0)])

                sage: matrix(ZZ,[[2, 0], [0, 1]]) * p
                Polygon(vertices=[(2, 0), (0, 1), (-2, -1)])

            """
            from flatsurf import Polygon

            if g in self.base_ring():
                from sage.all import MatrixSpace

                g = MatrixSpace(self.base_ring(), 2)(g)

            det = g.det()
            if det == 0:
                raise ValueError(
                    "Can not act on a polygon with matrix with zero determinant"
                )

            if det < 0:
                # Note that in this case we reverse the order
                vertices = [g * self.vertex(0)]
                for i in range(len(self.vertices()) - 1, 0, -1):
                    vertices.append(g * self.vertex(i))

                return Polygon(vertices=vertices, check=False)

            return Polygon(
                vertices=[g * v for v in self.vertices()],
                check=False,
                category=self.category(),
            )

        @cached_method
        def is_rational(self):
            r"""
            Return whether this is a rational polygon, i.e., all its
            :meth:`angles` are rational multiples of π.

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices = [(1, 0), (0, 1), (-1, -1)])
                sage: p.is_rational()
                False

            Note that determining rationality is somewhat costly. Once
            established, this refines the category of the triangle::

                sage: p = Polygon(vertices = [(0, 0), (1, 0), (0, 1)])
                sage: p.category()
                Category of convex simple euclidean polygons over Rational Field
                sage: p.is_rational()
                True
                sage: p.category()
                Category of rational convex simple euclidean polygons over Rational Field

            """
            for e in range(len(self.vertices())):
                u = self.edge(e)
                v = -self.edge((e - 1) % len(self.vertices()))

                cos = u.dot_product(v)
                sin = u[0] * v[1] - u[1] * v[0]

                from flatsurf.geometry.euclidean import is_cosine_sine_of_rational

                if not is_cosine_sine_of_rational(cos, sin, scaled=True):
                    return False

            self._refine_category_(self.category().Rational())

            return True

        def is_simple(self):
            r"""
            Return whether this is a simple polygon, i.e., without
            self-intersection.

            EXAMPLES:

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.is_simple()
                True

            """
            n = len(self.vertices())
            for i in range(n):
                ei = (self.vertex(i), self.vertex(i + 1))
                for j in range(i + 2, n + 1):
                    if (i - j) % n in [-1, 0, 1]:
                        continue

                    ej = (self.vertex(j), self.vertex(j + 1))

                    from flatsurf.geometry.euclidean import is_segment_intersecting

                    if is_segment_intersecting(ei, ej):
                        return False

            return True

        @abstract_method
        def vertices(self, marked_vertices=True):
            r"""
            Return the vertices of this polygon in counterclockwise order as
            vectors in the real plane.

            INPUT:

            - ``marked_vertices`` -- a boolean (default: ``True``); whether to
              include marked vertices that are not actually corners of the
              polygon.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.vertices()
                ((0, 0), (1, 0), (1, 1), (0, 1))

            """

        def vertex(self, i):
            r"""
            Return coordinates for the ``i``-th vertex of this polygon.

            EXAMPLES:

            The ``i`` wraps around if it is negative or exceeds the number of
            vertices in this polygon::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.vertex(-1)
                (0, 1)
                sage: s.vertex(0)
                (0, 0)
                sage: s.vertex(1)
                (1, 0)
                sage: s.vertex(2)
                (1, 1)
                sage: s.vertex(3)
                (0, 1)
                sage: s.vertex(4)
                (0, 0)

            """
            vertices = self.vertices()
            return vertices[i % len(vertices)]

        def edges(self):
            r"""
            Return the edges of this polygon as vectors in the plane going from
            one vertex to the next one.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.edges()
                [(1, 0), (0, 1), (-1, 0), (0, -1)]

            """
            return [self.edge(i) for i in range(len(self.vertices()))]

        def edge(self, i):
            r"""
            Return the vector going from vertex ``i`` to the following vertex
            in counter-clockwise order.

            EXAMPLES:

            Note that this wraps around if ``i`` is negative or exceeds the
            number of vertices::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.edge(-1)
                (0, -1)
                sage: s.edge(0)
                (1, 0)
                sage: s.edge(1)
                (0, 1)
                sage: s.edge(2)
                (-1, 0)
                sage: s.edge(3)
                (0, -1)
                sage: s.edge(4)
                (1, 0)

            """
            return self.vertex(i + 1) - self.vertex(i)

        def is_convex(self, strict=False):
            r"""
            Return whether this is a convex polygon.

            INPUT:

            - ``strict`` -- whether to check for strict convexity, i.e., a
              polygon with a π angle is not considered convex.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.is_convex()
                True
                sage: S.is_convex(strict=True)
                True

            """
            from flatsurf.geometry.euclidean import ccw

            for i in range(len(self.vertices())):
                consecutive_ccw = ccw(self.edge(i), self.edge(i + 1))
                if strict:
                    if consecutive_ccw <= 0:
                        return False
                else:
                    if consecutive_ccw < 0:
                        return False

            return True

        def _test_marked_vertices(self, **options):
            r"""
            Verify that :meth:`vertices` and :meth:`is_convex` are compatible.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S._test_marked_vertices()

            """
            tester = self._tester(**options)

            if self.is_convex():
                tester.assertEqual(
                    self.is_convex(strict=True),
                    self.vertices() == self.vertices(marked_vertices=False),
                )

        def is_degenerate(self):
            r"""
            Return whether this polygon is considered degenerate.

            This implements
            :meth:`flatsurf.geometry.categories.polygons.Polygons.ParentMethods.is_degenerate`.

            EXAMPLES:

            Polygons with zero area are considered degenerate::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (1, 0)], check=False)
                sage: p.is_degenerate()
                True

            Polygons with marked vertices are considered degenerate::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (4, 0), (2, 2)])
                sage: p.is_degenerate()
                True

            """
            if self.area() == 0:
                return True

            if self.vertices() != self.vertices(marked_vertices=False):
                return True

            return False

        def slopes(self, relative=False):
            r"""
            Return the slopes of this polygon as vectors in the plane.

            INPUT:

            - ``relative`` -- a boolean (default: ``False``); whether to return the
              slopes not as absolute vectors parallel to the edges but relative to
              the previous edge, i.e., after turning the previous edge to be
              parallel to the x axis.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: s = polygons.square()
                sage: s.slopes()
                [(1, 0), (0, 1), (-1, 0), (0, -1)]

                sage: s.slopes(relative=True)
                [(0, 1), (0, 1), (0, 1), (0, 1)]

            A polygon with a marked point::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (4, 0), (2, 2)])
                sage: p.slopes()
                [(2, 0), (2, 0), (-2, 2), (-2, -2)]
                sage: p.slopes(relative=True)
                [(-4, 4), (4, 0), (-4, 4), (0, 8)]

            """
            if not relative:
                return self.edges()

            edges = [
                (self.edge((e - 1) % len(self.vertices())), self.edge(e))
                for e in range(len(self.vertices()))
            ]

            cos = [u.dot_product(v) for (u, v) in edges]
            sin = [u[0] * v[1] - u[1] * v[0] for (u, v) in edges]

            from sage.all import vector

            return [vector((c, s)) for (c, s) in zip(cos, sin)]

        def erase_marked_vertices(self):
            r"""
            Return a copy of this polygon without marked vertices.

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (4, 0), (2, 2)])
                sage: p.erase_marked_vertices()
                Polygon(vertices=[(0, 0), (4, 0), (2, 2)])

            """
            from flatsurf import Polygon

            return Polygon(vertices=self.vertices(marked_vertices=False))

        def is_equilateral(self):
            r"""
            Return whether all sides of this polygon have the same length.

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)])
                sage: p.is_equilateral()
                True

            """
            return len({edge[0] ** 2 + edge[1] ** 2 for edge in self.edges()}) == 1

        def is_equiangular(self):
            r"""
            Return whether all sides of this polygon meet at the same angle.

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: p = Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)])
                sage: p.is_equiangular()
                True

            """
            slopes = self.slopes(relative=True)

            from flatsurf.geometry.euclidean import is_parallel

            return all(
                is_parallel(slopes[i - 1], slopes[i]) for i in range(len(slopes))
            )

        def plot(
            self,
            translation=None,
            polygon_options={},
            edge_options={},
            vertex_options={},
        ):
            r"""
            Return a plot of this polygon with the origin at ``translation``.

            EXAMPLES:

            .. jupyter-execute::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.plot()
                ...Graphics object consisting of 3 graphics primitives

            We can specify an explicit ``zorder`` to render edges and vertices on
            top of the axes which are rendered at z-order 3:

            .. jupyter-execute::

                sage: S.plot(edge_options={'zorder': 3}, vertex_options={'zorder': 3})
                ...Graphics object consisting of 3 graphics primitives

            We can control the colors, e.g., we can render transparent polygons,
            with red edges and blue vertices:

            .. jupyter-execute::

                sage: S.plot(polygon_options={'fill': None}, edge_options={'color': 'red'}, vertex_options={'color': 'blue'})
                ...Graphics object consisting of 3 graphics primitives

            """
            from sage.plot.point import point2d
            from sage.plot.line import line2d
            from sage.plot.polygon import polygon2d

            P = self.vertices(translation)

            polygon_options = {"alpha": 0.3, "zorder": 1, **polygon_options}
            edge_options = {"color": "orange", "zorder": 2, **edge_options}
            vertex_options = {"color": "red", "zorder": 2, **vertex_options}

            return (
                polygon2d(P, **polygon_options)
                + line2d(P + (P[0],), **edge_options)
                + point2d(P, **vertex_options)
            )

        def angles(self, numerical=None, assume_rational=None):
            r"""
            Return the list of angles of this polygon (divided by `2 \pi`).

            EXAMPLES::

                sage: from flatsurf import Polygon

                sage: T = Polygon(angles=[1, 2, 3])
                sage: [T.angle(i) for i in range(3)]
                [1/12, 1/6, 1/4]
                sage: T.angles()
                (1/12, 1/6, 1/4)
                sage: sum(T.angle(i) for i in range(3))
                1/2
            """
            if assume_rational is not None:
                import warnings

                warnings.warn(
                    "assume_rational has been deprecated as a keyword to angles() and will be removed from a future version of sage-flatsurf"
                )

            angles = tuple(
                self.angle(i, numerical=numerical) for i in range(len(self.vertices()))
            )

            if not numerical:
                self._refine_category_(self.category().WithAngles(angles))

            return angles

        def angle(self, e, numerical=None, assume_rational=None):
            r"""
            Return the angle at the beginning of the start point of the edge ``e``.

            EXAMPLES::

                sage: from flatsurf.geometry.polygon import polygons
                sage: polygons.square().angle(0)
                1/4
                sage: polygons.regular_ngon(8).angle(0)
                3/8

                sage: from flatsurf import Polygon
                sage: T = Polygon(vertices=[(0,0), (3,1), (1,5)])
                sage: [T.angle(i, numerical=True) for i in range(3)]  # abs tol 1e-13
                [0.16737532973071603, 0.22741638234956674, 0.10520828791971722]
                sage: sum(T.angle(i, numerical=True) for i in range(3))   # abs tol 1e-13
                0.5
            """
            if assume_rational is not None:
                import warnings

                warnings.warn(
                    "assume_rational has been deprecated as a keyword to angle() and will be removed from a future version of sage-flatsurf"
                )

            if numerical is None:
                numerical = not self.is_rational()

                if numerical:
                    import warnings

                    warnings.warn(
                        "the behavior of angle() has been changed in recent versions of sage-flatsurf; for non-rational polygons, numerical=True must be set explicitly to get a numerical approximation of the angle"
                    )

            from flatsurf.geometry.euclidean import angle

            return angle(
                self.edge(e),
                -self.edge((e - 1) % len(self.vertices())),
                numerical=numerical,
            )

        def area(self):
            r"""
            Return the area of this polygon.

            EXAMPLES::

                sage: from flatsurf.geometry.polygon import polygons
                sage: polygons.regular_ngon(8).area()
                2*a + 2
                sage: _ == 2*AA(2).sqrt() + 2
                True

                sage: AA(polygons.regular_ngon(11).area())
                9.36563990694544?

                sage: polygons.square().area()
                1
                sage: (2*polygons.square()).area()
                4
            """
            # Will use an area formula obtainable from Green's theorem. See for instance:
            # http://math.blogoverflow.com/2014/06/04/greens-theorem-and-area-of-polygons/
            total = self.base_ring().zero()
            for i in range(len(self.vertices())):
                total += (self.vertex(i)[0] + self.vertex(i + 1)[0]) * self.edge(i)[1]

            from sage.all import ZZ

            return total / ZZ(2)

        def centroid(self):
            r"""
            Return the coordinates of the centroid of this polygon.

            ALGORITHM:

            We use the customary formula of the centroid of polygons, see
            https://en.wikipedia.org/wiki/Centroid#Of_a_polygon

            EXAMPLES::

                sage: from flatsurf.geometry.polygon import polygons
                sage: P = polygons.regular_ngon(4)
                sage: P
                Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                sage: P.centroid()
                (1/2, 1/2)

                sage: P = polygons.regular_ngon(8); P
                Polygon(vertices=[(0, 0), (1, 0), (1/2*a + 1, 1/2*a), (1/2*a + 1, 1/2*a + 1), (1, a + 1), (0, a + 1), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a)])
                sage: P.centroid()
                (1/2, 1/2*a + 1/2)

                sage: P = polygons.regular_ngon(11)
                sage: C = P.centroid()
                sage: P = P.translate(-C)
                sage: P.centroid()
                (0, 0)

            """
            x, y = list(zip(*self.vertices()))
            nvertices = len(x)
            A = self.area()

            from sage.all import vector

            return vector(
                (
                    ~(6 * A)
                    * sum(
                        [
                            (x[i - 1] + x[i]) * (x[i - 1] * y[i] - x[i] * y[i - 1])
                            for i in range(nvertices)
                        ]
                    ),
                    ~(6 * A)
                    * sum(
                        [
                            (y[i - 1] + y[i]) * (x[i - 1] * y[i] - x[i] * y[i - 1])
                            for i in range(nvertices)
                        ]
                    ),
                )
            )

        def get_point_position(self, point, translation=None):
            r"""
            Return the combinatorial classification of a point and a polygon.

            INPUT:

            - ``point`` -- a point in the plane as a SageMath vector or pair of
              numbers

            OUTPUT:

            a :class:`.geometry.polygon.PolygonPosition` object

            ALGORITHM:

            We use a winding number algorithm see
            https://en.wikipedia.org/wiki/Point_in_polygon#Winding_number_algorithm

            EXAMPLES::

                sage: from flatsurf import polygons, Polygon
                sage: S = polygons.square()

                sage: S.get_point_position((1/2, 1/2))
                point positioned in interior of polygon
                sage: S.get_point_position((1, 0))
                point positioned on vertex 1 of polygon
                sage: S.get_point_position((1, 1/2))
                point positioned on interior of edge 1 of polygon
                sage: S.get_point_position((1, 3/2))
                point positioned outside polygon

                sage: p = Polygon(edges=[(1, 0), (1, 0), (1, 0), (0, 1), (-3, 0), (0, -1)])
                sage: p.get_point_position([10, 0])
                point positioned outside polygon
                sage: p.get_point_position([1/2, 0])
                point positioned on interior of edge 0 of polygon
                sage: p.get_point_position([3/2, 0])
                point positioned on interior of edge 1 of polygon
                sage: p.get_point_position([2,0])
                point positioned on vertex 2 of polygon
                sage: p.get_point_position([5/2, 0])
                point positioned on interior of edge 2 of polygon
                sage: p.get_point_position([5/2, 1/4])
                point positioned in interior of polygon

            """
            from sage.all import vector

            point = vector(point)

            if translation is not None:
                import warnings

                warnings.warn(
                    "the translation keyword argument to get_point_position() has been deprecated and will be removed in a future version of sage-flatsurf; shift the point instead with the - operator"
                )

                from sage.all import vector

                return self.get_point_position(point - vector(translation))

            from flatsurf.geometry.euclidean import ccw
            from flatsurf.geometry.polygon import PolygonPosition

            # Determine whether the point is a vertex of the polygon.
            for i, v in enumerate(self.vertices()):
                if point == v:
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=i)

            # Determine whether the point is on an edge of the polygon.
            for i, (v, e) in enumerate(zip(self.vertices(), self.edges())):
                if ccw(e, point - v) == 0:
                    # The point lies on the line through this edge.
                    if 0 < e.dot_product(point - v) < e.dot_product(e):
                        return PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)

            # Determine whether the point is inside or outside by computing the
            # winding number of the polygon.
            winding_number = 0

            for v, w in zip(self.vertices(), self.vertices()[1:] + self.vertices()[:1]):
                if v[1] < point[1] and w[1] >= point[1] and ccw(w - v, point - v) > 0:
                    winding_number += 1
                if v[1] >= point[1] and w[1] < point[1] and ccw(w - v, point - v) < 0:
                    winding_number -= 1

            if winding_number % 2:
                return PolygonPosition(PolygonPosition.INTERIOR)

            return PolygonPosition(PolygonPosition.OUTSIDE)

        def join(self, other, edge, other_edge):
            r"""
            Return the polygon obtained by gluing this polygon and ``other``
            along their ``edge`` and ``other_edge``, respectively.

            The polygons have to be such that the glued edges are identical but
            with opposite orientation.

            INPUT:

            - ``other`` -- a polygon over the same base ring as this polygon

            - ``edge`` -- an integer; the index of the edge of this polygon
              along which to glue

            - ``other_edge`` -- an integer; the index of the edge of ``other``
              along which to glue

            EXAMPLES::

                sage: from flatsurf import Polygon
                sage: P = Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
                sage: Q = Polygon(vertices=[(1, 0), (1, 1), (0, 1)])
                sage: P.join(Q, 1, 2)
                Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

                sage: P.join(P, 1, 1)
                Traceback (most recent call last):
                ...
                ValueError: glued edges must be identical with opposite orientation

                sage: Q = Polygon(vertices=[(0, 0), (0, 1), (-1, 1)])
                sage: P.join(Q, 1, 2)
                Traceback (most recent call last):
                ...
                ValueError: glued edges must be identical with opposite orientation

                sage: P.join(Q, 2, 0)
                Polygon(vertices=[(0, 0), (1, 0), (0, 1), (-1, 1)])

            Polygons cannot be joined if that would lead to a self-intersecting
            polygon::

                sage: P = Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2), (1, 1)])
                sage: Q = Polygon(vertices=[(0, 0), (1, 1), (2, 2), (0, 2)])
                sage: P.join(Q, 4, 0)
                Traceback (most recent call last):
                ...
                NotImplementedError: polygon self-intersects

            """
            if self.vertex(edge) != other.vertex(other_edge + 1) or self.vertex(
                edge + 1
            ) != other.vertex(other_edge):
                raise ValueError(
                    "glued edges must be identical with opposite orientation"
                )

            from flatsurf import Polygon

            return Polygon(
                base_ring=self.base_ring(),
                vertices=self.vertices()[:edge]
                + other.vertices()[other_edge + 1 :]
                + other.vertices()[:other_edge]
                + self.vertices()[edge + 1 :],
            )

    class Rational(CategoryWithAxiom_over_base_ring):
        r"""
        The category of rational Euclidean polygons.

        .. NOTE::

            This category must be defined here to make SageMath's test suite
            pass. Otherwise we get "The super categories of a category over
            base should be a category over base (or the related Bimodules) or a
            singleton category"; we did not investigate what exactly is going on
            here.

        """

    class SubcategoryMethods:
        @cached_method
        def module(self):
            r"""
            Return the free module of rank 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.module()
                doctest:warning
                ...
                UserWarning: module() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring()**2 instead
                Vector space of dimension 2 over Rational Field

            ::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: C = EuclideanPolygonsWithAngles(1, 2, 3)
                sage: C.module()
                Vector space of dimension 2 over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?

            """
            import warnings

            warnings.warn(
                "module() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring()**2 instead"
            )

            return FreeModule(self.base_ring(), 2)

        @cached_method
        def vector_space(self):
            r"""
            Return the vector space of dimension 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.vector_space()
                Vector space of dimension 2 over Rational Field

            ::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: C = EuclideanPolygonsWithAngles(1, 2, 3)
                sage: C.vector_space()
                doctest:warning
                ...
                UserWarning: vector_space() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring().fraction_field()**2 instead
                Vector space of dimension 2 over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?

            """
            import warnings

            warnings.warn(
                "vector_space() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring().fraction_field()**2 instead"
            )

            from sage.all import VectorSpace

            return VectorSpace(self.base_ring().fraction_field(), 2)

        def WithAngles(self, angles):
            r"""
            Return the subcategory of polygons with fixed ``angles``.

            INPUT:

            - ``angles`` -- a finite sequence of numbers, the inner angles of
              the polygon; the angles are automatically normalized to sum to
              `(n-2)π`.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import EuclideanPolygons
                sage: EuclideanPolygons(AA).WithAngles([1, 2, 3])
                Category of euclidean triangles with angles (1/12, 1/6, 1/4) over Algebraic Real Field

            """
            from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
                EuclideanPolygonsWithAngles,
            )

            angles = EuclideanPolygonsWithAngles._normalize_angles(angles)
            return EuclideanPolygonsWithAngles(self.base_ring(), angles) & self

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf import Polygons, ConvexPolygons

            sage: C = Polygons(QQ)
            sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
            doctest:warning
            ...
            UserWarning: Polygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead
            sage: p
            Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: C(p) is p
            False
            sage: C(p) == p
            True
            sage: C((1,0), (0,1), (-1, 1))
            Traceback (most recent call last):
            ...
            ValueError: the polygon does not close up

            sage: D = ConvexPolygons(QQbar)
            doctest:warning
            ...
            UserWarning: ConvexPolygons() has been deprecated and will be removed from a future version of sage-flatsurf; use Polygon() to create polygons.
            If you really need the category of convex polygons over a ring use EuclideanPolygons(ring).Simple().Convex() instead.
            sage: D(p)
            doctest:warning
            ...
            UserWarning: ConvexPolygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead
            Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: D(vertices=p.vertices())
            Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: D(edges=p.edges())
            Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
        """
        # We cannot have a __call__() in SubcategoryMethods so there is no good
        # way to support this in the category framework. Also, this code is
        # duplicated in several places and the Polygon() helper seems to be
        # much more versatile.
        import warnings

        warnings.warn(
            "Polygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead"
        )

        check = kwds.pop("check", True)

        from flatsurf.geometry.polygon import EuclideanPolygon

        if len(args) == 1 and isinstance(args[0], EuclideanPolygon):
            if args[0].category() is self:
                return args[0]
            vertices = [self.vector_space()(v) for v in args[0].vertices()]
            args = ()

        else:
            vertices = kwds.pop("vertices", None)
            edges = kwds.pop("edges", None)
            base_point = kwds.pop("base_point", (0, 0))

            if (vertices is None) and (edges is None):
                if len(args) == 1:
                    edges = args[0]
                elif args:
                    edges = args
                else:
                    raise ValueError(
                        "exactly one of 'vertices' or 'edges' must be provided"
                    )
            if kwds:
                raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

            if edges is not None:
                v = self.vector_space()(base_point)
                vertices = []
                for e in map(self.vector_space(), edges):
                    vertices.append(v)
                    v += e
                if v != vertices[0]:
                    raise ValueError("the polygon does not close up")

        from flatsurf.geometry.polygon import Polygon

        return Polygon(
            base_ring=self.base(), vertices=vertices, category=self, check=check
        )

    class Convex(CategoryWithAxiom_over_base_ring):
        r"""
        The subcategory of convex Euclidean polygons in the real plane.

        EXAMPLES:

        For historic reasons, there is the shortcut ``ConvexPolygons`` to get
        the Euclidean convex polygons::

            sage: from flatsurf import ConvexPolygons
            sage: C = ConvexPolygons(QQ)

            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: C is EuclideanPolygons(QQ).Convex().Simple()
            True

            sage: C(vertices=[(0,0), (2,0), (1,1)])
            Polygon(vertices=[(0, 0), (2, 0), (1, 1)])

            sage: C(edges=[(1,0), (0,1), (-1,0), (0,-1)])
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        This axiom can also be created over non-fields::

            sage: ConvexPolygons(ZZ)
            Category of convex simple euclidean polygons over Integer Ring

        TESTS::

            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: TestSuite(EuclideanPolygons(QQ).Convex()).run()
            sage: TestSuite(EuclideanPolygons(QQbar).Convex()).run()
            sage: TestSuite(EuclideanPolygons(ZZ).Convex()).run()

        """

        class ParentMethods:
            r"""
            Provides methods available to all convex Euclidean polygons in the
            real plane.

            If you want to add functionality to such polygons, you probably
            want to put it here.
            """

            def is_convex(self, strict=False):
                r"""
                Return whether this is a convex polygon.

                INPUT:

                - ``strict`` -- whether to check for strict convexity, i.e., a
                  polygon with a π angle is not considered convex.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: S = polygons.square()
                    sage: S.is_convex()
                    True
                    sage: S.is_convex(strict=True)
                    True

                """
                if not strict:
                    return True

                return EuclideanPolygons.ParentMethods.is_convex(self, strict=strict)

    class Simple(CategoryWithAxiom_over_base_ring):
        r"""
        The subcategory of Euclidean polygons without self-intersection in the
        real plane.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: EuclideanPolygons(QQ).Simple()
            Category of simple euclidean polygons over Rational Field

        """

        class ParentMethods:
            r"""
            Provides methods available to all simple Euclidean polygons.

            If you want to add functionality to all polygons, independent of
            implementation, you probably want to put it here.
            """

            def triangulation(self):
                r"""
                Return a list of pairs of indices of vertices that together with the boundary
                form a triangulation.

                EXAMPLES:

                We triangulate a non-convex polygon::

                    sage: from flatsurf import Polygon
                    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1), (0,2), (-1,2), (-1,1), (-2,1),
                    ....:                    (-2,0), (-1,0), (-1,-1), (0,-1)])
                    sage: P.triangulation()
                    [(0, 2), (2, 8), (3, 5), (6, 8), (8, 3), (3, 6), (9, 11), (0, 9), (2, 9)]

                TESTS::

                    sage: Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)]).triangulation()
                    [(0, 2)]

                    sage: quad = [(0,0), (1,-1), (0,1), (-1,-1)]
                    sage: for i in range(4):
                    ....:     Polygon(vertices=quad[i:] + quad[:i]).triangulation()
                    [(0, 2)]
                    [(1, 3)]
                    [(0, 2)]
                    [(1, 3)]

                    sage: poly = [(0,0),(1,1),(2,0),(3,1),(4,0),(4,2),
                    ....:     (-4,2),(-4,0),(-3,1),(-2,0),(-1,1)]
                    sage: Polygon(vertices=poly).triangulation()
                    [(1, 3), (3, 5), (5, 8), (6, 8), (8, 10), (10, 1), (1, 5), (5, 10)]
                    sage: for i in range(len(poly)):
                    ....:     Polygon(vertices=poly[i:] + poly[:i]).triangulation()
                    [(1, 3), (3, 5), (5, 8), (6, 8), (8, 10), (10, 1), (1, 5), (5, 10)]
                    [(0, 2), (2, 4), (4, 7), (5, 7), (7, 9), (9, 0), (0, 4), (4, 9)]
                    [(1, 3), (3, 6), (4, 6), (6, 8), (8, 10), (10, 1), (3, 8), (10, 3)]
                    [(0, 2), (2, 5), (3, 5), (5, 7), (7, 9), (9, 0), (2, 7), (9, 2)]
                    [(1, 4), (2, 4), (4, 6), (6, 8), (8, 10), (10, 1), (1, 6), (8, 1)]
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 9), (9, 0), (0, 5), (7, 0)]
                    [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10), (10, 2), (4, 10), (6, 10)]
                    [(1, 3), (3, 5), (5, 7), (7, 9), (9, 1), (10, 1), (3, 9), (5, 9)]
                    [(0, 2), (2, 4), (4, 6), (6, 8), (8, 0), (9, 0), (2, 8), (4, 8)]
                    [(1, 3), (3, 5), (5, 7), (7, 10), (8, 10), (10, 1), (1, 7), (3, 7)]
                    [(0, 2), (2, 4), (4, 6), (6, 9), (7, 9), (9, 0), (0, 6), (2, 6)]

                    sage: poly = [(0,0), (1,0), (2,0), (2,1), (2,2), (1,2), (0,2), (0,1)]
                    sage: Polygon(vertices=poly).triangulation()
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    sage: for i in range(len(poly)):
                    ....:     Polygon(vertices=poly[i:] + poly[:i]).triangulation()
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]

                    sage: poly = [(0,0), (1,2), (3,3), (1,4), (0,6), (-1,4), (-3,-3), (-1,2)]
                    sage: Polygon(vertices=poly).triangulation()
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    sage: for i in range(len(poly)):
                    ....:     Polygon(vertices=poly[i:] + poly[:i]).triangulation()
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]
                    [(0, 2), (3, 5), (5, 7), (7, 3), (0, 3)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]
                    [(0, 6), (1, 3), (3, 5), (5, 1), (6, 1)]
                    [(0, 2), (2, 4), (4, 6), (6, 0), (0, 4)]

                    sage: x = polygen(QQ)
                    sage: p = x^4 - 5*x^2 + 5
                    sage: r = AA.polynomial_root(p, RIF(1.17,1.18))
                    sage: K.<a> = NumberField(p, embedding=r)

                    sage: poly = [(1/2*a^2 - 3/2, 1/2*a),
                    ....:  (-a^3 + 2*a^2 + 2*a - 4, 0),
                    ....:  (1/2*a^2 - 3/2, -1/2*a),
                    ....:   (1/2*a^3 - a^2 - 1/2*a + 1, 1/2*a^2 - a),
                    ....:   (-1/2*a^2 + 1, 1/2*a^3 - 3/2*a),
                    ....:   (-1/2*a + 1, a^3 - 3/2*a^2 - 2*a + 5/2),
                    ....:   (1, 0),
                    ....:   (-1/2*a + 1, -a^3 + 3/2*a^2 + 2*a - 5/2),
                    ....:   (-1/2*a^2 + 1, -1/2*a^3 + 3/2*a),
                    ....:   (1/2*a^3 - a^2 - 1/2*a + 1, -1/2*a^2 + a)]
                    sage: Polygon(vertices=poly).triangulation()
                    [(0, 3), (1, 3), (3, 5), (5, 7), (7, 9), (9, 3), (3, 7)]

                    sage: z = QQbar.zeta(24)
                    sage: pts = [(1+i%2) * z**i for i in range(24)]
                    sage: pts = [vector(AA, (x.real(), x.imag())) for x in pts]
                    sage: Polygon(vertices=pts).triangulation()
                    [(0, 2), ..., (16, 0)]

                This is https://github.com/flatsurf/sage-flatsurf/issues/87 ::

                    sage: x = polygen(QQ)
                    sage: K.<c> = NumberField(x^2 - 3, embedding=AA(3).sqrt())

                    sage: Polygon(vertices=[(0, 0), (1, 0), (1/2*c + 1, -1/2), (c + 1, 0), (-3/2*c + 1, 5/2), (0, c - 2)]).triangulation()
                    [(0, 4), (1, 3), (4, 1)]

                """
                vertices = self.vertices()

                n = len(vertices)

                if n < 3:
                    raise ValueError

                if n == 3:
                    return []

                # NOTE: The algorithm is naive. We look at all possible chords between
                # the i-th and j-th vertices. If the chord does not intersect any edge
                # then we cut the polygon along this edge and call recursively
                # triangulate on the two pieces.
                for i in range(n - 1):
                    eiright = vertices[(i + 1) % n] - vertices[i]
                    eileft = vertices[(i - 1) % n] - vertices[i]
                    for j in range(i + 2, (n if i else n - 1)):
                        ejright = vertices[(j + 1) % n] - vertices[j]
                        ejleft = vertices[(j - 1) % n] - vertices[j]
                        chord = vertices[j] - vertices[i]

                        from flatsurf.geometry.euclidean import is_between

                        # check angles with neighbouring edges
                        if not (
                            is_between(eiright, eileft, chord)
                            and is_between(ejright, ejleft, -chord)
                        ):
                            continue

                        # check intersection with other edges
                        e = (vertices[i], vertices[j])
                        good = True
                        for k in range(n):
                            if k == i or k == j or k == (i - 1) % n or k == (j - 1) % n:
                                continue

                            f = (vertices[k], vertices[(k + 1) % n])

                            from flatsurf.geometry.euclidean import (
                                is_segment_intersecting,
                            )

                            res = is_segment_intersecting(e, f)

                            assert res != 1

                            if res == 2:
                                good = False
                                break

                        if good:
                            from flatsurf import Polygon

                            part0 = [
                                (s + i, t + i)
                                for s, t in Polygon(
                                    vertices=vertices[i : j + 1], check=False
                                ).triangulation()
                            ]
                            part1 = []
                            for s, t in Polygon(
                                vertices=vertices[j:] + vertices[: i + 1], check=False
                            ).triangulation():
                                if s < n - j:
                                    s += j
                                else:
                                    s -= n - j
                                if t < n - j:
                                    t += j
                                else:
                                    t -= n - j
                                part1.append((s, t))
                            return [(i, j)] + part0 + part1

                assert False

            def triangulate(self):
                r"""
                Return a triangulation of this polygon.

                Returns a pair consisting of a surface and a bidict. The
                surfaces consists of the triangles that form the triangulation,
                glued to the original polygon. The bidict provides the mapping
                from the edge index of the polygon to the (label, edge) of the
                corresponding unglued edge in the surface.

                ALGORITHM:

                We use a simple quadratic time ear-clipping algorithm. There
                are asymptotically faster algorithms out there. But we are
                usually dealing with a very small number of vertices, so the
                asymptotic behavior does not seem to be the limiting factor
                here.

                EXAMPLES:

                We triangulate a non-convex polygon::

                    sage: from flatsurf import Polygon
                    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1), (0,2), (-1,2), (-1,1), (-2,1),
                    ....:                    (-2,0), (-1,0), (-1,-1), (0,-1)])
                    sage: P.triangulate()
                    (Translation Surface with boundary built from 6 isosceles triangles and 4 triangles,
                     bidict({0: (0, 0), 3: (2, 0), 7: (6, 0), 1: (0, 1), 2: (1, 1), 4: (2, 1), 5: (4, 1), 6: (5, 1), 8: (6, 1), 9: (8, 1), 10: (9, 1), 11: (9, 2)}))

                TESTS::

                    sage: _ = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)]).triangulate()

                    sage: quad = [(0,0), (1,-1), (0,1), (-1,-1)]
                    sage: for i in range(4):
                    ....:     _ = Polygon(vertices=quad[i:] + quad[:i]).triangulate()

                    sage: poly = [(0,0),(1,1),(2,0),(3,1),(4,0),(4,2),
                    ....:     (-4,2),(-4,0),(-3,1),(-2,0),(-1,1)]
                    sage: _ = Polygon(vertices=poly).triangulate()
                    sage: for i in range(len(poly)):
                    ....:     _ = Polygon(vertices=poly[i:] + poly[:i]).triangulate()

                    sage: poly = [(0,0), (1,0), (2,0), (2,1), (2,2), (1,2), (0,2), (0,1)]
                    sage: _ = Polygon(vertices=poly).triangulate()
                    sage: for i in range(len(poly)):
                    ....:     _ = Polygon(vertices=poly[i:] + poly[:i]).triangulate()

                    sage: poly = [(0,0), (1,2), (3,3), (1,4), (0,6), (-1,4), (-3,-3), (-1,2)]
                    sage: _ = Polygon(vertices=poly).triangulate()
                    sage: for i in range(len(poly)):
                    ....:     _ = Polygon(vertices=poly[i:] + poly[:i]).triangulate()

                    sage: x = polygen(QQ)
                    sage: p = x^4 - 5*x^2 + 5
                    sage: r = AA.polynomial_root(p, RIF(1.17,1.18))
                    sage: K.<a> = NumberField(p, embedding=r)

                    sage: poly = [(1/2*a^2 - 3/2, 1/2*a),
                    ....:  (-a^3 + 2*a^2 + 2*a - 4, 0),
                    ....:  (1/2*a^2 - 3/2, -1/2*a),
                    ....:   (1/2*a^3 - a^2 - 1/2*a + 1, 1/2*a^2 - a),
                    ....:   (-1/2*a^2 + 1, 1/2*a^3 - 3/2*a),
                    ....:   (-1/2*a + 1, a^3 - 3/2*a^2 - 2*a + 5/2),
                    ....:   (1, 0),
                    ....:   (-1/2*a + 1, -a^3 + 3/2*a^2 + 2*a - 5/2),
                    ....:   (-1/2*a^2 + 1, -1/2*a^3 + 3/2*a),
                    ....:   (1/2*a^3 - a^2 - 1/2*a + 1, -1/2*a^2 + a)]
                    sage: _ = Polygon(vertices=poly).triangulate()

                    sage: z = QQbar.zeta(24)
                    sage: pts = [(1+i%2) * z**i for i in range(24)]
                    sage: pts = [vector(AA, (x.real(), x.imag())) for x in pts]
                    sage: _ = Polygon(vertices=pts).triangulate()

                This is https://github.com/flatsurf/sage-flatsurf/issues/87 ::

                    sage: x = polygen(QQ)
                    sage: K.<c> = NumberField(x^2 - 3, embedding=AA(3).sqrt())
                    sage: _ = Polygon(vertices=[(0, 0), (1, 0), (1/2*c + 1, -1/2), (c + 1, 0), (-3/2*c + 1, 5/2), (0, c - 2)]).triangulate()

                """
                from flatsurf import MutableOrientedSimilaritySurface

                triangulation = MutableOrientedSimilaritySurface(self.base_ring())

                vertices = list(self.vertices())
                nvertices = len(vertices)

                # The vertices of the polygon that have not been ear-clipped.
                untriangulated = list(range(len(vertices)))

                # Maps triples of polygon vertices to the label of the triangle in the resulting triangulation.
                triangles = {}

                def next_label():
                    return len(triangles)

                if len(untriangulated) < 3:
                    raise ValueError("cannot triangulate such a degenerate polygon")

                if len(untriangulated) == 3:
                    # No need to triangulate. We special-case so we get nicer labels.
                    label = triangulation.add_polygon(self)

                    from bidict import bidict

                    return triangulation, bidict(
                        {0: (label, 0), 1: (label, 1), 2: (label, 2)}
                    )

                while len(untriangulated) > 3:
                    for i in range(len(untriangulated)):
                        # We attempt to clip the j-th untriangulated vertex.
                        j = (i + 1) % len(untriangulated)
                        k = (j + 1) % len(untriangulated)

                        # a, b, c are the indexes of the vertices j-1, j, j+1 in the original polygon.
                        a = untriangulated[i]
                        b = untriangulated[j]
                        c = untriangulated[k]

                        if (
                            ccw(vertices[b] - vertices[a], vertices[c] - vertices[a])
                            <= 0
                        ):
                            # The triangle (a, b, c) has non-positive area.
                            continue

                        # Check that (a, b, c) form an ear, i.e., that there
                        # are no other vertices contained in the triangle.
                        if any(
                            ccw(vertices[b] - vertices[a], vertices[m] - vertices[a])
                            >= 0
                            and ccw(
                                vertices[c] - vertices[b], vertices[m] - vertices[b]
                            )
                            >= 0
                            and ccw(
                                vertices[a] - vertices[c], vertices[m] - vertices[c]
                            )
                            >= 0
                            for m in untriangulated
                            if m not in [a, b, c]
                        ):
                            continue

                        triangles[(a, b, c)] = next_label()
                        untriangulated = untriangulated[:j] + untriangulated[j + 1 :]
                        break
                    else:
                        assert False, "cannot triangulate this polygon"

                triangles[tuple(untriangulated)] = next_label()

                # Add triangles to triangulated surface.
                for triangle, label in triangles.items():
                    from flatsurf import Polygon

                    triangulation.add_polygon(
                        Polygon(
                            vertices=[
                                vertices[triangle[0]],
                                vertices[triangle[1]],
                                vertices[triangle[2]],
                            ]
                        ),
                        label=label,
                    )

                # Establish gluings between triangles
                edges = {
                    **{
                        triangle[:2]: (label, 0)
                        for (triangle, label) in triangles.items()
                    },
                    **{
                        triangle[1:]: (label, 1)
                        for (triangle, label) in triangles.items()
                    },
                    **{
                        (triangle[2], triangle[0]): (label, 2)
                        for (triangle, label) in triangles.items()
                    },
                }

                outer_edges = {}

                for a, b in edges:
                    glued = (b, a) in edges
                    assert not glued == (
                        b == (a + 1) % nvertices or a == (b + 1) % nvertices
                    )
                    if glued:
                        triangulation.glue(edges[(a, b)], edges[(b, a)])
                    else:
                        outer_edges[a] = edges[(a, b)]

                from bidict import bidict

                return triangulation, bidict(outer_edges)

        class Convex(CategoryWithAxiom_over_base_ring):
            r"""
            The subcategory of the simple convex Euclidean polygons.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import EuclideanPolygons
                sage: EuclideanPolygons(QQ).Simple().Convex()
                Category of convex simple euclidean polygons over Rational Field

            """

            def __call__(self, *args, **kwds):
                r"""
                TESTS::

                    sage: from flatsurf.geometry.categories import EuclideanPolygons

                    sage: C = EuclideanPolygons(QQ).Convex().Simple()
                    sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
                    doctest:warning
                    ...
                    UserWarning: ConvexPolygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead
                    sage: p
                    Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                    sage: C(p) is p
                    True
                    sage: C((1,0), (0,1), (-1, 1))
                    Traceback (most recent call last):
                    ...
                    ValueError: the polygon does not close up

                    sage: D = EuclideanPolygons(QQbar).Convex().Simple()
                    sage: D(p)
                    Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                    sage: D(vertices=p.vertices())
                    Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                    sage: D(edges=p.edges())
                    Polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])

                """
                # We cannot have a __call__() in SubcategoryMethods so there is no good
                # way to support this in the category framework. Also, this code is
                # duplicated in several places and the Polygon() helper seems to be
                # much more versatile.
                import warnings

                warnings.warn(
                    "ConvexPolygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead"
                )

                check = kwds.pop("check", True)

                from flatsurf.geometry.polygon import EuclideanPolygon

                if len(args) == 1 and isinstance(args[0], EuclideanPolygon):
                    if args[0].category() is self:
                        return args[0]

                    vertices = [self.vector_space()(v) for v in args[0].vertices()]
                    args = ()

                else:
                    vertices = kwds.pop("vertices", None)
                    edges = kwds.pop("edges", None)
                    base_point = kwds.pop("base_point", (0, 0))

                    if (vertices is None) and (edges is None):
                        if len(args) == 1:
                            edges = args[0]
                        elif args:
                            edges = args
                        else:
                            raise ValueError(
                                "exactly one of 'vertices' or 'edges' must be provided"
                            )
                    if kwds:
                        raise ValueError(
                            "invalid keyword {!r}".format(next(iter(kwds)))
                        )

                    if edges is not None:
                        v = (self.base_ring() ** 2)(base_point)
                        vertices = []
                        for e in map(self.base_ring() ** 2, edges):
                            vertices.append(v)
                            v += e
                        if v != vertices[0]:
                            raise ValueError("the polygon does not close up")

                from flatsurf.geometry.polygon import Polygon

                return Polygon(
                    base_ring=self.base(), vertices=vertices, category=self, check=check
                )

            class ParentMethods:
                r"""
                Provides methods available to all simple convex Euclidean polygons.

                If you want to add functionality to all polygons, independent of
                implementation, you probably want to put it here.
                """

                def find_separatrix(self, direction=None, start_vertex=0):
                    r"""
                    Returns a pair (v,same) where v is a vertex and same is a boolean.
                    The provided parameter "direction" should be a non-zero vector with
                    two entries, or by default direction=(0,1).

                    A separatrix is a ray leaving a vertex and entering the polygon.

                    The vertex v will have a separatrix leaving it which is parallel to
                    direction. The returned value "same" answers the question if this separatrix
                    points in the same direction as "direction". There is a boundary case:
                    we allow the separatrix to be an edge if and only if traveling along
                    the sepatrix from the vertex would travel in a counter-clockwise
                    direction about the polygon.

                    The vertex returned is uniquely defined from the above if the polygon
                    is a triangle. Otherwise, we return the first vertex with this property
                    obtained by inspecting starting at start_vertex (defaults to 0) and
                    then moving in the counter-clockwise direction.

                    EXAMPLES::

                        sage: from flatsurf import polygons
                        sage: p=polygons.square()
                        sage: print(p.find_separatrix())
                        (1, True)
                        sage: print(p.find_separatrix(start_vertex=2))
                        (3, False)

                    """
                    if direction is None:
                        direction = (self.base_ring() ** 2)(
                            (self.base_ring().zero(), self.base_ring().one())
                        )
                    else:
                        assert not direction.is_zero()
                    v = start_vertex
                    n = len(self.vertices())
                    for i in range(len(self.vertices())):
                        if (
                            ccw(self.edge(v), direction) >= 0
                            and ccw(self.edge(v + n - 1), direction) > 0
                        ):
                            return v, True
                        if (
                            ccw(self.edge(v), direction) <= 0
                            and ccw(self.edge(v + n - 1), direction) < 0
                        ):
                            return v, False
                        v = v + 1 % n
                    raise RuntimeError("Failed to find a separatrix")

                def contains_point(self, point, translation=None):
                    r"""
                    Return whether the point is within the polygon (after the polygon is possibly translated)
                    """
                    return self.get_point_position(
                        point, translation=translation
                    ).is_inside()

                def flow_to_exit(self, point, direction):
                    r"""
                    Flow a point in the direction of holonomy until the point leaves the
                    polygon.  Note that ValueErrors may be thrown if the point is not in the
                    polygon, or if it is on the boundary and the holonomy does not point
                    into the polygon.

                    INPUT:

                    - ``point`` -- a point in the closure of the polygon (as a vector)

                    - ``holonomy`` -- direction of motion (a vector of non-zero length)

                    OUTPUT:

                    - The point in the boundary of the polygon where the trajectory exits

                    - a PolygonPosition object representing the combinatorial position of the stopping point
                    """
                    from flatsurf.geometry.polygon import PolygonPosition

                    V = self.base_ring().fraction_field() ** 2
                    if direction == V.zero():
                        raise ValueError("Zero vector provided as direction.")
                    v0 = self.vertex(0)
                    for i in range(len(self.vertices())):
                        e = self.edge(i)
                        from sage.all import matrix

                        m = matrix([[e[0], -direction[0]], [e[1], -direction[1]]])
                        try:
                            ret = m.inverse() * (point - v0)
                            s = ret[0]
                            t = ret[1]
                            # What if the matrix is non-invertible?

                            # Answer: You'll get a ZeroDivisionError which means that the edge is parallel
                            # to the direction.

                            # s is location it intersects on edge, t is the portion of the direction to reach this intersection
                            if t > 0 and 0 <= s and s <= 1:
                                # The ray passes through edge i.
                                if s == 1:
                                    # exits through vertex i+1
                                    v0 = v0 + e
                                    return v0, PolygonPosition(
                                        PolygonPosition.VERTEX,
                                        vertex=(i + 1) % len(self.vertices()),
                                    )
                                if s == 0:
                                    # exits through vertex i
                                    return v0, PolygonPosition(
                                        PolygonPosition.VERTEX, vertex=i
                                    )
                                    # exits through vertex i
                                # exits through interior of edge i
                                prod = t * direction
                                return point + prod, PolygonPosition(
                                    PolygonPosition.EDGE_INTERIOR, edge=i
                                )
                        except ZeroDivisionError:
                            # Here we know the edge and the direction are parallel
                            if ccw(e, point - v0) == 0:
                                # In this case point lies on the edge.
                                # We need to work out which direction to move in.
                                from flatsurf.geometry.euclidean import is_parallel

                                if (point - v0).is_zero() or is_parallel(e, point - v0):
                                    # exits through vertex i+1
                                    return self.vertex(i + 1), PolygonPosition(
                                        PolygonPosition.VERTEX,
                                        vertex=(i + 1) % len(self.vertices()),
                                    )
                                else:
                                    # exits through vertex i
                                    return v0, PolygonPosition(
                                        PolygonPosition.VERTEX, vertex=i
                                    )
                            pass
                        v0 = v0 + e
                    # Our loop has terminated. This can mean one of several errors...
                    pos = self.get_point_position(point)
                    if pos.is_outside():
                        raise ValueError("Started with point outside polygon")
                    raise ValueError(
                        "Point on boundary of polygon and direction not pointed into the polygon."
                    )

                def flow_map(self, direction):
                    r"""
                    Return a polygonal map associated to the flow in ``direction`` in this
                    polygon.

                    EXAMPLES::

                        sage: from flatsurf.geometry.polygon import Polygon
                        sage: S = Polygon(vertices=[(0,0),(2,0),(2,2),(1,2),(0,2),(0,1)])
                        sage: S.flow_map((0,1))
                         Flow polygon map:
                          3 2
                          0
                         top lengths: [1, 1]
                         bot lengths: [2]
                        sage: S.flow_map((1,1))
                        Flow polygon map:
                         3 2 1
                         4 5 0
                        top lengths: [1, 1, 2]
                        bot lengths: [1, 1, 2]
                        sage: S.flow_map((-1,-1))
                        Flow polygon map:
                         0 5 4
                         1 2 3
                        top lengths: [2, 1, 1]
                        bot lengths: [2, 1, 1]

                        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
                        sage: S.flow_map((sqrt2,1))
                        Flow polygon map:
                         3 2 1
                         4 5 0
                        top lengths: [1, 1, 2*sqrt2]
                        bot lengths: [sqrt2, sqrt2, 2]
                    """
                    from sage.all import vector

                    direction = vector(direction)
                    DP = direction.parent()
                    P = self.base_ring().fraction_field() ** 2
                    if DP != P:
                        P = cm.common_parent(DP, P)
                        ring = P.base_ring()
                        direction = direction.change_ring(ring)
                    else:
                        ring = P.base_ring()

                    rt = None
                    lt = None
                    lb = None
                    rb = None

                    # first compute the transversal length of each edge
                    t = P([direction[1], -direction[0]])
                    lengths = [t.dot_product(e) for e in self.edges()]
                    n = len(lengths)
                    for i in range(n):
                        j = (i + 1) % len(lengths)
                        l0 = lengths[i]
                        l1 = lengths[j]
                        if l0 >= 0 and l1 < 0:
                            rt = j
                        if l0 > 0 and l1 <= 0:
                            rb = j
                        if l0 <= 0 and l1 > 0:
                            lb = j
                        if l0 < 0 and l1 >= 0:
                            lt = j

                    assert rt is not None
                    assert lt is not None
                    assert lb is not None
                    assert rb is not None

                    if rt < lt:
                        top_lengths = lengths[rt:lt]
                        top_labels = list(range(rt, lt))
                    else:
                        top_lengths = lengths[rt:] + lengths[:lt]
                        top_labels = list(range(rt, n)) + list(range(lt))
                    top_lengths = [-x for x in reversed(top_lengths)]
                    top_labels.reverse()

                    if lb < rb:
                        bot_lengths = lengths[lb:rb]
                        bot_labels = list(range(lb, rb))
                    else:
                        bot_lengths = lengths[lb:] + lengths[:rb]
                        bot_labels = list(range(lb, n)) + list(range(rb))

                    from flatsurf.geometry.interval_exchange_transformation import (
                        FlowPolygonMap,
                    )

                    return FlowPolygonMap(
                        ring, bot_labels, bot_lengths, top_labels, top_lengths
                    )

                def flow(self, point, holonomy, translation=None):
                    r"""
                    Flow a point in the direction of holonomy for the length of the
                    holonomy, or until the point leaves the polygon.  Note that ValueErrors
                    may be thrown if the point is not in the polygon, or if it is on the
                    boundary and the holonomy does not point into the polygon.

                    INPUT:

                    - ``point`` -- a point in the closure of the polygon (vector over the underlying base_ring())

                    - ``holonomy`` -- direction and magnitude of motion (vector over the underlying base_ring())

                    - ``translation`` -- optional translation to applied to the polygon (vector over the underlying base_ring())

                    OUTPUT:

                    - The point within the polygon where the motion stops (or leaves the polygon)

                    - The amount of holonomy left to flow

                    - a PolygonPosition object representing the combinatorial position of the stopping point

                    EXAMPLES::

                        sage: from flatsurf.geometry.polygon import polygons
                        sage: s = polygons.square()
                        sage: V = QQ**2
                        sage: p = V((1/2, 1/2))
                        sage: w = V((2, 0))
                        sage: s.flow(p, w)
                        ((1, 1/2), (3/2, 0), point positioned on interior of edge 1 of polygon)
                    """
                    from flatsurf.geometry.polygon import PolygonPosition

                    V = self.base_ring().fraction_field() ** 2
                    if holonomy == V.zero():
                        # not flowing at all!
                        return (
                            point,
                            V.zero(),
                            self.get_point_position(point, translation=translation),
                        )
                    if translation is None:
                        v0 = self.vertex(0)
                    else:
                        v0 = self.vertex(0) + translation
                    for i in range(len(self.vertices())):
                        e = self.edge(i)
                        from sage.all import matrix

                        m = matrix([[e[0], -holonomy[0]], [e[1], -holonomy[1]]])
                        try:
                            ret = m.inverse() * (point - v0)
                            s = ret[0]
                            t = ret[1]
                            # What if the matrix is non-invertible?

                            # s is location it intersects on edge, t is the portion of the holonomy to reach this intersection
                            if t > 0 and 0 <= s and s <= 1:
                                # The ray passes through edge i.
                                if t > 1:
                                    # the segment from point with the given holonomy stays within the polygon
                                    return (
                                        point + holonomy,
                                        V.zero(),
                                        PolygonPosition(PolygonPosition.INTERIOR),
                                    )
                                if s == 1:
                                    # exits through vertex i+1
                                    v0 = v0 + e
                                    return (
                                        v0,
                                        point + holonomy - v0,
                                        PolygonPosition(
                                            PolygonPosition.VERTEX,
                                            vertex=(i + 1) % len(self.vertices()),
                                        ),
                                    )
                                if s == 0:
                                    # exits through vertex i
                                    return (
                                        v0,
                                        point + holonomy - v0,
                                        PolygonPosition(
                                            PolygonPosition.VERTEX, vertex=i
                                        ),
                                    )
                                    # exits through vertex i
                                # exits through interior of edge i
                                prod = t * holonomy
                                return (
                                    point + prod,
                                    holonomy - prod,
                                    PolygonPosition(
                                        PolygonPosition.EDGE_INTERIOR, edge=i
                                    ),
                                )
                        except ZeroDivisionError:
                            # can safely ignore this error. It means that the edge and the holonomy are parallel.
                            pass
                        v0 = v0 + e
                    # Our loop has terminated. This can mean one of several errors...
                    pos = self.get_point_position(point, translation=translation)
                    if pos.is_outside():
                        raise ValueError("Started with point outside polygon")
                    raise ValueError(
                        "Point on boundary of polygon and holonomy not pointed into the polygon."
                    )

                def circumscribing_circle(self):
                    import warnings

                    warnings.warn(
                        "circumscribing_circle() has been deprecated and will be removed from a future version of sage-flatsurf; use circumscribed_circle() instead"
                    )

                    return self.circumscribed_circle()

                def circumscribed_circle(self):
                    r"""
                    Returns the circle which circumscribes this polygon.
                    Raises a ValueError if the polygon is not circumscribed by a circle.

                    EXAMPLES::

                        sage: from flatsurf import Polygon
                        sage: P = Polygon(vertices=[(0,0),(1,0),(2,1),(-1,1)])
                        sage: P.circumscribed_circle()
                        Circle((1/2, 3/2), 5/2)
                    """
                    from flatsurf.geometry.circle import circle_from_three_points

                    circle = circle_from_three_points(
                        self.vertex(0), self.vertex(1), self.vertex(2), self.base_ring()
                    )
                    for i in range(3, len(self.vertices())):
                        if not circle.point_position(self.vertex(i)) == 0:
                            raise ValueError(
                                "Vertex " + str(i) + " is not on the circle."
                            )
                    return circle

                def subdivide(self):
                    r"""
                    Return a list of triangles that partition this polygon.

                    For each edge of the polygon one triangle is created that joins this
                    edge to the
                    :meth:`~.EuclideanPolygons.ParentMethods.centroid` of
                    this polygon.

                    EXAMPLES::

                        sage: from flatsurf import polygons
                        sage: P = polygons.regular_ngon(3); P
                        Polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2*a)])
                        sage: print(P.subdivide())
                        [Polygon(vertices=[(0, 0), (1, 0), (1/2, 1/6*a)]),
                         Polygon(vertices=[(1, 0), (1/2, 1/2*a), (1/2, 1/6*a)]),
                         Polygon(vertices=[(1/2, 1/2*a), (0, 0), (1/2, 1/6*a)])]

                    ::

                        sage: P = polygons.regular_ngon(4)
                        sage: print(P.subdivide())
                        [Polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2)]),
                         Polygon(vertices=[(1, 0), (1, 1), (1/2, 1/2)]),
                         Polygon(vertices=[(1, 1), (0, 1), (1/2, 1/2)]),
                         Polygon(vertices=[(0, 1), (0, 0), (1/2, 1/2)])]

                    Sometimes alternating with :meth:`subdivide_edges` can produce a more
                    uniform subdivision::

                        sage: P = polygons.regular_ngon(4)
                        sage: print(P.subdivide_edges(2).subdivide())
                        [Polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2)]),
                         Polygon(vertices=[(1/2, 0), (1, 0), (1/2, 1/2)]),
                         Polygon(vertices=[(1, 0), (1, 1/2), (1/2, 1/2)]),
                         Polygon(vertices=[(1, 1/2), (1, 1), (1/2, 1/2)]),
                         Polygon(vertices=[(1, 1), (1/2, 1), (1/2, 1/2)]),
                         Polygon(vertices=[(1/2, 1), (0, 1), (1/2, 1/2)]),
                         Polygon(vertices=[(0, 1), (0, 1/2), (1/2, 1/2)]),
                         Polygon(vertices=[(0, 1/2), (0, 0), (1/2, 1/2)])]

                    """
                    vertices = self.vertices()
                    center = self.centroid()
                    from flatsurf import Polygon

                    return [
                        Polygon(
                            vertices=(
                                vertices[i],
                                vertices[(i + 1) % len(vertices)],
                                center,
                            ),
                        )
                        for i in range(len(vertices))
                    ]

                def subdivide_edges(self, parts=2):
                    r"""
                    Return a copy of this polygon whose edges have been split into
                    ``parts`` equal parts each.

                    INPUT:

                    - ``parts`` -- a positive integer (default: 2)

                    EXAMPLES::

                        sage: from flatsurf import polygons
                        sage: P = polygons.regular_ngon(3); P
                        Polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2*a)])
                        sage: P.subdivide_edges(1) == P
                        True
                        sage: P.subdivide_edges(2)
                        Polygon(vertices=[(0, 0), (1/2, 0), (1, 0), (3/4, 1/4*a), (1/2, 1/2*a), (1/4, 1/4*a)])
                        sage: P.subdivide_edges(3)
                        Polygon(vertices=[(0, 0), (1/3, 0), (2/3, 0), (1, 0), (5/6, 1/6*a), (2/3, 1/3*a), (1/2, 1/2*a), (1/3, 1/3*a), (1/6, 1/6*a)])

                    """
                    if parts < 1:
                        raise ValueError("parts must be a positive integer")

                    steps = [e / parts for e in self.edges()]
                    from flatsurf import Polygon

                    return Polygon(edges=[e for e in steps for p in range(parts)])

                def j_invariant(self):
                    r"""
                    Return the Kenyon-Smille J-invariant of this polygon.

                    The base ring of the polygon must be a number field.

                    The output is a triple ``(Jxx, Jyy, Jxy)`` that corresponds
                    respectively to the Sah-Arnoux-Fathi invariant of the vertical flow,
                    the Sah-Arnoux-Fathi invariant of the horizontal flow and the `xy`-matrix.

                    EXAMPLES::

                        sage: from flatsurf import polygons

                        sage: polygons.right_triangle(1/3,1).j_invariant()
                        (
                                  [0 0]
                        (0), (0), [1 0]
                        )

                    The regular 8-gon::

                        sage: polygons.regular_ngon(8).j_invariant()
                        (
                                  [2 2]
                        (0), (0), [2 1]
                        )

                        (
                                       [  0 3/2]
                        (1/2), (-1/2), [3/2   0]
                        )

                    Some extra debugging::

                        sage: K.<a> = NumberField(x^3 - 2, embedding=AA(2)**(1/3))
                        sage: ux = 1 + a + a**2
                        sage: uy = -2/3 + a
                        sage: vx = 1/5 - a**2
                        sage: vy = a + 7/13*a**2

                        sage: from flatsurf import Polygon
                        sage: p = Polygon(edges=[(ux, uy), (vx,vy), (-ux-vx,-uy-vy)], base_ring=K)
                        sage: Jxx, Jyy, Jxy = p.j_invariant()

                        sage: from flatsurf.geometry.categories import EuclideanPolygons
                        sage: EuclideanPolygons.Simple.Convex.ParentMethods._wedge_product(ux.vector(), vx.vector()) == Jxx
                        True
                        sage: EuclideanPolygons.Simple.Convex.ParentMethods._wedge_product(uy.vector(), vy.vector()) == Jyy
                        True

                    """
                    from sage.all import QQ, matrix

                    if self.base_ring() is QQ:
                        raise NotImplementedError

                    K = self.base_ring()
                    try:
                        V, from_V, to_V = K.vector_space()
                    except (AttributeError, ValueError):
                        raise ValueError(
                            "the surface needs to be define over a number field"
                        )

                    dim = K.degree()
                    M = K ** (dim * (dim - 1) // 2)
                    Jxx = Jyy = M.zero()
                    Jxy = matrix(K, dim)
                    vertices = list(self.vertices())
                    vertices.append(vertices[0])

                    for i in range(len(vertices) - 1):
                        a = to_V(vertices[i][0])
                        b = to_V(vertices[i][1])
                        c = to_V(vertices[i + 1][0])
                        d = to_V(vertices[i + 1][1])
                        Jxx += self._wedge_product(a, c)
                        Jyy += self._wedge_product(b, d)
                        Jxy += self._tensor_product(a, d)
                        Jxy -= self._tensor_product(c, b)

                    return (Jxx, Jyy, Jxy)

                @staticmethod
                def _wedge_product(v, w):
                    r"""
                    Return the wedge product of ``v`` and ``w``.

                    This is a helper method for :meth:`j_invariant`.

                    EXAMPLES::

                        sage: from flatsurf.geometry.categories import EuclideanPolygons

                        sage: EuclideanPolygons.Simple.Convex.ParentMethods._wedge_product(vector((1, 2)), vector((1, 2)))
                        (0)
                        sage: EuclideanPolygons.Simple.Convex.ParentMethods._wedge_product(vector((1, 2)), vector((2, 1)))
                        (-3)

                        sage: EuclideanPolygons.Simple.Convex.ParentMethods._wedge_product(vector((1, 2, 3)), vector((2, 3, 4)))
                        (-1, -2, -1)

                    """
                    d = len(v)

                    assert len(w) == d

                    R = v.base_ring()

                    from sage.all import free_module_element

                    return free_module_element(
                        R,
                        d * (d - 1) // 2,
                        [
                            (v[i] * w[j] - v[j] * w[i])
                            for i in range(d - 1)
                            for j in range(i + 1, d)
                        ],
                    )

                @staticmethod
                def _tensor_product(u, v):
                    r"""
                    Return the tensor product of ``u`` and ``v``.

                    This is a helper method for :meth:`j_invariant`.

                    EXAMPLES::

                        sage: from flatsurf.geometry.categories import EuclideanPolygons
                        sage: EuclideanPolygons.Simple.Convex.ParentMethods._tensor_product(vector((2, 3, 5)), vector((7, 11, 13)))
                        [14 21 35]
                        [22 33 55]
                        [26 39 65]

                    """
                    from sage.all import vector

                    u = vector(u)
                    v = vector(v)

                    d = len(u)
                    R = u.base_ring()

                    assert len(u) == len(v) and v.base_ring() == R
                    from sage.all import matrix

                    return matrix(
                        R, d, [u[i] * v[j] for j in range(d) for i in range(d)]
                    )

                def is_isometric(self, other, certificate=False):
                    r"""
                    Return whether ``self`` and ``other`` are isometric convex polygons via an orientation
                    preserving isometry.

                    If ``certificate`` is set to ``True`` return also a pair ``(index, rotation)``
                    of an integer ``index`` and a matrix ``rotation`` such that the given rotation
                    matrix identifies this polygon with the other and the edges 0 in this polygon
                    is mapped to the edge ``index`` in the other.

                    .. TODO::

                        Implement ``is_linearly_equivalent`` and ``is_similar``.

                    EXAMPLES::

                        sage: from flatsurf import Polygon, polygons
                        sage: S = polygons.square()
                        sage: S.is_isometric(S)
                        True
                        sage: U = matrix(2,[0,-1,1,0]) * S
                        sage: U.is_isometric(S)
                        True

                        sage: x = polygen(QQ)
                        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2)**(1/2))
                        sage: S = S.change_ring(K)
                        sage: U = matrix(2, [sqrt2/2, -sqrt2/2, sqrt2/2, sqrt2/2]) * S
                        sage: U.is_isometric(S)
                        True

                        sage: U2 = Polygon(edges=[(1,0), (sqrt2/2, sqrt2/2), (-1,0), (-sqrt2/2, -sqrt2/2)])
                        sage: U2.is_isometric(U)
                        False
                        sage: U2.is_isometric(U, certificate=True)
                        (False, None)

                        sage: S = Polygon(edges=[(1,0), (sqrt2/2, 3), (-2,3), (-sqrt2/2+1, -6)])
                        sage: T = Polygon(edges=[(sqrt2/2,3), (-2,3), (-sqrt2/2+1, -6), (1,0)])
                        sage: isometric, cert = S.is_isometric(T, certificate=True)
                        sage: assert isometric
                        sage: shift, rot = cert
                        sage: Polygon(edges=[rot * S.edge((k + shift) % 4) for k in range(4)]).translate(T.vertex(0)) == T
                        True


                        sage: T = (matrix(2, [sqrt2/2, -sqrt2/2, sqrt2/2, sqrt2/2]) * S).translate((3,2))
                        sage: isometric, cert = S.is_isometric(T, certificate=True)
                        sage: assert isometric
                        sage: shift, rot = cert
                        sage: Polygon(edges=[rot * S.edge(k + shift) for k in range(4)]).translate(T.vertex(0)) == T
                        True
                    """
                    from flatsurf.geometry.polygon import EuclideanPolygon

                    if not isinstance(other, EuclideanPolygon):
                        raise TypeError("other must be a polygon")

                    if not other.is_convex():
                        raise TypeError("other must be convex")

                    n = len(self.vertices())
                    if len(other.vertices()) != n:
                        return False
                    sedges = self.edges()
                    oedges = other.edges()

                    slengths = [x**2 + y**2 for x, y in sedges]
                    olengths = [x**2 + y**2 for x, y in oedges]
                    for i in range(n):
                        if slengths == olengths:
                            # we have a match of lengths after a shift by i
                            xs, ys = sedges[0]
                            xo, yo = oedges[0]
                            from sage.all import matrix

                            ms = matrix(2, [xs, -ys, ys, xs])
                            mo = matrix(2, [xo, -yo, yo, xo])
                            rot = mo * ~ms
                            assert rot.det() == 1 and (rot * rot.transpose()).is_one()
                            assert oedges[0] == rot * sedges[0]
                            if all(oedges[i] == rot * sedges[i] for i in range(1, n)):
                                return (
                                    (True, (0 if i == 0 else n - i, rot))
                                    if certificate
                                    else True
                                )
                        olengths.append(olengths.pop(0))
                        oedges.append(oedges.pop(0))
                    return (False, None) if certificate else False

                def is_translate(self, other, certificate=False):
                    r"""
                    Return whether ``other`` is a translate of ``self``.

                    EXAMPLES::

                        sage: from flatsurf import Polygon
                        sage: S = Polygon(vertices=[(0,0), (3,0), (1,1)])
                        sage: T1 = S.translate((2,3))
                        sage: S.is_translate(T1)
                        True
                        sage: T2 = Polygon(vertices=[(-1,1), (1,0), (2,1)])
                        sage: S.is_translate(T2)
                        False
                        sage: T3 = Polygon(vertices=[(0,0), (3,0), (2,1)])
                        sage: S.is_translate(T3)
                        False

                        sage: S.is_translate(T1, certificate=True)
                        (True, (0, 1))
                        sage: S.is_translate(T2, certificate=True)
                        (False, None)
                        sage: S.is_translate(T3, certificate=True)
                        (False, None)
                    """
                    if type(self) is not type(other):
                        raise TypeError

                    n = len(self.vertices())
                    if len(other.vertices()) != n:
                        return False
                    sedges = self.edges()
                    oedges = other.edges()
                    for i in range(n):
                        if sedges == oedges:
                            return (True, (i, 1)) if certificate else True
                        oedges.append(oedges.pop(0))
                    return (False, None) if certificate else False

                def is_half_translate(self, other, certificate=False):
                    r"""
                    Return whether ``other`` is a translate or half-translate of ``self``.

                    If ``certificate`` is set to ``True`` then return also a pair ``(orientation, index)``.

                    EXAMPLES::

                        sage: from flatsurf import Polygon
                        sage: S = Polygon(vertices=[(0,0), (3,0), (1,1)])
                        sage: T1 = S.translate((2,3))
                        sage: S.is_half_translate(T1)
                        True
                        sage: T2 = Polygon(vertices=[(-1,1), (1,0), (2,1)])
                        sage: S.is_half_translate(T2)
                        True
                        sage: T3 = Polygon(vertices=[(0,0), (3,0), (2,1)])
                        sage: S.is_half_translate(T3)
                        False

                        sage: S.is_half_translate(T1, certificate=True)
                        (True, (0, 1))
                        sage: half_translate, cert = S.is_half_translate(T2, certificate=True)
                        sage: assert half_translate
                        sage: shift, rot = cert
                        sage: Polygon(edges=[rot * S.edge(k + shift) for k in range(3)]).translate(T2.vertex(0)) == T2
                        True
                        sage: S.is_half_translate(T3, certificate=True)
                        (False, None)
                    """
                    if type(self) is not type(other):
                        raise TypeError

                    n = len(self.vertices())
                    if len(other.vertices()) != n:
                        return False

                    sedges = self.edges()
                    oedges = other.edges()
                    for i in range(n):
                        if sedges == oedges:
                            return (True, (i, 1)) if certificate else True
                        oedges.append(oedges.pop(0))

                    assert oedges == other.edges()
                    oedges = [-e for e in oedges]
                    for i in range(n):
                        if sedges == oedges:
                            return (
                                (True, (0 if i == 0 else n - i, -1))
                                if certificate
                                else True
                            )
                        oedges.append(oedges.pop(0))

                    return (False, None) if certificate else False
