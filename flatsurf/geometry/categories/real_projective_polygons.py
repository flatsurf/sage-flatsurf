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
from sage.all import FreeModule
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.categories.polygons import Polygons
from flatsurf.geometry.matrix_2x2 import wedge_product

from sage.structure.element import get_coercion_model
cm = get_coercion_model()


class RealProjectivePolygons(Category_over_base_ring):
    def super_categories(self):
        return [Polygons(self.base_ring())]

    class ParentMethods:
        def vector_space(self):
            r"""
            Return the vector space of dimension 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.vector_space()
                Vector space of dimension 2 over Rational Field

            """
            return self.category().vector_space()

        def module(self):
            r"""
            Return the free module of rank 2 in which this polygon embeds.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.module()
                Vector space of dimension 2 over Rational Field

            """
            return self.category().module()

        def field(self):
            import warnings
            warnings.warn("field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() instead")

            return self.base_ring()

    class SubcategoryMethods:
        @cached_method
        def module(self):
            r"""
            Return the free module of rank 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.module()
                Vector space of dimension 2 over Rational Field

            """
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

            """
            from sage.all import VectorSpace
            return VectorSpace(self.base_ring().fraction_field(), 2)

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf import Polygons, ConvexPolygons

            sage: C = Polygons(QQ)
            sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
            sage: p
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: C(p) is p
            True
            sage: C((1,0), (0,1), (-1, 1))
            Traceback (most recent call last):
            ...
            ValueError: the polygon does not close up

            sage: D = ConvexPolygons(QQbar)
            sage: D(p)
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: D(vertices=p.vertices())
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: D(edges=p.edges())
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
        """
        # TODO: We should deprecate this. It will conflict with more complex categories eventually.
        check = kwds.pop("check", True)

        from flatsurf.geometry.polygon import Polygon
        if len(args) == 1 and isinstance(args[0], Polygon):
            if args[0].category() is self:
                return args[0]
            vertices = map(self.vector_space(), args[0].vertices())
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

        from flatsurf.geometry.polygon import EuclideanPolygon
        return EuclideanPolygon(ring=self.base(), vertices=vertices, category=self, check=check)

    class Convex(CategoryWithAxiom_over_base_ring):
        def __call__(self, *args, **kwds):
            r"""
            TESTS::

                sage: from flatsurf import ConvexPolygons

                sage: C = ConvexPolygons(QQ)
                sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
                sage: p
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                sage: C(p) is p
                True
                sage: C((1,0), (0,1), (-1, 1))
                Traceback (most recent call last):
                ...
                ValueError: the polygon does not close up

                sage: D = ConvexPolygons(QQbar)
                sage: D(p)
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                sage: D(vertices=p.vertices())
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                sage: D(edges=p.edges())
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])

            """
            # TODO: We should deprecate this. It will conflict with more complex categories eventually.
            check = kwds.pop("check", True)

            from flatsurf.geometry.polygon import Polygon
            if len(args) == 1 and isinstance(args[0], Polygon):
                if args[0].category() is self:
                    return args[0]

                vertices = map(self.vector_space(), args[0].vertices())
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
                    v = self.module()(base_point)
                    vertices = []
                    for e in map(self.module(), edges):
                        vertices.append(v)
                        v += e
                    if v != vertices[0]:
                        raise ValueError("the polygon does not close up")

            from flatsurf.geometry.polygon import EuclideanPolygon
            return EuclideanPolygon(ring=self.base(), vertices=vertices, category=self, check=check)

        class ParentMethods:
            def _check(self):
                r"""
                TESTS::

                    sage: from flatsurf import *
                    sage: polygons(vertices=[(0,0),(1,0)])
                    Traceback (most recent call last):
                    ...
                    ValueError: a polygon should have more than two edges!
                    sage: polygons(vertices=[(0,0),(1,2),(0,1),(-1,2)])
                    Traceback (most recent call last):
                    ...
                    ValueError: not convex
                    sage: polygons(vertices=[(0,0),(1,0),(2,0)])
                    Traceback (most recent call last):
                    ...
                    ValueError: degenerate polygon
                """
                if self.num_edges() <= 2:
                    raise ValueError("a polygon should have more than two edges!")

                if not sum(self.edges()).is_zero():
                    raise ValueError("the sum over the edges do not sum up to 0")

                for i in range(self.num_edges()):
                    if self.edge(i).is_zero():
                        raise ValueError("zero edge")
                    if wedge_product(self.edge(i), self.edge(i + 1)) < 0:
                        raise ValueError("not convex")
                    from flatsurf.geometry.matrix_2x2 import is_opposite_direction
                    if is_opposite_direction(self.edge(i), self.edge(i + 1)):
                        raise ValueError("degenerate polygon")

                Polygons.ParentMethods._check(self)

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
                    direction = self.module()((self.base_ring().zero(), self.base_ring().one()))
                else:
                    assert not direction.is_zero()
                v = start_vertex
                n = self.num_edges()
                zero = self.base_ring().zero()
                for i in range(self.num_edges()):
                    if (
                        wedge_product(self.edge(v), direction) >= zero
                        and wedge_product(self.edge(v + n - 1), direction) > zero
                    ):
                        return v, True
                    if (
                        wedge_product(self.edge(v), direction) <= zero
                        and wedge_product(self.edge(v + n - 1), direction) < zero
                    ):
                        return v, False
                    v = v + 1 % n
                raise RuntimeError("Failed to find a separatrix")

            def contains_point(self, point, translation=None):
                r"""
                Return true if the point is within the polygon (after the polygon is possibly translated)
                """
                return self.get_point_position(point, translation=translation).is_inside()

            def get_point_position(self, point, translation=None):
                r"""
                Get a combinatorial position of a points position compared to the polygon

                INPUT:

                - ``point`` -- a point in the plane (vector over the underlying base_ring())

                - ``translation`` -- optional translation to applied to the polygon (vector over the underlying base_ring())

                OUTPUT:

                - a PolygonPosition object

                EXAMPLES::

                    sage: from flatsurf.geometry.polygon import polygons
                    sage: s = polygons.square()
                    sage: V = s.parent().vector_space()
                    sage: s.get_point_position(V((1/2,1/2)))
                    point positioned in interior of polygon
                    sage: s.get_point_position(V((1,0)))
                    point positioned on vertex 1 of polygon
                    sage: s.get_point_position(V((1,1/2)))
                    point positioned on interior of edge 1 of polygon
                    sage: s.get_point_position(V((1,3/2)))
                    point positioned outside polygon

                    sage: p=polygons(edges=[(1,0),(1,0),(1,0),(0,1),(-3,0),(0,-1)])
                    sage: V=p.vector_space()
                    sage: p.get_point_position(V([10,0]))
                    point positioned outside polygon
                    sage: p.get_point_position(V([1/2,0]))
                    point positioned on interior of edge 0 of polygon
                    sage: p.get_point_position(V([3/2,0]))
                    point positioned on interior of edge 1 of polygon
                    sage: p.get_point_position(V([2,0]))
                    point positioned on vertex 2 of polygon
                    sage: p.get_point_position(V([5/2,0]))
                    point positioned on interior of edge 2 of polygon
                    sage: p.get_point_position(V([5/2,1/4]))
                    point positioned in interior of polygon
                """
                from flatsurf.geometry.polygon import PolygonPosition

                if translation is None:
                    # Since we allow the initial vertex to be non-zero, this changed:
                    v1 = self.vertex(0)
                else:
                    # Since we allow the initial vertex to be non-zero, this changed:
                    v1 = translation + self.vertex(0)
                # Below, we only make use of edge vectors:
                for i in range(self.num_edges()):
                    v0 = v1
                    e = self.edge(i)
                    v1 = v0 + e
                    w = wedge_product(e, point - v0)
                    if w < 0:
                        return PolygonPosition(PolygonPosition.OUTSIDE)
                    if w == 0:
                        # Lies on the line through edge i!
                        dp1 = e.dot_product(point - v0)
                        if dp1 == 0:
                            return PolygonPosition(PolygonPosition.VERTEX, vertex=i)
                        dp2 = e.dot_product(e)
                        if 0 < dp1 and dp1 < dp2:
                            return PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)
                # Loop terminated (on inside of each edge)
                return PolygonPosition(PolygonPosition.INTERIOR)

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

                V = self.parent().vector_space()
                if direction == V.zero():
                    raise ValueError("Zero vector provided as direction.")
                v0 = self.vertex(0)
                for i in range(self.num_edges()):
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
                                    PolygonPosition.VERTEX, vertex=(i + 1) % self.num_edges()
                                )
                            if s == 0:
                                # exits through vertex i
                                return v0, PolygonPosition(PolygonPosition.VERTEX, vertex=i)
                                # exits through vertex i
                            # exits through interior of edge i
                            prod = t * direction
                            return point + prod, PolygonPosition(
                                PolygonPosition.EDGE_INTERIOR, edge=i
                            )
                    except ZeroDivisionError:
                        # Here we know the edge and the direction are parallel
                        if wedge_product(e, point - v0) == 0:
                            # In this case point lies on the edge.
                            # We need to work out which direction to move in.
                            from flatsurf.geometry.matrix_2x2 import is_same_direction
                            if (point - v0).is_zero() or is_same_direction(e, point - v0):
                                # exits through vertex i+1
                                return self.vertex(i + 1), PolygonPosition(
                                    PolygonPosition.VERTEX, vertex=(i + 1) % self.num_edges()
                                )
                            else:
                                # exits through vertex i
                                return v0, PolygonPosition(PolygonPosition.VERTEX, vertex=i)
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

                    sage: from flatsurf.geometry.polygon import polygons
                    sage: S = polygons(vertices=[(0,0),(2,0),(2,2),(1,2),(0,2),(0,1)])
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
                P = self.vector_space()
                if DP != P:
                    P = cm.common_parent(DP, P)
                    ring = P.base_ring()
                    direction = direction.change_ring(ring)
                else:
                    ring = P.base_ring()

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

                from flatsurf.geometry.interval_exchange_transformation import FlowPolygonMap

                return FlowPolygonMap(ring, bot_labels, bot_lengths, top_labels, top_lengths)

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
                    sage: V = s.parent().vector_space()
                    sage: p = V((1/2,1/2))
                    sage: w = V((2,0))
                    sage: s.flow(p,w)
                    ((1, 1/2), (3/2, 0), point positioned on interior of edge 1 of polygon)
                """
                from flatsurf.geometry.polygon import PolygonPosition

                V = self.parent().vector_space()
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
                for i in range(self.num_edges()):
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
                                        vertex=(i + 1) % self.num_edges(),
                                    ),
                                )
                            if s == 0:
                                # exits through vertex i
                                return (
                                    v0,
                                    point + holonomy - v0,
                                    PolygonPosition(PolygonPosition.VERTEX, vertex=i),
                                )
                                # exits through vertex i
                            # exits through interior of edge i
                            prod = t * holonomy
                            return (
                                point + prod,
                                holonomy - prod,
                                PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i),
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
                r"""
                Returns the circle which circumscribes this polygon.
                Raises a ValueError if the polygon is not circumscribed by a circle.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: P = polygons(vertices=[(0,0),(1,0),(2,1),(-1,1)])
                    sage: P.circumscribing_circle()
                    Circle((1/2, 3/2), 5/2)
                """
                from flatsurf.geometry.circle import circle_from_three_points

                circle = circle_from_three_points(
                    self.vertex(0), self.vertex(1), self.vertex(2), self.base_ring()
                )
                for i in range(3, self.num_edges()):
                    if not circle.point_position(self.vertex(i)) == 0:
                        raise ValueError("Vertex " + str(i) + " is not on the circle.")
                return circle

            def subdivide(self):
                r"""
                Return a list of triangles that partition this polygon.

                For each edge of the polygon one triangle is created that joins this
                edge to the :meth:`centroid <Polygon.centroid>` of this polygon.

                EXAMPLES::

                    sage: from flatsurf import polygons
                    sage: P = polygons.regular_ngon(3); P
                    polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2*a)])
                    sage: P.subdivide()
                    [polygon(vertices=[(0, 0), (1, 0), (1/2, 1/6*a)]),
                     polygon(vertices=[(1, 0), (1/2, 1/2*a), (1/2, 1/6*a)]),
                     polygon(vertices=[(1/2, 1/2*a), (0, 0), (1/2, 1/6*a)])]

                ::

                    sage: P = polygons.regular_ngon(4)
                    sage: P.subdivide()
                    [polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2)]),
                     polygon(vertices=[(1, 0), (1, 1), (1/2, 1/2)]),
                     polygon(vertices=[(1, 1), (0, 1), (1/2, 1/2)]),
                     polygon(vertices=[(0, 1), (0, 0), (1/2, 1/2)])]

                Sometimes alternating with :meth:`subdivide_edges` can produce a more
                uniform subdivision::

                    sage: P = polygons.regular_ngon(4)
                    sage: P.subdivide_edges(2).subdivide()
                    [polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2)]),
                     polygon(vertices=[(1/2, 0), (1, 0), (1/2, 1/2)]),
                     polygon(vertices=[(1, 0), (1, 1/2), (1/2, 1/2)]),
                     polygon(vertices=[(1, 1/2), (1, 1), (1/2, 1/2)]),
                     polygon(vertices=[(1, 1), (1/2, 1), (1/2, 1/2)]),
                     polygon(vertices=[(1/2, 1), (0, 1), (1/2, 1/2)]),
                     polygon(vertices=[(0, 1), (0, 1/2), (1/2, 1/2)]),
                     polygon(vertices=[(0, 1/2), (0, 0), (1/2, 1/2)])]

                """
                vertices = self.vertices()
                center = self.centroid()
                return [
                    self.parent()(
                        vertices=(vertices[i], vertices[(i + 1) % len(vertices)], center)
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
                    polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2*a)])
                    sage: P.subdivide_edges(1) == P
                    True
                    sage: P.subdivide_edges(2)
                    polygon(vertices=[(0, 0), (1/2, 0), (1, 0), (3/4, 1/4*a), (1/2, 1/2*a), (1/4, 1/4*a)])
                    sage: P.subdivide_edges(3)
                    polygon(vertices=[(0, 0), (1/3, 0), (2/3, 0), (1, 0), (5/6, 1/6*a), (2/3, 1/3*a), (1/2, 1/2*a), (1/3, 1/3*a), (1/6, 1/6*a)])

                """
                if parts < 1:
                    raise ValueError("parts must be a positive integer")

                steps = [e / parts for e in self.edges()]
                return self.parent()(edges=[e for e in steps for p in range(parts)])
