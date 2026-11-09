r"""
The category of hyperbolic cone surfaces.

A hyperbolic cone surface is a surface that can be built by gluing compact
polygons of the hyperbolic plane along their edges with isometries. The
metric is hyperbolic away from the vertices, which can be cone points of any
total angle.

Two geodesic segments of the hyperbolic plane are related by an isometry only
when they have the same length, so there is no hyperbolic counterpart of the
gluings of edges of different lengths by similarities: the edges glued to each
other have the same length and the gluing is the unique orientation preserving
isometry that sends one edge onto the other one with reverse orientation.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

EXAMPLES:

A sphere with three cone points, glued from a Lambert quadrilateral::

    sage: from flatsurf import HyperbolicPlane
    sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
    sage: H = HyperbolicPlane(QQ)
    sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
    sage: S = MutableOrientedHyperbolicSurface(QQ)
    sage: S.add_polygon(P)
    0
    sage: S.glue((0, 0), (0, 3))
    sage: S.glue((0, 1), (0, 2))
    sage: S.set_immutable()

    sage: from flatsurf.geometry.categories import HyperbolicConeSurfaces
    sage: S in HyperbolicConeSurfaces()
    True

"""

# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2026 Vincent Delecroix
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)


class HyperbolicConeSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing compact polygons of the
    hyperbolic plane with isometries on the edges.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import HyperbolicConeSurfaces
        sage: HyperbolicConeSurfaces()
        Category of hyperbolic cone surfaces
        sage: HyperbolicConeSurfaces().Oriented().FiniteType()
        Category of finite type oriented hyperbolic cone surfaces

    """

    def super_categories(self):
        r"""
        Return the categories that a hyperbolic cone surface is also always a
        member of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import HyperbolicConeSurfaces
            sage: HyperbolicConeSurfaces().super_categories()
            [Category of cone surfaces]

        """
        from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces

        return [ConeSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all hyperbolic cone surfaces.

        If you want to add functionality for such surfaces you most likely
        want to put it here.
        """

        def hyperbolic_plane(self):
            r"""
            Return the hyperbolic plane that the polygons of this surface live
            in.

            EXAMPLES::

                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.hyperbolic_plane()
                Hyperbolic Plane over Rational Field

            """
            from flatsurf.geometry.hyperbolic import HyperbolicPlane

            return HyperbolicPlane(self.base_ring())

        def _polygon_coordinates(self, label, point, ring=None):
            r"""
            Return ``point`` in the coordinates of the polygon ``label``,
            i.e., its coordinates in the Klein model.

            This is a helper method for
            :class:`~flatsurf.geometry.surface_objects.SurfacePoint` which
            turns the coordinates provided by the user into the immutable
            vectors that the polygons of this surface work with.

            INPUT:

            - ``label`` -- a polygon label

            - ``point`` -- a point of the hyperbolic plane, or its coordinates
              in the Klein model as a pair or a vector

            - ``ring`` -- a ring or ``None`` (default: ``None``); the ring
              of the coordinates, the base ring of this surface if ``None``

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
                0
                sage: p = S._polygon_coordinates(0, (1/4, 0))
                sage: p
                (1/4, 0)
                sage: p.is_immutable()
                True
                sage: S._polygon_coordinates(0, H(2*I))
                (0, 3/5)
                sage: S._polygon_coordinates(0, (1/4, 0), ring=AA).parent()
                Vector space of dimension 2 over Algebraic Real Field

            """
            from flatsurf.geometry.hyperbolic import HyperbolicPoint

            if ring is None:
                ring = self.base_ring()

            if isinstance(point, HyperbolicPoint):
                point = point.coordinates(model="klein")

            point = (ring**2)(point)
            point.set_immutable()
            return point

        def _polygon_point(self, label, coordinates):
            r"""
            Return the point of the hyperbolic plane with the Klein
            ``coordinates`` in the polygon ``label``.

            This is the inverse of :meth:`_polygon_coordinates`.

            INPUT:

            - ``label`` -- a polygon label

            - ``coordinates`` -- an immutable vector as returned by
              :meth:`_polygon_coordinates`

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
                0
                sage: S._polygon_point(0, vector((0, 3/5)))
                2*I

            """
            x, y = coordinates
            return self.hyperbolic_plane().point(x, y, model="klein", check=False)

        def some_elements(self):
            r"""
            Return some typical points of this surface (for testing): a point
            in the interior of each of the first polygons, then the vertices.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S.glue((0, 0), (0, 3))
                sage: S.glue((0, 1), (0, 2))
                sage: S.set_immutable()
                sage: points = list(S.some_elements())
                sage: points[0]
                Point (1/4, 1/4) of polygon 0
                sage: len(points)
                4

            For a non-convex polygon, the point is taken near a convex
            corner::

                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)]))
                0
                sage: S.an_element()
                Point (-1/6, 0) of polygon 0

            """
            from itertools import islice

            from sage.modules.free_module_element import vector

            for label in islice(self.labels(), 32):
                polygon = self.polygon(label)
                vertices = [
                    vector(v.coordinates(model="klein")) for v in polygon.vertices()
                ]
                n = len(vertices)

                if polygon.is_convex():
                    # Convexity in the Klein model is Euclidean convexity, so
                    # the average of the vertices is inside.
                    point = sum(vertices) / n
                else:
                    # The lowest of the leftmost vertices is a strictly convex
                    # corner, so the polygon contains a neighbourhood of that
                    # corner in the cone spanned by its two sides.
                    corner = min(vertices, key=tuple)
                    i = vertices.index(corner)
                    direction = (vertices[i - 1] - corner) + (
                        vertices[(i + 1) % n] - corner
                    )
                    t = self.base_ring()(1) / 3
                    while True:
                        point = corner + t * direction
                        if polygon.get_point_position(point).is_in_interior():
                            break
                        t /= 2

                yield self(label, point)  # pylint: disable=not-callable

            if self.is_finite_type():
                yield from self.vertices()

        def _describe_surface(self):
            r"""
            Return a string describing this kind of surface.

            This is a helper method for ``_repr_``.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S._describe_surface()
                'Hyperbolic Surface with boundary'
                sage: S.glue((0, 0), (0, 3))
                sage: S.glue((0, 1), (0, 2))
                sage: S.set_immutable()
                sage: S._describe_surface()
                'Genus 0 Hyperbolic Surface'

            ::

                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
                0
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
                1
                sage: S._describe_surface()
                'Disconnected Surface'

            """
            if not self.is_finite_type():
                return "Surface built from infinitely many polygons"

            if not self.is_connected():
                return "Disconnected Surface"

            description = "Hyperbolic Surface"

            # The genus is only provided by the categories of connected
            # surfaces without boundary, i.e., not for mutable surfaces.
            if self.genus is not NotImplemented:
                description = f"Genus {self.genus()} {description}"

            if self.is_with_boundary():
                description += " with boundary"

            return description

        def _mismatched_gluing(self):
            r"""
            Return the first pair of glued edges of this surface that do not
            have the same length, or ``None`` if there is none.

            For a surface of infinite type, only the edges of the first 32
            polygons are considered.

            OUTPUT:

            A pair of pairs (polygon label, edge index) or ``None``.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S.glue((0, 0), (0, 3))
                sage: S._mismatched_gluing() is None
                True
                sage: S.glue((0, 1), (0, 0))
                sage: S._mismatched_gluing()
                ((0, 0), (0, 1))

            """
            from itertools import islice

            labels = self.labels()
            if not self.is_finite_type():
                labels = islice(labels, 32)

            for label in labels:
                sides = self.polygon(label)._cosh2_sides()
                for edge, cosh2 in enumerate(sides):
                    cross = self.opposite_edge(label, edge)
                    if cross is None:
                        continue

                    cross_label, cross_edge = cross
                    if self.polygon(cross_label)._cosh2_sides()[cross_edge] != cosh2:
                        return (label, edge), cross

            return None

        def _test_edge_lengths(self, **options):
            r"""
            Verify that the edges glued to each other have the same length.

            The test is skipped for mutable surfaces, which may go through
            states with mismatched gluings while they are built.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S.glue((0, 0), (0, 1))
                sage: S._test_edge_lengths()
                sage: S.glue((0, 0), (0, 3))
                sage: S.glue((0, 1), (0, 2))
                sage: S.set_immutable()
                sage: S._test_edge_lengths()

            """
            if self.is_mutable():
                return

            tester = self._tester(**options)

            tester.assertIsNone(self._mismatched_gluing())

    class ElementMethods:
        r"""
        Provides methods available to all points of hyperbolic cone surfaces.

        If you want to add functionality for such points, you most likely
        want to put it here.
        """

        def edges(self):
            r"""
            Return the edges of the polygons that contain this point.

            Implements
            :meth:`.polygonal_surfaces.PolygonalSurfaces.ElementMethods.edges`
            for hyperbolic cone surfaces.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S.glue((0, 0), (0, 3))
                sage: S.glue((0, 1), (0, 2))
                sage: S.set_immutable()

                sage: S(0, (1/4, 1/4)).edges()
                set()
                sage: sorted(S(0, (1/4, 0)).edges())
                [(0, 0), (0, 3)]
                sage: sorted(S(0, 1).edges())
                [(0, 0), (0, 1), (0, 2), (0, 3)]

            """
            surface = self.parent()

            edges = set()
            for label, coordinates in self.representatives():
                polygon = surface.polygon(label)
                position = polygon.get_point_position(coordinates)

                if position.is_in_edge_interior():
                    edges.add((label, position.get_edge()))
                elif position.is_vertex():
                    vertex = position.get_vertex()
                    n = len(polygon.vertices())
                    edges.add((label, vertex))
                    edges.add((label, (vertex - 1) % n))

            return edges

    class Oriented(SurfaceCategoryWithAxiom):
        r"""
        The category of hyperbolic cone surfaces whose orientation is
        compatible with the orientation of the hyperbolic plane that their
        polygons live in.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import HyperbolicConeSurfaces
            sage: HyperbolicConeSurfaces().Oriented()
            Category of oriented hyperbolic cone surfaces

        """

        class ParentMethods:
            r"""
            Provides methods available to all oriented hyperbolic cone
            surfaces.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def edge_transformation(self, label, edge):
                r"""
                Return the orientation preserving isometry that sends the
                ``edge`` of the polygon ``label`` to the edge it is glued to.

                The start of the edge is sent to the end of the opposite edge
                and its end to the start of the opposite edge.

                INPUT:

                - ``label`` -- a polygon label

                - ``edge`` -- an edge index of that polygon

                OUTPUT:

                A :class:`~flatsurf.geometry.hyperbolic_isometry.HyperbolicIsometry`
                over the base ring of this surface.

                EXAMPLES::

                    sage: from flatsurf import HyperbolicPlane
                    sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                    sage: H = HyperbolicPlane(QQ)
                    sage: S = MutableOrientedHyperbolicSurface(QQ)
                    sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    0
                    sage: S.glue((0, 0), (0, 3))
                    sage: S.glue((0, 1), (0, 2))
                    sage: S.set_immutable()

                    sage: g = S.edge_transformation(0, 1)
                    sage: g
                    [   0    1    0]
                    [-5/3    0  4/3]
                    [-4/3    0  5/3]
                    sage: g * vector((1/2, 0)), g * vector((1/2, 1/2))
                    ((0, 1/2), (1/2, 1/2))

                The isometry of a self-glued edge is the half turn about its
                midpoint, which may need a quadratic extension::

                    sage: S = MutableOrientedHyperbolicSurface(QQ)
                    sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    0
                    sage: S.glue((0, 0), (0, 0))
                    sage: S.edge_transformation(0, 0)
                    Traceback (most recent call last):
                    ...
                    ValueError: isometry is not defined over Rational Field; adjoin a square root of 4/3

                TESTS::

                    sage: S.edge_transformation(0, 1)
                    Traceback (most recent call last):
                    ...
                    ValueError: edge 1 of polygon 0 is not glued

                """
                from flatsurf.geometry.hyperbolic_isometry import (
                    HyperbolicIsometryGroup,
                    segment_isometry,
                )

                cross = self.opposite_edge(label, edge)
                if cross is None:
                    raise ValueError(f"edge {edge} of polygon {label} is not glued")

                cross_label, cross_edge = cross

                polygon = self.polygon(label)
                cross_polygon = self.polygon(cross_label)

                K = self.base_ring()

                matrix = segment_isometry(
                    (polygon.vertex(edge), polygon.vertex(edge + 1)),
                    (
                        cross_polygon.vertex(cross_edge + 1),
                        cross_polygon.vertex(cross_edge),
                    ),
                    ring=K,
                )

                return HyperbolicIsometryGroup(K)(matrix)

            def edge_matrix(self, label, edge):
                r"""
                Return the matrix of
                `\mathrm{SO}(2, 1)^+` of the isometry that sends the ``edge``
                of the polygon ``label`` to the edge it is glued to.

                See :meth:`edge_transformation` for the isometry itself.

                EXAMPLES::

                    sage: from flatsurf import HyperbolicPlane
                    sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                    sage: H = HyperbolicPlane(QQ)
                    sage: S = MutableOrientedHyperbolicSurface(QQ)
                    sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    0
                    sage: S.glue((0, 0), (0, 3))
                    sage: S.glue((0, 1), (0, 2))
                    sage: S.set_immutable()

                    sage: S.edge_matrix(0, 0)
                    [ 0 -1  0]
                    [ 1  0  0]
                    [ 0  0  1]

                """
                return self.edge_transformation(label, edge).matrix()

            def _angle_vertex(self, vertex, numerical):
                r"""
                Return the total angle at ``vertex`` as a multiple of 2π.

                This is a helper method for
                :meth:`~flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.Oriented.ElementMethods.angle`;
                it sums the interior angles of the polygons at the corners of
                ``vertex``.

                INPUT:

                - ``vertex`` -- a vertex of this surface

                - ``numerical`` -- a boolean; whether to return a floating
                  point approximation, or an exact rational number, which
                  requires every corner angle to be a rational multiple of π
                  (see
                  :meth:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon.angle`)

                EXAMPLES::

                    sage: from flatsurf import HyperbolicPlane
                    sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                    sage: H = HyperbolicPlane(QQ)
                    sage: S = MutableOrientedHyperbolicSurface(QQ)
                    sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    0
                    sage: S.glue((0, 0), (0, 3))
                    sage: S.glue((0, 1), (0, 2))
                    sage: S.set_immutable()

                    sage: S._angle_vertex(S(0, 1), numerical=False)
                    1/2
                    sage: S._angle_vertex(S(0, 2), numerical=True)  # abs tol 1e-12
                    0.195913276015304
                    sage: S._angle_vertex(S(0, 2), numerical=False)
                    Traceback (most recent call last):
                    ...
                    ValueError: angle is not a rational multiple of 2π; set numerical=True to get a numerical approximation

                """
                total = 0
                for label, point in vertex.representatives():
                    polygon = self.polygon(label)
                    corner = polygon.get_point_position(point).get_vertex()
                    total += polygon.angle(corner, numerical=numerical)

                return total
