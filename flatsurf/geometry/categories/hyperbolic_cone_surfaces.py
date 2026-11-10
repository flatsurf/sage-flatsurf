r"""
The category of hyperbolic cone surfaces.

A hyperbolic cone surface is a surface that can be built by gluing compact
polygons of the hyperbolic plane along their edges with isometries. The
metric is hyperbolic away from the vertices, which can be cone points of any
total angle.

Two geodesic segments of the hyperbolic plane are related by an isometry only
when they have the same length, so there is no hyperbolic counterpart of the
gluings of edges of different lengths by similarities: the edges glued to each
other have the same length. On an oriented surface, the gluing is the unique
orientation preserving isometry that sends one edge onto the other one with
reverse orientation.

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

from flatsurf.cache import cached_surface_method
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
            in the interior of each of the first polygons, a point in the
            interior of each self-glued or unglued edge of these polygons,
            then the vertices.

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

            The unglued edges contribute a point each::

                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S.glue((0, 0), (0, 3))
                sage: S.set_immutable()
                sage: list(S.some_elements())
                [Point (1/4, 1/4) of polygon 0,
                 Point (1/2, 1/4) of polygon 0,
                 Point (1/4, 1/2) of polygon 0,
                 Vertex 0 of polygon 0,
                 Vertex 1 of polygon 0,
                 Vertex 2 of polygon 0]

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

            labels = list(islice(self.labels(), 32))

            for label in labels:
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

            for label in labels:
                polygon = self.polygon(label)
                vertices = [
                    vector(v.coordinates(model="klein")) for v in polygon.vertices()
                ]
                n = len(vertices)

                for e in range(n):
                    opposite = self.opposite_edge(label, e)
                    if opposite is None or opposite == (label, e):
                        # Geodesics are Euclidean segments in the Klein model.
                        yield self(  # pylint: disable=not-callable
                            label, (vertices[e] + vertices[(e + 1) % n]) / 2
                        )

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

                sage: from flatsurf.geometry.surface import MutablePolygonalSurface
                sage: MutablePolygonalSurface.glue(S, (0, 1), (0, 0))
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
            states with mismatched gluings while they are built, except for
            a
            :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`:
            its ``glue()`` refuses such gluings, so the property holds in
            every state.

            EXAMPLES::

                sage: from flatsurf import HyperbolicPlane
                sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(QQ)
                sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                0
                sage: S.glue((0, 0), (0, 3))
                sage: S.glue((0, 1), (0, 2))
                sage: S._test_edge_lengths()

            The test fails on a surface glued with the shared method, which
            does not check the lengths::

                sage: from flatsurf.geometry.surface import MutablePolygonalSurface
                sage: MutablePolygonalSurface.glue(S, (0, 1), (0, 0))
                sage: S._test_edge_lengths()
                Traceback (most recent call last):
                ...
                AssertionError: ((0, 0), (0, 1)) is not None

            """
            from flatsurf.geometry.hyperbolic_surface import (
                MutableOrientedHyperbolicSurface,
            )

            if self.is_mutable() and not isinstance(
                self, MutableOrientedHyperbolicSurface
            ):
                return

            tester = self._tester(**options)

            tester.assertIsNone(self._mismatched_gluing())

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

                TESTS::

                    sage: S = MutableOrientedHyperbolicSurface(QQ)
                    sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    0
                    sage: S.glue((0, 0), (0, 3))
                    sage: S.edge_transformation(0, 1)
                    Traceback (most recent call last):
                    ...
                    ValueError: edge 1 of polygon 0 is not glued

                """
                from flatsurf.geometry.hyperbolic_isometry import (
                    HyperbolicIsometryGroup,
                )

                cross = self.opposite_edge(label, edge)
                if cross is None:
                    raise ValueError(f"edge {edge} of polygon {label} is not glued")

                matrix = self._side_pairing_matrix(label, edge, *cross)

                return HyperbolicIsometryGroup(self.base_ring())(matrix)

            def _side_pairing_matrix(self, label, edge, cross_label, cross_edge):
                r"""
                Return the matrix of the orientation preserving isometry that
                sends the ``edge`` of the polygon ``label`` onto the edge
                ``cross_edge`` of ``cross_label`` with reverse orientation,
                whether or not these edges are glued.

                INPUT:

                - ``label``, ``edge`` -- a polygon label and an edge index

                - ``cross_label``, ``cross_edge`` -- a polygon label and an
                  edge index

                OUTPUT:

                A matrix over the base ring of this surface. A ``ValueError``
                is raised if the edges do not have the same length or if the
                isometry is not defined over the base ring.

                EXAMPLES::

                    sage: from flatsurf import HyperbolicPlane
                    sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                    sage: H = HyperbolicPlane(QQ)
                    sage: S = MutableOrientedHyperbolicSurface(QQ)
                    sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    0
                    sage: S._side_pairing_matrix(0, 0, 0, 3)
                    [ 0 -1  0]
                    [ 1  0  0]
                    [ 0  0  1]
                    sage: S._side_pairing_matrix(0, 0, 0, 1)
                    Traceback (most recent call last):
                    ...
                    ValueError: segments do not have the same length

                """
                from flatsurf.geometry.hyperbolic_isometry import segment_isometry

                polygon = self.polygon(label)
                cross_polygon = self.polygon(cross_label)

                return segment_isometry(
                    (polygon.vertex(edge), polygon.vertex(edge + 1)),
                    (
                        cross_polygon.vertex(cross_edge + 1),
                        cross_polygon.vertex(cross_edge),
                    ),
                    ring=self.base_ring(),
                )

            def edge_matrix(self, label, edge):
                r"""
                Return the matrix of
                `\mathrm{SO}(2, 1)^+` of the isometry that sends the ``edge``
                of the polygon ``label`` to the edge it is glued to.

                See :meth:`edge_transformation` for the isometry itself.

                INPUT:

                - ``label`` -- a polygon label

                - ``edge`` -- an edge index of that polygon

                OUTPUT:

                A `3 \times 3` matrix over the base ring of this surface.

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

            @cached_surface_method
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

                When the corners have irrational angles, the exact total of
                an interior vertex is read from the trace of the holonomy of
                its loop::

                    sage: T = MutableOrientedHyperbolicSurface(QQ)
                    sage: T.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
                    0
                    sage: T.add_polygon(H.polygon(vertices=[(1/2, 0), (1/2, 1/2), (0, 1/2)]))
                    1
                    sage: T.glue((0, 1), (1, 2)); T.glue((0, 0), (0, 2)); T.glue((1, 0), (1, 1))
                    sage: T.set_immutable()
                    sage: T(1, 0).angle()
                    1/2
                    sage: T(1, 1).angle()
                    Traceback (most recent call last):
                    ...
                    ValueError: angle is not a rational multiple of 2π; set numerical=True to get a numerical approximation

                """
                total = 0
                for label, point in vertex.representatives():
                    polygon = self.polygon(label)
                    corner = polygon.get_point_position(point).get_vertex()
                    try:
                        total += polygon.angle(corner, numerical=numerical)
                    except ValueError:
                        if numerical:
                            raise
                        from flatsurf.geometry.holonomy import _holonomy_angle

                        return _holonomy_angle(self, vertex)

                return total
