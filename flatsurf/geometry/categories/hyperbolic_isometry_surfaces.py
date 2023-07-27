r"""
The category of hyperbolic surfaces built from polygons glued with isometries.

This module provides shared functionality for all surfaces in sage-flatsurf
that are built from hyperbolic polygons that are glued by isometries, i.e.,
identified edges can be transformed into each other by application of an
isometry of the hyperbolic plane.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedHyperbolicSurface, HyperbolicPlane
    sage: H = HyperbolicPlane(QQ)
    sage: C = MutableOrientedHyperbolicSurface(H).category()

    sage: from flatsurf.geometry.categories import HyperbolicIsometrySurfaces
    sage: C.is_subcategory(HyperbolicIsometrySurfaces())
    True

We create a surface from scratch (using
:class:`~flatsurf.geometry.surface.MutableOrientedHyperbolicSurface`)::

    sage: S = MutableOrientedHyperbolicSurface(H)
    sage: S.add_polygon(H.polygon([
    ....:     H.vertical(3).left_half_space(),
    ....:     H.vertical(1).right_half_space(),
    ....:     H.half_circle(2, 4).left_half_space(),
    ....:     H.half_circle(2, 16).right_half_space(),
    ....: ]))
    0
    sage: S.add_polygon(H.polygon([
    ....:     H.vertical(3).left_half_space(),
    ....:     H.vertical(1).right_half_space(),
    ....:     H.half_circle(2, 4).left_half_space(),
    ....:     H.half_circle(2, 16).right_half_space(),
    ....: ]))
    1

    sage: S.polygon(0).edges()
    {{(x^2 + y^2) - 4*x = 0} ∩ {3*(x^2 + y^2) + 5*x - 17 ≥ 0} ∩ {(x^2 + y^2) + 13*x - 51 ≤ 0},
     {-x + 3 = 0} ∩ {10*(x^2 + y^2) - 39*x - 3 ≥ 0} ∩ {2*(x^2 + y^2) - 15*x - 3 ≤ 0},
     {-(x^2 + y^2) + 4*x + 12 = 0} ∩ {7*(x^2 + y^2) - 65*x + 27 ≥ 0} ∩ {9*(x^2 + y^2) - 221*x + 77 ≤ 0},
     {x - 1 = 0} ∩ {2*(x^2 + y^2) - 17*x - 15 ≤ 0} ∩ {2*(x^2 + y^2) - 5*x - 3 ≥ 0}}

    sage: S.polygon(1).edges()
    {{(x^2 + y^2) - 4*x = 0} ∩ {3*(x^2 + y^2) + 5*x - 17 ≥ 0} ∩ {(x^2 + y^2) + 13*x - 51 ≤ 0},
     {-x + 3 = 0} ∩ {10*(x^2 + y^2) - 39*x - 3 ≥ 0} ∩ {2*(x^2 + y^2) - 15*x - 3 ≤ 0},
     {-(x^2 + y^2) + 4*x + 12 = 0} ∩ {7*(x^2 + y^2) - 65*x + 27 ≥ 0} ∩ {9*(x^2 + y^2) - 221*x + 77 ≤ 0},
     {x - 1 = 0} ∩ {2*(x^2 + y^2) - 17*x - 15 ≤ 0} ∩ {2*(x^2 + y^2) - 5*x - 3 ≥ 0}}

    sage: S.glue((0, 0), (1, 0))
    sage: S.glue((0, 1), (1, 3))
    sage: S.glue((0, 2), (1, 2))
    sage: S.glue((0, 3), (1, 1))

    sage: S
    Hyperbolic Surface with 2 orbifold points built from 2 quadrilaterals

To perform a sanity check on the obtained surface, you can run its test
suite::

    sage: TestSuite(S).run()

If there are no errors reported, no consistency problems could be detected in
your surface.

Once you mark the surface as immutable, it gets more functionality, e.g.,
coming from its structure as a surface without boundary. This also adds more
tests to its test suite::

    sage: S.category()
    Category of finite type oriented hyperbolic isometry surfaces
    sage: S.set_immutable()
    sage: S.category()
    Category of connected without boundary finite type oriented hyperbolic isometry surfaces

    sage: S.genus()
    0

    sage: TestSuite(S).run()

"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
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
from flatsurf.geometry.categories.surface_category import SurfaceCategory


class HyperbolicIsometrySurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing hyperbolic polygons with
    isometries.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import HyperbolicIsometrySurfaces
        sage: HyperbolicIsometrySurfaces()
        Category of hyperbolic isometry surfaces

    """

    def super_categories(self):
        r"""
        Return the categories such surfaces are also automatically contained
        in, namely the category of surfaces built by gluing hyperbolic
        polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import HyperbolicIsometrySurfaces
            sage: C = HyperbolicIsometrySurfaces()
            sage: C.super_categories()
            [Category of hyperbolic polygonal surfaces]

        """
        from flatsurf.geometry.categories.hyperbolic_polygonal_surfaces import HyperbolicPolygonalSurfaces

        return [HyperbolicPolygonalSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces built by gluing hyperbolic
        polygons with isometries.

        If you want to add functionality to all such surfaces, you most likely
        want to put it here.
        """

        def cusps(self):
            r"""
            Return the set of cusps of this surface, i.e., the set of vertices
            that are at infinite points.

            EXAMPLES::

                sage: from flatsurf import MutableOrientedHyperbolicSurface, HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(H)

                sage: S.add_polygon(H.convex_hull(0, I + 1, I - 1))
                0

                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 1))

                sage: S
                Hyperbolic Surface with 1 cusp and with 1 orbifold point built from a degenerate triangle

                sage: S.cusps()
                {Vertex 0 of polygon 0}

            """
            return set(vertex for vertex in self.vertices() if next(iter(vertex.representatives()))[1].is_ideal())

        def orbifold_points(self):
            r"""
            Return the set of orbifold points of this surface with their
            respective order, i.e., for each point an n such that the total
            angle at the point is 2π/n.

            EXAMPLES::

                sage: from flatsurf import MutableOrientedHyperbolicSurface, HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(H)

                sage: S.add_polygon(H.convex_hull(0, I + 2, I - 2))
                0

                sage: S.glue((0, 1), (0, 1))
                sage: S.glue((0, 0), (0, 2))

                sage: S
                Hyperbolic Surface with 1 cusp and with 1 orbifold point built from a degenerate triangle

                sage: S.orbifold_points()
                {Point (2 + I → -2 + I) / 2 of polygon 0}

            .. SEEALSO::

                :meth:`ElementMethods.orbifold_order` to determine the order of the orbifold points

            """
            vertices = {vertex for vertex in self.vertices() if vertex.orbifold_order()}
            midpoints = {self(a[0], self.polygon(a[0]).edges()[a[1]].midpoint()) for (a, b) in self.gluings() if a == b}

            return vertices.union(midpoints)

        def _describe_surface(self):
            if not self.is_finite_type():
                return "Surface built from infinitely many polygons"

            if not self.is_connected():
                return "Disconnected surface"

            description = "Hyperbolic Surface"

            if self.genus is not NotImplemented:
                description = f"Genus {self.genus()} {description}"

            if self.is_with_boundary():
                description = f"{description} with boundary"
            else:
                cusps = self.cusps()
                orbifold_points = self.orbifold_points()

                if len(cusps) == 0:
                    cusps = ""
                elif len(cusps) == 1:
                    cusps = "with 1 cusp"
                else:
                    cusps = f"with {len(cusps)} cusps"

                if len(orbifold_points) == 0:
                    orbifold_points = ""
                elif len(orbifold_points) == 1:
                    orbifold_points = "with 1 orbifold point"
                else:
                    orbifold_points = f"with {len(orbifold_points)} orbifold points"

                description = " ".join(filter(None, [description, " and ".join(filter(None, [cusps, orbifold_points]))]))

            return description

        def edge_transformation(self, label, edge):
            opposite_label, opposite_edge = self.opposite_edge(label, edge)
            return self._hyperbolic_plane.isometry(self.polygon(label).edges()[edge], -self.polygon(opposite_label).edges()[opposite_edge])

    class ElementMethods:
        def orbifold_order(self):
            r"""
            Return the order of this point or ``None`` if this is not an
            orbifold point.

            INPUT:

            - ``point`` -- a point of this surface

            EXAMPLES::

                sage: from flatsurf import MutableOrientedHyperbolicSurface, HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(H)

                sage: S.add_polygon(H.convex_hull(0, I + 2, I - 2))
                0

                sage: S.glue((0, 1), (0, 1))
                sage: S.glue((0, 0), (0, 2))

                sage: P = S(0, S.polygon(0).edges()[1].midpoint())
                sage: P.orbifold_order()
                2

                sage: P = S(0, 0)
                sage: P.orbifold_order() is None
                True

            """
            surface = self.parent()

            label, p = self.representative()
            position = surface.polygon(label).get_point_position(p)

            if p.is_ideal():
                return None

            if position.is_in_interior():
                return None

            if position.is_in_edge_interior():
                edge = position.get_edge()
                # This is an orbifold iff it is on a self-glued edge
                if surface.opposite_edge(label, edge) == (label, edge):
                    return 2
                return None

            angle = self.angle(numerical=True)

            order = (1/angle).round()

            if order <= 1:
                return None

            # This does not work because we cannot decide whether an angle is rational. We could probably do something like https://mathproblems123.wordpress.com/2018/03/31/when-is-arccos-a-rational-multiple-of-pi/
            # angle = point.angle()
            # from sage.all import QQ
            # if angle not in QQ:
            #     return False
            # angle = QQ(angle)
            # return angle.numerator() == 1 and angle.denominator() > 1

            # We rearrange all the polygons attached to this point so that they
            # are glued by identities. Then we take the first and the last edge
            # of this gadget to get a hold of the total angle at the point.
            start_label, start_edge = label, position.get_edge()
            start_edge_geometry = surface.polygon(start_label).edges()[start_edge]
            end_label, end_edge = label, (start_edge - 1) % len(surface.polygon(start_label).vertices())
            end_edge_geometry = surface.polygon(end_label).edges()[end_edge]

            while True:
                label, edge = surface.opposite_edge(end_label, end_edge)
                if label == start_label and edge == start_edge:
                    break

                normalization = end_edge_geometry.parent().isometry(-end_edge_geometry, surface.polygon(label).edges()[edge])

                end_label, end_edge = label, (edge - 1) % len(surface.polygon(label).vertices())

                end_edge_geometry = surface.polygon(end_label).edges()[end_edge]
                end_edge_geometry = normalization * end_edge_geometry

            # Now we check if by arranging n copies of this widget at the vertex, we get a 2π angle.
            isometry = start_edge_geometry.parent().isometry(start_edge_geometry, -end_edge_geometry)
            identity = isometry ** order
            if identity == 1 or identity == -1:
                return order

            return None

        def angle(self, numerical=False):
            surface = self.parent()
            # TODO: Make ElementMethods (and make sure that the Euclidean version is as well. Adapt news!
            if numerical:
                from sage.all import RR
                angle = RR(0)
            else:
                from sage.all import ZZ
                angle = ZZ(0)

            if not self.is_vertex():
                label, point = self.representative()
                position = surface.polygon(label).get_point_position(point)
                if position.is_in_interior():
                    angle += 1
                elif surface.opposite_edge(label, position.get_edge()) == (label, position.get_edge()):
                    # point on self-glued edge
                    angle += ZZ(1) / 2
                else:
                    # point on non-self-glued edge
                    angle += 1

                return angle

            # TODO: Check that all edges at this vertex are glued.

            for label, vertex in self.vertices():
                angle += surface.polygon(label).angle(vertex, numerical=numerical)

            return angle
