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
    Hyperbolic Surface built from 2 quadrilaterals

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
    def super_categories(self):
        from flatsurf.geometry.categories.hyperbolic_polygonal_surfaces import HyperbolicPolygonalSurfaces

        return [HyperbolicPolygonalSurfaces()]

    class ParentMethods:
        def cusps(self):
            return set(vertex for vertex in self.vertices() if next(iter(vertex.representatives()))[1].is_ideal())

        def orbifold_points(self):
            orbifold_points = [vertex for vertex in self.vertices() if self.is_orbifold_point(vertex)]

            for (a, b) in self.gluings():
                if a == b:
                    orbifold_points.append(self(a[0], self.polygon(a[0]).edges()[a[1]].midpoint()))

            return orbifold_points

        def is_orbifold_point(self, point):
            label, p = point.representative()
            position = self.polygon(label).get_point_position(p)

            if p.is_ideal():
                return False

            if position.is_in_interior():
                return False

            if position.is_in_edge_interior():
                edge = position.edge()
                # This is an orbifold iff it is on a self-glued edge
                return self.opposite_edge(label, edge) == (label, edge)

            angle = point.angle(numerical=True)

            n = (1/angle).round()

            if n <= 1:
                return False

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
            start_edge_geometry = self.polygon(start_label).edges()[start_edge]
            end_label, end_edge = label, (start_edge - 1) % len(self.polygon(start_label).vertices())
            end_edge_geometry = self.polygon(end_label).edges()[end_edge]

            while True:
                label, edge = self.opposite_edge(end_label, end_edge)
                if label == start_label and edge == start_edge:
                    break

                normalization = end_edge_geometry.parent().isometry(-end_edge_geometry, self.polygon(label).edges()[edge])

                end_label, end_edge = label, (edge - 1) % len(self.polygon(label).vertices())

                end_edge_geometry = self.polygon(end_label).edges()[end_edge]
                end_edge_geometry = normalization * end_edge_geometry

            # Now we check if by arranging n copies of this widget at the vertex, we get a 2π angle.
            isometry = start_edge_geometry.parent().isometry(start_edge_geometry, -end_edge_geometry)
            identity = isometry ** n
            return identity == 1 or identity == -1

        def angle(self, point, numerical=False):
            if numerical:
                from sage.all import RR
                angle = RR(0)
            else:
                from sage.all import ZZ
                angle = ZZ(0)

            if not point.is_vertex():
                label, point = point.representative()
                position = self.polygon(label).get_point_position(point)
                if position.is_in_interior():
                    angle += 1
                elif self.opposite_edge(label, position.get_edge()) == (label, position.get_edge()):
                    # point on self-glued edge
                    angle += ZZ(1) / 2
                else:
                    # point on non-self-glued edge
                    angle += 1

                return angle

            # TODO: Check that all edges at this vertex are glued.

            for label, vertex in point.vertices():
                angle += self.polygon(label).angle(vertex, numerical=numerical)

            return angle

        def _describe_surface(self):
            if not self.is_finite_type():
                return "Surface built from infinitely many polygons"

            if not self.is_connected():
                return "Disconnected surface"

            description = "Hyperbolic Surface"

            if self.genus is not NotImplemented:
                description = f"Genus {self.genus()} {description}"

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
