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
    Surface built from 2 quadrilaterals

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
