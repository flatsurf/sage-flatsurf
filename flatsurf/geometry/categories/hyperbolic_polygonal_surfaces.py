r"""
The category of surfaces built by gluing hyperbolic polygons.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedHyperbolicSurface, HyperbolicPlane
    sage: H = HyperbolicPlane(QQ)
    sage: C = MutableOrientedHyperbolicSurface(H).category()

    sage: from flatsurf.geometry.categories import HyperbolicPolygonalSurfaces
    sage: C.is_subcategory(HyperbolicPolygonalSurfaces())
    True

"""
# ####################################################################
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import SurfaceCategory


class HyperbolicPolygonalSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing hyperbolic polygons.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import HyperbolicPolygonalSurfaces
        sage: HyperbolicPolygonalSurfaces()
        Category of hyperbolic polygonal surfaces

    """

    def super_categories(self):
        r"""
        Return the categories such surfaces are also automatically contained
        in, namely the category of surfaces built from polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import HyperbolicPolygonalSurfaces
            sage: C = HyperbolicPolygonalSurfaces()
            sage: C.super_categories()
            [Category of polygonal surfaces]

        """
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces

        return [PolygonalSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces that are built by gluing
        hyperboic polygons.

        If you want to provide methods to such surfaces, independent of the
        nature of these gluings, you should put them here.
        """

        def _an_element_(self):
            r"""
            Return a point in this surface, mostly used for automated testing.

            EXAMPLES::

                sage: from flatsurf import MutableOrientedHyperbolicSurface, HyperbolicPlane
                sage: H = HyperbolicPlane(QQ)
                sage: S = MutableOrientedHyperbolicSurface(H)

                sage: S.add_polygon(H.convex_hull(0, I + 2, I - 2))
                0

                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 1))

                sage: S.an_element()
                Vertex 0 of polygon 0

            """
            label = next(iter(self.labels()))
            polygon = self.polygon(label)

            return self(label, polygon.vertices()[0])
