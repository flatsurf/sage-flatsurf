r"""
The category of surfaces built by gluing Euclidean polygons.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: C = MutableOrientedSimilaritySurface(QQ).category()

    sage: from flatsurf.geometry.categories import EuclideanPolygonalSurfaces
    sage: C.is_subcategory(EuclideanPolygonalSurfaces())
    True

"""
# ####################################################################
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import SurfaceCategory


class EuclideanPolygonalSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing Euclidean polygons (or more
    generally, polygons in two-dimensional real space.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories import EuclideanPolygonalSurfaces
        sage: EuclideanPolygonalSurfaces()
        Category of euclidean polygonal surfaces

    """

    def super_categories(self):
        r"""
        The categories such surfaces are also automatically contained in,
        namely the category of surfaces built from polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygonalSurfaces
            sage: C = EuclideanPolygonalSurfaces()
            sage: C.super_categories()
            [Category of polygonal surfaces]

        """
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces

        return [PolygonalSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces that are built from polygons
        in the real plane.

        If you want to add functionality for such surfaces you most likely
        want to put it here.
        """

        def _an_element_(self):
            r"""
            Return a point on this surface.

            EXAMPLES::

                sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                sage: s = SimilaritySurfaceGenerators.example()
                sage: s.an_element()
                Point (4/3, -2/3) of polygon 0

            ::

                sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface

                sage: S = MutableOrientedSimilaritySurface(QQ)
                sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
                0
                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 3))

                sage: S.an_element()
                Point (1/2, 1/2) of polygon 0

            TESTS:

            Verify that this method works over non-fields (if 2 is
            invertible)::

              sage: from flatsurf import similarity_surfaces
              sage: from flatsurf import EuclideanPolygonsWithAngles
              sage: E = EuclideanPolygonsWithAngles((3, 3, 5))
              sage: from pyexactreal import ExactReals # optional: exactreal  # random output due to pkg_resources deprecation warnings in some contexts
              sage: R = ExactReals(E.base_ring()) # optional: exactreal
              sage: angles = (3, 3, 5)
              sage: slopes = EuclideanPolygonsWithAngles(*angles).slopes()
              sage: P = Polygon(angles=angles, edges=[R.random_element() * slopes[0]])  # optional: exactreal
              sage: S = similarity_surfaces.billiard(P) # optional: exactreal
              sage: S.an_element()  # optional: exactreal
              Point ((1/2 ~ 0.50000000)*ℝ(0.303644…), 0) of polygon 0

            """
            label = next(iter(self.labels()))
            polygon = self.polygon(label)

            from sage.categories.all import Fields

            # We use a point that can be constructed without problems on an
            # infinite surface.
            if polygon.is_convex() and self.base_ring() in Fields():
                coordinates = polygon.centroid()
            else:
                # Sometimes, this is not implemented because it requires the edge
                # transformation to be known, so we prefer the centroid.
                coordinates = polygon.edge(0) / 2
                coordinates.set_immutable()
            return self(label, coordinates)  # pylint: disable=not-callable
