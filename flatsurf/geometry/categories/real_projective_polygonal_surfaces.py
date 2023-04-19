r"""
The category of surfaces built by gluing Euclidean polygons.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict(QQ).category()

    sage: from flatsurf.geometry.categories import RealProjectivePolygonalSurfaces
    sage: C.is_subcategory(RealProjectivePolygonalSurfaces())
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


class RealProjectivePolygonalSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons (or more
    generally, polygons in two-dimensional real-projective space.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories import RealProjectivePolygonalSurfaces
        sage: RealProjectivePolygonalSurfaces()
        Category of real projective polygonal surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces
        return [PolygonalSurfaces()]
