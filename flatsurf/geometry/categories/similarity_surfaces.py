r"""
The category of similarity surfaces.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from Euclidean polygons that are glued by similarities, i.e., identified
edges can be transformed into each other by application of rotation and
homothety (scaling) and translation.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict().category()

    sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
    sage: C.is_subcategory(SimilaritySurfaces())
    True

The easiest way to construct a similarity surface is to use the pre-built
constructions from
:class:`flatsurf.geometry.similarity_surface_generators.SimilaritySurfaceGenerators`::

    sage: from flatsurf import polygons, similarity_surfaces
    sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
    sage: similarity_surfaces.self_glued_polygon(P)
    HalfTranslationSurface built from 1 polygon

The second way is to build a surface (using e.g. :class:`flatsurf.geometry.surface.Surface_list`)
and then use this surface as an argument for class:`SimilaritySurface`)::

    sage: from flatsurf.geometry.similarity_surface import SimilaritySurface
    sage: from flatsurf.geometry.surface import Surface_list
    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: Stop = Surface_list(QQ)
    sage: Stop.add_polygon(P)
    0
    sage: Stop.add_polygon(2*P)
    1
    sage: Stop.add_polygon(3*P)
    2
    sage: Stop.set_edge_pairing(0, 1, 1, 3)
    sage: Stop.set_edge_pairing(0, 0, 2, 2)
    sage: Stop.set_edge_pairing(0, 2, 2, 0)
    sage: Stop.set_edge_pairing(0, 3, 1, 1)
    sage: Stop.set_edge_pairing(1, 2, 2, 1)
    sage: Stop.set_edge_pairing(1, 0, 2, 3)
    sage: S = SimilaritySurface(Stop)
    sage: S
    SimilaritySurface built from 3 polygons

To perform a sanity check on the obtained surface, you can run its test
suite::

    sage: TestSuite(S).run()

In the following example, we build two broken surfaces and
check that the test suite fails as expected::

    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: Stop = Surface_list(QQ)
    sage: Stop.add_polygon(P)
    0
    sage: S = SimilaritySurface(Stop)
    sage: TestSuite(S).run()
    ...
      AssertionError: edge (0, 0) is not glued
      ------------------------------------------------------------
      The following tests failed: _test_gluings
    Failure in _test_underlying_surface
    The following tests failed: _test_underlying_surface

    sage: Stop.set_edge_pairing(0, 0, 0, 3)
    sage: Stop.set_edge_pairing(0, 1, 0, 3)
    sage: Stop.set_edge_pairing(0, 2, 0, 3)
    sage: S = SimilaritySurface(Stop)
    sage: TestSuite(S).run()
    ...
      AssertionError: edge gluing is not a pairing:
      (0, 0) -> (0, 3) -> (0, 2)
      ------------------------------------------------------------
      The following tests failed: _test_gluings
    Failure in _test_underlying_surface
    The following tests failed: _test_underlying_surface

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

from sage.categories.category import Category


class SimilaritySurfaces(Category):
    r"""
    The category of surfaces built from polygons with edges identified by
    similarities.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
        sage: SimilaritySurfaces()

    """

    def super_categories(self):
        from flatsurf.geometry.categories.real_projective_polygonal_surfaces import RealProjectivePolygonalSurfaces
        return [RealProjectivePolygonalSurfaces()]
