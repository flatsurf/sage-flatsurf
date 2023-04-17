r"""
The category of half-translation surfaces.

A half-translation surface is a surface built by gluing Euclidean polygons. The sides of the polygons can be glued with translations or half-translations (translation followed by a rotation of angle π.)

EXAMPLES:

We glue all the sides of a square to themselves. Since each gluing is just a
rotation of π, this is a half-translation surface::

    sage: from flatsurf import polygons, similarity_surfaces
    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = similarity_surfaces.self_glued_polygon(P)
    sage: S.set_immutable()

    sage: C = S.category()

    sage: from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
    sage: C.is_subcategory(HalfTranslationSurfaces())
    True

"""
# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2013-2019 Vincent Delecroix
#                      2013-2019 W. Patrick Hooper
#                           2023 Julian Rüth
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
from sage.categories.category_with_axiom import CategoryWithAxiom


class HalfTranslationSurfaces(Category):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations and half-translations (translations followed by rotations
    among an angle π.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
        sage: HalfTranslationSurfaces()

    """

    def super_categories(self):
        from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
        # TODO: We can be more limited here, something like HalfDilationSurfaces() & RationalConeSurfaces()
        return [SimilaritySurfaces()]

    class Oriented(CategoryWithAxiom):
        class ParentMethods:
            def stratum(self):
                r"""
                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                    sage: H = B.minimal_cover(cover_type="half-translation")
                    sage: H.stratum()
                    Q_1(3, -1^3)
                """
                angles = self.angles()
                if all(x.denominator() == 1 for x in angles):
                    raise NotImplementedError
                from surface_dynamics import QuadraticStratum

                return QuadraticStratum(*[2 * a - 2 for a in angles])
