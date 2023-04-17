r"""
The category of translation surfaces.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from Euclidean polygons whose glued edges can be transformed into each
other with translations.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.square_torus()
    sage: C = S.category()

    sage: from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
    sage: C.is_subcategory(TranslationSurfaces())
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
from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces


class TranslationSurfaces(Category):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
        sage: TranslationSurfaces()
        Category of translation surfaces

    """
    _base_category_class_and_axiom = (HalfTranslationSurfaces, 'PositiveDilation')

    def super_categories(self):
        # TODO: They satisfy the axiom of no self-gluing.

        return [HalfTranslationSurfaces()]

    class Orientable(CategoryWithAxiom):
        class ParentMethods:
            def stratum(self):
                r"""
                Return the stratum this surface belongs to.

                This uses the package ``surface-dynamics``
                (see http://www.labri.fr/perso/vdelecro/flatsurf_sage.html)

                EXAMPLES::

                    sage: import flatsurf.geometry.similarity_surface_generators as sfg
                    sage: sfg.translation_surfaces.octagon_and_squares().stratum()
                    H_3(4)
                """
                from surface_dynamics import AbelianStratum
                from sage.rings.integer_ring import ZZ

                return AbelianStratum([ZZ(a - 1) for a in self.angles()])
