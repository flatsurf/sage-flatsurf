r"""
The category of dilation surfaces.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from Euclidean polygons that are glued by translation followed by
homothety, i.e., application of a diagonal matrix.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import polygons, similarity_surfaces
    sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
    sage: S = similarity_surfaces.self_glued_polygon(P)
    sage: C = S.category()

    sage: from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
    sage: C.is_subcategory(DilationSurfaces())
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
from sage.categories.category_with_axiom import CategoryWithAxiom, all_axioms


class DilationSurfaces(Category):
    r"""
    The category of surfaces built from polygons with edges identified by
    translations and homothety.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
        sage: DilationSurfaces()
        Category of dilation surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
        return [SimilaritySurfaces()]

    class Positive(CategoryWithAxiom):
        r"""
        The axiom satisfied by dilation surfaces that use homothety with
        positive scaling.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: 'Positive' in S.category().axioms()
            True

        """

    class SubcategoryMethods:
        def Positive(self):
            r"""
            Return the subcategory of surfaces glued by positive dilation.

            EXAMPLES::

                sage: from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
                sage: C = DilationSurfaces()
                sage: C.Positive()
                Category of positive dilation surfaces

            """
            return self._with_axiom("Positive")


all_axioms += ('Positive',)
