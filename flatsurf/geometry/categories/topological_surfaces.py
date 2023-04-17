r"""
The category of topological surfaces.

This provides a base category for all the surfaces in sage-flatsurf.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict().category()

    sage: from flatsurf.geometry.categories.topological_surfaces import TopologicalSurfaces
    sage: C.is_subcategory(TopologicalSurfaces())
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

from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom, all_axioms
from sage.categories.topological_spaces import TopologicalSpaces


class TopologicalSurfaces(Category):
    r"""
    The category of topological surfaces.

    This category does not provide much functionality but just a common base
    for all the other categories defined in sage-flatsurf.

    In particular, this does not really require a topology since there is no
    general concept of open subsets of a surface in sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.topological_surfaces import TopologicalSurfaces
        sage: TopologicalSurfaces()

    """

    def super_categories(self):
        return [TopologicalSpaces()]

    class Oriented(CategoryWithAxiom):
        r"""
        The axiom satisfied by surfaces that can be oriented.

        As of 2023, all surfaces in sage-flatsurf satisfy this axiom.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: 'Oriented' in S.category().axioms()
            True

        """


all_axioms += ("Oriented",)
