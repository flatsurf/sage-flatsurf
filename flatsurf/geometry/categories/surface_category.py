r"""
Category types for surfaces.

This module provides alternative implementations for the SageMath types
``Category`` and ``CategoryWithAxiom``. It patches the ``_cmp_key`` of these
types to produce a more stable sorting of surface categories exactly in the
same way that SageMath does for its builtin categories.

Without this, the MRO of surfaces is session dependent and in particular we get
somewhat random printing of categories such as translation surfaces.

While we do not claim to actually understand all the details here, you can
consult ``c3_controlled.py`` in SageMath for all the (rather technical)
details.
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

import sage.categories.category
import sage.categories.category_with_axiom
from sage.misc.c3_controlled import _cmp_key

flags = {
    atom: 1 << (30 - i)
    for i, atom in enumerate(
        [
            "TopologicalSurfaces",
            "PolygonalSurfaces",
            "EuclideanPolygonalSurfaces",
            "SimilaritySurfaces",
            "ConeSurfaces",
            "DilationSurfaces",
            "HalfTranslationSurfaces",
            "TranslationSurfaces",
        ]
    )
}


class SurfaceCmpKey:
    def __get__(self, instance, cls):
        (flag, counter) = instance._cmp_key_vanilla
        flag |= flags.get(cls.__base__.__name__, 0)
        for cat in instance._super_categories:
            flag |= cat._cmp_key[0]
        instance._cmp_key = (flag, counter)
        return flag, counter


class SurfaceCategory(sage.categories.category.Category):
    _cmp_key_vanilla = _cmp_key
    _cmp_key = SurfaceCmpKey()


class SurfaceCategoryWithAxiom(sage.categories.category_with_axiom.CategoryWithAxiom):
    _cmp_key_vanilla = _cmp_key
    _cmp_key = SurfaceCmpKey()
