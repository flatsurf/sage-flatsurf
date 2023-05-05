r"""
Categories of Surfaces.

sage-flatsurf uses SageMath categories to distinguish different kinds of
surfaces such as hyperbolic surfaces, translation surfaces, …. See
https://doc.sagemath.org/html/en/reference/categories/sage/categories/primer.html
for a detailed introduction of categories in SageMath. In short, "categories" are
not so much mathematical categories but more similar to a normal class
hierarchy in Python; however, they extend the idea of a classical class
hierarchy by allowing us to dynamically change the category (and methods) of a
surface as we learn more about it.

Note that you normally don't have to create categories explicitly. Categories
are deduced automatically (at least for surfaces of finite type.) You should
think of categories as an implementation detail. As a user of sage-flatsurf,
you don't need to know about them. As a developer of sage-flatsurf, they
provide entry points to place your code; e.g., to add a method to all
translation surfaces, actually add a method to
:class:`TranslationSurfaces.ParentMethods`.

.. NOTE::

    Categories are deduced by calling methods such as :meth:`is_orientable`,
    :meth:`is_with_boundary`, :meth:`is_compact`, :meth:`is_connected`,
    :meth:`is_finite_type`, :meth:`is_cone_surface`,
    :meth:`is_dilation_surface`, :meth:`is_translation_surface`, and
    :meth:`is_rational`. There are default implementations for these for finite
    type surfaces. Once a surfaces has been found to be in a subcategory, these
    methods are replaced to simply return ``True`` instead of computing
    anything. If a class explicitly overrides these methods, then that category
    machinery cannot replace that method anymore when the category of the
    surface gets refined. Consequently, it can be beneficial for the override
    to shortcut the question by querying the category, e.g., ``is_rational``
    could start with ``if "Rational" in self.category().axioms(): return True``
    before actually performing any computation.

EXAMPLES::

A single square without any gluings::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: S = MutableOrientedSimilaritySurface(QQ)

    sage: from flatsurf import polygons
    sage: S.add_polygon(polygons.square(), label=0)
    0

This is considered to be a surface built from polygons with all gluings being
similarities (however there are none)::

    sage: S.category()
    Category of finite type oriented similarity surfaces

It does not really make sense to ask which stratum this surface belongs to::

    sage: S.stratum()
    Traceback (most recent call last):
    ...
    AttributeError: ... has no attribute 'stratum'

Once we add gluings, this turns into a square torus::

    sage: S.glue((0, 0), (0, 2))
    sage: S.glue((0, 1), (0, 3))

We signal to sage-flatsurf that we are done building this surface, and its
category gets refined::

    sage: S.set_immutable()
    sage: S.category()
    Category of connected without boundary finite type translation surfaces

Since this is now a translation surface, we can ask for its stratum again::

    sage: S.stratum()
    H_1(0)

"""
# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2021-2023 Julian Rüth
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

from flatsurf.geometry.categories.topological_surfaces import TopologicalSurfaces
from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces
from flatsurf.geometry.categories.real_projective_polygonal_surfaces import RealProjectivePolygonalSurfaces
from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces
from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
