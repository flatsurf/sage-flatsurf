r"""
Categories of Surfaces and Polygons.

sage-flatsurf uses SageMath categories to distinguish different kinds of
surfaces such as hyperbolic surfaces, translation surfaces, …. See
https://doc.sagemath.org/html/en/reference/categories/sage/categories/primer.html
for a detailed introduction of categories in SageMath. In short, "categories" are
not so much mathematical categories but more similar to a normal class
hierarchy in Python; however, they extend the idea of a classical class
hierarchy by allowing us to dynamically change the category (and methods) of a
surface as we learn more about it.

Note that you normally don't have to create categories explicitly. Categories
are deduced automatically (at least for surfaces of finite type). You should
think of categories as an implementation detail. As a user of sage-flatsurf,
you don't need to know about them. As a developer of sage-flatsurf, they
provide entry points to place your code; e.g., to add a method to all
translation surfaces, actually add a method to
:class:`translation_surfaces.TranslationSurfaces.ParentMethods`.

A similar but smaller hierarchy of categories exists for polygons, Euclidean
polygons, Hyperbolic polygons.

.. NOTE::

    Categories are deduced by calling methods such as
    :meth:`~topological_surfaces.TopologicalSurfaces.ParentMethods.is_orientable`,
    :meth:`~topological_surfaces.TopologicalSurfaces.ParentMethods.is_with_boundary`,
    :meth:`~topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`,
    :meth:`~topological_surfaces.TopologicalSurfaces.ParentMethods.is_connected`,
    :meth:`~polygonal_surfaces.PolygonalSurfaces.ParentMethods.is_finite_type`,
    :meth:`~similarity_surfaces.SimilaritySurfaces.ParentMethods.is_cone_surface`,
    :meth:`~similarity_surfaces.SimilaritySurfaces.ParentMethods.is_dilation_surface`,
    :meth:`~similarity_surfaces.SimilaritySurfaces.ParentMethods.is_translation_surface`,
    and
    :meth:`~similarity_surfaces.SimilaritySurfaces.ParentMethods.is_rational_surface`.
    There are default implementations for these for finite
    type surfaces. Once a surfaces has been found to be in a certain
    subcategory, these methods are replaced to simply return ``True`` instead
    of computing anything. If a class explicitly overrides these methods, then
    the category machinery cannot replace that method anymore when the category
    of the surface gets refined. Consequently, it can be beneficial for the
    override to shortcut the question by querying the category, e.g.,
    ``is_rational`` could start with ``if "Rational" in
    self.category().axioms(): return True`` before actually performing any
    computation.

EXAMPLES:

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
    AttributeError: ... has no attribute 'stratum'...

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

There are a number of workarounds if you want to compute things such as
``.stratum()`` on a mutable surface. For the sake of this demonstration, lets
make our surface mutable again::

    sage: S = MutableOrientedSimilaritySurface.from_surface(S)

The recommended approach is to create a copy of the surface, make it immutable
and then call the methods you need::

    sage: T = MutableOrientedSimilaritySurface.from_surface(S)
    sage: T.set_immutable()
    sage: T.stratum()
    H_1(0)

You might be worried about the performance implications but most of the time
that might be a `premature optimization
<https://en.wikipedia.org/wiki/Program_optimization#When_to_optimize>`_.

Often enough, the copy is actually not the problem but the time that is spent
to figure out that that copy is actually a translation surface. If you already
know that to be true, you can simplify things to::

    sage: from flatsurf.geometry.categories import TranslationSurfaces
    sage: T = MutableOrientedSimilaritySurface.from_surface(S, category=TranslationSurfaces().WithoutBoundary())
    sage: T.stratum()
    H_1(0)

You can also change the category of a mutable surface to provide all the
functionality that is available to surfaces in that category::

    sage: S._refine_category_(TranslationSurfaces().WithoutBoundary())
    sage: S.stratum()
    H_1(0)

Note however, that the category of a surface cannot be generalized anymore.
This surface is now at least a "translation surface", no matter what mutations
you make to it. (And the system will not check that it is indeed a translation
surface.) That approach might still be beneficial if you make lots of minor
changes to the surface, e.g., lots of edge flips, and want to query such
methods frequently.

Finally, we can try to call the
:meth:`~flatsurf.geometry.categories.translation_surfaces.TranslationSurfaces.FiniteType.WithoutBoundary.ParentMethods.stratum`
method directly but it might have dependencies on other methods that are not
available::

    sage: S = MutableOrientedSimilaritySurface.from_surface(S)

    sage: TranslationSurfaces.FiniteType.WithoutBoundary.ParentMethods.stratum(S)
    H_1(0)

While this works, this approach is quite brittle and might sometimes need a mix
with the above to work::

    sage: from flatsurf.geometry.categories import ConeSurfaces
    sage: S._refine_category_(ConeSurfaces().WithoutBoundary())
    sage: TranslationSurfaces.FiniteType.WithoutBoundary.ParentMethods.stratum(S)
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
from flatsurf.geometry.categories.euclidean_polygonal_surfaces import (
    EuclideanPolygonalSurfaces,
)
from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces
from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
from flatsurf.geometry.categories.half_translation_surfaces import (
    HalfTranslationSurfaces,
)
from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces

from flatsurf.geometry.categories.polygons import Polygons
from flatsurf.geometry.categories.euclidean_polygons import EuclideanPolygons
from flatsurf.geometry.categories.hyperbolic_polygons import HyperbolicPolygons
from flatsurf.geometry.categories.euclidean_polygons import EuclideanPolygons
from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
    EuclideanPolygonsWithAngles,
)
