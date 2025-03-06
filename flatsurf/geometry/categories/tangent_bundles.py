r"""
The category of tangent bundles over a surface.

This module provides the framework for shared functionality for all tangent
bundles and tangent vectors in sage-flatsurf. Actual functionality of tangent
bundle or vector over a surface is provided by ``TangentBundleMethods`` or
``TangentVectorMethods`` class in the category of that surface.

See
:class:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.TangentVectorMethods`
for an example.

EXAMPLES::

   sage: from flatsurf import translation_surfaces
   sage: S = translation_surfaces.square_torus()
   sage: t = S.tangent_bundle()
   sage: t.category()
   Category of tangent bundles over connected without boundary finite type translation surfaces

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2025 Julian Rüth
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
# ********************************************************************

from sage.categories.category import Category


class TangentBundles(Category):
    r"""
    The category of tangent bundles over a category of ``surfaces``.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
        sage: from flatsurf.geometry.categories.tangent_bundles import TangentBundles
        sage: S = SimilaritySurfaces()
        sage: C = TangentBundles(S)
        sage: C
        Category of tangent bundles over similarity surfaces

    TESTS::

        sage: TestSuite(C).run()

    """
    def __init__(self, surfaces):
        self._surfaces = surfaces
        super().__init__()

    def _repr_object_names(self):
        r"""
        Helper method to produce a printable representation of this category.

        EXAMPLES::

            sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
            sage: from flatsurf.geometry.categories.tangent_bundles import TangentBundles
            sage: S = SimilaritySurfaces()
            sage: C = TangentBundles(S)
            sage: C._repr_object_names()
            'tangent bundles over similarity surfaces'

        """
        return f"{super()._repr_object_names()} over {self._surfaces._repr_object_names()}"

    def super_categories(self):
        r"""
        Return the super categories of this category of tangent bundles.

        EXAMPLES::

            sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
            sage: from flatsurf.geometry.categories.tangent_bundles import TangentBundles
            sage: S = SimilaritySurfaces()
            sage: C = TangentBundles(S)
            sage: C.super_categories()
            [Category of tangent bundles over euclidean polygonal surfaces]
            sage: C = C.super_categories()[0]
            sage: C.super_categories()
            [Category of tangent bundles over polygonal surfaces]
            sage: C = C.super_categories()[0]
            sage: C.super_categories()
            [Category of tangent bundles over topological surfaces]
            sage: C = C.super_categories()[0]
            sage: C.super_categories()
            [Category of topological spaces]

        """
        from flatsurf.geometry.categories.surface_category import SurfaceCategory, SurfaceCategoryWithAxiom
        induced = [TangentBundles(surfaces) for surfaces in self._surfaces.super_categories() if isinstance(surfaces, (SurfaceCategory, SurfaceCategoryWithAxiom))]

        if induced:
            return induced

        # We cannot use the category of vector bundles because that requires a
        # fixed base ring. However, our surface categories do not have a fixed
        # base ring.
        from sage.categories.all import Sets
        return [Sets().Topological()]

    @property
    def ParentMethods(self):
        r"""
        Return the methods available to all tangent bundles in this category.

        These methods are obtained from the ``TangentBundleMethods`` in the
        corresponding category of surfaces.

        EXAMPLES::

            sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
            sage: from flatsurf.geometry.categories.tangent_bundles import TangentBundles
            sage: S = SimilaritySurfaces()
            sage: C = TangentBundles(S)
            sage: C.ParentMethods
            <class 'flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.TangentBundleMethods'>

        """
        return getattr(self._surfaces, "TangentBundleMethods", None)

    @property
    def ElementMethods(self):
        r"""
        Return the methods available to all tangent vectors for bundles in this
        category.

        These methods are obtained from the ``TangentVectorMethods`` in the
        corresponding category of surfaces.

        EXAMPLES::

            sage: from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
            sage: from flatsurf.geometry.categories.tangent_bundles import TangentBundles
            sage: S = SimilaritySurfaces()
            sage: C = TangentBundles(S)
            sage: C.ElementMethods
            <class 'flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.TangentVectorMethods'>

        """
        return getattr(self._surfaces, "TangentVectorMethods", None)
