r"""
The category of similarity surfaces.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from Euclidean polygons that are glued by similarities, i.e., identified
edges can be transformed into each other by application of rotation and
homothety (dilation) and translation.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict(QQ).category()

    sage: from flatsurf.geometry.categories import SimilaritySurfaces
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
#        Copyright (C) 2016-2020 Vincent Delecroix
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
from sage.misc.cachefunc import cached_method


class SimilaritySurfaces(Category):
    r"""
    The category of surfaces built from polygons with edges identified by
    similarities.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import SimilaritySurfaces
        sage: SimilaritySurfaces()
        Category of similarity surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.real_projective_polygonal_surfaces import RealProjectivePolygonalSurfaces
        return [RealProjectivePolygonalSurfaces()]

    class ParentMethods:
        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in by consulting
            how its edges are glued.

            The result of this method can be fed to ``_refine_category_`` to
            change the category of the surface (and enable functionality
            specific to the smaller classes of surfaces.)

            EXAMPLES::

                sage: from flatsurf import Surface_dict
                sage: S = Surface_dict(QQ)

                sage: from flatsurf import polygons
                sage: S.add_polygon(polygons.square(), label=0)
                0
                sage: S.refined_category()
                Category of compact connected finite type oriented orientable with boundary translation surfaces

                sage: S.set_edge_pairing(0, 0, 0, 2)
                sage: S.set_edge_pairing(0, 1, 0, 3)
                sage: S.refined_category()
                Category of compact connected finite type oriented orientable without boundary translation surfaces

            """
            from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces
            category = PolygonalSurfaces.ParentMethods.refined_category(self)

            # TODO: Force all surfaces to implement these methods and remove the try/except blocks.
            # TODO: Add rational
            # TODO: Document why is_… exists if there is in category.
            # TODO: Drop the NotImplemented path of the is_ for infinite type surfaces and rather put the implementation in the finite type axiom.
            try:
                cone_surface = self.is_cone_surface()
            except NotImplementedError:
                pass
            else:
                if cone_surface:
                    from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces
                    category &= ConeSurfaces()

            try:
                dilation_surface = self.is_dilation_surface()
            except NotImplementedError:
                pass
            else:
                if dilation_surface:
                    from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
                    category &= DilationSurfaces()

                    if self.is_dilation_surface(positive=True):
                        category &= DilationSurfaces().Positive()

                        try:
                            translation_surface = self.is_translation_surface()
                        except NotImplementedError:
                            pass
                        else:
                            if translation_surface:
                                from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
                                category &= TranslationSurfaces()
                    else:
                        try:
                            half_translation_surface = self.is_translation_surface(positive=False)
                        except NotImplementedError:
                            pass
                        else:
                            if half_translation_surface:
                                from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
                                category &= HalfTranslationSurfaces()

            return category

        def is_dilation_surface(self, positive=False):
            r"""
            Return whether this surface is a dilation surface, i.e., whether
            glued edges can be transformed into each other by translation
            followed by a dilation (multiplication by a diagonal matrix.)

            INPUT:

            - ``positive`` -- a boolean (default: ``False``); whether the
              entries of the diagonal matrix must be positive or are allowed to
              be negative.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_dilation_surface()
                True
                sage: S.is_dilation_surface(positive=True)
                False

            """
            if not self.is_finite():
                raise NotImplementedError("cannot decide whether this infinite type surface is a dilation surface")

            for label in self.label_iterator():
                for edge in range(self.polygon(label).num_edges()):
                    cross = self.opposite_edge(label, edge)

                    if cross is None:
                        continue

                    matrix = self.edge_matrix(label, edge)

                    if not matrix.is_diagonal():
                        return False

                    if positive:
                        if matrix[0][0] < 0 or matrix[1][1] < 0:
                            return False

            return True

        def is_translation_surface(self, positive=True):
            r"""
            Return whether this surface is a translation surface, i.e., glued
            edges can be transformed into each other by translations.

            INPUT:

            - ``positive`` -- a boolean (default: ``True``); whether the
              transformation must be a translation or is allowed to be a
              half-translation, i.e., a translation followed by a reflection in
              a point (equivalently, a rotation by π.)

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_translation_surface()
                False
                sage: S.is_translation_surface(False)
                True

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.is_translation_surface()
                True

            """
            if not self.is_finite():
                raise NotImplementedError("cannot decide whether this infinite type surface is a translation surface")

            for label in self.label_iterator():
                for edge in range(self.polygon(label).num_edges()):
                    cross = self.opposite_edge(label, edge)

                    if cross is None:
                        continue

                    matrix = self.edge_matrix(label, edge)

                    if not matrix.is_diagonal():
                        return False

                    if matrix[0][0] == 1 and matrix[1][1] == 1:
                        continue

                    if matrix[0][0] == -1 and matrix[1][1] == -1:
                        if not positive:
                            continue

                    return False

            return True

        def is_cone_surface(self):
            r"""
            Return whether this surface is a cone surface, i.e., glued edges
            can be transformed into each other with isometries.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_cone_surface()
                True

            """
            if not self.is_finite():
                raise NotImplementedError("cannot decide whether this infinite type surface is a cone surface")

            for label in self.label_iterator():
                for edge in range(self.polygon(label).num_edges()):
                    cross = self.opposite_edge(label, edge)

                    if cross is None:
                        continue

                    matrix = self.edge_matrix(label, edge)

                    if matrix * matrix.transpose() != 1:
                        return False

            return True

        def genus(self):
            r"""
            Return the genus of this surface.

            ALGORITHM:

            We use the angles around the vertices of the surface to compute the
            genus, see e.g. [Massart2021] p.17. It would probably be better to
            just compute the Euler characteristic directly from the polygon
            gluings here (and implement this on the level of polygonal
            surfaces.)

            EXAMPLES::

                sage: import flatsurf.geometry.similarity_surface_generators as sfg
                sage: sfg.translation_surfaces.octagon_and_squares().genus()
                3

                sage: from flatsurf import *
                sage: T = polygons.triangle(3,4,5)
                sage: B = similarity_surfaces.billiard(T)
                sage: B.genus()
                0
                sage: B.minimal_cover("translation").genus()
                3

            .. [Massart2021] \D. Massart. A short introduction to translation
            surfaces, Veech surfaces, and Teichműller dynamics.
            https://hal.science/hal-03300179/document

            """
            return sum(a - 1 for a in self.angles()) // 2 + 1

    class Oriented(CategoryWithAxiom):
        class ParentMethods:
            @cached_method
            def edge_matrix(self, p, e=None):
                r"""
                Returns the 2x2 matrix representing a similarity which when applied to the polygon with label `p`
                makes it so the edge `e` can be glued to its opposite edge by translation.

                If `e` is not provided, then `p` should be a pair consisting of a polygon label and an edge.

                EXAMPLES::

                    sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                    sage: s = SimilaritySurfaceGenerators.example()
                    sage: print(s.polygon(0))
                    Polygon: (0, 0), (2, -2), (2, 0)
                    sage: print(s.polygon(1))
                    Polygon: (0, 0), (2, 0), (1, 3)
                    sage: s.opposite_edge(0,0)
                    (1, 1)
                    sage: m = s.edge_matrix(0, 0)
                    sage: m
                    [   1  1/2]
                    [-1/2    1]
                    sage: m * vector((2,-2)) == -vector((-1, 3))
                    True
                """
                if e is None:
                    import warnings

                    warnings.warn("edge_matrix will now only take two arguments")
                    p, e = p
                u = self.polygon(p).edge(e)
                pp, ee = self.opposite_edge(p, e)
                v = self.polygon(pp).edge(ee)

                # note the orientation, it is -v and not v
                from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
                from sage.matrix.matrix_space import MatrixSpace

                return similarity_from_vectors(u, -v, MatrixSpace(self.base_ring(), 2))

    class Rational(CategoryWithAxiom):
        r"""
        The axiom satisfied by similarity surfaces with rational monodromy,
        i.e., a walk around a point leads to a turn by a rational multiple of
        2π.
        """

    class SubcategoryMethods:
        def Rational(self):
            return self._with_axiom("Rational")


all_axioms += ("Rational",)
