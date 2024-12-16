r"""
The category of similarity surfaces.

This module provides shared functionality for all surfaces in sage-flatsurf
that are built from Euclidean polygons that are glued by similarities, i.e.,
identified edges can be transformed into each other by application of rotation
and homothety (dilation) and translation.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: C = MutableOrientedSimilaritySurface(QQ).category()

    sage: from flatsurf.geometry.categories import SimilaritySurfaces
    sage: C.is_subcategory(SimilaritySurfaces())
    True

The easiest way to construct a similarity surface is to use the constructions
from
:class:`flatsurf.geometry.similarity_surface_generators.SimilaritySurfaceGenerators`::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
    sage: similarity_surfaces.self_glued_polygon(P)
    Half-Translation Surface in Q_0(0, -1^4) built from a quadrilateral

Another way is to build a surface from scratch (using e.g.
:class:`flatsurf.geometry.surface.MutableOrientedSimilaritySurface`)::

    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = MutableOrientedSimilaritySurface(QQ)
    sage: S.add_polygon(P)
    0
    sage: S.add_polygon(2*P)
    1
    sage: S.add_polygon(3*P)
    2
    sage: S.glue((0, 1), (1, 3))
    sage: S.glue((0, 0), (2, 2))
    sage: S.glue((0, 2), (2, 0))
    sage: S.glue((0, 3), (1, 1))
    sage: S.glue((1, 2), (2, 1))
    sage: S.glue((1, 0), (2, 3))
    sage: S
    Surface built from 3 squares

To perform a sanity check on the obtained surface, you can run its test
suite::

    sage: TestSuite(S).run()

If there are no errors reported, no consistency problems could be detected in
your surface.

Once you mark the surface as immutable, it gets more functionality, e.g.,
coming from its structure as a translation surface. This also adds more tests
to its test suite::

    sage: S.category()
    Category of finite type oriented similarity surfaces
    sage: S.set_immutable()
    sage: S.category()
    Category of connected without boundary finite type oriented rational similarity surfaces

    sage: TestSuite(S).run()

In the following example, we attempt to build a broken surface by gluing more
than two edges to each other; however, edges get unglued automatically::

    sage: S = MutableOrientedSimilaritySurface.from_surface(S)
    sage: S.glue((0, 0), (0, 3))
    sage: S.glue((0, 1), (0, 3))
    sage: S.glue((0, 2), (0, 3))

    sage: S.gluings()
    (((0, 2), (0, 3)), ((0, 3), (0, 2)), ((1, 0), (2, 3)), ((1, 2), (2, 1)), ((2, 1), (1, 2)), ((2, 3), (1, 0)))

    sage: S.set_immutable()
    sage: S.category()
    Category of with boundary finite type oriented rational similarity surfaces
    sage: TestSuite(S).run()

If we don't glue all the edges, we get a surface with boundary::

    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = MutableOrientedSimilaritySurface(QQ)
    sage: S.add_polygon(P)
    0
    sage: TestSuite(S).run()

"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2023-2024 Julian Rüth
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
# ****************************************************************************

from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)
from flatsurf.cache import cached_surface_method

from sage.categories.category_with_axiom import all_axioms
from sage.all import QQ, AA


class SimilaritySurfaces(SurfaceCategory):
    r"""
    The category of surfaces built from polygons with edges identified by
    similarities.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import SimilaritySurfaces
        sage: SimilaritySurfaces()
        Category of similarity surfaces

    """

    def super_categories(self):
        r"""
        Return the categories that a similarity surface is also a member of,
        namely the surfaces formed by Euclidean polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import SimilaritySurfaces
            sage: SimilaritySurfaces().super_categories()
            [Category of euclidean polygonal surfaces]

        """
        from flatsurf.geometry.categories.euclidean_polygonal_surfaces import (
            EuclideanPolygonalSurfaces,
        )

        return [EuclideanPolygonalSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces that are built from
        Euclidean polygons that are glued by similarities.

        If you want to add functionality for such surfaces you most likely
        want to put it here.
        """

        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in by consulting
            how its edges are glued.

            The result of this method can be fed to ``_refine_category_`` to
            change the category of the surface (and enable functionality
            specific to the smaller classes of surfaces).


            .. NOTE::

                If a surface cannot implement the various ``is_`` methods used in
                the implementation of this method (i.e., if any of them throws a
                ``NotImplementedError``,) then this method ``refined_category``
                must be overridden to skip that check. We don't want to actively
                catch a ``NotImplementedError`` and instead encourage authors
                to explicitly select the category their surfaces lives in.

            EXAMPLES::

                sage: from flatsurf import MutableOrientedSimilaritySurface
                sage: S = MutableOrientedSimilaritySurface(QQ)

                sage: from flatsurf import polygons
                sage: S.add_polygon(polygons.square(), label=0)
                0
                sage: S.refined_category()
                Category of connected with boundary finite type translation surfaces

                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 3))
                sage: S.refined_category()
                Category of connected without boundary finite type translation surfaces

            """
            from flatsurf.geometry.categories.polygonal_surfaces import (
                PolygonalSurfaces,
            )

            category = PolygonalSurfaces.ParentMethods.refined_category(self)

            if self.is_cone_surface():
                from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces

                category &= ConeSurfaces()

            if self.is_dilation_surface():
                from flatsurf.geometry.categories.dilation_surfaces import (
                    DilationSurfaces,
                )

                category &= DilationSurfaces()

                if self.is_dilation_surface(positive=True):
                    category &= DilationSurfaces().Positive()

                    if self.is_translation_surface():
                        from flatsurf.geometry.categories.translation_surfaces import (
                            TranslationSurfaces,
                        )

                        category &= TranslationSurfaces()
                elif self.is_translation_surface(positive=False):
                    from flatsurf.geometry.categories.half_translation_surfaces import (
                        HalfTranslationSurfaces,
                    )

                    category &= HalfTranslationSurfaces()

            if "Rational" not in category.axioms():
                if self.is_rational_surface():
                    category = category.Rational()

            return category

        def is_cone_surface(self):
            r"""
            Return whether this surface is a cone surface, i.e., glued edges
            can be transformed into each other with isometries.

            .. NOTE::

                This is a stronger requirement than the usual
                definition of a cone surface, see :mod:`.cone_surfaces` for
                details.

            .. NOTE::

                This method is used to determine whether this surface is in the
                category of :class:`~.cone_surfaces.ConeSurfaces`. Surfaces can
                override this method to perform specialized logic, see the note
                in :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_cone_surface()
                True

            """
            if self.is_translation_surface():
                return True

            from flatsurf.geometry.categories import ConeSurfaces

            return ConeSurfaces.ParentMethods._is_cone_surface(self)

        def is_dilation_surface(self, positive=False):
            r"""
            Return whether this surface is a dilation surface, i.e., whether
            glued edges can be transformed into each other by translation
            followed by a dilation (multiplication by a diagonal matrix).

            .. NOTE::

                This method is used to determine whether this surface is in the
                category of :class:`~.dilation_surfaces.DilationSurfaces` or
                :class:`~.dilation_surfaces.DilationSurfaces.Positive`.
                Surfaces can override this method to perform specialized logic,
                see the note in :mod:`~flatsurf.geometry.categories` for
                performance considerations.

            INPUT:

            - ``positive`` -- a boolean (default: ``False``); whether the
              entries of the diagonal matrix must be positive or are allowed to
              be negative.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_dilation_surface()
                True
                sage: S.is_dilation_surface(positive=True)
                False

            """
            if self.is_translation_surface(positive=positive):
                return True

            from flatsurf.geometry.categories import DilationSurfaces

            return DilationSurfaces.ParentMethods._is_dilation_surface(
                self, positive=positive
            )

        def is_translation_surface(self, positive=True):
            r"""
            Return whether this surface is a translation surface, i.e., glued
            edges can be transformed into each other by translations.

            This method must be implemented if this surface is a dilation surface.

            .. NOTE::

                This method is used to determine whether this surface is in the
                category of
                :class:`~.half_translation_surfaces.HalfTranslationSurfaces` or
                :class:`~.translation_surfaces.TranslationSurfaces`. Surfaces
                can override this method to perform specialized logic, see the
                note in :mod:`~flatsurf.geometry.categories` for performance
                considerations.

            INPUT:

            - ``positive`` -- a boolean (default: ``True``); whether the
              transformation must be a translation or is allowed to be a
              half-translation, i.e., a reflection about a point (equivalently,
              a rotation by π).

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
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
            from flatsurf.geometry.categories import TranslationSurfaces

            return TranslationSurfaces.ParentMethods._is_translation_surface(
                self, positive=positive
            )

        def is_rational_surface(self):
            r"""
            Return whether this surface is a rational surface, i.e., the
            rotational part of all gluings is a rational multiple of π.

            .. NOTE::

                This method is used to determine whether this surface satisfies
                the :class:`~.SimilaritySurfaces.Rational` axiom. Surfaces can
                override this method to perform specialized logic, see the note
                in :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_rational_surface()
                True

            """
            if self.is_dilation_surface(positive=False):
                return True

            return SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(self)

        def _mul_(self, matrix, switch_sides=True):
            r"""
            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: s = translation_surfaces.infinite_staircase()
                sage: s
                The infinite staircase
                sage: m=Matrix([[1,2],[0,1]])
                sage: s2=m*s
                sage: TestSuite(s2).run()
                sage: s2.polygon(0)
                Polygon(vertices=[(0, 0), (1, 0), (3, 1), (2, 1)])

            Testing multiplication by a matrix with negative determinant::

                sage: from flatsurf import dilation_surfaces
                sage: ds1 = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: ds1.polygon(0)
                Polygon(vertices=[(0, 0), (1/2, 0), (1, 1/3), (1, 1), (3/4, 1), (0, 4/5)])
                sage: m = matrix(QQ, [[0, 1], [1, 0]]) # maps (x,y) to (y, x)
                sage: ds2 = m*ds1
                sage: ds2.polygon(0)
                Polygon(vertices=[(0, 0), (4/5, 0), (1, 3/4), (1, 1), (1/3, 1), (0, 1/2)])
            """
            if not switch_sides:
                raise NotImplementedError

            from sage.structure.element import Matrix

            if not isinstance(matrix, Matrix):
                raise NotImplementedError("only implemented for matrices")
            if not matrix.dimensions != (2, 2):
                raise NotImplementedError("only implemented for 2x2 matrices")

            from flatsurf.geometry.lazy import GL2RImageSurface

            return GL2RImageSurface(self, matrix)

        def apply_matrix(self, m, in_place=None):
            r"""
            Apply the 2×2 matrix ``m`` to the polygons of this surface and
            return a morphism from this surface to the deformed surface.

            INPUT:

            - ``m`` -- a 2×2 matrix

            - ``in_place`` -- a boolean (default: ``True``); whether to modify
              this surface itself or return a modified copy of this surface
              instead.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

                sage: morphism.domain()
                Translation Surface in H_1(0) built from a square
                sage: morphism.codomain()
                Translation Surface in H_1(0) built from a rectangle

                sage: morphism.codomain().polygon(0)
                Polygon(vertices=[(0, 0), (2, 0), (2, 1), (0, 1)])

            """
            if in_place is None:
                import warnings

                warnings.warn(
                    "The defaults for apply_matrix() are going to change in a future version of sage-flatsurf; previously, apply_matrix() was performed in_place=True. "
                    "In a future version of sage-flatsurf the default is going to change to in_place=False. In the meantime, please pass in_place=True/False explicitly."
                )

                in_place = True

            if in_place:
                raise NotImplementedError(
                    "this surface does not support applying a GL(2,R) action in-place yet"
                )

            from flatsurf.geometry.lazy import GL2RImageSurface

            image = GL2RImageSurface(self, m)

            from flatsurf.geometry.morphism import GL2RMorphism

            return GL2RMorphism._create_morphism(self, image, m)

        def homology(
            self,
            k=1,
            coefficients=None,
            relative=None,
            implementation="generic",
            category=None,
        ):
            r"""
            Return the ``k``-th simplicial homology group of this surface.

            INPUT:

            - ``k`` -- an integer (default: ``1``)

            - ``coefficients`` -- a ring (default: the integer ring); consider
              the homology with coefficients in this ring

            - ``relative`` -- a set (default: the empty set); if non-empty, then
              relative homology with respect to this set is constructed.

            - ``implementation`` -- a string (default: ``"generic"``); the
              algorithm used to compute the homology groups. Currently only
              ``"generic"`` is supported, i.e., the groups are computed with
              the generic homology machinery from SageMath.

            - ``category`` -- a category; if not specified, a category for the
              homology group is chosen automatically depending on
              ``coefficients``.

            EXAMPLES::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: S.homology()
                H₁(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

                sage: S.homology(0)
                H₀(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

            """
            if self.is_mutable():
                raise ValueError("surface must be immutable to compute homology")

            from sage.all import ZZ

            k = ZZ(k)

            coefficients = coefficients or ZZ

            if category is None:
                from sage.categories.all import Modules

                category = Modules(coefficients)

            relative = frozenset(relative or {})

            return self._homology(
                k=k,
                coefficients=coefficients,
                relative=relative,
                implementation=implementation,
                category=category,
            )

        @cached_surface_method
        def _homology(self, k, coefficients, relative, implementation, category):
            r"""
            Return the ``k``-th homology group of this surface.

            This is a helper method for :meth:`homology` which passes
            normalized parameters to us. The approach chosen here is quite
            non-standard. Usually, one would not use a cached method but make
            :class:`SimplicialHomologyGroup` a unique representation. But this
            is not possible because equal surfaces can be non-identical so the
            homology of non-identical surfaces would be identical. We work
            around this issue by attaching the homology to the actual surface
            so a surface has a unique homology, but it is different from
            another equal surface's.

            A proper fix for the underlying problem could be to change
            :meth:`MutableOrientedSimilaritySurface.set_immutable` to return an
            immutable copy of the surface so that all immutable surfaces are
            unique.

            TESTS:

            Homology of a surface is unique::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: S.homology() is S.homology()
                True

            But non-identical surfaces have different homology::

                sage: from flatsurf import MutableOrientedSimilaritySurface
                sage: T = MutableOrientedSimilaritySurface.from_surface(S)
                sage: T.set_immutable()
                sage: S == T
                True
                sage: S is T
                False
                sage: S.homology() is T.homology()
                False

            """
            from flatsurf.geometry.homology import SimplicialHomologyGroup

            return SimplicialHomologyGroup(
                self, k, coefficients, relative, implementation, category
            )

        def cohomology(
            self,
            k=1,
            coefficients=None,
            relative=None,
            implementation="dual",
            category=None,
        ):
            r"""
            Return the ``k``-th simplicial cohomology group of this surface.

            INPUT:

            - ``k`` -- an integer (default: ``1``)

            - ``coefficients`` -- a ring (default: the reals); consider
              cohomology with coefficients in this ring

            - ``relative`` -- a set (default: the empty set); if non-empty, then
              relative cohomology with respect to this set is constructed.

            - ``implementation`` -- a string (default: ``"dual"``); the
              algorithm used to compute the cohomology groups. Currently only
              ``"dual"`` is supported, i.e., the groups are computed as duals
              of the generic homology groups from SageMath.

            - ``category`` -- a category; if not specified, a category for the
              cohomology group is chosen automatically depending on
              ``coefficients``.

            EXAMPLES::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: S.cohomology()
                H¹(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

            ::

                sage: S.cohomology(0)
                H⁰(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

            """
            if self.is_mutable():
                raise ValueError("surface must be immutable to compute cohomology")

            from sage.all import ZZ

            k = ZZ(k)

            from sage.all import RR

            coefficients = coefficients or RR

            relative = frozenset(relative or {})

            if category is None:
                from sage.categories.all import Modules

                category = Modules(coefficients)

            return self._cohomology(
                k=k,
                coefficients=coefficients,
                relative=relative,
                implementation=implementation,
                category=category,
            )

        @cached_surface_method
        def _cohomology(self, k, coefficients, relative, implementation, category):
            r"""
            Return the ``k``-th cohomology group of this surface.

            This is a helper method for :meth:`cohomology`, see
            :meth:`homology` for details.

            TESTS:

            Cohomology of a surface is unique::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: S.cohomology() is S.cohomology()
                True

            But non-identical surfaces have different cohomology::

                sage: from flatsurf import MutableOrientedSimilaritySurface
                sage: T = MutableOrientedSimilaritySurface.from_surface(S)
                sage: T.set_immutable()
                sage: S == T
                True
                sage: S is T
                False
                sage: S.cohomology() is T.cohomology()
                False

            """
            from flatsurf.geometry.cohomology import SimplicialCohomologyGroup

            return SimplicialCohomologyGroup(
                self, k, coefficients, relative, implementation, category
            )

    class Oriented(SurfaceCategoryWithAxiom):
        r"""
        The category of oriented surfaces built from Euclidean polygons that
        are glued by similarities with the orientation compatible with the
        orientation of the real plane that polygons are defined in.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import SimilaritySurfaces
            sage: SimilaritySurfaces().Oriented()
            Category of oriented similarity surfaces

        """

        class ElementMethods:
            r"""
            Provides methods available to all points of oriented similarity
            surfaces.

            If you want to add functionality for such surfaces, you most likely
            want to put it here.
            """

            def angle(self, numerical=False):
                r"""
                Return the total angle at this point as a multiple of 2π.

                INPUT:

                - ``numerical`` -- a boolean (default: ``False``); whether to
                  return a floating point approximation of the angle or its exact
                  value.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S(0, 0).angle()
                    1
                    sage: S(0, (1, 0)).angle()
                    1/2
                    sage: S(0, (1, 0)).angle(numerical=True)
                    0.500000000000000

                TESTS:

                For mutable surfaces, the angles are not cached::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: T = MutableOrientedSimilaritySurface.from_surface(S)
                    sage: T._refine_category_(S.category())  # make angle() available
                    sage: T(0, 0).angle()
                    1
                    sage: T(0, 0).angle() is T(0, 0).angle()
                    False

                ::

                    sage: T.glue((0, 0), (0, 2))
                    sage: T(0, 0).angle()
                    3/8

                .. SEEALSO::

                    :meth:`turns` for the number of full 2π turns at a point.

                """
                if self.is_vertex():
                    return self.parent()._angle_vertex(self, numerical=numerical)

                from sage.all import QQ, RR

                ring = RR if numerical else QQ

                if self.is_in_edge_interior():
                    label, coordinates = self.representative()
                    edge = (
                        self.parent()
                        .polygon(label)
                        .get_point_position(coordinates)
                        .get_edge()
                    )

                    opposite_edge = self.parent().opposite_edge(label, edge)
                    if opposite_edge is None or opposite_edge == (label, edge):
                        # The total angle at a self-glued or unglued edge is π.
                        return ring(0.5)

                return ring.one()

            def turns(self, start_edge=None, start_holonomy=None):
                r"""
                Return the number of total 2π turns at this vertex.

                Returns a triple ``(turns, start, end)`` where ``turns`` is an
                integer, and ``start`` and ``end`` are vectors such that start is
                the vector at which we begin counting the turns, and ``end`` is the
                vector at which we stop counting the turns.

                If the total angle at this point is a multiple of 2π, then
                ``start`` and ``end`` are parallel.

                If the surface has no boundary, then ``start`` and ``end`` are
                holonomy vectors of the same edge of the surface.

                INPUT:

                - ``start_edge`` -- an edge or ``None`` (default: ``None``); a
                  particular outgoing edge at which to start counting the
                  turns.

                - ``start_holonomy`` -- a vector or ``None`` (default: ``None``);
                  the holonomy that the first edge reported should use. If
                  ``None``, then use the holonomy that the polygon edge reports
                  normally.

                ALGORITHM:

                We count the half turns of angle π by walking around the vertex
                and counting how often the outgoing edges go between being
                clockwise and being counterclockwise from the ``start_edge``.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S(0, 0).turns()
                    (1, (0, -1), (0, -1))

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S(0, 0).turns()
                    (1, (-1, 1), (-1, 1))
                    sage: S(0, (1, 0)).turns()
                    Traceback (most recent call last):
                    ...
                    ValueError: point must be a vertex

                .. SEEALSO::

                    :meth:`angle` for the possibly non-integral number of
                    turns at this point.

                """
                if not self.is_vertex():
                    raise ValueError("point must be a vertex")

                edges = self.edges_ccw(
                    start_edge=start_edge, start_holonomy=start_holonomy
                )

                start = edges[0][1]

                half_turns = 0

                vectors = [vector for (edge, vector) in edges[1::2]]

                for vector, previous in zip(vectors, [start] + vectors):
                    from flatsurf.geometry.euclidean import ccw, is_parallel

                    # This is a vertex at a concave vertex of a polygon (with
                    # angle < 2π). So at least one π turn must have happened.
                    if ccw(previous, vector) < 0:
                        half_turns += 1

                    # Determine whether the turn from start to vector is
                    # counterclockwise (0) or clockwise (1). If they are
                    # orthogonal, we count parallel vectors as counterclockwise
                    # and anti-parallel vectors as clockwise.
                    ccw_start_vector = (
                        ccw(start, vector) or (1 if is_parallel(start, vector) else -1)
                    ) > 0
                    ccw_start_previous = half_turns % 2 == 0

                    if ccw_start_vector != ccw_start_previous:
                        half_turns += 1

                return (half_turns // 2, start, vector)

            def edges(self):
                r"""
                Return the edges of the polygons that contain this point.

                Implements
                :meth:`.polygonal_surfaces.PolygonalSurfaces.ElementMethods.edges`
                for similarity surfaces.

                """
                return {edge for (edge, holonomy) in self.edges_ccw()}

            def edges_ccw(self, start_edge=None, start_holonomy=None):
                r"""
                Return the edges of the polygons that are crossed when walking
                around this point in counterclockwise direction.

                Each edge is reported together with its holonomy vector in the
                coordinate system of the ``start_edge``.

                Each edge is reported "twice", once when leaving a polygon, and
                once when entering a polygon.

                INPUT:

                - ``start_edge`` -- an edge or ``None`` (default: ``None``); a
                  particular outgoing edge at which to start the walk.

                - ``start_holonomy`` -- a vector or ``None`` (default: ``None``);
                  the holonomy that the first edge reported should use. If
                  ``None``, then use the holonomy that the polygon edge reports
                  normally.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S(0, 0).edges_ccw()
                    [((0, 2), (-1, 1)),
                     ((0, 1), (1, -4)),
                     ((0, 1), (1, -4)),
                     ((0, 0), (2, 0)),
                     ((0, 0), (2, 0)),
                     ((0, 3), (0, 5)),
                     ((0, 3), (0, 5)),
                     ((0, 2), (-1, 1))]
                    sage: S(0, 0).edges_ccw(start_edge=(0, 0))
                    [((0, 0), (2, 0)),
                     ((0, 3), (0, 5)),
                     ((0, 3), (0, 5)),
                     ((0, 2), (-1, 1)),
                     ((0, 2), (-1, 1)),
                     ((0, 1), (1, -4)),
                     ((0, 1), (1, -4)),
                     ((0, 0), (2, 0))]
                    sage: S(0, 0).edges_ccw(start_holonomy=(1, 0))
                    [((0, 2), (1, 0)),
                     ((0, 1), (-5/2, 3/2)),
                     ((0, 1), (-5/2, 3/2)),
                     ((0, 0), (-1, -1)),
                     ((0, 0), (-1, -1)),
                     ((0, 3), (5/2, -5/2)),
                     ((0, 3), (5/2, -5/2)),
                     ((0, 2), (1, 0))]

                .. SEEALSO::

                    :meth:`edges` for the edges containing this point in no
                    particular order.

                """
                S = self.parent()

                label, coordinates = self.representative()

                position = S.polygon(label).get_point_position(coordinates)

                if position.is_in_interior():
                    return []

                if position.is_in_edge_interior():
                    return self._edges_ccw_edge(
                        start_edge=start_edge, start_holonomy=start_holonomy
                    )

                return self._edges_ccw_vertex(
                    start_edge=start_edge, start_holonomy=start_holonomy
                )

            def _test_edges_ccw(self, **options):
                r"""
                Verify that :meth:`edges_ccw` is implemented correctly.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S(0, 0)._test_edges_ccw()

                """
                tester = self._tester(**options)

                edges = self.edges_ccw()

                tester.assertTrue(len(edges) % 2 == 0)

                for e, f in zip(edges[1::2], edges[2::2]):
                    tester.assertEqual(e[1], f[1])

                for e, f in zip(edges[::2], edges[1::2]):
                    from flatsurf.geometry.euclidean import is_parallel

                    tester.assertFalse(is_parallel(e[1], f[1]))

            def _edges_ccw_edge(self, start_edge=None, start_holonomy=None):
                r"""
                Return the edges that are crossed when walking around this
                point on an edge.

                This is a helper method for :meth:`edges_ccw`.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S(0, (1/2, 0))._edges_ccw_edge()
                    [((0, 2), (-1, 0)), ((0, 0), (1, 0))]

                """
                S = self.parent()

                label, coordinates = self.representative()

                if start_edge is None:
                    position = S.polygon(label).get_point_position(coordinates)
                    start_edge = label, position.get_edge()

                if start_holonomy is None:
                    start_holonomy = S.polygon(start_edge[0]).edge(start_edge[1])
                else:
                    start_holonomy = (S.base_ring() ** 2)(start_holonomy)

                opposite = S.opposite_edge(*start_edge)
                if opposite is None or opposite == start_edge:
                    return [(start_edge, start_holonomy)]

                return [(start_edge, start_holonomy), (opposite, -start_holonomy)]

            def _edges_ccw_vertex(self, start_edge=None, start_holonomy=None):
                r"""
                Return the edges that are crossed when walking around this
                point at a vertex.

                This is a helper method for :meth:`edges_ccw`.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S(0, 0)._edges_ccw_vertex()
                    [((0, 3), (0, -1)),
                     ((0, 2), (1, 0)),
                     ((0, 0), (1, 0)),
                     ((0, 3), (0, 1)),
                     ((0, 1), (0, 1)),
                     ((0, 0), (-1, 0)),
                     ((0, 2), (-1, 0)),
                     ((0, 1), (0, -1))]

                """
                S = self.parent()

                if start_edge is None:
                    label, coordinates = self.representative()

                    initial = (
                        label,
                        S.polygon(label).get_point_position(coordinates).get_vertex(),
                    )
                    start_edge = initial

                    while True:
                        opposite = S.opposite_edge(*start_edge)
                        if opposite is None:
                            break
                        start_edge = opposite[0], (opposite[1] + 1) % len(
                            S.polygon(opposite[0]).vertices()
                        )
                        if start_edge == initial:
                            break

                if start_holonomy is None:
                    start_holonomy = S.polygon(start_edge[0]).edge(start_edge[1])
                else:
                    start_holonomy = (S.base_ring() ** 2)(start_holonomy)

                from flatsurf.geometry.similarity import similarity_from_vectors

                similarity = similarity_from_vectors(
                    S.polygon(start_edge[0]).edge(start_edge[1]), start_holonomy
                )

                edges = [(start_edge, start_holonomy)]

                while True:
                    exit_edge = (
                        edges[-1][0][0],
                        (edges[-1][0][1] - 1)
                        % len(S.polygon(edges[-1][0][0]).vertices()),
                    )
                    exit_holonomy = -similarity * S.polygon(exit_edge[0]).edge(
                        exit_edge[1]
                    )
                    edges.append((exit_edge, exit_holonomy))

                    enter_edge = S.opposite_edge(*exit_edge)
                    if enter_edge is None:
                        if S.opposite_edge(*start_edge) is not None:
                            raise ValueError(
                                "start_edge must be set to the most clockwise edge at this vertex on the boundary"
                            )
                        break

                    similarity *= S.edge_transformation(*enter_edge).derivative()

                    if enter_edge == start_edge:
                        break

                    enter_holonomy = similarity * S.polygon(enter_edge[0]).edge(
                        enter_edge[1]
                    )
                    edges.append((enter_edge, enter_holonomy))

                return edges

        class ParentMethods:
            r"""
            Provides methods available to all oriented surfaces that are built
            from Euclidean polygons that are glued by similarities.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def some_elements(self):
                r"""
                Return some typical points of this surface (for testing).

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: list(S.some_elements())
                    [Point (1/2, 1/2) of polygon 0, Vertex 0 of polygon 0]

                """
                from itertools import islice

                for label in islice(self.labels(), 32):
                    polygon = self.polygon(label)

                    from sage.categories.all import Fields

                    if polygon.is_convex() and self.base_ring() in Fields():
                        point = self(  # pylint: disable=not-callable
                            label, polygon.centroid()
                        )
                    else:
                        point = self(  # pylint: disable=not-callable
                            label, polygon.vertex(0) + polygon.edge(0) / 2
                        )

                    yield point

                if self.is_finite_type():
                    yield from self.vertices()

            def angles(self, numerical=False, return_adjacent_edges=None):
                r"""
                Return the angles around the vertices of the surface as
                multiples of 2π.

                INPUT:

                - ``numerical`` -- a boolean (default: ``False``);
                  whether the angles are returned as floating point
                  numbers or as exact algebraic numbers.

                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces, translation_surfaces
                    sage: T = polygons.triangle(3, 4, 5)
                    sage: S = similarity_surfaces.billiard(T)
                    sage: S.angles()
                    [1/4, 5/12, 1/3]
                    sage: S.angles(numerical=True)   # abs tol 1e-14
                    [0.25, 0.4166666666666667, 0.33333333333333337]

                    sage: S.angles(return_adjacent_edges=True)
                    doctest:warning
                    ...
                    UserWarning: return_adjacent_edges has been deprecated as a keyword argument to angles() and will be removed in a future version of sage-flatsurf; use angle() and edges_ccw()[::2] on each vertex instead.
                    [(1/4, [(0, 0), (1, 0)]), (5/12, [(0, 2), (1, 1)]), (1/3, [(0, 1), (1, 2)])]

                ::

                    sage: translation_surfaces.regular_octagon().angles()
                    [3]
                    sage: S = translation_surfaces.veech_2n_gon(5)
                    sage: S.angles()
                    [2, 2]
                    sage: S.angles(numerical=True)
                    [2.0, 2.0]
                    sage: S.angles(return_adjacent_edges=True) # random output
                    [(2, [(0, 1), (0, 5), (0, 9), (0, 3), (0, 7)]),
                     (2, [(0, 0), (0, 4), (0, 8), (0, 2), (0, 6)])]
                    sage: S.angles(numerical=True, return_adjacent_edges=True) # random output
                    [(2.0, [(0, 1), (0, 5), (0, 9), (0, 3), (0, 7)]),
                     (2.0, [(0, 0), (0, 4), (0, 8), (0, 2), (0, 6)])]

                    sage: translation_surfaces.veech_2n_gon(6).angles()
                    [5]
                    sage: translation_surfaces.veech_double_n_gon(5).angles()
                    [3]
                    sage: translation_surfaces.cathedral(1, 1).angles()
                    [3, 3, 3]

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                    sage: H = B.minimal_cover(cover_type="half-translation")
                    sage: S = B.minimal_cover(cover_type="translation")
                    sage: H.angles()
                    [1/2, 1/2, 1/2, 5/2]
                    sage: S.angles()
                    [1, 1, 5, 1]

                    sage: H.angles(return_adjacent_edges=True)
                    [(1/2, [...]), (1/2, [...]), (1/2, [...]), (5/2, [...])]
                    sage: S.angles(return_adjacent_edges=True)
                    [(1, [...]), (1, [...]), (5, [...]), (1, [...])]

                For self-glued edges, no angle is reported for the
                "vertex" at the midpoint of the edge::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.angles()
                    [1]

                Non-convex examples::

                    sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: L = Polygon(vertices=[(0,0),(1,0),(2,0),(2,1),(1,1),(1,2),(0,2),(0,1)])
                    sage: S.add_polygon(L)
                    0
                    sage: S.glue((0, 0), (0, 5))
                    sage: S.glue((0, 1), (0, 3))
                    sage: S.glue((0, 2), (0, 7))
                    sage: S.glue((0, 4), (0, 6))
                    sage: S.set_immutable()
                    sage: S.angles()
                    [3]

                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: P = Polygon(vertices=[(0,0),(1,0),(2,0),(2,1),(3,1),(3,2),(2,2),(1,2),(1,1),(0,1)])
                    sage: S.add_polygon(P)
                    0
                    sage: S.glue((0, 0), (0, 8))
                    sage: S.glue((0, 1), (0, 6))
                    sage: S.glue((0, 2), (0, 9))
                    sage: S.glue((0, 3), (0, 5))
                    sage: S.glue((0, 4), (0, 7))
                    sage: S.set_immutable()
                    sage: S.angles()
                    [2, 2]

                Angles can also be determined for surfaces with boundary::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: P = Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 2)])
                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(P)
                    0
                    sage: S.set_immutable()
                    sage: S.angles()
                    [1/4, 1/4, 3/8, 1/8]

                """
                if return_adjacent_edges is not None:
                    if return_adjacent_edges:
                        import warnings

                        warnings.warn(
                            "return_adjacent_edges has been deprecated as a keyword argument to angles() and will be removed in a future version of sage-flatsurf; use angle() and edges_ccw()[::2] on each vertex instead."
                        )
                        return [
                            (
                                p.angle(numerical=numerical),
                                [e for (e, h) in p.edges_ccw()[::2]],
                            )
                            for p in self.vertices()
                        ]

                    import warnings

                    warnings.warn(
                        "return_adjacent_edges has been deprecated as a keyword argument to angles() and will be removed in a future version of sage-flatsurf; omit the argument to not include the edges instead."
                    )

                return [p.angle(numerical=numerical) for p in self.vertices()]

            @cached_surface_method
            def _angle_vertex(self, vertex, numerical):
                r"""
                Return the total angle at this vertex in multiples of 2π.

                INPUT:

                - ``numerical`` -- a boolean; whether to return the angle as a
                  floating point number.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S._angle_vertex(S(0, 0), numerical=False)
                    1

                """
                turns, start, end = vertex.turns()

                from flatsurf.geometry.euclidean import angle

                return turns + angle(start, end, numerical=numerical)

            @cached_surface_method
            def edge_matrix(self, p, e=None):
                r"""
                Returns the 2x2 matrix representing a similarity which when
                applied to the polygon with label `p` makes it so the edge `e`
                can be glued to its opposite edge by translation.

                If `e` is not provided, then `p` should be a pair consisting of
                a polygon label and an edge.

                EXAMPLES::

                    sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                    sage: s = SimilaritySurfaceGenerators.example()
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
                    sage: s.polygon(1)
                    Polygon(vertices=[(0, 0), (2, 0), (1, 3)])
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

                    warnings.warn(
                        "passing only a single tuple argument to edge_matrix() has been deprecated and will be deprecated in a future version of sage-flatsurf; pass the label and edge index as separate arguments instead"
                    )
                    p, e = p

                u = self.polygon(p).edge(e)
                pp, ee = self.opposite_edge(p, e)
                v = self.polygon(pp).edge(ee)

                # note the orientation, it is -v and not v
                from flatsurf.geometry.similarity import similarity_from_vectors
                from sage.matrix.matrix_space import MatrixSpace

                return similarity_from_vectors(u, -v, MatrixSpace(self.base_ring(), 2))

            def underlying_surface(self):
                r"""
                Return this surface.

                EXAMPLES::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.underlying_surface() is S
                    doctest:warning
                    ...
                    UserWarning: underlying_surface() has been deprecated and will be removed in a future version of sage-flatsurf; this function has no effect anymore since there is no distinction between a surface and its underlying surface anymore
                    True

                """
                import warnings

                warnings.warn(
                    "underlying_surface() has been deprecated and will be removed in a future version of sage-flatsurf; this function has no effect anymore since there is no distinction between a surface and its underlying surface anymore"
                )

                return self

            def edge_transformation(self, p, e):
                r"""
                Return the similarity bringing the provided edge to the opposite edge.

                EXAMPLES::

                    sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                    sage: s = SimilaritySurfaceGenerators.example()
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
                    sage: s.polygon(1)
                    Polygon(vertices=[(0, 0), (2, 0), (1, 3)])
                    sage: s.opposite_edge(0,0)
                    (1, 1)
                    sage: g = s.edge_transformation(0,0)
                    sage: g((0,0))
                    (1, 3)
                    sage: g((2,-2))
                    (2, 0)

                """
                from flatsurf.geometry.similarity import SimilarityGroup

                G = SimilarityGroup(self.base_ring())
                q = self.polygon(p)
                a = q.vertex(e)
                b = q.vertex(e + 1)
                # This is the similarity carrying the origin to a and (1,0) to b:
                g = G(b[0] - a[0], b[1] - a[1], a[0], a[1])

                pp, ee = self.opposite_edge(p, e)
                qq = self.polygon(pp)
                # Be careful here: opposite vertices are identified
                aa = qq.vertex(ee + 1)
                bb = qq.vertex(ee)
                # This is the similarity carrying the origin to aa and (1,0) to bb:
                gg = G(bb[0] - aa[0], bb[1] - aa[1], aa[0], aa[1])

                # This is the similarity carrying (a,b) to (aa,bb):
                return gg / g

            def set_vertex_zero(self, label, v, in_place=False):
                r"""
                Applies a combinatorial rotation to the polygon with the provided label.

                This makes what is currently vertex v of this polygon vertex 0. In other words,
                what is currently vertex (or edge) e will now become vertex (e-v)%n where
                n is the number of sides of the polygon.

                For the updated polygons, the polygons will be translated so that vertex
                0 is the origin.

                EXAMPLES:

                Example with polygon glued to another polygon::

                    sage: from flatsurf import translation_surfaces
                    sage: s = translation_surfaces.veech_double_n_gon(4)
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                    sage: [s.opposite_edge(0,i) for i in range(4)]
                    [(1, 0), (1, 1), (1, 2), (1, 3)]
                    sage: ss = s.set_vertex_zero(0,1)
                    sage: ss.polygon(0)
                    Polygon(vertices=[(0, 0), (0, 1), (-1, 1), (-1, 0)])
                    sage: [ss.opposite_edge(0,i) for i in range(4)]
                    [(1, 1), (1, 2), (1, 3), (1, 0)]
                    sage: TestSuite(ss).run()

                Example with polygon glued to self::

                    sage: s = translation_surfaces.veech_2n_gon(2)
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                    sage: [s.opposite_edge(0,i) for i in range(4)]
                    [(0, 2), (0, 3), (0, 0), (0, 1)]
                    sage: ss = s.set_vertex_zero(0,3)
                    sage: ss.polygon(0)
                    Polygon(vertices=[(0, 0), (0, -1), (1, -1), (1, 0)])
                    sage: [ss.opposite_edge(0,i) for i in range(4)]
                    [(0, 2), (0, 3), (0, 0), (0, 1)]
                    sage: TestSuite(ss).run()

                """
                if in_place:
                    raise NotImplementedError(
                        "this surface does not support set_vertex_zero(mutable=True)"
                    )

                from flatsurf.geometry.surface import (
                    MutableOrientedSimilaritySurface,
                )

                s = MutableOrientedSimilaritySurface.from_surface(self)
                s.set_vertex_zero(label, v, in_place=True)
                s.set_immutable()
                return s

            def relabel(self, relabeling=None, in_place=False):
                r"""
                Return a surface whose polygons have been relabeled according
                to ``relabeling``.

                INPUT:

                - ``relabeling`` -- a dict, a callable, or ``None`` (default:
                  ``None``); the mapping from labels of this surface to labels
                  of the relabeled surface. If ``None``, then relabel to the
                  non-negative integers.

                - ``in_place`` -- a boolean (default: ``False``); whether to
                  modify this surface or return a relabeled copy instead.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.veech_double_n_gon(5)
                    sage: SS = S.relabel({0: 1, 1: 2})
                    sage: SS
                    Translation Surface in H_2(2) built from 2 regular pentagons

                    sage: SS.root()
                    1

                    sage: SS.opposite_edge(1, 0)
                    (2, 0)

                    sage: TestSuite(SS).run()

                The relabeling can also be a callable::

                    sage: SSS = SS.relabel(lambda label: label -1)
                    sage: SSS == S
                    True

                However, this is only supported on finite-type surfaces::

                    sage: S = translation_surfaces.infinite_staircase()
                    sage: S.labels()
                    (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

                    sage: S.relabel(lambda label: -label)
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: cannot relabel a surface with an infinite number of polygons with a callable relabeling yet

                ::

                    sage: S.relabel({1: 'X'})
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: cannot relabel a surface with an infinite number of polygons with a dict relabeling yet

                For infinite surfaces, we only support relabeling to the
                non-negative integers:

                    sage: S.relabel().labels()
                    (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, …)

                """
                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement relabel(in_place=True) yet"
                    )

                if not self.is_finite_type():
                    if callable(relabeling):
                        # We do not support this because the resulting surface would not be picklable.
                        # If you need this, implement a subclass of the LazyRelabeledSurface.
                        raise NotImplementedError(
                            "cannot relabel a surface with an infinite number of polygons with a callable relabeling yet"
                        )

                    if relabeling is not None:
                        # It would not be difficult to add this to
                        # LazyRelabeledSurface but so far we did not need it.
                        raise NotImplementedError(
                            "cannot relabel a surface with an infinite number of polygons with a dict relabeling yet"
                        )

                    from flatsurf.geometry.lazy import LazyRelabeledSurface

                    return LazyRelabeledSurface(self)

                from flatsurf.geometry.surface import (
                    MutableOrientedSimilaritySurface,
                )

                S = MutableOrientedSimilaritySurface.from_surface(self)
                S = S.relabel(relabeling=relabeling, in_place=True)
                S.set_immutable()

                return S

            def copy(
                self,
                relabel=False,
                mutable=False,
                lazy=None,
                new_field=None,
                optimal_number_field=False,
            ):
                r"""
                Returns a copy of this surface. The method takes several flags to modify how the copy is taken.

                If relabel is True, then instead of returning an exact copy, it returns a copy indexed by the
                non-negative integers. This uses the Surface_list implementation. If relabel is False (default),
                then we return an exact copy. The returned surface uses the Surface_dict implementation.

                The mutability flag returns if the resulting surface should be mutable or not. By default, the
                resulting surface will not be mutable.

                If lazy is True, then the surface is copied by reference. This is the only type of copy
                possible for infinite surfaces. The parameter defaults to False for finite surfaces, and
                defaults to True for infinite surfaces.

                The new_field parameter can be used to place the vertices in a larger field than the basefield
                for the original surface.

                The optimal_number_field option can be used to find a best NumberField containing the
                (necessarily finite) surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: ss=translation_surfaces.ward(3)
                    sage: ss.is_mutable()
                    False
                    sage: s=ss.copy(mutable=True)
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead.
                    sage: s.is_mutable()
                    True
                    sage: TestSuite(s).run()
                    sage: s == ss
                    False

                Changing the base field::

                    sage: s=translation_surfaces.veech_double_n_gon(5)
                    sage: ss=s.copy(mutable=False,new_field=AA)
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead.
                    Use set_immutable() to make the resulting surface immutable. Use change_ring() to change the field over which the surface is defined.
                    sage: TestSuite(ss).run()
                    sage: ss.base_ring()
                    Algebraic Real Field

                Optimization of number field::

                    sage: s = translation_surfaces.arnoux_yoccoz(3)
                    sage: ss = s.copy(new_field=AA).copy(optimal_number_field=True)
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead.
                    Use set_immutable() to make the resulting surface immutable. Use change_ring() to change the field over which the surface is defined.
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead.
                    Use set_immutable() to make the resulting surface immutable. There is currently no replacement for optimal number field.
                    If you are relying on this features, let the authors of sage-flatsurf know and we will try to make it available again.
                    sage: TestSuite(ss).run()
                    sage: ss.base_ring().discriminant()
                    -44
                """
                message = "copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead."

                if not mutable:
                    message += (
                        " Use set_immutable() to make the resulting surface immutable."
                    )

                if relabel:
                    message += " Use relabel({old: new for (new, old) in enumerate(surface.labels())}) for integer labels."

                if not self.is_finite_type():
                    message += " However, there is no immediate replacement for lazy copying of infinite surfaces. Have a look at the implementation of flatsurf.geometry.delaunay.LazyMutableSurface and adapt it to your needs."

                if new_field is not None:
                    message += " Use change_ring() to change the field over which the surface is defined."

                if optimal_number_field:
                    message += " There is currently no replacement for optimal number field. If you are relying on this features, let the authors of sage-flatsurf know and we will try to make it available again."

                import warnings

                warnings.warn(message)

                category = self.category()
                s = None  # This will be the surface we copy. (Likely we will set s=self below.)
                if new_field is not None and optimal_number_field:
                    raise ValueError(
                        "You can not set a new_field and also set optimal_number_field=True."
                    )
                if optimal_number_field is True:
                    if not self.is_finite_type():
                        raise NotImplementedError(
                            "can only optimize_number_field for a finite surface"
                        )
                    if lazy:
                        raise NotImplementedError(
                            "lazy copying is unavailable when optimize_number_field=True"
                        )
                    coordinates_AA = []
                    for label, p in zip(self.labels(), self.polygons()):
                        for e in p.edges():
                            coordinates_AA.append(AA(e[0]))
                            coordinates_AA.append(AA(e[1]))
                    from sage.rings.qqbar import number_field_elements_from_algebraics

                    field, coordinates_NF, hom = number_field_elements_from_algebraics(
                        coordinates_AA, minimal=True
                    )
                    if field is QQ:
                        new_field = QQ
                        # We pretend new_field = QQ was passed as a parameter.
                        # It will now get picked up by the "if new_field is not None:" line below.
                    else:
                        # Unfortunately field doesn't come with an real embedding (which is given by hom!)
                        # So, we make a copy of the field, and add the embedding.
                        from sage.all import NumberField

                        field2 = NumberField(
                            field.polynomial(), name="a", embedding=hom(field.gen())
                        )
                        # The following converts from field to field2:
                        hom2 = field.hom(im_gens=[field2.gen()])

                        from flatsurf.geometry.surface import (
                            MutableOrientedSimilaritySurface,
                        )

                        ss = MutableOrientedSimilaritySurface(field2)
                        index = 0

                        from flatsurf import Polygon

                        for label, p in zip(self.labels(), self.polygons()):
                            new_edges = []
                            for i in range(len(p.vertices())):
                                new_edges.append(
                                    (
                                        hom2(coordinates_NF[index]),
                                        hom2(coordinates_NF[index + 1]),
                                    )
                                )
                                index += 2
                            pp = Polygon(edges=new_edges, base_ring=field2)
                            ss.add_polygon(pp, label=label)
                        ss.set_roots(self.roots())
                        for (l1, e1), (l2, e2) in self.gluings():
                            ss.glue((l1, e1), (l2, e2))
                        s = ss
                        if not relabel:
                            if not mutable:
                                s.set_immutable()
                            return s
                        # Otherwise we are supposed to relabel. We will make a relabeled copy of s below.
                if new_field is not None:
                    s = self.change_ring(new_field)
                if s is None:
                    s = self
                if s.is_finite_type():
                    if relabel:
                        from flatsurf.geometry.surface import Surface_list

                        return Surface_list(
                            surface=s,
                            copy=not lazy,
                            mutable=mutable,
                            category=category,
                            deprecation_warning=False,
                        )
                    else:
                        from flatsurf.geometry.surface import Surface_dict

                        return Surface_dict(
                            surface=s,
                            copy=not lazy,
                            mutable=mutable,
                            category=category,
                            deprecation_warning=False,
                        )
                else:
                    if lazy is False:
                        raise ValueError(
                            "Only lazy copying available for infinite surfaces."
                        )
                    if self.is_mutable():
                        raise ValueError(
                            "An infinite surface can only be copied if it is immutable."
                        )
                    if relabel:
                        from flatsurf.geometry.surface import Surface_list

                        return Surface_list(
                            surface=s,
                            copy=False,
                            mutable=mutable,
                            category=category,
                            deprecation_warning=False,
                        )
                    else:
                        from flatsurf.geometry.surface import Surface_dict

                        return Surface_dict(
                            surface=s,
                            copy=False,
                            mutable=mutable,
                            category=category,
                            deprecation_warning=False,
                        )

            def change_ring(self, ring):
                r"""
                Return a copy of this surface whose polygons are defined over
                ``ring``.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.veech_2n_gon(4)
                    sage: T = S.change_ring(AA)
                    sage: T.base_ring()
                    Algebraic Real Field

                This method creates a copy even if the base ring is unchanged::

                    sage: S.change_ring(S.base_ring()) is S
                    False

                Note that the resulting surface might not be functional if the
                polygon cannot be represented in the ring::

                    sage: S.change_ring(QQ).polygon(0)
                    Traceback (most recent call last):
                    ...
                    TypeError: Unable to coerce 1/2*a + 1 to a rational

                """
                from flatsurf.geometry.surface import BaseRingChangedSurface

                return BaseRingChangedSurface(self, ring)

            def triangle_flip(
                self, label, edge, in_place=False, test=None, direction=None
            ):
                r"""
                Flips the diagonal of the quadrilateral formed by two triangles
                glued together along the provided ``edge`` of the polygon with
                ``label``.

                The diagonal is flipped by turning it counterclockwise (from
                the point of view of the polygon with ``label``.) The ``label``
                is given to the triangle which has the edge that was previously
                the edge preceding ``edge`` in the polygon. The diagonal does
                not change its index in the triangles.

                In other words,

                .. code-block:: text

                   +<----------+       +<----------+
                   |\          ^       |          /^
                   | \      L' |       | L'      / |
                   |  \e'      |       |        /  |
                   |   \|      |       |       /e  |
                   |e+1 \      |       |      /-   |
                   |     \     |       |     /     |
                   |      \    |       |   -/      |
                   |       \   |       | e'/   L   |
                   |  L    |\  |       |  /        |
                   |        e\ |       | /      e-1|
                   v e-1      \|       v/  e+1     |
                   +---------->+       +---------->+

                where ``L`` is ``label``, ``L'`` is the label of the opposite
                polygon, ``e`` is ``edge`` and ``e'`` is the edge index of the
                diagonal in ``L'``.

                INPUT:

                - ``label`` -- a polygon label

                - ``edge`` -- an int; the index of an edge of the triangle with
                  ``label``

                - ``in_place`` -- a bool (default: ``False``) - whether this
                  surface is modified or a modified copy is produced instead.

                - ``test`` -- a bool (default: ``False``) - whether not to
                  perform the flip but return whether the flip is possible.

                - ``direction`` -- ignored

                EXAMPLES::

                    sage: from flatsurf import similarity_surfaces, MutableOrientedSimilaritySurface

                    sage: S = similarity_surfaces.right_angle_triangle(1, 1)
                    sage: S.triangle_flip(0, 0)
                    Traceback (most recent call last):
                    ...
                    ValueError: cannot flip this edge because the surrounding quadrilateral is not strictly convex

                    sage: S.triangle_flip(0, 1)
                    Genus 0 Rational Cone Surface built from 2 isosceles triangles

                    sage: S.triangle_flip(0, 2)
                    Traceback (most recent call last):
                    ...
                    ValueError: cannot flip this edge because the surrounding quadrilateral is not strictly convex

                We can flip edges of self-glued polygons::

                    sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface, similarity_surfaces

                    sage: P = Polygon(edges=[(2, 0), (-1, 3), (-1, -3)])
                    sage: T = MutableOrientedSimilaritySurface.from_surface(similarity_surfaces.self_glued_polygon(P))

                    sage: T.triangle_flip(0, 1, in_place=True)
                    Half-Translation Surface built from a triangle

                    sage: T.set_immutable()

                    sage: from flatsurf.geometry.categories import DilationSurfaces
                    sage: T in DilationSurfaces()
                    True
                    sage: T.labels()
                    (0,)
                    sage: T.polygons()
                    (Polygon(vertices=[(0, 0), (-1, -3), (2, 0)]),)
                    sage: T.gluings()
                    (((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)))

                TESTS::

                    sage: TestSuite(S).run()
                    sage: TestSuite(T).run()

                The deprecated ``test`` keyword argument::

                    sage: S = similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
                    sage: S.polygon(0)
                    Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
                    sage: S.triangle_flip(0, 0, test=True)
                    doctest:warning
                    ...
                    UserWarning: the test keyword has been deprecated in triangle_flip and will be removed in a future version of sage-flatsurf; is is_convex(strict=True) instead.
                    False
                    sage: S.is_convex(0, 0, strict=True)
                    False

                    sage: S.triangle_flip(0, 1, test=True)
                    True
                    sage: S.triangle_flip(0, 2, test=True)
                    False

                """
                if test is not None:
                    if not test:
                        import warnings

                        warnings.warn(
                            "the test keyword has been deprecated in triangle_flip and will be removed in a future version of sage-flatsurf; it should not be passed in anymore"
                        )
                    else:
                        import warnings

                        warnings.warn(
                            "the test keyword has been deprecated in triangle_flip and will be removed in a future version of sage-flatsurf; is is_convex(strict=True) instead."
                        )
                        if len(self.polygon(label).vertices()) != 3:
                            return False
                        if (
                            len(
                                self.polygon(
                                    self.opposite_edge(label, edge)[0]
                                ).vertices()
                            )
                            != 3
                        ):
                            return False
                        return self.is_convex(label, edge, strict=True)

                if direction is not None:
                    import warnings

                    warnings.warn(
                        "the direction keyword has been removed from triangle_flip(); the diagonal of the quadrilateral is now always turned counterclockwise"
                    )

                if in_place:
                    raise NotImplementedError(
                        "this surface does not support triangle_flip(in_place=True) yet"
                    )

                from flatsurf.geometry.surface import (
                    MutableOrientedSimilaritySurface,
                )

                s = MutableOrientedSimilaritySurface.from_surface(self)
                s.triangle_flip(label, edge, in_place=True)
                s.set_immutable()

                return s

            def is_convex(self, label, edge, strict=False):
                r"""
                Return whether the polygon with ``label`` is convex after it
                has been glued with the polygon across its ``edge``.

                INPUT:

                - ``label`` -- a label of a polygon in this surface

                - ``edge`` -- the index of an edge of the polygon with ``label``

                - ``strict`` -- a boolean (default: ``False``); whether to
                  determine that the polygon is strictly convex, i.e., none of
                  its inner angles greater or equal than π.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
                    sage: S.is_convex(0, 0)
                    True
                    sage: S.is_convex(0, 0, strict=True)
                    False

                    sage: S.is_convex(2, 0)
                    True

                """
                opposite_label, opposite_edge = self.opposite_edge(label, edge)

                p1 = self.polygon(label)

                from flatsurf import Polygon

                sim = self.edge_transformation(opposite_label, opposite_edge)
                p2 = self.polygon(opposite_label)
                p2 = Polygon(
                    vertices=[sim(v) for v in p2.vertices()], base_ring=p1.base_ring()
                )

                p = p1.join(p2, edge, opposite_edge)

                return p.is_convex(strict=strict)

            def random_flip(self, repeat=1, in_place=False):
                r"""
                Perform random edge flip on a triangulated surface.

                INPUT:

                - ``repeat`` -- integer (default ``1``). The number of edge flip to perform.

                - ``in_place`` -- whether the transformation is done in place.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: ss = translation_surfaces.ward(3).triangulate()
                    sage: ss.random_flip(15)  # random
                    Translation Surface in H_1(0^3) built from 6 triangles

                """
                if not self.is_triangulated():
                    raise ValueError("random_flip only works for triangulated surfaces")
                if not in_place:
                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    self = MutableOrientedSimilaritySurface.from_surface(self)
                labels = list(self.labels())
                i = 0
                from sage.misc.prandom import choice

                while i < repeat:
                    p1 = choice(labels)
                    e1 = choice(range(3))
                    if not self.is_convex(p1, e1, strict=True):
                        continue
                    self.triangle_flip(p1, e1, in_place=True)
                    i += 1
                if not in_place:
                    self.set_immutable()
                return self

            def join_polygons(self, p1, e1, test=False, in_place=False):
                r"""
                Join polygons across the provided edge (p1,e1). By default,
                it returns the surface obtained by joining the two polygons
                together. It raises a ValueError if gluing the two polygons
                together results in a non-convex polygon. This is done to the
                current surface if in_place is True, and otherwise a mutable
                copy is made and then modified.

                If test is True then instead of changing the surface, it just
                checks to see if the change would be successful and returns
                True if successful or False if not.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
                    sage: ss = translation_surfaces.ward(3)
                    sage: s = MutableOrientedSimilaritySurface.from_surface(ss)
                    sage: s.join_polygons(0,0, in_place=True)
                    Translation Surface built from an equilateral triangle and a pentagon with 2 marked vertices
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (0, 2*a), (-1, a)])
                    sage: s.join_polygons(0,4, in_place=True)
                    Translation Surface built from a rhombus
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (1, 3*a), (0, 2*a), (-1, a)])

                TESTS::

                    sage: from flatsurf.geometry.categories import TranslationSurfaces
                    sage: s.set_immutable()
                    sage: s in TranslationSurfaces()
                    True

                """
                if test:
                    in_place = False

                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement join_polygons(in_place=True) yet"
                    )

                if not test:
                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    s = MutableOrientedSimilaritySurface.from_surface(self)
                    s.join_polygons(p1=p1, e1=e1, test=False, in_place=True)
                    s.set_immutable()
                    return s

                poly1 = self.polygon(p1)
                p2, e2 = self.opposite_edge(p1, e1)
                poly2 = self.polygon(p2)

                if p1 == p2:
                    return False

                t = self.edge_transformation(p2, e2)
                dt = t.derivative()
                es = []
                for i in range(e1):
                    es.append(poly1.edge(i))
                ne = len(poly2.vertices())
                for i in range(1, ne):
                    ee = (e2 + i) % ne
                    es.append(dt * poly2.edge(ee))
                for i in range(e1 + 1, len(poly1.vertices())):
                    es.append(poly1.edge(i))

                try:
                    from flatsurf import Polygon

                    Polygon(edges=es, base_ring=self.base_ring())
                except (ValueError, TypeError):
                    return False

                # Gluing would be successful
                return True

            def subdivide_polygon(self, p, v1, v2, test=False, new_label=None):
                r"""
                Cut the polygon with label p along the diagonal joining vertex
                v1 to vertex v2. This cuts p into two polygons, one will keep the same
                label. The other will get a new label, which can be provided
                via new_label. Otherwise a default new label will be provided.
                If test=False, then the surface will be changed (in place). If
                test=True, then it just checks to see if the change would be successful

                The convention is that the resulting subdivided polygon which has an oriented
                edge going from the original vertex v1 to vertex v2 will keep the label p.
                The other polygon will get a new label.

                The change will be done in place.
                """
                if not test:
                    raise NotImplementedError(
                        "this surface does not implement subdivide_polygon(test=False) yet"
                    )

                poly = self.polygon(p)
                ne = len(poly.vertices())
                if v1 < 0 or v2 < 0 or v1 >= ne or v2 >= ne:
                    return False
                if abs(v1 - v2) <= 1 or abs(v1 - v2) >= ne - 1:
                    return False

                return True

            def singularity(self, label, v, limit=None):
                r"""
                Represents the Singularity associated to the v-th vertex of the polygon
                with label ``label``.

                If the surface is infinite, the limit can be set. In this case the
                construction of the singularity is successful if the sequence of
                vertices hit by passing through edges closes up in limit or less steps.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: s = translation_surfaces.square_torus()
                    sage: pc = s.minimal_cover(cover_type="planar")
                    sage: pc.singularity(pc.root(), 0)
                    doctest:warning
                    ...
                    UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
                    Vertex 0 of polygon (0, (x, y) |-> (x, y))
                    sage: pc.singularity(pc.root(), 0, limit=1)
                    Traceback (most recent call last):
                    ...
                    ValueError: number of edges at singularity exceeds limit

                """
                from flatsurf.geometry.surface_objects import Singularity

                return Singularity(self, label, v, limit)

            def point(self, label, point, ring=None, limit=None):
                r"""
                Return a point in this surface.

                INPUT:

                - ``label`` - label of the polygon

                - ``point`` - coordinates of the point inside the polygon or
                  the index of the vertex of the polygon.

                - ``ring`` (optional) - a ring for the coordinates

                - ``limit`` (optional) - undocumented (only relevant if the point
                  corresponds to a singularity in an infinite surface)

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: s = translation_surfaces.square_torus()
                    sage: pc = s.minimal_cover(cover_type="planar")
                    sage: pc.point(pc.root(), (0, 0))
                    Vertex 0 of polygon (0, (x, y) |-> (x, y))
                    sage: pc.point(pc.root(), 0)
                    Vertex 0 of polygon (0, (x, y) |-> (x, y))
                    sage: pc.point(pc.root(), 1)
                    Vertex 0 of polygon (0, (x, y) |-> (x + 1, y))
                    sage: pc.point(pc.root(), (1, 1))
                    Vertex 0 of polygon (0, (x, y) |-> (x + 1, y + 1))
                    sage: z = pc.point(pc.root(),(sqrt(2)-1,sqrt(3)-1),ring=AA)
                    doctest:warning
                    ...
                    UserWarning: the ring parameter is deprecated and will be removed in a future version of sage-flatsurf; define the surface over a larger ring instead so that this points' coordinates live in the base ring
                    sage: next(iter(z.coordinates(next(iter(z.labels()))))).parent()
                    Vector space of dimension 2 over Algebraic Real Field

                ::

                    sage: s = translation_surfaces.cathedral(2, 3)
                    sage: s.point(0, 0)
                    Vertex 0 of polygon 0
                    sage: s.point(0, (0, 0))
                    Vertex 0 of polygon 0
                    sage: s.point(0, (1, 1))
                    Point (1, 0) of polygon 0
                    sage: s.point(0, 1)
                    Vertex 0 of polygon 1

                """
                # pylint: disable-next=not-callable
                return self(label, point, limit=limit, ring=ring)

            def surface_point(self, *args, **kwargs):
                r"""
                Return a point in this surface.

                This is an alias for :meth:`point`.
                """
                import warnings

                warnings.warn(
                    "surface_point() has been deprecated and will be removed in a future version of sage-flatsurf; use point() instead"
                )

                return self.point(*args, **kwargs)

            def minimal_cover(self, cover_type="translation"):
                r"""
                Return the minimal translation or half-translation cover of the surface.

                Cover type may be either "translation", "half-translation" or "planar".

                The minimal planar cover of a surface S is the smallest cover C so that
                the developing map from the universal cover U to the plane induces a
                well defined map from C to the plane. This is an infinite translation
                surface that is naturally a branched cover of the plane.

                EXAMPLES::

                    sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
                    sage: s = MutableOrientedSimilaritySurface(QQ)
                    sage: square = polygons.square(base_ring=QQ)
                    sage: s.add_polygon(square)
                    0
                    sage: s.glue((0,0), (0,1))
                    sage: s.glue((0,2) ,(0,3))
                    sage: cs = s
                    sage: ts = cs.minimal_cover(cover_type="translation")
                    sage: ts
                    Minimal Translation Cover of Rational Cone Surface built from a square
                    sage: from flatsurf.geometry.categories import TranslationSurfaces
                    sage: ts in TranslationSurfaces()
                    True
                    sage: hts = cs.minimal_cover(cover_type="half-translation")
                    sage: hts
                    Minimal Half-Translation Cover of Genus 0 Rational Cone Surface built from a square
                    sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
                    sage: hts in HalfTranslationSurfaces()
                    True
                    sage: TestSuite(hts).run()
                    sage: ps = cs.minimal_cover(cover_type="planar")
                    sage: ps
                    Minimal Planar Cover of Genus 0 Rational Cone Surface built from a square
                    sage: ps in TranslationSurfaces()
                    True
                    sage: TestSuite(ps).run()

                    sage: from flatsurf import similarity_surfaces
                    sage: S = similarity_surfaces.example()
                    sage: T = S.minimal_cover(cover_type="translation")
                    sage: T
                    Minimal Translation Cover of Genus 1 Surface built from 2 isosceles triangles
                    sage: T in TranslationSurfaces()
                    True
                    sage: T.polygon(T.root())
                    Polygon(vertices=[(0, 0), (2, -2), (2, 0)])

                """
                if cover_type == "translation":
                    from flatsurf.geometry.minimal_cover import MinimalTranslationCover

                    return MinimalTranslationCover(self)

                if cover_type == "half-translation":
                    from flatsurf.geometry.minimal_cover import (
                        MinimalHalfTranslationCover,
                    )

                    return MinimalHalfTranslationCover(self)

                if cover_type == "planar":
                    from flatsurf.geometry.minimal_cover import MinimalPlanarCover

                    return MinimalPlanarCover(self)

                raise ValueError("Provided cover_type is not supported.")

            def vector_space(self):
                r"""
                Return the vector space in which self naturally embeds.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S.vector_space()
                    doctest:warning
                    ...
                    UserWarning: vector_space() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring()**2 or base_ring().fraction_field()**2 instead
                    Vector space of dimension 2 over Rational Field

                    sage: S.base_ring()**2
                    Vector space of dimension 2 over Rational Field

                """
                import warnings

                warnings.warn(
                    "vector_space() has been deprecated and will be removed in a future version of sage-flatsurf; use base_ring()**2 or base_ring().fraction_field()**2 instead"
                )

                from sage.modules.free_module import VectorSpace

                return VectorSpace(self.base_ring(), 2)

            @cached_surface_method(key=lambda self, ring: ring or self.base_ring())
            def tangent_bundle(self, ring=None):
                r"""
                Return the tangent bundle

                INPUT:

                - ``ring`` -- an optional field (defaults to the coordinate field of the
                  surface)
                """
                if ring is None:
                    ring = self.base_ring()

                if self.is_mutable():
                    raise NotImplementedError(
                        "cannot compute the tangent bundle of a mutable surface"
                    )

                from flatsurf.geometry.tangent_bundle import (
                    SimilaritySurfaceTangentBundle,
                )

                return SimilaritySurfaceTangentBundle(self, ring)

            def tangent_vector(self, lab, p, v, ring=None):
                r"""
                Return a tangent vector.

                INPUT:

                - ``lab`` -- label of a polygon

                - ``p`` -- coordinates of a point in the polygon

                - ``v`` -- coordinates of a vector in R^2

                EXAMPLES::

                    sage: from flatsurf.geometry.chamanara import chamanara_surface
                    sage: S = chamanara_surface(1/2)
                    sage: S.tangent_vector(S.root(), (1/2,1/2), (1,1))
                    SimilaritySurfaceTangentVector in polygon 1 based at (1/2, -3/2) with vector (1, 1)
                    sage: K.<sqrt2> = QuadraticField(2)
                    sage: S.tangent_vector(S.root(), (1/2,1/2), (1,sqrt2), ring=K)
                    SimilaritySurfaceTangentVector in polygon 1 based at (1/2, -3/2) with vector (1, sqrt2)
                """
                from sage.all import vector

                p = vector(p)
                v = vector(v)

                if p.parent().dimension() != 2 or v.parent().dimension() != 2:
                    raise ValueError(
                        "p (={!r}) and v (={!v}) should have two coordinates"
                    )

                return self.tangent_bundle(ring=ring)(lab, p, v)

            def triangulation_mapping(self):
                r"""
                Return a ``SurfaceMapping`` triangulating the surface
                or ``None`` if the surface is already triangulated.
                """
                from flatsurf.geometry.mappings import triangulation_mapping

                return triangulation_mapping(self)

            def triangulate(self, in_place=False, label=None, relabel=None):
                r"""
                Return a triangulated version of this surface. (This may be mutable
                or not depending on the input.)

                INPUT:

                - ``in_place`` -- a boolean (default: ``False``); whether to
                  modify this surface or return a triangulated copy instead.

                - ``label`` -- a label or ``None`` (default: ``None``); if set,
                  then only the polygon with that label is triangulated and all
                  other polygons are unchanged. Otherwise, all polygons are
                  triangulated.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
                    sage: S.triangulate()
                    Triangulation morphism:
                      From: Translation Surface in H_2(2) built from 3 squares
                      To:   Triangulation of Translation Surface in H_2(2) built from 3 squares

                A non-strictly convex example that caused trouble at some point:

                    sage: from flatsurf import similarity_surfaces, Polygon
                    sage: P = Polygon(edges=[(1, 1), (-3, -1), (1, 0), (1, 0)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S
                    Half-Translation Surface in Q_0(0, -1^4) built from a triangle

                    sage: T = S.triangulate().codomain()
                    sage: len(T.polygon((0, 0)).vertices())
                    3

                The surface returned is explicitly a triangulation of the
                original surface. Use
                :meth:`flatsurf.geometry.surface.MutableOrientedSimilaritySurface.from_surface`
                or :meth:`relabel` to get a primitive surface::

                    sage: T
                    Triangulation of Half-Translation Surface in Q_0(0, -1^4) built from a triangle
                    sage: T.relabel()
                    Half-Translation Surface in Q_0(0, -1^4) built from 2 triangles

                ::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: T = MutableOrientedSimilaritySurface.from_surface(T)
                    sage: T.set_immutable()
                    sage: T
                    Half-Translation Surface in Q_0(0, -1^4) built from 2 triangles

                We can also only triangulate part of a surface, namely a single
                polygon::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.cathedral(1, 1)
                    sage: S
                    Translation Surface in H_4(2^3) built from 2 squares, a hexagon with 4 marked vertices and an octagon
                    sage: S.triangulate(label=1).codomain().relabel()
                    Translation Surface in H_4(2^3) built from 2 isosceles triangles, 5 triangles, a right triangle, 2 squares and an octagon

                """
                if relabel is not None:
                    import warnings

                    warnings.warn(
                        "the relabel keyword argument of triangulate() is ignored, it has been deprecated and will be removed in a future version of sage-flatsurf"
                    )

                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement triangulate(in_place=True) yet"
                    )

                if self.is_mutable():
                    from flatsurf import MutableOrientedSimilaritySurface

                    s = MutableOrientedSimilaritySurface.from_surface(self)
                    s.set_immutable()
                    triangulation = s.triangulate(label=label, relabel=relabel)

                    from flatsurf.geometry.morphism import NamedUnknownMorphism

                    return NamedUnknownMorphism._create_morphism(
                        self, triangulation.codomain(), "Triangulation"
                    )

                if label is not None:
                    label = {label}

                from flatsurf.geometry.morphism import (
                    TriangulationMorphism_LazyTriangulatedSurface,
                )
                from flatsurf.geometry.lazy import LazyTriangulatedSurface

                return TriangulationMorphism_LazyTriangulatedSurface._create_morphism(
                    self, LazyTriangulatedSurface(self, labels=label)
                )

            def _delaunay_edge_needs_flip(self, p1, e1):
                r"""
                Return whether edge ``e1`` of polygon ``p1`` should be flipped
                to get closer to a Delaunay triangulated surface.
                """
                p2, e2 = self.opposite_edge(p1, e1)
                poly1 = self.polygon(p1)
                poly2 = self.polygon(p2)
                if len(poly1.vertices()) != 3 or len(poly2.vertices()) != 3:
                    raise ValueError("Edge must be adjacent to two triangles.")
                from flatsurf.geometry.similarity import similarity_from_vectors

                sim1 = similarity_from_vectors(poly1.edge(e1 + 2), -poly1.edge(e1 + 1))
                sim2 = similarity_from_vectors(poly2.edge(e2 + 2), -poly2.edge(e2 + 1))
                sim = sim1 * sim2
                return sim[1][0] < 0

            def _delaunay_edge_needs_join(self, p1, e1):
                r"""
                Return whether edge ``e1`` of polygon ``p1`` should be
                eliminated and the polygons attached to it joined to get closer
                to a Delaunay cell decomposition.
                """
                p2, e2 = self.opposite_edge(p1, e1)
                poly1 = self.polygon(p1)
                poly2 = self.polygon(p2)
                from flatsurf.geometry.similarity import similarity_from_vectors

                sim1 = similarity_from_vectors(
                    poly1.vertex(e1) - poly1.vertex(e1 + 2), -poly1.edge(e1 + 1)
                )
                sim2 = similarity_from_vectors(
                    poly2.vertex(e2) - poly2.vertex(e2 + 2), -poly2.edge(e2 + 1)
                )
                sim = sim1 * sim2

                return sim[1][0] == 0

            def is_delaunay_triangulated(self, limit=None):
                r"""
                Return whether the surface is triangulated and the
                triangulation is Delaunay.

                INPUT:

                - ``limit`` -- an integer or ``None`` (default: ``None``);
                  check only ``limit`` many edges if set

                """
                if limit is not None:
                    import warnings

                    warnings.warn(
                        "limit has been deprecated as a keyword argument for is_delaunay_triangulated() and will be removed from a future version of sage-flatsurf; "
                        "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                    )

                if not self.is_finite_type() and limit is None:
                    raise NotImplementedError(
                        "a limit must be set for infinite surfaces."
                    )

                count = 0

                for (l1, e1), (l2, e2) in self.gluings():
                    if limit is not None and count >= limit:
                        break
                    count += 1
                    if len(self.polygon(l1).vertices()) != 3:
                        return False
                    if len(self.polygon(l2).vertices()) != 3:
                        return False
                    if self._delaunay_edge_needs_flip(l1, e1):
                        return False

                return True

            def is_delaunay_decomposed(self, limit=None):
                r"""
                Return if the decomposition of the surface into polygons is Delaunay.

                INPUT:

                - ``limit`` -- an integer or ``None`` (default: ``None``);
                  check only ``limit`` many polygons if set

                """
                if limit is not None:
                    import warnings

                    warnings.warn(
                        "limit has been deprecated as a keyword argument for is_delaunay_decomposed() and will be removed from a future version of sage-flatsurf; "
                        "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                    )

                if not self.is_finite_type() and limit is None:
                    raise NotImplementedError(
                        "a limit must be set for infinite surfaces."
                    )

                count = 0

                for l1, p1 in zip(self.labels(), self.polygons()):
                    if limit is not None and count >= limit:
                        break

                    count += 1

                    try:
                        c1 = p1.circumscribed_circle()
                    except ValueError:
                        # p1 is not circumscribed
                        return False

                    for e1 in range(len(p1.vertices())):
                        c2 = self.edge_transformation(l1, e1) * c1
                        l2, e2 = self.opposite_edge(l1, e1)
                        if c2.point_position(self.polygon(l2).vertex(e2 + 2)) != -1:
                            # The circumscribed circle developed into the adjacent polygon
                            # contains a vertex in its interior or boundary.
                            return False

                return True

            def delaunay_triangulation(
                self,
                triangulated=None,
                in_place=False,
                direction=None,
                relabel=None,
            ):
                r"""
                Return a Delaunay triangulation of this surface.

                INPUT:

                - ``triangulated`` - ``None``; ignored

                - ``in_place`` - boolean (default: ``False``); whether this
                  ssurface is modified and returned or just a copy of the
                  surface

                - ``direction`` - ``None``; ignored

                - ``relabel`` - ``None``; ignored

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces

                    sage: m = matrix([[2, 1], [1, 1]])
                    sage: s = m * translation_surfaces.infinite_staircase()
                    sage: ss = s.delaunay_triangulation()
                    doctest:warning
                    ...
                    UserWarning: delaunay_triangulation() is deprecated and will be removed in a future version of sage-flatsurf; use delaunay_triangulate().codomain() instead
                    sage: ss.root()
                    (0, 0)
                    sage: ss.polygon((0, 0))
                    Polygon(vertices=[(0, 0), (-1, 0), (-1, -1)])
                    sage: TestSuite(ss).run()
                    sage: ss.is_delaunay_triangulated()
                    True

                """
                import warnings

                warnings.warn(
                    "delaunay_triangulation() is deprecated and will be removed in a future version of sage-flatsurf; use delaunay_triangulate().codomain() instead"
                )

                if triangulated is not None:
                    import warnings

                    warnings.warn(
                        "the triangulated keyword has been deprecated as a parameter for delaunay_triangulation() and will be removed in a future version of sage-flatsurf"
                    )

                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement delaunay_triangulation(in_place=True) yet"
                    )

                if direction is not None:
                    import warnings

                    warnings.warn(
                        "the direction parameter of delaunay_triangulation() has been removed since it did not work correctly in previous versions of sage-flatsurf"
                    )

                if relabel is not None:
                    if relabel:
                        raise NotImplementedError(
                            "the relabel keyword has been removed from delaunay_triangulation(); use relabel() to use integer labels instead"
                        )
                    else:
                        import warnings

                        warnings.warn(
                            "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_triangulation()"
                        )

                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement delaunay_triangulation(in_place=True) yet"
                    )

                if relabel is not None:
                    if relabel:
                        raise NotImplementedError(
                            "the relabel keyword has been removed from delaunay_triangulation(); use relabel() to use integer labels instead"
                        )
                    else:
                        import warnings

                        warnings.warn(
                            "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_triangulation()"
                        )

                return self.delaunay_triangulate().codomain()

            def delaunay_triangulate(self):
                r"""
                Return a Delaunay triangulation of this surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces

                    sage: M = matrix([[2, 1], [1, 1]])
                    sage: S = M * translation_surfaces.infinite_staircase()
                    sage: S = S.delaunay_triangulate().codomain()

                    sage: S.root()
                    (0, 0)
                    sage: S.polygon((0, 0))
                    Polygon(vertices=[(0, 0), (-1, 0), (-1, -1)])

                    sage: S.is_delaunay_triangulated()
                    True

                TESTS::

                    sage: TestSuite(S).run()

                """
                triangulation = self.triangulate()

                from flatsurf.geometry.lazy import (
                    LazyDelaunayTriangulatedSurface,
                )

                codomain = LazyDelaunayTriangulatedSurface(
                    triangulation.codomain(), category=self.category()
                )

                from flatsurf.geometry.morphism import DelaunayTriangulationMorphism

                delaunay = DelaunayTriangulationMorphism._create_morphism(
                    triangulation.codomain(), codomain
                )

                from flatsurf.geometry.morphism import NamedFactorizationMorphism

                return NamedFactorizationMorphism._create_morphism(
                    self, codomain, "Delaunay triangulation", delaunay * triangulation
                )

            def delaunay_decomposition(
                self,
                triangulated=None,
                delaunay_triangulated=None,
                in_place=False,
                direction=None,
                relabel=None,
            ):
                r"""
                Return the Delaunay Decomposition of this surface.

                INPUT:

                - ``triangulated`` - ``None``; ignored

                - ``delaunay_triangulated`` - ``None``; ignored

                - ``in_place`` - boolean (default: ``False``); whether this
                  surface is modified and returned or just a copy of the
                  surface

                - ``direction`` - ``None``; ignored

                - ``relabel`` - ``None``; ignored

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces, Polygon, similarity_surfaces
                    sage: s0 = translation_surfaces.octagon_and_squares()
                    sage: a = s0.base_ring().gens()[0]
                    sage: m = Matrix([[1,2+a],[0,1]])
                    sage: s = m*s0
                    sage: s = s.triangulate().codomain()
                    sage: ss = s.delaunay_decomposition()
                    doctest:warning
                    ...
                    UserWarning: delaunay_decomposition() has been deprecated and will be removed in a future version of sage-flatsurf; use delaunay_decompose().codomain() instead

                    sage: len(ss.polygons())
                    3

                    sage: p = Polygon(edges=[(4,0),(-2,1),(-2,-1)])
                    sage: s0 = similarity_surfaces.self_glued_polygon(p)
                    sage: s = s0.delaunay_decomposition()
                    sage: TestSuite(s).run()

                    sage: m = matrix([[2,1],[1,1]])
                    sage: s = m*translation_surfaces.infinite_staircase()
                    sage: ss = s.delaunay_decomposition()
                    doctest:warning
                    ...
                    UserWarning: delaunay_decomposition() has been deprecated and will be removed in a future version of sage-flatsurf; use delaunay_decompose().codomain() instead
                    sage: ss.root()
                    (0, 0)
                    sage: ss.polygon(ss.root())
                    Polygon(vertices=[(0, 0), (-1, 0), (-1, -1), (0, -1)])
                    sage: TestSuite(ss).run()
                    sage: ss.is_delaunay_decomposed()
                    True

                """
                import warnings

                warnings.warn(
                    "delaunay_decomposition() has been deprecated and will be removed in a future version of sage-flatsurf; use delaunay_decompose().codomain() instead"
                )

                if triangulated is not None:
                    import warnings

                    warnings.warn(
                        "the triangulated keyword has been deprecated as a parameter for delaunay_decomposition() and will be removed in a future version of sage-flatsurf"
                    )

                if delaunay_triangulated is not None:
                    import warnings

                    warnings.warn(
                        "the delaunay_triangulated keyword has been deprecated as a parameter for delaunay_decomposition() and will be removed in a future version of sage-flatsurf"
                    )

                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement delaunay_decomposition(in_place=True) yet"
                    )

                if direction is not None:
                    import warnings

                    warnings.warn(
                        "the direction parameter of delaunay_decomposition() has been removed since it did not work correctly in previous versions of sage-flatsurf"
                    )

                if relabel is not None:
                    if relabel:
                        raise NotImplementedError(
                            "the relabel keyword has been removed from delaunay_decomposition(); use relabel() to use integer labels instead"
                        )
                    else:
                        import warnings

                        warnings.warn(
                            "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_decomposition()"
                        )

                return self.delaunay_decompose().codomain()

            def delaunay_decompose(self, codomain=None):
                r"""
                Return a Delaunay decomposition of this surface, i.e., a
                representation of this surface such that the circumscribed disk
                of each polygon contains exactly the vertices of that polygon.

                ALGORITHM:

                The Delaunay decomposition is obtained by removing ambiguous
                edges from a Delaunay triangulation. An edge is ambiguous if
                after its removal the circumscribed disk of each polygon
                (still) contains no vertices in its interior.

                INPUT:

                - ``codomain`` -- a Delaunay decomposed surface or ``None``
                  (default: ``None``); if present, the morphism returned goes
                  from this surface to ``codomain``.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.octagon_and_squares()
                    sage: a = S.base_ring().gens()[0]
                    sage: M = Matrix([[1, 2 + a], [0, 1]])
                    sage: S = M * S
                    sage: S = S.triangulate().codomain()

                    sage: S = S.delaunay_decompose().codomain()
                    sage: S.is_delaunay_decomposed()
                    True
                    sage: S
                    Delaunay cell decomposition of Triangulation of Translation Surface in H_3(4) built from 2 quadrilaterals and an octagon

                ::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(edges=[(4, 0), (-2, 1), (-2, -1)])
                    sage: T = similarity_surfaces.self_glued_polygon(P)

                    sage: T = T.delaunay_decompose().codomain()
                    sage: T.is_delaunay_decomposed()
                    True
                    sage: T
                    Delaunay cell decomposition of Half-Translation Surface in Q_0(-1^4) built from an isosceles triangle


                ::

                    sage: M = matrix([[2, 1], [1, 1]])
                    sage: U = M * translation_surfaces.infinite_staircase()

                    sage: U = U.delaunay_decompose().codomain()
                    sage: U.is_delaunay_decomposed()
                    True

                    sage: U.root()
                    (0, 0)
                    sage: U.polygon((0, 0))
                    Polygon(vertices=[(0, 0), (-1, 0), (-1, -1), (0, -1)])

                TESTS::

                    sage: TestSuite(S).run()
                    sage: TestSuite(T).run()
                    sage: TestSuite(U).run()

                """
                if codomain is not None:
                    from sage.all import End

                    decomposition = (
                        End(self).identity()
                        if self.is_delaunay_decomposed()
                        else self.delaunay_decompose()
                    )

                    from flatsurf.geometry.morphism import (
                        DelaunayDecompositionIsomorphism,
                    )

                    isomorphism = DelaunayDecompositionIsomorphism._create_morphism(
                        decomposition.codomain(), codomain
                    )

                    # This will raise an error if no such isomorphism exists.
                    isomorphism._factorization()

                    from flatsurf.geometry.morphism import NamedFactorizationMorphism

                    return NamedFactorizationMorphism._create_morphism(
                        self,
                        codomain,
                        "Delaunay decomposition",
                        isomorphism * decomposition,
                    )

                delaunay_triangulation = self.delaunay_triangulate()

                from flatsurf.geometry.lazy import LazyDelaunaySurface

                decomposed = LazyDelaunaySurface(
                    delaunay_triangulation.codomain(), category=self.category()
                )

                from flatsurf.geometry.morphism import DelaunayDecompositionMorphism

                decomposition = DelaunayDecompositionMorphism._create_morphism(
                    delaunay_triangulation.codomain(), decomposed
                )

                from flatsurf.geometry.morphism import NamedFactorizationMorphism

                return NamedFactorizationMorphism._create_morphism(
                    self,
                    decomposed,
                    "Delaunay decomposition",
                    decomposition * delaunay_triangulation,
                )

            def saddle_connections(
                self,
                squared_length_bound,
                initial_label=None,
                initial_vertex=None,
                sc_list=None,
                check=False,
            ):
                r"""
                Returns a list of saddle connections on the surface whose length squared is less than or equal to squared_length_bound.
                The length of a saddle connection is measured using holonomy from polygon in which the trajectory starts.

                If initial_label and initial_vertex are not provided, we return all saddle connections satisfying the bound condition.

                If initial_label and initial_vertex are provided, it only provides saddle connections emanating from the corresponding
                vertex of a polygon. If only initial_label is provided, the added saddle connections will only emanate from the
                corresponding polygon.

                If sc_list is provided the found saddle connections are appended to this list and the resulting list is returned.

                If check==True it uses the checks in the SaddleConnection class to sanity check our results.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: s = translation_surfaces.square_torus()
                    sage: sc_list = s.saddle_connections(13, check=True)
                    sage: len(sc_list)
                    32

                TESTS:

                Verify that #255 has been resolved::

                    sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon

                    sage: P = Polygon(vertices=[(0, 0), (64, 0), (1, 1)])
                    sage: Q = Polygon(vertices=[(0, 0), (64, 0), (32, 1)])

                    sage: S = MutableOrientedSimilaritySurface(QQ)

                    sage: S.add_polygon(P)
                    0
                    sage: S.add_polygon(Q)
                    1

                    sage: S.glue((1, 2), (1, 2))
                    sage: S.glue((1, 1), (1, 1))
                    sage: S.glue((0, 2), (0, 2))
                    sage: S.glue((0, 0), (0, 0))
                    sage: S.glue((0, 1), (1, 0))

                    sage: S.set_immutable()

                    sage: len(S.saddle_connections(63**2, initial_label=0, initial_vertex=0))
                    19

                """
                if squared_length_bound <= 0:
                    raise ValueError

                if sc_list is None:
                    sc_list = []
                if initial_label is None:
                    if not self.is_finite_type():
                        raise NotImplementedError
                    if initial_vertex is not None:
                        raise ValueError(
                            "when initial_label is not provided, then initial_vertex must not be provided either"
                        )
                    for label in self.labels():
                        self.saddle_connections(
                            squared_length_bound, initial_label=label, sc_list=sc_list
                        )
                    return sc_list
                if initial_vertex is None:
                    for vertex in range(len(self.polygon(initial_label).vertices())):
                        self.saddle_connections(
                            squared_length_bound,
                            initial_label=initial_label,
                            initial_vertex=vertex,
                            sc_list=sc_list,
                        )
                    return sc_list

                # Now we have a specified initial_label and initial_vertex
                from flatsurf.geometry.similarity import SimilarityGroup

                SG = SimilarityGroup(self.base_ring())
                start_data = (initial_label, initial_vertex)
                from flatsurf.geometry.circle import Circle

                circle = Circle(
                    (0, 0),
                    squared_length_bound,
                    base_ring=self.base_ring(),
                )
                p = self.polygon(initial_label)
                v = p.vertex(initial_vertex)
                last_sim = SG(-v[0], -v[1])

                # First check the edge eminating rightward from the start_vertex.
                e = p.edge(initial_vertex)
                if e[0] ** 2 + e[1] ** 2 <= squared_length_bound:
                    from flatsurf.geometry.surface_objects import SaddleConnection

                    sc_list.append(SaddleConnection(self, start_data, e))

                # Represents the bounds of the beam of trajectories we are sending out.
                wedge = (
                    last_sim(p.vertex((initial_vertex + 1) % len(p.vertices()))),
                    last_sim(
                        p.vertex(
                            (initial_vertex + len(p.vertices()) - 1) % len(p.vertices())
                        )
                    ),
                )

                # This will collect the data we need for a depth first search.
                chain = [
                    (
                        last_sim,
                        initial_label,
                        wedge,
                        [
                            (initial_vertex + len(p.vertices()) - i) % len(p.vertices())
                            for i in range(2, len(p.vertices()))
                        ],
                    )
                ]

                while len(chain) > 0:
                    # Should verts really be edges?
                    sim, label, wedge, verts = chain[-1]
                    if len(verts) == 0:
                        chain.pop()
                        continue
                    vert = verts.pop()
                    p = self.polygon(label)
                    # First check the vertex
                    vert_position = sim(p.vertex(vert))
                    from flatsurf.geometry.euclidean import ccw

                    if (
                        ccw(wedge[0], vert_position) > 0
                        and ccw(vert_position, wedge[1]) > 0
                        and vert_position[0] ** 2 + vert_position[1] ** 2
                        <= squared_length_bound
                    ):
                        from flatsurf.geometry.surface_objects import SaddleConnection

                        sc_list.append(
                            SaddleConnection(
                                self,
                                start_data,
                                vert_position,
                                end_data=(label, vert),
                                end_direction=~sim.derivative() * -vert_position,
                                holonomy=vert_position,
                                end_holonomy=~sim.derivative() * -vert_position,
                                check=check,
                            )
                        )
                    # Now check if we should develop across the edge
                    vert_position2 = sim(p.vertex((vert + 1) % len(p.vertices())))
                    if (
                        ccw(vert_position, vert_position2) > 0
                        and ccw(wedge[0], vert_position2) > 0
                        and ccw(vert_position, wedge[1]) > 0
                        and circle.line_segment_position(vert_position, vert_position2)
                        == 1
                    ):
                        if ccw(wedge[0], vert_position) > 0:
                            # First in new_wedge should be vert_position
                            if ccw(vert_position2, wedge[1]) > 0:
                                new_wedge = (vert_position, vert_position2)
                            else:
                                new_wedge = (vert_position, wedge[1])
                        else:
                            if ccw(vert_position2, wedge[1]) > 0:
                                new_wedge = (wedge[0], vert_position2)
                            else:
                                new_wedge = wedge
                        new_label, new_edge = self.opposite_edge(label, vert)
                        new_sim = sim * ~self.edge_transformation(label, vert)
                        p = self.polygon(new_label)
                        chain.append(
                            (
                                new_sim,
                                new_label,
                                new_wedge,
                                [
                                    (new_edge + len(p.vertices()) - i)
                                    % len(p.vertices())
                                    for i in range(1, len(p.vertices()))
                                ],
                            )
                        )
                return sc_list

            def ramified_cover(self, d, data):
                r"""
                Return a ramified cover of this surface.

                INPUT:

                - ``d`` - integer (the degree of the cover)

                - ``data`` - a dictionary which to a pair ``(label, edge_num)`` associates a permutation
                  of {1,...,d}

                EXAMPLES:

                The L-shape origami::

                    sage: import flatsurf
                    sage: T = flatsurf.translation_surfaces.square_torus()
                    sage: T.ramified_cover(3, {(0,0): '(1,2)', (0,1): '(1,3)'})
                    Translation Surface in H_2(2) built from 3 squares
                    sage: O = T.ramified_cover(3, {(0,0): '(1,2)', (0,1): '(1,3)'})
                    sage: O.stratum()
                    H_2(2)

                TESTS::

                    sage: import flatsurf
                    sage: T = flatsurf.translation_surfaces.square_torus()
                    sage: T.ramified_cover(3, {(0,0): '(1,2)', (0,2): '(1,3)'})
                    Traceback (most recent call last):
                    ...
                    ValueError: inconsistent covering data

                """
                from sage.groups.perm_gps.permgroup_named import SymmetricGroup

                G = SymmetricGroup(d)
                for k in data:
                    data[k] = G(data[k])
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                cover = MutableOrientedSimilaritySurface(self.base_ring())
                edges = set(self.edges())
                cover_labels = {}
                for i in range(1, d + 1):
                    for lab in self.labels():
                        cover_labels[(lab, i)] = cover.add_polygon(self.polygon(lab))
                while edges:
                    lab, e = elab = edges.pop()
                    llab, ee = eelab = self.opposite_edge(lab, e)
                    edges.remove(eelab)
                    if elab in data:
                        if eelab in data:
                            if not (data[elab] * data[eelab]).is_one():
                                raise ValueError("inconsistent covering data")
                        s = data[elab]
                    elif eelab in data:
                        s = ~data[eelab]
                    else:
                        s = G.one()

                    for i in range(1, d + 1):
                        p0 = cover_labels[(lab, i)]
                        p1 = cover_labels[(lab, s(i))]
                        cover.glue((p0, e), (p1, ee))
                cover.set_immutable()
                return cover

            def subdivide(self):
                r"""
                Return a copy of this surface whose polygons have been partitioned into
                smaller triangles with
                :meth:`~.euclidean_polygons.EuclideanPolygons.Simple.Convex.ParentMethods.subdivide`.

                EXAMPLES:

                A surface consisting of a single triangle::

                    sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon

                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(Polygon(edges=[(1, 0), (0, 1), (-1, -1)]), label="Δ")
                    'Δ'

                Subdivision of this surface yields a surface with three triangles::

                    sage: T = S.subdivide()
                    sage: T.labels()
                    (('Δ', 0), ('Δ', 1), ('Δ', 2))

                Note that the new labels are old labels plus an index. We verify that
                the triangles are glued correctly::

                    sage: list(T.gluings())
                    [((('Δ', 0), 1), (('Δ', 1), 2)),
                     ((('Δ', 0), 2), (('Δ', 2), 1)),
                     ((('Δ', 1), 1), (('Δ', 2), 2)),
                     ((('Δ', 1), 2), (('Δ', 0), 1)),
                     ((('Δ', 2), 1), (('Δ', 0), 2)),
                     ((('Δ', 2), 2), (('Δ', 1), 1))]

                If we add another polygon to the original surface and glue things, we
                can see how existing gluings are preserved when subdividing::

                    sage: S.add_polygon(Polygon(edges=[(1, 0), (0, 1), (-1, 0), (0, -1)]), label='□')
                    '□'

                    sage: S.glue(("Δ", 0), ("□", 2))
                    sage: S.glue(("□", 1), ("□", 3))

                    sage: T = S.subdivide()

                    sage: T.labels()
                    (('Δ', 0), ('□', 2), ('Δ', 1), ('Δ', 2), ('□', 3), ('□', 1), ('□', 0))
                    sage: list(sorted(T.gluings()))
                    [((('Δ', 0), 0), (('□', 2), 0)),
                     ((('Δ', 0), 1), (('Δ', 1), 2)),
                     ((('Δ', 0), 2), (('Δ', 2), 1)),
                     ((('Δ', 1), 1), (('Δ', 2), 2)),
                     ((('Δ', 1), 2), (('Δ', 0), 1)),
                     ((('Δ', 2), 1), (('Δ', 0), 2)),
                     ((('Δ', 2), 2), (('Δ', 1), 1)),
                     ((('□', 0), 1), (('□', 1), 2)),
                     ((('□', 0), 2), (('□', 3), 1)),
                     ((('□', 1), 0), (('□', 3), 0)),
                     ((('□', 1), 1), (('□', 2), 2)),
                     ((('□', 1), 2), (('□', 0), 1)),
                     ((('□', 2), 0), (('Δ', 0), 0)),
                     ((('□', 2), 1), (('□', 3), 2)),
                     ((('□', 2), 2), (('□', 1), 1)),
                     ((('□', 3), 0), (('□', 1), 0)),
                     ((('□', 3), 1), (('□', 0), 2)),
                     ((('□', 3), 2), (('□', 2), 1))]

                """
                labels = list(self.labels())
                polygons = [self.polygon(label) for label in labels]

                subdivisions = [p.subdivide() for p in polygons]

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                surface = MutableOrientedSimilaritySurface(self.base())

                # Add subdivided polygons
                for s, subdivision in enumerate(subdivisions):
                    label = labels[s]
                    for p, polygon in enumerate(subdivision):
                        surface.add_polygon(polygon, label=(label, p))

                surface.set_roots((label, 0) for label in self.roots())

                # Add gluings between subdivided polygons
                for s, subdivision in enumerate(subdivisions):
                    label = labels[s]
                    for p in range(len(subdivision)):
                        surface.glue(
                            ((label, p), 1), ((label, (p + 1) % len(subdivision)), 2)
                        )

                        # Add gluing from original surface
                        opposite = self.opposite_edge(label, p)
                        if opposite is not None:
                            surface.glue(((label, p), 0), (opposite, 0))

                surface.set_immutable()

                return surface

            def subdivide_edges(self, parts=2):
                r"""
                Return a copy of this surface whose edges have been split into
                ``parts`` equal pieces each.

                INPUT:

                - ``parts`` -- a positive integer (default: 2)

                EXAMPLES:

                A surface consisting of a single triangle::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: from flatsurf.geometry.polygon import Polygon

                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(Polygon(edges=[(1, 0), (0, 1), (-1, -1)]), label="Δ")
                    'Δ'

                Subdividing this triangle yields a triangle with marked points along
                the edges::

                    sage: T = S.subdivide_edges()

                If we add another polygon to the original surface and glue them, we
                can see how existing gluings are preserved when subdividing::

                    sage: S.add_polygon(Polygon(edges=[(1, 0), (0, 1), (-1, 0), (0, -1)]), label='□')
                    '□'

                    sage: S.glue(("Δ", 0), ("□", 2))
                    sage: S.glue(("□", 1), ("□", 3))

                    sage: T = S.subdivide_edges().codomain()
                    sage: list(sorted(T.gluings()))
                    [(('Δ', 0), ('□', 5)),
                     (('Δ', 1), ('□', 4)),
                     (('□', 2), ('□', 7)),
                     (('□', 3), ('□', 6)),
                     (('□', 4), ('Δ', 1)),
                     (('□', 5), ('Δ', 0)),
                     (('□', 6), ('□', 3)),
                     (('□', 7), ('□', 2))]

                """
                labels = list(self.labels())
                polygons = [self.polygon(label) for label in labels]

                subdivideds = [p.subdivide_edges(parts=parts) for p in polygons]

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                surface = MutableOrientedSimilaritySurface(self.base())

                # Add subdivided polygons
                for s, subdivided in enumerate(subdivideds):
                    surface.add_polygon(subdivided, label=labels[s])

                surface.set_roots(self.roots())

                # Reestablish gluings between polygons
                for label, polygon, subdivided in zip(labels, polygons, subdivideds):
                    for e in range(len(polygon.vertices())):
                        opposite = self.opposite_edge(label, e)
                        if opposite is not None:
                            for p in range(parts):
                                surface.glue(
                                    (label, e * parts + p),
                                    (
                                        opposite[0],
                                        opposite[1] * parts + (parts - p - 1),
                                    ),
                                )

                surface.set_immutable()

                from flatsurf.geometry.morphism import SubdivideEdgesMorphism

                return SubdivideEdgesMorphism._create_morphism(self, surface, parts)

    class Rational(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by similarity surfaces where all similarities that
        describe how edges are glued only use rational rotations, i.e.,
        rotations by a rational multiple of π.

        Note that this differs slightly from the usual definition of
        "rational". Normally, a surface would be rational if it can be
        described using only such similarities. Here we require that the
        similarities used are actually of that kind.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: "Rational" in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all surfaces built from Euclidean
            polygons glued by similarities that have rational monodromy, i.e.,
            `monodromy
            <https://en.wikipedia.org/wiki/(G,X)-manifold#Monodromy>`_ gives
            similarities whose rotational part has finite order.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            @staticmethod
            def _is_rational_surface(surface, limit=None):
                r"""
                Return whether the gluings of this surface lead to a rational
                surface, i.e., whether all gluings use similarities whose
                rotational part uses only a rational multiple of π as a
                rotation.

                This is a helper method for
                :meth:`flatsurf.geometry.categories.similarity_surfaces.ParentMethods.is_rational_surface`.

                INPUT:

                - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
                  the first ``limit`` polygons are checked

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.infinite_staircase()

                    sage: from flatsurf.geometry.categories import SimilaritySurfaces
                    sage: SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(S, limit=8)
                    doctest:warning
                    ...
                    UserWarning: limit has been deprecated as a keyword argument for _is_rational_surface() and will be removed from a future version of sage-flatsurf; ...
                    True

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(MutableOrientedSimilaritySurface.from_surface(S, labels=S.labels()[:8]))
                    True

                """
                if limit is not None:
                    import warnings

                    warnings.warn(
                        "limit has been deprecated as a keyword argument for _is_rational_surface() and will be removed from a future version of sage-flatsurf; "
                        "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                    )

                if "Oriented" not in surface.category().axioms():
                    raise NotImplementedError

                labels = surface.labels()

                if limit is not None:
                    from itertools import islice

                    labels = islice(labels, limit)

                checked = set()

                for label in labels:
                    for edge in range(len(surface.polygon(label).vertices())):
                        cross = surface.opposite_edge(label, edge)

                        if cross is None:
                            continue

                        if cross in checked:
                            continue

                        checked.add((label, edge))

                        # We do not call self.edge_matrix() since the surface might
                        # have overridden this (just returning the identity matrix e.g.)
                        # and we want to deduce the matrix from the attached polygon
                        # edges instead.
                        matrix = SimilaritySurfaces.Oriented.ParentMethods.edge_matrix.f(  # pylint: disable=no-member
                            surface, label, edge
                        )

                        if matrix.is_diagonal():
                            continue

                        a = AA(matrix[0, 0])
                        b = AA(matrix[1, 0])
                        q = (a**2 + b**2).sqrt()

                        from flatsurf.geometry.euclidean import (
                            is_cosine_sine_of_rational,
                        )

                        if not is_cosine_sine_of_rational(a / q, b / q):
                            return False

                return True

            def is_rational_surface(self):
                r"""
                Return whether all edges of this surface are glued with
                similarities whose rotational part is by a rational multiple of
                π, i.e., return ``True`` since this is a rational surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.infinite_staircase()
                    sage: S.is_rational_surface()
                    True

                """
                return True

            def _test_rational_surface(self, **options):
                r"""
                Verify that this is a rational similarity surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S._test_rational_surface()

                """
                tester = self._tester(**options)

                limit = None

                if not self.is_finite_type():
                    from flatsurf import MutableOrientedSimilaritySurface

                    self = MutableOrientedSimilaritySurface.from_surface(
                        self, labels=self.labels()[:32]
                    )

                tester.assertTrue(
                    SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(
                        self, limit=limit
                    )
                )

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The category of surfaces built by gluing a finite number of Euclidean
        polygons with similarities.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
            sage: S = similarity_surfaces.self_glued_polygon(P)
            sage: "FiniteType" in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all surfaces that are built from
            finitely many polygons in the real plane glued with similarities.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def num_singularities(self):
                r"""
                EXAMPLES::

                    sage: from flatsurf import translation_surfaces

                    sage: translation_surfaces.regular_octagon().num_singularities()
                    doctest:warning
                    ...
                    UserWarning: num_singularities() has been deprecated and will be removed in a future version of sage-flatsurf; use len(vertices()) instead
                    1

                    sage: S = SymmetricGroup(4)
                    sage: r = S('(1,2)(3,4)')
                    sage: u = S('(2,3)')
                    sage: translation_surfaces.origami(r,u).num_singularities()
                    2

                    sage: S = SymmetricGroup(8)
                    sage: r = S('(1,2,3,4,5,6,7,8)')
                    sage: u = S('(1,8,5,4)(2,3)(6,7)')
                    sage: translation_surfaces.origami(r,u).num_singularities()
                    4
                """
                import warnings

                warnings.warn(
                    "num_singularities() has been deprecated and will be removed in a future version of sage-flatsurf; use len(vertices()) instead"
                )

                return len(self.vertices())

            def _test_eq_surface(self, **options):
                r"""
                Verify that this surface follows our standards for equality of
                surfaces.

                We want two surfaces to compare equal (`S == T`) iff they are
                virtually indistinguishable; so without a lot of non-Pythonic
                effort, you should not be able to tell them apart. They have
                (virtually) the same type, are made from equally labeled
                polygons with indistinguishable coordinates and equal gluings.
                Any other data that was used when creating them should be
                indistinguishable. They might of course live at different
                memory addresses have differences in their internal caches and
                representation but everything user-facing should be the same.

                People often want `==` to mean that the two surfaces are
                isomorphic in some more-or-less strong sense. Such a notion for
                `==` always leads to trouble down the road. The operator `==`
                is used to identify surfaces in caches and identify surfaces in
                sets. Sometimes "are isomorphic" is a good notion in such cases
                but most of the time "are indistinguishable" is the much safer
                default. Also, "are isomorphic" is often costly or, e.g. in the case
                of infinite surfaces, not even decidable.

                Currently, we do treat two surfaces as equal even if they
                differ by category because categories can presently be changed
                for immutable surfaces (as more properties of the surface are
                found.)

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S._test_eq_surface()

                :meta public:

                """
                tester = self._tester(**options)

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                copy = MutableOrientedSimilaritySurface.from_surface(self)
                if not self.is_mutable():
                    copy.set_immutable()

                if isinstance(self, MutableOrientedSimilaritySurface):
                    tester.assertEqual(self, copy)
                    tester.assertFalse(self != copy)
                else:
                    tester.assertNotEqual(self, copy)
                    tester.assertTrue(self != copy)

        class Oriented(SurfaceCategoryWithAxiom):
            r"""
            The category of surfaces built from finitely many Euclidean
            polygons glued with singularities with an orientation that is
            compatible with the embedding that the polygons inherit from the
            real plane.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: "Oriented" in S.category().axioms()
                True

            """

            class ParentMethods:
                r"""
                Provides methods available to all surfaces that are built from
                finitely many Euclidean polygons that are glued by similarities
                and are oriented with the natural orientation of the polygons
                in the real plane.

                If you want to add functionality for such surfaces you most likely
                want to put it here.
                """

                def reposition_polygons(self, in_place=False, relabel=None):
                    r"""
                    We choose a maximal tree in the dual graph of the decomposition into
                    polygons, and ensure that the gluings between two polygons joined by
                    an edge in this tree is by translation.

                    This guarantees that the group generated by the edge identifications is
                    minimal among representations of the surface. In particular, if for instance
                    you have a translation surface which is anot representable as a translation
                    surface (because polygons are presented with rotations) then after this
                    change it will be representable as a translation surface.
                    """
                    if in_place:
                        raise NotImplementedError(
                            "this surface does not implement reposition_polygons(in_place=True) yet"
                        )

                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    s = MutableOrientedSimilaritySurface.from_surface(self)
                    s.reposition_polygons(in_place=True, relabel=relabel)
                    s.set_immutable()
                    return s

                def standardize_polygons(self, in_place=False):
                    r"""
                    Return a surface with each polygon replaced with a new
                    polygon which differs by translation and reindexing. The
                    new polygon will have the property that vertex zero is the
                    origin, and all vertices lie either in the upper half
                    plane, or on the x-axis with non-negative x-coordinate.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: s=translation_surfaces.veech_double_n_gon(4)
                        sage: s.polygon(1)
                        Polygon(vertices=[(0, 0), (-1, 0), (-1, -1), (0, -1)])
                        sage: [s.opposite_edge(0,i) for i in range(4)]
                        [(1, 0), (1, 1), (1, 2), (1, 3)]
                        sage: ss=s.standardize_polygons()
                        sage: ss.polygon(1)
                        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                        sage: [ss.opposite_edge(0,i) for i in range(4)]
                        [(1, 2), (1, 3), (1, 0), (1, 1)]
                        sage: TestSuite(ss).run()

                    """
                    if in_place:
                        raise NotImplementedError(
                            "cannot standardize polygons in_place anymore on this surface; use in_place=False to create a copy of the surface with standardized polygons"
                        )

                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    S = MutableOrientedSimilaritySurface.from_surface(
                        self, category=self.category()
                    )
                    S.standardize_polygons(in_place=True)
                    S.set_immutable()
                    return S

                def fundamental_group(self, base_label=None):
                    r"""
                    Return the fundamental group of this surface.
                    """
                    if base_label is None:
                        base_label = self.root()

                    from flatsurf.geometry.fundamental_group import FundamentalGroup

                    return FundamentalGroup(self, base_label)

    class SubcategoryMethods:
        def Rational(self):
            r"""
            Return the subcategory of surfaces with rational monodromy, see
            :class:`~.SimilaritySurfaces.Rational`.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import SimilaritySurfaces
                sage: C = SimilaritySurfaces()
                sage: C.Rational()
                Category of rational similarity surfaces

            """
            return self._with_axiom("Rational")


# Currently, there is no "Rational" axiom in SageMath so we make it known to
# the category framework.
all_axioms += ("Rational",)
