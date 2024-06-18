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
from sage.categories.category_with_axiom import all_axioms
from sage.misc.cachefunc import cached_method, cached_in_parent_method
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
            specific to the smaller classes of surfaces.)


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
            followed by a dilation (multiplication by a diagonal matrix.)

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
              half-translation, i.e., a translation followed by a reflection in
              a point (equivalently, a rotation by π.)

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
            Apply the 2×2 ``matrix`` to the polygons of this surface.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: s = translation_surfaces.infinite_staircase()
                sage: s
                The infinite staircase
                sage: m = matrix([[1,2],[0,1]])
                sage: s2 = m * s
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

            return self.apply_matrix(matrix).codomain()

        def harmonic_differentials(
            self, error, cell_decomposition, check=True, category=None
        ):
            if self.is_mutable():
                raise ValueError(
                    "surface must be immutable to compute harmonic differentials"
                )

            from sage.all import RR

            coefficients = RR

            if category is None:
                from sage.categories.all import Modules

                category = Modules(coefficients)

            return self._harmonic_differentials(
                error=error,
                cell_decomposition=cell_decomposition,
                check=check,
                category=category,
            )

        @cached_method
        def _harmonic_differentials(self, error, cell_decomposition, check, category):
            from flatsurf.geometry.harmonic_differentials import (
                HarmonicDifferentialSpace,
            )

            return HarmonicDifferentialSpace(
                self, error, cell_decomposition, check, category
            )

        def homology(
            self,
            k=1,
            coefficients=None,
            generators="edge",
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

            - ``generators`` -- a string (default: ``"edge"``); how the
              generators of homology are represented. Currently only ``"edge"``
              is implemented, i.e., the generators are written as formal sums
              of half edges.

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
                generators=generators,
                relative=relative,
                implementation=implementation,
                category=category,
            )

        @cached_method
        def _homology(
            self, k, coefficients, generators, relative, implementation, category
        ):
            r"""
            Return the ``k``-th homology group of this surface.

            This is a helper method for :meth:`homology`. We cannot make
            :class:`SimplicialHomologyGroup` a unique representation because
            equal surfaces can be non-identical so the homology of
            non-identical surfaces could be identical. We work around this
            issue by attaching the homology to the actual surface so a surface
            has a unique homology, but it is different from another equal
            surface's.

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
                self, k, coefficients, generators, relative, implementation, category
            )

        def cohomology(
            self, k=1, coefficients=None, implementation="dual", category=None
        ):
            r"""
            Return the ``k``-th simplicial cohomology group of this surface.

            INPUT:

            - ``k`` -- an integer (default: ``1``)

            - ``coefficients`` -- a ring (default: the reals); consider
              cohomology with coefficients in this ring

            - ``implementation`` -- a string (default: ``"dual"``); the
              algorithm used to compute the cohomology groups. Currently only
              ``"dual"`` is supported, i.e., the groups are computed as duals
              of the generic homology groups from SageMath.

            - ``category`` -- a category; if not specified, a category for the
              homology group is chosen automatically depending on
              ``coefficients``.

            EXAMPLES::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: S.cohomology()
                H¹(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon; Real Field with 53 bits of precision)

            ::

                sage: S.cohomology(0)
                H⁰(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon; Real Field with 53 bits of precision)

            """
            if self.is_mutable():
                raise ValueError("surface must be immutable to compute cohomology")

            from sage.all import ZZ

            k = ZZ(k)

            from sage.all import RR

            coefficients = coefficients or RR

            if category is None:
                from sage.categories.all import Modules

                category = Modules(coefficients)

            return self._cohomology(
                k=k,
                coefficients=coefficients,
                implementation=implementation,
                category=category,
            )

        @cached_method
        def _cohomology(self, k, coefficients, implementation, category):
            r"""
            Return the ``k``-th cohomology group of this surface.

            This is a helper method for :meth:`cohomology`. We cannot make
            :class:`SimplicialCohomologyGroup` a unique representation because
            equal surfaces can be non-identical so the cohomology of
            non-identical surfaces could be identical. We work around this
            issue by attaching the cohomology to the actual surface so a
            surface has a unique cohomology, but it is different from another
            equal surface's.

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
                self, k, coefficients, implementation, category
            )

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
                    "The defaults for apply_matrix() are going to change in a future version of sage-flatsurf; previously, apply_matrix() was performed in_place=True. In a future version of sage-flatsurf the default is going to change to in_place=False. In the meantime, please pass in_place=True/False explicitly."
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
                    sage: T.glue((0, 0), (0, 2))
                    sage: T(0, 0).angle()
                    3/8

                .. SEEALSO::

                    :meth:`turns` for the number of full 2π turns at a point.

                """
                if self.is_vertex():
                    angle = self._angle_vertex(numerical=numerical)

                    if self.parent().is_mutable():
                        self._angle_vertex.clear_cache()  # pylint: disable=no-member

                    return angle

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

            @cached_in_parent_method
            def _angle_vertex(self, numerical):
                r"""
                Return the total angle at this vertex in multiples of 2π.

                INPUT:

                - ``numerical`` -- a boolean; whether to return the angle as a
                  floating point number.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S(0, 0)._angle_vertex(numerical=False)
                    1

                """
                turns, start, end = self.turns()

                from flatsurf.geometry.euclidean import angle

                return turns + angle(start, end, numerical=numerical)

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
                    # angle <2π.) So at least one π turn must have happened.
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

            @cached_method
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

            @cached_method
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
                Return a morphism to a surface whose polygons have been
                relabeled according to ``relabeling``.

                INPUT:

                - ``relabeling`` -- a dict or a callable mapping all labels of
                  this surface to new (unique) labels

                - ``in_place`` -- a boolean (default: ``False``); whether to
                  mutate this surface or return a morphism to an independent copy.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.veech_double_n_gon(5)
                    sage: relabeling = S.relabel({0: 1, 1: 2})
                    sage: SS = relabeling.codomain()
                    sage: SS
                    Translation Surface in H_2(2) built from 2 regular pentagons
                    sage: SS.root()
                    1
                    sage: SS.opposite_edge(1, 0)
                    (2, 0)

                    sage: TestSuite(SS).run()

                The relabeling can also be a callable::

                    sage: unrelabeling = SS.relabel(lambda label: label -1)
                    sage: SSS = unrelabeling.codomain()
                    sage: SSS == S
                    True

                """
                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement relabel(in_place=True) yet"
                    )

                from flatsurf.geometry.surface import (
                    MutableOrientedSimilaritySurface,
                )

                s = MutableOrientedSimilaritySurface.from_surface(self)
                morphism = s.relabel(relabeling=relabeling, in_place=True)
                s.set_immutable()

                return morphism.change(domain=self, codomain=s, check=False)

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
                    message += " However, there is no immediate replacement for lazy copying of infinite surfaces. Have a look at the implementation of flatsurf.geometry.lazy.LazyMutableSurface and adapt it to your needs."

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

                """
                from flatsurf.geometry.surface import BaseRingChangedSurface

                return BaseRingChangedSurface(self, ring)

            def triangle_flip(
                self, label, edge, in_place=False, test=None, direction=None
            ):
                r"""
                Flip the diagonal of the quadrilateral formed by two triangles
                glued together along the ``edge`` of the polygon with
                ``label``.

                ALGORITHM:

                This can be broken into two steps: join along the diagonal
                ``edge`` to form a convex quadilateral, then cut along the
                other diagonal.

                The labels of the new triangles will be the labels of the old
                triangles, assigned such that the flip turns the diagonal
                counterclockwise. To make this precise suppose that triangle
                with the ``label`` has the diagonal ``edge`` followed (in a
                counterclockwise walk of the triangle edges) by the two edges
                ``c`` and ``d``, and the triangle across the ``edge`` has the
                label ``opposite_label`` with the diagonal ``opposite_edge``
                followed by the edges ``a`` and ``b``. Then the new triangles
                are going to be labeled ``label`` and ``opposite_label`` and
                their diagonals are going to be ``edge`` and ``opposite_edge``
                respectively, and the choice is made such that the ``edge`` is
                followed by the edges that used to be called ``b`` and ``c``.
                Note that there is no choice here to be made if ``edge`` is
                glued to itself.

                INPUT:

                - ``label`` - label of a polygon

                - ``edge`` - an integer; the edge of a polygon joining two
                  triangles such that the quadrilateral formed after gluing
                  these triangles is strictly convex

                - ``in_place`` (boolean) - If True do the flip to the current surface
                  which must be mutable. In this case the updated surface will be
                  returned.  Otherwise a mutable copy is made and then an edge is
                  flipped, which is then returned.

                - ``test`` (boolean) - If True we don't actually flip, and we return
                  True or False depending on whether or not the flip would be
                  successful. (deprecated)

                - ``direction`` (ignored and deprecated)

                OUTPUT:

                A morphism from the original surface to the surface with the
                flipped edge.

                EXAMPLES::

                    sage: from flatsurf import similarity_surfaces, MutableOrientedSimilaritySurface, Polygon

                    sage: s = similarity_surfaces.right_angle_triangle(1, 1)
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
                    sage: s.triangle_flip(0, 0, test=True)
                    doctest:warning
                    ...
                    UserWarning: triangle_flip(test=True) has been deprecated and will be removed in a future version of sage-flatsurf. Use is_convex(strict=True) instead.
                    False
                    sage: s.triangle_flip(0, 1, test=True)
                    True
                    sage: s.triangle_flip(0, 2, test=True)
                    False

                    sage: s = similarity_surfaces.right_angle_triangle(1, 1)
                    sage: s = MutableOrientedSimilaritySurface.from_surface(s)
                    sage: s.triangle_flip(0, 0, in_place=True)
                    Traceback (most recent call last):
                    ...
                    ValueError: polygon has zero area

                    sage: s.triangle_flip(0, 1, in_place=True)
                    Triangle Flip morphism:
                      From: Unknown Surface
                      To:   Unknown Surface
                    sage: s
                    Rational Cone Surface built from 2 isosceles triangles
                    sage: s.polygon(0)
                    Polygon(vertices=[(0, 0), (0, 1), (-1, 0)])
                    sage: s.polygon(1)
                    Polygon(vertices=[(-1, 1), (-1, 0), (0, 1)])
                    sage: s.gluings()
                    (((0, 0), (0, 2)), ((0, 1), (1, 1)), ((0, 2), (0, 0)), ((1, 0), (1, 2)), ((1, 1), (0, 1)), ((1, 2), (1, 0)))
                    sage: s.triangle_flip(0, 2, in_place=True)
                    Traceback (most recent call last):
                    ...
                    ValueError: polygon has zero area

                    sage: p = Polygon(edges=[(2,0),(-1,3),(-1,-3)])
                    sage: s = similarity_surfaces.self_glued_polygon(p)
                    sage: s = MutableOrientedSimilaritySurface.from_surface(s)
                    sage: s.triangle_flip(0, 1, in_place=True)
                    Triangle Flip morphism:
                      From: Unknown Surface
                      To:   Unknown Surface
                    sage: s
                    Half-Translation Surface built from a triangle

                    sage: s.set_immutable()

                    sage: from flatsurf.geometry.categories import DilationSurfaces
                    sage: s in DilationSurfaces()
                    True
                    sage: s.labels()
                    (0,)
                    sage: s.polygons()
                    (Polygon(vertices=[(0, 0), (1, 3), (-2, 0)]),)
                    sage: s.gluings()
                    (((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)))
                    sage: TestSuite(s).run()

                """
                if test is not None:
                    if test:
                        import warnings

                        warnings.warn(
                            "triangle_flip(test=True) has been deprecated and will be removed in a future version of sage-flatsurf. Use is_convex(strict=True) instead."
                        )
                        return self.is_convex(label, edge, strict=True)
                    else:
                        import warnings

                        warnings.warn(
                            "the test keyword argument of triangle_flip() has been deprecated and will be removed in a future version of sage-flatsurf."
                        )

                if direction is not None:
                    direction = None

                    import warnings

                    warnings.warn(
                        "the direction keyword argument of triangle_flip() has been deprecated and will be removed in a future version of sage-flatsurf. The direction keyword argument is ignored. The flip is always performed such that the flipped edge rotates counterclockwise."
                    )

                if in_place:
                    raise NotImplementedError(
                        "this surface does not support triangle_flip(in_place=True) yet"
                    )

                from flatsurf.geometry.surface import (
                    MutableOrientedSimilaritySurface,
                )

                s = MutableOrientedSimilaritySurface.from_surface(self)
                morphism = s.triangle_flip(
                    label=label,
                    edge=edge,
                    in_place=True,
                    test=test,
                    direction=direction,
                )
                s.set_immutable()

                from flatsurf.geometry.morphism import TriangleFlipMorphism

                return TriangleFlipMorphism(self, s, morphism._edge_map)

            def is_convex(self, label, edge, strict=False):
                r"""
                Return whether the polygon obtained by joining the polygon with
                ``label`` and the one acress ``edge`` is convex.
                """
                opposite_label, opposite_edge = self.opposite_edge(label, edge)

                if opposite_label == label:
                    raise ValueError("edge must be joining two different polygons")

                p1 = self.polygon(label)

                from flatsurf import Polygon

                sim = self.edge_transformation(opposite_label, opposite_edge)
                p2 = self.polygon(opposite_label)
                p2 = Polygon(
                    vertices=[sim(v) for v in p2.vertices()], base_ring=p1.base_ring()
                )

                p = p1.join(p2, edge, opposite_edge)

                return p.is_convex(strict=strict)

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
                    Traceback (most recent call last):
                    ...
                    ValueError: edge must be joining two different polygons

                """
                opposite_label, opposite_edge = self.opposite_edge(label, edge)

                if opposite_label == label:
                    raise ValueError("edge must be joining two different polygons")

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
                    sage: ss = translation_surfaces.ward(3).triangulate().codomain()
                    sage: ss.random_flip(15)
                    Translation Surface in H_1(0^3) built from ...

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

            @cached_method(key=lambda self, ring: ring or self.base_ring())
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
                    SimilaritySurfaceTangentVector in polygon (1, -1, 0) based at (1/2, -3/2) with vector (1, 1)
                    sage: K.<sqrt2> = QuadraticField(2)
                    sage: S.tangent_vector(S.root(), (1/2,1/2), (1,sqrt2), ring=K)
                    SimilaritySurfaceTangentVector in polygon (1, -1, 0) based at (1/2, -3/2) with vector (1, sqrt2)
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
                Return a morphism triangulating the surface or
                ``None`` if the surface is already triangulated.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
                    sage: S.triangulation_mapping()
                    doctest:warning
                    ...
                    UserWarning: triangulation_mapping() has been deprecated and will be removed in a future version of sage-flatsurf; use triangulate() instead
                    Triangulation morphism:
                      From: Translation Surface in H_2(2) built from 3 squares
                      To:   Triangulation of Translation Surface in H_2(2) built from 3 squares

                """
                import warnings

                warnings.warn(
                    "triangulation_mapping() has been deprecated and will be removed in a future version of sage-flatsurf; use triangulate() instead"
                )

                return self.triangulate()

            # TODO: Rename to triangulation()
            def triangulate(self, in_place=False, label=None, relabel=None):
                r"""
                Return a morphism to a triangulated version of this surface.

                If label=None (the default) all polygons are triangulated. Otherwise,
                label should be a polygon label. In this case, just this polygon
                is split into triangles.

                This is done in place if in_place is True (defaults to False).

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: s=translation_surfaces.mcmullen_L(1,1,1,1)
                    sage: ss=s.triangulate().codomain()
                    sage: gs=ss.graphical_surface()
                    sage: gs.make_all_visible()
                    sage: gs
                    Graphical representation of Triangulation of Translation Surface in H_2(2) built from 3 squares

                A non-strictly convex example that caused trouble:

                    sage: from flatsurf import similarity_surfaces, Polygon
                    sage: s=similarity_surfaces.self_glued_polygon(Polygon(edges=[(1,1),(-3,-1),(1,0),(1,0)]))
                    sage: s=s.triangulate().codomain()
                    sage: len(s.polygon((0, 0)).vertices())
                    3
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

                    return MutableOrientedSimilaritySurface.from_surface(
                        self
                    ).triangulate(in_place=True, label=label)

                labels = {label} if label is not None else None

                from flatsurf.geometry.morphism import TriangulationMorphism
                from flatsurf.geometry.lazy import LazyTriangulatedSurface

                return TriangulationMorphism._create_morphism(
                    self, LazyTriangulatedSurface(self, labels=labels)
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

            # TODO: Should we rename this to is_delaunay? And remove is_delaunay_triangulated?
            def is_delaunay_decomposed(self, limit=None):
                r"""
                Return if the decomposition of the surface into polygons is Delaunay.

                INPUT:

                - ``limit`` -- an integer or ``None`` (default: ``None``);
                  check only ``limit`` many polygons if set

                """
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
                        c1 = p1.circumscribing_circle()
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

            # TODO: In the long run, what should the naming strategy be?
            # Should it say delaunay_triangulation() or delaunay_triangulate()?
            # The latter sounds like an in-place operation. So maybe
            # delaunay_triangulation() is better.
            # Should there be any in-place operations at all? For anything that
            # is linear-time or slower, it makes not a lot of sense really. If
            # we want speed, we should call out to libflatsurf.
            # TODO: Maybe the codomain should be a
            # MutableOrientedSimilaritySurface that is immutable when
            # this is finite type. (Here and for all the other functions
            # returning a lazy surface at the moment.)
            def delaunay_triangulation(
                self,
                triangulated=False,
                in_place=False,
                direction=None,
                relabel=None,
            ):
                r"""
                Returns a Delaunay triangulation of a surface, or make some
                triangle flips to get closer to the Delaunay decomposition.

                INPUT:

                - ``triangulated`` (boolean) - If true, the algorithm assumes the
                  surface is already triangulated. It does this without verification.

                - ``in_place`` (boolean) - If true, the triangulating and the
                  triangle flips are done in place.  Otherwise, a mutable copy of the
                  surface is made.

                - ``direction`` (ignored and deprecated)

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces

                    sage: m = matrix([[2,1],[1,1]])
                    sage: s = m*translation_surfaces.infinite_staircase()
                    sage: ss = s.delaunay_triangulation().codomain()
                    sage: ss.root()
                    (0, 0)
                    sage: ss.polygon((0, 0))
                    Polygon(vertices=[(0, 0), (-1, 0), (-1, -1)])
                    sage: TestSuite(ss).run()
                    sage: ss.is_delaunay_triangulated(limit=10)
                    True
                """
                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement delaunay_triangulation(in_place=True) yet"
                    )

                if direction is not None:
                    direction = None
                    import warnings

                    warnings.warn(
                        "the direction keyword argument for delaunay_triangulation() has been deprecated and will be removed in a future version of sage-flatsurf; the keyword argument has no effect"
                    )

                if relabel is not None:
                    if relabel:
                        raise NotImplementedError(
                            "the relabel keyword has been removed from delaunay_triangulation(); use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                        )
                    else:
                        import warnings

                        warnings.warn(
                            "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_triangulation()"
                        )

                morphism = self.triangulate()

                from flatsurf.geometry.lazy import (
                    LazyDelaunayTriangulatedSurface,
                )

                s = LazyDelaunayTriangulatedSurface(
                    morphism.codomain(), direction=direction, category=self.category()
                )

                from flatsurf.geometry.morphism import DelaunayTriangulationMorphism

                return (
                    DelaunayTriangulationMorphism._create_morphism(
                        morphism.codomain(), s
                    )
                    * morphism
                )

            def delaunay_decomposition(
                self,
                triangulated=False,
                delaunay_triangulated=False,
                in_place=False,
                direction=None,
                relabel=None,
            ):
                r"""
                Return the Delaunay Decomposition of this surface.

                INPUT:

                - ``triangulated`` (boolean) - If true, the algorithm assumes the
                  surface is already triangulated. It does this without verification.

                - ``delaunay_triangulated`` (boolean) - If true, the algorithm assumes
                  the surface is already delaunay_triangulated. It does this without
                  verification.

                - ``in_place`` (boolean) - If true, the triangulating and the triangle
                  flips are done in place. Otherwise, a mutable copy of the surface is
                  made.

                - ``direction`` - (ignored and deprecated)

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces, Polygon, similarity_surfaces
                    sage: s0 = translation_surfaces.octagon_and_squares()
                    sage: a = s0.base_ring().gens()[0]
                    sage: m = Matrix([[1,2+a],[0,1]])
                    sage: s = m*s0
                    sage: s = s.triangulate().codomain()
                    sage: ss = s.delaunay_decomposition(triangulated=True).codomain()
                    sage: len(ss.polygons())
                    3

                    sage: p = Polygon(edges=[(4,0),(-2,1),(-2,-1)])
                    sage: s0 = similarity_surfaces.self_glued_polygon(p)
                    sage: s = s0.delaunay_decomposition().codomain()
                    sage: TestSuite(s).run()

                    sage: m = matrix([[2,1],[1,1]])
                    sage: s = m*translation_surfaces.infinite_staircase()
                    sage: ss = s.delaunay_decomposition().codomain()
                    sage: ss.root()
                    (0, 0)
                    sage: ss.polygon(ss.root())
                    Polygon(vertices=[(0, 0), (-1, 0), (-1, -1), (0, -1)])
                    sage: TestSuite(ss).run()
                    sage: ss.is_delaunay_decomposed(limit=10)
                    True

                """
                if in_place:
                    raise NotImplementedError(
                        "this surface does not implement delaunay_decomposition(in_place=True) yet"
                    )

                if direction is not None:
                    direction = None
                    import warnings

                    warnings.warn(
                        "the direction keyword argument for delaunay_decomposition() has been deprecated and will be removed in a future version of sage-flatsurf; the keyword argument has no effect"
                    )

                if relabel is not None:
                    if relabel:
                        raise NotImplementedError(
                            "the relabel keyword has been removed from delaunay_decomposition(); use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                        )
                    else:
                        import warnings

                        warnings.warn(
                            "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_decomposition()"
                        )

                morphism = self.delaunay_triangulation()

                from flatsurf.geometry.lazy import LazyDelaunaySurface

                s = LazyDelaunaySurface(
                    morphism.codomain(), direction=direction, category=self.category()
                )
                from flatsurf.geometry.morphism import DelaunayDecompositionMorphism

                return (
                    DelaunayDecompositionMorphism._create_morphism(
                        morphism.codomain(), s
                    )
                    * morphism
                )

            def _saddle_connections_unbounded(
                self, initial_label, initial_vertex, algorithm
            ):
                r"""
                Enumerate all saddle connections in this surface ordered by length.

                This is a helper method for :meth:`saddle_connections`.
                """

                def squared_length(v):
                    return v[0] ** 2 + v[1] ** 2

                # Enumerate all saddle connections by length
                connections = set()
                shortest_edge = min(
                    squared_length(self.polygon(label).edge(edge))
                    for (label, edge) in self.edges()
                )

                length_bound = shortest_edge
                while True:
                    more_connections = [
                        connection
                        for connection in self.saddle_connections(
                            squared_length_bound=length_bound,
                            initial_label=initial_label,
                            initial_vertex=initial_vertex,
                            algorithm=algorithm,
                        )
                        if connection not in connections
                    ]
                    for connection in sorted(
                        more_connections,
                        key=lambda connection: squared_length(connection.holonomy()),
                    ):
                        connections.add(connection)
                        yield connection

                    length_bound *= 2

            def _saddle_connections_generic_cone_bounded(
                self, squared_length_bound, source, incoming_edge, similarity, cone
            ):
                r"""
                Enumerate the saddle connections of length at most square root
                of ``squared_length_bound`` and which are strictly inside the
                ``cone``.

                This is a helper method for :meth:`saddle_connections`.

                ALGORITHM:

                We check for each vertex of the polygon if it is contained in
                the cone. If it is, it leads to a saddle connection (unless
                it's hidden or too far away.)

                Then we recursively propagate the cone across each edge it hits
                into the neighboring polygons.

                INPUT:

                - ``squared_length_bound`` -- a number, the square of the
                  length up to which saddle connections should be considered;
                  the length of saddle connections is determined using their
                  holonomy vector written in their source polygon.

                - ``source`` -- a pair of a polygon label and a vertex
                  index; the vertex where all saddle connections enumerated by
                  this method start.

                - ``incoming_edge`` -- a pair of a polygon label and an edge
                  index; the ``cone`` is crossing over this ``incoming_edge``
                  (after transforming that polygon with the inverse of the
                   ``similarity``.)

                - ``similarity`` -- a similarity of the plane; describes a
                  translation, rotation, and dilation of the ``incoming_edge``
                  polygon to position it relative to the ``cone``

                - ``cone`` -- an open cone in the plane which bounds the
                  holonomy of the saddle connections.

                """
                assert not cone.is_empty()

                label = incoming_edge[0]
                polygon = similarity(self.polygon(label))

                incoming_edge_segment = (
                    polygon.vertex(incoming_edge[1]),
                    polygon.vertex(incoming_edge[1] + 1),
                )
                origin = (polygon.base_ring() ** 2).zero()

                from flatsurf.geometry.cone import Cones

                if polygon.vertex(incoming_edge[1]):
                    incoming_edge_cone = Cones(self.base_ring())(
                        polygon.vertex(incoming_edge[1] + 1),
                        polygon.vertex(incoming_edge[1]),
                    )
                else:
                    incoming_edge_cone = Cones(self.base_ring())(
                        polygon.vertex(incoming_edge[1] + 1),
                        polygon.vertex(incoming_edge[1] - 1),
                    )

                if not cone.is_subset(incoming_edge_cone):
                    raise ValueError(
                        "cone must be contained in the cone formed by the incoming edge"
                    )

                from flatsurf import EuclideanPlane

                bounding_circle = EuclideanPlane(self.base_ring()).circle(
                    (0, 0),
                    radius_squared=squared_length_bound,
                )

                # Each vertex that is contained in the cone's interior yields a
                # saddle connection (if it is "behind" the incoming edge and
                # not hidden by some other edge; these conditions are only
                # possible for polygons that are not strictly convex.)
                for v, vertex in enumerate(polygon.vertices()):
                    if cone.contains_point(vertex):
                        if v == incoming_edge[1]:
                            continue

                        from flatsurf.geometry.euclidean import (
                            time_on_ray,
                            ray_segment_intersection,
                        )

                        vertex_time_on_ray = time_on_ray(
                            origin,
                            vertex,
                            ray_segment_intersection(
                                origin, vertex, incoming_edge_segment
                            ),
                        )
                        if vertex_time_on_ray[0] > vertex_time_on_ray[1]:
                            assert not polygon.is_convex()
                            # The cone hits the vertex before entering the polygon.
                            continue

                        exit = polygon.flow_to_exit(vertex, -vertex)

                        # TODO: This is probably very inefficient. It would be enough to check if this point is on the incoming_edge.
                        exit = polygon.get_point_position(exit)

                        if exit.is_vertex() and exit.get_vertex() != incoming_edge[1]:
                            # Another vertex hides this vertex.
                            continue

                        if (
                            exit.is_in_edge_interior()
                            and exit.get_edge() != incoming_edge[1]
                        ):
                            # Another edge hides this vertex.
                            continue

                        # The vertex is not hidden by some other vertex or
                        # edge. This is a saddle connection.
                        holonomy = vertex
                        if bounding_circle.point_position(holonomy) >= 0:
                            # The saddle connection is within the squared_length_bound.
                            from flatsurf.geometry.saddle_connection import (
                                SaddleConnection,
                            )

                            yield SaddleConnection(
                                surface=self,
                                start=source,
                                end=(label, v),
                                holonomy=vertex,
                                end_holonomy=~similarity.derivative() * vertex,
                            )

                # We need to propagate the cone across edges to neighboring
                # polygons. For this, we split the cone into smaller subcones,
                # such that each subcone contains no vertex in its interior.
                cone_space = cone.parent()
                ray_space = cone_space.rays()
                vertex_directions = (
                    [cone.start()]
                    + cone.sorted_rays(
                        [
                            ray_space(vertex)
                            for v, vertex in enumerate(polygon.vertices())
                            if v != incoming_edge[1]
                            and cone.contains_point(vertex)
                            and ray_space(vertex) not in [cone.start(), cone.end()]
                        ]
                    )
                    + [cone.end()]
                )

                subcones = [
                    cone_space(v, w)
                    for (v, w) in zip(vertex_directions, vertex_directions[1:])
                ]
                # Now we propagate each subcone across the first edge it hits
                # after crossing over the incoming edge.
                for subcone in subcones:
                    ray = subcone.a_ray()
                    from flatsurf.geometry.euclidean import ray_segment_intersection

                    start = ray_segment_intersection(
                        origin, ray.vector(), incoming_edge_segment
                    )
                    exit = polygon.flow_to_exit(start, ray.vector())
                    exit = polygon.get_point_position(exit)
                    assert exit.is_in_edge_interior()
                    outgoing_edge = exit.get_edge()
                    if (
                        bounding_circle.line_segment_position(
                            polygon.vertex(outgoing_edge),
                            polygon.vertex(outgoing_edge + 1),
                        )
                        != 1
                    ):
                        # No part of the edge is inside the
                        # squared_length_bound, search ends here.
                        continue

                    opposite_edge = self.opposite_edge(label, outgoing_edge)
                    if opposite_edge is None:
                        # Unglued edge. Search ends here.
                        continue

                    # Recurse
                    yield from self._saddle_connections_generic_cone_bounded(
                        squared_length_bound,
                        source,
                        opposite_edge,
                        similarity * self.edge_transformation(*opposite_edge),
                        subcone,
                    )

            def _saddle_connections_generic_from_vertex_bounded(
                self, squared_length_bound, source
            ):
                r"""
                Enumerate all the saddle connections up to length
                ``squared_length_bound`` which start at ``source``.

                This is a helper method for :meth:`saddle_connections`.

                ALGORITHM:

                We consider saddle connections that come from the edges
                adjacent to the vertex of ``source`` and then use
                :meth:`_saddle_connections_generic_cone_bounded` to enumerate the
                saddle connections in the open cone formed by these edges.

                INPUT:

                - ``squared_length_bound`` -- a number, the square of the
                  length up to which saddle connections should be considered;
                  the length of saddle connections is determined using their
                  holonomy vector written in their source polygon.

                - ``source`` -- a pair consisting of a polygon label and a vertex
                  index. The saddle connection starts at that vertex and is
                  contained in the closed cone that is formed by the two edges
                  adjacent to the vertex.

                """
                polygon = self.polygon(source[0])

                from flatsurf.geometry.saddle_connection import SaddleConnection

                for connection in [
                    SaddleConnection.from_half_edge(self, source[0], source[1]),
                    -SaddleConnection.from_half_edge(
                        self, source[0], (source[1] - 1) % len(polygon.edges())
                    ),
                ]:
                    if connection.length_squared() <= squared_length_bound:
                        yield connection

                from flatsurf.geometry.similarity import SimilarityGroup

                G = SimilarityGroup(self.base_ring())
                similarity = G.translation(*-polygon.vertex(source[1]))

                from flatsurf.geometry.cone import Cones

                cone = Cones(polygon.base_ring())(
                    polygon.edge(source[1]), -polygon.edge(source[1] - 1)
                )

                yield from self._saddle_connections_generic_cone_bounded(
                    squared_length_bound, source, source, similarity, cone
                )

            def saddle_connections(
                self,
                squared_length_bound=None,
                initial_label=None,
                initial_vertex=None,
                algorithm=None,
            ):
                r"""
                Return the saddle connections on this surface whose length
                squared is at most ``squared_length_bound`` (ordered by
                length.)

                The length of a saddle connection is measured using holonomy
                from the polygon in which the trajectory starts.

                If no ``squared_length_bound`` is given, all saddle connections
                are enumerated (ordered by length.)

                If ``initial_label`` and ``initial_vertex`` are provided, only
                saddle connections are returned which emanate from the
                corresponding vertex of a polygon (and only pointing into the
                polygon or along the edges adjacent to that vertex.)

                If only ``initial_label`` is provided, the saddle connections
                will only emanate from vertices of the corresponding polygon.

                If only ``initial_vertex`` is provided, the saddle connections
                will only emanate from that vertex.

                EXAMPLES:

                Return the connections of length up to square root of 5::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: connections = S.saddle_connections(5)
                    sage: list(connections)
                    [Saddle connection (0, -1) from vertex 3 of polygon 0 to vertex 1 of polygon 0,
                     Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 2 of polygon 0,
                     Saddle connection (0, 1) from vertex 1 of polygon 0 to vertex 3 of polygon 0,
                     Saddle connection (-1, 0) from vertex 2 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (1, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0,
                     Saddle connection (1, -1) from vertex 3 of polygon 0 to vertex 1 of polygon 0,
                     Saddle connection (-1, 1) from vertex 1 of polygon 0 to vertex 3 of polygon 0,
                     Saddle connection (-1, -1) from vertex 2 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (-1, 2) from vertex 1 of polygon 0 to vertex 3 of polygon 0,
                     Saddle connection (2, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0,
                     Saddle connection (1, -2) from vertex 3 of polygon 0 to vertex 1 of polygon 0,
                     Saddle connection (-2, 1) from vertex 1 of polygon 0 to vertex 3 of polygon 0,
                     Saddle connection (-2, -1) from vertex 2 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (-1, -2) from vertex 2 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (2, -1) from vertex 3 of polygon 0 to vertex 1 of polygon 0,
                     Saddle connection (1, 2) from vertex 0 of polygon 0 to vertex 2 of polygon 0]

                We get the same result if we take the first 16 saddle
                connections without a length bound::

                    sage: from itertools import islice
                    sage: set(connections) == set(islice(S.saddle_connections(), 16))
                    True

                While enumerating saddle connections without a bound is not
                asymptotically slower than enumerating with a bound, in the
                current implementation it is quite a bit slower in practice in
                particular if the bound is small.

                TESTS:

                Verify that saddle connections are enumerated correctly when
                there are unglued edges::

                    sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (0, 1)]))
                    0
                    sage: S.set_immutable()
                    sage: len(list(S.saddle_connections(10)))
                    6

                Verify that saddle connections are enumerated correctly when
                there are self-glued edges::

                    sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (0, 1)]))
                    0
                    sage: S.glue((0, 0), (0, 0))
                    sage: S.glue((0, 1), (0, 1))
                    sage: S.glue((0, 2), (0, 2))
                    sage: S.set_immutable()

                    sage: from itertools import islice
                    sage: list(islice(S.saddle_connections(), 8))
                    [Saddle connection (0, -1) from vertex 2 of polygon 0 to vertex 2 of polygon 0,
                     Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (1, 1) from vertex 0 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (-1, 1) from vertex 1 of polygon 0 to vertex 1 of polygon 0,
                     Saddle connection (2, 1) from vertex 0 of polygon 0 to vertex 0 of polygon 0,
                     Saddle connection (1, -2) from vertex 2 of polygon 0 to vertex 2 of polygon 0,
                     Saddle connection (-2, 1) from vertex 1 of polygon 0 to vertex 1 of polygon 0,
                     Saddle connection (1, 2) from vertex 0 of polygon 0 to vertex 0 of polygon 0]

                    sage: len(list(S.saddle_connections(1)))
                    2

                    sage: len(list(S.saddle_connections(1, initial_label=0, initial_vertex=1)))
                    1

                We can also enumerate saddle connections on surfaces that are
                built from non-convex polygons such as this L shaped polygon::

                    sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
                    sage: L = MutableOrientedSimilaritySurface(QQ)
                    sage: L.add_polygon(Polygon(vertices=[(0, 0), (3, 0), (7, 0), (7, 2), (3, 2), (3, 3), (0, 3), (0, 2)]))
                    0
                    sage: L.glue((0, 0), (0, 5))
                    sage: L.glue((0, 1), (0, 3))
                    sage: L.glue((0, 2), (0, 7))
                    sage: L.glue((0, 4), (0, 6))
                    sage: L.set_immutable()

                    sage: connections = L.saddle_connections(128)
                    sage: len(connections)
                    164

                Note that on translation surfaces, enumerating saddle
                connections with the (default) ``"pyflatsurf"`` algorithm is
                usually much faster than the ``"generic"`` algorithm::

                    sage: connections = L.saddle_connections(128)
                    sage: len(connections)
                    164

                    sage: connections = L.saddle_connections(128, algorithm="generic")
                    sage: len(connections)
                    164

                """
                # TODO: Add benchmarks of the generic "cone" algorithm against
                # the pyflatsurf algorithm. Also benchmark how much slower this
                # is now since we are supporting much more complicated
                # geometries.

                if squared_length_bound is not None and squared_length_bound < 0:
                    raise ValueError("length bound must be non-negative")

                if algorithm is None:
                    algorithm = "generic"

                if algorithm == "generic":
                    return self._saddle_connections_generic(
                        squared_length_bound, initial_label, initial_vertex
                    )

                raise NotImplementedError(
                    "cannot enumerate saddle connections with this algorithm yet"
                )

            def _saddle_connections_generic(
                self, squared_length_bound, initial_label, initial_vertex
            ):
                if squared_length_bound is None:
                    # Enumerate all (usually infinitely many) saddle connections.
                    return self._saddle_connections_unbounded(
                        initial_label=initial_label,
                        initial_vertex=initial_vertex,
                        algorithm="generic",
                    )

                connections = []

                if initial_label is None and initial_vertex is None:
                    if not self.is_finite_type():
                        raise NotImplementedError(
                            "cannot enumerate saddle connections on surfaces that are built from inifinitely many polygons yet"
                        )

                    initial = [(label, range(len(self.polygon(label).vertices()))) for label in self.labels()]
                elif initial_label is None:
                    representatives = [(label, self.polygon(label).get_point_position(coordinates).get_vertex()) for (label, coordinates) in initial_vertex.representatives()]
                    labels = {label for (label, _) in representatives}
                    initial = [(label, [vertex for (lbl, vertex) in representatives if lbl == label]) for label in labels]
                elif initial_vertex is None:
                    initial = [(initial_label, range(len(self.polygon(initial_label).vertices())))]
                else:
                    initial = [(initial_label, [initial_vertex])]

                for label, vertices in initial:
                    for vertex in vertices:
                        connections.extend(
                            self._saddle_connections_generic_from_vertex_bounded(
                                squared_length_bound=squared_length_bound,
                                source=(label, vertex),
                            )
                        )

                # The connections might contain duplicates because each glued
                # edge can show up in two different polygons. Note that we
                # cannot just change _saddle_connections_from_vertex_bounded()
                # to only consider the edge clockwise from the vertex since we
                # would then miss saddle connections in surfaces with
                # self-glued and unglued edges.
                connections = set(connections)

                return sorted(
                    connections, key=lambda connection: connection.length_squared()
                )

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
                # TODO: Returns a morphism actually.
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

                    sage: T = S.subdivide().codomain()
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

                    sage: T = S.subdivide().codomain()

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
                return self.insert_marked_points(
                    *[
                        self(label, self.polygon(label).centroid())
                        for label in self.labels()
                    ]
                )

            def insert_marked_points(self, *points):
                if self.is_mutable():
                    from flatsurf import MutableOrientedSimilaritySurface

                    copy = MutableOrientedSimilaritySurface.from_surface(self)
                    copy.set_immutable()

                    codomain = copy.insert_marked_points(*points).codomain()

                    # Since the domain is mutable, the morphism is not
                    # functional but only has a .codomain().
                    from flatsurf.geometry.morphism import SurfaceMorphism

                    return SurfaceMorphism._create_morphism(None, codomain)

                from flatsurf.geometry.morphism import IdentityMorphism

                morphism = IdentityMorphism._create_morphism(self)

                for p in points:
                    if p.is_vertex():
                        raise ValueError("cannot insert marked points at vertices")

                edge_points = [p for p in points if p.is_in_edge_interior()]
                face_points = [p for p in points if p not in edge_points]

                assert len(edge_points) + len(face_points) == len(points)

                if edge_points:
                    insert_morphism = morphism.codomain()._insert_marked_points_edges(
                        *edge_points
                    )
                    morphism = insert_morphism * morphism
                    face_points = [insert_morphism(point) for point in face_points]
                if face_points:
                    insert_morphism = morphism.codomain()._insert_marked_points_faces(
                        *face_points
                    )
                    morphism = insert_morphism * morphism

                return morphism

            def _insert_marked_points_edges(self, *points):
                assert (
                    points
                ), "_insert_marked_points_edges must be called with some points to insert"

                from flatsurf.geometry.euclidean import time_on_ray

                points = {
                    label: {
                        # TODO: Sort vertices along edge
                        edge: sorted(
                            [
                                coordinates
                                for point in points
                                for (lbl, coordinates) in point.representatives()
                                if lbl == label
                                and self.polygon(label)
                                .get_point_position(coordinates)
                                .get_edge()
                                == edge
                            ],
                            key=lambda coordinates: time_on_ray(
                                self.polygon(label).vertex(edge),
                                self.polygon(label).edge(edge),
                                coordinates,
                            ),
                        )
                        for edge in range(len(self.polygon(label).edges()))
                    }
                    for label in self.labels()
                }

                from flatsurf import MutableOrientedSimilaritySurface

                surface = MutableOrientedSimilaritySurface(self.base_ring())

                # Add polygons to surface with marked point
                for label in self.labels():
                    vertices = []
                    for v, vertex in enumerate(self.polygon(label).vertices()):
                        vertices.append(vertex)
                        vertices.extend(points[label][v])

                    from flatsurf import Polygon

                    surface.add_polygon(Polygon(vertices=vertices), label=label)

                # TODO: Make static on InsertMarkedPointsOnEdgeMorphism
                def edgenum(label, edge, section):
                    edgenum = 0
                    for e in range(edge):
                        edgenum += len(points[label][e]) + 1

                    edgenum += section

                    return edgenum

                # Glue polygons in surface.
                for label in self.labels():
                    for edge in range(len(self.polygon(label).edges())):
                        opposite = self.opposite_edge(label, edge)
                        if opposite is None:
                            continue

                        opposite_label, opposite_edge = opposite
                        for e in range(len(points[label][edge]) + 1):
                            surface.glue(
                                (label, edgenum(label, edge, e)),
                                (
                                    opposite_label,
                                    edgenum(opposite_label, opposite_edge + 1, -1 - e),
                                ),
                            )

                surface.set_immutable()

                from flatsurf.geometry.morphism import InsertMarkedPointsOnEdgeMorphism

                return InsertMarkedPointsOnEdgeMorphism._create_morphism(
                    self, surface, points
                )

            def _insert_marked_points_faces(self, *points):
                # Recursively insert points by only inserting at most one point
                # in each face at a time.
                first_point = {}
                more_points = {}

                for point in points:
                    label, coordinates = point.representative()
                    if label not in first_point:
                        first_point[label] = point.coordinates(label)[0]
                    else:
                        more_points[label] = more_points.get(label, []) + [point]

                assert (
                    first_point
                ), "_insert_marked_points_faces must be called with some points to insert"

                def is_subdivided(label):
                    return label in first_point

                subdivisions = {
                    label: self.polygon(label).subdivide(first_point[label])
                    if is_subdivided(label)
                    else [self.polygon(label)]
                    for label in self.labels()
                }

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                surface = MutableOrientedSimilaritySurface(self.base())

                # Add subdivided polygons
                for label in self.labels():
                    if is_subdivided(label):
                        for p, polygon in enumerate(subdivisions[label]):
                            surface.add_polygon(polygon, label=(label, p))
                    else:
                        surface.add_polygon(self.polygon(label), label=label)

                surface.set_roots(
                    (label, 0) if is_subdivided(label) else label
                    for label in self.roots()
                )

                # Establish gluings
                for label in self.labels():
                    for e in range(len(self.polygon(label).vertices())):
                        # Reestablish the original gluings
                        opposite = self.opposite_edge(label, e)

                        if opposite is not None:
                            opposite_label, opposite_edge = opposite
                            surface.glue(
                                (
                                    (label, e) if is_subdivided(label) else label,
                                    0 if is_subdivided(label) else e,
                                ),
                                (
                                    (opposite_label, opposite_edge)
                                    if is_subdivided(opposite_label)
                                    else opposite_label,
                                    0
                                    if is_subdivided(opposite_label)
                                    else opposite_edge,
                                ),
                            )

                        # Glue subdivided polygons internally
                        if is_subdivided(label):
                            surface.glue(
                                ((label, e), 1),
                                (
                                    (
                                        label,
                                        (e + 1) % len(self.polygon(label).vertices()),
                                    ),
                                    2,
                                ),
                            )

                surface.set_immutable()

                from flatsurf.geometry.morphism import InsertMarkedPointsInFaceMorphism

                insert_first_point = InsertMarkedPointsInFaceMorphism._create_morphism(
                    self, surface, subdivisions
                )

                if more_points:
                    from itertools import chain

                    more_points = list(chain.from_iterable(more_points.values()))
                    more_points = [insert_first_point(p) for p in more_points]
                    return (
                        insert_first_point.codomain().insert_marked_points(*more_points)
                        * insert_first_point
                    )

                return insert_first_point

            def subdivide_edges(self, parts=2):
                r"""
                # TODO: Returns a morphism actually.
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

                    sage: T = S.subdivide_edges().codomain()

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

                from collections.abc import Iterable

                if not isinstance(parts, Iterable):
                    parts = tuple(1 / parts for k in range(parts))
                if sum(parts) != 1:
                    raise ValueError

                subdivideds = [
                    p.subdivide_edges(parts=[parts for _ in p.edges()])
                    for p in polygons
                ]

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
                            for p in range(len(parts)):
                                surface.glue(
                                    (label, e * len(parts) + p),
                                    (
                                        opposite[0],
                                        opposite[1] * len(parts) + (len(parts) - p - 1),
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
                    True

                """
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
                    limit = 32

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
                    Return a morphism to a surface with each polygon replaced
                    with a new polygon which differs by translation and
                    reindexing. The new polygon will have the property that
                    vertex zero is the origin, and each vertex lies in the
                    upper half plane or on the x-axis with non-negative
                    x-coordinate.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: s=translation_surfaces.veech_double_n_gon(4)
                        sage: s.polygon(1)
                        Polygon(vertices=[(0, 0), (-1, 0), (-1, -1), (0, -1)])
                        sage: [s.opposite_edge(0,i) for i in range(4)]
                        [(1, 0), (1, 1), (1, 2), (1, 3)]
                        sage: ss = s.standardize_polygons().codomain()
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
                    morphism = S.standardize_polygons(in_place=True)
                    S.set_immutable()
                    return morphism.change(domain=self, codomain=S)

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
