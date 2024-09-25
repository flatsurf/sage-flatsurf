r"""
The category of dilation surfaces.

This module provides shared functionality for all surfaces in sage-flatsurf
that are built from Euclidean polygons that are glued by translation followed
by homothety, i.e., application of a diagonal matrix.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
    sage: S = similarity_surfaces.self_glued_polygon(P)
    sage: C = S.category()

    sage: from flatsurf.geometry.categories import DilationSurfaces
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

from sage.categories.category_with_axiom import all_axioms

from flatsurf.cache import cached_surface_method
from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)


class DilationSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built from polygons with edges identified by
    translations and homothety.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import DilationSurfaces
        sage: DilationSurfaces()
        Category of dilation surfaces

    """

    def super_categories(self):
        r"""
        Return the categories that a dilation surfaces is also always a member
        of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import DilationSurfaces
            sage: DilationSurfaces().super_categories()
            [Category of rational similarity surfaces]

        """
        from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces

        return [SimilaritySurfaces().Rational()]

    class Positive(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by dilation surfaces that use homothety with
        positive scaling.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: 'Positive' in S.category().axioms()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all positive dilation surfaces.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def is_dilation_surface(self, positive=False):
                r"""
                Return whether this surface is a dilation surface.

                See
                :meth:`.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_dilation_surface`
                for details.
                """
                return True

            def _test_positive_dilation_surface(self, **options):
                r"""
                Verify that this is a positive dilation surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S.set_immutable()
                    sage: S._test_positive_dilation_surface()

                """
                tester = self._tester(**options)

                limit = None

                if not self.is_finite_type():
                    from flatsurf import MutableOrientedSimilaritySurface

                    self = MutableOrientedSimilaritySurface.from_surface(
                        self, labels=self.labels()[:32]
                    )

                tester.assertTrue(
                    DilationSurfaces.ParentMethods._is_dilation_surface(
                        self, positive=True, limit=limit
                    )
                )

    class SubcategoryMethods:
        def Positive(self):
            r"""
            Return the subcategory of surfaces glued by positive dilation.

            EXAMPLES::

                sage: from flatsurf.geometry.categories import DilationSurfaces
                sage: C = DilationSurfaces()
                sage: C.Positive()
                Category of positive dilation surfaces

            """
            return self._with_axiom("Positive")

    class ParentMethods:
        r"""
        Provides methods available to all dilation surfaces.

        If you want to add functionality for such surfaces you most likely want
        to put it here.
        """

        def is_dilation_surface(self, positive=False):
            r"""
            Return whether this surface is a dilation surface.

            See :meth:`.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_dilation_surface`
            for details.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()

                sage: S.is_dilation_surface(positive=True)
                True
                sage: S.is_dilation_surface(positive=False)
                True

            """
            if not positive:
                return True

            # We do not know whether this surface is a positive dilation
            # surface or not so we have to rely on the generic implementation
            # of this.
            # pylint: disable-next=bad-super-call
            return super(DilationSurfaces().parent_class, self).is_dilation_surface(
                positive=positive
            )

        @staticmethod
        def _is_dilation_surface(surface, positive=False, limit=None):
            r"""
            Return whether ``surface`` is a dilation surface by checking how
            its polygons are glued.

            This is a helper method for
            :meth:`~.similarity_surfaces.ParentMethods.is_dilation_surface`.

            INPUT:

            - ``surface`` -- an oriented similarity surface

            - ``positive`` -- a boolean (default: ``False``); whether the
              entries of the diagonal matrix must be positive or are allowed to
              be negative.

            - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
              the first ``limit`` polygons are checked

            EXAMPLES::

                sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
                sage: S = translation_surfaces.infinite_staircase()

                sage: from flatsurf.geometry.categories import DilationSurfaces
                sage: DilationSurfaces.ParentMethods._is_dilation_surface(S, limit=8)
                doctest:warning
                ...
                UserWarning: limit has been deprecated as a keyword argument for _is_dilation_surface() and will be removed from a future version of sage-flatsurf; ...
                True
                sage: DilationSurfaces.ParentMethods._is_dilation_surface(MutableOrientedSimilaritySurface.from_surface(S, labels=S.labels()[:8]))
                True

            ::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

                sage: DilationSurfaces.ParentMethods._is_dilation_surface(S, positive=True)
                False
                sage: DilationSurfaces.ParentMethods._is_dilation_surface(S)
                True

            """
            if "Oriented" not in surface.category().axioms():
                raise NotImplementedError(
                    "cannot decide whether a non-oriented surface is dilation surface yet"
                )

            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for _is_dilation_surface() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

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
                    from flatsurf.geometry.categories import SimilaritySurfaces

                    matrix = SimilaritySurfaces.Oriented.ParentMethods.edge_matrix.f(  # pylint: disable=no-member
                        surface, label, edge
                    )

                    if not matrix.is_diagonal():
                        return False

                    if positive:
                        if matrix[0][0] < 0 or matrix[1][1] < 0:
                            return False

            return True

        def is_veering_triangulated(self, certificate=False):
            r"""
            Return whether this dilation surface is given by a veering triangulation.

            A triangulation is *veering* if none of its triangles are made of
            three edges of the same slope (positive or negative).

            EXAMPLES::

                sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon

                sage: s = MutableOrientedSimilaritySurface(QQ)
                sage: s.add_polygon(Polygon(vertices=[(0,0), (1,1), (-1,2)]))
                0
                sage: s.add_polygon(Polygon(vertices=[(0,0), (-1,2), (-2,1)]))
                1
                sage: s.glue((0, 0), (1, 1))
                sage: s.glue((0, 1), (1, 2))
                sage: s.glue((0, 2), (1, 0))
                sage: s.set_immutable()
                sage: s.is_veering_triangulated()
                True

                sage: m = matrix(ZZ, 2, 2, [5, 3, 3, 2])
                sage: (m * s).is_veering_triangulated()
                False

            """
            from flatsurf.geometry.euclidean import slope

            for label in self.labels():
                p = self.polygon(label)
                edges = p.edges()
                if len(edges) != 3:
                    return (False, label) if certificate else False
                s0, s1, s2 = map(slope, edges)
                if s0 == s1 == s2:
                    return (False, label) if certificate else False
            return (True, None) if certificate else True

        def _delaunay_edge_needs_flip_Linfinity(self, p1, e1, p2, e2):
            r"""
            Return whether the provided edge which bounds two triangles should be flipped
            to get closer to the L-infinity Delaunay decomposition.

            The return code is either ``0``: no flip needed, ``2``: flip needed to make
            it veering, ``1``: flip needed to make it L-infinity.

            TESTS::

                sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
                sage: s = MutableOrientedSimilaritySurface(QQ)
                sage: s.add_polygon(Polygon(vertices=[(0,0), (1,0), (0,1)]))
                0
                sage: s.add_polygon(Polygon(vertices=[(1,1), (0,1), (1,0)]))
                1
                sage: s.glue((0, 0), (1, 0))
                sage: s.glue((0, 1), (1, 1))
                sage: s.glue((0, 2), (1, 2))
                sage: s.set_immutable()
                sage: [s._delaunay_edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [0, 0, 0]

                sage: ss = matrix(2, [1,1,0,1]) * s
                sage: [ss._delaunay_edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [0, 0, 0]
                sage: ss = matrix(2, [1,0,1,1]) * s
                sage: [ss._delaunay_edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [0, 0, 0]

                sage: ss = matrix(2, [1,2,0,1]) * s
                sage: [ss._delaunay_edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [0, 0, 2]

                sage: ss = matrix(2, [1,0,2,1]) * s
                sage: [ss._delaunay_edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [1, 0, 0]

            """
            assert self.opposite_edge(p1, e1) == (p2, e2), "not opposite edges"

            # triangles
            poly1 = self.polygon(p1)
            poly2 = self.polygon(p2)
            if len(poly1.vertices()) != 3 or len(poly2.vertices()) != 3:
                raise ValueError("edge must be adjacent to two triangles")

            edge1 = poly1.edge(e1)
            edge1L = poly1.edge(e1 - 1)
            edge1R = poly1.edge(e1 + 1)
            edge2 = poly2.edge(e2)
            edge2L = poly2.edge(e2 - 1)
            edge2R = poly2.edge(e2 + 1)

            # if one of the triangle is monochromatic and the edge is the largest
            sim = self.edge_transformation(p2, e2)
            m = sim.derivative()  # matrix carrying p2 to p1
            if not m.is_one():
                edge2 = m * edge2
                edge2L = m * edge2L
                edge2R = m * edge2R

            # convexity check of the quadrilateral
            from flatsurf.geometry.euclidean import ccw, slope

            if ccw(edge2L, edge1R) <= 0 or ccw(edge1L, edge2R) <= 0:
                return False

            s = slope(edge1)
            s1L = slope(edge1L)
            s1R = slope(edge1R)
            s2L = slope(edge2L)
            s2R = slope(edge2R)
            x = abs(edge1[0])
            y = abs(edge1[1])
            edge1new = edge2L + edge1R
            xnew = abs(edge1new[0])
            ynew = abs(edge1new[1])
            monochromatic1 = s == s1L == s1R
            monochromatic2 = s == s2L == s2R

            # 1. flip to get closer to a veering triangulation
            if monochromatic1 and monochromatic2:
                # two monochromatic triangles
                m = max(x, y)
                mnew = max(xnew, ynew)
                snew = slope(edge1new)
                killed_monochromatic = snew != s
                needs_flip = killed_monochromatic or mnew < m
                return 2 if needs_flip else 0
            if monochromatic1:
                # only the first triangle is monochromatic
                needs_flip = (x == abs(edge1L[0]) + abs(edge1R[0])) and (
                    y == abs(edge1L[1]) + abs(edge1R[1])
                )
                return 2 if needs_flip else 0
            if monochromatic2:
                # only the second triangle is monochromatic
                needs_flip = (x == abs(edge2L[0]) + abs(edge2R[0])) and (
                    y == abs(edge2L[1]) + abs(edge2R[1])
                )
                return 2 if needs_flip else 0

            # 2. flip to get closer to the l-infinity delaunay
            if s1L == s2L == 1 and s1R == s2R == -1:
                # veering backward flip
                assert x == abs(edge1L[0]) + abs(edge1R[0])
                assert x == abs(edge2L[0]) + abs(edge2R[0])
                return 1 if ynew < x else 0
            if s1L == s2L == -1 and s1R == s2R == 1:
                # veering forward flip
                assert y == abs(edge1L[1]) + abs(edge1R[1])
                assert y == abs(edge2L[1]) + abs(edge2R[1])
                return 1 if xnew < y else 0

            return 0

        def _test_dilation_surface(self, **options):
            r"""
            Verify that this is a dilation surface.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_dilation_surface()

            """
            tester = self._tester(**options)

            limit = None

            if not self.is_finite_type():
                from flatsurf import MutableOrientedSimilaritySurface

                self = MutableOrientedSimilaritySurface.from_surface(
                    self, labels=self.labels()[:32]
                )

            tester.assertTrue(
                DilationSurfaces.ParentMethods._is_dilation_surface(
                    self, positive=False
                )
            )

        @cached_surface_method
        def affine_automorphism_group(self):
            r"""
            Return the group of affine automorphisms of this surface, i.e., the
            group of homeomorphisms that can be locally expressed as affine
            transformations.

            EXAMPLES::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: A = S.affine_automorphism_group(); A
                AffineAutomorphismGroup(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

            TESTS:

            This group is uniquely attached to a surface::

                sage: A is S.affine_automorphism_group()
                True

            """
            if self.is_mutable():
                raise NotImplementedError(
                    "affine automorphism group only implemented for immutable surfaces"
                )

            from flatsurf.geometry.veech_group import AffineAutomorphismGroup_generic

            return AffineAutomorphismGroup_generic(self)

        @cached_surface_method
        def veech_group(self):
            r"""
            Return the Veech group of this surface, i.e., the group of matrices
            that fix the vertices of this surface.

            EXAMPLES::

                sage: from flatsurf import dilation_surfaces
                sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: V = S.veech_group(); V
                VeechGroup(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

            TESTS:

            This group is uniquely attached to a surface::

                sage: V is S.veech_group()
                True

            """
            if self.is_mutable():
                raise NotImplementedError(
                    "affine automorphism group only implemented for immutable surfaces"
                )

            from flatsurf.geometry.veech_group import VeechGroup_generic

            return VeechGroup_generic(self)

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The category of dilation surfaces built from a finite number of polygons.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
            sage: S = similarity_surfaces.self_glued_polygon(P)

            sage: from flatsurf.geometry.categories import DilationSurfaces
            sage: S in DilationSurfaces().FiniteType()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all dilation surfaces built from
            finitely many polygons.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def l_infinity_delaunay_triangulation(
                self, triangulated=None, in_place=None, limit=None, direction=None
            ):
                r"""
                Return an L-infinity Delaunay triangulation of a surface, or make
                some triangle flips to get closer to the Delaunay decomposition.

                INPUT:

                - ``triangulated`` -- deprecated and ignored.

                - ``in_place`` -- deprecated  and must be ``None`` (the default);
                  otherwise an error is produced

                - ``limit`` -- optional (positive integer) If provided, then at most ``limit``
                    many diagonal flips will be done.

                - ``direction`` -- optional (vector). Used to determine labels when a
                  pair of triangles is flipped. Each triangle has a unique separatrix
                  which points in the provided direction or its negation. As such a
                  vector determines a sign for each triangle.  A pair of adjacent
                  triangles have opposite signs. Labels are chosen so that this sign is
                  preserved (as a function of labels).

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: s0 = translation_surfaces.veech_double_n_gon(5)
                    sage: field = s0.base_ring()
                    sage: a = field.gen()
                    sage: m = matrix(field, 2, [2,a,1,1])
                    sage: for _ in range(5): assert s0.triangulate().codomain().random_flip(5).l_infinity_delaunay_triangulation().is_veering_triangulated()

                    sage: s = m*s0
                    sage: s = s.l_infinity_delaunay_triangulation()
                    sage: assert s.is_veering_triangulated()
                    sage: TestSuite(s).run()

                    sage: s = (m**2)*s0
                    sage: s = s.l_infinity_delaunay_triangulation()
                    sage: assert s.is_veering_triangulated()
                    sage: TestSuite(s).run()

                    sage: s = (m**3)*s0
                    sage: s = s.l_infinity_delaunay_triangulation()
                    sage: assert s.is_veering_triangulated()
                    sage: TestSuite(s).run()

                The octagon which has horizontal and vertical edges::

                    sage: t0 = translation_surfaces.regular_octagon()
                    sage: for _ in range(5): assert t0.triangulate().codomain().random_flip(5).l_infinity_delaunay_triangulation().is_veering_triangulated()
                    sage: r = matrix(t0.base_ring(), [
                    ....:    [ sqrt(2)/2, -sqrt(2)/2 ],
                    ....:    [ sqrt(2)/2, sqrt(2)/2 ]])
                    sage: p = matrix(t0.base_ring(), [
                    ....:    [ 1, 2+2*sqrt(2) ],
                    ....:    [ 0, 1 ]])

                    sage: t = (r * p * r * t0).l_infinity_delaunay_triangulation()
                    sage: systole_count = 0
                    sage: for l, e in t.edges():
                    ....:     v = t.polygon(l).edge(e)
                    ....:     systole_count += v[0]**2 + v[1]**2 == 1
                    sage: assert t.is_veering_triangulated() and systole_count == 8, (systole_count, {lab: t.polygon(lab) for lab in t.labels()}, t.gluings())

                    sage: t = (r**4 * p * r**5 * p**2 * r * t0).l_infinity_delaunay_triangulation()
                    sage: systole_count = 0
                    sage: for l, e in t.edges():
                    ....:     v = t.polygon(l).edge(e)
                    ....:     systole_count += v[0]**2 + v[1]**2 == 1
                    sage: assert t.is_veering_triangulated() and systole_count == 8, (systole_count, {lab: t.polygon(lab) for lab in t.labels()}, t.gluings())

                TESTS:

                Verify that deprecated keywords do not cause errors::

                    sage: s.l_infinity_delaunay_triangulation(triangulated=True)  # long time (.8s)
                    doctest:warning
                    ...
                    UserWarning: The triangulated keyword of l_infinity_delaunay_triangulation() has been deprecated and will be removed from a future version of sage-flatsurf. The keyword has no effect anymore.
                    Translation Surface in H_2(2) built from 6 triangles
                    sage: s.l_infinity_delaunay_triangulation(triangulated=False)  # long time (.9s)
                    doctest:warning
                    ...
                    UserWarning: The triangulated keyword of l_infinity_delaunay_triangulation() has been deprecated and will be removed from a future version of sage-flatsurf. The keyword has no effect anymore.
                    Translation Surface in H_2(2) built from 6 triangles

                ::

                    sage: s.l_infinity_delaunay_triangulation(in_place=True)
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: The in_place keyword for l_infinity_delaunay_triangulation() is not supported anymore. It did not work correctly in previous versions of sage-flatsurf and will be fully removed in a future version of sage-flatsurf.

                """
                if triangulated is not None:
                    import warnings

                    warnings.warn(
                        "The triangulated keyword of l_infinity_delaunay_triangulation() has been deprecated and will be removed from a future version of sage-flatsurf. The keyword has no effect anymore."
                    )
                if in_place is not None:
                    raise NotImplementedError(
                        "The in_place keyword for l_infinity_delaunay_triangulation() is not supported anymore. It did not work correctly in previous versions of sage-flatsurf and will be fully removed in a future version of sage-flatsurf."
                    )

                if direction is not None:
                    import warnings

                    warnings.warn(
                        "the direction parameter of l_infinity_delaunay_triangulation() has been removed since it did not work correctly in previous versions of sage-flatsurf"
                    )

                return self.veering_triangulation(l_infinity=True, limit=limit)

            def veering_triangulation(
                self, l_infinity=False, limit=None, direction=None
            ):
                r"""
                Return a veering triangulated surface, or make some triangle
                flips to get closer to the Delaunay decomposition.

                INPUT:

                - ``l_infinity`` -- optional (boolean, default ``False``).
                  Whether to return a L^oo-Delaunay triangulation or any
                  veering triangulation.

                - ``limit`` -- optional (positive integer) If provided, then at most ``limit``
                    many diagonal flips will be done.

                - ``direction`` -- optional (vector). Used to determine labels when a
                  pair of triangles is flipped. Each triangle has a unique separatrix
                  which points in the provided direction or its negation. As such a
                  vector determines a sign for each triangle.  A pair of adjacent
                  triangles have opposite signs. Labels are chosen so that this sign is
                  preserved (as a function of labels).

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: s0 = translation_surfaces.veech_double_n_gon(5)
                    sage: field = s0.base_ring()
                    sage: a = field.gen()
                    sage: m = matrix(field, 2, [2,a,1,1])
                    sage: for _ in range(5): assert s0.triangulate().codomain().random_flip(5).veering_triangulation().is_veering_triangulated()

                    sage: s = m*s0
                    sage: s = s.veering_triangulation()
                    sage: assert s.is_veering_triangulated()
                    sage: TestSuite(s).run()

                    sage: s = (m**2)*s0
                    sage: s = s.veering_triangulation()
                    sage: assert s.is_veering_triangulated()
                    sage: TestSuite(s).run()

                    sage: s = (m**3)*s0
                    sage: s = s.veering_triangulation()
                    sage: assert s.is_veering_triangulated()
                    sage: TestSuite(s).run()

                The octagon which has horizontal and vertical edges::

                    sage: t0 = translation_surfaces.regular_octagon().triangulate().codomain()
                    sage: t0.is_veering_triangulated()
                    False
                    sage: t0.veering_triangulation().is_veering_triangulated()
                    True
                    sage: for _ in range(5): assert t0.random_flip(5).veering_triangulation().is_veering_triangulated()

                    sage: r = matrix(t0.base_ring(), [
                    ....:    [ sqrt(2)/2, -sqrt(2)/2 ],
                    ....:    [ sqrt(2)/2, sqrt(2)/2 ]])
                    sage: p = matrix(t0.base_ring(), [
                    ....:    [ 1, 2+2*sqrt(2) ],
                    ....:    [ 0, 1 ]])
                    sage: t = (r * p * r * t0).veering_triangulation()
                    sage: assert t.is_veering_triangulated()
                    sage: t = (r**4 * p * r**5 * p**2 * r * t0).veering_triangulation()
                    sage: assert t.is_veering_triangulated()
                """
                if direction is not None:
                    import warnings

                    warnings.warn(
                        "the direction parameter of veering_triangulation() has been removed since it did not work correctly in previous versions of sage-flatsurf"
                    )

                self = self.triangulate().codomain()

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                self = MutableOrientedSimilaritySurface.from_surface(
                    self, category=DilationSurfaces()
                )

                flip_bound = 1 if l_infinity else 2

                triangles = set(self.labels())
                if limit is None:
                    limit = -1
                else:
                    limit = int(limit)
                while triangles and limit:
                    p1 = triangles.pop()
                    for e1 in range(3):
                        p2, e2 = self.opposite_edge(p1, e1)
                        needs_flip = self._delaunay_edge_needs_flip_Linfinity(
                            p1, e1, p2, e2
                        )
                        if needs_flip >= flip_bound:
                            self.triangle_flip(p1, e1, in_place=True)
                            triangles.add(p1)
                            triangles.add(p2)
                            limit -= 1
                            if not limit:
                                break
                self.set_immutable()
                assert self.is_triangulated()
                return self


# Currently, there is no "Positive" axiom in SageMath so we make it known to
# the category framework.
all_axioms += ("Positive",)
