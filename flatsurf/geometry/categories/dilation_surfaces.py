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

from flatsurf.geometry.categories.surface_category import SurfaceCategory, SurfaceCategoryWithAxiom
from sage.categories.category_with_axiom import all_axioms


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
        from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
        return [SimilaritySurfaces()]

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

                if not self.is_finite():
                    limit = 32

                tester.assertTrue(DilationSurfaces.ParentMethods._is_dilation_surface(self, positive=True, limit=limit))

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
        @staticmethod
        def _is_dilation_surface(surface, positive=False, limit=None):
            r"""
            Return whether ``surface`` is a dilation surface by checking how
            its polygons are glued.

            This is a helper method for
            :meth:`flatsurf.geometry.categories.simililarity_surfaces.ParentMethods.is_dilation_surface`.

            INPUT:

            - ``surface`` -- an oriented similarity surface

            - ``positive`` -- a boolean (default: ``False``); whether the
              entries of the diagonal matrix must be positive or are allowed to
              be negative.

            - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
              the first ``limit`` polygons are checked

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()

                sage: from flatsurf.geometry.categories import DilationSurfaces
                sage: DilationSurfaces.ParentMethods._is_dilation_surface(S, limit=8)
                True

            ::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons((2, 0),(-1, 3),(-1, -3))
                sage: S = similarity_surfaces.self_glued_polygon(P)

                sage: DilationSurfaces.ParentMethods._is_dilation_surface(S, positive=True)
                False
                sage: DilationSurfaces.ParentMethods._is_dilation_surface(S)
                True

            """
            if 'Oriented' not in surface.category().axioms():
                raise NotImplementedError

            labels = surface.labels()

            if limit is not None:
                from itertools import islice
                labels = islice(labels, limit)

            for label in labels:
                for edge in range(surface.polygon(label).num_edges()):
                    cross = surface.opposite_edge(label, edge)

                    if cross is None:
                        continue

                    # We do not call self.edge_matrix() since the surface might
                    # have overriden this (just returning the identity matrix e.g.)
                    # and we want to deduce the matrix from the attached polygon
                    # edges instead.
                    from flatsurf.geometry.categories import SimilaritySurfaces
                    matrix = SimilaritySurfaces.Oriented.ParentMethods.edge_matrix.f(surface, label, edge)

                    if not matrix.is_diagonal():
                        return False

                    if positive:
                        if matrix[0][0] < 0 or matrix[1][1] < 0:
                            return False

            return True

        def apply_matrix(self, m, in_place=True, mapping=False):
            r"""
            Carry out the GL(2,R) action of m on this surface and return the result.

            If in_place=True, then this is done in place and changes the surface.
            This can only be carried out if the surface is finite and mutable.

            If mapping=True, then we return a GL2RMapping between this surface and its image.
            In this case in_place must be False.

            If in_place=False, then a copy is made before the deformation.

            TESTS::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: T = S.apply_matrix(matrix([[1, 0], [0, 1]]), in_place=False)

            """
            if mapping is True:
                if in_place:
                    raise NotImplementedError(
                        "can not modify in place and return a mapping"
                    )
                from flatsurf.geometry.dilation_surface import GL2RMapping
                return GL2RMapping(self, m)
            if not in_place:
                if self.is_finite():
                    from sage.structure.element import get_coercion_model

                    cm = get_coercion_model()
                    field = cm.common_parent(self.base_ring(), m.base_ring())
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    s = MutableOrientedSimilaritySurface.from_surface(self.change_ring(field))
                    return s.apply_matrix(m)
                else:
                    return m * self
            else:
                # Make sure m is in the right state
                from sage.matrix.constructor import Matrix

                m = Matrix(self.base_ring(), 2, 2, m)
                if m.det() == self.base_ring().zero():
                    raise ValueError("can not deform by degenerate matrix")
                if not self.is_finite():
                    raise NotImplementedError(
                        "in-place GL(2,R) action only works for finite surfaces"
                    )
                us = self.underlying_surface()
                if not us.is_mutable():
                    raise ValueError("in-place changes only work for mutable surfaces")
                for label in self.labels():
                    us.change_polygon(label, m * self.polygon(label))
                if m.det() < self.base_ring().zero():
                    # Polygons were all reversed orientation. Need to redo gluings.

                    # First pass record new gluings in a dictionary.
                    new_glue = {}
                    seen_labels = set()
                    for p1 in self.labels():
                        n1 = self.polygon(p1).num_edges()
                        for e1 in range(n1):
                            p2, e2 = self.opposite_edge(p1, e1)
                            n2 = self.polygon(p2).num_edges()
                            if p2 in seen_labels:
                                pass
                            elif p1 == p2 and e1 > e2:
                                pass
                            else:
                                new_glue[(p1, n1 - 1 - e1)] = (p2, n2 - 1 - e2)
                        seen_labels.add(p1)
                    # Second pass: reassign gluings
                    for (p1, e1), (p2, e2) in new_glue.items():
                        us.change_edge_gluing(p1, e1, p2, e2)
                return self

        def _edge_needs_flip_Linfinity(self, p1, e1, p2, e2):
            r"""
            Check whether the provided edge which bounds two triangles should be flipped
            to get closer to the L-infinity Delaunay decomposition.

            TESTS::

                sage: from flatsurf import *
                sage: s = Surface_list(base_ring=QQ)
                sage: t1 = polygons((1,0),(-1,1),(0,-1))
                sage: t2 = polygons((0,1),(-1,0),(1,-1))
                sage: s.add_polygon(polygons(vertices=[(0,0), (1,0), (0,1)]))
                0
                sage: s.add_polygon(polygons(vertices=[(1,1), (0,1), (1,0)]))
                1
                sage: s.change_polygon_gluings(0, [(1,0), (1,1), (1,2)])
                sage: s.set_immutable()
                sage: [s._edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [False, False, False]

                sage: ss = matrix(2, [1,1,0,1]) * s
                sage: [ss._edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [False, False, False]
                sage: ss = matrix(2, [1,0,1,1]) * s
                sage: [ss._edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [False, False, False]

                sage: ss = matrix(2, [1,2,0,1]) * s
                sage: [ss._edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [False, False, True]

                sage: ss = matrix(2, [1,0,2,1]) * s
                sage: [ss._edge_needs_flip_Linfinity(0, i, 1, i) for i in range(3)]
                [True, False, False]
            """
            # safety check for now
            assert self.opposite_edge(p1, e1) == (p2, e2), "not opposite edges"

            # triangles
            poly1 = self.polygon(p1)
            poly2 = self.polygon(p2)
            if poly1.num_edges() != 3 or poly2.num_edges() != 3:
                raise ValueError("edge must be adjacent to two triangles")

            edge1 = poly1.edge(e1)
            edge1L = poly1.edge(e1 - 1)
            edge1R = poly1.edge(e1 + 1)
            edge2 = poly2.edge(e2)
            edge2L = poly2.edge(e2 - 1)
            edge2R = poly2.edge(e2 + 1)

            sim = self.edge_transformation(p2, e2)
            m = sim.derivative()  # matrix carrying p2 to p1
            if not m.is_one():
                edge2 = m * edge2
                edge2L = m * edge2L
                edge2R = m * edge2R

            # convexity check of the quadrilateral
            from flatsurf.geometry.polygon import wedge_product

            if wedge_product(edge2L, edge1R) <= 0 or wedge_product(edge1L, edge2R) <= 0:
                return False

            # compare the norms
            new_edge = edge2L + edge1R
            n1 = max(abs(edge1[0]), abs(edge1[1]))
            n = max(abs(new_edge[0]), abs(new_edge[1]))
            return n < n1

        def l_infinity_delaunay_triangulation(
            self, triangulated=None, in_place=None, limit=None, direction=None
        ):
            r"""
            Returns a L-infinity Delaunay triangulation of a surface, or make some
            triangle flips to get closer to the Delaunay decomposition.

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

                sage: from flatsurf import *
                sage: s0 = translation_surfaces.veech_double_n_gon(5)
                sage: field = s0.base_ring()
                sage: a = field.gen()
                sage: m = matrix(field, 2, [2,a,1,1])

                sage: s = m*s0
                sage: s = s.l_infinity_delaunay_triangulation()
                sage: TestSuite(s).run()

                sage: s = (m**2)*s0
                sage: s = s.l_infinity_delaunay_triangulation()
                sage: TestSuite(s).run()

                sage: s = (m**3)*s0
                sage: s = s.l_infinity_delaunay_triangulation()
                sage: TestSuite(s).run()

            TESTS:

            Verify that deprecated keywords do not cause errors::

                sage: s.l_infinity_delaunay_triangulation(triangulated=True)
                doctest:warning
                ...
                UserWarning: The triangulated keyword of l_infinity_delaunay_triangulation() has been deprecated and will be removed from a future version of sage-flatsurf. The keyword has no effect anymore.
                Surface built from 6 polygons
                sage: s.l_infinity_delaunay_triangulation(triangulated=False)
                Surface built from 6 polygons

            ::

                sage: s.l_infinity_delaunay_triangulation(in_place=True)
                Traceback (most recent call last):
                ...
                NotImplementedError: The in_place keyword for l_infinity_delaunay_triangulation() is not supported anymore. It did not work correctly in previous versions of sage-flatsurf and will be fully removed in a future version of sage-flatsurf.

            """
            if not self.is_finite():
                raise NotImplementedError(
                    "no L-infinity Delaunay implemented for infinite surfaces"
                )

            if triangulated is not None:
                import warnings
                warnings.warn("The triangulated keyword of l_infinity_delaunay_triangulation() has been deprecated and will be removed from a future version of sage-flatsurf. The keyword has no effect anymore.")
            if in_place is not None:
                raise NotImplementedError("The in_place keyword for l_infinity_delaunay_triangulation() is not supported anymore. It did not work correctly in previous versions of sage-flatsurf and will be fully removed in a future version of sage-flatsurf.")

            self = self.triangulate()

            from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
            self = MutableOrientedSimilaritySurface.from_surface(self)

            if direction is None:
                base_ring = self.base_ring()
                direction = self.vector_space()((base_ring.zero(), base_ring.one()))

            if direction.is_zero():
                raise ValueError("direction must be non-zero")

            triangles = set(self.labels())
            if limit is None:
                limit = -1
            else:
                limit = int(limit)
            while triangles and limit:
                p1 = triangles.pop()
                for e1 in range(3):
                    p2, e2 = self.opposite_edge(p1, e1)
                    if self._edge_needs_flip_Linfinity(p1, e1, p2, e2):
                        self.triangle_flip(p1, e1, in_place=True, direction=direction)
                        triangles.add(p1)
                        triangles.add(p2)
                        limit -= 1
            return self

        def _test_dilation_surface(self, **options):
            r"""
            Verify that this is a dilation surface.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S._test_dilation_surface()

            """
            tester = self._tester(**options)

            limit = None

            if not self.is_finite():
                limit = 32

            tester.assertTrue(DilationSurfaces.ParentMethods._is_dilation_surface(self, positive=False, limit=limit))


all_axioms += ('Positive',)
