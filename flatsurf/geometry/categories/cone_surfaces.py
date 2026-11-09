r"""
The categories of cone surfaces and of flat cone surfaces.

A cone surface is a surface with a metric of constant curvature that has cone
points, in any of the geometries (Euclidean or hyperbolic). The category
:class:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces` is what flat
and hyperbolic cone surfaces have in common, e.g., that a surface built from
finitely many polygons has an area.

A flat cone surface is a surface that can be built by gluing Euclidean polygons
along their edges such that the matrix describing `monodromy
<https://en.wikipedia.org/wiki/(G,X)-manifold#Monodromy>`_ along a closed path
is an isometry; that matrix is given by multiplying the individual matrices
that describe how to transition between pairs of glued edges, see
:meth:`~.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.edge_matrix`.

In sage-flatsurf, we restrict flat cone surfaces slightly by requiring that a
flat cone surface is given by polygons such that each edge matrix is an
isometry. This is the category
:class:`~flatsurf.geometry.categories.cone_surfaces.FlatConeSurfaces`.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES:

We glue the sides of a square with a rotation of π/2. Since each gluing is just
a rotation, this is a flat cone surface::

    sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = MutableOrientedSimilaritySurface(QQ)
    sage: S.add_polygon(P, label=0)
    0
    sage: S.glue((0, 0), (0, 1))
    sage: S.glue((0, 2), (0, 3))
    sage: S.set_immutable()

    sage: C = S.category()

    sage: from flatsurf.geometry.categories import FlatConeSurfaces
    sage: C.is_subcategory(FlatConeSurfaces())
    True

"""

# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2013-2019 Vincent Delecroix
#                      2013-2019 W. Patrick Hooper
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)


class ConeSurfaces(SurfaceCategory):
    r"""
    The category of surfaces with a metric of constant curvature that has cone
    points, in any geometry.

    This is the category that flat and hyperbolic cone surfaces have in common.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import ConeSurfaces
        sage: ConeSurfaces()
        Category of cone surfaces
        sage: ConeSurfaces().super_categories()
        [Category of polygonal surfaces]

    The category of flat cone surfaces is a subcategory::

        sage: from flatsurf.geometry.categories import FlatConeSurfaces
        sage: FlatConeSurfaces().is_subcategory(ConeSurfaces())
        True

    """

    def super_categories(self):
        r"""
        Return the categories that a cone surface is also always a member of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import ConeSurfaces
            sage: ConeSurfaces().super_categories()
            [Category of polygonal surfaces]

        """
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces

        return [PolygonalSurfaces()]

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The category of cone surfaces built from finitely many polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import ConeSurfaces
            sage: ConeSurfaces().FiniteType()
            Category of finite type cone surfaces

        """

        class ParentMethods:
            r"""
            Provides methods available to all cone surfaces built from finitely
            many polygons.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def area(self):
                r"""
                Return the area of this surface.

                The unit depends on the geometry: for flat polygons this is an
                element of the base ring; for hyperbolic polygons, the area is
                measured in units of 2π, so that for a hyperbolic surface it
                is `-\chi(S, \beta)`.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.area()
                    3

                """
                return sum(p.area() for p in self.polygons())

        class Oriented(SurfaceCategoryWithAxiom):
            r"""
            The category of oriented cone surfaces built from finitely many
            polygons.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: from flatsurf.geometry.categories import ConeSurfaces
                sage: translation_surfaces.square_torus() in ConeSurfaces().FiniteType().Oriented()
                True

            """

            class ParentMethods:
                r"""
                Provides methods available to all oriented cone surfaces built
                from finitely many polygons.

                If you want to add functionality for such surfaces you most
                likely want to put it here.
                """

                def holonomy_representation(
                    self, base_label=None, punctured_vertices=True
                ):
                    r"""
                    Return the holonomy representation of this surface.

                    It sends a closed walk crossing the edges from polygon to
                    polygon to the isometry obtained by developing the
                    polygons along the walk, see
                    :mod:`~flatsurf.geometry.holonomy`.

                    INPUT:

                    - ``base_label`` -- a polygon label or ``None`` (default:
                      ``None``); the polygon at which the walks start and
                      end, the root of the surface if ``None``

                    - ``punctured_vertices`` -- a boolean or a list of
                      vertices of this surface (default: ``True``); the
                      vertices at which the surface is punctured, ``True``
                      standing for all of them. Every other interior vertex
                      must have an angle that is an integer multiple of `2\pi`

                    OUTPUT:

                    A :class:`~flatsurf.geometry.holonomy.HolonomyRepresentation`
                    from
                    :meth:`~flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.FiniteType.Oriented.ParentMethods.fundamental_group`
                    to the group of the ``edge_transformation`` of this
                    surface: similarities for a flat surface, isometries of
                    the hyperbolic plane for a hyperbolic one.

                    EXAMPLES:

                    On a translation surface, the holonomy consists of
                    translations::

                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.square_torus()
                        sage: phi = S.holonomy_representation()
                        sage: [phi(g) for g in phi.domain().gens()]
                        [(x, y) |-> (x, y - 1), (x, y) |-> (x + 1, y)]

                    The vertex of the torus is regular, so it can be filled::

                        sage: S.holonomy_representation(punctured_vertices=False)
                        Holonomy morphism:
                          From: Fundamental group of Translation Surface in H_1(0) built from a square based at polygon 0
                          To:   Similarity group over Rational Field

                    On a hyperbolic surface, the holonomy is a group of
                    isometries of the hyperbolic plane::

                        sage: from flatsurf import HyperbolicPlane
                        sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                        sage: H = HyperbolicPlane(QQ)
                        sage: S = MutableOrientedHyperbolicSurface(QQ)
                        sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                        0
                        sage: S.glue((0, 0), (0, 3))
                        sage: S.glue((0, 1), (0, 2))
                        sage: S.set_immutable()
                        sage: S.holonomy_representation().image()
                        Matrix group over Rational Field with 2 generators (
                        [ 0  1  0]  [   0 -5/3  4/3]
                        [-1  0  0]  [   1    0    0]
                        [ 0  0  1], [   0 -4/3  5/3]
                        )

                    A vertex whose angle is not an integer multiple of `2\pi`
                    cannot be filled::

                        sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                        sage: S = MutableOrientedSimilaritySurface(QQ)
                        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
                        0
                        sage: S.glue((0, 0), (0, 1))
                        sage: S.glue((0, 2), (0, 3))
                        sage: S.set_immutable()
                        sage: S.holonomy_representation(punctured_vertices=False)
                        Traceback (most recent call last):
                        ...
                        ValueError: cannot fill the vertex Vertex 0 of polygon 0 whose angle 1/2 is not an integer multiple of 2π

                    """
                    from flatsurf.geometry.holonomy import holonomy_representation

                    return holonomy_representation(
                        self,
                        base_label=base_label,
                        punctured_vertices=punctured_vertices,
                    )


class FlatConeSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    isometries on the edges.

    See :mod:`~flatsurf.geometry.categories.cone_surfaces` and
    :meth:`~ParentMethods.is_cone_surface` on how this differs slightly from
    the customary definition of a cone surface.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import FlatConeSurfaces
        sage: FlatConeSurfaces()
        Category of flat cone surfaces

    """

    def super_categories(self):
        r"""
        Return the categories that a flat cone surface is also always a member
        of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import FlatConeSurfaces
            sage: FlatConeSurfaces().super_categories()
            [Category of similarity surfaces, Category of cone surfaces]

        """
        from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces

        return [SimilaritySurfaces(), ConeSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all cone surfaces.

        If you want to add functionality for such surfaces you most likely want
        to put it here.
        """

        def is_cone_surface(self):
            r"""
            Return whether this surface is a cone surface, i.e., whether its
            edges are glued by isometries.

            .. NOTE::

                This is a stronger requirement than the usual definition of a
                cone surface, see
                :mod:`~flatsurf.geometry.categories.cone_surfaces` for details.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.is_cone_surface()
                True

            """
            return True

        @staticmethod
        def _is_cone_surface(surface, limit=None):
            r"""
            Return whether ``surface`` is a cone surface by checking how its
            polygons are glued.

            This is a helper method for :meth:`is_cone_surface` and
            :meth:`_test_cone_surface`.

            INPUT:

            - ``surface`` -- an oriented similarity surface

            - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
              the first ``limit`` polygons are checked

            EXAMPLES::

                sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
                sage: S = translation_surfaces.infinite_staircase()

                sage: from flatsurf.geometry.categories import FlatConeSurfaces
                sage: FlatConeSurfaces.ParentMethods._is_cone_surface(S, limit=8)
                doctest:warning
                ...
                UserWarning: limit has been deprecated as a keyword argument for _is_cone_surface() and will be removed from a future version of sage-flatsurf; ...
                True
                sage: FlatConeSurfaces.ParentMethods._is_cone_surface(MutableOrientedSimilaritySurface.from_surface(S, labels=S.labels()[:8]))
                True

            ::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

                sage: FlatConeSurfaces.ParentMethods._is_cone_surface(S)
                True

            """
            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for _is_cone_surface() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

            if "Oriented" not in surface.category().axioms():
                raise NotImplementedError(
                    "cannot check whether a non-oriented surface is a cone surface yet"
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

                    if matrix * matrix.transpose() != 1:
                        return False

            return True

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The category of cone surfaces built from finitely many polygons.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
            sage: S = similarity_surfaces.self_glued_polygon(P)

            sage: from flatsurf.geometry.categories import FlatConeSurfaces
            sage: S in FlatConeSurfaces().FiniteType()
            True

        """

        class Oriented(SurfaceCategoryWithAxiom):
            r"""
            The category of oriented cone surfaces, i.e., orientable cone surfaces
            whose orientation can be chosen to be compatible with the embedding of
            its polygons in the real plane.

            EXAMPLES::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

                sage: from flatsurf.geometry.categories import FlatConeSurfaces
                sage: S in FlatConeSurfaces().Oriented()
                True

            """

            class ParentMethods:
                r"""
                Provides methods available to all oriented cone surfaces.

                If you want to add functionality for such surfaces you most likely
                want to put it here.
                """

                def _test_cone_surface(self, **options):
                    r"""
                    Verify that this is a cone surface.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.square_torus()
                        sage: S.set_immutable()
                        sage: S._test_cone_surface()

                    """
                    tester = self._tester(**options)

                    limit = None

                    if not self.is_finite_type():
                        limit = 32

                    tester.assertTrue(
                        FlatConeSurfaces.ParentMethods._is_cone_surface(
                            self, limit=limit
                        )
                    )

            class WithoutBoundary(SurfaceCategoryWithAxiom):
                r"""
                The category of oriented cone surfaces without boundary.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)

                    sage: from flatsurf.geometry.categories import FlatConeSurfaces
                    sage: S in FlatConeSurfaces().Oriented().WithoutBoundary()
                    True

                """

                class Connected(SurfaceCategoryWithAxiom):
                    r"""
                    The category of oriented connected cone surfaces without boundary.

                    EXAMPLES::

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                        sage: S = similarity_surfaces.self_glued_polygon(P)

                        sage: from flatsurf.geometry.categories import FlatConeSurfaces
                        sage: S in FlatConeSurfaces().Oriented().Connected()
                        True

                    """

                    class ParentMethods:
                        r"""
                        Provides methods available to all oriented connected
                        cone surfaces without boundary.

                        If you want to add functionality for such surfaces you
                        most likely want to put it here.
                        """

                        def _test_genus(self, **options):
                            r"""
                            Verify that the genus is compatible with the angles of the
                            singularities.

                            ALGORITHM:

                            We use the angles around the vertices of the surface to compute the
                            genus, see e.g. [Massart2021] p.17 and compare this
                            to the genus computed directly from the polygon
                            gluings.

                            EXAMPLES::

                                sage: from flatsurf import translation_surfaces
                                sage: translation_surfaces.octagon_and_squares()._test_genus()

                            .. [Massart2021] D. Massart. A short introduction to translation
                                surfaces, Veech surfaces, and Teichműller dynamics.
                                https://hal.science/hal-03300179/document

                            """
                            tester = self._tester(**options)

                            for edge in self.edges():
                                if self.opposite_edge(*edge) == edge:
                                    # The genus formula below is wrong when
                                    # there is a non-materialized vertex on an
                                    # edge.
                                    return

                            tester.assertAlmostEqual(
                                self.genus(),
                                sum(a - 1 for a in self.angles(numerical=True)) / 2.0
                                + 1,
                            )
