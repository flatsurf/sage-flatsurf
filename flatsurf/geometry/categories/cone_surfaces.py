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

                def orbifold_points(self):
                    r"""
                    Return the orbifold points of this surface with their
                    orders.

                    An orbifold point of order `m \ge 2` is a vertex not on
                    the boundary whose total angle is `1/m` (in multiples of
                    `2\pi`), or the midpoint of an edge glued to itself, of
                    order 2. Every one of them has coordinates over the base
                    ring: a hyperbolic surface refuses the self-gluing of an
                    edge whose half-turn, and so whose midpoint, is not over
                    the base ring.

                    OUTPUT:

                    A list of pairs ``(p, m)``, ``p`` a point of this surface
                    and ``m`` an integer: first the vertices, in the order of
                    their first corner in
                    :meth:`~flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.edges`,
                    then the midpoints of the self-glued edges, in the same
                    order. A vertex whose angle is not a rational multiple of
                    `2\pi` is no orbifold point.

                    EXAMPLES:

                    The quotient of the square by its quarter turn, glued
                    from its half, a sphere with two points of order 4 and
                    one of order 2::

                        sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                        sage: S = MutableOrientedSimilaritySurface(QQ)
                        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1)]))
                        0
                        sage: S.glue((0, 2), (0, 2))
                        sage: S.glue((0, 0), (0, 1))
                        sage: S.set_immutable()
                        sage: S.orbifold_points()
                        [(Vertex 0 of polygon 0, 4),
                         (Vertex 1 of polygon 0, 4),
                         (Point (1/2, 1/2) of polygon 0, 2)]

                    A vertex on the boundary is no orbifold point::

                        sage: S = MutableOrientedSimilaritySurface(QQ)
                        sage: S.add_polygon(Polygon(vertices=[(1/2, 1/2), (0, 0), (1, 0)]))
                        0
                        sage: S.glue((0, 0), (0, 2))
                        sage: S.set_immutable()
                        sage: S.orbifold_points()
                        [(Vertex 0 of polygon 0, 4)]

                    A polygon glued to itself along each of its edges::

                        sage: from flatsurf import similarity_surfaces
                        sage: S = similarity_surfaces.self_glued_polygon(Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]))
                        sage: S.orbifold_points()
                        [(Point (1, 0) of polygon 0, 2),
                         (Point (3/2, 2) of polygon 0, 2),
                         (Point (1/2, 9/2) of polygon 0, 2),
                         (Point (0, 5/2) of polygon 0, 2)]

                    The regular octagon has a single vertex, of angle `3`::

                        sage: from flatsurf import translation_surfaces
                        sage: translation_surfaces.regular_octagon().orbifold_points()
                        []

                    A hyperbolic orbifold of signature `(0; 2, 8, 8)`, the
                    quotient of the regular octagon of angles `\pi/4` by its
                    rotation of order 8, its third side glued to itself::

                        sage: from flatsurf import HyperbolicPlane
                        sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                        sage: K.<r> = NumberField(x^4 + 32*x^2 - 32, embedding=AA(12*AA(2).sqrt() - 16).sqrt())
                        sage: s2 = r^2/12 + 4/3
                        sage: P = HyperbolicPlane(K).polygon(vertices=[(0, 0), (r, 0), (r / s2, r / s2)])
                        sage: S = MutableOrientedHyperbolicSurface(K)
                        sage: S.add_polygon(P)
                        0
                        sage: S.glue((0, 0), (0, 2))
                        sage: S.glue((0, 1), (0, 1))
                        sage: S.set_immutable()
                        sage: S.orbifold_points()
                        [(Vertex 0 of polygon 0, 8),
                         (Vertex 1 of polygon 0, 8),
                         (Point (1/48*r^3 + 5/6*r, 1/48*r^3 + 1/3*r) of polygon 0, 2)]

                    A vertex of the Lambert quadrilateral has an angle that is
                    not a rational multiple of `2\pi`::

                        sage: S = MutableOrientedHyperbolicSurface(QQ)
                        sage: S.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                        0
                        sage: S.glue((0, 0), (0, 3))
                        sage: S.glue((0, 1), (0, 2))
                        sage: S.set_immutable()
                        sage: S.orbifold_points()
                        [(Vertex 0 of polygon 0, 4), (Vertex 1 of polygon 0, 2)]

                    TESTS:

                    A flat vertex whose angle is not a rational multiple of
                    `2\pi`::

                        sage: S = similarity_surfaces.pillowcase(Polygon(vertices=[(0, 0), (1, 0), (1/3, 1/2)]))
                        sage: S.orbifold_points()
                        []

                    """
                    from sage.all import ZZ

                    boundary = set()
                    for label, edge in self.edges():
                        if self.opposite_edge(label, edge) is None:
                            n = len(self.polygon(label).vertices())
                            # pylint: disable-next=not-callable
                            boundary.add(self(label, edge))
                            # pylint: disable-next=not-callable
                            boundary.add(self(label, (edge + 1) % n))

                    points = []

                    seen = set()
                    for label, corner in self.edges():
                        vertex = self(label, corner)  # pylint: disable=not-callable
                        if vertex in seen:
                            continue
                        seen.add(vertex)
                        if vertex in boundary:
                            continue
                        order = self._orbifold_order(vertex)
                        if order is not None and order >= 2:
                            points.append((vertex, order))

                    for label, edge in self.edges():
                        if self.opposite_edge(label, edge) == (label, edge):
                            midpoint = self._edge_midpoint(label, edge)
                            # pylint: disable-next=not-callable
                            points.append((self(label, midpoint), ZZ(2)))

                    return points

                def _orbifold_order(self, vertex):
                    r"""
                    Return the positive integer `m` such that the angle of
                    ``vertex`` is `1/m` (in multiples of `2\pi`), or ``None``
                    if there is none.

                    This is a helper method for :meth:`orbifold_points` and
                    :meth:`is_orbifold`; a regular vertex has order 1.

                    EXAMPLES::

                        sage: from flatsurf import polygons, similarity_surfaces
                        sage: S = similarity_surfaces.pillowcase(polygons.square())
                        sage: S._orbifold_order(S(0, 0))
                        2

                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.square_torus()
                        sage: S._orbifold_order(S(0, 0))
                        1
                        sage: S = translation_surfaces.regular_octagon()
                        sage: S._orbifold_order(S(0, 0)) is None
                        True

                    """
                    from sage.all import QQ, ZZ

                    try:
                        angle = QQ(vertex.angle())
                    except (ValueError, NotImplementedError):
                        # The exact angle of a hyperbolic vertex raises a
                        # ValueError, the one of a flat vertex a
                        # NotImplementedError, when it is not a rational
                        # multiple of 2π.
                        return None

                    if angle.numerator() == 1:
                        return ZZ(angle.denominator())
                    return None

                def is_orbifold(self):
                    r"""
                    Return whether this surface is an orbifold, i.e., whether
                    the angle of each of its vertices is `1` or `1/m` (in
                    multiples of `2\pi`) for an integer `m \ge 2`.

                    The midpoints of the edges glued to themselves are
                    orbifold points of order 2 and always qualify.

                    EXAMPLES::

                        sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                        sage: S = MutableOrientedSimilaritySurface(QQ)
                        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1)]))
                        0
                        sage: S.glue((0, 2), (0, 2))
                        sage: S.glue((0, 0), (0, 1))
                        sage: S.set_immutable()
                        sage: S.is_orbifold()
                        True

                        sage: from flatsurf import similarity_surfaces
                        sage: similarity_surfaces.self_glued_polygon(Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])).is_orbifold()
                        True

                        sage: from flatsurf import translation_surfaces
                        sage: translation_surfaces.regular_octagon().is_orbifold()
                        False

                    A hyperbolic sphere with a cone point of angle neither
                    `1` nor `1/m`::

                        sage: from flatsurf import HyperbolicPlane
                        sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                        sage: S = MutableOrientedHyperbolicSurface(QQ)
                        sage: S.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                        0
                        sage: S.glue((0, 0), (0, 3))
                        sage: S.glue((0, 1), (0, 2))
                        sage: S.set_immutable()
                        sage: S.is_orbifold()
                        False

                    Orbifolds with boundary are not supported::

                        sage: S = MutableOrientedSimilaritySurface(QQ)
                        sage: S.add_polygon(Polygon(vertices=[(1/2, 1/2), (0, 0), (1, 0)]))
                        0
                        sage: S.glue((0, 0), (0, 2))
                        sage: S.set_immutable()
                        sage: S.is_orbifold()
                        Traceback (most recent call last):
                        ...
                        NotImplementedError: orbifolds with boundary are not supported yet

                    """
                    if self.is_with_boundary():
                        raise NotImplementedError(
                            "orbifolds with boundary are not supported yet"
                        )

                    return all(
                        self._orbifold_order(vertex) is not None
                        for vertex in self.vertices()
                    )

            class Connected(SurfaceCategoryWithAxiom):
                r"""
                The category of oriented connected cone surfaces built from
                finitely many polygons.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: from flatsurf.geometry.categories import ConeSurfaces
                    sage: translation_surfaces.square_torus() in ConeSurfaces().FiniteType().Oriented().Connected()
                    True

                """

                class ParentMethods:
                    r"""
                    Provides methods available to all oriented connected cone
                    surfaces built from finitely many polygons.

                    If you want to add functionality for such surfaces you
                    most likely want to put it here.
                    """

                    def signature(self):
                        r"""
                        Return the signature of this surface: its genus, the
                        orders of its cone points, its number of cusps and the
                        traces of the holonomy of its boundary components.

                        OUTPUT:

                        A :class:`~flatsurf.geometry.signature.SurfaceSignature`.
                        Its cone points are the vertices not on the boundary
                        whose angle `2\pi (1 + k)` is not `2\pi`, of order
                        `k`, and the midpoints of the edges glued to
                        themselves, of order `-1/2`. Its boundary traces are
                        the traces in `\mathrm{SO}(2, 1)` of the holonomy of
                        the
                        :meth:`~flatsurf.geometry.fundamental_group.FundamentalGroup.boundary_loops`
                        on a hyperbolic surface, ``None`` on a flat one.

                        A ``NotImplementedError`` is raised if the angle of an
                        interior vertex is not a rational multiple of `2\pi`;
                        a ``ValueError`` if the angle of a vertex on the
                        boundary is not `\pi`, i.e., the boundary is not
                        geodesic.

                        ALGORITHM:

                        The genus `g` is read from the Euler characteristic
                        `2 - 2g - b` of the surface, where `b` is the number
                        of boundary components; it counts the midpoint of an
                        edge glued to itself as a vertex.

                        EXAMPLES:

                        A translation surface in `H(2)`::

                            sage: from flatsurf import translation_surfaces
                            sage: translation_surfaces.regular_octagon().signature()
                            SurfaceSignature(genus=2, cone_orders=(2,))

                        A polygon glued to itself along each of its edges is
                        a sphere with four cone points of angle `\pi`, the
                        midpoints of the edges::

                            sage: from flatsurf import Polygon, similarity_surfaces
                            sage: S = similarity_surfaces.self_glued_polygon(Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)]))
                            sage: S.signature()
                            SurfaceSignature(genus=0, cone_orders=(-1/2, -1/2, -1/2, -1/2))

                        A flat cylinder has two boundary components without a
                        trace::

                            sage: from flatsurf import MutableOrientedSimilaritySurface
                            sage: S = MutableOrientedSimilaritySurface(QQ)
                            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
                            0
                            sage: S.glue((0, 0), (0, 2))
                            sage: S.set_immutable()
                            sage: S.signature()
                            SurfaceSignature(genus=0, cone_orders=(), boundary_traces=(None, None))

                        A pair of pants, glued from two right-angled regular
                        hexagons: its three boundary components have length
                        `\operatorname{arccosh} 7`, of trace `1 + 2 \cdot 7`::

                            sage: from flatsurf import HyperbolicPlane
                            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
                            sage: rho = AA(6).sqrt() / 3
                            sage: h = AA(2).sqrt() / 2
                            sage: P = HyperbolicPlane(AA).polygon(vertices=[(rho, 0), (rho/2, h), (-rho/2, h), (-rho, 0), (-rho/2, -h), (rho/2, -h)])
                            sage: S = MutableOrientedHyperbolicSurface(AA)
                            sage: S.add_polygon(P)
                            0
                            sage: S.add_polygon(P)
                            1
                            sage: S.glue((0, 0), (1, 0))
                            sage: S.glue((0, 2), (1, 4))
                            sage: S.glue((0, 4), (1, 2))
                            sage: S.set_immutable(minimal_field=True)
                            sage: s = S.signature()
                            sage: s
                            SurfaceSignature(genus=0, cone_orders=(), boundary_traces=(15, 15, 15))
                            sage: s.area() == S.area()
                            True

                        A hyperbolic orbifold of signature `(0; 2, 8, 8)`, the
                        quotient of the regular octagon of angles `\pi/4` by
                        its rotation of order 8::

                            sage: K.<r> = NumberField(x^4 + 32*x^2 - 32, embedding=AA(12*AA(2).sqrt() - 16).sqrt())
                            sage: s2 = r^2/12 + 4/3
                            sage: P = HyperbolicPlane(K).polygon(vertices=[(0, 0), (r, 0), (r / s2, r / s2)])
                            sage: S = MutableOrientedHyperbolicSurface(K)
                            sage: S.add_polygon(P)
                            0
                            sage: S.glue((0, 0), (0, 2))
                            sage: S.glue((0, 1), (0, 1))
                            sage: S.set_immutable()
                            sage: S.signature()
                            SurfaceSignature(genus=0, cone_orders=(-7/8, -7/8, -1/2))
                            sage: S.signature().katok()
                            (0, (2, 8, 8), 0)

                        A vertex of the Lambert quadrilateral has an angle that
                        is not a rational multiple of `2\pi`::

                            sage: S = MutableOrientedHyperbolicSurface(QQ)
                            sage: S.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
                            0
                            sage: S.glue((0, 0), (0, 3))
                            sage: S.glue((0, 1), (0, 2))
                            sage: S.set_immutable()
                            sage: S.signature()
                            Traceback (most recent call last):
                            ...
                            NotImplementedError: the angle of Vertex 2 of polygon 0 is not a rational multiple of 2π

                        A corner on the boundary::

                            sage: S = MutableOrientedSimilaritySurface(QQ)
                            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1)]))
                            0
                            sage: S.set_immutable()
                            sage: S.signature()
                            Traceback (most recent call last):
                            ...
                            ValueError: the boundary of the surface is not geodesic at Vertex 0 of polygon 0, of angle 1/8

                        """
                        from sage.all import QQ, ZZ

                        from flatsurf.geometry.signature import SurfaceSignature

                        components = self.boundary_components()

                        # Each vertex on the boundary is the start of an unglued
                        # edge.
                        boundary = {
                            self(label, edge)  # pylint: disable=not-callable
                            for component in components
                            for label, edge in component
                        }

                        orders = []
                        seen = set()
                        for label, corner in self.edges():
                            vertex = self(label, corner)  # pylint: disable=not-callable
                            if vertex in seen:
                                continue
                            seen.add(vertex)

                            try:
                                angle = QQ(vertex.angle())
                            except (ValueError, NotImplementedError) as error:
                                # The exact angle of a hyperbolic vertex raises
                                # a ValueError, the one of a flat vertex a
                                # NotImplementedError, when it is not a
                                # rational multiple of 2π.
                                if vertex in boundary:
                                    raise ValueError(
                                        f"the boundary of the surface is not geodesic at {vertex}"
                                    ) from error
                                raise NotImplementedError(
                                    f"the angle of {vertex} is not a rational multiple of 2π"
                                ) from error

                            if vertex in boundary:
                                if angle != QQ((1, 2)):
                                    raise ValueError(
                                        f"the boundary of the surface is not geodesic at {vertex}, of angle {angle}"
                                    )
                            elif angle != 1:
                                orders.append(angle - 1)

                        for label, edge in self.edges():
                            if self.opposite_edge(label, edge) == (label, edge):
                                orders.append(-QQ((1, 2)))

                        genus = ZZ(
                            (2 - len(components) - self.euler_characteristic()) / 2
                        )

                        loops = self.fundamental_group(
                            punctured_vertices=True
                        ).boundary_loops()

                        return SurfaceSignature(
                            genus,
                            orders,
                            boundary_traces=[
                                self._boundary_trace(loop) for loop in loops
                            ],
                        )

                    def _test_gauss_bonnet(self, **options):
                        r"""
                        Verify the Gauss-Bonnet theorem for the
                        :meth:`signature` of this surface: the Euler
                        characteristic `\chi(S, \beta)` is zero on a flat
                        surface, and minus the area, in units of `2\pi`, on a
                        hyperbolic one.

                        The test is skipped when the signature is not
                        defined.

                        EXAMPLES::

                            sage: from flatsurf import translation_surfaces
                            sage: translation_surfaces.regular_octagon()._test_gauss_bonnet()

                        """
                        tester = self._tester(**options)

                        try:
                            signature = self.signature()
                        except (NotImplementedError, ValueError):
                            return

                        chi = signature.euler_characteristic()

                        from flatsurf.geometry.categories import (
                            FlatConeSurfaces,
                            HyperbolicConeSurfaces,
                        )

                        if self in HyperbolicConeSurfaces():
                            if all(
                                polygon.is_rational() for polygon in self.polygons()
                            ):
                                tester.assertEqual(self.area(), -chi)
                            else:
                                # The area of a polygon whose corners are not
                                # rational multiples of π is only known
                                # numerically.
                                tester.assertAlmostEqual(
                                    sum(
                                        polygon.area(numerical=True)
                                        for polygon in self.polygons()
                                    ),
                                    float(-chi),
                                )
                        elif self in FlatConeSurfaces():
                            tester.assertEqual(chi, 0)


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

                def _boundary_trace(self, loop):
                    r"""
                    Return the invariant of the holonomy of the boundary
                    ``loop`` recorded by the signature, ``None`` for a flat
                    surface.

                    This is a helper method for
                    :meth:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces.FiniteType.Oriented.Connected.ParentMethods.signature`.
                    The holonomy of a loop parallel to a flat boundary is a
                    translation or a rotation, which carries no invariant of
                    the length of the boundary.

                    INPUT:

                    - ``loop`` -- an element of the fundamental group of this
                      surface

                    EXAMPLES::

                        sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
                        sage: S = MutableOrientedSimilaritySurface(QQ)
                        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
                        0
                        sage: S.glue((0, 0), (0, 2))
                        sage: S.set_immutable()
                        sage: loop = S.fundamental_group().boundary_loops()[0]
                        sage: S._boundary_trace(loop) is None
                        True

                    """
                    return None

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
