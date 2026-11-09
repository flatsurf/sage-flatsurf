r"""
Surfaces built by gluing compact polygons of the hyperbolic plane.

The polygons are
:class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon` objects, whose
vertices are finite points of the hyperbolic plane. Two edges are glued by the
orientation preserving isometry that sends one onto the other one with reverse
orientation, which requires the two edges to have the same length; this is
checked when the surface is made immutable.

Points of the surface are given by a polygon label and a point of the
hyperbolic plane, or its coordinates in the Klein model.

EXAMPLES:

We glue the sides of a Lambert quadrilateral, which has three right angles,
to a sphere with three cone points::

    sage: from flatsurf import HyperbolicPlane
    sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
    sage: H = HyperbolicPlane(QQ)
    sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
    sage: S = MutableOrientedHyperbolicSurface(QQ)
    sage: S.add_polygon(P)
    0
    sage: S.glue((0, 0), (0, 3))
    sage: S.glue((0, 1), (0, 2))
    sage: S.set_immutable()
    sage: S
    Genus 0 Hyperbolic Surface built from a quadrilateral

    sage: S.category()
    Category of connected without boundary finite type oriented hyperbolic cone surfaces

Points are given by a polygon label and Klein coordinates, a point of the
hyperbolic plane, or the index of a vertex of the polygon; they print with
their Klein coordinates::

    sage: S(0, (1/4, 1/4))
    Point (1/4, 1/4) of polygon 0
    sage: S(0, (1/4, 1/4)).coordinates(0)
    ((1/4, 1/4),)
    sage: S(0, H.point(1/4, 1/4, model="klein")) == S(0, (1/4, 1/4))
    True
    sage: S(0, 0)
    Vertex 0 of polygon 0
    sage: S(0, 0).coordinates(0)
    (I,)
    sage: S(0, (1/2, 0)) == S(0, (0, 1/2))
    True
    sage: set(S(0, 1).coordinates(0)) == {P.vertex(1), P.vertex(3)}
    True
    sage: len(S.vertices())
    3

The gluings are isometries of the hyperbolic plane; the one of the edge 0
sends its start to the end of the edge 3, and its end to the start of the
edge 3::

    sage: g = S.edge_transformation(0, 0)
    sage: g
    [ 0 -1  0]
    [ 1  0  0]
    [ 0  0  1]
    sage: g(P.vertex(0)), g(P.vertex(1))
    (I, (0, 1/2))

The three vertices are cone points of angles π/2, π and about 0.39π; the
last one is not a rational multiple of π::

    sage: S(0, 0).angle()
    1/4
    sage: S(0, 1).angle()
    1/2
    sage: S(0, 2).angle()
    Traceback (most recent call last):
    ...
    ValueError: angle is not a rational multiple of 2π; set numerical=True to get a numerical approximation
    sage: S(0, 2).angle(numerical=True)  # abs tol 1e-12
    0.195913276015304
    sage: sorted(S.angles(numerical=True))  # abs tol 1e-12
    [0.195913276015304, 0.25, 0.5]

Glued edges must have the same length::

    sage: T = MutableOrientedHyperbolicSurface(QQ)
    sage: T.add_polygon(P)
    0
    sage: T.glue((0, 0), (0, 1))
    sage: T.set_immutable()
    Traceback (most recent call last):
    ...
    ValueError: edge 0 of polygon 0 and edge 1 of polygon 0 are glued but do not have the same length

An edge glued to itself is folded at its midpoint by a half turn, which is
not always defined over the base ring::

    sage: T.glue((0, 0), (0, 0))
    sage: T.edge_transformation(0, 0)
    Traceback (most recent call last):
    ...
    ValueError: isometry is not defined over Rational Field; adjoin a square root of 4/3

TESTS::

    sage: TestSuite(S).run()

"""

######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2026 Vincent Delecroix
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
######################################################################

from flatsurf.geometry.surface import MutablePolygonalSurface, Surface_base


class OrientedHyperbolicSurface(Surface_base):
    r"""
    Base class for surfaces built from compact polygons of the hyperbolic
    plane that are glued with orientation preserving isometries.

    INPUT:

    - ``base`` -- a ring, the base ring of the hyperbolic plane of the
      polygons

    - ``category`` -- a category or ``None`` (default: ``None``); refined to
      a subcategory of oriented hyperbolic cone surfaces

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface, OrientedHyperbolicSurface
        sage: S = MutableOrientedHyperbolicSurface(QQ)
        sage: isinstance(S, OrientedHyperbolicSurface)
        True

    TESTS::

        sage: MutableOrientedHyperbolicSurface(0)
        Traceback (most recent call last):
        ...
        TypeError: base ring must be a ring

    """

    def __init__(self, base, category=None):
        from sage.categories.all import Rings

        if base not in Rings():
            raise TypeError("base ring must be a ring")

        from flatsurf.geometry.categories import HyperbolicConeSurfaces

        if category is None:
            category = HyperbolicConeSurfaces().Oriented()

        category &= HyperbolicConeSurfaces().Oriented()

        super().__init__(base, category=category)


class MutableOrientedHyperbolicSurface(
    OrientedHyperbolicSurface, MutablePolygonalSurface
):
    r"""
    A surface built from compact polygons of the hyperbolic plane glued by
    orientation preserving isometries.

    INPUT:

    - ``base`` -- an ordered field, the base ring of the hyperbolic plane of
      the polygons

    - ``category`` -- a category or ``None`` (default: ``None``); refined to
      a subcategory of finite type oriented hyperbolic cone surfaces

    EXAMPLES:

    The regular octagon with angles π/4, with opposite sides glued, is a
    closed surface of genus 2 whose vertices form a single regular point::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
        sage: K.<r> = NumberField(x^4 + 32*x^2 - 32, embedding=AA(12*AA(2).sqrt() - 16).sqrt())
        sage: s2 = r^2/12 + 4/3
        sage: cos = [1, s2/2, 0, -s2/2, -1, -s2/2, 0, s2/2]
        sage: sin = [0, s2/2, 1, s2/2, 0, -s2/2, -1, -s2/2]
        sage: P = HyperbolicPlane(K).polygon(vertices=[(r * c, r * s) for c, s in zip(cos, sin)])
        sage: S = MutableOrientedHyperbolicSurface(K)
        sage: S.add_polygon(P)
        0
        sage: for i in range(4):
        ....:     S.glue((0, i), (0, i + 4))
        sage: S.set_immutable()
        sage: S
        Genus 2 Hyperbolic Surface built from a regular octagon
        sage: S.category()
        Category of connected without boundary finite type oriented hyperbolic cone surfaces

        sage: S.vertices()
        {Vertex 0 of polygon 0}
        sage: S(0, 0).angle()
        1
        sage: S.area()
        2

    A point in the interior of the polygon, and a point in the interior of
    the edge 0, which has a second representative on the edge 4::

        sage: S(0, (r/2, 0))
        Point (1/2*r, 0) of polygon 0
        sage: m = (vector(P.vertex(0).coordinates(model="klein")) + vector(P.vertex(1).coordinates(model="klein"))) / 2
        sage: p = S(0, m)
        sage: len(p.coordinates(0))
        2

    The gluing of the edge 0 to the edge 4 is a hyperbolic isometry; it
    sends the representative on the edge 0 to the one on the edge 4 and the
    vertex 0 to the vertex 5::

        sage: g = S.edge_transformation(0, 0)
        sage: g.trace()
        2/3*r^2 + 65/3
        sage: P.get_point_position(g * m)
        point positioned on interior of edge 4 of polygon
        sage: H = S.hyperbolic_plane()
        sage: H.point(*(g * m), model="klein") in p.coordinates(0)
        True
        sage: g * vector((r, 0)) == vector(P.vertex(5).coordinates(model="klein"))
        True

    TESTS::

        sage: TestSuite(S).run()

    """

    def __init__(self, base, category=None):
        from flatsurf.geometry.categories import HyperbolicConeSurfaces

        if category is None:
            category = HyperbolicConeSurfaces().Oriented().FiniteType()

        category &= HyperbolicConeSurfaces().Oriented().FiniteType()

        super().__init__(base, category=category)

    def add_polygon(self, polygon, *, label=None):
        r"""
        Add an unglued polygon to this surface and return its label.

        INPUT:

        - ``polygon`` -- a
          :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`

        - ``label`` -- a hashable identifier or ``None`` (default: ``None``);
          if ``None`` an integer identifier is automatically selected

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: H = HyperbolicPlane(QQ)
            sage: S = MutableOrientedHyperbolicSurface(AA)
            sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            0
            sage: S.polygon(0).base_ring()
            Algebraic Real Field

        Euclidean polygons cannot be added::

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            Traceback (most recent call last):
            ...
            TypeError: polygon must be a hyperbolic polygon given by its vertices

        """
        self._check_polygon(polygon)

        return super().add_polygon(polygon, label=label)

    def replace_polygon(self, label, polygon):
        r"""
        Replace the polygon ``label`` with ``polygon`` while keeping its
        gluings intact.

        INPUT:

        - ``label`` -- a polygon label of this surface

        - ``polygon`` -- a
          :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
          with the same number of vertices

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: H = HyperbolicPlane(QQ)
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            0
            sage: S.glue((0, 1), (0, 1))
            sage: S.replace_polygon(0, H.polygon(vertices=[(0, 0), (1/3, 0), (0, 1/3)]))
            sage: S.polygon(0)
            Hyperbolic polygon with vertices [I, (1/3, 0), (0, 1/3)]
            sage: S.opposite_edge(0, 1)
            (0, 1)

        TESTS::

            sage: from flatsurf import polygons
            sage: S.replace_polygon(0, polygons.square())
            Traceback (most recent call last):
            ...
            TypeError: polygon must be a hyperbolic polygon given by its vertices

        """
        self._check_polygon(polygon)

        return super().replace_polygon(label, polygon.change_ring(self.base_ring()))

    @staticmethod
    def _check_polygon(polygon):
        r"""
        Raise a ``TypeError`` if ``polygon`` cannot be a polygon of this
        surface.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane, polygons
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: MutableOrientedHyperbolicSurface._check_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            sage: MutableOrientedHyperbolicSurface._check_polygon(polygons.square())
            Traceback (most recent call last):
            ...
            TypeError: polygon must be a hyperbolic polygon given by its vertices

        """
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

        if not isinstance(polygon, HyperbolicPolygon):
            raise TypeError(
                "polygon must be a hyperbolic polygon given by its vertices"
            )

    def set_immutable(self):
        r"""
        Make this surface immutable, after checking that the edges glued to
        each other have the same length.

        Any mutation attempts from now on will be an error.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: H = HyperbolicPlane(QQ)
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S.add_polygon(H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
            0
            sage: S.glue((0, 0), (0, 1))
            sage: S.set_immutable()
            Traceback (most recent call last):
            ...
            ValueError: edge 0 of polygon 0 and edge 1 of polygon 0 are glued but do not have the same length
            sage: S.is_mutable()
            True

            sage: S.glue((0, 0), (0, 3))
            sage: S.glue((0, 1), (0, 2))
            sage: S.set_immutable()
            sage: S.is_mutable()
            False
            sage: S.category()
            Category of connected without boundary finite type oriented hyperbolic cone surfaces

        """
        if self._mutable:
            mismatch = self._mismatched_gluing()
            if mismatch is not None:
                (label, edge), (cross_label, cross_edge) = mismatch
                raise ValueError(
                    f"edge {edge} of polygon {label} and edge {cross_edge} of polygon {cross_label} are glued but do not have the same length"
                )

        super().set_immutable()
