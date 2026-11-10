r"""
Surfaces built by gluing compact polygons of the hyperbolic plane.

The polygons are
:class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon` objects, whose
vertices are finite points of the hyperbolic plane. Two edges are glued by the
orientation preserving isometry that sends one onto the other one with reverse
orientation, which requires the two edges to have the same length; this is
checked by ``glue()``, which also checks that the isometry exists over the base
ring of the surface.

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
    Traceback (most recent call last):
    ...
    ValueError: cannot glue edge 0 of polygon 0 to edge 1 of polygon 0: segments do not have the same length

An edge glued to itself is folded at its midpoint by a half turn, which is
not always defined over the base ring::

    sage: T.glue((0, 0), (0, 0))
    Traceback (most recent call last):
    ...
    ValueError: cannot glue edge 0 of polygon 0 to edge 0 of polygon 0: isometry is not defined over Rational Field; adjoin a square root of 4/3

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

    - ``base`` -- a ring, the base ring of the hyperbolic plane of the
      polygons

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

    def glue(self, x, y):
        r"""
        Glue the edge ``x`` to the edge ``y``, after checking that the
        isometry that sends one onto the other exists over the base ring.

        Gluings of ``x`` and ``y`` that already exist are removed, unless the
        new gluing is refused, in which case the surface is unchanged.

        INPUT:

        - ``x`` -- a pair (polygon label, edge index)

        - ``y`` -- a pair (polygon label, edge index)

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S.add_polygon(P)
            0
            sage: S.glue((0, 0), (0, 3))

        Glued edges must have the same length::

            sage: S.glue((0, 0), (0, 1))
            Traceback (most recent call last):
            ...
            ValueError: cannot glue edge 0 of polygon 0 to edge 1 of polygon 0: segments do not have the same length
            sage: S.opposite_edge(0, 0)
            (0, 3)

        An edge glued to itself is folded at its midpoint by a half turn,
        which is not always defined over the base ring::

            sage: S.glue((0, 0), (0, 0))
            Traceback (most recent call last):
            ...
            ValueError: cannot glue edge 0 of polygon 0 to edge 0 of polygon 0: isometry is not defined over Rational Field; adjoin a square root of 4/3

        TESTS::

            sage: S.glue((0, 0), (1, 0))
            Traceback (most recent call last):
            ...
            ValueError: there is no polygon 1 in this surface
            sage: S.glue((0, 7), (0, 0))
            Traceback (most recent call last):
            ...
            ValueError: polygon 0 has no edge 7

        """
        if not self._mutable:
            raise Exception(
                "cannot modify immutable surface; create a mutable copy first"
            )

        for label, edge in (x, y):
            if label not in self._polygons:
                raise ValueError(f"there is no polygon {label} in this surface")
            if not 0 <= edge < len(self._polygons[label].vertices()):
                raise ValueError(f"polygon {label} has no edge {edge}")

        self._check_side_pairing(x, y)

        super().glue(x, y)

    def _check_side_pairing(self, x, y):
        r"""
        Raise a ``ValueError`` if the edges ``x`` and ``y`` cannot be glued
        over the base ring.

        INPUT:

        - ``x``, ``y`` -- pairs (polygon label, edge index)

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]))
            0
            sage: S._check_side_pairing((0, 0), (0, 3))
            sage: S._check_side_pairing((0, 0), (0, 1))
            Traceback (most recent call last):
            ...
            ValueError: cannot glue edge 0 of polygon 0 to edge 1 of polygon 0: segments do not have the same length

        """
        try:
            self._side_pairing_matrix(*x, *y)
        except ValueError as error:
            raise ValueError(
                f"cannot glue edge {x[1]} of polygon {x[0]} to edge {y[1]} of polygon {y[0]}: {error}"
            ) from error

    def replace_polygon(self, label, polygon):
        r"""
        Replace the polygon ``label`` with ``polygon`` while keeping its
        gluings intact.

        The replacement is refused, and the surface left unchanged, if one of
        the gluings of the polygon is not defined over the base ring for
        ``polygon``.

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

        A replacement for which a gluing does not exist is refused, and the
        surface is left unchanged::

            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: T = MutableOrientedHyperbolicSurface(QQ)
            sage: T.add_polygon(P)
            0
            sage: T.glue((0, 0), (0, 3))
            sage: T.glue((0, 1), (0, 2))
            sage: T.replace_polygon(0, H.polygon(vertices=[(0, 0), (1/3, 0), (1/2, 1/2), (0, 1/2)]))
            Traceback (most recent call last):
            ...
            ValueError: cannot replace polygon 0: cannot glue edge 0 of polygon 0 to edge 3 of polygon 0: segments do not have the same length
            sage: T.polygon(0) == P
            True

        An immutable surface refuses the replacement::

            sage: S.set_immutable()
            sage: S.replace_polygon(0, H.polygon(vertices=[(0, 0), (1/3, 0), (0, 1/3)]))
            Traceback (most recent call last):
            ...
            Exception: cannot modify immutable surface; create a mutable copy first

        """
        if not self._mutable:
            raise Exception(
                "cannot modify immutable surface; create a mutable copy first"
            )

        self._check_polygon(polygon)

        polygon = polygon.change_ring(self.base_ring())

        if len(polygon.vertices()) != len(self.polygon(label).vertices()):
            # The shared method refuses the replacement with its message.
            return super().replace_polygon(label, polygon)

        # The side pairings are computed from the polygons of the surface:
        # put the new polygon in place and restore the old one if one of its
        # gluings does not exist.
        old = self._polygons[label]
        self._polygons[label] = polygon
        try:
            for edge in range(len(polygon.vertices())):
                cross = self.opposite_edge(label, edge)
                if cross is None or (cross[0] == label and cross[1] < edge):
                    continue
                try:
                    self._check_side_pairing((label, edge), cross)
                except ValueError as error:
                    raise ValueError(
                        f"cannot replace polygon {label}: {error}"
                    ) from error
        except ValueError:
            self._polygons[label] = old
            raise

        self._polygons[label] = old
        super().replace_polygon(label, polygon)

    def _check_polygon(self, polygon):
        r"""
        Raise an error if ``polygon`` cannot be a polygon of this surface.

        A ``TypeError`` is raised for a polygon that is not hyperbolic, and a
        ``ValueError`` for a polygon whose base ring does not coerce into the
        base ring of this surface.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane, polygons
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S._check_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            sage: S._check_polygon(polygons.square())
            Traceback (most recent call last):
            ...
            TypeError: polygon must be a hyperbolic polygon given by its vertices
            sage: S._check_polygon(HyperbolicPlane(AA).polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            Traceback (most recent call last):
            ...
            ValueError: polygon is defined over Algebraic Real Field, which does not coerce into Rational Field

        """
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

        if not isinstance(polygon, HyperbolicPolygon):
            raise TypeError(
                "polygon must be a hyperbolic polygon given by its vertices"
            )

        if not self.base_ring().has_coerce_map_from(polygon.base_ring()):
            raise ValueError(
                f"polygon is defined over {polygon.base_ring()}, which does not coerce into {self.base_ring()}"
            )
