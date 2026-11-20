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
      polygons, or a
      :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane`, the plane of
      the polygons, whose quadratic form may differ from `x^2 + y^2 - z^2`

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

        from flatsurf.geometry.hyperbolic import HyperbolicPlane

        plane = None
        if isinstance(base, HyperbolicPlane):
            plane = base
            base = plane.base_ring()

        if base not in Rings():
            raise TypeError("base ring must be a ring")

        from flatsurf.geometry.categories import HyperbolicConeSurfaces

        if category is None:
            category = HyperbolicConeSurfaces().Oriented()

        category &= HyperbolicConeSurfaces().Oriented()

        # A plane of the form x^2 + y^2 - z^2 is not stored, so that a surface
        # given by its plane equals the surface given by its base ring. The
        # plane is not built from a ring, which may have no hyperbolic plane,
        # such as CC.
        self._hyperbolic_plane = None if plane is None or plane._standard else plane

        super().__init__(base, category=category)

    def hyperbolic_plane(self):
        r"""
        Return the hyperbolic plane that the polygons of this surface live
        in.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: MutableOrientedHyperbolicSurface(QQ).hyperbolic_plane()
            Hyperbolic Plane over Rational Field

        A sphere glued from a triangle of the plane with the form `x^2 + 3y^2 -
        z^2`, whose Klein model is the disk `x^2 + 3y^2 < 1`; its edge 1 is
        folded at its midpoint, an orbifold point of order 2::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ, form=[1, 3, -1])
            sage: C = H.polygon(vertices=[(1/2, 0), (-1/4, 1/4), (-1/4, -1/4)]).cyclic_quotient(3)
            sage: C
            Hyperbolic polygon with vertices [(0, 0), (1/2, 0), (-1/4, 1/4)]
            sage: S = MutableOrientedHyperbolicSurface(H)
            sage: S.add_polygon(C)
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 1))
            sage: S.set_immutable()
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field with form x^2 + 3*y^2 - z^2
            sage: S.edge_transformation(0, 0)
            [-1/2 -3/2    0]
            [ 1/2 -1/2    0]
            [   0    0    1]
            sage: S.edge_transformation(0, 0).parent()
            Group of isometries of the hyperbolic plane over Rational Field with form x^2 + 3*y^2 - z^2
            sage: S(0, 0).angle()
            1/3
            sage: S(0, 1).angle(numerical=True)  # abs tol 1e-12
            0.147583617650433
            sage: S.holonomy_representation().invariant_trace_field()
            Rational Field
            sage: TestSuite(S).run()

        TESTS:

        The plane is kept by pickling and by a change of ring::

            sage: loads(dumps(S)).hyperbolic_plane() is H
            True
            sage: loads(dumps(S)) == S
            True
            sage: S.change_ring(AA).hyperbolic_plane()
            Hyperbolic Plane over Algebraic Real Field with form x^2 + 3*y^2 - z^2

        """
        if self._hyperbolic_plane is not None:
            return self._hyperbolic_plane

        from flatsurf.geometry.hyperbolic import HyperbolicPlane

        return HyperbolicPlane(self.base_ring())


class MutableOrientedHyperbolicSurface(
    OrientedHyperbolicSurface, MutablePolygonalSurface
):
    r"""
    A surface built from compact polygons of the hyperbolic plane glued by
    orientation preserving isometries.

    INPUT:

    - ``base`` -- a ring, the base ring of the hyperbolic plane of the
      polygons, or a
      :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane`, the plane of
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

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        Surfaces in hyperbolic planes of different quadratic forms differ,
        also when they have no polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S == MutableOrientedHyperbolicSurface(HyperbolicPlane(QQ))
            True
            sage: S == MutableOrientedHyperbolicSurface(HyperbolicPlane(QQ, form=[1, 2, -3]))
            False
            sage: S != MutableOrientedHyperbolicSurface(HyperbolicPlane(QQ, form=[1, 2, -3]))
            True

        """
        if not super().__eq__(other):
            return False

        return self._hyperbolic_plane == getattr(other, "_hyperbolic_plane", None)

    # Defining __eq__ drops the inherited __hash__, which is compatible with
    # this equality: it does not depend on the plane.
    __hash__ = MutablePolygonalSurface.__hash__

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
        # put the new polygon in place for the checks and restore the old one
        # whatever their outcome, so that any exception leaves the surface
        # unchanged.
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
        finally:
            self._polygons[label] = old

        super().replace_polygon(label, polygon)

    def set_immutable(self, minimal_field=False):
        r"""
        Make this surface immutable.

        INPUT:

        - ``minimal_field`` -- a boolean (default: ``False``); whether to
          first move the surface to the number field generated by the Klein
          coordinates of its vertices, by its side pairings and by the
          coefficients of the quadratic form of its hyperbolic plane, which
          needs the base ring to coerce into ``AA``. The base ring is kept if
          it is a number field, or ``QQ``, of the degree of that field.

        EXAMPLES:

        A surface built over ``AA`` is moved to its minimal field::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: K.<r> = NumberField(x^4 + 32*x^2 - 32, embedding=AA(12*AA(2).sqrt() - 16).sqrt())
            sage: s2 = r^2/12 + 4/3
            sage: cos = [1, s2/2, 0, -s2/2, -1, -s2/2, 0, s2/2]
            sage: sin = [0, s2/2, 1, s2/2, 0, -s2/2, -1, -s2/2]
            sage: P = HyperbolicPlane(K).polygon(vertices=[(r * c, r * s) for c, s in zip(cos, sin)])
            sage: S = MutableOrientedHyperbolicSurface(AA)
            sage: S.add_polygon(P)
            0
            sage: for i in range(4):
            ....:     S.glue((0, i), (0, i + 4))
            sage: S.set_immutable(minimal_field=True)
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 2 with a = 1.189207115002722?
            sage: S.polygon(0).base_ring() is S.base_ring()
            True
            sage: S.edge_transformation(0, 0).trace()
            8*a^2 + 11
            sage: TestSuite(S).run()

        A base ring that is a number field of the right degree is kept::

            sage: S = MutableOrientedHyperbolicSurface(K)
            sage: S.add_polygon(P)
            0
            sage: for i in range(4):
            ....:     S.glue((0, i), (0, i + 4))
            sage: S.set_immutable(minimal_field=True)
            sage: S.base_ring() is K
            True
            sage: TestSuite(S).run()

        Over ``AA``, a surface defined over ``QQ`` moves to ``QQ``::

            sage: H = HyperbolicPlane(QQ)
            sage: Q = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: S = MutableOrientedHyperbolicSurface(AA)
            sage: S.add_polygon(Q)
            0
            sage: S.glue((0, 1), (0, 2))
            sage: S.set_immutable(minimal_field=True)
            sage: S.base_ring()
            Rational Field
            sage: S.polygon(0)
            Hyperbolic polygon with vertices [I, (1/2, 0), (1/2, 1/2), (0, 1/2)]
            sage: TestSuite(S).run()

        An edge glued to itself needs a square root::

            sage: S = MutableOrientedHyperbolicSurface(AA)
            sage: S.add_polygon(Q)
            0
            sage: S.glue((0, 1), (0, 2))
            sage: S.glue((0, 0), (0, 0))
            sage: S.set_immutable(minimal_field=True)
            sage: S.base_ring()
            Number Field in a with defining polynomial y^2 - 3 with a = 1.732050807568878?
            sage: S.polygon(0).base_ring() is S.base_ring()
            True
            sage: TestSuite(S).run()

        The plane moves with the polygons; the octagon above, stretched into
        the disk `x^2 + 2y^2 < 3`::

            sage: H = HyperbolicPlane(AA, form=[1, 2, -3])
            sage: O = H.polygon(vertices=[(AA(3).sqrt() * AA(v.coordinates(model="klein")[0]), AA(3/2).sqrt() * AA(v.coordinates(model="klein")[1])) for v in P.vertices()])
            sage: S = MutableOrientedHyperbolicSurface(H)
            sage: S.add_polygon(O)
            0
            sage: for i in range(4):
            ....:     S.glue((0, i), (0, i + 4))
            sage: S.edge_transformation(0, 0).trace()
            22.31370849898476?
            sage: S.set_immutable(minimal_field=True)
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 18 with a = 2.059767143907118?
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Number Field in a with defining polynomial y^4 - 18 with a = 2.059767143907118? with form x^2 + 2*y^2 - 3*z^2
            sage: S.edge_transformation(0, 0).trace()
            8/3*a^2 + 11
            sage: TestSuite(S).run()

        TESTS::

            sage: S.set_immutable(minimal_field=True)
            Traceback (most recent call last):
            ...
            Exception: cannot modify immutable surface; create a mutable copy first

        A base ring that does not embed into ``AA`` is refused::

            sage: S = MutableOrientedHyperbolicSurface(CC)
            sage: S.set_immutable(minimal_field=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: the minimal field is only computed for base rings that embed into AA

        """
        if not minimal_field:
            return super().set_immutable()

        if not self._mutable:
            raise Exception(
                "cannot modify immutable surface; create a mutable copy first"
            )

        from sage.all import AA, QQ

        if not AA.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError(
                "the minimal field is only computed for base rings that embed into AA"
            )

        elements = []
        if self._hyperbolic_plane is not None:
            form = self._hyperbolic_plane.quadratic_form()
            elements.extend(AA(form[i, i]) for i in range(3))
        for label, polygon in self._polygons.items():
            for vertex in polygon.vertices():
                elements.extend(AA(c) for c in vertex.coordinates(model="klein"))

        # Each gluing is seen from both of its edges; one of them suffices.
        done = set()
        for label, polygon in self._polygons.items():
            for edge in range(len(polygon.vertices())):
                cross = self.opposite_edge(label, edge)
                if cross is None or (label, edge) in done:
                    continue
                done.add((label, edge))
                done.add(cross)
                matrix = self._side_pairing_matrix(label, edge, *cross)
                elements.extend(AA(c) for c in matrix.list())

        from flatsurf.geometry.subfield import number_field_elements_from_algebraics

        field = number_field_elements_from_algebraics(elements)[0]

        from sage.rings.number_field.number_field_base import NumberField

        base = self.base_ring()
        if not (
            (base is QQ or isinstance(base, NumberField))
            and base.degree() == field.degree()
        ):
            self._base = field
            if self._hyperbolic_plane is not None:
                self._hyperbolic_plane = self._hyperbolic_plane.change_ring(
                    AA
                ).change_ring(field)
            for label in list(self._polygons):
                self._polygons[label] = (
                    self._polygons[label].change_ring(AA).change_ring(field)
                )

        super().set_immutable()

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

        TESTS:

        A polygon over a smaller ring is added over the ring of the surface,
        and a Euclidean polygon is refused by :meth:`add_polygon`::

            sage: A = MutableOrientedHyperbolicSurface(AA)
            sage: A.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            0
            sage: A.polygon(0).base_ring()
            Algebraic Real Field
            sage: A.add_polygon(polygons.square())
            Traceback (most recent call last):
            ...
            TypeError: polygon must be a hyperbolic polygon given by its vertices

        The polygon must lie in a hyperbolic plane of the quadratic form of the
        plane of the surface::

            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: P = H3.polygon(vertices=[(0, 0), (1/2, -1/2), (3/4, 0), (1/2, 1/2)])
            sage: S.add_polygon(P)
            Traceback (most recent call last):
            ...
            ValueError: polygon lies in a hyperbolic plane with form x^2 + 2*y^2 - 3*z^2, the surface in a hyperbolic plane with form x^2 + y^2 - z^2
            sage: T = MutableOrientedHyperbolicSurface(H3)
            sage: T.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            Traceback (most recent call last):
            ...
            ValueError: polygon lies in a hyperbolic plane with form x^2 + y^2 - z^2, the surface in a hyperbolic plane with form x^2 + 2*y^2 - 3*z^2
            sage: T.add_polygon(P)
            0

        A polygon over a smaller ring with the same form is accepted::

            sage: U = MutableOrientedHyperbolicSurface(H3.change_ring(AA))
            sage: U.add_polygon(P)
            0
            sage: U.polygon(0).parent() is H3.change_ring(AA)
            True

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

        plane = polygon.parent()
        if self._hyperbolic_plane is None:
            same_form = plane._standard
        else:
            ring = self.base_ring()
            same_form = plane._diagonal_form_in(
                ring
            ) == self._hyperbolic_plane._diagonal_form_in(ring)
        if not same_form:
            from flatsurf.geometry.hyperbolic_isometry import _minkowski_form
            from flatsurf.geometry.hyperbolic_polygon import _form_repr

            form = _minkowski_form(self.base_ring())
            if self._hyperbolic_plane is not None:
                form = self._hyperbolic_plane.quadratic_form()
            raise ValueError(
                f"polygon lies in a hyperbolic plane with form {_form_repr(plane.quadratic_form())}, the surface in a hyperbolic plane with form {_form_repr(form)}"
            )
