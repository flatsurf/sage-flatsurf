r"""
Generators of hyperbolic polygons and surfaces.

The polygons are built exactly from their interior angles, given as
multiples of `2\pi`, and from the hyperbolic cosine of some of their sides.
Each polygon has its vertex 0 at the origin of the Klein model and its edge
0 on the positive `x`-axis, with its vertices in counterclockwise order.

EXAMPLES::

    sage: from flatsurf import hyperbolic_polygons
    sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 145/24)
    sage: P.base_ring()
    Rational Field
    sage: P.angles()
    (1/4, 1/4, 1/4, 1/4, 1/4, 1/4)

"""

# ****************************************************************************
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
# ****************************************************************************

from sage.rings.all import AA, QQ, QQbar, ZZ
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix, identity_matrix
from sage.modules.free_module_element import vector


def _angle(angle, name="angle"):
    r"""
    Return ``angle`` as a rational number, the multiple of `2\pi` of an
    angle.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _angle
        sage: _angle(1/4)
        1/4
        sage: _angle(sqrt(2))
        Traceback (most recent call last):
        ...
        ValueError: angle must be a rational multiple of 2π

    """
    try:
        return QQ(angle)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a rational multiple of 2π")


def _cos_sin(angle):
    r"""
    Return the cosine and the sine of `2\pi` ``angle`` in ``AA``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _cos_sin
        sage: _cos_sin(1/8)
        (0.7071067811865475?, 0.7071067811865475?)

    """
    angle = QQ(angle)
    z = QQbar.zeta(angle.denominator()) ** angle.numerator()
    return AA(z.real()), AA(z.imag())


def _from_cosh(cosh, name="cosh"):
    r"""
    Return the hyperbolic cosine and sine in ``AA`` of the side of
    hyperbolic cosine ``cosh``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _from_cosh
        sage: _from_cosh(5/4)
        (5/4, 3/4)
        sage: _from_cosh(1)
        Traceback (most recent call last):
        ...
        ValueError: cosh must be larger than 1

    """
    try:
        cosh = AA(cosh)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a real algebraic number")
    if not cosh > 1:
        raise ValueError(f"{name} must be larger than 1")
    return cosh, (cosh * cosh - 1).sqrt()


def _from_sinh(sinh):
    r"""
    Return the hyperbolic cosine and sine in ``AA`` of the side of positive
    hyperbolic sine ``sinh``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _from_sinh
        sage: _from_sinh(AA(3/4))
        (5/4, 3/4)

    """
    return (sinh * sinh + 1).sqrt(), sinh


def _walk(sides, angles):
    r"""
    Return the Klein coordinates in ``AA`` of the vertices of the polygon
    reached by a walk from the origin along the positive `x`-axis.

    INPUT:

    - ``sides`` -- a list of pairs ``(cosh, sinh)`` of the sides 0, 1, ...
      of the walk

    - ``angles`` -- a list of pairs ``(cos, sin)`` of the interior angles at
      the vertices 1, 2, ..., one fewer than ``sides``

    ALGORITHM:

    A frame of the hyperboloid model starts at the origin, pointing along
    the positive `x`-axis. A side of hyperbolic cosine `c` and sine `s`
    multiplies the frame on the right by the translation `((c, 0, s), (0, 1,
    0), (s, 0, c))` along its first axis, and a vertex of interior angle
    `\alpha` by the rotation by `\pi - \alpha`, a left turn. Each vertex is
    the last column of its frame, divided by its last entry.

    EXAMPLES:

    A walk along two sides of `\cosh = 5/4` with a right angle between
    them::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _walk
        sage: _walk([(AA(5/4), AA(3/4))] * 2, [(AA(0), AA(1))])
        [(0, 0), (3/5, 0), (3/5, 12/25)]

    """
    frame = identity_matrix(AA, 3)
    vertices = [(AA.zero(), AA.zero())]
    for i, (cosh, sinh) in enumerate(sides):
        frame = frame * matrix(AA, [[cosh, 0, sinh], [0, 1, 0], [sinh, 0, cosh]])
        vertex = frame.column(2)
        vertices.append((vertex[0] / vertex[2], vertex[1] / vertex[2]))
        if i < len(angles):
            cos, sin = angles[i]
            frame = frame * matrix(AA, [[-cos, -sin, 0], [sin, -cos, 0], [0, 0, 1]])
    return vertices


def _polygon(vertices, ring):
    r"""
    Return the hyperbolic polygon with the Klein coordinates ``vertices`` in
    ``AA``, over ``ring``.

    INPUT:

    - ``vertices`` -- a list of pairs of elements of ``AA``

    - ``ring`` -- a ring or ``None``; when ``None``, the field of
      :func:`~flatsurf.geometry.subfield.number_field_elements_from_algebraics`
      of the coordinates, ``QQ`` when they are rational; otherwise the
      coordinates are converted into ``ring``

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _polygon
        sage: vertices = [(AA(0), AA(0)), (AA(1/2), AA(0)), (AA(0), AA(2).sqrt() / 2)]
        sage: _polygon(vertices, None).base_ring()
        Number Field in a with defining polynomial y^2 - 2 with a = 1.414213562373095?
        sage: _polygon(vertices, AA).base_ring()
        Algebraic Real Field
        sage: _polygon(vertices, QQ)
        Traceback (most recent call last):
        ...
        ValueError: the coordinates of the polygon do not lie in Rational Field

    """
    from flatsurf.geometry.hyperbolic import HyperbolicPlane

    coordinates = [c for vertex in vertices for c in vertex]

    if ring is None:
        from flatsurf.geometry.subfield import number_field_elements_from_algebraics

        ring, coordinates = number_field_elements_from_algebraics(coordinates)
    elif ring is not AA:
        try:
            coordinates = [ring(c) for c in coordinates]
        except (TypeError, ValueError):
            raise ValueError(f"the coordinates of the polygon do not lie in {ring}")

    return HyperbolicPlane(ring).polygon(
        vertices=list(zip(coordinates[::2], coordinates[1::2]))
    )


_RIGHT_ANGLE = (AA.zero(), AA.one())


class HyperbolicPolygonsConstructor:
    r"""
    Constructors of compact hyperbolic polygons from their angles and the
    hyperbolic cosine of some of their sides.

    Angles are multiples of `2\pi`, so that `1/4` is a right angle. Each
    polygon has its vertex 0 at the origin of the Klein model and its edge 0
    on the positive `x`-axis. With ``ring=None``, the polygon is defined
    over the number field generated by the Klein coordinates of its
    vertices, ``QQ`` when they are rational; ``ring=AA`` keeps the
    coordinates in ``AA``; any other ring receives the coordinates by
    conversion.

    EXAMPLES::

        sage: from flatsurf import hyperbolic_polygons
        sage: hyperbolic_polygons
        hyperbolic polygons

    """

    def __repr__(self):
        r"""
        Return a printable representation of this constructor.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_polygons
            sage: hyperbolic_polygons
            hyperbolic polygons

        """
        return "hyperbolic polygons"

    @staticmethod
    def triangle(a, b, c, ring=None):
        r"""
        Return the hyperbolic triangle with interior angles ``a``, ``b`` and
        ``c`` at its vertices 0, 1 and 2.

        INPUT:

        - ``a``, ``b``, ``c`` -- positive rational numbers, the angles as
          multiples of `2\pi`, of sum less than `1/2`

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the triangle, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin and edge 0 on the positive `x`-axis.

        ALGORITHM:

        The side `c` opposite the angle `\gamma` has `\cosh c = (\cos\alpha
        \cos\beta + \cos\gamma) / (\sin\alpha \sin\beta)`, the second
        cosine rule of Katok, Fuchsian groups (1992), Theorem 1.5.2, p. 16.

        The field of the coordinates depends on the vertex at the origin,
        chosen by the order of the angles.

        EXAMPLES:

        A triangle with a right angle at the origin::

            sage: from flatsurf import hyperbolic_polygons
            sage: T = hyperbolic_polygons.triangle(1/4, 1/16, 1/16)
            sage: T.base_ring()
            Number Field in a with defining polynomial y^4 - 2*y^2 - 1 with a = ...
            sage: T.angles()
            (1/4, 1/16, 1/16)
            sage: T.area()
            1/8

        The triangle of the triangle group `(2, 3, 7)`, of degree 6::

            sage: T = hyperbolic_polygons.triangle(1/4, 1/6, 1/14)  # long time
            sage: T.base_ring().polynomial()  # long time
            y^6 + y^4 - 2*y^2 - 1
            sage: T.area()  # long time
            1/84

        The coordinates can be kept in ``AA``::

            sage: T = hyperbolic_polygons.triangle(QQ((1, 4)), QQ((1, 16)), QQ((1, 16)), ring=AA)
            sage: T.base_ring()
            Algebraic Real Field
            sage: T.angles()
            (1/4, 1/16, 1/16)

        A triangle whose angles sum to `1/2` is Euclidean::

            sage: hyperbolic_polygons.triangle(1/4, 1/8, 1/8)
            Traceback (most recent call last):
            ...
            ValueError: the angles of a hyperbolic triangle must sum to less than 1/2

        TESTS::

            sage: hyperbolic_polygons.triangle(0, 1/8, 1/8)
            Traceback (most recent call last):
            ...
            ValueError: the angles of a hyperbolic triangle must be positive
            sage: hyperbolic_polygons.triangle(1/4, 1/16, 1/16, ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field

        """
        a, b, c = _angle(a), _angle(b), _angle(c)
        if not (a > 0 and b > 0 and c > 0):
            raise ValueError("the angles of a hyperbolic triangle must be positive")
        if not a + b + c < QQ((1, 2)):
            raise ValueError(
                "the angles of a hyperbolic triangle must sum to less than 1/2"
            )

        cos_a, sin_a = _cos_sin(a)
        cos_b, sin_b = _cos_sin(b)
        cos_c, sin_c = _cos_sin(c)

        side0 = _from_cosh((cos_a * cos_b + cos_c) / (sin_a * sin_b))
        side1 = _from_cosh((cos_b * cos_c + cos_a) / (sin_b * sin_c))

        return _polygon(_walk([side0, side1], [(cos_b, sin_b)]), ring)

    @staticmethod
    def regular_ngon(n, angle, ring=None):
        r"""
        Return the regular hyperbolic polygon with ``n`` vertices and
        interior angle ``angle``.

        INPUT:

        - ``n`` -- an integer, at least 3

        - ``angle`` -- a rational number, the interior angle as a multiple of
          `2\pi`, in `(0, 1/2 - 1/n)`

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the polygon, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin and edge 0 on the positive `x`-axis.

        ALGORITHM:

        The triangles from the centre to the sides have angles `2\pi/n` at
        the centre and `\theta/2` at the vertices, with `\theta = 2\pi`
        ``angle``. The second cosine rule of Katok, Fuchsian groups (1992),
        Theorem 1.5.2, p. 16, gives the side `\cosh s = (\cos^2(\theta/2) +
        \cos(2\pi/n)) / \sin^2(\theta/2)`; the octagon of angle `1/8` is
        the one of Katok, Section 4.3, Example C, p. 95.

        EXAMPLES:

        The octagon of the genus two surface of Katok::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.regular_ngon(8, 1/8)
            sage: P.base_ring()
            Number Field in a with defining polynomial y^4 - 2*y^2 - 1 with a = ...
            sage: P.angles()
            (1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8)
            sage: P.area()
            2

        Python integers are accepted::

            sage: P = hyperbolic_polygons.regular_ngon(int(5), QQ((1, 5)), ring=AA)
            sage: P.angles()
            (1/5, 1/5, 1/5, 1/5, 1/5)

        A square with right angles is Euclidean::

            sage: hyperbolic_polygons.regular_ngon(4, 1/4)
            Traceback (most recent call last):
            ...
            ValueError: the angle of a regular hyperbolic polygon with 4 vertices must be in (0, 1/4)

        TESTS::

            sage: hyperbolic_polygons.regular_ngon(2, 1/8)
            Traceback (most recent call last):
            ...
            ValueError: a polygon has at least 3 vertices

        """
        n = ZZ(n)
        if n < 3:
            raise ValueError("a polygon has at least 3 vertices")
        angle = _angle(angle)
        bound = QQ((1, 2)) - QQ((1, n))
        if not 0 < angle < bound:
            raise ValueError(
                f"the angle of a regular hyperbolic polygon with {n} vertices must be in (0, {bound})"
            )

        cos_centre, _ = _cos_sin(QQ((1, n)))
        cos_half, sin_half = _cos_sin(angle / 2)
        side = _from_cosh((cos_half**2 + cos_centre) / sin_half**2)

        turn = _cos_sin(angle)
        return _polygon(_walk([side] * (n - 1), [turn] * (n - 2)), ring)

    @staticmethod
    def right_angled_hexagon(c0, c1, c2, ring=None):
        r"""
        Return the right-angled hyperbolic hexagon whose edges 1, 3 and 5
        have hyperbolic cosine ``c0``, ``c1`` and ``c2``.

        The edges 0, 2 and 4 follow: edge 0 is opposite edge 3, edge 2 is
        opposite edge 5, and edge 4 is opposite edge 1.

        INPUT:

        - ``c0``, ``c1``, ``c2`` -- real algebraic numbers larger than 1

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the hexagon, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin, edge 0 on the positive `x`-axis, and
        edge 5 ending at vertex 0.

        ALGORITHM:

        Such a hexagon exists for any three lengths of its alternate sides,
        Buser, Geometry and spectra of compact Riemann surfaces (1992),
        Lemma 1.7.1. The side `b_k` opposite the side `a_k` has `\cosh b_k =
        (\cosh a_k + \cosh a_i \cosh a_j) / (\sinh a_i \sinh a_j)`, Buser,
        Chapter 2.

        EXAMPLES:

        A hexagon with rational vertices::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 145/24)
            sage: P.base_ring()
            Rational Field
            sage: [v.coordinates(model="klein") for v in P.vertices()]
            [(0, 0), (56/65, 0), (56/65, 99/325), (840/4121, 4015/4121), (1344/9425, 143/145), (0, 143/145)]
            sage: P.angles()
            (1/4, 1/4, 1/4, 1/4, 1/4, 1/4)
            sage: P.area()
            1/2

        The three alternate sides equal::

            sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 5/4)
            sage: P.base_ring()
            Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?

        Python integers are accepted::

            sage: P = hyperbolic_polygons.right_angled_hexagon(int(2), int(2), int(2))
            sage: P.angles()
            (1/4, 1/4, 1/4, 1/4, 1/4, 1/4)

        Sides must have a hyperbolic cosine larger than 1::

            sage: hyperbolic_polygons.right_angled_hexagon(1, 5/4, 5/4)
            Traceback (most recent call last):
            ...
            ValueError: c0 must be larger than 1

        """
        h0 = _from_cosh(c0, "c0")
        h1 = _from_cosh(c1, "c1")
        h2 = _from_cosh(c2, "c2")

        def opposite(k, i, j):
            return _from_cosh((k[0] + i[0] * j[0]) / (i[1] * j[1]))

        sides = [
            opposite(h1, h0, h2),
            h0,
            opposite(h2, h0, h1),
            h1,
            opposite(h0, h1, h2),
        ]
        return _polygon(_walk(sides, [_RIGHT_ANGLE] * 4), ring)

    @staticmethod
    def right_angled_pentagon(angle, c0, c1, ring=None):
        r"""
        Return the hyperbolic pentagon with interior angle ``angle`` at
        vertex 1, right angles at the other vertices, and edges 2 and 4 of
        hyperbolic cosine ``c0`` and ``c1``.

        The edges 0 and 1 meet at vertex 1; edge 3 is opposite it.

        INPUT:

        - ``angle`` -- a rational number in `(0, 1/2)`, the angle as a
          multiple of `2\pi`

        - ``c0``, ``c1`` -- real algebraic numbers larger than 1

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the pentagon, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin, edge 0 on the positive `x`-axis, and
        edge 4 ending at vertex 0.

        ALGORITHM:

        The formula for right-angled hexagons of Buser, Geometry and spectra
        of compact Riemann surfaces (1992), Chapter 2, with a side replaced
        by the angle `\varphi = 2\pi` ``angle``: the edge 3 has `\cosh s_3 =
        (\cos\varphi + \cosh s_2 \cosh s_4) / (\sinh s_2 \sinh s_4)`, and the
        edges 1 and 0 have `\sinh s_1 = (\cosh s_4 + \cosh s_2 \cos\varphi)
        / (\sinh s_2 \sin\varphi)` and symmetrically. The pentagon exists if
        and only if these two hyperbolic sines are positive; for `\varphi =
        \pi/2` this is `\cosh s_4 = \sinh s_1 \sinh s_2`.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_pentagon(3/8, 5/4, 5/4)
            sage: P.base_ring().polynomial()
            y^4 - 56*y^2 - 72*y + 414
            sage: P.angles()
            (1/4, 3/8, 1/4, 1/4, 1/4)
            sage: P.area()
            1/8

        A pentagon of degree 6::

            sage: P = hyperbolic_polygons.right_angled_pentagon(1/14, 5/4, 5/4)
            sage: P.base_ring().degree()
            6
            sage: P.area()
            3/7

        Python integers are accepted::

            sage: P = hyperbolic_polygons.right_angled_pentagon(QQ((1, 6)), int(2), int(2), ring=AA)
            sage: P.angles()
            (1/4, 1/6, 1/4, 1/4, 1/4)

        An obtuse angle requires sides of comparable lengths::

            sage: hyperbolic_polygons.right_angled_pentagon(3/8, 145/24, 5/4)
            Traceback (most recent call last):
            ...
            ValueError: no right-angled pentagon with these parameters: c1 + c0 cos(2π angle) must be positive

        TESTS::

            sage: hyperbolic_polygons.right_angled_pentagon(1/2, 2, 2)
            Traceback (most recent call last):
            ...
            ValueError: the angle must be in (0, 1/2)
            sage: hyperbolic_polygons.right_angled_pentagon(1/3, 2, 5/4).angles()
            (1/4, 1/3, 1/4, 1/4, 1/4)
            sage: hyperbolic_polygons.right_angled_pentagon(1/3, 3, 5/4)
            Traceback (most recent call last):
            ...
            ValueError: no right-angled pentagon with these parameters: c1 + c0 cos(2π angle) must be positive
            sage: hyperbolic_polygons.right_angled_pentagon(1/3, 5/4, 3)
            Traceback (most recent call last):
            ...
            ValueError: no right-angled pentagon with these parameters: c0 + c1 cos(2π angle) must be positive

        """
        angle = _angle(angle)
        if not 0 < angle < QQ((1, 2)):
            raise ValueError("the angle must be in (0, 1/2)")
        h0 = _from_cosh(c0, "c0")
        h1 = _from_cosh(c1, "c1")
        cos, sin = _cos_sin(angle)

        if not h1[0] + h0[0] * cos > 0:
            raise ValueError(
                "no right-angled pentagon with these parameters: c1 + c0 cos(2π angle) must be positive"
            )
        if not h0[0] + h1[0] * cos > 0:
            raise ValueError(
                "no right-angled pentagon with these parameters: c0 + c1 cos(2π angle) must be positive"
            )

        side0 = _from_sinh((h0[0] + h1[0] * cos) / (h1[1] * sin))
        side1 = _from_sinh((h1[0] + h0[0] * cos) / (h0[1] * sin))
        side3 = _from_cosh((cos + h0[0] * h1[0]) / (h0[1] * h1[1]))

        return _polygon(
            _walk([side0, side1, h0, side3], [(cos, sin), _RIGHT_ANGLE, _RIGHT_ANGLE]),
            ring,
        )

    @staticmethod
    def right_angled_quadrilateral(angle0, angle1, c, ring=None):
        r"""
        Return the hyperbolic quadrilateral with interior angles ``angle0``
        and ``angle1`` at vertices 1 and 2, right angles at vertices 3 and
        0, and edge 3 of hyperbolic cosine ``c``.

        Edge 1, between the vertices 1 and 2, is opposite edge 3.

        INPUT:

        - ``angle0``, ``angle1`` -- positive rational numbers, the angles as
          multiples of `2\pi`, of sum less than `1/2`

        - ``c`` -- a real algebraic number larger than 1

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the quadrilateral, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin, edge 0 on the positive `x`-axis, and
        edge 3 ending at vertex 0.

        ALGORITHM:

        The formula for right-angled hexagons of Buser, Geometry and spectra
        of compact Riemann surfaces (1992), Chapter 2, with two sides
        replaced by the angles `\varphi_0, \varphi_1`: the edge 0, at the
        vertex of `\varphi_0`, has `\sinh h_0 = (\cos\varphi_1 + \cosh a
        \cos\varphi_0) / (\sinh a \sin\varphi_0)`, and symmetrically for
        the edge 2. The quadrilateral exists if and only if these two
        hyperbolic sines are positive; for `\varphi_1 = \pi/2` this is the
        relation `\cos\varphi_0 = \sinh a \sinh h_2` of a Lambert
        quadrilateral.

        Vertex 1 is reached along edge 0, and the vertices 3 and 2 along
        the edges 3 and 2 walked from the origin in the mirror image under
        `(x, y) \mapsto (y, x)`, so that the edge 1 is never needed.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_quadrilateral(1/4, 1/8, 5/4)
            sage: P.base_ring()
            Number Field in a with defining polynomial y^2 - 34 with a = 5.830951894845300?
            sage: P.angles()
            (1/4, 1/4, 1/8, 1/4)
            sage: P.area()
            1/8

        Python integers are accepted::

            sage: P = hyperbolic_polygons.right_angled_quadrilateral(QQ((1, 6)), QQ((1, 6)), int(2), ring=AA)
            sage: P.angles()
            (1/4, 1/6, 1/6, 1/4)

        An obtuse angle requires a short edge 3::

            sage: hyperbolic_polygons.right_angled_quadrilateral(3/8, 1/16, 145/24)
            Traceback (most recent call last):
            ...
            ValueError: no right-angled quadrilateral with these parameters: cos(2π angle1) + c cos(2π angle0) must be positive

        TESTS::

            sage: hyperbolic_polygons.right_angled_quadrilateral(3/8, 1/16, 5/4).angles()
            (1/4, 3/8, 1/16, 1/4)
            sage: hyperbolic_polygons.right_angled_quadrilateral(1/16, 3/8, 145/24)
            Traceback (most recent call last):
            ...
            ValueError: no right-angled quadrilateral with these parameters: cos(2π angle0) + c cos(2π angle1) must be positive
            sage: hyperbolic_polygons.right_angled_quadrilateral(1/4, 1/4, 2)
            Traceback (most recent call last):
            ...
            ValueError: the two angles must be positive and sum to less than 1/2
            sage: hyperbolic_polygons.right_angled_quadrilateral(1/4, 1/8, 1)
            Traceback (most recent call last):
            ...
            ValueError: c must be larger than 1

        """
        angle0 = _angle(angle0, "angle0")
        angle1 = _angle(angle1, "angle1")
        if not (angle0 > 0 and angle1 > 0 and angle0 + angle1 < QQ((1, 2))):
            raise ValueError("the two angles must be positive and sum to less than 1/2")
        h = _from_cosh(c, "c")
        cos0, sin0 = _cos_sin(angle0)
        cos1, sin1 = _cos_sin(angle1)

        if not cos1 + h[0] * cos0 > 0:
            raise ValueError(
                "no right-angled quadrilateral with these parameters: cos(2π angle1) + c cos(2π angle0) must be positive"
            )
        if not cos0 + h[0] * cos1 > 0:
            raise ValueError(
                "no right-angled quadrilateral with these parameters: cos(2π angle0) + c cos(2π angle1) must be positive"
            )

        side0 = _from_sinh((cos1 + h[0] * cos0) / (h[1] * sin0))
        side2 = _from_sinh((cos0 + h[0] * cos1) / (h[1] * sin1))

        forward = _walk([side0], [])
        mirror = _walk([h, side2], [_RIGHT_ANGLE])
        vertices = forward + [mirror[2][::-1], mirror[1][::-1]]
        return _polygon(vertices, ring)


hyperbolic_polygons = HyperbolicPolygonsConstructor()


def _reflection(p, q):
    r"""
    Return the matrix of the reflection in the hyperbolic line through the
    points of Klein coordinates ``p`` and ``q``, acting on the Klein lifts
    `(x, y, 1)`.

    ALGORITHM:

    With `J = \mathrm{diag}(1, 1, -1)` and `n = J (\hat p \times \hat q)`
    the Minkowski normal of the plane of the lifts `\hat p, \hat q`, the
    reflection is `R = I - 2 n (J n)^T / \langle n, n \rangle`, where
    `\langle n, n \rangle = n^T J n` is positive. Its coefficients are
    rational in the coordinates of ``p`` and ``q``.

    EXAMPLES:

    The reflection in the edge 1 of the triangle `(0, 0), (1/2, 0), (1/3,
    1/2)`, which sends the origin to `(36/49, 12/49)`::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _reflection
        sage: R = _reflection((1/2, 0), (1/3, 1/2))
        sage: R
        [-41/31 -24/31  36/31]
        [-24/31  23/31  12/31]
        [-36/31 -12/31  49/31]
        sage: R.determinant(), (R^2).is_one()
        (-1, True)
        sage: R * vector((1/2, 0, 1))
        (1/2, 0, 1)

    """
    J = diagonal_matrix([1, 1, -1])
    n = J * vector((p[0], p[1], 1)).cross_product(vector((q[0], q[1], 1)))
    return identity_matrix(3) - 2 * matrix(n).transpose() * matrix(J * n) / (n * J * n)


def _optimize_placement(vertices):
    r"""
    Return the image of the polygon with the Klein coordinates ``vertices``
    by an isometry that makes its coordinates rational, or ``vertices`` when
    no such isometry is found.

    INPUT:

    - ``vertices`` -- a list of pairs of elements of ``AA``, the Klein
      coordinates of the vertices of a polygon

    OUTPUT:

    A list of pairs of elements of ``AA`` of the same length, the image of
    the vertices, in the same order, by an orientation preserving isometry;
    ``vertices`` itself when their coordinates are rational, generate a
    field of degree other than 2, or when no image is rational.

    ALGORITHM:

    When the coordinates generate `\QQ(\sqrt{d})`, with `d` the squarefree
    part of the discriminant, the point `(x_0, 0)` with `x_0 = (q^2 -
    dp^2)/(q^2 + dp^2)` has `1 - x_0^2 \in d \QQ^2`. For each vertex `i`
    and its neighbours `j = i + 1, i - 1`, in this order, the isometry that
    sends the vertex `i` to `(x_0, 0)` and the direction of the vertex `j`
    to the positive `x`-axis is applied; the first rational image is
    returned.

    The start point `x_0` is the one of smallest height, then the smallest,
    among `1 \le p, q \le 40`. Any two start points with `1 - x_0^2 \in d
    \QQ^2` differ by a boost along the `x`-axis with rational `\cosh` and
    `\sinh`, so the choice of `x_0` only changes the size of the rational
    coordinates, not whether there are some; `p = q = 1` is always a
    candidate.

    EXAMPLES:

    The right-angled hexagon whose alternate sides have `\cosh = 5/4`, the
    half of a pair of pants with three cuffs of trace `21/4`, has its
    coordinates in `\QQ(\sqrt{6})` when its vertex 0 is at the origin, and
    rational ones once moved::

        sage: from flatsurf import hyperbolic_polygons
        sage: from flatsurf.geometry.hyperbolic_surface_generators import _optimize_placement
        sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 5/4, ring=AA)
        sage: vertices = [v.coordinates(model="klein") for v in P.vertices()]
        sage: _optimize_placement(vertices)
        [(-1/5, 0), (-1/5, -24/25), (-5/61, -60/61), (23/65, -12/13), (5/11, -48/55), (5/11, 0)]

    Rational coordinates are kept::

        sage: vertices = [(AA(0), AA(0)), (AA(1/2), AA(0)), (AA(0), AA(1/2))]
        sage: _optimize_placement(vertices) is vertices
        True

    """
    from flatsurf.geometry.hyperbolic_isometry import tangent_vector_isometry
    from flatsurf.geometry.subfield import number_field_elements_from_algebraics

    field, _ = number_field_elements_from_algebraics(
        [c for vertex in vertices for c in vertex]
    )
    if field is QQ or field.degree() != 2:
        return vertices

    d = field.discriminant().squarefree_part()

    def height_and_value(x):
        return max(abs(x.numerator()), x.denominator()), x

    x0 = min(
        (
            QQ((q * q - d * p * p, q * q + d * p * p))
            for p in range(1, 41)
            for q in range(1, 41)
        ),
        key=height_and_value,
    )

    n = len(vertices)
    lifts = [vector(AA, (x, y, 1)) for x, y in vertices]
    for i in range(n):
        for j in ((i + 1) % n, (i - 1) % n):
            direction = (
                vertices[j][0] - vertices[i][0],
                vertices[j][1] - vertices[i][1],
            )
            M = tangent_vector_isometry(
                (vertices[i], direction), ((x0, 0), (1, 0)), ring=AA
            )
            images = []
            for lift in lifts:
                x, y, z = M * lift
                images.append((x / z, y / z))
            if all(c in QQ for image in images for c in image):
                return [(AA(QQ(x)), AA(QQ(y))) for x, y in images]

    return vertices


def _pants_polygon(boundary_traces, cone_orders, optimize_placement):
    r"""
    Return the Klein coordinates in ``AA`` of the front polygon of the pair
    of pants of :meth:`HyperbolicSurfaceGenerators.pants`, and its seams.

    INPUT:

    - ``boundary_traces``, ``cone_orders``, ``optimize_placement`` -- as in
      :meth:`HyperbolicSurfaceGenerators.pants`

    OUTPUT:

    A pair ``(vertices, seams)``: the list of the Klein coordinates of the
    vertices of the front polygon, and the list of its edges that are
    seams. The other edges are the half-cuffs, in the order of
    ``boundary_traces``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _pants_polygon
        sage: _pants_polygon((21/4, 21/4), (-1/2,), True)
        ([(-15/19, 0),
          (-15/19, -10/19),
          (-225/421, -340/421),
          (-9/25, -544/625),
          (-9/25, 0)],
         [0, 1, 3])

    """
    boundary_traces = tuple(boundary_traces)
    cone_orders = tuple(cone_orders)
    if len(boundary_traces) + len(cone_orders) != 3:
        raise ValueError(
            f"a pair of pants has three cuffs and cone points in total, not {len(boundary_traces) + len(cone_orders)}"
        )

    orders = []
    for k in cone_orders:
        try:
            order = QQ(k)
        except (TypeError, ValueError):
            raise ValueError(f"the cone order {k} must be a rational number")
        if order <= -1:
            raise ValueError(f"the cone order {order} is not larger than -1")
        if order == 0:
            raise ValueError("a cone order must be nonzero")
        if order > 0:
            raise NotImplementedError(
                f"cone points of angle larger than 2π are not implemented, the cone order {order} is positive"
            )
        orders.append(order)

    cosh = []
    for t in boundary_traces:
        try:
            trace = AA(t)
        except (TypeError, ValueError):
            raise ValueError(f"the boundary trace {t} must be a real algebraic number")
        if not trace > 3:
            raise ValueError(f"the boundary trace {t} must be larger than 3")
        cosh.append((trace + 1).sqrt() / 2)

    divisor_euler_characteristic = 2 + sum(orders) - len(cosh)
    if not divisor_euler_characteristic < 0:
        raise ValueError(
            f"the divisor Euler characteristic 2 + sum(cone_orders) - len(boundary_traces) must be negative, not {divisor_euler_characteristic}"
        )

    angles = [(1 + k) / 2 for k in orders]

    # The primitives test the inequalities below; the ones not listed
    # cannot fail once the cuffs, the cone orders and the Euler
    # characteristic are checked.
    if len(orders) == 0:
        constructor, arguments, seams = (
            hyperbolic_polygons.right_angled_hexagon,
            cosh,
            [0, 2, 4],
        )
        inequalities = {}
    elif len(orders) == 1:
        constructor, arguments, seams = (
            hyperbolic_polygons.right_angled_pentagon,
            [angles[0]] + cosh,
            [0, 1, 3],
        )
        slots = "with k0 the cone order of slot 0 and cs = √(ts + 1)/2 for the boundary trace ts of slot s"
        inequalities = {
            "c1 + c0 cos(2π angle)": f"c2 + c1 cos(π(1 + k0)) must be positive, {slots}",
            "c0 + c1 cos(2π angle)": f"c1 + c2 cos(π(1 + k0)) must be positive, {slots}",
        }
    elif len(orders) == 2:
        constructor, arguments, seams = (
            hyperbolic_polygons.right_angled_quadrilateral,
            angles + cosh,
            [0, 1, 2],
        )
        slots = "with k0, k1 the cone orders of slots 0, 1 and c2 = √(t2 + 1)/2 for the boundary trace t2 of slot 2"
        inequalities = {
            "cos(2π angle1) + c cos(2π angle0)": f"cos(π(1 + k1)) + c2 cos(π(1 + k0)) must be positive, {slots}",
            "cos(2π angle0) + c cos(2π angle1)": f"cos(π(1 + k0)) + c2 cos(π(1 + k1)) must be positive, {slots}",
        }
    else:
        constructor, arguments, seams = (
            hyperbolic_polygons.triangle,
            [angles[2], angles[0], angles[1]],
            [0, 1, 2],
        )
        inequalities = {}

    try:
        P = constructor(*arguments, ring=AA)
    except ValueError as error:
        for inequality, message in inequalities.items():
            if inequality in str(error):
                raise ValueError(
                    f"no pair of pants with these parameters: {message}"
                ) from error
        raise

    vertices = [v.coordinates(model="klein") for v in P.vertices()]
    if optimize_placement:
        vertices = _optimize_placement(vertices)

    return vertices, seams


def _mirror(P, glued_edges, mirror_edge):
    r"""
    Return the mirror image of the hyperbolic polygon ``P`` and the gluings
    of :meth:`HyperbolicSurfaceGenerators.pillowcase`.

    INPUT:

    - ``P``, ``glued_edges``, ``mirror_edge`` -- as in
      :meth:`HyperbolicSurfaceGenerators.pillowcase`

    OUTPUT:

    A pair ``(Q, gluings)``: the reflection `Q` of ``P`` in the line of the
    mirror edge, a polygon over the base ring of ``P``, and the list of
    the pairs of edges ``(label, edge)`` to glue, ``P`` having the label 0
    and `Q` the label 1.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_surface_generators import _mirror
        sage: P = HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/3, 1/2)])
        sage: Q, gluings = _mirror(P, None, 1)
        sage: [v.coordinates(model="klein") for v in Q.vertices()]
        [(36/49, 12/49), (1/3, 1/2), (1/2, 0)]
        sage: gluings
        [((0, 0), (1, 2)), ((0, 1), (1, 1)), ((0, 2), (1, 0))]

    """
    from flatsurf.geometry.similarity_surface_generators import _pillowcase_gluings

    vertices = [v.coordinates(model="klein") for v in P.vertices()]
    m, permutation, gluings = _pillowcase_gluings(
        len(vertices), glued_edges, mirror_edge
    )

    reflection = _reflection(vertices[m], vertices[(m + 1) % len(vertices)])

    mirrored = []
    for j in permutation:
        x, y, z = reflection * vector((vertices[j][0], vertices[j][1], 1))
        mirrored.append((x / z, y / z))
    return P.parent().polygon(vertices=mirrored), gluings


class HyperbolicSurfaceGenerators:
    r"""
    Constructors of hyperbolic surfaces from hyperbolic polygons.

    EXAMPLES::

        sage: from flatsurf import hyperbolic_surfaces
        sage: hyperbolic_surfaces
        hyperbolic surfaces

    """

    def __repr__(self):
        r"""
        Return a printable representation of this constructor.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_surfaces
            sage: hyperbolic_surfaces
            hyperbolic surfaces

        """
        return "hyperbolic surfaces"

    @staticmethod
    def pillowcase(P, glued_edges=None, mirror_edge=None):
        r"""
        Return the surface built from the hyperbolic polygon ``P`` and its
        mirror image, glued along ``glued_edges``.

        The polygon ``P`` has label 0 and its reflection `Q` in the line of
        the edge `m` = ``mirror_edge`` has label 1. With `n` the number of
        vertices of ``P``, the vertex `j` of `Q` is the reflection of the
        vertex `2m + 1 - j \bmod n` of ``P``, so that the edge `j` of `Q` is
        the reflection of the edge `2m - j \bmod n` of ``P``. The edge `e`
        of ``P`` is glued to the edge `2m - e \bmod n` of `Q` for every `e`
        in ``glued_edges``, by the reflection in the line of `e` followed by
        the reflection in the line of `m`; the edge `m` is glued to itself
        in the plane, by the identity. The other edges of ``P`` and of `Q`
        form the boundary of the surface.

        INPUT:

        - ``P`` -- a
          :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`,
          convex or not

        - ``glued_edges`` -- an iterable of edges of ``P`` or ``None``
          (default: ``None``) for all the edges of ``P``

        - ``mirror_edge`` -- an edge in ``glued_edges`` or ``None`` (default:
          ``None``) for the first of ``glued_edges``, the edge `0` when
          ``glued_edges`` is ``None``

        OUTPUT:

        An immutable
        :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`
        over the base ring of ``P``: a reflection in a line through two
        points has rational coefficients in their Klein coordinates.

        EXAMPLES:

        A triangle reflected in its edge 1, a line that does not pass through
        the origin::

            sage: from flatsurf import HyperbolicPlane, hyperbolic_surfaces
            sage: P = HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/3, 1/2)])
            sage: S = hyperbolic_surfaces.pillowcase(P, mirror_edge=1)
            sage: S
            Genus 0 Hyperbolic Surface built from 2 triangles
            sage: S.base_ring()
            Rational Field
            sage: [v.coordinates(model="klein") for v in S.polygon(1).vertices()]
            [(36/49, 12/49), (1/3, 1/2), (1/2, 0)]

        Each gluing is the reflection in the line of its edge followed by
        the reflection in the line of the edge 1; the edge 1 is glued by the
        identity::

            sage: from flatsurf.geometry.hyperbolic_surface_generators import _reflection
            sage: vertices = [v.coordinates(model="klein") for v in P.vertices()]
            sage: R = [_reflection(vertices[e], vertices[(e + 1) % 3]) for e in range(3)]
            sage: S.edge_transformation(0, 0)
            [-41/31  24/31  36/31]
            [-24/31 -23/31  12/31]
            [-36/31  12/31  49/31]
            sage: S.edge_transformation(0, 1)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: S.edge_transformation(0, 2)
            [ -83/403 -612/403    36/31]
            [ 396/403 -173/403    12/31]
            [  36/403 -492/403    49/31]
            sage: all(S.edge_transformation(0, e).matrix() == R[1] * R[e] for e in range(3))
            True
            sage: S.holonomy_representation().invariant_trace_field()
            Rational Field
            sage: TestSuite(S).run()

        A non-convex quadrilateral is not triangulated; its reflex vertex
        is a cone point of angle `3/2`::

            sage: P = HyperbolicPlane(QQ).polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
            sage: S = hyperbolic_surfaces.pillowcase(P)
            sage: S
            Genus 0 Hyperbolic Surface built from 2 non-convex quadrilaterals
            sage: S(0, 2).angle()
            3/2
            sage: TestSuite(S).run()

        A right-angled hexagon glued along its alternate edges is a pair of
        pants::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 145/24)
            sage: S = hyperbolic_surfaces.pillowcase(P, glued_edges=[0, 2, 4])
            sage: S
            Hyperbolic Surface with boundary built from 2 right-angled hexagons
            sage: S.base_ring()
            Rational Field
            sage: S.boundary_components()
            [((0, 1), (1, 5)), ((0, 3), (1, 3)), ((0, 5), (1, 1))]
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(21/4, 21/4, 20881/144))
            sage: TestSuite(S).run()

        The mirror edge must be glued, and the glued edges distinct edges of
        ``P``::

            sage: hyperbolic_surfaces.pillowcase(P, glued_edges=[])
            Traceback (most recent call last):
            ...
            ValueError: glued_edges must not be empty, the two polygons would be disjoint
            sage: hyperbolic_surfaces.pillowcase(P, glued_edges=[0, 2], mirror_edge=1)
            Traceback (most recent call last):
            ...
            ValueError: mirror_edge 1 is not one of the glued edges
            sage: hyperbolic_surfaces.pillowcase(P, glued_edges=[0, 2, 0])
            Traceback (most recent call last):
            ...
            ValueError: edge 0 is repeated in glued_edges
            sage: hyperbolic_surfaces.pillowcase(P, glued_edges=[6])
            Traceback (most recent call last):
            ...
            ValueError: 6 is not an edge of a polygon with 6 vertices

        A Euclidean polygon is refused, see
        :meth:`~flatsurf.geometry.similarity_surface_generators.SimilaritySurfaceGenerators.pillowcase`
        instead::

            sage: from flatsurf import polygons
            sage: hyperbolic_surfaces.pillowcase(polygons.square())
            Traceback (most recent call last):
            ...
            TypeError: P must be a hyperbolic polygon

        """
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon
        from flatsurf.geometry.hyperbolic_surface import (
            MutableOrientedHyperbolicSurface,
        )

        if not isinstance(P, HyperbolicPolygon):
            raise TypeError("P must be a hyperbolic polygon")

        Q, gluings = _mirror(P, glued_edges, mirror_edge)

        surface = MutableOrientedHyperbolicSurface(P.parent().base_ring())
        surface.add_polygon(P, label=0)
        surface.add_polygon(Q, label=1)
        for x, y in gluings:
            surface.glue(x, y)
        surface.set_immutable()
        return surface

    @staticmethod
    def pants(boundary_traces=(), cone_orders=(), ring=None, optimize_placement=True):
        r"""
        Return the hyperbolic pair of pants with geodesic boundary components
        of traces ``boundary_traces`` and cone points of orders
        ``cone_orders``.

        A pair of pants has three slots, each a geodesic boundary component,
        a cuff, or a cone point. The cone points take the first slots, in the
        order of ``cone_orders``, and the cuffs the last ones, in the order
        of ``boundary_traces``. The pants is the
        :meth:`pillowcase` of its front polygon glued along its seams, the
        mirror edge being the seam between the slots 2 and 0. The front
        polygon has the edges, counterclockwise, seam 20, half-cuff 0, seam
        01, half-cuff 1, seam 12, half-cuff 2, all its angles right; a cone
        point in the slot `s` replaces the half-cuff `s` by a vertex of half
        its cone angle.

        INPUT:

        - ``boundary_traces`` -- a sequence of real algebraic numbers larger
          than 3 (default: ``()``); the traces in `\mathrm{SO}(2, 1)` of the
          holonomy of the cuffs. A cuff of trace `t` has length `\ell` with
          `\cosh(\ell / 2) = \sqrt{t + 1} / 2`.

        - ``cone_orders`` -- a sequence of rational numbers in `(-1, 0)`
          (default: ``()``); the orders `k` of the cone points, of angle
          `2\pi(1 + k)`

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the surface: with ``None``, the number field generated by the
          Klein coordinates of the front polygon, ``QQ`` when they are
          rational; ``AA`` keeps the coordinates in ``AA``; any other ring
          receives the coordinates by conversion

        - ``optimize_placement`` -- a boolean (default: ``True``); whether to
          move the front polygon by an isometry that makes its coordinates
          rational when they generate a quadratic field, see below.
          Otherwise the front polygon has its vertex 0 at the origin and its
          edge 0 on the positive `x`-axis.

        OUTPUT:

        An immutable
        :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`
        with the labels 0, the front polygon, and 1, its mirror image. Its
        :meth:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces.FiniteType.Oriented.Connected.ParentMethods.signature`
        has the genus 0, the cone orders ``cone_orders`` and the boundary
        traces ``boundary_traces``.

        ALGORITHM:

        The front polygon is
        :meth:`~HyperbolicPolygonsConstructor.right_angled_hexagon` without
        cone point,
        :meth:`~HyperbolicPolygonsConstructor.right_angled_pentagon` with
        one, :meth:`~HyperbolicPolygonsConstructor.right_angled_quadrilateral`
        with two, and :meth:`~HyperbolicPolygonsConstructor.triangle` with
        three, built over ``AA``.

        When its coordinates generate `\QQ(\sqrt{d})` and
        ``optimize_placement`` is set, the polygon is moved by the isometry
        that sends one of its vertices to a point `(x_0, 0)` with `1 -
        x_0^2 \in d \QQ^2`, and the direction of a neighbouring vertex to
        the positive `x`-axis, for the first vertex and neighbour whose image
        has rational coordinates; the polygon is kept when there is none.

        EXAMPLES:

        The pants with three cuffs of trace `21/4` is defined over `\QQ`
        once its front hexagon is moved off the origin::

            sage: from flatsurf import hyperbolic_surfaces
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4))
            sage: S
            Hyperbolic Surface with boundary built from 2 right-angled hexagons
            sage: S.base_ring()
            Rational Field
            sage: S.polygon(0).vertices()
            ((-1/5, 0),
             (-1/5, -24/25),
             (-5/61, -60/61),
             (23/65, -12/13),
             (5/11, -48/55),
             (5/11, 0))
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(21/4, 21/4, 21/4))
            sage: S.area()
            1
            sage: TestSuite(S).run()

        With the vertex 0 of the hexagon at the origin, the coordinates need
        `\sqrt{6}`::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), optimize_placement=False)
            sage: S.base_ring()
            Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(21/4, 21/4, 21/4))
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), ring=QQ, optimize_placement=False)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field

        A cone point of angle `\pi` and two cuffs, from a pentagon::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(QQ(21)/4, QQ(21)/4), cone_orders=(QQ(-1)/2,))
            sage: S
            Hyperbolic Surface with boundary built from 2 right-angled pentagons
            sage: S.base_ring()
            Rational Field
            sage: S.polygon(0).vertices()
            ((-15/19, 0),
             -15/29 + 6/29*I,
             (-225/421, -340/421),
             (-9/25, -544/625),
             (-9/25, 0))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-1/2,), boundary_traces=(21/4, 21/4))
            sage: S.area()
            1/2

        Two cone points and a cuff, from a quadrilateral::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4,), cone_orders=(-1/2, -2/3))
            sage: S
            Hyperbolic Surface with boundary built from 2 quadrilaterals
            sage: S.base_ring()
            Rational Field
            sage: S.polygon(0).vertices()
            ((-6/7, 0), -2/3 + 1/3*I, (-9/17, -10/17), (-9/17, 0))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-2/3, -1/2), boundary_traces=(21/4,))
            sage: S.area()
            1/6

        Three cone points, from a triangle, give a closed orbifold of genus
        0::

            sage: S = hyperbolic_surfaces.pants(cone_orders=(-1/2, -3/4, -7/8))
            sage: S
            Genus 0 Hyperbolic Surface built from 2 right triangles
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 2 with a = 1.189207115002722?
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-7/8, -3/4, -1/2))
            sage: S.area()
            1/8

        A cone angle larger than `\pi` requires cuffs of comparable lengths::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(-1/4,))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-1/4,), boundary_traces=(21/4, 21/4))
            sage: S.area()
            1/4
            sage: hyperbolic_surfaces.pants(boundary_traces=(20881/144, 21/4), cone_orders=(-1/4,))
            Traceback (most recent call last):
            ...
            ValueError: no pair of pants with these parameters: c2 + c1 cos(π(1 + k0)) must be positive, with k0 the cone order of slot 0 and cs = √(ts + 1)/2 for the boundary trace ts of slot s

        Two cone points of angle `\pi` and a cuff would give a flat
        surface::

            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4,), cone_orders=(-1/2, -1/2))
            Traceback (most recent call last):
            ...
            ValueError: the divisor Euler characteristic 2 + sum(cone_orders) - len(boundary_traces) must be negative, not 0

        A pair of pants has three slots, the traces of its cuffs are larger
        than 3, and its cone orders are in `(-1, 0)`::

            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4))
            Traceback (most recent call last):
            ...
            ValueError: a pair of pants has three cuffs and cone points in total, not 2
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 3))
            Traceback (most recent call last):
            ...
            ValueError: the boundary trace 3 must be larger than 3
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(-1,))
            Traceback (most recent call last):
            ...
            ValueError: the cone order -1 is not larger than -1
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(0,))
            Traceback (most recent call last):
            ...
            ValueError: a cone order must be nonzero
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(1,))
            Traceback (most recent call last):
            ...
            NotImplementedError: cone points of angle larger than 2π are not implemented, the cone order 1 is positive

        TESTS:

        The slots of the cuffs and of the cone points, in their input
        order::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(int(4), int(4), int(4)))
            sage: S.base_ring().polynomial()
            y^4 - 2*y^2 - 4
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(4, 4, 4))
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 20881/144))
            sage: S.base_ring()
            Rational Field
            sage: [v.coordinates(model="klein") for v in S.polygon(0).vertices()]
            [(0, 0), (56/65, 0), (56/65, 99/325), (840/4121, 4015/4121), (1344/9425, 143/145), (0, 143/145)]
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(-2/3,)).base_ring()
            Rational Field
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4,), cone_orders=(-2/3, -2/3)).base_ring()
            Rational Field
            sage: hyperbolic_surfaces.pants(boundary_traces=(20881/144,), cone_orders=(-1/4, -7/8))
            Traceback (most recent call last):
            ...
            ValueError: no pair of pants with these parameters: cos(π(1 + k1)) + c2 cos(π(1 + k0)) must be positive, with k0, k1 the cone orders of slots 0, 1 and c2 = √(t2 + 1)/2 for the boundary trace t2 of slot 2
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(1/2,))
            Traceback (most recent call last):
            ...
            NotImplementedError: cone points of angle larger than 2π are not implemented, the cone order 1/2 is positive
            sage: hyperbolic_surfaces.pants(cone_orders=(-1/2, -2/3, -5/6))
            Traceback (most recent call last):
            ...
            ValueError: the divisor Euler characteristic 2 + sum(cone_orders) - len(boundary_traces) must be negative, not 0
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(sqrt(2),))
            Traceback (most recent call last):
            ...
            ValueError: the cone order sqrt(2) must be a rational number
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, I))
            Traceback (most recent call last):
            ...
            ValueError: the boundary trace I must be a real algebraic number
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), ring=AA)
            sage: S.base_ring()
            Algebraic Real Field

        """
        vertices, seams = _pants_polygon(
            boundary_traces, cone_orders, optimize_placement
        )
        return HyperbolicSurfaceGenerators.pillowcase(
            _polygon(vertices, ring), glued_edges=seams
        )

    @staticmethod
    def regular_octagon(gluing="opposite", ring=None):
        r"""
        Return the hyperbolic surface of genus 2 built from the regular
        octagon with interior angles `\pi/4`.

        INPUT:

        - ``gluing`` -- a string (default: ``"opposite"``); ``"opposite"``
          glues each edge `e` to the edge `e + 4`, ``"katok"`` glues the
          edges `0, 1, 4, 5` to the edges `2, 3, 6, 7`, the identification
          `a_1 b_1 a_1^{-1} b_1^{-1} a_2 b_2 a_2^{-1} b_2^{-1}` of Katok,
          Fuchsian groups (1992), Section 4.3, Example C, p. 95

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the octagon, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        An immutable
        :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`
        with the label 0, the octagon of
        :meth:`~HyperbolicPolygonsConstructor.regular_ngon`, with its vertex
        0 at the origin and its edge 0 on the positive `x`-axis.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_surfaces
            sage: S = hyperbolic_surfaces.regular_octagon()
            sage: S
            Genus 2 Hyperbolic Surface built from a regular octagon
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 2*y^2 - 1 with a = ...
            sage: S.signature()
            SurfaceSignature(g=2)
            sage: S.area()
            2
            sage: S.holonomy_representation().invariant_trace_field()
            Number Field in ... with defining polynomial x^2 - 2 with ...
            sage: TestSuite(S).run()

        The identification of Katok::

            sage: S = hyperbolic_surfaces.regular_octagon(gluing="katok")
            sage: S
            Genus 2 Hyperbolic Surface built from a regular octagon
            sage: S.signature()
            SurfaceSignature(g=2)
            sage: S.area()
            2
            sage: S.holonomy_representation().invariant_trace_field()
            Number Field in ... with defining polynomial x^2 - 2 with ...
            sage: TestSuite(S).run()

        TESTS::

            sage: hyperbolic_surfaces.regular_octagon(gluing="adjacent")
            Traceback (most recent call last):
            ...
            ValueError: gluing must be "opposite" or "katok", not 'adjacent'
            sage: hyperbolic_surfaces.regular_octagon(ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field

        """
        from flatsurf.geometry.hyperbolic_surface import (
            MutableOrientedHyperbolicSurface,
        )

        if gluing == "opposite":
            pairs = [(0, 4), (1, 5), (2, 6), (3, 7)]
        elif gluing == "katok":
            pairs = [(0, 2), (1, 3), (4, 6), (5, 7)]
        else:
            raise ValueError(f'gluing must be "opposite" or "katok", not {gluing!r}')

        P = hyperbolic_polygons.regular_ngon(8, QQ((1, 8)), ring=ring)

        surface = MutableOrientedHyperbolicSurface(P.base_ring())
        surface.add_polygon(P, label=0)
        for e, f in pairs:
            surface.glue((0, e), (0, f))
        surface.set_immutable()
        return surface

    @staticmethod
    def from_signature(signature, cuff_traces=None, ring=None, optimize_placement=True):
        r"""
        Return a compact hyperbolic surface of signature ``signature``, glued
        from pairs of pants.

        INPUT:

        - ``signature`` -- a
          :class:`~flatsurf.geometry.signature.SurfaceSignature` without
          cusps, whose cone orders are in `(-1, 0)`, whose boundary traces
          are larger than 3, and whose geometry is hyperbolic

        - ``cuff_traces`` -- a sequence of real algebraic numbers larger than
          3 or ``None`` (default: ``None``); the traces in `\mathrm{SO}(2,
          1)` of the interior cuffs, in the order given below. With
          ``None``, every interior cuff has the trace `21/4`.

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the surface: with ``None``, the number field generated by the
          Klein coordinates of all the polygons, ``QQ`` when they are
          rational; ``AA`` keeps the coordinates in ``AA``; any other ring
          receives the coordinates by conversion

        - ``optimize_placement`` -- a boolean (default: ``True``); whether to
          move each polygon to coordinates of small degree, as in
          :meth:`pants` and :meth:`~FuchsianGroupGenerators.triangle`

        OUTPUT:

        An immutable
        :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`
        whose
        :meth:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces.FiniteType.Oriented.Connected.ParentMethods.signature`
        is ``signature``. The pair of pants `i` of the decomposition below
        is made of the polygons with the labels ``(i, 0)``, its front
        polygon, and ``(i, 1)``, its mirror image, as in :meth:`pants`.

        EXAMPLES:

        A closed surface of genus 2, from two pairs of pants::

            sage: from flatsurf import hyperbolic_surfaces
            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature(2))
            sage: S
            Genus 2 Hyperbolic Surface built from 4 right-angled hexagons
            sage: S.signature()
            SurfaceSignature(g=2)
            sage: S.base_ring()
            Rational Field
            sage: S.area()
            2
            sage: TestSuite(S).run()

        Orbifolds, given in Katok's form; the cone points of angle `\pi`
        come in pairs as folds::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(1, (2,)))
            sage: S.signature()
            SurfaceSignature(g=1, cone_orders=(-1/2,))
            sage: S.base_ring()
            Rational Field
            sage: S.area()
            1/2
            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (2, 2, 2, 3)))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-2/3, -1/2, -1/2, -1/2))
            sage: S.base_ring()
            Rational Field
            sage: S.area()
            1/6

        Three cone points of angles `2\pi/m` give the pillowcase of the
        triangle of the triangle group, placed by the search of
        :meth:`~FuchsianGroupGenerators.triangle`::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (3, 3, 4)))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-3/4, -2/3, -2/3))
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 2 with a = 1.189207115002722?
            sage: S.area()
            1/12

        ::

            sage: # long time
            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (2, 3, 7)))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-6/7, -2/3, -1/2))
            sage: S.base_ring().polynomial()
            y^6 + y^4 - 2*y^2 - 1
            sage: S.area()
            1/42

        With the default traces, the polygons of the orbifold of signature
        `(2; 2)` are rational but not its side pairings: the surface is
        glued over ``AA`` and moved to the field of its coordinates and
        side pairings, which is its invariant trace field. Traces `99`
        give a surface over ``QQ``::

            sage: s = SurfaceSignature.from_orbifold_orders(2, (2,))
            sage: S = hyperbolic_surfaces.from_signature(s)
            sage: S.signature()
            SurfaceSignature(g=2, cone_orders=(-1/2,))
            sage: S.base_ring()
            Number Field in a with defining polynomial y^2 - 51 with a = 7.141428428542850?
            sage: S.area()
            5/2
            sage: S.holonomy_representation().invariant_trace_field()
            Number Field in a with defining polynomial y^2 - 51 with a = 7.141428428542850?
            sage: S = hyperbolic_surfaces.from_signature(s, cuff_traces=(int(99),) * 4)
            sage: S.signature()
            SurfaceSignature(g=2, cone_orders=(-1/2,))
            sage: S.base_ring()
            Rational Field
            sage: S.area()
            5/2

        A surface with geodesic boundary::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature(1, boundary_traces=(QQ(21)/4,)))
            sage: S.signature()
            SurfaceSignature(g=1, boundary_traces=(21/4,))
            sage: S.base_ring()
            Rational Field
            sage: S.area()
            1

        Cusps, cone angles larger than `2\pi`, boundary components without
        a trace larger than 3, and signatures that are not hyperbolic are
        refused::

            sage: hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (2, 3), cusps=1))
            Traceback (most recent call last):
            ...
            NotImplementedError: surfaces with cusps are not implemented, the signature has 1 cusps
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2, cone_orders=(1/2,)))
            Traceback (most recent call last):
            ...
            NotImplementedError: cone points of angle larger than 2π are not implemented, the cone order 1/2 is positive
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(1, boundary_traces=(3,)))
            Traceback (most recent call last):
            ...
            ValueError: the boundary trace 3 must be larger than 3
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (2, 3, 6)))
            Traceback (most recent call last):
            ...
            ValueError: the signature SurfaceSignature(g=0, cone_orders=(-5/6, -2/3, -1/2)) is euclidean, not hyperbolic

        A pair of pants of the decomposition may not exist, while a surface
        of the signature does::

            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(0, cone_orders=(-1/4,), boundary_traces=(21/4, 20881/144)))
            Traceback (most recent call last):
            ...
            NotImplementedError: a pair of pants of the decomposition into pants of this signature does not exist, the surface may exist without such a decomposition: no pair of pants with these parameters: ...

        TESTS:

        The labels of the closed surface of genus 3, from four pairs of
        pants::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature(3))
            sage: list(S.labels())
            [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
            sage: S.signature()
            SurfaceSignature(g=3)

        The number of cuff traces, their values, and an explicit ring
        that does not hold the side pairings::

            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(21/4,))
            Traceback (most recent call last):
            ...
            ValueError: cuff_traces must have 3 traces, one per interior cuff, not 1
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(21/4, 21/4, 2))
            Traceback (most recent call last):
            ...
            ValueError: the cuff trace 2 must be larger than 3
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(2, (2,)), ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: cannot glue edge 5 of polygon (0, 0) to edge 2 of polygon (1, 0): isometry is not defined over Rational Field; adjoin a square root of 390625/55539
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), ring=AA).base_ring()
            Algebraic Real Field
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(1, boundary_traces=(None,)))
            Traceback (most recent call last):
            ...
            ValueError: the boundary traces of the signature must be given, not None
            sage: hyperbolic_surfaces.from_signature((2,))
            Traceback (most recent call last):
            ...
            TypeError: signature must be a SurfaceSignature

        ALGORITHM:

        The legs of the signature are its boundary components, its folds,
        one for each pair of cone points of angle `\pi`, and its other cone
        points. A fold is a cuff whose front half-cuff is glued to its own
        back half-cuff, the edge `e` of the polygon ``(i, 0)`` to the edge
        `-e \bmod n` of the polygon ``(i, 1)``; its two seam feet become cone
        points of angle `\pi`. The legs are ordered boundary components
        first, then folds, then the cone points sorted by angle, paired
        smallest with largest, `k_0, k_{c-1}, k_1, k_{c-2}, \dots`, so that
        the pants at the ends of the chain below, the only ones with two
        legs, receive them first.

        With `L_1, \dots, L_n` the legs and `g` the genus, the `2g - 2 + n`
        pairs of pants are glued along a chain:

        - `g = 0`: the pants `(L_1, L_2, c_1)`, then `(c_{i-1}, L_{i+1},
          c_i)` for `i = 2, \dots, n - 3`, then `(c_{n-3}, L_{n-1}, L_n)`;
          the single pants `(L_1, L_2, L_3)` when `n = 3`;

        - `g \ge 1`: the pants `(\ell, \ell, c)` whose first two cuffs, a
          loop, are glued to each other, or `(\ell, \ell, L_1)` when `g = 1`
          and `n = 1`; then one pants `(c, L_i, c')` per leg, but for the
          last two when `g = 1`; then, for each handle `2, \dots, g - 1`, two
          pants `(c, d, d')` and `(d, d', c')` glued along two cuffs; then
          `(c, \ell', \ell')` when `g \ge 2`, `(c, L_{n-1}, L_n)` when `g =
          1`.

        A pair of pants is the one of :meth:`pants`, its cone points in its
        first slots, its cuffs in its last ones. Two cuffs are glued without
        twist: gluing the slot `a` of the pants `i` to the slot `b` of the
        pants `j` glues the front half-cuffs ``((i, 0), e_a)`` and ``((j, 0),
        e_b)`` and the back ones ``((i, 1), -e_a mod n_i)`` and ``((j, 1),
        -e_b mod n_j)``, so that the reflection exchanging every front and
        back polygon is an isometry of the surface.

        The interior cuffs are the loops, the links `c`, the doubled cuffs
        `d, d'` and the folds. Each has the trace `21/4` by default, for
        which the pants have rational coordinates once placed, or the trace
        of ``cuff_traces``, given in the order in which the interior cuffs
        first appear along the pants of the chain and their slots: on `(2;
        2)` the loop `\ell`, the links `c_1, c_2`, the loop `\ell'`.

        With ``ring=None``, the polygons are built over the number field of
        all their Klein coordinates. The side pairings between two pants
        may need square roots that this field does not hold: when a gluing
        is refused, the surface is glued over ``AA`` and moved to the field
        of its coordinates and side pairings by
        :meth:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface.set_immutable`.
        With an explicit ``ring``, the refusal is raised.

        Three cone points of angles `2\pi/m_i`, without boundary, give the
        pillowcase of the triangle of
        :meth:`~FuchsianGroupGenerators.triangle`.

        """
        from flatsurf.geometry.hyperbolic import HyperbolicPlane
        from flatsurf.geometry.hyperbolic_surface import (
            MutableOrientedHyperbolicSurface,
        )
        from flatsurf.geometry.signature import SurfaceSignature

        if not isinstance(signature, SurfaceSignature):
            raise TypeError("signature must be a SurfaceSignature")

        def check_trace(t, name):
            try:
                trace = AA(t)
            except (TypeError, ValueError):
                raise ValueError(f"the {name} {t} must be a real algebraic number")
            if not trace > 3:
                raise ValueError(f"the {name} {t} must be larger than 3")

        if signature.cusps():
            raise NotImplementedError(
                f"surfaces with cusps are not implemented, the signature has {signature.cusps()} cusps"
            )
        for k in signature.cone_orders():
            if k > 0:
                raise NotImplementedError(
                    f"cone points of angle larger than 2π are not implemented, the cone order {k} is positive"
                )
        for t in signature.boundary_traces():
            if t is None:
                raise ValueError(
                    "the boundary traces of the signature must be given, not None"
                )
            check_trace(t, "boundary trace")
        geometry = signature.geometry()
        if geometry != "hyperbolic":
            raise ValueError(
                f"the signature {signature!r} is {geometry}, not hyperbolic"
            )

        genus = signature.genus()
        cone_orders = signature.cone_orders()
        boundary_traces = signature.boundary_traces()

        triangle = (
            genus == 0
            and not boundary_traces
            and len(cone_orders) == 3
            and all((1 + k).numerator() == 1 for k in cone_orders)
        )

        # A leg is ("boundary", trace), ("fold", None) or ("cone", order).
        half = QQ((-1, 2))
        cones = list(cone_orders)
        folds = 0 if triangle else cones.count(half) // 2
        for _ in range(2 * folds):
            cones.remove(half)
        paired = []
        while cones:
            paired.append(cones.pop(0))
            if cones:
                paired.append(cones.pop())
        legs = (
            [("leg", ("boundary", t)) for t in boundary_traces]
            + [("leg", ("fold", None))] * folds
            + [("leg", ("cone", k)) for k in paired]
        )

        # Each slot of a pants is a leg or ("cuff", name), each name in two
        # slots.
        n = len(legs)
        if genus == 0 and n == 3:
            chain = [legs]
        elif genus == 0:
            chain = [[legs[0], legs[1], ("cuff", ("link", 1))]]
            for i in range(2, n - 2):
                chain.append(
                    [("cuff", ("link", i - 1)), legs[i], ("cuff", ("link", i))]
                )
            chain.append([("cuff", ("link", n - 3)), legs[n - 2], legs[n - 1]])
        elif genus == 1 and n == 1:
            loop = ("cuff", ("loop", 0))
            chain = [[loop, loop, legs[0]]]
        else:
            loop = ("cuff", ("loop", 0))
            chain = [[loop, loop, ("cuff", ("link", 0))]]
            link = 0
            for leg in legs if genus >= 2 else legs[: n - 2]:
                chain.append(
                    [("cuff", ("link", link)), leg, ("cuff", ("link", link + 1))]
                )
                link += 1
            for handle in range(2, genus):
                d, e = ("cuff", ("handle", handle, 0)), ("cuff", ("handle", handle, 1))
                chain.append([("cuff", ("link", link)), d, e])
                chain.append([d, e, ("cuff", ("link", link + 1))])
                link += 1
            if genus >= 2:
                loop = ("cuff", ("loop", 1))
                chain.append([("cuff", ("link", link)), loop, loop])
            else:
                chain.append([("cuff", ("link", link)), legs[n - 2], legs[n - 1]])

        # The interior cuffs in the order they first appear; a fold is
        # identified by its slot.
        interior = []
        for i, slots in enumerate(chain):
            for s, slot in enumerate(slots):
                if slot[0] == "cuff":
                    key = slot[1]
                elif slot[1][0] == "fold":
                    key = ("fold", i, s)
                else:
                    continue
                if key not in interior:
                    interior.append(key)

        if cuff_traces is None:
            cuff_traces = [QQ((21, 4))] * len(interior)
        else:
            cuff_traces = list(cuff_traces)
            if len(cuff_traces) != len(interior):
                raise ValueError(
                    f"cuff_traces must have {len(interior)} traces, one per interior cuff, not {len(cuff_traces)}"
                )
            for t in cuff_traces:
                check_trace(t, "cuff trace")
        traces = dict(zip(interior, cuff_traces))

        # The half-cuff edge of the front polygon of each slot that is not a
        # cone point.
        half_cuffs = {}
        if triangle:
            periods = signature.orbifold_orders()[1]
            polygons = [_triangle_group_polygon(*periods, ring, optimize_placement)]
            seams = [None]
        else:
            fronts = []
            seams = []
            for i, slots in enumerate(chain):
                orders = []
                cuffs = []
                pants_traces = []
                for s, slot in enumerate(slots):
                    if slot[0] == "cuff":
                        pants_traces.append(traces[slot[1]])
                    elif slot[1][0] == "boundary":
                        pants_traces.append(slot[1][1])
                    elif slot[1][0] == "fold":
                        pants_traces.append(traces[("fold", i, s)])
                    else:
                        orders.append(slot[1][1])
                        continue
                    cuffs.append(s)
                try:
                    vertices, pants_seams = _pants_polygon(
                        pants_traces, orders, optimize_placement
                    )
                except ValueError as error:
                    raise NotImplementedError(
                        f"a pair of pants of the decomposition into pants of this signature does not exist, the surface may exist without such a decomposition: {error}"
                    ) from error
                edges = [e for e in range(len(vertices)) if e not in pants_seams]
                for s, e in zip(cuffs, edges):
                    half_cuffs[i, s] = e
                fronts.append(vertices)
                seams.append(pants_seams)

            if ring is None:
                from flatsurf.geometry.subfield import (
                    number_field_elements_from_algebraics,
                )

                field, coordinates = number_field_elements_from_algebraics(
                    [c for vertices in fronts for vertex in vertices for c in vertex]
                )
                H = HyperbolicPlane(field)
                coordinates = iter(coordinates)
                polygons = [
                    H.polygon(
                        vertices=[
                            (next(coordinates), next(coordinates)) for _ in vertices
                        ]
                    )
                    for vertices in fronts
                ]
            else:
                polygons = [_polygon(vertices, ring) for vertices in fronts]

        labelled = {}
        gluings = []
        for i, P in enumerate(polygons):
            Q, pants_gluings = _mirror(P, seams[i], None)
            labelled[i, 0] = P
            labelled[i, 1] = Q
            gluings.extend(
                (((i, a), e), ((i, b), f)) for (a, e), (b, f) in pants_gluings
            )

        ends = {}
        for i, slots in enumerate(chain):
            size = len(polygons[i].vertices())
            for s, slot in enumerate(slots):
                if (i, s) not in half_cuffs:
                    continue
                e = half_cuffs[i, s]
                if slot[0] == "cuff":
                    ends.setdefault(slot[1], []).append((i, e, size))
                elif slot[1][0] == "fold":
                    gluings.append((((i, 0), e), ((i, 1), -e % size)))
        for (i, a, n_i), (j, b, n_j) in ends.values():
            gluings.append((((i, 0), a), ((j, 0), b)))
            gluings.append((((i, 1), -a % n_i), ((j, 1), -b % n_j)))

        def assemble(base_ring):
            surface = MutableOrientedHyperbolicSurface(base_ring)
            for label, polygon in labelled.items():
                surface.add_polygon(polygon, label=label)
            for x, y in gluings:
                surface.glue(x, y)
            return surface

        try:
            surface = assemble(polygons[0].parent().base_ring())
        except ValueError as error:
            if ring is not None or "isometry is not defined over" not in str(error):
                raise
            surface = assemble(AA)
            surface.set_immutable(minimal_field=True)
        else:
            surface.set_immutable()
        return surface


hyperbolic_surfaces = HyperbolicSurfaceGenerators()


def _triangle_group_degree_bound(m1, m2, m3):
    r"""
    Return a lower bound on the degree of the real fields that contain the
    Klein coordinates of a triangle of the triangle group with periods
    ``m1``, ``m2``, ``m3``.

    INPUT:

    - ``m1``, ``m2``, ``m3`` -- integers at least 2, with `1/m_1 + 1/m_2 +
      1/m_3 < 1`

    OUTPUT:

    The integer `[k_0 : \QQ]` when the quaternion algebra `A_0` is a matrix
    algebra over `k_0`, `2 [k_0 : \QQ]` otherwise.

    ALGORITHM:

    With `c_i = \cos(\pi/m_i)`, the invariant trace field of the group `\Gamma`
    is `k_0 = \QQ(c_1^2, c_2^2, c_3^2, c_1 c_2 c_3)`, Takeuchi, Arithmetic
    triangle groups (1977), Propositions 4 and 5. For two rotations `g, h`
    of periods `m_i, m_j` other than 2, of traces `2 c_i, 2 c_j` in
    `\mathrm{SL}_2(\RR)`, its quaternion algebra is `A_0 = (a, b | k_0)` with
    `a = \mathrm{tr}^2 g (\mathrm{tr}^2 g - 4) = 4 c_i^2 (4 c_i^2 - 4)` and `b
    = \mathrm{tr}^2 g \mathrm{tr}^2 h (\mathrm{tr} [g, h] - 2) = 16 c_i^2
    c_j^2 \cdot 4 \Delta`, `\Delta = c_1^2 + c_2^2 + c_3^2 + 2 c_1 c_2 c_3 -
    1`, by Maclachlan and Reid, The Arithmetic of Hyperbolic 3-Manifolds
    (2003), Theorem 3.6.2 and the Fricke identity. A hyperbolic triangle
    group has at most one period 2.

    If the Klein coordinates of a triangle lie in a real field `F`, the
    pairings of its pillowcase lie over `F`, so `\Gamma` lies in the
    orthogonal group of the isotropic form `x^2 + y^2 - z^2` over `F`, which
    is `\mathrm{PGL}_2(F)`. The lifts to `\mathrm{SL}_2(\RR)` of the squares
    of `\Gamma` then lie in a conjugate of `M_2(F)`, so `F` contains `k_0`
    and splits `A_0`: `[F : \QQ] \ge 2 [k_0 : \QQ]` when `A_0` is a division
    algebra.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _triangle_group_degree_bound
        sage: _triangle_group_degree_bound(2, 3, 7)
        6
        sage: _triangle_group_degree_bound(3, 4, 4)
        2
        sage: _triangle_group_degree_bound(int(2), int(4), int(5))
        4

    """
    from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
    from flatsurf.geometry.subfield import number_field_elements_from_algebraics

    periods = (m1, m2, m3)
    c1, c2, c3 = (_cos_sin(QQ((1, 2 * m)))[0] for m in periods)
    k0, (s1, s2, s3, p) = number_field_elements_from_algebraics(
        [c1 * c1, c2 * c2, c3 * c3, c1 * c2 * c3]
    )
    squares = (s1, s2, s3)
    i, j = [k for k in range(3) if periods[k] != 2][:2]

    delta = s1 + s2 + s3 + 2 * p - 1
    a = 4 * squares[i] * (4 * squares[i] - 4)
    b = 16 * squares[i] * squares[j] * 4 * delta

    degree = k0.degree()
    if QuaternionAlgebra(k0, a, b).is_matrix_ring():
        return degree
    return 2 * degree


def _triangle_placements(angles, order):
    r"""
    Return an iterator over placements of the triangle with the interior
    angles ``angles`` at its vertices 0, 1, 2.

    INPUT:

    - ``angles`` -- a list of three positive rational numbers, the angles as
      multiples of `2\pi`, of sum less than `1/2`

    - ``order`` -- a sequence of vertices of the triangle

    OUTPUT:

    An iterator over lists of three pairs of elements of ``AA``, the Klein
    coordinates of the vertices 0, 1, 2 of the triangle, counterclockwise.
    For each vertex of ``order`` in turn, three placements with that vertex
    at the origin: the edge to the next vertex on the positive `x`-axis, the
    edge to the previous vertex on it, and the bisector of its angle on it.

    EXAMPLES:

    The triangle with a right angle at the vertex 0, placed at that vertex::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _triangle_placements
        sage: placements = list(_triangle_placements([1/4, 1/12, 1/12], [0]))
        sage: placements[0]
        [(0, 0), (0.8164965809277260?, 0), (0, 0.816496580927726?)]
        sage: placements[1]
        [(0, 0), (0, -0.8164965809277260?), (0.816496580927726?, 0)]
        sage: placements[2]
        [(0, 0),
         (0.5773502691896258?, -0.5773502691896258?),
         (0.5773502691896258?, 0.5773502691896258?)]

    """
    for i in order:
        T = hyperbolic_polygons.triangle(
            angles[i], angles[(i + 1) % 3], angles[(i + 2) % 3], ring=AA
        )
        vertices = [v.coordinates(model="klein") for v in T.vertices()]
        for rotation in (0, -angles[i], -angles[i] / 2):
            cos, sin = _cos_sin(rotation)
            rotated = [(cos * x - sin * y, sin * x + cos * y) for x, y in vertices]
            yield [rotated[(j - i) % 3] for j in range(3)]


def _triangle_group_polygon(m1, m2, m3, ring, optimize_placement):
    r"""
    Return the triangle of :meth:`FuchsianGroupGenerators.triangle`, with
    the vertices of periods ``m1``, ``m3``, ``m2`` at its vertices 0, 1, 2.

    INPUT:

    - ``m1``, ``m2``, ``m3``, ``ring``, ``optimize_placement`` -- as in
      :meth:`FuchsianGroupGenerators.triangle`, the group being hyperbolic

    OUTPUT:

    A hyperbolic triangle over ``ring``, over the field of its Klein
    coordinates when ``ring`` is ``None``, placed by the search of
    :meth:`FuchsianGroupGenerators.triangle` when ``optimize_placement`` is
    set.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _triangle_group_polygon
        sage: T = _triangle_group_polygon(3, 4, 4, None, True)
        sage: T.base_ring()
        Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?
        sage: T.angles()
        (1/6, 1/8, 1/8)

    """
    m1, m2, m3 = ZZ(m1), ZZ(m2), ZZ(m3)
    angles = [QQ((1, 2 * m)) for m in (m1, m3, m2)]

    if not optimize_placement:
        return hyperbolic_polygons.triangle(*angles, ring=ring)

    from flatsurf.geometry.hyperbolic import HyperbolicPlane
    from flatsurf.geometry.subfield import number_field_elements_from_algebraics

    # The input periods m1, m2, m3 are at the vertices 0, 2, 1.
    order = [
        vertex
        for _, vertex in sorted(
            zip((m1, m2, m3), (0, 2, 1)), key=lambda pv: (pv[0] % 2, -pv[0])
        )
    ]

    # No placement over QQ is sought as _optimize_placement does: when
    # k_0 is QQ the group is arithmetic and A_0 a division algebra
    # (Takeuchi, Theorem 1), so the bound is 2.
    bound = _triangle_group_degree_bound(m1, m2, m3)

    candidates = []
    for vertices in _triangle_placements(angles, order):
        field, coordinates = number_field_elements_from_algebraics(
            [c for vertex in vertices for c in vertex]
        )
        candidates.append(
            (field.degree(), len(candidates), field, coordinates, vertices)
        )
        if field.degree() <= bound:
            break
    _, _, field, coordinates, vertices = min(
        candidates, key=lambda candidate: candidate[:2]
    )

    if ring is None:
        # The minimal field is already computed; _polygon would compute it
        # again.
        return HyperbolicPlane(field).polygon(
            vertices=list(zip(coordinates[::2], coordinates[1::2]))
        )
    return _polygon(vertices, ring)


class FuchsianGroupGenerators:
    r"""
    Constructors of Fuchsian groups, as groups of matrices of `\mathrm{SO}(2,
    1)^+` acting on the Klein model of the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf import fuchsian_groups
        sage: fuchsian_groups
        fuchsian groups

    """

    def __repr__(self):
        r"""
        Return a printable representation of this constructor.

        EXAMPLES::

            sage: from flatsurf import fuchsian_groups
            sage: fuchsian_groups
            fuchsian groups

        """
        return "fuchsian groups"

    @staticmethod
    def triangle(m1, m2, m3, ring=None, optimize_placement=True):
        r"""
        Return the cocompact triangle group with periods ``m1``, ``m2`` and
        ``m3``.

        The group is generated by the rotations `g_1, g_2, g_3` by the
        angles `2\pi/m_1, 2\pi/m_2, 2\pi/m_3` about the vertices of a
        triangle with interior angles `\pi/m_i`, and has the presentation
        `g_1^{m_1} = g_2^{m_2} = g_3^{m_3} = g_1 g_2 g_3 = 1` of Katok,
        Fuchsian groups (1992), Section 4.3, p. 99. The vertices of periods
        `m_1, m_3, m_2` are the vertices 0, 1, 2 of the triangle,
        counterclockwise.

        INPUT:

        - ``m1``, ``m2``, ``m3`` -- integers at least 2, with `1/m_1 + 1/m_2
          + 1/m_3 < 1`

        - ``ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the matrices: with ``None``, the number field generated by the
          Klein coordinates of the triangle; ``AA`` keeps the coordinates in
          ``AA``, which makes the construction slow; any other ring receives
          the coordinates by conversion

        - ``optimize_placement`` -- a boolean (default: ``True``); whether to
          search for a placement of the triangle whose coordinates generate a
          field of small degree, see below. Otherwise the vertex of period
          ``m1`` is at the origin and its edge to the vertex of period
          ``m3`` on the positive `x`-axis.

        OUTPUT:

        A matrix group generated by the `3 \times 3` matrices of `g_1, g_2,
        g_3` in `\mathrm{SO}(2, 1)^+`. They need not lift to
        `\mathrm{PSL}_2` over the base ring.

        ALGORITHM:

        The generators are the images by the
        :meth:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces.FiniteType.Oriented.ParentMethods.holonomy_representation`
        of the loops around the vertices of the
        :meth:`~HyperbolicSurfaceGenerators.pillowcase` of the triangle.
        The loops turn counterclockwise, and their images multiply to the
        identity in the clockwise order of the vertices, `0, 2, 1`.

        The base ring is the field of the Klein coordinates of the
        triangle, which depends on its placement. With
        ``optimize_placement``, the vertices are tried by decreasing even
        period, then by decreasing odd period, in the order of the input on
        ties; at each vertex, the vertex is at the origin and the edge to the
        next vertex, the edge to the previous vertex, or the bisector of its
        angle on the positive `x`-axis. The search stops at the first
        placement whose field has a degree at most a lower bound, otherwise
        it keeps the first placement of smallest degree.

        The bound is `[k_0 : \QQ]` or `2 [k_0 : \QQ]`, with `k_0` the
        invariant trace field of Takeuchi, Arithmetic triangle groups
        (1977), Propositions 4 and 5, the second when the quaternion algebra
        `A_0` of the group is a division algebra over `k_0`: a real field
        `F` that contains the coordinates of a triangle puts the group in
        `\mathrm{PGL}_2(F)`, so it contains `k_0` and splits `A_0`.

        EXAMPLES:

        The triangle group `(3, 4, 4)`::

            sage: from flatsurf import fuchsian_groups
            sage: G = fuchsian_groups.triangle(3, 4, 4)
            sage: G
            Matrix group over Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178? with 3 generators (
            [  -7/2   -3/2  3/2*a]  [ 0 -1  0]  [  -3/2   -1/2  1/2*a]
            [   3/2   -1/2 -1/2*a]  [ 1  0  0]  [   7/2   -3/2 -3/2*a]
            [-3/2*a -1/2*a      4], [ 0  0  1], [-3/2*a  1/2*a      4]
            )
            sage: G.base_ring()
            Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?
            sage: g1, g2, g3 = G.gens()
            sage: (g1^3).is_one(), (g2^4).is_one(), (g3^4).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True

        The triangle group `(2, 4, 6)`, given by Python integers::

            sage: G = fuchsian_groups.triangle(int(2), int(4), int(6))
            sage: G.base_ring()
            Number Field in a with defining polynomial y^2 - 3 with a = 1.732050807568878?
            sage: g1, g2, g3 = G.gens()
            sage: (g1^2).is_one(), (g2^4).is_one(), (g3^6).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True

        The triangle group `(2, 3, 7)` is defined over a field of degree
        6, the lower bound, twice the degree of its invariant trace field::

            sage: # long time
            sage: G = fuchsian_groups.triangle(2, 3, 7)
            sage: G.base_ring().polynomial()
            y^6 + y^4 - 2*y^2 - 1
            sage: g1, g2, g3 = G.gens()
            sage: (g1^2).is_one(), (g2^3).is_one(), (g3^7).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True

        Its triangle has the vertex of period 2 at the origin and the edge
        to the vertex of period 7 on the positive `x`-axis; the invariant
        trace field of its pillowcase is `\QQ(\cos(2\pi/7))`::

            sage: # long time
            sage: from flatsurf import hyperbolic_polygons, hyperbolic_surfaces
            sage: T = hyperbolic_polygons.triangle(1/4, 1/14, 1/6, ring=G.base_ring())
            sage: hyperbolic_surfaces.pillowcase(T).holonomy_representation().invariant_trace_field()
            Number Field in ... with defining polynomial x^3 - 7*x - 7 with ...

        With the vertex of period `m_1` at the origin, the field may be
        larger::

            sage: G = fuchsian_groups.triangle(3, 4, 4, optimize_placement=False)
            sage: G.base_ring().degree()
            4

        The coordinates can be kept in ``AA``::

            sage: G = fuchsian_groups.triangle(3, 4, 4, ring=AA)
            sage: G.base_ring()
            Algebraic Real Field
            sage: g1, g2, g3 = G.gens()
            sage: (g1 * g2 * g3).is_one()
            True

        The periods are integers at least 2, and the group is hyperbolic::

            sage: fuchsian_groups.triangle(1, 4, 4)
            Traceback (most recent call last):
            ...
            ValueError: the period 1 is not an integer at least 2
            sage: fuchsian_groups.triangle(2, 3, 6)
            Traceback (most recent call last):
            ...
            ValueError: the triangle group (2, 3, 6) is euclidean, not hyperbolic: 1/m1 + 1/m2 + 1/m3 must be less than 1
            sage: fuchsian_groups.triangle(2, 3, 5)
            Traceback (most recent call last):
            ...
            ValueError: the triangle group (2, 3, 5) is spherical, not hyperbolic: 1/m1 + 1/m2 + 1/m3 must be less than 1

        TESTS::

            sage: fuchsian_groups.triangle(3, 4, 4, ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field
            sage: fuchsian_groups.triangle(5/2, 4, 4)
            Traceback (most recent call last):
            ...
            ValueError: the period 5/2 is not an integer at least 2

        """
        from sage.groups.matrix_gps.finitely_generated import MatrixGroup
        from flatsurf.geometry.signature import SurfaceSignature

        geometry = SurfaceSignature.from_orbifold_orders(0, (m1, m2, m3)).geometry()
        if geometry != "hyperbolic":
            raise ValueError(
                f"the triangle group ({m1}, {m2}, {m3}) is {geometry}, not hyperbolic: 1/m1 + 1/m2 + 1/m3 must be less than 1"
            )
        T = _triangle_group_polygon(m1, m2, m3, ring, optimize_placement)

        S = HyperbolicSurfaceGenerators.pillowcase(T)
        phi = S.holonomy_representation()
        loops = phi.domain().vertex_loops()
        g1, g3, g2 = (phi(loops[S(0, i)]).matrix() for i in range(3))
        return MatrixGroup([g1, g2, g3])


fuchsian_groups = FuchsianGroupGenerators()
