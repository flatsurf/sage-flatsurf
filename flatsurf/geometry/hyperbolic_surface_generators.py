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
            SurfaceSignature(genus=0, cone_orders=(), boundary_traces=(21/4, 21/4, 20881/144))
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
        from flatsurf.geometry.similarity_surface_generators import (
            _pillowcase_gluings,
        )

        if not isinstance(P, HyperbolicPolygon):
            raise TypeError("P must be a hyperbolic polygon")

        vertices = [v.coordinates(model="klein") for v in P.vertices()]
        m, permutation, gluings = _pillowcase_gluings(
            len(vertices), glued_edges, mirror_edge
        )

        H = P.parent()
        reflection = _reflection(vertices[m], vertices[(m + 1) % len(vertices)])

        mirrored = []
        for j in permutation:
            x, y, z = reflection * vector((vertices[j][0], vertices[j][1], 1))
            mirrored.append((x / z, y / z))
        Q = H.polygon(vertices=mirrored)

        surface = MutableOrientedHyperbolicSurface(H.base_ring())
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
            SurfaceSignature(genus=0, cone_orders=(), boundary_traces=(21/4, 21/4, 21/4))
            sage: S.area()
            1
            sage: TestSuite(S).run()

        With the vertex 0 of the hexagon at the origin, the coordinates need
        `\sqrt{6}`::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), optimize_placement=False)
            sage: S.base_ring()
            Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?
            sage: S.signature()
            SurfaceSignature(genus=0, cone_orders=(), boundary_traces=(21/4, 21/4, 21/4))
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
            SurfaceSignature(genus=0, cone_orders=(-1/2,), boundary_traces=(21/4, 21/4))
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
            SurfaceSignature(genus=0, cone_orders=(-2/3, -1/2), boundary_traces=(21/4,))
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
            SurfaceSignature(genus=0, cone_orders=(-7/8, -3/4, -1/2))
            sage: S.area()
            1/8

        A cone angle larger than `\pi` requires cuffs of comparable lengths::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4), cone_orders=(-1/4,))
            sage: S.signature()
            SurfaceSignature(genus=0, cone_orders=(-1/4,), boundary_traces=(21/4, 21/4))
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
            ValueError: the Euler characteristic 2 + sum(cone_orders) - len(boundary_traces) must be negative, not 0

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
            SurfaceSignature(genus=0, cone_orders=(), boundary_traces=(4, 4, 4))
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
            ValueError: the Euler characteristic 2 + sum(cone_orders) - len(boundary_traces) must be negative, not 0
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
                raise ValueError(
                    f"the boundary trace {t} must be a real algebraic number"
                )
            if not trace > 3:
                raise ValueError(f"the boundary trace {t} must be larger than 3")
            cosh.append((trace + 1).sqrt() / 2)

        euler_characteristic = 2 + sum(orders) - len(cosh)
        if not euler_characteristic < 0:
            raise ValueError(
                f"the Euler characteristic 2 + sum(cone_orders) - len(boundary_traces) must be negative, not {euler_characteristic}"
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
            SurfaceSignature(genus=2, cone_orders=())
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
            SurfaceSignature(genus=2, cone_orders=())
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


hyperbolic_surfaces = HyperbolicSurfaceGenerators()
