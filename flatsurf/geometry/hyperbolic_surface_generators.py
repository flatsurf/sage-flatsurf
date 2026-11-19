r"""
Generators of hyperbolic polygons and surfaces.

The polygons are built exactly from their interior angles, given as
multiples of `2\pi`, and from the hyperbolic cosine of some of their sides.
Each polygon has its vertex 0 at the origin of the Klein model and its edge
0 on the positive `x`-axis, with its vertices in counterclockwise order.
Every polygon except the regular one lies in the hyperbolic plane of a
quadratic form of its own.

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

from flatsurf.geometry.subfield import _check_exact


def _angle(angle, name="angle"):
    r"""
    Return ``angle`` as a rational number, the multiple of `2\pi` of an
    angle.

    An angle must coerce into `\QQ`: a float is turned into the simplest
    rational close to it, which is wrong after arithmetic and selects another
    number field.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _angle
        sage: _angle(1/4)
        1/4
        sage: _angle(sqrt(2))
        Traceback (most recent call last):
        ...
        ValueError: angle must be a rational multiple of 2π
        sage: _angle(0.125)
        Traceback (most recent call last):
        ...
        ValueError: angle must be a rational multiple of 2π

    """
    try:
        return QQ.coerce(angle)
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


def _cos(angle):
    r"""
    Return the cosine of `2\pi` ``angle`` in ``AA``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _cos
        sage: _cos(1/8)
        0.7071067811865475?
        sage: _cos(1/3)
        -0.500000000000000?

    """
    angle = QQ(angle)
    z = QQbar.zeta(angle.denominator()) ** angle.numerator()
    return AA(z.real())


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


def _walk(sides, angles):
    r"""
    Return the Klein coordinates in ``AA`` of the vertices of the polygon
    reached by a walk from the origin along the positive `x`-axis.

    INPUT:

    - ``sides`` -- a list of pairs ``(cosh, sinh)`` of the sides 0, 1, ...
      of the walk

    - ``angles`` -- a list of pairs ``(cos, sin)`` of the interior angles at
      the vertices 1, 2, ..., one fewer than ``sides``

    EXAMPLES:

    A walk along two sides of `\cosh = 5/4` with a right angle between
    them::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _walk
        sage: _walk([(AA(5/4), AA(3/4))] * 2, [(AA(0), AA(1))])
        [(0, 0), (3/5, 0), (3/5, 12/25)]

    ALGORITHM:

    A frame of the hyperboloid model starts at the origin, pointing along
    the positive `x`-axis. A side of hyperbolic cosine `c` and sine `s`
    multiplies the frame on the right by the translation `((c, 0, s), (0, 1,
    0), (s, 0, c))` along its first axis, and a vertex of interior angle
    `\alpha` by the rotation by `\pi - \alpha`, a left turn. Each vertex is
    the last column of its frame, divided by its last entry.

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


def _polygon(vertices, base_ring):
    r"""
    Return the hyperbolic polygon with the Klein coordinates ``vertices`` in
    ``AA``, over ``base_ring``.

    INPUT:

    - ``vertices`` -- a list of pairs of elements of ``AA``

    - ``base_ring`` -- a ring or ``None``; when ``None``, the field of
      :func:`~flatsurf.geometry.subfield.number_field_elements_from_algebraics`
      of the coordinates, ``QQ`` when they are rational; otherwise the
      coordinates are converted into ``base_ring``

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

    if base_ring is None:
        from flatsurf.geometry.subfield import number_field_elements_from_algebraics

        base_ring, coordinates = number_field_elements_from_algebraics(coordinates)
    elif base_ring is not AA:
        try:
            coordinates = [base_ring(c) for c in coordinates]
        except (TypeError, ValueError):
            raise ValueError(
                f"the coordinates of the polygon do not lie in {base_ring}"
            ) from None

    return HyperbolicPlane(base_ring).polygon(
        vertices=list(zip(coordinates[::2], coordinates[1::2]))
    )


def _vinberg_triangle(angles):
    r"""
    Return the field, the quadratic form and the Klein coordinates of the
    hyperbolic triangle with the interior angles ``angles`` in the basis of
    Vinberg.

    INPUT:

    - ``angles`` -- three positive rational numbers, the angles at the
      vertices 0, 1, 2 as multiples of `2\pi`, of sum less than `1/2`

    OUTPUT:

    A triple ``(K, (b, c), vertices)``: ``K`` is the field `k_0 = \QQ(u_0,
    u_1, u_2, p)` with `u_k = 4 \cos^2 \alpha_k` and `p = 8 \cos \alpha_0
    \cos \alpha_1 \cos \alpha_2`, ``QQ`` or a number field embedded in the
    reals; the triangle lies in the hyperbolic plane of the form `x^2 + b
    y^2 - c z^2`, with `b` and `c` positive elements of ``K``, squarefree
    integers when ``K`` is ``QQ``; ``vertices`` is the list of the Klein
    coordinates in ``K`` of the vertices 0, 1, 2,
    counterclockwise, the vertex 0 at the origin and the vertex 1 on the
    positive `x`-axis.

    EXAMPLES:

    The triangle of the triangle group `(3, 4, 4)` is over `\QQ`, with a
    form of squarefree integers::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _vinberg_triangle
        sage: K, (b, c), vertices = _vinberg_triangle([1/6, 1/8, 1/8])
        sage: K
        Rational Field
        sage: b, c
        (3, 6)
        sage: vertices
        [(0, 0), (2, 0), (1, 1)]

    The triangle of the triangle group `(2, 3, 7)` is over the cubic
    invariant trace field of the group::

        sage: K, (b, c), vertices = _vinberg_triangle([1/4, 1/14, 1/6])
        sage: K
        Number Field in a with defining polynomial x^3 - 7*x - 7 with a = 3.048917339522306?
        sage: b, c
        (a^2 - a - 3, 4*a + 4)
        sage: vertices
        [(0, 0), (2, 0), (0, -6*a^2 + 8*a + 32)]

    ALGORITHM:

    The triangle is the polygon of :func:`_seam_data` with three cone points,
    the vertex `k` of the triangle in the slot `k - 1 \bmod 3`, so that its
    edge 0 is on the mirror 1 of the slots. The root mirror is the one
    opposite the right angle if there is one, otherwise the side opposite
    the vertex 1 of the triangle, the mirror 0 of the slots.

    """
    from flatsurf.geometry.subfield import _triangle_invariant_trace_field

    K, u, p = _triangle_invariant_trace_field(angles)

    right = [k for k in range(3) if u[k] == 0]
    root = (right[0] if right else 1) - 1
    form, vertices = _seam_data(K, (u[1], u[2], u[0]), p, root % 3)
    return K, form, vertices


def _seam_field(cosines):
    r"""
    Return the field `\QQ(u_0, u_1, u_2, p)` of the generalized triangle
    whose slots have the signed cosines ``cosines``, with these elements.

    INPUT:

    - ``cosines`` -- three elements `c_0, c_1, c_2` of ``AA``; `c_k = \cos
      \theta_k` for a cone point of polygon angle `\theta_k`, `c_k = \cosh(\ell_k
      / 2)` for a cuff of length `\ell_k`

    OUTPUT:

    A triple ``(K, u, p)``: ``K`` is ``QQ`` or a number field embedded in the
    reals, ``u`` the tuple of the `u_k = 4 c_k^2` and ``p`` the element `8 c_0
    c_1 c_2` of ``K``.

    EXAMPLES:

    Three cuffs of trace `5`, of `\cosh = \sqrt{6}/2`::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _seam_field
        sage: c = AA(6).sqrt() / 2
        sage: _seam_field([c, c, c])
        (Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?,
         (6, 6, 6),
         6*a)
        sage: _seam_field([AA(5/4), AA(5/4), AA(5/4)])
        (Rational Field, (25/4, 25/4, 25/4), 125/8)

    """
    from flatsurf.geometry.subfield import number_field_elements_from_algebraics

    u = [4 * c * c for c in cosines]
    p = 8 * cosines[0] * cosines[1] * cosines[2]
    K, elements = number_field_elements_from_algebraics(u + [p])
    return K, tuple(elements[:3]), elements[3]


def _seam_data(K, u, p, root=None):
    r"""
    Return the quadratic form and the Klein coordinates of the front polygon
    of the generalized triangle with the data ``u``, ``p`` over ``K``.

    A generalized triangle has three slots `k = 0, 1, 2`, each a cone point or
    a cuff, and three mirrors, the mirror `k` opposite the slot `k`. The
    mirrors `i` and `j` meet at the cone point of the slot `k` with the angle
    `\theta_k`, or are disjoint, their common perpendicular being the half
    of the cuff of the slot `k`, of length `\ell_k/2`. The front polygon has
    the edges, counterclockwise, mirror 1, slot 0, mirror 2, slot 1, mirror
    0, slot 2, a cone point being a vertex instead of an edge.

    INPUT:

    - ``K`` -- ``QQ`` or a number field embedded in the reals

    - ``u`` -- three elements `u_k = 4 c_k^2` of ``K``, with `c_k = \cos
      \theta_k` for a cone point and `c_k = \cosh(\ell_k / 2)` for a cuff, so
      that the slot `k` is a cuff when `u_k > 4`

    - ``p`` -- the element `8 c_0 c_1 c_2` of ``K``

    - ``root`` -- ``0``, ``1``, ``2`` or ``None`` (default: ``None``); the root
      mirror, by default the mirror opposite a right angle if there is one,
      the mirror 1 otherwise

    OUTPUT:

    A pair ``((b, c), vertices)``: the polygon lies in the hyperbolic plane of
    the form `x^2 + b y^2 - c z^2`, with `b` and `c` positive elements of
    ``K``, squarefree integers when ``K`` is ``QQ``; ``vertices`` is the list
    of the Klein coordinates in ``K`` of its vertices, counterclockwise, the
    vertex 0 at the origin and the vertex 1 on the positive `x`-axis.

    EXAMPLES:

    The right-angled hexagon of the pair of pants with three cuffs of trace
    `21/4`, `c_k = 5/4`::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _seam_data
        sage: _seam_data(QQ, (25/4, 25/4, 25/4), 125/8)
        ((6, 6),
         [(0, 0), (12/5, 0), (12/5, 3/25), (60/29, 15/29), (48/25, 3/5), (0, 3/5)])

    A cone point of angle `\pi/2` in the slot 0 gives a pentagon::

        sage: _seam_data(QQ, (0, 25/4, 25/4), 0)
        ((34, 34), [(0, 0), (5, 0), (5, 15/34), (544/125, 3/5), (0, 3/5)])

    ALGORITHM:

    The unit normals `e_k` of the mirrors have the Gram matrix `B_{kk} = 1`,
    `B_{ij} = -c_k` for `\{i, j, k\} = \{0, 1, 2\}`, of signature `(2, 1)`.
    Vinberg's rescaling of the normals takes the root mirror `r`, such that
    `B(e_r, e_i)` is nonzero for both other mirrors, and sets `f_r = e_r`, `f_i
    = 2 B(e_r, e_i) e_i`. The Gram matrix `G` of the `f` has `G_{rr} = 1`,
    `G_{ii} = u_k` and `G_{ri} = u_k / 2`, with `k` the slot between the
    mirrors `r` and `i`, and `G_{ij} = -p/2` for `i, j \neq r`: its entries are
    in ``K``. In the basis `f`, the cone point of the slot `k`, orthogonal to
    the two mirrors through it, is the column `v_k` of `G^{-1}`; for a cuff,
    `v_k` is spacelike, the pole of the line of the cuff, and the corner of
    the mirror `a` and the cuff `k`, orthogonal to `f_a` and to `v_k`, is the
    cross product of the row `a` of `G` and of the `k`-th vector of the
    canonical basis; its Klein coordinates do not depend on its sign.

    The basis `b_1, b_2 = f_1, o` with `o` the vertex 0, and `b_1` the pole
    `v_2` of the cuff of the slot 2, or `f_0 - G(f_1, f_0) / G(f_1, f_1) f_1`
    when the slot 2 is a cone point, is orthogonal for `G`, so that the form
    is `\mathrm{diag}(1, b, -c)` in it once divided by `G(b_1, b_1)`. The
    vertex 1, on the mirror 1, is on the `x`-axis, and the signs of `b_1` and
    `b_2` put it on its positive half and the vertex 2 above it.

    Over `\QQ`, with `\lambda^2 = c / c'` and `\mu^2 = \lambda^2 b' / b`, `b'`
    and `c'` the squarefree parts of `b` and `c`, the Klein coordinates `(x /
    \lambda, y / \mu)` are those of the polygon in the plane of the form `x^2
    + b' y^2 - c' z^2`.

    """
    if root is None:
        # f_i vanishes when the slot between the mirrors r and i is a right
        # angle, so the root mirror is opposite the right angle.
        right = [k for k in range(3) if u[k] == 0]
        root = right[0] if right else 1

    G = matrix(K, 3, 3)
    for i in range(3):
        for j in range(3):
            if i == j:
                G[i, i] = 1 if i == root else u[3 - root - i]
            elif root in (i, j):
                # The mirrors i and j meet at the slot 3 - i - j.
                G[i, j] = u[3 - i - j] / 2
            else:
                G[i, j] = -p / 2

    poles = G.inverse().columns()
    cuff = [u[k] > 4 for k in range(3)]

    corners = []
    for a, k in [(1, 2), (1, 0), (2, 0), (2, 1), (0, 1), (0, 2)]:
        if cuff[k]:
            v = G.row(a).cross_product(vector(K, [int(i == k) for i in range(3)]))
        elif (a, k) in [(0, 2), (2, 0), (0, 1)]:
            # the second corner of a cone point, the same vertex as the first
            continue
        else:
            v = poles[k]
        corners.append(v)

    f0 = vector(K, (1, 0, 0))
    f1 = vector(K, (0, 1, 0))
    if cuff[2]:
        b1 = poles[2]
    else:
        b1 = f0 - (f1 * G * f0) / (f1 * G * f1) * f1
    W = matrix(K, [b1, f1, corners[0]]).transpose()
    norm = b1 * G * b1
    b = (f1 * G * f1) / norm
    c = -(corners[0] * G * corners[0]) / norm

    W_inverse = W.inverse()
    coordinates = []
    for v in corners:
        x, y, z = W_inverse * v
        coordinates.append((x / z, y / z))

    sx = 1 if coordinates[1][0] > 0 else -1
    sy = 1 if coordinates[2][1] > 0 else -1
    if K is QQ:
        b_squarefree, c_squarefree = b.squarefree_part(), c.squarefree_part()
        lambda_ = (c / c_squarefree).sqrt()
        mu = (lambda_**2 * b_squarefree / b).sqrt()
        sx, sy = sx / lambda_, sy / mu
        b, c = QQ(b_squarefree), QQ(c_squarefree)
    return (b, c), [(sx * x, sy * y) for x, y in coordinates]


def _seam_polygon(cosines, base_ring):
    r"""
    Return the front polygon of the generalized triangle whose slots have the
    signed cosines ``cosines``, in the plane of its form.

    INPUT:

    - ``cosines`` -- as in :func:`_seam_field`

    - ``base_ring`` -- as in :func:`_plane_polygon`

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _seam_polygon
        sage: _seam_polygon([AA(5/4), AA(5/4), AA(5/4)], None).parent()
        Hyperbolic Plane over Rational Field with form x^2 + 6*y^2 - 6*z^2

    """
    K, u, p = _seam_field(cosines)
    form, vertices = _seam_data(K, u, p)
    return _plane_polygon(K, form, vertices, base_ring)


def _plane_polygon(K, form, vertices, base_ring):
    r"""
    Return the polygon with the Klein coordinates ``vertices`` in the
    hyperbolic plane of the form `x^2 + b y^2 - c z^2`, over ``base_ring``.

    INPUT:

    - ``K`` -- ``QQ`` or a number field embedded in the reals

    - ``form`` -- a pair ``(b, c)`` of positive elements of ``K``

    - ``vertices`` -- a list of pairs of elements of ``K``

    - ``base_ring`` -- a field or ``None``; with ``None``, ``K``, otherwise an
      exact field into which ``K`` converts along the embeddings in the
      reals, see :func:`_field_embedding`

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _plane_polygon
        sage: P = _plane_polygon(QQ, (3, 6), [(0, 0), (2, 0), (1, 1)], None)
        sage: P.parent()
        Hyperbolic Plane over Rational Field with form x^2 + 3*y^2 - 6*z^2
        sage: _plane_polygon(QQ, (3, 6), [(0, 0), (2, 0), (1, 1)], AA).parent()
        Hyperbolic Plane over Algebraic Real Field with form x^2 + 3*y^2 - 6*z^2

    """
    from flatsurf.geometry.hyperbolic import HyperbolicPlane

    if base_ring is None:
        base_ring = K
    convert = _field_embedding(K, base_ring)

    H = HyperbolicPlane(base_ring, form=[1, convert(form[0]), -convert(form[1])])
    return H.polygon(vertices=[(convert(x), convert(y)) for x, y in vertices])


def _field_embedding(K, ring):
    r"""
    Return a map from the field ``K`` into ``ring`` that is compatible with
    their embeddings in the reals.

    INPUT:

    - ``K`` -- ``QQ`` or a number field embedded in the real algebraic
      numbers

    - ``ring`` -- an exact ring

    OUTPUT:

    ``ring`` itself, as a conversion, when ``K`` coerces into it, otherwise
    a field morphism into a number field ``ring`` embedded in the real
    algebraic numbers.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _field_embedding
        sage: K.<a> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: _field_embedding(QQ, K)
        Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?
        sage: _field_embedding(K, AA)(a)
        1.414213562373095?
        sage: L.<s> = NumberField(x^4 - 2*x^2 - 1, embedding=AA.polynomial_root(AA['x']([-1, 0, -2, 0, 1]), RIF(1.5, 1.6)))
        sage: _field_embedding(K, L)(a)
        s^2 - 1
        sage: _field_embedding(K, QQ)
        Traceback (most recent call last):
        ...
        ValueError: the coordinates of the polygon do not lie in Rational Field

    """
    if ring.has_coerce_map_from(K):
        return ring

    from sage.rings.number_field.number_field_base import NumberField

    if (
        K is not QQ
        and isinstance(ring, NumberField)
        and ring.coerce_embedding() is not None
        and ring.coerce_embedding().codomain() is K.coerce_embedding().codomain()
    ):
        from flatsurf.geometry.subfield import is_embedded_subfield

        embedded, morphism = is_embedded_subfield(K, ring, certificate=True)
        if embedded:
            return morphism

    raise ValueError(f"the coordinates of the polygon do not lie in {ring}")


class HyperbolicPolygonsConstructor:
    r"""
    Constructors of compact hyperbolic polygons from their angles and the
    hyperbolic cosine of some of their sides.

    Angles are multiples of `2\pi`, so that `1/4` is a right angle. Each
    polygon has its vertex 0 at the origin of the Klein model and its edge 0
    on the positive `x`-axis. Every polygon but the regular one lies in the
    hyperbolic plane of a quadratic form `x^2 + b y^2 - c z^2` of its own,
    over the field `\QQ(u_0, u_1, u_2, p)` of :meth:`triangle` and
    :meth:`right_angled_hexagon`, with `b` and `c` squarefree integers over
    `\QQ`; ``base_ring`` is then an exact field into which that field
    converts along the embeddings in the reals, such as ``AA``, which
    receives the form and the coordinates. The regular polygon lies in the
    plane of the form `x^2 + y^2 - z^2`: with ``base_ring=None``, over the
    number field generated by the Klein coordinates of its vertices, ``QQ``
    when they are rational; ``base_ring=AA`` keeps the coordinates in
    ``AA``; any other ring receives the coordinates by conversion.

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
    def triangle(a, b, c, base_ring=None):
        r"""
        Return the hyperbolic triangle with interior angles ``a``, ``b`` and
        ``c`` at its vertices 0, 1 and 2.

        INPUT:

        - ``a``, ``b``, ``c`` -- positive rational numbers, the angles as
          multiples of `2\pi`, of sum less than `1/2`

        - ``base_ring`` -- a field or ``None`` (default: ``None``); with
          ``None``, the field `k_0` below; otherwise an exact field into
          which `k_0` converts along the embeddings in the reals, such as
          ``AA``, which receives the form and the coordinates

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin and edge 0 on the positive `x`-axis, in
        the hyperbolic plane of the form `x^2 + b y^2 - c z^2` over the
        field `k_0 = \QQ(u_0, u_1, u_2, p)`, with `u_k = 4 \cos^2 \alpha_k`
        for the angles `\alpha_0, \alpha_1, \alpha_2` and `p = 8 \cos \alpha_0
        \cos \alpha_1 \cos \alpha_2`. For the angles `\pi/m_k`, `k_0` is the
        invariant trace field of the triangle group of periods `(m_0, m_1,
        m_2)`.

        The plane of the form `x^2 + y^2 - z^2` may need a field of twice
        the degree: :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.polygon`
        with ``angles`` builds the triangle there over a field of choice,
        and :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.standard_form_isometry`
        of the plane of the triangle carries it there.

        EXAMPLES:

        A triangle with a right angle at the origin::

            sage: from flatsurf import hyperbolic_polygons
            sage: T = hyperbolic_polygons.triangle(1/4, 1/16, 1/16)
            sage: T.parent()
            Hyperbolic Plane over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? with form ...
            sage: T.angles()
            (1/4, 1/16, 1/16)
            sage: T.area()
            1/8

        The triangle of the triangle group `(2, 3, 7)`, over the cubic
        invariant trace field of the group::

            sage: T = hyperbolic_polygons.triangle(1/4, 1/6, 1/14)
            sage: T.base_ring().polynomial()
            x^3 - 7*x - 7
            sage: T.area()
            1/84

        The isomorphism to the plane of `x^2 + y^2 - z^2` adjoins a square
        root::

            sage: f = T.parent().standard_form_isometry(algorithm="extension")
            sage: f(T).parent()
            Hyperbolic Plane over Number Field in s with defining polynomial ... with s = ...
            sage: f(T).base_ring().degree()
            6
            sage: f(T).angles()
            (1/4, 1/6, 1/14)

        The coordinates can be kept in ``AA``::

            sage: T = hyperbolic_polygons.triangle(QQ((1, 4)), QQ((1, 16)), QQ((1, 16)), base_ring=AA)
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
            sage: hyperbolic_polygons.triangle(1/4, 1/16, 1/16, base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field

        A float angle is refused::

            sage: hyperbolic_polygons.triangle(1/6, 0.125, 0.125)
            Traceback (most recent call last):
            ...
            ValueError: angle must be a rational multiple of 2π

        The ring must be exact::

            sage: hyperbolic_polygons.triangle(1/8, 1/8, 1/8, base_ring=RR)
            Traceback (most recent call last):
            ...
            ValueError: the ring must be exact, not Real Field with 53 bits of precision

        A field that contains `k_0` is accepted, with or without a coercion
        from `k_0`::

            sage: K.<a> = NumberField(x^2 - 6, embedding=AA(6).sqrt())
            sage: T = hyperbolic_polygons.triangle(1/4, 1/12, 1/12, base_ring=K)
            sage: T.vertex(1).coordinates(model="klein")
            (2, 0)
            sage: T.parent()
            Hyperbolic Plane over Number Field in a with defining polynomial x^2 - 6 with a = 2.449489742783178? with form x^2 + y^2 - 6*z^2
            sage: L.<s> = NumberField(x^4 - 2*x^2 - 1, embedding=AA.polynomial_root(AA['x']([-1, 0, -2, 0, 1]), RIF(1.5, 1.6)))
            sage: T = hyperbolic_polygons.triangle(1/4, 1/16, 1/16, base_ring=L)
            sage: T.base_ring() is L
            True
            sage: T.angles()
            (1/4, 1/16, 1/16)

        A base ring embedded in the real lazy field::

            sage: x = polygen(QQ, 'x')
            sage: L = NumberField(x^4 - 2, 'r', embedding=1.18)
            sage: L.coerce_embedding().codomain()
            Real Lazy Field
            sage: hyperbolic_polygons.triangle(1/6, 1/6, 1/8, base_ring=L).parent()
            Hyperbolic Plane over Number Field in r with defining polynomial x^4 - 2 with r = 1.189207115002722? with form x^2 + 4/3*y^2 - r^2*z^2

        ALGORITHM:

        The unit normals `e_i` of the sides of the triangle, the side `i`
        opposite the vertex `i`, have the Gram matrix `B` with the entries
        `1` on the diagonal and `-\cos \alpha_k` for the two sides that meet
        at the vertex `k`. By Davis, The geometry and topology of Coxeter
        groups (2008), (6.15), p. 97, and Proposition 6.8.4, p. 98, a
        triangle of the hyperbolic plane with the angles `\alpha_k` exists
        if and only if this matrix has the signature `(2, 1)` and its
        principal submatrices are positive definite, which holds for
        positive angles of sum less than `\pi`. Vinberg's rescaling of the
        normals, `f_r = e_r` and `f_i = 2 B(e_r, e_i) e_i` for a root side
        `r`, the side opposite the right angle if there is one, the side 1
        otherwise, turns this Gram matrix into one over `k_0`, which is then
        diagonalised at the vertex 0.

        """
        _check_exact(base_ring)
        a, b, c = _angle(a), _angle(b), _angle(c)
        if not (a > 0 and b > 0 and c > 0):
            raise ValueError("the angles of a hyperbolic triangle must be positive")
        if not a + b + c < QQ((1, 2)):
            raise ValueError(
                "the angles of a hyperbolic triangle must sum to less than 1/2"
            )

        K, form, vertices = _vinberg_triangle([a, b, c])
        return _plane_polygon(K, form, vertices, base_ring)

    @staticmethod
    def regular_ngon(n, angle, base_ring=None):
        r"""
        Return the regular hyperbolic polygon with ``n`` vertices and
        interior angle ``angle``.

        INPUT:

        - ``n`` -- an integer, at least 3

        - ``angle`` -- a rational number, the interior angle as a multiple of
          `2\pi`, in `(0, 1/2 - 1/n)`

        - ``base_ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the polygon, see :class:`HyperbolicPolygonsConstructor`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin and edge 0 on the positive `x`-axis.

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

            sage: P = hyperbolic_polygons.regular_ngon(int(5), QQ((1, 5)), base_ring=AA)
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
            sage: hyperbolic_polygons.regular_ngon(QQ((11, 2)), 1/8)
            Traceback (most recent call last):
            ...
            ValueError: the number of vertices must be an integer
            sage: hyperbolic_polygons.regular_ngon(5.5, 1/8)
            Traceback (most recent call last):
            ...
            ValueError: the number of vertices must be an integer
            sage: hyperbolic_polygons.regular_ngon(8, 1/16, base_ring=RR)
            Traceback (most recent call last):
            ...
            ValueError: the ring must be exact, not Real Field with 53 bits of precision

        ALGORITHM:

        The triangles from the centre to the sides have angles `2\pi/n` at
        the centre and `\theta/2` at the vertices, with `\theta = 2\pi`
        ``angle``. The second cosine rule of Katok, Fuchsian groups (1992),
        Theorem 1.5.2, p. 16, gives the side `\cosh s = (\cos^2(\theta/2) +
        \cos(2\pi/n)) / \sin^2(\theta/2)`; the octagon of angle `1/8` is
        the one of Katok, Section 4.3, Example C, p. 95.

        """
        _check_exact(base_ring)
        try:
            n = ZZ(n)
        except (TypeError, ValueError):
            raise ValueError("the number of vertices must be an integer") from None
        if n < 3:
            raise ValueError("a polygon has at least 3 vertices")
        angle = _angle(angle)
        bound = QQ((1, 2)) - QQ((1, n))
        if not 0 < angle < bound:
            raise ValueError(
                f"the angle of a regular hyperbolic polygon with {n} vertices must be in (0, {bound})"
            )

        cos_centre = _cos(QQ((1, n)))
        cos_half, sin_half = _cos_sin(angle / 2)
        side = _from_cosh((cos_half**2 + cos_centre) / sin_half**2)

        turn = _cos_sin(angle)
        return _polygon(_walk([side] * (n - 1), [turn] * (n - 2)), base_ring)

    @staticmethod
    def right_angled_hexagon(c0, c1, c2, base_ring=None):
        r"""
        Return the right-angled hyperbolic hexagon whose edges 1, 3 and 5
        have hyperbolic cosine ``c0``, ``c1`` and ``c2``.

        The edges 0, 2 and 4 follow: edge 0 is opposite edge 3, edge 2 is
        opposite edge 5, and edge 4 is opposite edge 1.

        INPUT:

        - ``c0``, ``c1``, ``c2`` -- real algebraic numbers larger than 1

        - ``base_ring`` -- a field or ``None`` (default: ``None``); with
          ``None``, the field `K` below, otherwise an exact field into which
          `K` converts along the embeddings in the reals, such as ``AA``,
          which receives the form and the coordinates

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin, edge 0 on the positive `x`-axis, and
        edge 5 ending at vertex 0, in the hyperbolic plane of the form `x^2
        + b y^2 - c z^2` over the field `K = \QQ(u_0, u_1, u_2, p)`, with
        `u_k = 4 c_k^2` and `p = 8 c_0 c_1 c_2`; over `\QQ`, `b` and `c` are
        squarefree integers.

        EXAMPLES:

        A hexagon over `\QQ`::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 145/24)
            sage: P.parent()
            Hyperbolic Plane over Rational Field
            sage: [v.coordinates(model="klein") for v in P.vertices()]
            [(0, 0), (56/65, 0), (56/65, 99/325), (840/4121, 4015/4121), (1344/9425, 143/145), (0, 143/145)]
            sage: P.angles()
            (1/4, 1/4, 1/4, 1/4, 1/4, 1/4)
            sage: P.area()
            1/2

        The three alternate sides equal, over `\QQ` where `u_k = 25/4` and `p
        = 125/8`::

            sage: P = hyperbolic_polygons.right_angled_hexagon(5/4, 5/4, 5/4)
            sage: P.parent()
            Hyperbolic Plane over Rational Field with form x^2 + 6*y^2 - 6*z^2
            sage: P.vertices()
            ((0, 0), (12/5, 0), (12/5, 3/25), (60/29, 15/29), (48/25, 3/5), (0, 3/5))
            sage: P._cosh2_sides()
            (25, 25/16, 25, 25/16, 25, 25/16)

        The plane of the form `x^2 + y^2 - z^2` may need a larger field,
        see :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.standard_form_isometry`::

            sage: f = P.parent().standard_form_isometry()
            sage: f(P).base_ring()
            Rational Field

        Python integers are accepted::

            sage: P = hyperbolic_polygons.right_angled_hexagon(int(2), int(2), int(2))
            sage: P.angles()
            (1/4, 1/4, 1/4, 1/4, 1/4, 1/4)

        Sides must have a hyperbolic cosine larger than 1::

            sage: hyperbolic_polygons.right_angled_hexagon(1, 5/4, 5/4)
            Traceback (most recent call last):
            ...
            ValueError: c0 must be larger than 1

        ALGORITHM:

        Such a hexagon exists, and is unique, for any three lengths of its
        alternate sides, Buser, Geometry and spectra of compact Riemann
        surfaces (1992), Theorem 2.4.2. It is the front polygon of the pair
        of pants whose cuffs have the half-lengths of the edges 1, 3, 5, the
        slots 0, 1, 2: the lines of the edges 4, 0, 2, the mirrors 0, 1, 2,
        have unit normals `e_k` with the Gram matrix `B_{kk} = 1`, `B_{ij} =
        -c_k`, the mirrors `i` and `j` being disjoint with the half-cuff `k`
        as common perpendicular. Vinberg's rescaling of these normals gives a
        Gram matrix over `K`, in which the vertices and an orthogonal basis
        at the vertex 0 are computed, see :meth:`triangle`.

        """
        _check_exact(base_ring)
        h0 = _from_cosh(c0, "c0")
        h1 = _from_cosh(c1, "c1")
        h2 = _from_cosh(c2, "c2")

        return _seam_polygon([h0[0], h1[0], h2[0]], base_ring)

    @staticmethod
    def right_angled_pentagon(angle, c0, c1, base_ring=None):
        r"""
        Return the hyperbolic pentagon with interior angle ``angle`` at
        vertex 1, right angles at the other vertices, and edges 2 and 4 of
        hyperbolic cosine ``c0`` and ``c1``.

        The edges 0 and 1 meet at vertex 1; edge 3 is opposite it.

        INPUT:

        - ``angle`` -- a rational number in `(0, 1/2)`, the angle as a
          multiple of `2\pi`

        - ``c0``, ``c1`` -- real algebraic numbers larger than 1

        - ``base_ring`` -- a field or ``None`` (default: ``None``), as in
          :meth:`right_angled_hexagon`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin, edge 0 on the positive `x`-axis, and
        edge 4 ending at vertex 0, in the hyperbolic plane of the form `x^2
        + b y^2 - c z^2` over the field `K = \QQ(u_0, u_1, u_2, p)` of
        :meth:`right_angled_hexagon`, with `c_0` the cosine of the angle
        ``angle`` and `c_1, c_2` the arguments ``c0``, ``c1``.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_pentagon(3/8, 5/4, 5/4)
            sage: P.base_ring().polynomial()
            y^2 - 2
            sage: P.angles()
            (1/4, 3/8, 1/4, 1/4, 1/4)
            sage: P.area()
            1/8

        The Klein coordinates of the vertices of a pentagon over `\QQ`::

            sage: P = hyperbolic_polygons.right_angled_pentagon(1/6, 5/4, 5/4)
            sage: P.parent()
            Hyperbolic Plane over Rational Field with form x^2 + 7*y^2 - 7*z^2
            sage: [v.coordinates(model="klein") for v in P.vertices()]
            [(0, 0), (5/2, 0), (70/31, 15/31), (112/55, 3/5), (0, 3/5)]

        A pentagon over a cubic field::

            sage: P = hyperbolic_polygons.right_angled_pentagon(1/14, 5/4, 5/4)
            sage: P.base_ring().degree()
            3
            sage: P.area()
            3/7

        Python integers are accepted::

            sage: P = hyperbolic_polygons.right_angled_pentagon(QQ((1, 6)), int(2), int(2), base_ring=AA)
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

        ALGORITHM:

        The pentagon with the angle `\varphi = 2\pi` ``angle`` at the vertex
        1 and right angles elsewhere is the generalized triangle of Buser,
        Geometry and spectra of compact Riemann surfaces (1992), Example
        2.2.7 (second figure of Figure 2.2.4), p. 37. Its relation (v), with
        the edges 2 and 4 for the sides `\alpha, \beta` and the edge 3 for
        `c`, is `\cos\varphi = \sinh s_2 \sinh s_4 \cosh s_3 - \cosh s_2
        \cosh s_4`, so that `\cosh s_3 = (\cos\varphi + \cosh s_2 \cosh s_4)
        / (\sinh s_2 \sinh s_4)`. The book does not state the edges 1 and 0
        in this form; they are derived here from the same relations of the
        example: `\sinh s_1 = (\cosh s_4 + \cosh s_2 \cos\varphi) / (\sinh
        s_2 \sin\varphi)` and symmetrically. The pentagon exists if and only
        if these two hyperbolic sines are positive. For `\varphi = \pi/2` the
        pentagon is the right-angled one, with `\cosh s_3 = \sinh s_0 \sinh
        s_1`, Theorem 2.3.4 (i), p. 39.

        The pentagon is then built as the hexagon of
        :meth:`right_angled_hexagon` whose slot 0 is a cone point: the
        mirrors 1 and 2, the lines of the edges 0 and 1, meet at the vertex 1
        with the angle `\varphi`, and `B_{12} = -\cos\varphi`.

        """
        _check_exact(base_ring)
        angle = _angle(angle)
        if not 0 < angle < QQ((1, 2)):
            raise ValueError("the angle must be in (0, 1/2)")
        h0 = _from_cosh(c0, "c0")
        h1 = _from_cosh(c1, "c1")
        cos = _cos(angle)

        if not h1[0] + h0[0] * cos > 0:
            raise ValueError(
                "no right-angled pentagon with these parameters: c1 + c0 cos(2π angle) must be positive"
            )
        if not h0[0] + h1[0] * cos > 0:
            raise ValueError(
                "no right-angled pentagon with these parameters: c0 + c1 cos(2π angle) must be positive"
            )

        return _seam_polygon([cos, h0[0], h1[0]], base_ring)

    @staticmethod
    def right_angled_quadrilateral(angle0, angle1, c, base_ring=None):
        r"""
        Return the hyperbolic quadrilateral with interior angles ``angle0``
        and ``angle1`` at vertices 1 and 2, right angles at vertices 3 and
        0, and edge 3 of hyperbolic cosine ``c``.

        Edge 1, between the vertices 1 and 2, is opposite edge 3.

        INPUT:

        - ``angle0``, ``angle1`` -- positive rational numbers, the angles as
          multiples of `2\pi`, of sum less than `1/2`

        - ``c`` -- a real algebraic number larger than 1

        - ``base_ring`` -- a field or ``None`` (default: ``None``), as in
          :meth:`right_angled_hexagon`

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`
        with vertex 0 at the origin, edge 0 on the positive `x`-axis, and
        edge 3 ending at vertex 0, in the hyperbolic plane of the form `x^2
        + b y^2 - c z^2` over the field `K = \QQ(u_0, u_1, u_2, p)` of
        :meth:`right_angled_hexagon`, with `c_0` and `c_1` the cosines of
        the angles ``angle0`` and ``angle1`` and `c_2` the argument ``c``.

        EXAMPLES::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.right_angled_quadrilateral(1/4, 1/8, 5/4)
            sage: P.parent()
            Hyperbolic Plane over Rational Field with form x^2 + 34*y^2 - 34*z^2
            sage: [v.coordinates(model="klein") for v in P.vertices()]
            [(0, 0), (4, 0), (4, 3/5), (0, 3/5)]
            sage: P.angles()
            (1/4, 1/4, 1/8, 1/4)
            sage: P.area()
            1/8

        Python integers are accepted::

            sage: P = hyperbolic_polygons.right_angled_quadrilateral(QQ((1, 6)), QQ((1, 6)), int(2), base_ring=AA)
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

        ALGORITHM:

        Buser, Geometry and spectra of compact Riemann surfaces (1992),
        states the relations of a quadrilateral with three right angles,
        Theorem 2.3.1, and of one with two opposite angles that are not
        right, Example 2.2.7 (i) and (ii), p. 37; the two angles
        `\varphi_0, \varphi_1` here are adjacent. The book does not state the
        relations below; they are derived here as the analogue of Theorem
        2.4.1 (iii), p. 40, for the hexagon with two sides replaced by
        angles, from Theorem 2.2.6, p. 36, the identity on the matrices of
        the sides and angles from which Buser obtains all of them. The edge
        0, at the vertex of `\varphi_0`, has `\sinh h_0 = (\cos\varphi_1 +
        \cosh a \cos\varphi_0) / (\sinh a \sin\varphi_0)`, and symmetrically
        for the edge 2. The quadrilateral exists if and only if these two
        hyperbolic sines are positive. For `\varphi_1 = \pi/2` this is
        Theorem 2.3.1 (vi), `\sinh h_0 = \coth a \cot\varphi_0`, and Theorem
        2.3.1 (i), `\cos\varphi_0 = \sinh a \sinh h_2`, for the quadrilateral
        with three right angles.

        The quadrilateral is then built as the hexagon of
        :meth:`right_angled_hexagon` whose slots 0 and 1 are cone points, at
        the vertices 1 and 2.

        """
        _check_exact(base_ring)
        angle0 = _angle(angle0, "angle0")
        angle1 = _angle(angle1, "angle1")
        if not (angle0 > 0 and angle1 > 0 and angle0 + angle1 < QQ((1, 2))):
            raise ValueError("the two angles must be positive and sum to less than 1/2")
        h = _from_cosh(c, "c")
        cos0 = _cos(angle0)
        cos1 = _cos(angle1)

        if not cos1 + h[0] * cos0 > 0:
            raise ValueError(
                "no right-angled quadrilateral with these parameters: cos(2π angle1) + c cos(2π angle0) must be positive"
            )
        if not cos0 + h[0] * cos1 > 0:
            raise ValueError(
                "no right-angled quadrilateral with these parameters: cos(2π angle0) + c cos(2π angle1) must be positive"
            )

        return _seam_polygon([cos0, cos1, h[0]], base_ring)


hyperbolic_polygons = HyperbolicPolygonsConstructor()


def _pants_polygon(boundary_traces, cone_orders):
    r"""
    Return the field, the quadratic form and the Klein coordinates of the
    front polygon of the pair of pants of
    :meth:`HyperbolicSurfaceGenerators.pants`, and its seams.

    INPUT:

    - ``boundary_traces``, ``cone_orders`` -- as in
      :meth:`HyperbolicSurfaceGenerators.pants`

    OUTPUT:

    A tuple ``(K, (b, c), vertices, seams)``: the front polygon lies in the
    hyperbolic plane of the form `x^2 + b y^2 - c z^2` over the field ``K``
    of its slots, ``vertices`` is the list of the Klein coordinates in ``K``
    of its vertices, and ``seams`` the list of its edges that are seams. The
    other edges are the half-cuffs, in the order of ``boundary_traces``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _pants_polygon
        sage: _pants_polygon((21/4, 21/4), (-1/2,))
        (Rational Field,
         (34, 34),
         [(0, 0), (5, 0), (5, 15/34), (544/125, 3/5), (0, 3/5)],
         [0, 1, 3])
        sage: _pants_polygon((), (-1/2, -3/4, -7/8))
        (Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?,
         (8*a + 12, a),
         [(0, 0), (1, 0), (1, -a + 3/2)],
         [0, 1, 2])

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
            raise ValueError(f"the cone order {k} must be a rational number") from None
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
            ) from None
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
        K, form, vertices = _vinberg_triangle([angles[2], angles[0], angles[1]])
        return K, form, vertices, [0, 1, 2]

    try:
        P = constructor(*arguments)
    except ValueError as error:
        for inequality, message in inequalities.items():
            if inequality in str(error):
                raise ValueError(
                    f"no pair of pants with these parameters: {message}"
                ) from error
        raise

    form = P.parent().quadratic_form()
    vertices = [v.coordinates(model="klein") for v in P.vertices()]
    return P.base_ring(), (form[1, 1], -form[2, 2]), vertices, seams


def _mirror(P, glued_edges, mirror_edge):
    r"""
    Return the mirror image of the hyperbolic polygon ``P`` and the gluings
    of :meth:`HyperbolicSurfaceGenerators.pillowcase`.

    INPUT:

    - ``P``, ``glued_edges``, ``mirror_edge`` -- as in
      :meth:`HyperbolicSurfaceGenerators.pillowcase`

    OUTPUT:

    A pair ``(Q, gluings)``: the polygon `Q`, over the base ring of ``P``,
    and the list of the pairs of edges ``(label, edge)`` to glue, ``P``
    having the label 0 and `Q` the label 1. The polygon `Q` is the image of
    ``P`` under the reflection in the line of the mirror edge `m` (see
    :meth:`~flatsurf.geometry.categories.polygons.Polygons.ParentMethods.apply_transformation`)
    with its vertex 0 moved to the reflection of the vertex `2m + 1 \bmod n`
    of ``P``.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_surface_generators import _mirror
        sage: P = HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (1/2, 0), (1/3, 1/2)])
        sage: Q, gluings = _mirror(P, None, 1)
        sage: [v.coordinates(model="klein") for v in Q.vertices()]
        [(36/49, 12/49), (1/3, 1/2), (1/2, 0)]
        sage: gluings
        [((0, 0), (1, 2)), ((0, 1), (1, 1)), ((0, 2), (1, 0))]

    The reflection preserves the quadratic form of the plane of ``P``::

        sage: H = HyperbolicPlane(QQ, form=[1, 3, -1])
        sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (-1/4, 1/4)])
        sage: Q, _ = _mirror(P, None, 1)
        sage: [v.coordinates(model="klein") for v in Q.vertices()]
        [(4/17, 4/17), (-1/4, 1/4), (1/2, 0)]
        sage: P._cosh2_sides(), Q._cosh2_sides()
        ((4/3, 9/4, 4/3), (4/3, 9/4, 4/3))

    """
    from flatsurf.geometry.similarity_surface_generators import _pillowcase_gluings
    from flatsurf.geometry.hyperbolic_isometry import reflection

    n = len(P.vertices())
    m, gluings = _pillowcase_gluings(n, glued_edges, mirror_edge)

    image = P.apply_transformation(reflection(P.vertex(m), P.vertex((m + 1) % n)))
    Q = P.parent().polygon(
        vertices=[image.vertex((j - 2 * m - 1) % n) for j in range(n)]
    )
    return Q, gluings


def _standard_surface(surface, half_cuffs):
    r"""
    Return the image of ``surface`` in the hyperbolic plane of the form `x^2
    + y^2 - z^2`.

    INPUT:

    - ``surface`` -- an immutable hyperbolic surface in the plane of a form
      `x^2 + b y^2 - c z^2` over ``QQ`` or a number field `K`

    - ``half_cuffs`` -- a list of pairs ``(label, edge)`` of edges of
      ``surface`` whose lines may have their ends over `K`

    OUTPUT:

    The image of ``surface`` by the
    :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.standard_form_isometry`
    of its plane: without isotropic vector when `c` or `c/b` is a square in
    `K`, otherwise with the isotropic vector of an end of the line of the
    first of ``half_cuffs`` whose ends are over `K`, otherwise with
    ``algorithm="extension"``.

    EXAMPLES:

    In the plane of `x^2 + 2 y^2 - 3 z^2`, the line of the edge 2 of a
    triangle has its ends `(\pm 1, \pm 1)` over `\QQ`, the line of its edge
    0 has the ends `(\pm\sqrt{3}, 0)`::

        sage: from flatsurf import HyperbolicPlane, hyperbolic_surfaces
        sage: from flatsurf.geometry.hyperbolic_surface_generators import _standard_surface
        sage: H = HyperbolicPlane(QQ, form=[1, 2, -3])
        sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2)])
        sage: S = hyperbolic_surfaces.pillowcase(P)
        sage: _standard_surface(S, [(0, 0), (0, 2)]).base_ring()
        Rational Field
        sage: _standard_surface(S, [(0, 0)]).base_ring()
        Number Field in s with defining polynomial x^2 - 3 with s = 1.732050807568878?

    ALGORITHM:

    The ends of the line through the points of Klein lifts `P` and `Q` are
    `P + \lambda Q` with `\langle P, P\rangle + 2 \lambda \langle P,
    Q\rangle + \lambda^2 \langle Q, Q\rangle = 0`, over `K` exactly when
    `\langle P, Q\rangle^2 - \langle P, P\rangle \langle Q, Q\rangle` is a
    square in `K`.

    """
    H = surface.hyperbolic_plane()
    Q = H.quadratic_form()
    b, c = Q[1, 1], -Q[2, 2]
    if c.is_square() or (c / b).is_square():
        return H.standard_form_isometry()(surface)

    for label, edge in half_cuffs:
        vertices = surface.polygon(label).vertices()
        P, R = (
            vector(vertices[e % len(vertices)].coordinates(model="klein") + (1,))
            for e in (edge, edge + 1)
        )
        pp, pr, rr = P * Q * P, P * Q * R, R * Q * R
        discriminant = pr**2 - pp * rr
        if discriminant.is_square():
            end = P + (discriminant.sqrt() - pr) / rr * R
            return H.standard_form_isometry(isotropic_vector=end)(surface)

    return H.standard_form_isometry(algorithm="extension")(surface)


def _adjoin_root(K, f, approximation):
    r"""
    Return a field that contains ``K`` and the real root of ``f`` enclosed
    by ``approximation``.

    INPUT:

    - ``K`` -- ``QQ`` or a number field embedded in ``AA``

    - ``f`` -- a squarefree polynomial over ``K``

    - ``approximation`` -- a function that returns, for a precision, a real
      interval of this precision that contains the root; the intervals
      shrink to the root as the precision grows

    OUTPUT:

    A triple ``(L, embed, root)``: ``L`` is ``K`` when the root lies in
    ``K``, otherwise an absolute number field embedded in ``AA``; ``embed``
    is a function from ``K`` to ``L`` compatible with the embeddings, and
    ``root`` the root in ``L``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _adjoin_root
        sage: x = polygen(QQ)
        sage: L, embed, root = _adjoin_root(QQ, x^2 - 6, lambda prec: RealIntervalField(prec)(6).sqrt())
        sage: L, root
        (Number Field in a with defining polynomial x^2 - 6 with a = 2.449489742783178?, a)
        sage: M, embed, r = _adjoin_root(L, x^2 - (root + 3), lambda prec: RealIntervalField(prec)(root + 3).sqrt())
        sage: M
        Number Field in a with defining polynomial x^4 - 6*x^2 + 3 with a = 2.334414218338977?
        sage: embed(root), r^2 == embed(root) + 3
        (a^2 - 3, True)

    A root that lies in ``K``, here `-\sqrt{6}`, gives ``K`` back::

        sage: M, embed, r = _adjoin_root(L, x^2 - 6, lambda prec: -RealIntervalField(prec)(6).sqrt())
        sage: M is L, r
        (True, -a)

    TESTS:

    An approximation that encloses no root of ``f``::

        sage: _adjoin_root(QQ, x^2 - 6, lambda prec: RealIntervalField(prec)(3))
        Traceback (most recent call last):
        ...
        ValueError: no root of x^2 - 6 lies in the approximation

    ALGORITHM:

    The factors of ``f`` over ``K`` are evaluated in interval arithmetic on
    ``approximation``, with the precision doubled until a single one may
    vanish there, the factor of the root; when none may vanish, ``f`` has
    no root there. The extension of ``K`` by this
    factor is turned into an absolute field, embedded at the real root of
    its polynomial that encloses the value of its generator, computed from
    ``approximation`` and the embedding of ``K``.

    """
    from sage.rings.number_field.number_field import NumberField
    from sage.rings.real_mpfi import RealIntervalField

    def evaluate(coefficients, value, prec):
        R = RealIntervalField(prec)
        return sum(R(c) * value**k for k, c in enumerate(coefficients))

    factors = [h for h, _ in K["x"](f).factor()]
    prec = 64
    while True:
        value = approximation(prec)
        candidates = [
            h for h in factors if evaluate(h.list(), value, prec).contains_zero()
        ]
        if not candidates:
            # The interval evaluation of the factor of a root of f contains
            # zero at every precision.
            raise ValueError(f"no root of {f} lies in the approximation")
        if len(candidates) == 1:
            break
        prec *= 2
    h = candidates[0]

    if h.degree() == 1:
        return K, K, -h[0] / h[1]

    if K is QQ:
        # the generator of the absolute field is the root itself
        absolute, generator, root = h, [0, 1], [0, 1]

        def to_absolute(c):
            return c

    else:
        relative = K.extension(h, "b")
        A = relative.absolute_field("a")
        absolute_to_relative, relative_to_absolute = A.structure()
        absolute = A.polynomial()
        generator = absolute_to_relative(A.gen()).list()
        root = relative_to_absolute(relative.gen()).polynomial()

        def to_absolute(c):
            return relative_to_absolute(relative(c)).polynomial()

    roots = absolute.roots(AA, multiplicities=False)
    while True:
        value = evaluate(generator, approximation(prec), prec)
        close = [r for r in roots if RealIntervalField(prec)(r).overlaps(value)]
        if len(close) == 1:
            break
        prec *= 2
    L = NumberField(absolute, "a", embedding=close[0])

    def embed(c):
        return L(to_absolute(c))

    return L, embed, L(root)


def _composite(K, L):
    r"""
    Return the composite of the fields ``K`` and ``L`` in ``AA``.

    INPUT:

    - ``K``, ``L`` -- ``QQ`` or number fields embedded in ``AA``

    OUTPUT:

    A triple ``(M, embed_K, embed_L)``: ``M`` is ``K`` when it contains
    ``L``, ``L`` when ``K`` is ``QQ``, otherwise a number field embedded in
    ``AA``; ``embed_K`` and ``embed_L`` are functions from ``K`` and ``L``
    to ``M`` compatible with the embeddings.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _composite
        sage: x = polygen(QQ)
        sage: K.<a> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: L.<b> = NumberField(x^2 - 8, embedding=-AA(8).sqrt())
        sage: M, embed_K, embed_L = _composite(K, L)
        sage: M is K, embed_L(b)
        (True, -2*a)
        sage: L.<b> = NumberField(x^2 - 3, embedding=AA(3).sqrt())
        sage: M, embed_K, embed_L = _composite(K, L)
        sage: M.degree(), AA(embed_K(a)) == AA(2).sqrt(), AA(embed_L(b)) == AA(3).sqrt()
        (4, True, True)

    """
    from sage.rings.real_mpfi import RealIntervalField

    if L is QQ:
        return K, K, K
    if K is QQ:
        return L, L, L

    M, embed, root = _adjoin_root(
        K, L.polynomial(), lambda prec: RealIntervalField(prec)(L.gen())
    )
    return M, embed, lambda c: c.polynomial()(root)


def _coordinate_field(K, elements):
    r"""
    Return the field generated by ``elements`` with a positive generator,
    and the elements in it.

    INPUT:

    - ``K`` -- ``QQ`` or a number field embedded in ``AA``

    - ``elements`` -- a list of elements of ``K``

    OUTPUT:

    A pair ``(F, elements)``: ``F`` is ``QQ`` or a number field embedded in
    ``AA`` whose generator is positive, the subfield of ``K`` generated by
    ``elements``, whose defining polynomial is reduced as in ``QQbar``;
    ``elements`` is the list of the elements in ``F``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _coordinate_field
        sage: x = polygen(QQ)
        sage: K.<a> = NumberField(x^4 - 6*x^2 + 3, embedding=AA.polynomial_root(x^4 - 6*x^2 + 3, RIF(2.3, 2.4)))
        sage: _coordinate_field(K, [a^2 - 3, 1/2])
        (Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?, [a, 1/2])
        sage: _coordinate_field(K, [3 - a^2])
        (Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?, [-a])
        sage: _coordinate_field(K, [QQ(1), QQ(2)])
        (Rational Field, [1, 2])

    The field of `(1 - \sqrt{5}) / 2` is generated by the root `(1 +
    \sqrt{5}) / 2` of the same polynomial::

        sage: K.<a> = NumberField(x^2 - x - 1, embedding=(1 - AA(5).sqrt()) / 2)
        sage: _coordinate_field(K, [a])
        (Number Field in a with defining polynomial y^2 - y - 1 with a = 1.618033988749895?, [-a + 1])

    ALGORITHM:

    The generator `\gamma` of the subfield, scaled to an algebraic integer,
    is reduced by ``do_polred`` of ``QQbar`` (PARI's ``polredbest`` when
    the coefficients are below its cost threshold, otherwise unchanged)
    into `\delta` of polynomial `g`. When
    `\delta` is negative, the generator is a positive root of `g` in
    `\QQ(\delta)` if there is one, which keeps the polynomial, otherwise
    `-\delta`. The elements are expressed in the powers of the generator by
    linear algebra over `\QQ`.

    """
    from sage.matrix.constructor import matrix
    from sage.rings.number_field.number_field import NumberField
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.qqbar import do_polred
    from sage.rings.real_mpfi import RealIntervalField

    from flatsurf.geometry.subfield import subfield_from_elements

    if all(c in QQ for c in elements):
        return QQ, [QQ(c) for c in elements]

    _, _, inclusion = subfield_from_elements(K, elements, polred=False)
    gamma = inclusion(inclusion.domain().gen())
    p = gamma.minpoly()
    d = p.degree()
    denominator = ZZ.one()
    for c in p.list():
        denominator = denominator.lcm(c.denominator())
    fwd, _, g = do_polred(p(p.parent().gen() / denominator) * denominator**d)
    delta = fwd(denominator * gamma)

    if delta < 0:
        F = NumberField(g, "c")
        positive = [
            r.polynomial()(delta)
            for r in g.roots(F, multiplicities=False)
            if r.polynomial()(delta) > 0
        ]
        if positive:
            delta = positive[0]
        else:
            delta = -delta
            g = g(-g.parent().gen()) * (-1) ** d

    powers = matrix(QQ, [(delta**k).vector() for k in range(d)])
    coordinates = powers.solve_left(matrix(QQ, [K(c).vector() for c in elements]))

    roots = g.roots(AA, multiplicities=False)
    prec = 64
    while True:
        close = [
            r
            for r in roots
            if RealIntervalField(prec)(r).overlaps(RealIntervalField(prec)(delta))
        ]
        if len(close) == 1:
            break
        prec *= 2
    F = NumberField(PolynomialRing(QQ, "y")(g.list()), "a", embedding=close[0])
    return F, [F(list(row)) for row in coordinates.rows()]


def _develop(K, chain, half_cuffs, placed, i, plane, form, vertices):
    r"""
    Return the Klein coordinates of the front polygon of the pants ``i`` of
    :meth:`HyperbolicSurfaceGenerators.from_signature`, developed into the
    seam plane of the pants 0, and the field that holds them.

    INPUT:

    - ``K`` -- ``QQ`` or a number field embedded in ``AA`` that contains
      the coordinates and the forms below

    - ``chain``, ``half_cuffs`` -- the slots of the pants and the half-cuff
      edges of their slots, as in
      :meth:`HyperbolicSurfaceGenerators.from_signature`

    - ``placed`` -- the list of the Klein coordinates of the front polygons
      of the pants before ``i``, in the seam plane of the pants 0

    - ``i`` -- the index of the pants

    - ``plane`` -- the pair ``(b, c)`` of the seam plane of the pants 0

    - ``form`` -- the pair ``(b, c)`` of the seam plane of the pants ``i``

    - ``vertices`` -- the Klein coordinates of the front polygon of the
      pants ``i`` in its seam plane

    OUTPUT:

    A triple ``(L, embed, developed)``: ``L`` is ``K`` or an extension of
    ``K`` by a square root, embedded in ``AA``, ``embed`` the function from
    ``K`` to ``L``, and ``developed`` the list of the Klein coordinates in
    ``L`` of the developed polygon.

    EXAMPLES:

    The two pants of the closed surface of genus 2 share the cuff of their
    slots 2 and 0, the edge 5 of the first hexagon and the edge 1 of the
    second one::

        sage: from flatsurf.geometry.hyperbolic_surface_generators import _develop, _pants_polygon
        sage: chain = [[("cuff", "l"), ("cuff", "l"), ("cuff", "c")], [("cuff", "c"), ("cuff", "m"), ("cuff", "m")]]
        sage: half_cuffs = {(0, 0): 1, (0, 1): 3, (0, 2): 5, (1, 0): 1, (1, 1): 3, (1, 2): 5}
        sage: K, form, vertices, _ = _pants_polygon((21/4, 21/4, 21/4), ())
        sage: _, _, first = _develop(K, chain, half_cuffs, [], 0, form, form, vertices)
        sage: L, _, second = _develop(K, chain, half_cuffs, [first], 1, form, form, vertices)
        sage: L
        Rational Field
        sage: second[1] == first[0] and second[2] == first[5]
        True

    ALGORITHM:

    The first slot of the pants ``i`` is a cuff shared with an earlier
    pants. The isometry `M` sends the half-cuff edge `pq` of the pants
    ``i`` onto the matching edge `p'q'` of the earlier front polygon,
    reversed, from the plane of the form `Q_i` to the one of `Q_0`, with
    `M^T Q_0 M = \lambda Q_i`. On the Klein lifts, `M` sends `p` to `p'`,
    the tangent `t = q - (\langle p, q\rangle / \langle p, p\rangle) p` to
    `\mu t'`, and the normal `n = Q_i^{-1} (p \times q)` of the line `pq` to
    `\nu n'`. Then `\lambda = \langle p', p'\rangle / \langle p,
    p\rangle`, and `\mu = \kappa / \kappa'`, with `\kappa = \langle p, q
    \rangle / \langle p, p\rangle`, sends `q` to a multiple of `q'`, the two
    edges having the same length. The scale `\nu = \sqrt{\lambda \langle n,
    n\rangle / \langle n', n'\rangle}` is positive, so that `M` sends the
    side of the polygon ``i`` to the other side of the edge `p'q'`; it is
    the only square root, adjoined to ``K`` when it is not there.

    """
    from sage.rings.real_mpfi import RealIntervalField

    if i == 0:
        return K, K, vertices

    j, s = next(
        (j, s) for j in range(i) for s in range(3) if chain[j][s] == chain[i][0]
    )
    target = placed[j]
    a = half_cuffs[j, s]
    e = half_cuffs[i, 0]

    def lift(point):
        return vector(K, (point[0], point[1], 1))

    def frame(Q, p, q):
        t = q - (p * Q * q) / (p * Q * p) * p
        return t, Q.inverse() * p.cross_product(q)

    Qi = diagonal_matrix(K, [1, form[0], -form[1]])
    Q0 = diagonal_matrix(K, [1, plane[0], -plane[1]])
    p, q = lift(vertices[e]), lift(vertices[(e + 1) % len(vertices)])
    p0, q0 = lift(target[(a + 1) % len(target)]), lift(target[a])
    t, n = frame(Qi, p, q)
    t0, n0 = frame(Q0, p0, q0)

    scale = (p0 * Q0 * p0) / (p * Qi * p)
    mu = ((p * Qi * q) / (p * Qi * p)) / ((p0 * Q0 * q0) / (p0 * Q0 * p0))
    nu2 = scale * (n * Qi * n) / (n0 * Q0 * n0)

    x = K["x"].gen()
    L, embed, nu = _adjoin_root(
        K, x**2 - nu2, lambda prec: RealIntervalField(prec)(nu2).sqrt()
    )

    def column(v, factor=1):
        return [factor * embed(c) for c in v]

    source = matrix(L, [column(p), column(t), column(n)]).transpose()
    image = matrix(L, [column(p0), column(mu * t0), column(n0, nu)]).transpose()
    M = image * source.inverse()

    developed = []
    for point in vertices:
        X, Y, Z = M * vector(L, (embed(point[0]), embed(point[1]), 1))
        developed.append((X / Z, Y / Z))
    return L, embed, developed


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
        the reflection of the edge `2m - j \bmod n` of ``P``. The polygon
        `Q` is the image of ``P`` under the reflection, as given by
        :meth:`~flatsurf.geometry.categories.polygons.Polygons.ParentMethods.apply_transformation`,
        with its vertex 0 moved to the reflection of the vertex `2m + 1 \bmod
        n`, so that the edge `m` of `Q` is the edge `m` of ``P``. The edge `e`
        of ``P`` is glued to the edge `2m - e \bmod n` of `Q` for every `e`
        in ``glued_edges``, by the reflection in the line of `e` followed by
        the reflection in the line of `m`; the edge `m` of ``P`` is glued to
        the edge `m` of `Q`, which is the same segment of the plane, by the
        identity. The other edges of ``P`` and of `Q` form the boundary of
        the surface.

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
        in the hyperbolic plane of ``P``, over its base ring: a reflection in
        a line through two points has rational coefficients in their Klein
        coordinates and in the coefficients of the quadratic form.

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

            sage: from flatsurf.geometry.hyperbolic_isometry import reflection
            sage: vertices = [v.coordinates(model="klein") for v in P.vertices()]
            sage: R = [reflection(vertices[e], vertices[(e + 1) % 3]) for e in range(3)]
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

        A triangle of the plane with the form `x^2 + 3y^2 - z^2`, reflected
        in its edge 1 with that form::

            sage: H = HyperbolicPlane(QQ, form=[1, 3, -1])
            sage: C = H.polygon(vertices=[(0, 0), (1/2, 0), (-1/4, 1/4)])
            sage: T = hyperbolic_surfaces.pillowcase(C, mirror_edge=1)
            sage: T.hyperbolic_plane()
            Hyperbolic Plane over Rational Field with form x^2 + 3*y^2 - z^2
            sage: [v.coordinates(model="klein") for v in T.polygon(1).vertices()]
            [(4/17, 4/17), (-1/4, 1/4), (1/2, 0)]
            sage: sorted(T.angles(numerical=True))  # abs tol 1e-12
            [0.147583617650433, 0.147583617650433, 0.666666666666667]
            sage: T.holonomy_representation().invariant_trace_field()
            Rational Field
            sage: TestSuite(T).run()

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

        surface = MutableOrientedHyperbolicSurface(P.parent())
        surface.add_polygon(P, label=0)
        surface.add_polygon(Q, label=1)
        for x, y in gluings:
            surface.glue(x, y)
        surface.set_immutable()
        return surface

    @staticmethod
    def pants(boundary_traces=(), cone_orders=(), base_ring=None, standard=False):
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

        The front polygon, hence its base field and its coordinates, depends
        on the order of ``cone_orders`` and of ``boundary_traces``, which
        decides what slot each takes. The surfaces of two orders are
        isometric, and the signature of the surface does not depend on the
        order.

        INPUT:

        - ``boundary_traces`` -- a sequence of real algebraic numbers larger
          than 3 (default: ``()``); the traces in `\mathrm{SO}(2, 1)` of the
          holonomy of the cuffs. A cuff of trace `t` has length `\ell` with
          `\cosh(\ell / 2) = \sqrt{t + 1} / 2`.

        - ``cone_orders`` -- a sequence of rational numbers in `(-1, 0)`
          (default: ``()``); the orders `k` of the cone points, of angle
          `2\pi(1 + k)`

        - ``base_ring`` -- a field or ``None`` (default: ``None``); with
          ``None``, the field `K` below, otherwise an exact field into which
          `K` converts along the embeddings in the reals, such as ``AA``,
          which receives the form and the coordinates

        - ``standard`` -- a boolean (default: ``False``); whether to carry
          the surface to the plane of the form `x^2 + y^2 - z^2`, over `K`
          when that plane has an isometry to it over `K`, otherwise over a
          quadratic extension of `K`

        OUTPUT:

        An immutable
        :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`
        with the labels 0, the front polygon, and 1, its mirror image. Its
        :meth:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces.FiniteType.Oriented.Connected.ParentMethods.signature`
        has the genus 0, the cone orders ``cone_orders`` and the boundary
        traces ``boundary_traces``. The surface lies in the hyperbolic plane
        of a form `x^2 + b y^2 - c z^2` of its own, its seam plane, over the
        field `K = \QQ(u_0, u_1, u_2, p)` of its slots, with `u_k = t_k + 1`
        for a cuff of trace `t_k`, see the algorithm below; over `\QQ`, `b`
        and `c` are squarefree integers. The front polygon has its vertex 0
        at the origin and its edge 0 on the positive `x`-axis.

        EXAMPLES:

        The pants with three cuffs of trace `21/4` is defined over `\QQ`, in
        its seam plane::

            sage: from flatsurf import hyperbolic_surfaces
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4))
            sage: S
            Hyperbolic Surface with boundary built from 2 right-angled hexagons
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field with form x^2 + 6*y^2 - 6*z^2
            sage: S.polygon(0).vertices()
            ((0, 0), (12/5, 0), (12/5, 3/25), (60/29, 15/29), (48/25, 3/5), (0, 3/5))
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(21/4, 21/4, 21/4))
            sage: S.area()
            1
            sage: TestSuite(S).run()

        In the plane of the form `x^2 + y^2 - z^2`, still over `\QQ` since
        `c/b` is a square::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), standard=True)
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field
            sage: S.polygon(0).vertices()
            ((-23/25, 0),
             (-23/25, -48/125),
             (-125/139, -60/139),
             (-73/95, -12/19),
             (-5/7, -24/35),
             (-5/7, 0))
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(21/4, 21/4, 21/4))

        The cuffs of trace `5` give a seam plane over `\QQ(\sqrt{6})` that
        has no isometry to the plane of `x^2 + y^2 - z^2` over that field;
        ``standard=True`` adjoins a square root::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(5, 5, 5))
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178? with form ...
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(5, 5, 5), standard=True)
            sage: S.base_ring()
            Number Field in s with defining polynomial ... with s = ...
            sage: S.base_ring().degree()
            4
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(5, 5, 5))

        A cone point of angle `\pi` and two cuffs, from a pentagon::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(QQ(21)/4, QQ(21)/4), cone_orders=(QQ(-1)/2,))
            sage: S
            Hyperbolic Surface with boundary built from 2 right-angled pentagons
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field with form x^2 + 34*y^2 - 34*z^2
            sage: S.polygon(0).vertices()
            ((0, 0), (5, 0), (5, 15/34), (544/125, 3/5), (0, 3/5))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-1/2,), boundary_traces=(21/4, 21/4))
            sage: S.area()
            1/2

        Two cone points and a cuff, from a quadrilateral::

            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4,), cone_orders=(-1/2, -2/3))
            sage: S
            Hyperbolic Surface with boundary built from 2 quadrilaterals
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field with form x^2 + 13*y^2 - 13*z^2
            sage: S.polygon(0).vertices()
            ((0, 0), (2, 0), (2, 3/5), (0, 3/5))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-2/3, -1/2), boundary_traces=(21/4,))
            sage: S.area()
            1/6

        Three cone points, from a triangle, give a closed orbifold of genus
        0, in the plane of the invariant form of the triangle::

            sage: S = hyperbolic_surfaces.pants(cone_orders=(-1/2, -3/4, -7/8))
            sage: S
            Genus 0 Hyperbolic Surface built from 2 right triangles
            sage: S.base_ring().polynomial()
            x^2 - 2
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? with form ...
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
            y^2 - y - 1
            sage: S.signature()
            SurfaceSignature(g=0, boundary_traces=(4, 4, 4))
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 20881/144))
            sage: S.base_ring()
            Rational Field
            sage: [v.coordinates(model="klein") for v in S.polygon(0).vertices()]
            [(0, 0), (56/65, 0), (56/65, 99/325), (840/4121, 4015/4121), (1344/9425, 143/145), (0, 143/145)]

        The field `K` of the slots, and the way to the plane of the form
        `x^2 + y^2 - z^2`, the isotropic vector `(0, \sqrt{c/b}, 1)` when
        `c/b` is a square, otherwise an end of the line of a half-cuff when
        one is over `K`, otherwise a quadratic extension::

            sage: def fields(traces, orders=()):
            ....:     S = hyperbolic_surfaces.pants(boundary_traces=traces, cone_orders=orders)
            ....:     T = hyperbolic_surfaces.pants(boundary_traces=traces, cone_orders=orders, standard=True)
            ....:     return S.base_ring(), S.area(), T.base_ring(), T.signature() == S.signature()
            sage: fields((4, 5, 21/4))
            (Number Field in a with defining polynomial y^2 - 30 with a = ...,
             1,
             Number Field in a with defining polynomial y^2 - 30 with a = ...,
             True)
            sage: fields((21/4, 4, 5))
            (Number Field in a with defining polynomial y^2 - 30 with a = ...,
             1,
             Number Field in a with defining polynomial y^2 - 30 with a = ...,
             True)
            sage: fields((21/4, 21/4), (-1/2,))
            (Rational Field, 1/2, Rational Field, True)
            sage: K, area, L, _ = fields((4, 5), (-2/3,))
            sage: K.polynomial(), area, L.degree()
            (y^2 - 30, 2/3, 4)
            sage: fields((4,), (-2/3, -4/5))
            (Number Field in a with defining polynomial y^2 - y - 1 with a = ...,
             7/15,
             Number Field in a with defining polynomial y^2 - y - 1 with a = ...,
             True)
            sage: fields((8, 8, 8))
            (Rational Field,
             1,
             Number Field in s with defining polynomial x^2 - 2 with s = 1.414213562373095?,
             True)
            sage: K, area, L, _ = fields((15, 15, 15))
            sage: K, area, L.degree()
            (Rational Field, 1, 2)
            sage: fields((), (-1/2, -3/4, -7/8))[0].polynomial()
            x^2 - 2
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
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), base_ring=AA)
            sage: S.base_ring()
            Algebraic Real Field
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Algebraic Real Field with form x^2 + 6*y^2 - 6*z^2
            sage: S = hyperbolic_surfaces.pants(boundary_traces=(5, 5, 5), base_ring=AA, standard=True)
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Algebraic Real Field
            sage: hyperbolic_surfaces.pants(boundary_traces=(5, 5, 5), base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field
            sage: hyperbolic_surfaces.pants(boundary_traces=(21/4, 21/4, 21/4), base_ring=RR)
            Traceback (most recent call last):
            ...
            ValueError: the ring must be exact, not Real Field with 53 bits of precision

        ALGORITHM:

        The front polygon is
        :meth:`~HyperbolicPolygonsConstructor.right_angled_hexagon` without
        cone point,
        :meth:`~HyperbolicPolygonsConstructor.right_angled_pentagon` with
        one, and :meth:`~HyperbolicPolygonsConstructor.right_angled_quadrilateral`
        with two. With three cone points, the front polygon is
        :meth:`~HyperbolicPolygonsConstructor.triangle`.

        The mirror `k`, the seam opposite the slot `k`, has the unit normal
        `e_k`; the mirrors `i` and `j` meet at the cone point of the slot `k`
        and of order `\kappa_k` with the angle `\theta_k = \pi(1 +
        \kappa_k)`, or are disjoint with the
        half-cuff of the slot `k` as common perpendicular, of length `\ell_k
        / 2`. The Gram matrix of the `e_k` has the entries `1` on the
        diagonal and `-c_k` for the mirrors `i`, `j`, with `c_k = \cos
        \theta_k` or `c_k = \cosh(\ell_k / 2) = \sqrt{t_k + 1} / 2`. Vinberg's
        rescaling of the normals turns it into a Gram matrix `G` over the
        field of `u_k = 4 c_k^2`, `t_k + 1` for a cuff, and of `p = 8 c_0 c_1
        c_2`. In the rescaled basis `f`, a cone point is a column `v_k` of
        `G^{-1}`; for a cuff, that column is the pole of the line of the
        cuff, and the corner of the mirror `a` and the cuff `k` is the cross
        product of the row `a` of `G` and of the `k`-th vector of the
        canonical basis. In the orthogonal basis of the normal `f_1` of the
        seam 20, of a vector orthogonal to it in the plane of the cuff or the
        cone point of the slot 2, and of the vertex 0, the form is diagonal.

        With ``standard``, the isometry to the plane of `x^2 + y^2 - z^2`
        is the one of
        :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.standard_form_isometry`,
        found from the square roots of `c` or `c/b` when they are in `K`,
        otherwise from an ideal end of the line of a half-cuff when one has
        its coordinates in `K`, otherwise over a quadratic extension of `K`.

        """
        _check_exact(base_ring)
        boundary_traces = tuple(boundary_traces)
        K, form, vertices, seams = _pants_polygon(boundary_traces, cone_orders)
        if standard:
            P = _plane_polygon(K, form, vertices, None)
            S = HyperbolicSurfaceGenerators.pillowcase(P, glued_edges=seams)
            n = len(vertices)
            half_cuffs = [e for e in range(n) if e not in seams]
            S = _standard_surface(
                S,
                [(0, e) for e in half_cuffs]
                + [(1, e) for e in sorted(-e % n for e in half_cuffs)],
            )
            return S if base_ring is None else S.change_ring(base_ring)
        P = _plane_polygon(K, form, vertices, base_ring)
        return HyperbolicSurfaceGenerators.pillowcase(P, glued_edges=seams)

    @staticmethod
    def regular_octagon(gluing="opposite", base_ring=None):
        r"""
        Return the hyperbolic surface of genus 2 built from the regular
        octagon with interior angles `\pi/4`.

        INPUT:

        - ``gluing`` -- a string (default: ``"opposite"``); ``"opposite"``
          glues each edge `e` to the edge `e + 4`, ``"katok"`` glues the
          edges `0, 1, 4, 5` to the edges `2, 3, 6, 7`, the identification
          `a_1 b_1 a_1^{-1} b_1^{-1} a_2 b_2 a_2^{-1} b_2^{-1}` of Katok,
          Fuchsian groups (1992), Section 4.3, Example C, p. 95

        - ``base_ring`` -- a ring or ``None`` (default: ``None``), the base ring
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
            sage: hyperbolic_surfaces.regular_octagon(base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field

        """
        from flatsurf.geometry.hyperbolic_surface import (
            MutableOrientedHyperbolicSurface,
        )

        _check_exact(base_ring)
        if gluing == "opposite":
            pairs = [(0, 4), (1, 5), (2, 6), (3, 7)]
        elif gluing == "katok":
            pairs = [(0, 2), (1, 3), (4, 6), (5, 7)]
        else:
            raise ValueError(f'gluing must be "opposite" or "katok", not {gluing!r}')

        P = hyperbolic_polygons.regular_ngon(8, QQ((1, 8)), base_ring=base_ring)

        surface = MutableOrientedHyperbolicSurface(P.base_ring())
        surface.add_polygon(P, label=0)
        for e, f in pairs:
            surface.glue((0, e), (0, f))
        surface.set_immutable()
        return surface

    @staticmethod
    def from_signature(signature, cuff_traces=None, base_ring=None, standard=False):
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

        - ``base_ring`` -- a ring or ``None`` (default: ``None``), the base ring
          of the surface: with ``None``, the number field generated by the
          Klein coordinates of all the polygons and by the form of their
          plane, ``QQ`` when they are rational; ``AA`` keeps the coordinates
          in ``AA``; any other ring receives the coordinates and the form by
          conversion

        - ``standard`` -- a boolean (default: ``False``); whether to carry
          the surface to the plane of the form `x^2 + y^2 - z^2`, as in
          :meth:`pants`

        OUTPUT:

        An immutable
        :class:`~flatsurf.geometry.hyperbolic_surface.MutableOrientedHyperbolicSurface`
        whose
        :meth:`~flatsurf.geometry.categories.cone_surfaces.ConeSurfaces.FiniteType.Oriented.Connected.ParentMethods.signature`
        is ``signature``. The pair of pants `i` of the decomposition below
        is made of the polygons with the labels ``(i, 0)``, its front
        polygon, and ``(i, 1)``, its mirror image, as in :meth:`pants`. The
        surface lies in the seam plane of the pants 0 of :meth:`pants`.

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

        The surface lies in the seam plane of its first pants; the second
        front polygon is glued to the first one by the identity::

            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field with form x^2 + 6*y^2 - 6*z^2
            sage: S.edge_transformation((0, 0), 5).matrix().is_one()
            True
            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), standard=True)
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Rational Field
            sage: S.signature()
            SurfaceSignature(g=2)

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

        Three cone points give the pillowcase of a triangle, in the plane of
        its invariant form, see :meth:`pants`; for the angles `2\pi/m`, the
        triangle of the triangle group::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (3, 3, 4)))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-3/4, -2/3, -2/3))
            sage: S.base_ring().polynomial()
            x^2 - 2
            sage: S.area()
            1/12

        ::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (2, 3, 7)))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-6/7, -2/3, -1/2))
            sage: S.base_ring().polynomial()
            x^3 - 7*x - 7
            sage: S.area()
            1/42

        Cone angles that are not of the form `2\pi/m`::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature(0, cone_orders=(-1/3, -5/6, -7/8)))
            sage: S.signature()
            SurfaceSignature(g=0, cone_orders=(-7/8, -5/6, -1/3))
            sage: S.area()
            1/24

        With the default traces, the orbifold of signature `(2; 2)` is over
        its invariant trace field, a quadratic field. Traces `99` give a
        surface over ``QQ``::

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

        Closed surfaces of genus 2 with other traces, the second pants
        developed across the cuff of trace `c` in ``cuff_traces=(l, c,
        l')``; the base ring is the invariant trace field::

            sage: def check(traces):
            ....:     S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=traces)
            ....:     assert S.edge_transformation((0, 0), 5).matrix().is_one()
            ....:     assert S.signature() == SurfaceSignature(2) and S.area() == 2
            ....:     K = S.base_ring()
            ....:     assert K.degree() == S.holonomy_representation().invariant_trace_field().degree()
            ....:     TestSuite(S).run()
            ....:     return K
            sage: check((21/4, 21/4, 21/4))
            Rational Field
            sage: check((4, 4, 4)).polynomial()
            y^2 - y - 1

        With the cuff traces `4, 5, 21/4`, the two seam planes are over
        `\QQ(\sqrt{6})` and the development across the cuff of trace `5`
        adjoins a square root::

            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(4, 5, 21/4))
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 552*y^2 - 1200*y + 54696 with a = 9.696518361151262?
            sage: S.holonomy_representation().invariant_trace_field()
            Number Field in a with defining polynomial y^4 - 552*y^2 - 1200*y + 54696 with a = 9.696518361151262?
            sage: S.edge_transformation((0, 0), 5).matrix().is_one()
            True
            sage: TestSuite(S).run()

        The number of cuff traces, their values, and an explicit ring
        that does not hold the coordinates::

            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(21/4,))
            Traceback (most recent call last):
            ...
            ValueError: cuff_traces must have 3 traces, one per interior cuff, not 1
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(21/4, 21/4, 2))
            Traceback (most recent call last):
            ...
            ValueError: the cuff trace 2 must be larger than 3
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(2, (2,)), base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), base_ring=AA).base_ring()
            Algebraic Real Field
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), base_ring=AA, standard=True).hyperbolic_plane()
            Hyperbolic Plane over Algebraic Real Field
            sage: S = hyperbolic_surfaces.from_signature(SurfaceSignature.from_orbifold_orders(0, (3, 3, 4)), standard=True)
            sage: S.hyperbolic_plane()
            Hyperbolic Plane over Number Field in s with defining polynomial ... with s = ...
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(1, boundary_traces=(None,)))
            Traceback (most recent call last):
            ...
            ValueError: the boundary traces of the signature must be given, not None
            sage: hyperbolic_surfaces.from_signature((2,))
            Traceback (most recent call last):
            ...
            TypeError: signature must be a SurfaceSignature
            sage: hyperbolic_surfaces.from_signature(SurfaceSignature(2), base_ring=RR)
            Traceback (most recent call last):
            ...
            ValueError: the ring must be exact, not Real Field with 53 bits of precision

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
        which the pants of three cuffs is over `\QQ` in its seam plane, or
        the trace of ``cuff_traces``, given in the order in which the
        interior cuffs first appear along the pants of the chain and their
        slots: on `(2; 2)` the loop `\ell`, the links `c_1, c_2`, the loop
        `\ell'`.

        Each front polygon is built in the seam plane of its pants, see
        :meth:`pants`. The pants 0 stays there; the first slot of each
        other pants of the chain is a cuff shared with an earlier pants,
        across which its front polygon is developed: the isometry that sends
        its half-cuff edge onto the matching edge of the front polygon
        already placed, reversed, carries it into the seam plane of the
        pants 0, so that the two polygons lie on either side of the cuff and
        are glued by the identity. The coordinates are computed in a number
        field embedded in the reals, the composite of the fields of the seam
        planes, to which each development adjoins the single square root it
        needs: the isometry sends the lift of a vertex to the lift of its
        image, the tangent of the edge to a multiple in the field of its
        image, and the normal of the edge to the square root of a ratio of
        their norms times its image.

        The other cuffs, the loops and the second cuff shared by the two
        pants of a handle, join two front polygons already in the seam plane
        of the pants 0, of form `Q_0`. With `\langle u, v\rangle = u^T Q_0
        v`, the isometry that glues the edge `[A, B]` to the edge `[C, D]`,
        `A` to `D` and `B` to `C`, sends the lifts `A` to `\alpha D` and `B`
        to `\beta C`, and the normal of `[A, B]` to `\alpha\beta` times the
        normal of `[D, C]`, with `\alpha^2 = \langle A, A\rangle / \langle
        D, D\rangle` and `\alpha\beta = \langle A, B\rangle / \langle D,
        C\rangle` in the field of the coordinates. The side pairing is over
        a field that holds the coordinates if and only if it holds the
        positive square root `\alpha`, which is adjoined when it is not
        there.

        With ``base_ring=None``, the polygons are built over the subfield
        generated by all their Klein coordinates, by the form of the seam
        plane of the pants 0 and by the square roots `\alpha` of the side
        pairings, its polynomial reduced as in ``QQbar`` (by PARI's
        ``polredbest`` when its coefficients are small enough) and its
        generator positive, so that every side pairing is defined over it.
        An explicit ``base_ring`` that does not hold these numbers is
        refused.

        Genus 0 with three cone points and no boundary gives the single
        pants of three cone points, the pillowcase of a triangle of
        :meth:`HyperbolicPolygonsConstructor.triangle` in the plane of its
        invariant form, over its field or over ``base_ring``.

        """
        from sage.rings.real_mpfi import RealIntervalField

        from flatsurf.geometry.hyperbolic import HyperbolicPlane
        from flatsurf.geometry.hyperbolic_surface import (
            MutableOrientedHyperbolicSurface,
        )
        from flatsurf.geometry.signature import SurfaceSignature

        _check_exact(base_ring)
        if not isinstance(signature, SurfaceSignature):
            raise TypeError("signature must be a SurfaceSignature")

        def check_trace(t, name):
            try:
                trace = AA(t)
            except (TypeError, ValueError):
                raise ValueError(
                    f"the {name} {t} must be a real algebraic number"
                ) from None
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

        triangle = genus == 0 and not boundary_traces and len(cone_orders) == 3

        # A leg is stored as ("leg", kind), with kind ("boundary", trace),
        # ("fold", None) or ("cone", order); a slot of a pants is a leg or a
        # ("cuff", name).
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

        # Each cuff name is in two slots.
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

        # The surface is built over its own field before it is carried to
        # the standard plane.
        ring = None if standard else base_ring

        # The half-cuff edge of the front polygon of each slot that is not a
        # cone point, and the pairs (i, e) of a front polygon and its
        # half-cuff edge at the two ends of each interior cuff.
        half_cuffs = {}
        cuff_ends = {}
        if triangle:
            angles = [(1 + slot[1][1]) / 2 for slot in chain[0]]
            polygons = [
                hyperbolic_polygons.triangle(
                    angles[2], angles[0], angles[1], base_ring=ring
                )
            ]
            seams = [[0, 1, 2]]
        else:
            # A field that holds the seam fields of the pants met so far,
            # the coordinates of the front polygons placed so far and the
            # form of the seam plane of the pants 0.
            # TODO: Compute each development in the composite of the fields
            # it needs only, or pass to the subfield of the coordinates
            # placed so far before taking the next composite. The composite
            # doubles its degree with each new quadratic seam field or
            # square root, and every later operation runs in it, though the
            # coordinates often generate a much smaller subfield.
            field = QQ
            fronts = []
            seams = []
            plane = None

            def move(embed):
                return [[(embed(x), embed(y)) for x, y in front] for front in fronts], [
                    embed(t) for t in plane
                ]

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
                    K, form, vertices, pants_seams = _pants_polygon(
                        pants_traces, orders
                    )
                except ValueError as error:
                    raise NotImplementedError(
                        f"a pair of pants of the decomposition into pants of this signature does not exist, the surface may exist without such a decomposition: {error}"
                    ) from error
                edges = [e for e in range(len(vertices)) if e not in pants_seams]
                for s, e in zip(cuffs, edges):
                    half_cuffs[i, s] = e
                field, embed, embed_K = _composite(field, K)
                form = [embed_K(t) for t in form]
                vertices = [(embed_K(x), embed_K(y)) for x, y in vertices]
                if i == 0:
                    plane = form
                else:
                    fronts, plane = move(embed)
                field, embed, front = _develop(
                    field, chain, half_cuffs, fronts, i, plane, form, vertices
                )
                fronts, plane = move(embed)
                fronts.append(front)
                seams.append(pants_seams)

            # Each interior cuff joins the half-cuff edges of two front
            # polygons.
            for (i, s), e in half_cuffs.items():
                if chain[i][s][0] == "cuff":
                    cuff_ends.setdefault(chain[i][s][1], []).append((i, e))

            # The isometry of the side pairing of the edge [A, B] of the
            # front polygon i with the edge [C, D] of the front polygon j
            # sends A to alpha D and B to beta C on the Klein lifts, with
            # alpha^2 = <A, A> / <D, D> and alpha beta = <A, B> / <D, C>,
            # and the normal of [A, B] to alpha beta times the normal of
            # [D, C]. It is over the field if and only if alpha is. The root
            # alpha is positive: both norms are negative, the points being
            # inside the disk, and the isometry preserves the future cone.
            alphas = []
            for (i, a), (j, b) in cuff_ends.values():
                Q0 = diagonal_matrix(field, [1, plane[0], -plane[1]])
                A = vector(field, (*fronts[i][a], 1))
                D = vector(field, (*fronts[j][(b + 1) % len(fronts[j])], 1))
                ratio = (A * Q0 * A) / (D * Q0 * D)
                field, embed, alpha = _adjoin_root(
                    field,
                    field["x"].gen() ** 2 - ratio,
                    lambda prec, ratio=ratio: RealIntervalField(prec)(ratio).sqrt(),
                )
                fronts, plane = move(embed)
                alphas = [embed(c) for c in alphas] + [alpha]

            coordinates = (
                [c for vertices in fronts for vertex in vertices for c in vertex]
                + plane
                + alphas
            )
            field, coordinates = _coordinate_field(field, coordinates)
            if ring is not None:
                convert = _field_embedding(field, ring)
                field = ring
                coordinates = [convert(c) for c in coordinates]
            coordinates = coordinates[: len(coordinates) - len(alphas)]
            H = HyperbolicPlane(field, form=[1, coordinates[-2], -coordinates[-1]])
            coordinates = iter(coordinates)
            polygons = [
                H.polygon(
                    vertices=[(next(coordinates), next(coordinates)) for _ in vertices]
                )
                for vertices in fronts
            ]

        labelled = {}
        gluings = []
        for i, P in enumerate(polygons):
            Q, pants_gluings = _mirror(P, seams[i], None)
            labelled[i, 0] = P
            labelled[i, 1] = Q
            gluings.extend(
                (((i, a), e), ((i, b), f)) for (a, e), (b, f) in pants_gluings
            )

        for (i, s), e in half_cuffs.items():
            slot = chain[i][s]
            if slot[0] == "leg" and slot[1][0] == "fold":
                size = len(polygons[i].vertices())
                gluings.append((((i, 0), e), ((i, 1), -e % size)))
        for (i, a), (j, b) in cuff_ends.values():
            n_i, n_j = len(polygons[i].vertices()), len(polygons[j].vertices())
            gluings.append((((i, 0), a), ((j, 0), b)))
            gluings.append((((i, 1), -a % n_i), ((j, 1), -b % n_j)))

        surface = MutableOrientedHyperbolicSurface(polygons[0].parent())
        for label, polygon in labelled.items():
            surface.add_polygon(polygon, label=label)
        for x, y in gluings:
            surface.glue(x, y)
        surface.set_immutable()

        if not standard:
            return surface
        edges = []
        for (i, _), e in half_cuffs.items():
            size = len(polygons[i].vertices())
            edges.extend([((i, 0), e), ((i, 1), -e % size)])
        surface = _standard_surface(surface, sorted(edges))
        return surface if base_ring is None else surface.change_ring(base_ring)


hyperbolic_surfaces = HyperbolicSurfaceGenerators()


class FuchsianGroupGenerators:
    r"""
    Constructors of Fuchsian groups, as subgroups of the group of all the
    isometries of a hyperbolic plane, matrices of `\mathrm{O}(Q)^+` acting on
    its Klein model.

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
    def triangle(
        m1, m2, m3, base_ring=None, standard=False, orientation_preserving=True
    ):
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

        - ``base_ring`` -- a field or ``None`` (default: ``None``); with
          ``None``, the field of the group described below; otherwise an
          exact field into which that field converts along the embeddings in
          the reals, such as ``AA``, into which the group is converted

        - ``standard`` -- a boolean (default: ``False``); whether to carry
          the group to the plane of the form `x^2 + y^2 - z^2`

        - ``orientation_preserving`` -- a boolean (default: ``True``); with
          ``False``, the group generated by the reflections in the sides of
          the triangle instead

        OUTPUT:

        A :class:`~flatsurf.geometry.hyperbolic_isometry.HyperbolicIsometrySubgroup`
        generated by `g_1, g_2, g_3`, of the group of all the isometries of
        the plane (`\mathrm{O}(Q)^+`) of
        :meth:`HyperbolicPolygonsConstructor.triangle` with the angles
        `1/(2 m_1), 1/(2 m_3), 1/(2 m_2)`, over its field `k_0`, the
        invariant trace field of the group. The vertex of period `m_1` is at
        the origin and its edge to the vertex of period `m_3` on the
        positive `x`-axis. The matrices of the generators are reached by
        :meth:`~flatsurf.geometry.subgroup.FinitelyGeneratedSubgroup.matrix_group`.

        With ``standard``, the image of that group by the
        :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.standard_form_isometry`
        of its plane with ``algorithm="extension"``, over `k_0` or a
        quadratic extension of `k_0`.

        With ``orientation_preserving=False``, the group is the extended
        triangle group `W` generated by the reflections `r_0, r_1, r_2` in
        the sides of the triangle. It contains reflections, so it is an NEC
        group (non-Euclidean crystallographic), not a Fuchsian group, and its
        :meth:`~flatsurf.geometry.subgroup.FinitelyGeneratedSubgroup.orientation_preserving_subgroup`
        is the triangle group, with other generators than `g_1, g_2, g_3`.

        EXAMPLES:

        The triangle group `(3, 4, 4)`::

            sage: from flatsurf import fuchsian_groups
            sage: G = fuchsian_groups.triangle(3, 4, 4)
            sage: G
            Subgroup with 3 generators of Group of isometries of the hyperbolic plane over Rational Field with form x^2 + 3*y^2 - 6*z^2
            sage: G.gens()
            ([-1/2 -3/2    0]
             [ 1/2 -1/2    0]
             [   0    0    1],
             [-1/2 -9/2    6]
             [ 1/2 -3/2    2]
             [   0   -2    3],
             [-1/2 -3/2    0]
             [ 3/2 -3/2    4]
             [   1   -1    3])
            sage: G.base_ring()
            Rational Field
            sage: g1, g2, g3 = G.gens()
            sage: (g1^3).is_one(), (g2^4).is_one(), (g3^4).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True
            sage: G.matrix_group()
            Matrix group over Rational Field with 3 generators (
            ...
            )

        The triangle group `(2, 4, 6)`, given by Python integers::

            sage: G = fuchsian_groups.triangle(int(2), int(4), int(6))
            sage: G.base_ring()
            Rational Field
            sage: g1, g2, g3 = G.gens()
            sage: (g1^2).is_one(), (g2^4).is_one(), (g3^6).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True

        The triangle group `(2, 3, 7)` is defined over its cubic invariant
        trace field::

            sage: G = fuchsian_groups.triangle(2, 3, 7)
            sage: G.base_ring().polynomial()
            x^3 - 7*x - 7
            sage: g1, g2, g3 = G.gens()
            sage: (g1^2).is_one(), (g2^3).is_one(), (g3^7).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True
            sage: G.invariant_trace_field().polynomial()
            x^3 - 7*x - 7

        On the plane of the form `x^2 + y^2 - z^2`, it needs a field of
        degree 6::

            sage: G = fuchsian_groups.triangle(2, 3, 7, standard=True)
            sage: G.hyperbolic_plane()
            Hyperbolic Plane over Number Field in s with defining polynomial ... with s = ...
            sage: G.base_ring().degree()
            6
            sage: g1, g2, g3 = G.gens()
            sage: (g1^2).is_one(), (g2^3).is_one(), (g3^7).is_one()
            (True, True, True)
            sage: (g1 * g2 * g3).is_one()
            True

        The matrices can be converted into ``AA``::

            sage: G = fuchsian_groups.triangle(3, 4, 4, base_ring=AA)
            sage: G.base_ring()
            Algebraic Real Field
            sage: g1, g2, g3 = G.gens()
            sage: (g1 * g2 * g3).is_one()
            True

        With ``orientation_preserving=False``, the extended triangle group
        generated by the three reflections in the sides of the triangle::

            sage: W = fuchsian_groups.triangle(2, 3, 7, orientation_preserving=False)
            sage: W.ngens()
            3
            sage: W.is_orientable()
            False
            sage: all(g.is_reflection() for g in W.gens())
            True
            sage: W.invariant_trace_field().polynomial()
            x^3 - 7*x - 7

        On the plane of the form `x^2 + y^2 - z^2`, the reflections have
        determinant `-1`::

            sage: W = fuchsian_groups.triangle(3, 4, 4, standard=True, orientation_preserving=False)
            sage: W.base_ring().degree()
            2
            sage: [g.det() for g in W.gens()]
            [-1, -1, -1]

        The periods are integers at least 2, and the group is hyperbolic::

            sage: fuchsian_groups.triangle(1, 4, 4)
            Traceback (most recent call last):
            ...
            ValueError: the period 1 is not an integer at least 2
            sage: fuchsian_groups.triangle(2, 3, 6)
            Traceback (most recent call last):
            ...
            ValueError: the triangle group (2, 3, 6) is euclidean, not hyperbolic: 1/m1 + 1/m2 + 1/m3 must be less than 1; use crystallographic_groups.triangle
            sage: fuchsian_groups.triangle(2, 3, 5)
            Traceback (most recent call last):
            ...
            ValueError: the triangle group (2, 3, 5) is spherical, not hyperbolic: 1/m1 + 1/m2 + 1/m3 must be less than 1

        TESTS::

            sage: fuchsian_groups.triangle(2, 3, 7, base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field
            sage: fuchsian_groups.triangle(3, 4, 4, standard=True, base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates of the polygon do not lie in Rational Field
            sage: fuchsian_groups.triangle(5/2, 4, 4)
            Traceback (most recent call last):
            ...
            ValueError: the period 5/2 is not an integer at least 2
            sage: fuchsian_groups.triangle(3, 4, 4, base_ring=RR)
            Traceback (most recent call last):
            ...
            ValueError: the ring must be exact, not Real Field with 53 bits of precision

        The generators are the images by the holonomy of the pillowcase of
        the triangle of the counterclockwise loops around its vertices::

            sage: from flatsurf import hyperbolic_polygons, hyperbolic_surfaces
            sage: G = fuchsian_groups.triangle(3, 4, 4)
            sage: T = hyperbolic_polygons.triangle(1/6, 1/8, 1/8)
            sage: S = hyperbolic_surfaces.pillowcase(T)
            sage: phi = S.holonomy_representation()
            sage: loops = phi.domain().vertex_loops()
            sage: g1, g3, g2 = (phi(loops[S(0, i)]) for i in range(3))
            sage: [g.matrix() for g in G.gens()] == [g.matrix() for g in (g1, g2, g3)]
            True

        ALGORITHM:

        Let `r_0, r_1, r_2` be the reflections in the sides of the triangle,
        the side `i` going from the vertex `i` to the vertex `i + 1`. The
        generators are the words `r_2 r_0, r_1 r_2, r_2 r_1 r_0 r_2`. They
        are the images by the holonomy of the pillowcase of the triangle of
        the counterclockwise loops around its vertices 0, 2, 1, and they
        multiply to the identity in this order. With
        ``orientation_preserving=False``, the generators are `r_0, r_1, r_2`.

        The triangle lies in the plane of its invariant form, over the
        invariant trace field `k_0` of the group of Takeuchi, Arithmetic
        triangle groups (1977), Propositions 4 and 5, so that the reflections
        in its sides and the generators are over `k_0`. The plane of the
        form `x^2 + y^2 - z^2` needs a quadratic extension of `k_0` when the
        quaternion algebra of the group is a division algebra over `k_0`.

        """
        from flatsurf.geometry.hyperbolic_isometry import (
            HyperbolicIsometryGroup,
            reflection,
        )
        from flatsurf.geometry.signature import SurfaceSignature
        from flatsurf.geometry.subgroup import _triangle_group

        _check_exact(base_ring)
        geometry = SurfaceSignature.from_orbifold_orders(0, (m1, m2, m3)).geometry()
        if geometry != "hyperbolic":
            hint = (
                "; use crystallographic_groups.triangle"
                if geometry == "euclidean"
                else ""
            )
            raise ValueError(
                f"the triangle group ({m1}, {m2}, {m3}) is {geometry}, not hyperbolic: 1/m1 + 1/m2 + 1/m3 must be less than 1{hint}"
            )
        T = hyperbolic_polygons.triangle(*(QQ((1, 2 * ZZ(m))) for m in (m1, m3, m2)))

        ambient = T._transformation_group()
        V = T.vertices()
        reflections = [ambient(reflection(V[i], V[(i + 1) % 3])) for i in range(3)]
        G = _triangle_group(ambient, reflections, orientation_preserving)

        if standard:
            G = T.parent().standard_form_isometry(algorithm="extension")(G)

        if base_ring is not None:
            convert = _field_embedding(G.base_ring(), base_ring)
            form = G.quadratic_form()
            ambient = HyperbolicIsometryGroup(
                base_ring, form=[convert(form[i, i]) for i in range(3)]
            )
            G = ambient.subgroup([g.matrix().apply_map(convert) for g in G.gens()])

        return G


fuchsian_groups = FuchsianGroupGenerators()
