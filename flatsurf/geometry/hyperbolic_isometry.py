r"""
Orientation preserving isometries of the hyperbolic plane.

The orientation preserving isometries of the hyperbolic plane form the group
`\mathrm{SO}(2, 1)^+` of the 3×3 matrices `M` with `M^T J M = J` for `J =
\mathrm{diag}(1, 1, -1)`, determinant `1` and `M_{2,2} > 0`. Such a matrix
acts on the hyperboloid model, and on the Klein model by sending the point
`(x, y)` to `(X/Z, Y/Z)` where `(X, Y, Z) = M (x, y, 1)`.

The :class:`HyperbolicIsometryGroup` over a field `K` is the group of those
matrices with entries in `K`. The functions :func:`tangent_vector_isometry`
and :func:`segment_isometry` build such matrices in closed form.

EXAMPLES:

The isometry gluing the bottom side of a Lambert quadrilateral to its left
side is a quarter turn about the origin::

    sage: from flatsurf import HyperbolicPlane
    sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup, segment_isometry
    sage: H = HyperbolicPlane(QQ)
    sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
    sage: M = segment_isometry(P.edge_segment(0), -P.edge_segment(3))
    sage: M
    [ 0 -1  0]
    [ 1  0  0]
    [ 0  0  1]
    sage: G = HyperbolicIsometryGroup(QQ)
    sage: g = G(M)
    sage: g * vector((1/2, 0))
    (0, 1/2)
    sage: g(P)
    Hyperbolic polygon with vertices [I, (0, 1/2), (-1/2, 1/2), (-1/2, 0)]

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

import operator

from sage.categories.action import Action
from sage.groups.group import Group
from sage.misc.cachefunc import cached_method
from sage.structure.element import Matrix, MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation

from flatsurf.geometry.hyperbolic_polygon import _cosh2, _lift, _minkowski


def _minkowski_form(ring):
    r"""
    Return the matrix `J = \mathrm{diag}(1, 1, -1)` of the Minkowski form
    over ``ring``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import _minkowski_form
        sage: _minkowski_form(QQ)
        [ 1  0  0]
        [ 0  1  0]
        [ 0  0 -1]

    """
    from sage.matrix.constructor import diagonal_matrix

    return diagonal_matrix(ring, [1, 1, -1])


def _positive_sqrt(x):
    r"""
    Return the positive square root of the positive element ``x`` in its
    parent, or ``None`` if ``x`` is not a square there.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import _positive_sqrt
        sage: _positive_sqrt(QQ(9/4))
        3/2
        sage: _positive_sqrt(QQ(9/2)) is None
        True
        sage: _positive_sqrt(AA(9/2))
        2.121320343559643?

    """
    try:
        root = x.sqrt(extend=False)
    except ValueError:
        return None

    if root < 0:
        root = -root

    return root


def _not_defined_over(ring, square):
    r"""
    Return the error raised when an isometry needs the square root of
    ``square``, which is not in ``ring``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import _not_defined_over
        sage: _not_defined_over(QQ, 3/4)
        ValueError('isometry is not defined over Rational Field; adjoin a square root of 3/4')

    """
    return ValueError(
        f"isometry is not defined over {ring}; adjoin a square root of {square}"
    )


def _klein_vector(point):
    r"""
    Return the Klein coordinates of ``point`` as a vector.

    INPUT:

    - ``point`` -- a finite point of a hyperbolic plane, or a pair or a vector
      of Klein coordinates

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_isometry import _klein_vector
        sage: H = HyperbolicPlane(QQ)
        sage: _klein_vector(H(2*I))
        (0, 3/5)
        sage: _klein_vector((1/2, 0))
        (1/2, 0)
        sage: _klein_vector(H(0))
        Traceback (most recent call last):
        ...
        ValueError: 0 is not a finite point of the hyperbolic plane

    """
    from sage.modules.free_module_element import vector

    from flatsurf.geometry.hyperbolic import HyperbolicPoint

    if isinstance(point, HyperbolicPoint):
        if not point.is_finite():
            raise ValueError(f"{point} is not a finite point of the hyperbolic plane")
        return vector(point.parent().base_ring(), point.coordinates(model="klein"))

    point = vector(point)
    if len(point) != 2:
        raise ValueError("Klein coordinates must be a pair")

    return point


def _common_ring(vectors, ring):
    r"""
    Return ``ring``, or, when it is ``None``, the fraction field of the
    pushout of the base rings of ``vectors``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import _common_ring
        sage: _common_ring([vector((1, 0)), vector(AA, (0, 1))], None)
        Algebraic Real Field
        sage: _common_ring([vector((1, 0))], None)
        Rational Field
        sage: _common_ring([vector((1, 0))], AA)
        Algebraic Real Field

    """
    if ring is not None:
        return ring

    from sage.categories.pushout import pushout

    ring = vectors[0].base_ring()
    for v in vectors[1:]:
        ring = pushout(ring, v.base_ring())

    return ring.fraction_field()


def tangent_vector_isometry(u, v, ring=None):
    r"""
    Return the matrix of the orientation preserving isometry of the
    hyperbolic plane that sends the tangent vector ``u`` to the tangent
    vector ``v``.

    INPUT:

    - ``u``, ``v`` -- tangent vectors, each a pair ``(point, direction)``
      with ``point`` a finite point of the hyperbolic plane, or a pair or a
      vector of Klein coordinates, and ``direction`` a nonzero vector of the
      Klein chart at ``point``, taken up to positive scaling

    - ``ring`` -- a field or ``None`` (default: ``None``); the field of the
      result, by default the fraction field of the pushout of the rings of
      the coordinates

    OUTPUT:

    An immutable 3×3 matrix of `\mathrm{SO}(2, 1)^+` over ``ring``. The group
    `\mathrm{SO}(2, 1)^+` acts simply transitively on the unit tangent
    bundle; the Klein chart preserves the orientation, so that a ray in the
    chart is a unit tangent vector: the matrix is unique.

    ALGORITHM:

    The point `p` lifts to `P = (p, 1)` on the cone of timelike vectors, and
    the direction `d` to the tangent vector `T = D - \langle P, D\rangle /
    \langle P, P\rangle P` of the hyperboloid, the Minkowski projection of
    `D = (d, 0)` onto `P^\perp`. The vector `N = J (P \times T)` satisfies
    `\langle N, X\rangle = \det(P, T, X)`, so that `(P, T, N)` is an
    orthogonal basis with `\langle N, N\rangle = \det(P, T, N) = -\langle P,
    P\rangle \langle T, T\rangle`. With the same construction at the image,
    the isometry is

    .. MATH::

        M = (\lambda P', \mu T', \lambda\mu N') (P, T, N)^{-1}

    with `\lambda^2 = \langle P, P\rangle / \langle P', P'\rangle` and
    `\mu^2 = \langle T, T\rangle / \langle T', T'\rangle`, `\lambda, \mu >
    0`; its determinant is `\lambda^2 \mu^2 \langle N', N'\rangle / \langle
    N, N\rangle = 1`. The inverse of the orthogonal basis is
    `\mathrm{diag}(1/\langle P, P\rangle, 1/\langle T, T\rangle, 1/\langle
    N, N\rangle) (P, T, N)^T J`. Since `\lambda` is the last entry of `M P`
    and `\mu T' = M T`, the matrix is defined over ``ring`` if and only if
    `\lambda` and `\mu` are.

    EXAMPLES:

    A quarter turn about the origin; the directions are taken up to positive
    scaling::

        sage: from flatsurf.geometry.hyperbolic_isometry import tangent_vector_isometry
        sage: M = tangent_vector_isometry(((0, 0), (1, 0)), ((0, 0), (0, 1)))
        sage: M
        [ 0 -1  0]
        [ 1  0  0]
        [ 0  0  1]
        sage: tangent_vector_isometry(((0, 0), (2, 0)), ((0, 0), (0, 1))) == M
        True
        sage: J = diagonal_matrix([1, 1, -1])
        sage: M.transpose() * J * M == J
        True

    The boost along the horizontal axis that sends `(1/2, 0)` to the origin
    needs `\sqrt{3}`::

        sage: tangent_vector_isometry(((1/2, 0), (1, 0)), ((0, 0), (1, 0)))
        Traceback (most recent call last):
        ...
        ValueError: isometry is not defined over Rational Field; adjoin a square root of 3/4
        sage: M = tangent_vector_isometry(((1/2, 0), (1, 0)), ((0, 0), (1, 0)), ring=AA)
        sage: M * vector(AA, (1/2, 0, 1)) == vector(AA, (0, 0, sqrt(AA(3/4))))
        True

    Points of the hyperbolic plane are accepted::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane(QQ)
        sage: tangent_vector_isometry((H(I), (1, 0)), (H(I), (-1, 0)))
        [-1  0  0]
        [ 0 -1  0]
        [ 0  0  1]

    TESTS::

        sage: tangent_vector_isometry(((0, 0), (0, 0)), ((0, 0), (0, 1)))
        Traceback (most recent call last):
        ...
        ValueError: direction must be nonzero
        sage: tangent_vector_isometry(((1, 0), (0, 1)), ((0, 0), (0, 1)))
        Traceback (most recent call last):
        ...
        ValueError: (1, 0) is not a finite point of the hyperbolic plane

    """
    from sage.matrix.constructor import diagonal_matrix, matrix
    from sage.modules.free_module_element import vector

    p, d = u
    q, e = v

    p = _klein_vector(p)
    q = _klein_vector(q)
    d = vector(d)
    e = vector(e)

    ring = _common_ring([p, q, d, e], ring)

    p, q, d, e = (w.change_ring(ring) for w in (p, q, d, e))

    for point in (p, q):
        if not point[0] ** 2 + point[1] ** 2 < 1:
            raise ValueError(f"{point} is not a finite point of the hyperbolic plane")

    for direction in (d, e):
        if len(direction) != 2:
            raise ValueError("direction must be a vector of dimension two")
        if direction.is_zero():
            raise ValueError("direction must be nonzero")

    def frame(point, direction):
        P = vector(ring, _lift(point))
        D = vector(ring, (direction[0], direction[1], 0))
        T = D - (_minkowski(P, D) / _minkowski(P, P)) * P
        cross = P.cross_product(T)
        N = vector(ring, (cross[0], cross[1], -cross[2]))
        return P, T, N

    P, T, N = frame(p, d)
    P_, T_, N_ = frame(q, e)

    lambda2 = _minkowski(P, P) / _minkowski(P_, P_)
    lam = _positive_sqrt(lambda2)
    if lam is None:
        raise _not_defined_over(ring, lambda2)

    mu2 = _minkowski(T, T) / _minkowski(T_, T_)
    mu = _positive_sqrt(mu2)
    if mu is None:
        raise _not_defined_over(ring, mu2)

    source = matrix(ring, [P, T, N]).transpose()
    target = matrix(ring, [lam * P_, mu * T_, lam * mu * N_]).transpose()

    inverse = (
        diagonal_matrix(
            ring,
            [1 / _minkowski(P, P), 1 / _minkowski(T, T), 1 / _minkowski(N, N)],
        )
        * source.transpose()
        * _minkowski_form(ring)
    )

    M = target * inverse
    M.set_immutable()
    return M


def _segment_endpoints(segment):
    r"""
    Return the start and the end of ``segment``.

    INPUT:

    - ``segment`` -- an oriented segment of the hyperbolic plane, or a pair of
      points

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_isometry import _segment_endpoints
        sage: H = HyperbolicPlane(QQ)
        sage: _segment_endpoints(H(I).segment(H(2*I)))
        (I, 2*I)
        sage: _segment_endpoints(((0, 0), (1/2, 0)))
        ((0, 0), (1/2, 0))

    """
    from flatsurf.geometry.hyperbolic import HyperbolicOrientedSegment

    if isinstance(segment, HyperbolicOrientedSegment):
        return segment.start(), segment.end()

    start, end = segment
    return start, end


def segment_isometry(source, target, ring=None):
    r"""
    Return the matrix of the orientation preserving isometry of the
    hyperbolic plane that sends the start of ``source`` to the start of
    ``target`` and the end of ``source`` to the end of ``target``.

    INPUT:

    - ``source``, ``target`` -- oriented segments of the hyperbolic plane with
      finite endpoints, or pairs of finite points or of Klein coordinates

    - ``ring`` -- a field or ``None`` (default: ``None``); as in
      :func:`tangent_vector_isometry`

    OUTPUT:

    An immutable 3×3 matrix of `\mathrm{SO}(2, 1)^+` over ``ring``.

    ALGORITHM:

    The two segments must have the same length, compared through
    `\cosh^2 d(p, q) = \langle P, Q\rangle^2 / (\langle P, P\rangle \langle
    Q, Q\rangle)`. The isometry is then the one of
    :func:`tangent_vector_isometry` sending the tangent vector at the start
    of ``source`` pointing to its end to the one of ``target``: geodesics are
    straight lines in the Klein model.

    EXAMPLES:

    The side pairings of a Lambert quadrilateral; the pairing of the side `i`
    to the side `j` sends the start of `i` to the end of `j`::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_isometry import segment_isometry
        sage: H = HyperbolicPlane(QQ)
        sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
        sage: J = diagonal_matrix([1, 1, -1])
        sage: M = segment_isometry(P.edge_segment(0), -P.edge_segment(3))
        sage: M
        [ 0 -1  0]
        [ 1  0  0]
        [ 0  0  1]
        sage: M.transpose() * J * M == J
        True
        sage: M = segment_isometry(P.edge_segment(1), -P.edge_segment(2))
        sage: M
        [   0    1    0]
        [-5/3    0  4/3]
        [-4/3    0  5/3]
        sage: M.trace()
        5/3
        sage: M.transpose() * J * M == J
        True
        sage: M * vector((1/2, 0, 1)), M * vector((1/2, 1/2, 1))
        ((0, 1/2, 1), (1/2, 1/2, 1))

    The side pairings of a non-convex quadrilateral::

        sage: P = H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
        sage: M = segment_isometry(P.edge_segment(0), -P.edge_segment(3))
        sage: M
        [ 13/19 -16/19  -8/19]
        [ 16/19  11/19  -4/19]
        [ -8/19   4/19  21/19]
        sage: M.transpose() * J * M == J
        True
        sage: M * vector((0, -1/2, 1)), M * vector((1/2, 1/2, 1))
        ((0, -1/2, 1), (-1/2, 1/2, 1))
        sage: M = segment_isometry(P.edge_segment(1), -P.edge_segment(2))
        sage: M
        [ 0 -1  0]
        [ 1  0  0]
        [ 0  0  1]
        sage: M * vector((1/2, 1/2, 1))
        (-1/2, 1/2, 1)

    A side glued to itself with the reverse orientation; the half turn about
    its midpoint needs `\sqrt{3}`::

        sage: segment_isometry(((0, 0), (1/2, 0)), ((1/2, 0), (0, 0)))
        Traceback (most recent call last):
        ...
        ValueError: isometry is not defined over Rational Field; adjoin a square root of 4/3
        sage: L.<r> = NumberField(x^2 - 3, embedding=AA(3).sqrt())
        sage: M = segment_isometry(((0, 0), (1/2, 0)), ((1/2, 0), (0, 0)), ring=L)
        sage: M
        [-2/3*r      0  1/3*r]
        [     0     -1      0]
        [-1/3*r      0  2/3*r]
        sage: M.trace()
        -1
        sage: M.transpose() * J * M == J
        True
        sage: M * vector((0, 0, 1)), M * vector((1/2, 0, 1))
        ((1/3*r, 0, 2/3*r), (0, 0, 1/2*r))

    TESTS::

        sage: segment_isometry(((0, 0), (1/2, 0)), ((0, 0), (0, 1/3)))
        Traceback (most recent call last):
        ...
        ValueError: segments do not have the same length
        sage: segment_isometry(((0, 0), (0, 0)), ((0, 0), (0, 0)))
        Traceback (most recent call last):
        ...
        ValueError: direction must be nonzero
        sage: segment_isometry(((1, 0), (0, 0)), ((1, 0), (0, 0)))
        Traceback (most recent call last):
        ...
        ValueError: (1, 0) is not a finite point of the hyperbolic plane

    """
    a, b = (_klein_vector(point) for point in _segment_endpoints(source))
    c, d = (_klein_vector(point) for point in _segment_endpoints(target))

    ring = _common_ring([a, b, c, d], ring)
    a, b, c, d = (w.change_ring(ring) for w in (a, b, c, d))

    # The Minkowski norm of the lift of a point of the boundary vanishes, so
    # the lengths can only be compared for points inside the disk.
    for point in (a, b, c, d):
        if not point[0] ** 2 + point[1] ** 2 < 1:
            raise ValueError(f"{point} is not a finite point of the hyperbolic plane")

    if _cosh2(a, b) != _cosh2(c, d):
        raise ValueError("segments do not have the same length")

    return tangent_vector_isometry((a, b - a), (c, d - c), ring=ring)


class HyperbolicIsometry(MultiplicativeGroupElement):
    r"""
    An orientation preserving isometry of the hyperbolic plane, given by a
    matrix of `\mathrm{SO}(2, 1)^+`.

    INPUT:

    - ``parent`` -- a :class:`HyperbolicIsometryGroup`

    - ``matrix`` -- a 3×3 matrix of `\mathrm{SO}(2, 1)^+` over the base ring
      of ``parent``, immutable; not checked: isometries are built by calling
      the group, which checks its input

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
        sage: G = HyperbolicIsometryGroup(QQ)
        sage: g = G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
        sage: g
        [ 0 -1  0]
        [ 1  0  0]
        [ 0  0  1]

    TESTS::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometry
        sage: isinstance(g, HyperbolicIsometry)
        True
        sage: TestSuite(g).run()

    """

    def __init__(self, parent, matrix):
        MultiplicativeGroupElement.__init__(self, parent)
        self._matrix = matrix

    def __reduce__(self):
        r"""
        Return a picklable representation of this isometry.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: g = HyperbolicIsometryGroup(QQ).an_element()
            sage: loads(dumps(g)) == g
            True

        """
        return self.parent(), (self._matrix,)

    def matrix(self):
        r"""
        Return the matrix of `\mathrm{SO}(2, 1)^+` of this isometry.

        The matrix is immutable; it acts on the Klein model as in
        :meth:`~flatsurf.geometry.hyperbolic.HyperbolicConvexSet.apply_isometry`
        with ``model="klein"``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: g = HyperbolicIsometryGroup(QQ).an_element()
            sage: g.matrix()
            [ 7 -4  8]
            [ 4 -1  4]
            [ 8 -4  9]
            sage: g.matrix().is_immutable()
            True

        """
        return self._matrix

    def trace(self):
        r"""
        Return the trace of the matrix of this isometry.

        The trace is `1 + 2\cos\theta` for a rotation of angle `\theta`, `3`
        for a parabolic isometry and `1 + 2\cosh\ell` for a hyperbolic
        isometry of translation length `\ell`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).an_element().trace()
            15

        """
        return self._matrix.trace()

    def _repr_(self):
        r"""
        Return a printable representation of this isometry, its matrix.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).one()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        return repr(self._matrix)

    def _mul_(self, other):
        r"""
        Return the composition of this isometry with ``other``, ``other``
        applied first.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: g = G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
            sage: g * g
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

        """
        M = self._matrix * other._matrix
        M.set_immutable()
        return self.parent().element_class(self.parent(), M)

    def __invert__(self):
        r"""
        Return the inverse of this isometry, `J M^T J`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: g = G.an_element()
            sage: ~g
            [ 7  4 -8]
            [-4 -1  4]
            [-8 -4  9]
            sage: ~g * g == G.one()
            True

        """
        J = _minkowski_form(self.parent().base_ring())
        M = J * self._matrix.transpose() * J
        M.set_immutable()
        return self.parent().element_class(self.parent(), M)

    def __eq__(self, other):
        r"""
        Return whether ``other`` is an isometry of the same group with the
        same matrix.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: G.one() == G.one()
            True
            sage: G.one() == G.an_element()
            False
            sage: G.one() == HyperbolicIsometryGroup(AA).one()
            False

        """
        if not isinstance(other, HyperbolicIsometry):
            return False

        return self.parent() is other.parent() and self._matrix == other._matrix

    def __ne__(self, other):
        r"""
        Return whether this isometry differs from ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: G.one() != G.an_element()
            True

        """
        return not self == other

    def __hash__(self):
        r"""
        Return a hash value for this isometry, compatible with ``==``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: hash(G.an_element()) == hash(G.an_element())
            True

        """
        return hash(self._matrix)

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this isometry.

        INPUT:

        - ``x`` -- a vector or a pair of Klein coordinates, a subset of the
          hyperbolic plane such as a point, or a
          :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`

        OUTPUT:

        For Klein coordinates, an immutable vector as for ``self * x``;
        otherwise a subset of the same hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: H = HyperbolicPlane(QQ)
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: g = G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
            sage: g((1/2, 0))
            (0, 1/2)
            sage: g(H(2*I))
            -3/5 + 4/5*I
            sage: g(H(I).segment(H(2*I)))
            {-(x^2 + y^2) + 1 = 0} ∩ {x ≤ 0} ∩ {3*(x^2 + y^2) + 10*x + 3 ≥ 0}
            sage: g(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            Hyperbolic polygon with vertices [I, (0, 1/2), (-1/2, 0)]

        """
        from flatsurf.geometry.hyperbolic import HyperbolicConvexSet
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

        if isinstance(x, HyperbolicPolygon):
            return x.apply_isometry(self)

        if isinstance(x, HyperbolicConvexSet):
            return x.apply_isometry(self._matrix, model="klein")

        from sage.modules.free_module_element import FreeModuleElement, vector

        if not isinstance(x, FreeModuleElement):
            x = vector(x)

        return self * x

    def psl2_matrix(self, ring=None):
        r"""
        Return a matrix of `\mathrm{SL}(2, R)` that acts on the upper half
        plane as this isometry.

        INPUT:

        - ``ring`` -- a ring or ``None`` (default: ``None``); the ring of the
          result, by default the base ring

        OUTPUT:

        An immutable 2×2 matrix `A` of determinant `1` whose first nonzero
        entry is positive, such that
        :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane._isometry_gl2_to_sim12`
        sends it to the matrix of this isometry; so `A` acts as this
        isometry in
        :meth:`~flatsurf.geometry.hyperbolic.HyperbolicConvexSet.apply_isometry`
        with ``model="half_plane"``.

        ALGORITHM:

        The matrix of this isometry is quadratic in the entries `a, b, c, d`
        of `A`: its entries give `a^2, b^2, c^2, d^2` and the products of
        two of them. As `ad - bc = 1`, either `a` or `b` is nonzero. If
        `s = a^2` is nonzero, then `A = (a^2, ab; ac, ad)/a` with `a =
        \sqrt{s}`; otherwise `A = (ab, b^2; bc, bd)/b` with `b = \sqrt{s}`
        for `s = b^2`. So `A` is defined over ``ring`` if and only if `s` is
        a square there. Note that `\mathrm{tr}(M) = \mathrm{tr}(A)^2 - 1`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: H = HyperbolicPlane(QQ)
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: A = matrix([[2, 3], [1, 2]])
            sage: g = G(H._isometry_gl2_to_sim12(A))
            sage: g
            [ 7 -4  8]
            [ 4 -1  4]
            [ 8 -4  9]
            sage: g.psl2_matrix()
            [2 3]
            [1 2]

        A boost of translation length `\ell` with `\cosh\ell = 5/4` needs
        `\sqrt{2}`::

            sage: g = G(matrix([[5/4, 0, 3/4], [0, 1, 0], [3/4, 0, 5/4]]))
            sage: g.psl2_matrix()
            Traceback (most recent call last):
            ...
            ValueError: isometry is not defined over Rational Field; adjoin a square root of 9/8
            sage: A = g.psl2_matrix(AA)
            sage: A
            [ 1.060660171779822? 0.3535533905932738?]
            [0.3535533905932738?  1.060660171779822?]
            sage: H._isometry_gl2_to_sim12(A) == g.matrix()
            True
            sage: A.trace()^2 - 1 == g.trace()
            True

        A rotation of angle `\theta` about `i` is a rotation of angle
        `\theta/2` in `\mathrm{SL}(2, \RR)`; the quarter turn needs
        `\sqrt{2}`::

            sage: g = G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
            sage: g.psl2_matrix()
            Traceback (most recent call last):
            ...
            ValueError: isometry is not defined over Rational Field; adjoin a square root of 1/2

        TESTS:

        The first entry may vanish::

            sage: g = G(matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]))
            sage: A = g.psl2_matrix()
            sage: A
            [ 0  1]
            [-1  0]
            sage: H._isometry_gl2_to_sim12(A) == g.matrix()
            True

        The products of the entries may not be in ``ring`` at all::

            sage: L.<r> = NumberField(x^2 - 3, embedding=AA(3).sqrt())
            sage: g = HyperbolicIsometryGroup(L)(matrix([[2, 0, r], [0, 1, 0], [r, 0, 2]]))
            sage: g.psl2_matrix(QQ)
            Traceback (most recent call last):
            ...
            ValueError: isometry is not defined over Rational Field

        """
        from sage.matrix.constructor import matrix

        K = self.parent().base_ring()
        if ring is None:
            ring = K

        m00, m01, m02, m10, m11, m12, m20, m21, m22 = self._matrix.list()
        two = K(2)

        a2 = (m11 + m12 + m21 + m22) / two
        if a2:
            square = a2
            entries = [a2, (m10 + m20) / two, (m01 + m02) / two, (m00 + 1) / two]
        else:
            square = (m12 + m22 - m11 - m21) / two
            entries = [K.zero(), square, (m00 - 1) / two, (m02 - m01) / two]

        try:
            square_ = ring(square)
            entries = [ring(entry) for entry in entries]
        except (TypeError, ValueError):
            raise ValueError(f"isometry is not defined over {ring}") from None

        root = _positive_sqrt(square_)
        if root is None:
            raise _not_defined_over(ring, square)

        A = matrix(ring, 2, [entry / root for entry in entries])
        A.set_immutable()
        return A


class HyperbolicIsometryAction(Action):
    r"""
    The action of a :class:`HyperbolicIsometryGroup` on the Klein
    coordinates of the points of the hyperbolic plane, by multiplication
    from the left.

    The image of `(x, y)` is `(X/Z, Y/Z)` for `(X, Y, Z) = M (x, y, 1)`; it
    lives in the free module of rank two over the pushout of the base ring of
    the group and the base ring of the vector.

    INPUT:

    - ``G`` -- a :class:`HyperbolicIsometryGroup`

    - ``M`` -- a free module of rank two

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryAction, HyperbolicIsometryGroup
        sage: G = HyperbolicIsometryGroup(QQ)
        sage: A = HyperbolicIsometryAction(G, ZZ^2)
        sage: A
        Left action by Group of orientation preserving isometries of the hyperbolic plane over Rational Field on Ambient free module of rank 2 over the principal ideal domain Integer Ring
        sage: A.codomain()
        Vector space of dimension 2 over Rational Field

    """

    def __init__(self, G, M):
        from sage.categories.pushout import pushout
        from sage.modules.free_module import FreeModule

        self._codomain = FreeModule(
            pushout(G.base_ring(), M.base_ring()).fraction_field(), 2
        )
        super().__init__(G, M, True, operator.mul)

    def codomain(self):
        r"""
        Return the parent of the images of this action.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryAction, HyperbolicIsometryGroup
            sage: HyperbolicIsometryAction(HyperbolicIsometryGroup(QQ), AA^2).codomain()
            Vector space of dimension 2 over Algebraic Real Field

        """
        return self._codomain

    def _act_(self, g, v):
        r"""
        Return the image of the vector ``v`` under the isometry ``g``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: g = G.an_element()
            sage: w = g * vector(AA, (0, 0))
            sage: w
            (8/9, 4/9)
            sage: w.parent()
            Vector space of dimension 2 over Algebraic Real Field
            sage: w.is_immutable()
            True

        TESTS:

        A point outside the disk may be sent to infinity in the Klein chart::

            sage: g = G(matrix([[5/4, 0, 3/4], [0, 1, 0], [3/4, 0, 5/4]]))
            sage: g * vector((-5/3, 0))
            Traceback (most recent call last):
            ...
            ValueError: (-5/3, 0) is sent outside of the Klein chart

        """
        from sage.modules.free_module_element import vector

        R = self._codomain.base_ring()
        M = g.matrix()
        x, y = v
        X = M[0, 0] * x + M[0, 1] * y + M[0, 2]
        Y = M[1, 0] * x + M[1, 1] * y + M[1, 2]
        Z = M[2, 0] * x + M[2, 1] * y + M[2, 2]
        if not Z:
            raise ValueError(f"{v} is sent outside of the Klein chart")
        return vector(R, [X / Z, Y / Z], immutable=True)


class HyperbolicIsometryPlaneAction(Action):
    r"""
    The action of a :class:`HyperbolicIsometryGroup` on the points of a
    hyperbolic plane, by multiplication from the left.

    The other subsets of the hyperbolic plane (segments, polygons, ...) are
    parents, not elements of the plane, so that ``*`` does not reach this
    action for them; their image under an isometry ``g`` is ``g(x)``.

    INPUT:

    - ``G`` -- a :class:`HyperbolicIsometryGroup`

    - ``H`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane` whose
      base ring has a coercion from the base ring of ``G``

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup, HyperbolicIsometryPlaneAction
        sage: G = HyperbolicIsometryGroup(QQ)
        sage: H = HyperbolicPlane(QQ)
        sage: HyperbolicIsometryPlaneAction(G, H)
        Left action by Group of orientation preserving isometries of the hyperbolic plane over Rational Field on Hyperbolic Plane over Rational Field

    """

    def __init__(self, G, H):
        super().__init__(G, H, True, operator.mul)

    def _act_(self, g, x):
        r"""
        Return the image of the point ``x`` of the hyperbolic plane under the
        isometry ``g``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: H = HyperbolicPlane(QQ)
            sage: g = G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
            sage: g * H(2*I)
            -3/5 + 4/5*I

        """
        return g(x)


class HyperbolicIsometryGroup(UniqueRepresentation, Group):
    r"""
    The group of orientation preserving isometries of the hyperbolic plane
    whose matrices in `\mathrm{SO}(2, 1)^+` have entries in ``base_ring``.

    INPUT:

    - ``base_ring`` -- an ordered field

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
        sage: G = HyperbolicIsometryGroup(QQ)
        sage: G
        Group of orientation preserving isometries of the hyperbolic plane over Rational Field

    TESTS::

        sage: TestSuite(G).run()
        sage: TestSuite(HyperbolicIsometryGroup(AA)).run()

    """

    Element = HyperbolicIsometry

    def __init__(self, base_ring):
        from sage.categories.groups import Groups

        self._base_ring = base_ring
        Group.__init__(self, category=Groups().Infinite())

    def base_ring(self):
        r"""
        Return the field of the entries of the matrices of this group.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).base_ring()
            Rational Field

        """
        return self._base_ring

    def _repr_(self):
        r"""
        Return a printable representation of this group.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(AA)
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field

        """
        return f"Group of orientation preserving isometries of the hyperbolic plane over {self._base_ring}"

    @cached_method
    def one(self):
        r"""
        Return the identity of this group.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).one().is_one()
            True

        """
        from sage.matrix.matrix_space import MatrixSpace

        M = MatrixSpace(self._base_ring, 3).one()
        return self.element_class(self, M)

    def _an_element_(self):
        r"""
        Return an isometry of this group, the image of the matrix `(2, 3; 1,
        2)` acting on the upper half plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).an_element()
            [ 7 -4  8]
            [ 4 -1  4]
            [ 8 -4  9]

        """
        return self([[7, -4, 8], [4, -1, 4], [8, -4, 9]])

    def some_elements(self):
        r"""
        Return some isometries of this group: the identity, an elliptic, a
        parabolic and two hyperbolic isometries.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: [g.trace() for g in HyperbolicIsometryGroup(QQ).some_elements()]
            [3, 1, 3, 15, 7/2]

        """
        from sage.all import QQ

        return [
            self.one(),
            self([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            self(
                [[1, -1, 1], [1, QQ((1, 2)), QQ((1, 2))], [1, -QQ((1, 2)), QQ((3, 2))]]
            ),
            self.an_element(),
            self([[QQ((5, 4)), 0, QQ((3, 4))], [0, 1, 0], [QQ((3, 4)), 0, QQ((5, 4))]]),
        ]

    def _coerce_map_from_(self, other):
        r"""
        Return whether there is a coercion from ``other`` to this group: the
        groups of isometries over a field that coerces into the base ring.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: HyperbolicIsometryGroup(AA).has_coerce_map_from(G)
            True
            sage: G.has_coerce_map_from(HyperbolicIsometryGroup(AA))
            False
            sage: (HyperbolicIsometryGroup(AA).one() * G.an_element()).parent()
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field

        """
        if isinstance(other, HyperbolicIsometryGroup):
            return self._base_ring.has_coerce_map_from(other.base_ring())

        return False

    def _element_constructor_(self, x):
        r"""
        Return the isometry given by ``x``.

        INPUT:

        - ``x`` -- an isometry, a 3×3 matrix of `\mathrm{SO}(2, 1)^+`, or a
          2×2 matrix of positive determinant acting on the upper half plane
          as in
          :meth:`~flatsurf.geometry.hyperbolic.HyperbolicConvexSet.apply_isometry`

        ALGORITHM:

        A 2×2 matrix `A` goes through
        :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane._isometry_gl2_to_sim12`,
        which is quadratic in the entries of `A`; normalizing `A` to
        determinant `1` divides its image by `\det A`, so no root is needed.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: G(matrix([[2, 0], [0, 1]]))
            [   1    0    0]
            [   0  5/4  3/4]
            [   0  3/4  5/4]
            sage: G(matrix([[1, 0], [0, -1]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an orientation preserving isometry
            sage: G(matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(2, 1)+
            sage: G(matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(2, 1)+
            sage: G(matrix([[1, 2, 0], [0, 1, 0], [0, 0, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(2, 1)+
            sage: G(matrix([[1, 0], [0, 1], [0, 0]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix must be 2×2 or 3×3

        """
        from sage.matrix.constructor import matrix

        K = self._base_ring

        if isinstance(x, HyperbolicIsometry):
            M = x.matrix().change_ring(K)
            M.set_immutable()
            return self.element_class(self, M)

        if not isinstance(x, Matrix):
            x = matrix(K, x)

        if x.dimensions() == (2, 2):
            det = x.determinant()
            if not det > 0:
                raise ValueError(
                    "matrix does not define an orientation preserving isometry"
                )

            from flatsurf.geometry.hyperbolic import HyperbolicPlane

            M = HyperbolicPlane(K)._isometry_gl2_to_sim12(x.change_ring(K)) / K(det)
        elif x.dimensions() == (3, 3):
            M = x.change_ring(K)
            J = _minkowski_form(K)
            if M.transpose() * J * M != J or M.determinant() != 1 or not M[2, 2] > 0:
                raise ValueError("matrix is not in SO(2, 1)+")
        else:
            raise ValueError("matrix must be 2×2 or 3×3")

        M = M.change_ring(K)
        M.set_immutable()
        return self.element_class(self, M)

    def _get_action_(self, S, op, self_on_left):
        r"""
        Return the action of this group on ``S`` with respect to ``op``.

        This group acts by multiplication from the left on the free modules
        of rank two whose base ring has a common parent with the base ring of
        this group, see :class:`HyperbolicIsometryAction`, and on the
        hyperbolic planes over a field into which the base ring coerces, see
        :class:`HyperbolicIsometryPlaneAction`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: g = G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
            sage: g * vector((1/2, 0))
            (0, 1/2)
            sage: g * HyperbolicPlane(AA)(2*I)
            -3/5 + 4/5*I

        Only multiplication from the left is an action::

            sage: vector((1, 0)) * g
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Ambient free module of rank 2 over the principal ideal domain Integer Ring' and 'Group of orientation preserving isometries of the hyperbolic plane over Rational Field'

        """
        if self_on_left and op is operator.mul:
            from sage.modules.free_module import FreeModule_generic

            from flatsurf.geometry.hyperbolic import HyperbolicPlane

            if isinstance(S, FreeModule_generic) and S.degree() == 2:
                from sage.categories.pushout import pushout
                from sage.structure.coerce_exceptions import CoercionException

                try:
                    pushout(self._base_ring, S.base_ring())
                except (CoercionException, TypeError):
                    return None

                return HyperbolicIsometryAction(self, S)

            if isinstance(S, HyperbolicPlane):
                if S.base_ring().has_coerce_map_from(self._base_ring):
                    return HyperbolicIsometryPlaneAction(self, S)
                return None

        return super()._get_action_(S, op, self_on_left)
