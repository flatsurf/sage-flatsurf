r"""
Orientation preserving isometries of the hyperbolic plane.

The orientation preserving isometries of the hyperbolic plane form the group
`\mathrm{SO}(2, 1)^+` of the 3×3 matrices `M` with `M^T J M = J` for `J =
\mathrm{diag}(1, 1, -1)`, determinant `1` and `M_{2,2} > 0`. Such a matrix
acts on the hyperboloid model, and on the Klein model by sending the point
`(x, y)` to `(X/Z, Y/Z)` where `(X, Y, Z) = M (x, y, 1)`.

The :class:`HyperbolicIsometryGroup` over a field `K` is the group of those
matrices with entries in `K`. For a hyperbolic plane with another quadratic
form `Q = \mathrm{diag}(1, b, -c)`, whose Klein model is the disk `x^2 + b
y^2 < c`, the group is `\mathrm{SO}(Q)^+`, the matrices with `M^T Q M = Q`,
determinant `1` and `M_{2,2} > 0`. The functions :func:`tangent_vector_isometry`
and :func:`segment_isometry` build such matrices in closed form, and
:func:`segment_midpoint` gives the midpoint of a segment.

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

The group of the form `x^2 + 2y^2 - 3z^2` contains the boost that sends `(-1,
0)` to `(1, 0)`::

    sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
    sage: G3
    Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2
    sage: g = G3(matrix([[2, 0, 3], [0, 1, 0], [1, 0, 2]]))
    sage: g * vector((-1, 0))
    (1, 0)

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
from sage.categories.pushout import ConstructionFunctor
from sage.groups.group import Group
from sage.misc.cachefunc import cached_method
from sage.structure.element import Matrix, MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation

from flatsurf.geometry.hyperbolic_polygon import _cosh2, _form_repr, _lift, _minkowski


def _minkowski_form(ring, form=None):
    r"""
    Return the immutable matrix of the quadratic form over ``ring``.

    INPUT:

    - ``ring`` -- a ring

    - ``form`` -- ``None``, a 3×3 diagonal matrix or the list of its three
      diagonal coefficients (default: ``None``); the diagonal must be `(a, b,
      -c)` with `a, b, c` positive. ``None`` is `\mathrm{diag}(1, 1, -1)`.

    OUTPUT:

    The matrix of the form divided by its entry `(0, 0)`: a positive scaling
    changes neither the cone of the form nor its isometries.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import _minkowski_form
        sage: _minkowski_form(QQ)
        [ 1  0  0]
        [ 0  1  0]
        [ 0  0 -1]
        sage: _minkowski_form(QQ, [2, 4, -6])
        [ 1  0  0]
        [ 0  2  0]
        [ 0  0 -3]
        sage: _minkowski_form(QQ, diagonal_matrix([2, 4, -6])) == _minkowski_form(QQ, [1, 2, -3])
        True

    TESTS::

        sage: _minkowski_form(QQ, [1, 1, 1])
        Traceback (most recent call last):
        ...
        ValueError: the form must be diag(a, b, -c) with a, b, c positive
        sage: _minkowski_form(QQ, matrix([[1, 1, 0], [1, 2, 0], [0, 0, -3]]))
        Traceback (most recent call last):
        ...
        ValueError: the form must be diagonal

    """
    from sage.matrix.constructor import diagonal_matrix

    if form is None:
        diagonal = [1, 1, -1]
    elif isinstance(form, (list, tuple)):
        diagonal = list(form)
    else:
        if not form.is_diagonal():
            raise ValueError("the form must be diagonal")
        diagonal = [form[i, i] for i in range(3)]

    if len(diagonal) != 3 or not (
        diagonal[0] > 0 and diagonal[1] > 0 and diagonal[2] < 0
    ):
        raise ValueError("the form must be diag(a, b, -c) with a, b, c positive")

    matrix = diagonal_matrix(ring, [ring(c) / ring(diagonal[0]) for c in diagonal])
    matrix.set_immutable()
    return matrix


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


def _points_form(points, form, ring):
    r"""
    Return the matrix over ``ring`` of the quadratic form of a computation on
    ``points``: ``form`` if given, otherwise the form of the hyperbolic plane
    of the points.

    INPUT:

    - ``points`` -- a list of points of hyperbolic planes, or of Klein
      coordinates, which carry no form

    - ``form`` -- a quadratic form as in :func:`tangent_vector_isometry`, or
      ``None``; without it and without points of a hyperbolic plane, the form
      is `x^2 + y^2 - z^2`

    - ``ring`` -- the field of the computation

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_isometry import _minkowski_form, _points_form
        sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
        sage: p = H3.point(0, 0, model="klein")
        sage: _points_form([p, (1, 0)], None, QQ)
        [ 1  0  0]
        [ 0  2  0]
        [ 0  0 -3]
        sage: _points_form([(0, 0), (1, 0)], None, QQ)
        [ 1  0  0]
        [ 0  1  0]
        [ 0  0 -1]
        sage: _points_form([p], [2, 4, -6], AA) == _minkowski_form(AA, [1, 2, -3])
        True

    TESTS::

        sage: _points_form([p, HyperbolicPlane(QQ)(I)], None, QQ)
        Traceback (most recent call last):
        ...
        ValueError: the points lie in hyperbolic planes of different forms
        sage: _points_form([p], [1, 1, -1], QQ)
        Traceback (most recent call last):
        ...
        ValueError: the form differs from the form of the hyperbolic plane of the points

    """
    from flatsurf.geometry.hyperbolic import HyperbolicPoint

    forms = []
    for point in points:
        if isinstance(point, HyperbolicPoint):
            Q = _minkowski_form(ring, point.parent().quadratic_form())
            if Q not in forms:
                forms.append(Q)

    if len(forms) > 1:
        raise ValueError("the points lie in hyperbolic planes of different forms")

    if form is None:
        return forms[0] if forms else _minkowski_form(ring)

    Q = _minkowski_form(ring, form)
    if forms and forms[0] != Q:
        raise ValueError(
            "the form differs from the form of the hyperbolic plane of the points"
        )

    return Q


def _form_reflection(n, Q):
    r"""
    Return the matrix `1 - 2 n n^T Q / \langle n, n\rangle` of the reflection
    of the quadratic form ``Q`` in the vector ``n`` with `\langle n, n\rangle
    = n^T Q n > 0`.

    The reflection fixes the plane `n^\perp`, whose trace on the Klein model
    is the mirror, a geodesic since `n` is spacelike, and preserves the form;
    its determinant is `-1`.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import _form_reflection
        sage: Q = diagonal_matrix(QQ, [1, 2, -3])
        sage: R = _form_reflection(vector(QQ, (1, 1, 0)), Q)
        sage: R
        [ 1/3 -4/3    0]
        [-2/3 -1/3    0]
        [   0    0    1]
        sage: R.transpose() * Q * R == Q, R.determinant()
        (True, -1)

    """
    from sage.matrix.constructor import identity_matrix, matrix

    N = matrix(Q.base_ring(), [n])
    return identity_matrix(Q.base_ring(), 3) - 2 * N.transpose() * N * Q / _minkowski(
        n, n, Q
    )


def tangent_vector_isometry(u, v, ring=None, form=None):
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

    - ``form`` -- a diagonal 3×3 matrix `\mathrm{diag}(a, b, -c)` with `a, b,
      c` positive, or the list of its diagonal, or ``None`` (default:
      ``None``) for the form of the hyperbolic plane of the points, `x^2 +
      y^2 - z^2` for Klein coordinates; the quadratic form whose cone is the
      hyperbolic plane. The points of hyperbolic planes must all have this
      form.

    OUTPUT:

    An immutable 3×3 matrix of determinant `1` over ``ring`` that preserves
    the form, `M^T Q M = Q`; for the default form it lies in `\mathrm{SO}(2,
    1)^+`. The group acts simply transitively on the unit tangent
    bundle; the Klein chart preserves the orientation, so that a ray in the
    chart is a unit tangent vector: the matrix is unique.

    ALGORITHM:

    Write `Q` for the matrix of the form and `\langle u, v\rangle = u^T Q v`.
    The point `p` lifts to `P = (p, 1)` on the cone of timelike vectors, and
    the direction `d` to the tangent vector `T = D - \langle P, D\rangle /
    \langle P, P\rangle P` of the hyperboloid, the projection of
    `D = (d, 0)` onto `P^\perp`. The vector `N = Q^{-1} (P \times T)` satisfies
    `\langle N, X\rangle = \det(P, T, X)`, so that `(P, T, N)` is an
    orthogonal basis with `\langle N, N\rangle = \langle P, P\rangle \langle
    T, T\rangle / \det Q`. With the same construction at the image,
    the isometry is

    .. MATH::

        M = (\lambda P', \mu T', \lambda\mu N') (P, T, N)^{-1}

    with `\lambda^2 = \langle P, P\rangle / \langle P', P'\rangle` and
    `\mu^2 = \langle T, T\rangle / \langle T', T'\rangle`, `\lambda, \mu >
    0`; its determinant is `\lambda^2 \mu^2 \langle N', N'\rangle / \langle
    N, N\rangle = 1`. The inverse of the orthogonal basis is
    `\mathrm{diag}(1/\langle P, P\rangle, 1/\langle T, T\rangle, 1/\langle
    N, N\rangle) (P, T, N)^T Q`. Since `\lambda` is the last entry of `M P`
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

    The form can be changed; the hyperbolic plane is then the disk `x^2 + 2
    y^2 < 3`::

        sage: Q = diagonal_matrix([1, 2, -3])
        sage: M = tangent_vector_isometry(((0, 0), (1, 0)), ((0, 0), (1, 2)), form=Q)
        sage: M
        [ 1/3 -4/3    0]
        [ 2/3  1/3    0]
        [   0    0    1]
        sage: M.transpose() * Q * M == Q
        True
        sage: M.determinant()
        1
        sage: tangent_vector_isometry(((-1, 0), (1, 0)), ((1, 0), (1, 0)), form=Q)
        [2 0 3]
        [0 1 0]
        [1 0 2]
        sage: tangent_vector_isometry(((1, 1), (1, 0)), ((0, 0), (1, 0)), form=Q)
        Traceback (most recent call last):
        ...
        ValueError: (1, 1) is not a finite point of the hyperbolic plane

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
    points = [p, q]

    p = _klein_vector(p)
    q = _klein_vector(q)
    d = vector(d)
    e = vector(e)

    ring = _common_ring([p, q, d, e], ring)

    p, q, d, e = (w.change_ring(ring) for w in (p, q, d, e))

    Q = _points_form(points, form, ring)

    for point in (p, q):
        if not _minkowski(_lift(point), _lift(point), Q) < 0:
            raise ValueError(f"{point} is not a finite point of the hyperbolic plane")

    for direction in (d, e):
        if len(direction) != 2:
            raise ValueError("direction must be a vector of dimension two")
        if direction.is_zero():
            raise ValueError("direction must be nonzero")

    Q_inverse = Q.inverse()

    def frame(point, direction):
        P = vector(ring, _lift(point))
        D = vector(ring, (direction[0], direction[1], 0))
        T = D - (_minkowski(P, D, Q) / _minkowski(P, P, Q)) * P
        N = Q_inverse * P.cross_product(T)
        return P, T, N

    P, T, N = frame(p, d)
    P_, T_, N_ = frame(q, e)

    lambda2 = _minkowski(P, P, Q) / _minkowski(P_, P_, Q)
    lam = _positive_sqrt(lambda2)
    if lam is None:
        raise _not_defined_over(ring, lambda2)

    mu2 = _minkowski(T, T, Q) / _minkowski(T_, T_, Q)
    mu = _positive_sqrt(mu2)
    if mu is None:
        raise _not_defined_over(ring, mu2)

    source = matrix(ring, [P, T, N]).transpose()
    target = matrix(ring, [lam * P_, mu * T_, lam * mu * N_]).transpose()

    inverse = (
        diagonal_matrix(
            ring,
            [
                1 / _minkowski(P, P, Q),
                1 / _minkowski(T, T, Q),
                1 / _minkowski(N, N, Q),
            ],
        )
        * source.transpose()
        * Q
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


def segment_isometry(source, target, ring=None, form=None):
    r"""
    Return the matrix of the orientation preserving isometry of the
    hyperbolic plane that sends the start of ``source`` to the start of
    ``target`` and the end of ``source`` to the end of ``target``.

    INPUT:

    - ``source``, ``target`` -- oriented segments of the hyperbolic plane with
      finite endpoints, or pairs of finite points or of Klein coordinates

    - ``ring`` -- a field or ``None`` (default: ``None``); as in
      :func:`tangent_vector_isometry`

    - ``form`` -- a quadratic form or ``None`` (default: ``None``); as in
      :func:`tangent_vector_isometry`

    OUTPUT:

    An immutable 3×3 matrix of `\mathrm{SO}(2, 1)^+` over ``ring``, or of the
    group of the form.

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

    With another form, the rotation of the corresponding example of
    :func:`tangent_vector_isometry`::

        sage: Q = diagonal_matrix([1, 2, -3])
        sage: segment_isometry(((0, 0), (1, 0)), ((0, 0), (1/3, 2/3)), form=Q)
        [ 1/3 -4/3    0]
        [ 2/3  1/3    0]
        [   0    0    1]

    Points of a hyperbolic plane bring its form::

        sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
        sage: p0 = H3.point(0, 0, model="klein")
        sage: p1 = H3.point(1, 0, model="klein")
        sage: p2 = H3.point(1/3, 2/3, model="klein")
        sage: segment_isometry((p0, p1), (p0, p2))
        [ 1/3 -4/3    0]
        [ 2/3  1/3    0]
        [   0    0    1]
        sage: segment_isometry((p0, H3.point(1/2, 0, model="klein")), (p0, H3.point(0, 1/2, model="klein")))
        Traceback (most recent call last):
        ...
        ValueError: segments do not have the same length

    TESTS:

    The points must have the form of the computation::

        sage: segment_isometry((p0, p1), (H(I), H(2*I)))
        Traceback (most recent call last):
        ...
        ValueError: the points lie in hyperbolic planes of different forms
        sage: segment_isometry((p0, p1), (p0, p2), form=[1, 1, -1])
        Traceback (most recent call last):
        ...
        ValueError: the form differs from the form of the hyperbolic plane of the points

    Other refusals::

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
    points = list(_segment_endpoints(source)) + list(_segment_endpoints(target))
    a, b, c, d = (_klein_vector(point) for point in points)

    ring = _common_ring([a, b, c, d], ring)
    a, b, c, d = (w.change_ring(ring) for w in (a, b, c, d))

    Q = _points_form(points, form, ring)

    # The norm of the lift of a point of the boundary vanishes, so the lengths
    # can only be compared for points inside the disk.
    for point in (a, b, c, d):
        if not _minkowski(_lift(point), _lift(point), Q) < 0:
            raise ValueError(f"{point} is not a finite point of the hyperbolic plane")

    if _cosh2(a, b, Q) != _cosh2(c, d, Q):
        raise ValueError("segments do not have the same length")

    return tangent_vector_isometry((a, b - a), (c, d - c), ring=ring, form=Q)


def segment_midpoint(a, b, ring=None, form=None):
    r"""
    Return the Klein coordinates of the hyperbolic midpoint of the segment
    from ``a`` to ``b``.

    INPUT:

    - ``a``, ``b`` -- finite points of the hyperbolic plane, or pairs or
      vectors of Klein coordinates

    - ``ring`` -- a field or ``None`` (default: ``None``); the field of the
      result, by default the fraction field of the pushout of the rings of
      the coordinates

    - ``form`` -- a quadratic form or ``None`` (default: ``None``); as in
      :func:`tangent_vector_isometry`

    OUTPUT:

    An immutable vector of length two over ``ring``. A :class:`ValueError` is
    raised when the midpoint is not defined over ``ring``, exactly when the
    half turn about it is not (see :func:`segment_isometry`).

    ALGORITHM:

    With `\langle u, v\rangle = u^T Q v` for the matrix `Q` of the form, the
    lifts `A = (a, 1)` and `B = (b, 1)` are rescaled to `A/\sqrt{\langle A,
    A\rangle}` and `B/\sqrt{\langle B, B\rangle}`, which lie on the
    hyperboloid, and their sum is the lift of the midpoint. With
    `\lambda^2 = \langle A, A\rangle / \langle B, B\rangle` and `\lambda >
    0`, the midpoint is `(a + \lambda b) / (1 + \lambda)`.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import segment_midpoint
        sage: segment_midpoint((1/2, 0), (0, 1/2))
        (1/4, 1/4)
        sage: segment_midpoint((0, 0), (1/2, 0))
        Traceback (most recent call last):
        ...
        ValueError: isometry is not defined over Rational Field; adjoin a square root of 4/3
        sage: L.<s3> = NumberField(x^2 - 3, embedding=AA(3).sqrt())
        sage: segment_midpoint((0, 0), (1/2, 0), ring=L)
        (-s3 + 2, 0)
        sage: segment_midpoint((0, 0), (1/2, 0), ring=AA)
        (0.2679491924311227?, 0)

    With another form::

        sage: Q = diagonal_matrix([1, 2, -3])
        sage: segment_midpoint((-1, 0), (-1/3, -1), form=Q)
        (-3/5, -3/5)
        sage: segment_midpoint((0, 0), (1, 0), form=Q)
        Traceback (most recent call last):
        ...
        ValueError: isometry is not defined over Rational Field; adjoin a square root of 3/2

    """
    points = [a, b]
    a = _klein_vector(a)
    b = _klein_vector(b)

    ring = _common_ring([a, b], ring)
    a, b = (w.change_ring(ring) for w in (a, b))

    Q = _points_form(points, form, ring)

    for point in (a, b):
        if not _minkowski(_lift(point), _lift(point), Q) < 0:
            raise ValueError(f"{point} is not a finite point of the hyperbolic plane")

    lambda2 = _minkowski(_lift(a), _lift(a), Q) / _minkowski(_lift(b), _lift(b), Q)
    lam = _positive_sqrt(lambda2)
    if lam is None:
        raise _not_defined_over(ring, lambda2)

    midpoint = (a + lam * b) / (1 + lam)
    midpoint.set_immutable()
    return midpoint


class HyperbolicIsometry(MultiplicativeGroupElement):
    r"""
    An orientation preserving isometry of the hyperbolic plane, given by a
    matrix of `\mathrm{SO}(2, 1)^+`, or of `\mathrm{SO}(Q)^+` for the
    quadratic form `Q` of its group.

    INPUT:

    - ``parent`` -- a :class:`HyperbolicIsometryGroup`

    - ``matrix`` -- a 3×3 matrix of the group of ``parent`` over its base
      ring, immutable; not checked: isometries are built by calling the
      group, which checks its input

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

    def _fixed_vector(self):
        r"""
        Return a vector spanning the kernel of `M - 1`, where `M` is the
        matrix of this isometry.

        OUTPUT:

        A vector over the base ring, not normalised; its last coordinate
        vanishes when `M` is a hyperbolic isometry whose axis is a diameter.
        A :class:`ValueError` is raised for the identity.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: G(matrix([[0, -1], [1, 1]]))._fixed_vector()
            (1, 0, -2)
            sage: G(matrix([[1, 1], [0, 1]]))._fixed_vector()
            (0, 1, 1)
            sage: B = matrix([[5/4, 0, 3/4], [0, 1, 0], [3/4, 0, 5/4]])
            sage: G(B)._fixed_vector()
            (0, 1, 0)
            sage: G(matrix([[2, 1], [1, 1]]))._fixed_vector()
            (1, -2, 0)
            sage: G.one()._fixed_vector()
            Traceback (most recent call last):
            ...
            ValueError: every point is fixed by the identity

        """
        kernel = (self._matrix - 1).right_kernel()
        if kernel.dimension() == 3:
            raise ValueError("every point is fixed by the identity")
        assert kernel.dimension() == 1
        return kernel.basis()[0]

    def fixed_point(self):
        r"""
        Return the fixed point of this isometry in the projective plane of the
        Klein model.

        OUTPUT:

        A point of the hyperbolic plane over the base ring of this isometry:
        the finite fixed point of an elliptic isometry, or the ideal fixed
        point of a parabolic isometry. The fixed point of a hyperbolic
        isometry is not in the closed disk (it is the pole of its axis) and
        cannot be represented yet: a :class:`NotImplementedError` is raised.
        A :class:`ValueError` is raised for the identity.

        ALGORITHM:

        A fixed point of the action on the projective plane is a line of
        fixed vectors, the kernel of `M - 1`, spanned by `v`. The sign of the
        Minkowski norm `\langle v, v\rangle` is negative for a finite point,
        zero for an ideal point and positive outside of the closed disk.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: rho = G(matrix([[0, -1], [1, 1]]))
            sage: rho.trace()
            0
            sage: rho.fixed_point()
            (-1/2, 0)
            sage: G(matrix([[1, 1], [0, 1]])).fixed_point()
            ∞

        The fixed point is the midpoint of a segment for the half turn
        exchanging its endpoints::

            sage: from flatsurf.geometry.hyperbolic_isometry import segment_isometry, segment_midpoint
            sage: s = ((1/2, 0), (0, 1/2))
            sage: g = G(segment_isometry(s, (s[1], s[0])))
            sage: g.fixed_point()
            (1/4, 1/4)
            sage: segment_midpoint(*s)
            (1/4, 1/4)

        TESTS::

            sage: B = matrix([[5/4, 0, 3/4], [0, 1, 0], [3/4, 0, 5/4]])
            sage: G(B).fixed_point()
            Traceback (most recent call last):
            ...
            NotImplementedError: the fixed point of a hyperbolic isometry is outside of the closed disk and cannot be represented yet
            sage: G(matrix([[2, 1], [1, 1]])).fixed_point()
            Traceback (most recent call last):
            ...
            NotImplementedError: the fixed point of a hyperbolic isometry is outside of the closed disk and cannot be represented yet
            sage: G.one().fixed_point()
            Traceback (most recent call last):
            ...
            ValueError: every point is fixed by the identity

        With another form, the point lies in the hyperbolic plane of that
        form::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: p = G3([[-1/3, 4/3, 0], [-2/3, -1/3, 0], [0, 0, 1]]).fixed_point()
            sage: p
            (0, 0)
            sage: p.parent()
            Hyperbolic Plane over Rational Field with form x^2 + 2*y^2 - 3*z^2

        """
        v = self._fixed_vector()
        norm = _minkowski(v, v, self.parent()._form)
        if norm > 0:
            raise NotImplementedError(
                "the fixed point of a hyperbolic isometry is outside of the closed disk and cannot be represented yet"
            )

        H = self.parent().hyperbolic_plane()
        return H.point(v[0] / v[2], v[1] / v[2], model="klein")

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
        Return the inverse of this isometry, `Q^{-1} M^T Q`.

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

        With another form::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: ~G3(matrix([[2, 0, 3], [0, 1, 0], [1, 0, 2]]))
            [ 2  0 -3]
            [ 0  1  0]
            [-1  0  2]

        """
        Q = self.parent().quadratic_form()
        M = Q.inverse() * self._matrix.transpose() * Q
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
        otherwise a subset of the hyperbolic plane over the pushout of the
        base ring of this isometry and the base ring of ``x``.

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

        An isometry over a larger ring moves a set over a smaller ring to the
        larger ring::

            sage: L = NumberField(x^2 - 3, "s3", embedding=AA(3).sqrt())
            sage: s3 = L.gen()
            sage: GL = HyperbolicIsometryGroup(L)
            sage: g = GL(matrix(L, [[s3/2, -1/2, 0], [1/2, s3/2, 0], [0, 0, 1]]))
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: Q = g(P)
            sage: Q
            Hyperbolic polygon with vertices [I, 1/3*s3 + 2/3*s3*I, -1/13*s3 - 4/13 + (8/13*s3 + 6/13)*I]
            sage: Q.parent().base_ring()
            Number Field in s3 with defining polynomial x^2 - 3 with s3 = 1.732050807568878?
            sage: Q == g(P.change_ring(L))
            True
            sage: g(H(2*I)).parent().base_ring() is L
            True
            sage: GL.one()(P).parent().base_ring() is L
            True

        Two fields embedded in the algebraic reals meet there::

            sage: M = NumberField(x^2 - 2, "s2", embedding=AA(2).sqrt())
            sage: s2 = M.gen()
            sage: PM = HyperbolicPlane(M).polygon(vertices=[(0, 0), (s2/4, 0), (0, s2/4)])
            sage: R = g(PM)
            sage: R
            Hyperbolic polygon with vertices [I, 0.3719357993129545? + 1.136282636010346?*I, -0.2547898295538365? + 1.348221051175973?*I]
            sage: R.parent().base_ring()
            Algebraic Real Field
            sage: R == HyperbolicIsometryGroup(AA)(g)(PM.change_ring(AA))
            True

        TESTS:

        A point cannot be multiplied by an isometry on the right::

            sage: H(I) * G.one()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Hyperbolic Plane over Rational Field' and 'Group of orientation preserving isometries of the hyperbolic plane over Rational Field'

        Fields embedded in the reals of Sage have no common parent::

            sage: L = NumberField(x^2 - 3, "s3", embedding=RR(1.7))
            sage: M = NumberField(x^2 - 2, "s2", embedding=RR(1.4))
            sage: s3 = L.gen()
            sage: g = HyperbolicIsometryGroup(L)(matrix(L, [[s3/2, -1/2, 0], [1/2, s3/2, 0], [0, 0, 1]]))
            sage: p = HyperbolicPlane(M)(I)
            sage: g(p)
            Traceback (most recent call last):
            ...
            TypeError: no common base ring for the isometry over Number Field in s3 with defining polynomial x^2 - 3 with s3 = 1.732050807568878? and the set over Number Field in s2 with defining polynomial x^2 - 2 with s2 = 1.414213562373095?

        The set must lie in a hyperbolic plane of the form of the isometry::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: B = G3(matrix([[2, 0, 3], [0, 1, 0], [1, 0, 2]]))
            sage: B(H3.point(0, 0, model="klein"))
            (3/2, 0)
            sage: B(H(I))
            Traceback (most recent call last):
            ...
            TypeError: the isometry preserves x^2 + 2*y^2 - 3*z^2 and the plane has the form x^2 + y^2 - z^2

        .. SEEALSO::

            :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.isomorphism`
            to map a set or an isometry to a hyperbolic plane of another form

        """
        from flatsurf.geometry.hyperbolic import HyperbolicConvexSet
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

        if isinstance(x, (HyperbolicConvexSet, HyperbolicPolygon)):
            from sage.categories.pushout import pushout
            from sage.structure.coerce_exceptions import CoercionException

            K_g = self._matrix.base_ring()
            K_x = x.parent().base_ring()
            try:
                R = pushout(K_g, K_x)
            except (CoercionException, TypeError) as error:
                raise TypeError(
                    f"no common base ring for the isometry over {K_g} and the set over {K_x}"
                ) from error

            self.parent()._require_form_of(x.parent())

            if R is not K_x and R != K_x:
                x = x.change_ring(R)

            g = self
            if R is not K_g and R != K_g:
                g = HyperbolicIsometryGroup(R, form=self.parent()._form)(self)

            if isinstance(x, HyperbolicPolygon):
                return x.apply_isometry(g)

            return x.apply_isometry(g._matrix, model="klein")

        from sage.modules.free_module_element import FreeModuleElement, vector

        if not isinstance(x, FreeModuleElement):
            x = vector(x)

        return self * x

    def _half_lift(self):
        r"""
        Return a lift of this isometry to the spin group of the quadratic
        form, up to a scalar whose square lies in the base ring, and that
        square.

        OUTPUT:

        For the form `x^2 + y^2 - z^2`, a pair ``(B, s)`` of a 2×2 matrix
        ``B`` over the base ring and a nonzero element ``s`` of the base ring
        such that ``B = aA`` and ``s = a^2``, where `A` is the matrix of
        :meth:`psl2_matrix` and `a` a nonzero scalar. In particular
        ``B.det() == s``.

        For the form `x^2 + b y^2 - c z^2`, a pair ``(x, s)`` of a nonzero
        element ``x`` of the even Clifford algebra of the form, the
        quaternion algebra `(-b, bc)` over the base ring, and its reduced
        norm ``s``. The element `x / \sqrt{s}` of reduced norm `1` acts on the
        vectors of the form by the conjugation `v \mapsto x v x^{-1}` as this
        isometry; the lifts of the isometries `g` and `h` multiply to a
        nonzero scalar times the lift of `gh`. Its reduced trace `t`
        satisfies `t^2 / s = \mathrm{tr}(M) + 1` for the matrix `M` of this
        isometry, as the trace of `A` on the standard form.

        ALGORITHM:

        For the form `x^2 + y^2 - z^2`, the matrix of this isometry is
        quadratic in the entries `a, b, c, d` of `A`: its entries give `a^2,
        b^2, c^2, d^2` and the products of two of them. As `ad - bc = 1`,
        either `a` or `b` is nonzero. If `s = a^2` is nonzero, then `aA =
        (a^2, ab; ac, ad)`; otherwise `bA = (ab, b^2; bc, bd)` and `s = b^2`.
        No square root is taken.

        For another form, with `e_1, e_2, e_3` the basis of the vectors, the
        generators of the even Clifford algebra are `i = e_1 e_2`, `j = e_2
        e_3` and `k = ij = b e_1 e_3`. The element

        .. MATH::

            x = (1 + \mathrm{tr} M) + (M_{01}/b - M_{10}) i
            - (M_{21}/b + M_{12}/c) j - (M_{20}/b + M_{02}/(bc)) k

        (indices from `0`) is a multiple of the lift, by the classical
        formula of a lift as `1` plus the sum of the products of the dual
        basis vectors with their images by `M`, read in the basis `1, i, j,
        k`. It vanishes exactly for a half turn, when `\mathrm{tr} M = -1`
        and `QM` is symmetric: the lift of the half turn about the point `v`
        is then the product of `v` with the volume element, `v_0 j - v_1 k -
        c v_2 i`, for a nonzero vector `v` of the kernel of `M - 1`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: g = G(matrix([[5/4, 0, 3/4], [0, 1, 0], [3/4, 0, 5/4]]))
            sage: B, s = g._half_lift()
            sage: B
            [9/8 3/8]
            [3/8 9/8]
            sage: s
            9/8
            sage: B.det() == s
            True
            sage: B / AA(s).sqrt() == g.psl2_matrix(AA)
            True

        With the form `x^2 + 2y^2 - 3z^2`, the lifts of the identity, of the
        half turn about the origin, of a rotation about the origin and of a
        hyperbolic isometry::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: G3.one()._half_lift()
            (4, 16)
            sage: G3(diagonal_matrix([-1, -1, 1]))._half_lift()
            (-3*i, 18)
            sage: G3(matrix([[-1/3, 4/3, 0], [-2/3, -1/3, 0], [0, 0, 1]]))._half_lift()
            (4/3 + 4/3*i, 16/3)
            sage: h = G3(matrix([[19/13, 0, -24/13], [0, 1, 0], [-8/13, 0, 19/13]]))
            sage: x, s = h._half_lift()
            sage: x, s
            (64/13 + 8/13*k, 256/13)
            sage: x.parent()
            Quaternion Algebra (-2, 6) with base ring Rational Field
            sage: x.reduced_trace()^2 / s == h.trace() + 1
            True

        """
        from sage.matrix.constructor import matrix

        K = self.parent().base_ring()

        if self.parent()._form is not None:
            from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra

            Q = self.parent().quadratic_form()
            b, c = Q[1, 1], -Q[2, 2]
            i, j, k = QuaternionAlgebra(K, -b, b * c).gens()
            M = self._matrix
            t = 1 + M.trace()
            x = (
                t
                + (M[0, 1] / b - M[1, 0]) * i
                - (M[2, 1] / b + M[1, 2] / c) * j
                - (M[2, 0] / b + M[0, 2] / (b * c)) * k
            )
            if not t:
                v = (M - 1).right_kernel().basis()[0]
                x = v[0] * j - v[1] * k - c * v[2] * i
            return x, x.reduced_norm()

        m00, m01, m02, m10, m11, m12, m20, m21, m22 = self._matrix.list()
        two = K(2)

        a2 = (m11 + m12 + m21 + m22) / two
        if a2:
            square = a2
            entries = [a2, (m10 + m20) / two, (m01 + m02) / two, (m00 + 1) / two]
        else:
            square = (m12 + m22 - m11 - m21) / two
            entries = [K.zero(), square, (m00 - 1) / two, (m02 - m01) / two]
        return matrix(K, 2, entries), square

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

        The upper half plane model needs the form `x^2 + y^2 - z^2`::

            sage: HyperbolicIsometryGroup(QQ, form=[1, 2, -3]).one().psl2_matrix()
            Traceback (most recent call last):
            ...
            ValueError: the half-plane model needs the form x^2 + y^2 - z^2; see standard_form_isometry()

        """
        from sage.matrix.constructor import matrix

        if self.parent()._form is not None:
            self.parent().hyperbolic_plane()._require_standard_form()

        K = self.parent().base_ring()
        if ring is None:
            ring = K

        half, square = self._half_lift()
        entries = half.list()

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
      base ring has a pushout with the base ring of ``G``, with the quadratic
      form of ``G``

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup, HyperbolicIsometryPlaneAction
        sage: G = HyperbolicIsometryGroup(QQ)
        sage: H = HyperbolicPlane(QQ)
        sage: HyperbolicIsometryPlaneAction(G, H)
        Left action by Group of orientation preserving isometries of the hyperbolic plane over Rational Field on Hyperbolic Plane over Rational Field

    """

    def __init__(self, G, H):
        from sage.categories.pushout import pushout

        from flatsurf.geometry.hyperbolic import HyperbolicPlane

        self._codomain = HyperbolicPlane(
            pushout(G.base_ring(), H.base_ring()), form=G._form
        )
        super().__init__(G, H, True, operator.mul)

    def codomain(self):
        r"""
        Return the parent of the images of this action: the hyperbolic plane
        over the pushout of the base rings of the group and of the plane,
        with the quadratic form of the group.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup, HyperbolicIsometryPlaneAction
            sage: A = HyperbolicIsometryPlaneAction(HyperbolicIsometryGroup(QQ), HyperbolicPlane(AA))
            sage: A.codomain()
            Hyperbolic Plane over Algebraic Real Field
            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: HyperbolicIsometryPlaneAction(G3, HyperbolicPlane(AA, form=[1, 2, -3])).codomain()
            Hyperbolic Plane over Algebraic Real Field with form x^2 + 2*y^2 - 3*z^2

        """
        return self._codomain

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


class HyperbolicIsometryGroupFunctor(ConstructionFunctor):
    r"""
    The construction functor of :class:`HyperbolicIsometryGroup`, from the
    fields to the groups.

    It lets the coercion system find the pushout of two groups of isometries
    over different fields, which is the group over the pushout of the fields.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
        sage: F, K = HyperbolicIsometryGroup(QQ).construction()
        sage: F
        HyperbolicIsometryGroupFunctor
        sage: F(AA)
        Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field

    The functor of a group with another quadratic form carries the form::

        sage: F, K = HyperbolicIsometryGroup(QQ, form=[1, 2, -3]).construction()
        sage: F
        HyperbolicIsometryGroupFunctor(x^2 + 2*y^2 - 3*z^2)
        sage: F(AA)
        Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field with form x^2 + 2*y^2 - 3*z^2

    """

    rank = 10

    def __init__(self, form=None):
        from sage.categories.fields import Fields
        from sage.categories.groups import Groups

        # None for the standard form, otherwise the normalised matrix of the
        # form over the base ring of the group.
        self._form = form
        super().__init__(Fields(), Groups())

    def _apply_functor(self, R):
        r"""
        Return the group of isometries over ``R``, with the quadratic form of
        this functor.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: F = HyperbolicIsometryGroup(QQ).construction()[0]
            sage: F(QQ) is HyperbolicIsometryGroup(QQ)
            True
            sage: F = HyperbolicIsometryGroup(QQ, form=[1, 2, -3]).construction()[0]
            sage: F(AA) is HyperbolicIsometryGroup(AA, form=[1, 2, -3])
            True

        """
        if self._form is None:
            return HyperbolicIsometryGroup(R)

        return HyperbolicIsometryGroup(R, form=self._form.diagonal())

    def __eq__(self, other):
        r"""
        Return whether ``other`` is such a functor with the same quadratic
        form.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: F = HyperbolicIsometryGroup(QQ).construction()[0]
            sage: F == HyperbolicIsometryGroup(AA).construction()[0]
            True
            sage: F == 1
            False
            sage: F3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3]).construction()[0]
            sage: F3 == F
            False
            sage: F3 == HyperbolicIsometryGroup(AA, form=[1, 2, -3]).construction()[0]
            True

        """
        if not isinstance(other, HyperbolicIsometryGroupFunctor):
            return False

        if self._form is None or other._form is None:
            return self._form is None and other._form is None

        return self._form == other._form

    def __ne__(self, other):
        r"""
        Return whether ``other`` is not such a functor.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: F = HyperbolicIsometryGroup(QQ).construction()[0]
            sage: F != HyperbolicIsometryGroup(AA).construction()[0]
            False

        """
        return not self == other

    def __hash__(self):
        r"""
        Return a hash value compatible with ``==``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: F = HyperbolicIsometryGroup(QQ).construction()[0]
            sage: hash(F) == hash(HyperbolicIsometryGroup(AA).construction()[0])
            True

        """
        return hash(HyperbolicIsometryGroupFunctor)

    def _repr_(self):
        r"""
        Return a printable representation of this functor.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).construction()[0]
            HyperbolicIsometryGroupFunctor

        """
        if self._form is None:
            return "HyperbolicIsometryGroupFunctor"

        return f"HyperbolicIsometryGroupFunctor({_form_repr(self._form)})"


class HyperbolicIsometryGroup(UniqueRepresentation, Group):
    r"""
    The group of orientation preserving isometries of the hyperbolic plane
    whose matrices in `\mathrm{SO}(Q)^+` have entries in ``base_ring``.

    `\mathrm{SO}(Q)^+` is the group of the matrices `M` with `M^T Q M = Q`,
    determinant `1` and `M_{2,2} > 0`, for the matrix `Q` of the quadratic
    form, `\mathrm{SO}(2, 1)^+` for the default form `x^2 + y^2 - z^2`.

    INPUT:

    - ``base_ring`` -- an ordered field

    - ``form`` -- a diagonal 3×3 matrix `\mathrm{diag}(a, b, -c)` with `a, b,
      c` positive, or the list of its diagonal, or ``None`` (default:
      ``None``) for `x^2 + y^2 - z^2`; normalised as in
      :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane`

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
        sage: G = HyperbolicIsometryGroup(QQ)
        sage: G
        Group of orientation preserving isometries of the hyperbolic plane over Rational Field

    The group of the form `x^2 + 2y^2 - 3z^2`, whose Klein model is the disk
    `x^2 + 2y^2 < 3`::

        sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
        sage: G3
        Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2
        sage: G3.hyperbolic_plane()
        Hyperbolic Plane over Rational Field with form x^2 + 2*y^2 - 3*z^2
        sage: [g.trace() for g in G3.some_elements()]
        [3, -1, 1/3, 51/13]

    TESTS::

        sage: TestSuite(G).run()
        sage: TestSuite(HyperbolicIsometryGroup(AA)).run()
        sage: TestSuite(G3).run()
        sage: TestSuite(HyperbolicIsometryGroup(QQ, form=[1, 3, -1])).run()

    The form is normalised::

        sage: HyperbolicIsometryGroup(QQ, form=[2, 2, -2]) is G
        True
        sage: HyperbolicIsometryGroup(QQ, form=[2, 4, -6]) is G3
        True
        sage: loads(dumps(G3)) is G3
        True

    """

    Element = HyperbolicIsometry

    @staticmethod
    def __classcall__(cls, base_ring, form=None):
        r"""
        Create the group with a normalised form, ``None`` for the form `x^2 +
        y^2 - z^2` so that the standard group keeps its key.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ, form=diagonal_matrix([1, 2, -3])) is HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            True

        """
        if form is not None:
            if isinstance(form, (list, tuple)):
                form = [base_ring(c) for c in form]
            form = _minkowski_form(base_ring, form)
            if form == _minkowski_form(base_ring):
                form = None

        if form is None:
            return super().__classcall__(cls, base_ring)

        return super().__classcall__(cls, base_ring, form)

    def __init__(self, base_ring, form=None):
        from sage.categories.groups import Groups

        self._base_ring = base_ring
        # None for the standard form, so that the helpers of hyperbolic_polygon
        # evaluate it without any division.
        self._form = form
        self._quadratic_form = _minkowski_form(base_ring, form)
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

    def quadratic_form(self):
        r"""
        Return the matrix of the quadratic form preserved by the matrices of
        this group.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).quadratic_form()
            [ 1  0  0]
            [ 0  1  0]
            [ 0  0 -1]
            sage: HyperbolicIsometryGroup(QQ, form=[2, 4, -6]).quadratic_form()
            [ 1  0  0]
            [ 0  2  0]
            [ 0  0 -3]

        """
        return self._quadratic_form

    def hyperbolic_plane(self):
        r"""
        Return the hyperbolic plane over the base ring with the quadratic form
        of this group, on which this group acts.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).hyperbolic_plane()
            Hyperbolic Plane over Rational Field
            sage: HyperbolicIsometryGroup(AA, form=[1, 2, -3]).hyperbolic_plane()
            Hyperbolic Plane over Algebraic Real Field with form x^2 + 2*y^2 - 3*z^2

        """
        from flatsurf.geometry.hyperbolic import HyperbolicPlane

        return HyperbolicPlane(self._base_ring, form=self._form)

    def _require_form_of(self, plane):
        r"""
        Raise an error if the hyperbolic ``plane`` does not have the quadratic
        form of this group, on a base ring that has a pushout with the base
        ring of this group.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: G3._require_form_of(HyperbolicPlane(AA, form=[1, 2, -3]))
            sage: G3._require_form_of(HyperbolicPlane(QQ))
            Traceback (most recent call last):
            ...
            TypeError: the isometry preserves x^2 + 2*y^2 - 3*z^2 and the plane has the form x^2 + y^2 - z^2

        """
        if self._form is None and plane._form is None:
            return

        if self._form is None or plane._form is None or self._form != plane._form:
            raise TypeError(
                f"the isometry preserves {_form_repr(self._quadratic_form)} and the plane has the form {_form_repr(plane.quadratic_form())}"
            )

    def _repr_(self):
        r"""
        Return a printable representation of this group.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(AA)
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field
            sage: HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2

        """
        if self._form is not None:
            return f"Group of orientation preserving isometries of the hyperbolic plane over {self._base_ring} with form {_form_repr(self._quadratic_form)}"

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

        With another form `\mathrm{diag}(1, b, -c)`, the hyperbolic isometry
        of :meth:`some_elements`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ).an_element()
            [ 7 -4  8]
            [ 4 -1  4]
            [ 8 -4  9]
            sage: HyperbolicIsometryGroup(QQ, form=[1, 2, -3]).an_element()
            [ 19/13      0 -24/13]
            [     0      1      0]
            [ -8/13      0  19/13]

        """
        if self._form is not None:
            return self.some_elements()[-1]

        return self([[7, -4, 8], [4, -1, 4], [8, -4, 9]])

    def some_elements(self):
        r"""
        Return some isometries of this group: the identity, an elliptic, a
        parabolic and two hyperbolic isometries.

        With another form `\mathrm{diag}(1, b, -c)`: the identity, the half
        turn about the origin, a rotation and a hyperbolic isometry. An
        anisotropic form has no ideal point over the base ring, so that there
        is no parabolic isometry in general.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: [g.trace() for g in HyperbolicIsometryGroup(QQ).some_elements()]
            [3, 1, 3, 15, 7/2]
            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: G3.some_elements()[2]
            [-1/3  4/3    0]
            [-2/3 -1/3    0]
            [   0    0    1]
            sage: [g.trace() for g in G3.some_elements()]
            [3, -1, 1/3, 51/13]

        ALGORITHM:

        With another form, the rotation and the hyperbolic isometry are
        products of two reflections `R_n = 1 - 2 n n^T Q / \langle n,
        n\rangle` in spacelike vectors `n`: `R_{e_1} R_{(1, 1, 0)}`, whose
        mirrors meet at the origin, and `R_{e_1} R_{(1, 0, x_0/c)}` with `x_0
        = c/(1 + c)`, whose mirrors are the disjoint lines `x = 0` and `x =
        x_0`, since `x_0^2 < c`. Both are defined over every base ring.

        """
        from sage.all import QQ

        if self._form is not None:
            from sage.modules.free_module_element import vector

            K = self._base_ring
            Q = self._quadratic_form
            c = -Q[2, 2]
            x0 = c / (1 + c)
            mirror = _form_reflection(vector(K, (1, 0, 0)), Q)
            return [
                self.one(),
                self([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),
                self(mirror * _form_reflection(vector(K, (1, 1, 0)), Q)),
                self(mirror * _form_reflection(vector(K, (1, 0, x0 / c)), Q)),
            ]

        return [
            self.one(),
            self([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            self(
                [[1, -1, 1], [1, QQ((1, 2)), QQ((1, 2))], [1, -QQ((1, 2)), QQ((3, 2))]]
            ),
            self.an_element(),
            self([[QQ((5, 4)), 0, QQ((3, 4))], [0, 1, 0], [QQ((3, 4)), 0, QQ((5, 4))]]),
        ]

    def construction(self):
        r"""
        Return the construction of this group: the functor that sends a field
        to its group of isometries, and the base ring.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: F, K = G.construction()
            sage: F, K
            (HyperbolicIsometryGroupFunctor, Rational Field)
            sage: F(K) is G
            True

        The pushout of the groups over two number fields is the group over
        the algebraic reals, and their elements multiply there::

            sage: from sage.categories.pushout import pushout
            sage: L = NumberField(x^2 - 3, "s3", embedding=AA(3).sqrt())
            sage: M = NumberField(x^2 - 2, "s2", embedding=AA(2).sqrt())
            sage: G_L = HyperbolicIsometryGroup(L)
            sage: G_M = HyperbolicIsometryGroup(M)
            sage: pushout(G_L, G_M)
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field
            sage: (G_L.one() * G_M.one()).parent()
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field

        The functor carries the quadratic form::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: F, K = G3.construction()
            sage: F(K) is G3
            True

        """
        return HyperbolicIsometryGroupFunctor(self._form), self._base_ring

    def _pushout_(self, other):
        r"""
        Return the group of isometries over the pushout of the base rings of
        this group and of ``other``, or ``None`` if ``other`` is not a group
        of isometries or if the two groups do not have the same quadratic
        form.

        The base rings are joined directly: the generic pushout of the two
        constructions applies the functors of the fields in an order it
        cannot decide when neither field coerces into the other.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: HyperbolicIsometryGroup(QQ)._pushout_(HyperbolicIsometryGroup(AA))
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field
            sage: HyperbolicIsometryGroup(QQ)._pushout_(QQ) is None
            True

        With another form::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: G3._pushout_(HyperbolicIsometryGroup(AA, form=[1, 2, -3]))
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field with form x^2 + 2*y^2 - 3*z^2
            sage: (G3.one() * HyperbolicIsometryGroup(AA, form=[1, 2, -3]).one()).parent()
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field with form x^2 + 2*y^2 - 3*z^2
            sage: G3._pushout_(HyperbolicIsometryGroup(QQ)) is None
            True

        Isometries of different forms do not multiply::

            sage: G3.one() * HyperbolicIsometryGroup(QQ).one()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: ...

        """
        if not isinstance(other, HyperbolicIsometryGroup):
            return None

        from sage.categories.pushout import pushout

        R = pushout(self._base_ring, other.base_ring())

        if self._form is None and other._form is None:
            return HyperbolicIsometryGroup(R)

        if self._form is None or other._form is None:
            return None

        form = self.hyperbolic_plane()._diagonal_form_in(R)
        if form is None or form != other.hyperbolic_plane()._diagonal_form_in(R):
            return None

        return HyperbolicIsometryGroup(R, form=form)

    def _coerce_map_from_(self, other):
        r"""
        Return whether there is a coercion from ``other`` to this group: the
        groups of isometries with the same quadratic form over a field that
        coerces into the base ring.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: G = HyperbolicIsometryGroup(QQ)
            sage: HyperbolicIsometryGroup(AA).has_coerce_map_from(G)
            True
            sage: G.has_coerce_map_from(HyperbolicIsometryGroup(AA))
            False
            sage: (HyperbolicIsometryGroup(AA).one() * G.an_element()).parent()
            Group of orientation preserving isometries of the hyperbolic plane over Algebraic Real Field

        With another form::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: HyperbolicIsometryGroup(AA, form=[1, 2, -3]).has_coerce_map_from(G3)
            True
            sage: HyperbolicIsometryGroup(AA).has_coerce_map_from(G3)
            False

        """
        if isinstance(other, HyperbolicIsometryGroup):
            if not self._base_ring.has_coerce_map_from(other.base_ring()):
                return False

            if self._form is None or other._form is None:
                return self._form is None and other._form is None

            return other.quadratic_form().change_ring(self._base_ring) == self._form

        return False

    def _element_constructor_(self, x):
        r"""
        Return the isometry given by ``x``.

        INPUT:

        - ``x`` -- an isometry, a 3×3 matrix of `\mathrm{SO}(2, 1)^+`, or a
          2×2 matrix of positive determinant acting on the upper half plane
          as in
          :meth:`~flatsurf.geometry.hyperbolic.HyperbolicConvexSet.apply_isometry`;
          a matrix may be given as the list (or tuple) of its rows

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
            sage: G([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            [ 0 -1  0]
            [ 1  0  0]
            [ 0  0  1]

        Something that is neither an isometry nor a matrix is refused::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: G(SimilarityGroup(QQ)((1, 2)))
            Traceback (most recent call last):
            ...
            ValueError: element in Similarity group over Rational Field cannot be used to create element in Group of orientation preserving isometries of the hyperbolic plane over Rational Field

        With another form, the matrix must preserve that form; a 2×2 matrix
        acts on the upper half plane, which needs the form `x^2 + y^2 -
        z^2`::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: B = matrix([[2, 0, 3], [0, 1, 0], [1, 0, 2]])
            sage: G3(B).trace()
            5
            sage: G(B)
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(2, 1)+
            sage: G3([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(x^2 + 2*y^2 - 3*z^2)+
            sage: G3(matrix([[2, 0], [0, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: the half-plane model needs the form x^2 + y^2 - z^2; see standard_form_isometry()

        An isometry of another form is refused::

            sage: G3(G.an_element())
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(x^2 + 2*y^2 - 3*z^2)+
            sage: G(G3(B))
            Traceback (most recent call last):
            ...
            ValueError: matrix is not in SO(2, 1)+

        """
        from sage.matrix.constructor import matrix

        K = self._base_ring

        if isinstance(x, HyperbolicIsometry):
            if self._form is None and x.parent()._form is None:
                M = x.matrix().change_ring(K)
                M.set_immutable()
                return self.element_class(self, M)

            x = x.matrix()

        if not isinstance(x, Matrix):
            if not isinstance(x, (list, tuple)):
                from sage.structure.element import parent

                raise ValueError(
                    "element in {} cannot be used to create element in {}".format(
                        parent(x), self
                    )
                )

            x = matrix(K, x)

        if x.dimensions() == (2, 2):
            if self._form is not None:
                self.hyperbolic_plane()._require_standard_form()

            det = x.determinant()
            if not det > 0:
                raise ValueError(
                    "matrix does not define an orientation preserving isometry"
                )

            from flatsurf.geometry.hyperbolic import HyperbolicPlane

            M = HyperbolicPlane(K)._isometry_gl2_to_sim12(x.change_ring(K)) / K(det)
        elif x.dimensions() == (3, 3):
            M = x.change_ring(K)
            J = self._quadratic_form
            if M.transpose() * J * M != J or M.determinant() != 1 or not M[2, 2] > 0:
                if self._form is None:
                    raise ValueError("matrix is not in SO(2, 1)+")
                raise ValueError(f"matrix is not in SO({_form_repr(J)})+")
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
        hyperbolic planes over a field that has a pushout with the base ring,
        see :class:`HyperbolicIsometryPlaneAction`.

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

        A point over a smaller field is moved to the field of the isometry::

            sage: L = NumberField(x^2 - 3, "s3", embedding=AA(3).sqrt())
            sage: s3 = L.gen()
            sage: g = HyperbolicIsometryGroup(L)(matrix(L, [[s3/2, -1/2, 0], [1/2, s3/2, 0], [0, 0, 1]]))
            sage: p = g * HyperbolicPlane(QQ)(2*I)
            sage: p.parent().base_ring() is L
            True

        TESTS:

        Fields embedded in the reals of Sage have no common parent::

            sage: L = NumberField(x^2 - 3, "s3", embedding=RR(1.7))
            sage: M = NumberField(x^2 - 2, "s2", embedding=RR(1.4))
            sage: s3 = L.gen()
            sage: g = HyperbolicIsometryGroup(L)(matrix(L, [[s3/2, -1/2, 0], [1/2, s3/2, 0], [0, 0, 1]]))
            sage: g * HyperbolicPlane(M)(I)
            Traceback (most recent call last):
            ...
            TypeError: ...
            sage: HyperbolicIsometryGroup(M).one() * g
            Traceback (most recent call last):
            ...
            TypeError: ...

        A group with another form acts on the planes of that form only::

            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: B = G3(matrix([[2, 0, 3], [0, 1, 0], [1, 0, 2]]))
            sage: B * HyperbolicPlane(QQ, form=[1, 2, -3]).point(0, 0, model="klein")
            (3/2, 0)
            sage: B * HyperbolicPlane(QQ)(I)
            Traceback (most recent call last):
            ...
            TypeError: the isometry preserves x^2 + 2*y^2 - 3*z^2 and the plane has the form x^2 + y^2 - z^2

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
                from sage.categories.pushout import pushout
                from sage.structure.coerce_exceptions import CoercionException

                try:
                    pushout(self._base_ring, S.base_ring())
                except (CoercionException, TypeError) as error:
                    # Falling back to the generic action would reach the
                    # point's own action, which does not accept isometries.
                    raise TypeError(
                        f"no common base ring for the isometry over {self._base_ring} and the set over {S.base_ring()}"
                    ) from error

                self._require_form_of(S)

                return HyperbolicIsometryPlaneAction(self, S)

        return super()._get_action_(S, op, self_on_left)
