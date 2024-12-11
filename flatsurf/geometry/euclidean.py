r"""
A loose collection of tools for Euclidean geometry in the plane.

.. SEEALSO::

    :mod:`flatsurf.geometry.circle` for everything specific to circles in the plane

"""

######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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


def is_cosine_sine_of_rational(cos, sin, scaled=False):
    r"""
    Check whether the given pair is a cosine and sine of a same rational angle.

    INPUT:

    - ``cos`` -- a number

    - ``sin`` -- a number

    - ``scaled`` -- a boolean (default: ``False``); whether to allow ``cos``
      and ``sin`` to be scaled by the same positive algebraic number

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_cosine_sine_of_rational

        sage: c = s = AA(sqrt(2))/2
        sage: is_cosine_sine_of_rational(c, s)
        True

        sage: c = AA(sqrt(3))/2
        sage: s = AA(1/2)
        sage: is_cosine_sine_of_rational(c, s)
        True

        sage: c = AA(sqrt(5)/2)
        sage: s = (1 - c**2).sqrt()
        sage: c**2 + s**2
        1.000000000000000?
        sage: is_cosine_sine_of_rational(c, s)
        False

        sage: c = (AA(sqrt(5)) + 1)/4
        sage: s = (1 - c**2).sqrt()
        sage: is_cosine_sine_of_rational(c, s)
        True

        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: is_cosine_sine_of_rational(K.zero(), -K.one())
        True

    TESTS::

        sage: from pyexactreal import ExactReals # optional: pyexactreal  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: R = ExactReals() # optional: pyexactreal
        sage: is_cosine_sine_of_rational(R.one(), R.zero()) # optional: pyexactreal
        True

    """
    from sage.all import AA

    # We cannot check in AA due to https://github.com/flatsurf/exact-real/issues/172
    # We just trust that non-algebraic elements won't allow conversion to AA.
    # if cos not in AA:
    #     return False
    # if sin not in AA:
    #     return False

    if not scaled:
        if cos**2 + sin**2 != 1:
            return False

    try:
        cos = AA(cos)
    except ValueError:
        # This is a replacement for the "in AA" checked disabled above.
        return False

    cos = cos.as_number_field_element(embedded=True)
    # We need an explicit conversion to the number field due to https://github.com/sagemath/sage/issues/35613
    cos = cos[0](cos[1])

    try:
        sin = AA(sin)
    except ValueError:
        # This is a replacement for the "in AA" checked disabled above.
        return False
    sin = sin.as_number_field_element(embedded=True)
    # We need an explicit conversion to the number field due to https://github.com/sagemath/sage/issues/35613
    sin = sin[0](sin[1])

    from sage.all import ComplexBallField

    CBF = ComplexBallField(53)

    x = CBF(cos) + CBF.gen(0) * CBF(sin)
    xN = x

    # Suppose that (cos, sin) are indeed sine and cosine of a rational angle.
    # Then x = cos + I*sin generates a cyclotomic field C and for some N we
    # have x^N = ±1. Since C is contained in the compositum of K=Q(cos) and
    # L=Q(i*sin) and Q(cos) and Q(sin) are both contained in C, the degree of C
    # is bounded from above by twice (accounting for the imaginary unit) the
    # degrees of K and L. The degree of C is the totient of N which is bounded
    # from below by n / (e^γ loglog n + 3 / loglog n) [cf. wikipedia].
    degree_bound = 2 * cos.minpoly().degree() * sin.minpoly().degree()

    from itertools import count

    for n in count(2):
        xN *= x

        c = xN.real()
        s = xN.imag()

        if xN.real().contains_zero() or xN.imag().contains_zero():
            c, s = cos, sin
            for i in range(n - 1):
                c, s = c * cos - s * sin, s * cos + c * sin

            if c == 0 or s == 0:
                return True

            CBF = ComplexBallField(CBF.precision() * 2)
            x = CBF(cos) + CBF.gen(0) * CBF(sin)
            xN = x**n

        from math import log

        if n / (2.0 * log(log(n)) + 3 / log(log(n))) > 2 * degree_bound:
            return False


def acos(cos_angle, numerical=False):
    r"""
    Return the arccosine of ``cos_angle`` as a multiple of 2π, i.e., as a value
    between 0 and 1/2.

    INPUT:

    - ``cos_angle`` -- a floating point number, the cosine of an angle

    - ``numerical`` -- a boolean (default: ``False``); whether to return a
      numerical approximation of the arccosine or try to reconstruct an exact
      rational value for the arccosine (in radians.)

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import acos

        sage: acos(1)
        0
        sage: acos(.5)
        1/6
        sage: acos(0)
        1/4
        sage: acos(-.5)
        1/3
        sage: acos(-1)
        1/2

        sage: acos(.25)
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot recover a rational angle from these numerical results
        sage: acos(.25, numerical=True)
        0.2097846883724169

    """
    import math

    angle = math.acos(cos_angle) / (2 * math.pi)

    assert 0 <= angle <= 0.5

    if numerical:
        return angle

    # fast and dirty way using floating point approximation
    from sage.all import RR

    angle_rat = RR(angle).nearby_rational(0.00000001)
    if angle_rat.denominator() > 256:
        raise NotImplementedError(
            "cannot recover a rational angle from these numerical results"
        )
    return angle_rat


def angle(u, v, numerical=False):
    r"""
    Return the angle between the vectors ``u`` and ``v`` divided by `2 \pi`.

    INPUT:

    - ``u``, ``v`` - vectors

    - ``numerical`` - boolean (default: ``False``), whether to return floating
      point numbers

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import angle

    As the implementation is dirty, we at least check that it works for all
    denominator up to 20::

        sage: u = vector((AA(1),AA(0)))
        sage: for n in xsrange(1,20):       # long time  (1.5s)
        ....:     for k in xsrange(1,n):
        ....:         v = vector((AA(cos(2*k*pi/n)), AA(sin(2*k*pi/n))))
        ....:         assert angle(u,v) == k/n

    The numerical version (working over floating point numbers)::

        sage: import math
        sage: u = (1, 0)
        sage: for n in xsrange(1,20):
        ....:     for k in xsrange(1,n):
        ....:         a = 2 * k * math.pi / n
        ....:         v = (math.cos(a), math.sin(a))
        ....:         assert abs(angle(u,v,numerical=True) * 2 * math.pi - a) < 1.e-10

    If the angle is not rational, then the method returns an element in the real
    lazy field::

        sage: v = vector((AA(sqrt(2)), AA(sqrt(3))))
        sage: a = angle(u, v)
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot recover a rational angle from these numerical results
        sage: a = angle(u, v, numerical=True)
        sage: a    # abs tol 1e-14
        0.14102355421224375
        sage: exp(2*pi.n()*CC(0,1)*a)
        0.632455532033676 + 0.774596669241483*I
        sage: v / v.norm()
        (0.6324555320336758?, 0.774596669241484?)

    """
    import math

    u0 = float(u[0])
    u1 = float(u[1])
    v0 = float(v[0])
    v1 = float(v[1])

    cos_uv = (u0 * v0 + u1 * v1) / math.sqrt((u0 * u0 + u1 * u1) * (v0 * v0 + v1 * v1))
    if cos_uv < -1.0:
        assert cos_uv > -1.0000001
        cos_uv = -1.0
    elif cos_uv > 1.0:
        assert cos_uv < 1.0000001
        cos_uv = 1.0

    angle = acos(cos_uv, numerical=numerical)
    return 1 - angle if u0 * v1 - u1 * v0 < 0 else angle

    # a neater way is provided below by working only with number fields
    # but this method is slower...
    # sqnorm_u = u[0]*u[0] + u[1]*u[1]
    # sqnorm_v = v[0]*v[0] + v[1]*v[1]
    #
    # if sqnorm_u != sqnorm_v:
    #    # we need to take a square root in order that u and v have the
    #    # same norm
    #    u = (1 / AA(sqnorm_u)).sqrt() * u.change_ring(AA)
    #    v = (1 / AA(sqnorm_v)).sqrt() * v.change_ring(AA)
    #    sqnorm_u = AA.one()
    #    sqnorm_v = AA.one()
    #
    # cos_uv = (u[0]*v[0] + u[1]*v[1]) / sqnorm_u
    # sin_uv = (u[0]*v[1] - u[1]*v[0]) / sqnorm_u


def ccw(v, w):
    r"""
    Return a positive number if the turn from ``v`` to ``w`` is
    counterclockwise, a negative number if it is clockwise, and zero if the two
    vectors are collinear.

    .. NOTE::

        This function is sometimes also referred to as the wedge product or
        simply the determinant. We chose the more customary name ``ccw`` from
        computational geometry here.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import ccw
        sage: ccw((1, 0), (0, 1))
        1
        sage: ccw((1, 0), (-1, 0))
        0
        sage: ccw((1, 0), (0, -1))
        -1
        sage: ccw((1, 0), (1, 0))
        0

    """
    return v[0] * w[1] - v[1] * w[0]


def is_parallel(v, w):
    r"""
    Return whether the vectors ``v`` and ``w`` are parallel (but not
    anti-parallel.)

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_parallel
        sage: is_parallel((0, 1), (0, 1))
        True
        sage: is_parallel((0, 1), (0, 2))
        True
        sage: is_parallel((0, 1), (0, -2))
        False
        sage: is_parallel((0, 1), (0, 0))
        False
        sage: is_parallel((0, 1), (1, 0))
        False

    TESTS::

        sage: V = QQ**2

        sage: is_parallel(V((0,1)), V((0,2)))
        True
        sage: is_parallel(V((1,-1)), V((2,-2)))
        True
        sage: is_parallel(V((4,-2)), V((2,-1)))
        True
        sage: is_parallel(V((1,2)), V((2,4)))
        True
        sage: is_parallel(V((0,2)), V((0,1)))
        True

        sage: is_parallel(V((1,1)), V((1,2)))
        False
        sage: is_parallel(V((1,2)), V((2,1)))
        False
        sage: is_parallel(V((1,2)), V((1,-2)))
        False
        sage: is_parallel(V((1,2)), V((-1,-2)))
        False
        sage: is_parallel(V((2,-1)), V((-2,1)))
        False

    """
    if ccw(v, w) != 0:
        # vectors are not collinear
        return False

    return v[0] * w[0] + v[1] * w[1] > 0


def is_anti_parallel(v, w):
    r"""
    Return whether the vectors ``v`` and ``w`` are anti-parallel, i.e., whether
    ``v`` and ``-w`` are parallel.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_anti_parallel
        sage: V = QQ**2

        sage: is_anti_parallel(V((0,1)), V((0,-2)))
        True
        sage: is_anti_parallel(V((1,-1)), V((-2,2)))
        True
        sage: is_anti_parallel(V((4,-2)), V((-2,1)))
        True
        sage: is_anti_parallel(V((-1,-2)), V((2,4)))
        True

        sage: is_anti_parallel(V((1,1)), V((1,2)))
        False
        sage: is_anti_parallel(V((1,2)), V((2,1)))
        False
        sage: is_anti_parallel(V((0,2)), V((0,1)))
        False
        sage: is_anti_parallel(V((1,2)), V((1,-2)))
        False
        sage: is_anti_parallel(V((1,2)), V((-1,2)))
        False
        sage: is_anti_parallel(V((2,-1)), V((-2,-1)))
        False

    """
    return is_parallel(v, -w)


def line_intersection(p1, p2, q1, q2):
    r"""
    Return the point of intersection between the line joining p1 to p2
    and the line joining q1 to q2. If the lines do not have a single point of
    intersection, we return None. Here p1, p2, q1 and q2 should be vectors in
    the plane.
    """
    if ccw(p2 - p1, q2 - q1) == 0:
        return None

    # Since the wedge product is non-zero, the following is invertible:
    from sage.all import matrix

    m = matrix([[p2[0] - p1[0], q1[0] - q2[0]], [p2[1] - p1[1], q1[1] - q2[1]]])
    return p1 + (m.inverse() * (q1 - p1))[0] * (p2 - p1)


def is_segment_intersecting(e1, e2):
    r"""
    Return whether the segments ``e1`` and ``e2`` intersect.

    OUTPUT:

    - ``0`` - do not intersect
    - ``1`` - one endpoint in common
    - ``2`` - non-trivial intersection

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_segment_intersecting
        sage: is_segment_intersecting(((0,0),(1,0)),((0,1),(0,3)))
        0
        sage: is_segment_intersecting(((0,0),(1,0)),((0,0),(0,3)))
        1
        sage: is_segment_intersecting(((0,0),(1,0)),((0,-1),(0,3)))
        2
        sage: is_segment_intersecting(((-1,-1),(1,1)),((0,0),(2,2)))
        2
        sage: is_segment_intersecting(((-1,-1),(1,1)),((1,1),(2,2)))
        1

    """
    if e1[0] == e1[1] or e2[0] == e2[1]:
        raise ValueError("degenerate segments")

    elts = [e[i][j] for e in (e1, e2) for i in (0, 1) for j in (0, 1)]

    from sage.structure.element import get_coercion_model

    cm = get_coercion_model()

    base_ring = cm.common_parent(*elts)
    if isinstance(base_ring, type):
        from sage.structure.coerce import py_scalar_parent

        base_ring = py_scalar_parent(base_ring)

    from sage.all import matrix

    m = matrix(base_ring, 3)
    xs1, ys1 = map(base_ring, e1[0])
    xt1, yt1 = map(base_ring, e1[1])
    xs2, ys2 = map(base_ring, e2[0])
    xt2, yt2 = map(base_ring, e2[1])

    m[0] = [xs1, ys1, 1]
    m[1] = [xt1, yt1, 1]
    m[2] = [xs2, ys2, 1]
    s0 = m.det()
    m[2] = [xt2, yt2, 1]
    s1 = m.det()
    if (s0 > 0 and s1 > 0) or (s0 < 0 and s1 < 0):
        # e2 stands on one side of the line generated by e1
        return 0

    m[0] = [xs2, ys2, 1]
    m[1] = [xt2, yt2, 1]
    m[2] = [xs1, ys1, 1]
    s2 = m.det()
    m[2] = [xt1, yt1, 1]
    s3 = m.det()
    if (s2 > 0 and s3 > 0) or (s2 < 0 and s3 < 0):
        # e1 stands on one side of the line generated by e2
        return 0

    if s0 == 0 and s1 == 0:
        assert s2 == 0 and s3 == 0
        if xt1 < xs1 or (xt1 == xs1 and yt1 < ys1):
            xs1, xt1 = xt1, xs1
            ys1, yt1 = yt1, ys1
        if xt2 < xs2 or (xt2 == xs2 and yt2 < ys2):
            xs2, xt2 = xt2, xs2
            ys2, yt2 = yt2, ys2

        if xs1 == xt1 == xs2 == xt2:
            xs1, xt1, xs2, xt2 = ys1, yt1, ys2, yt2

        assert xs1 < xt1 and xs2 < xt2, (xs1, xt1, xs2, xt2)

        if (xs2 > xt1) or (xt2 < xs1):
            return 0  # no intersection
        elif (xs2 == xt1) or (xt2 == xs1):
            return 1  # one endpoint in common
        else:
            assert (
                xs1 <= xs2 < xt1
                or xs1 < xt2 <= xt1
                or (xs2 < xs1 and xt2 > xt1)
                or (xs2 > xs1 and xt2 < xt1)
            ), (xs1, xt1, xs2, xt2)
            return 2  # one dimensional

    elif s0 == 0 or s1 == 0:
        # treat alignment here
        if s2 == 0 or s3 == 0:
            return 1  # one endpoint in common
        else:
            return 2  # intersection in the middle

    return 2  # middle intersection


def is_between(e0, e1, f):
    r"""
    Check whether the vector ``f`` is strictly in the sector formed by the vectors
    ``e0`` and ``e1`` (in counter-clockwise order).

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_between
        sage: is_between((1, 0), (1, 1), (2, 1))
        True

        sage: from itertools import product
        sage: vecs = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
        sage: for (i, vi), (j, vj), (k, vk) in product(enumerate(vecs), repeat=3):
        ....:     assert is_between(vi, vj, vk) == ((i == j and i != k) or i < k < j or k < j < i or j < i < k), ((i, vi), (j, vj), (k, vk))
    """
    if e0[0] * e1[1] > e1[0] * e0[1]:
        # positive determinant
        # [ e0[0] e1[0] ]^-1 = [ e1[1] -e1[0] ]
        # [ e0[1] e1[1] ]      [-e0[1]  e0[0] ]
        # f[0] * e1[1] - e1[0] * f[1] > 0
        # - f[0] * e0[1] + e0[0] * f[1] > 0
        return e1[1] * f[0] > e1[0] * f[1] and e0[0] * f[1] > e0[1] * f[0]
    elif e0[0] * e1[1] == e1[0] * e0[1]:
        # aligned vector
        if e0[0] * e1[0] >= 0 and e0[1] * e1[1] >= 0:
            return f[0] * e0[1] != f[1] * e0[0] or f[0] * e0[0] < 0 or f[1] * e0[1] < 0
        else:
            return e0[0] * f[1] > e0[1] * f[0]
    else:
        # negative determinant
        # [ e1[0] e0[0] ]^-1 = [ e0[1] -e0[0] ]
        # [ e1[1] e0[1] ]      [-e1[1]  e1[0] ]
        # f[0] * e0[1] - e0[0] * f[1] > 0
        # - f[0] * e1[1] + e1[0] * f[1] > 0
        return e0[1] * f[0] < e0[0] * f[1] or e1[0] * f[1] < e1[1] * f[0]


def solve(x, u, y, v):
    r"""
    Return (a,b) so that: x + au = y + bv

    INPUT:

    - ``x``, ``u``, ``y``, ``v`` -- two dimensional vectors

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import solve
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: V = VectorSpace(K,2)
        sage: x = V((1,-sqrt2))
        sage: y = V((1,1))
        sage: a = V((0,1))
        sage: b = V((-sqrt2, sqrt2+1))
        sage: u = V((0,1))
        sage: v = V((-sqrt2, sqrt2+1))
        sage: a, b = solve(x,u,y,v)
        sage: x + a*u == y + b*v
        True

        sage: u = V((1,1))
        sage: v = V((1,sqrt2))
        sage: a, b = solve(x,u,y,v)
        sage: x + a*u == y + b*v
        True

    """
    d = -u[0] * v[1] + u[1] * v[0]
    if d.is_zero():
        raise ValueError("parallel vectors")
    a = v[1] * (x[0] - y[0]) + v[0] * (y[1] - x[1])
    b = u[1] * (x[0] - y[0]) + u[0] * (y[1] - x[1])
    return (a / d, b / d)


def projectivization(x, y, signed=True, denominator=None):
    r"""
    Return a simplified version of the projective coordinate [x: y].

    If ``signed`` (the default), the second coordinate is made non-negative;
    otherwise the coordinates keep their signs.

    If ``denominator`` is ``False``, returns [x/y: 1] up to sign. Otherwise,
    the returned coordinates have no denominator and no non-unit gcd.

    TESTS::

        sage: from flatsurf.geometry.euclidean import projectivization

        sage: projectivization(2/3, -3/5, signed=True, denominator=True)
        (10, -9)
        sage: projectivization(2/3, -3/5, signed=False, denominator=True)
        (-10, 9)
        sage: projectivization(2/3, -3/5, signed=True, denominator=False)
        (10/9, -1)
        sage: projectivization(2/3, -3/5, signed=False, denominator=False)
        (-10/9, 1)

        sage: projectivization(-1/2, 0, signed=True, denominator=True)
        (-1, 0)
        sage: projectivization(-1/2, 0, signed=False, denominator=True)
        (1, 0)
        sage: projectivization(-1/2, 0, signed=True, denominator=False)
        (-1, 0)
        sage: projectivization(-1/2, 0, signed=False, denominator=False)
        (1, 0)

    """
    from sage.all import Sequence

    parent = Sequence([x, y]).universe()
    if y:
        z = x / y
        if denominator is True or (denominator is None and hasattr(z, "denominator")):
            d = parent(z.denominator())
        else:
            d = parent(1)
        if signed and y < 0:
            d *= -1
        return (z * d, d)
    elif signed and x < 0:
        return (parent(-1), parent(0))
    else:
        return (parent(1), parent(0))


def slope(a, rotate=1):
    r"""
    Return either ``1`` (positive slope) or ``-1`` (negative slope).

    If ``rotate`` is set to 1 then consider the edge as if it was rotated counterclockwise
    infinitesimally.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import slope
        sage: slope((1, 1))
        1
        sage: slope((-1, 1))
        -1
        sage: slope((-1, -1))
        1
        sage: slope((1, -1))
        -1

        sage: slope((1, 0))
        1
        sage: slope((0, 1))
        -1
        sage: slope((-1, 0))
        1
        sage: slope((0, -1))
        -1

        sage: slope((1, 0), rotate=-1)
        -1
        sage: slope((0, 1), rotate=-1)
        1
        sage: slope((-1, 0), rotate=-1)
        -1
        sage: slope((0, -1), rotate=-1)
        1

        sage: slope((1, 0), rotate=0)
        0
        sage: slope((0, 1), rotate=0)
        0
        sage: slope((-1, 0), rotate=0)
        0
        sage: slope((0, -1), rotate=0)
        0

        sage: slope((0, 0))
        Traceback (most recent call last):
        ...
        ValueError: zero vector
    """
    x, y = a
    if not x and not y:
        raise ValueError("zero vector")
    if (x > 0 and y > 0) or (x < 0 and y < 0):
        return 1
    elif (x > 0 and y < 0) or (x < 0 and y > 0):
        return -1
    if rotate == 0:
        return 0
    if rotate == 1:
        return 1 if x else -1
    if rotate == -1:
        return 1 if y else -1
    raise ValueError("invalid argument rotate={}".format(rotate))


from sage.all import UniqueRepresentation, Parent
from sage.structure.element import Element
from sage.misc.cachefunc import cached_method


class Cone(Element):
    # This is an open cone.
    def __init__(self, parent, start, end):
        super().__init__(parent)
        self._start = parent.rays()(start)
        self._end = parent.rays()(end)

    @cached_method
    def is_empty(self):
        return self._start == self._end

    def is_convex(self):
        from flatsurf.geometry.euclidean import ccw

        return ccw(self._start.vector(), self._end.vector()) >= 0

    # TODO: Rename to contains_cone() [arguments are reversed!]
    def is_subset(self, other):
        r"""
        Return whether this cone is contained in ``other``.
        """
        if not isinstance(other, Cone):
            raise NotImplementedError

        if other.parent() is not self.parent():
            raise NotImplementedError

        if self.is_empty():
            return True

        if other.is_empty():
            return False

        from flatsurf.geometry.euclidean import ccw

        start_ccw = ccw(self._start.vector(), other._start.vector())
        end_ccw = ccw(self._end.vector(), other._end.vector())

        # Check whether the interior of this cone is contained in other.
        if not self.is_convex():
            if other.is_convex():
                return False

            # This non-convex cone C is contained in the other non-convex cone D,
            # iff for the complements we have D^c ⊆ C^c.
            return other.complement().is_subset(self.complement())

        if other.is_convex():
            if start_ccw > 0:
                return False
            if end_ccw < 0:
                return False

            return True

        raise NotImplementedError

    def complement(self):
        r"""
        Return the maximal cone contained in the complement of this cone, i.e.,
        the complement of this cone with the boundaries (except for the origin)
        missing.

        """
        if self.is_empty():
            raise NotImplementedError

        return self.parent()(self._end, self._start)

    def contains_ray(self, ray):
        if ray.parent() is not self.parent().rays():
            raise NotImplementedError

        if self.is_empty():
            return False

        from flatsurf.geometry.euclidean import ccw

        ccw_from_start = ccw(self._start.vector(), ray.vector())
        ccw_to_end = ccw(ray.vector(), self._end.vector())

        if self.is_convex():
            if ccw_from_start <= 0:
                return False

            if ccw_to_end <= 0:
                return False

            return True

        return (
            not self.complement().contains_ray(ray)
            and ray != self._start
            and ray != self._end
        )

    def sorted_rays(self, rays):
        class Key:
            def __init__(self, cone, ray):
                self._cone = cone
                self._ray = ray

            def __lt__(self, rhs):
                # TODO: Make sure all code paths are tested.
                from flatsurf.geometry.euclidean import ccw

                if not self._cone.is_convex():
                    start_to_self = ccw(self._cone._start.vector(), self._ray.vector())
                    start_to_rhs = ccw(self._cone._start.vector(), rhs._ray.vector())
                    if start_to_self > 0 and start_to_rhs > 0:
                        return ccw(self._ray.vector(), rhs._ray.vector()) > 0

                    end_to_self = ccw(self._cone._end.vector(), self._ray.vector())
                    end_to_rhs = ccw(self._cone._end.vector(), rhs._ray.vector())
                    if end_to_self < 0 and end_to_rhs < 0:
                        return ccw(self._ray.vector(), rhs._ray.vector()) > 0

                    if start_to_self > 0 and end_to_rhs < 0:
                        return True
                    if end_to_self < 0 and start_to_rhs > 0:
                        return False

                    raise NotImplementedError
                return ccw(self._ray.vector(), rhs._ray.vector()) > 0

        rays = sorted(rays, key=lambda ray: Key(self, ray))

        from itertools import groupby

        return [ray for ray, _ in groupby(rays)]

    def a_ray(self):
        if self.is_empty():
            raise TypeError

        if self.is_convex():
            return self.parent().rays()((self._start.vector() + self._end.vector()) / 2)

        raise NotImplementedError

    def start(self):
        return self._start

    def end(self):
        return self._end

    def contains_point(self, p):
        if self.is_empty():
            return False
        if p.is_zero():
            return True
        return self.contains_ray(self.parent().rays()(p))

    def _repr_(self):
        if self.is_empty():
            return "Empty cone"
        return f"Open cone between {self._start} and {self._end}"


class Cones(UniqueRepresentation, Parent):
    Element = Cone

    def __init__(self, base_ring, category=None):
        from sage.categories.all import Sets

        super().__init__(base_ring, category=category or Sets())

    @cached_method
    def rays(self):
        from flatsurf.geometry.ray import Rays

        return Rays(self.base_ring())


r"""
Geometry with rays in the Euclidean plane.

EXAMPLES::

    sage: from flatsurf.geometry.ray import Rays
    sage: R = Rays(QQ)
    sage: R((1, 0))
    Ray towards (1, 0)

    sage: R((1, 0)) == R((2, 0))
    True

    sage: R((1, 0)) == R((-2, 0))
    False

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2024 Julian Rüth
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
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method


class Ray(Element):
    r"""
    A ray in the Euclidean plane.

    EXAMPLES::

        sage: from flatsurf.geometry.ray import Rays
        sage: R = Rays(QQ)
        sage: r = R((1, 0)); r
        Ray towards (1, 0)

        sage: R((0, 0))
        Traceback (most recent call last):
        ...
        ValueError: direction must not be the zero vector

    TESTS::

        sage: from flatsurf.geometry.ray import Ray
        sage: isinstance(r, Ray)
        True
        sage: TestSuite(r).run()

    """

    def __init__(self, parent, direction):
        super().__init__(parent)

        direction = parent.ambient_space()(direction)
        direction.set_immutable()

        if direction.is_zero():
            raise ValueError("direction must not be the zero vector")

        self._direction = direction

    def vector(self):
        r"""
        Return a vector in this ray.

        EXAMPLES::

            sage: from flatsurf.geometry.ray import Rays
            sage: R = Rays(QQ)
            sage: r = R((1, 0)); r
            Ray towards (1, 0)
            sage: r.vector()
            (1, 0)

        """
        return self._direction

    def _neg_(self):
        return self.parent()(-self._direction)

    def _repr_(self):
        return f"Ray towards {self._direction}"

    def _richcmp_(self, other, op):
        r"""
        Return how this ray compares to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.ray import Rays
            sage: R = Rays(QQ)

            sage: R((0, 1)) == R((0, 2))
            True
            sage: R((1, 0)) == R((2, 0))
            True
            sage: R((1, 1)) == R((2, 2))
            True
            sage: R((1, 1)) == R((2, 1))
            False
            sage: R((0, 1)) == R((0, -2))
            False
            sage: R((1, 1)) == R((-2, -2))
            False

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            from sage.all import sgn

            return (
                self._direction[0] * other._direction[1]
                == other._direction[0] * self._direction[1]
                and sgn(self._direction[0]) == sgn(other._direction[0])
                and sgn(self._direction[1]) == sgn(other._direction[1])
            )


class Rays(UniqueRepresentation, Parent):
    r"""
    The space of rays from the origin in the Euclidean plane.

    EXAMPLES::

        sage: from flatsurf.geometry.ray import Rays
        sage: R = Rays(QQ)
        sage: R
        Rays in Vector space of dimension 2 over Rational Field

    TESTS::

        sage: isinstance(R, Rays)
        True
        sage: TestSuite(R).run()

    """
    Element = Ray

    def __init__(self, base_ring, category=None):
        from sage.categories.all import Sets, Rings

        if base_ring not in Rings():
            raise TypeError("base ring must be a ring")

        super().__init__(base=base_ring, category=category or Sets())

    def _an_element_(self):
        return self((1, 0))

    def some_element(self):
        return [self(v) for v in self.ambient_space().some_elements() if v]

    @cached_method
    def ambient_space(self):
        r"""
        Return the ambient Euclidean space containing these rays.

        EXAMPLES::

            sage: from flatsurf.geometry.ray import Rays
            sage: Rays(QQ).ambient_space()
            Vector space of dimension 2 over Rational Field

        """
        return self.base_ring() ** 2

    def _repr_(self):
        return f"Rays in {self.ambient_space()}"


def time_on_ray(p, direction, q):
    if direction[0]:
        dim = 0
    else:
        dim = 1

    delta = q[dim] - p[dim]
    length = direction[dim]
    if length < 0:
        delta *= -1
        length *= -1

    return delta, length


def ray_segment_intersection(p, direction, segment):
    r"""
    Return the intersection of the ray from ``p`` in ``direction`` with
    ``segment``.

    If the segment and the ray intersect in a point, return that point as a
    vector.

    If the segment and the ray overlap in a segment, return the end points of
    that segment (in order.)

    If the segment and the ray do not intersect, return ``None``.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import ray_segment_intersection
        sage: V = QQ**2

        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((1, -1)), V((1, 1))))
        (1, 0)
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((0, 0)), V((1, 0))))
        ((0, 0), (1, 0))
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((1, 0)), V((2, 0))))
        ((1, 0), (2, 0))
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((-1, 0)), V((1, 0))))
        ((0, 0), (1, 0))
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((-1, -1)), V((-1, 1))))

    TESTS::

        sage: ray_segment_intersection(V((0, 0)), V((5, 1)), (V((3, 2)), V((3, 3))))

    """
    intersection = line_intersection(*(p, p + direction), *segment)

    if intersection is None:
        # ray and segment are parallel.
        if ccw(direction, segment[0] - p) != 0:
            # ray and segment are not on the same line
            return None

        t0, length = time_on_ray(p, direction, segment[0])
        t1, _ = time_on_ray(p, direction, segment[1])

        if t1 < t0:
            t0, t1 = t1, t0

        if t1 < 0:
            return None

        if t1 == 0:
            return p

        t0 /= length
        t1 /= length

        if t0 < 0:
            return (p, p + t1 * direction)

        return (p + t0 * direction, p + t1 * direction)

    if time_on_ray(p, direction, intersection)[0] < 0:
        return None

    if ccw(segment[0] - p, direction) * ccw(segment[1] - p, direction) > 0:
        return None

    return intersection


def circle_from_three_points(p, q, r, base_ring=None):
    r"""
    Construct a circle from three points on the circle.
    """
    if base_ring is None:
        base_ring = p.base_ring()
    V2 = VectorSpace(base_ring.fraction_field(), 2)
    V3 = VectorSpace(base_ring.fraction_field(), 3)

    v1 = V3((p[0] + q[0], p[1] + q[1], 2))
    v2 = V3((p[1] - q[1], q[0] - p[0], 0))
    line1 = v1.cross_product(v2)
    v1 = V3((p[0] + r[0], p[1] + r[1], 2))
    v2 = V3((p[1] - r[1], r[0] - p[0], 0))
    line2 = v1.cross_product(v2)
    center_3 = line1.cross_product(line2)
    if center_3[2].is_zero():
        raise ValueError("The three points lie on a line.")
    center = V2((center_3[0] / center_3[2], center_3[1] / center_3[2]))
    return Circle(center, (p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)


class Circle:
    def __init__(self, center, radius_squared, base_ring=None):
        r"""
        Construct a circle from a Vector representing the center, and the
        radius squared.
        """
        if base_ring is None:
            self._base_ring = radius_squared.parent()
        else:
            self._base_ring = base_ring

        # for calculations:
        self._V2 = VectorSpace(self._base_ring, 2)
        self._V3 = VectorSpace(self._base_ring, 3)

        self._center = self._V2(center)
        self._center.set_immutable()
        self._radius_squared = self._base_ring(radius_squared)

    def center(self):
        r"""
        Return the center of the circle as a vector.
        """
        return self._center

    def radius_squared(self):
        r"""
        Return the square of the radius of the circle.
        """
        return self._radius_squared

    def point_position(self, point):
        r"""
        Return 1 if point lies in the circle, 0 if the point lies on the circle,
        and -1 if the point lies outide the circle.
        """
        value = (
            (point[0] - self._center[0]) ** 2
            + (point[1] - self._center[1]) ** 2
            - self._radius_squared
        )
        if value > self._base_ring.zero():
            return -1
        if value < self._base_ring.zero():
            return 1
        return 0

    def closest_point_on_line(self, point, direction_vector):
        r"""
        Consider the line through the provided point in the given direction.
        Return the closest point on this line to the center of the circle.
        """
        cc = self._V3((self._center[0], self._center[1], self._base_ring.one()))
        # point at infinite orthogonal to direction_vector:
        dd = self._V3(
            (direction_vector[1], -direction_vector[0], self._base_ring.zero())
        )
        l1 = cc.cross_product(dd)

        pp = self._V3((point[0], point[1], self._base_ring.one()))
        # direction_vector pushed to infinity
        ee = self._V3(
            (direction_vector[0], direction_vector[1], self._base_ring.zero())
        )
        l2 = pp.cross_product(ee)

        # This is the point we want to return
        rr = l1.cross_product(l2)
        try:
            return self._V2((rr[0] / rr[2], rr[1] / rr[2]))
        except ZeroDivisionError:
            raise ValueError(
                "Division by zero error. Perhaps direction is zero. "
                + "point="
                + str(point)
                + " direction="
                + str(direction_vector)
                + " circle="
                + str(self)
            )

    def line_position(self, point, direction_vector):
        r"""
        Consider the line through the provided point in the given direction.
        We return 1 if the line passes through the circle, 0 if it is tangent
        to the circle and -1 if the line does not intersect the circle.
        """
        return self.point_position(self.closest_point_on_line(point, direction_vector))

    def line_segment_position(self, p, q):
        r"""
        Consider the open line segment pq.We return 1 if the line segment
        enters the interior of the circle, zero if it touches the circle
        tangentially (at a point in the interior of the segment) and
        and -1 if it does not touch the circle or its interior.
        """
        if self.point_position(p) == 1:
            return 1
        if self.point_position(q) == 1:
            return 1
        r = self.closest_point_on_line(p, q - p)
        pos = self.point_position(r)
        if pos == -1:
            return -1
        # This checks if r lies in the interior of pq
        if p[0] == q[0]:
            if (p[1] < r[1] and r[1] < q[1]) or (p[1] > r[1] and r[1] > q[1]):
                return pos
        elif (p[0] < r[0] and r[0] < q[0]) or (p[0] > r[0] and r[0] > q[0]):
            return pos
        # It does not lie in the interior.
        return -1

    def tangent_vector(self, point):
        r"""
        Return a vector based at the provided point (which must lie on the circle)
        which is tangent to the circle and points in the counter-clockwise
        direction.

        EXAMPLES::

            sage: from flatsurf.geometry.circle import Circle
            sage: c=Circle(vector((0,0)), 2, base_ring=QQ)
            sage: c.tangent_vector(vector((1,1)))
            (-1, 1)
        """
        if not self.point_position(point) == 0:
            raise ValueError("point not on circle.")
        return vector((self._center[1] - point[1], point[0] - self._center[0]))

    def other_intersection(self, p, v):
        r"""
        Consider a point p on the circle and a vector v. Let L be the line
        through p in direction v. Then L intersects the circle at another
        point q. This method returns q.

        Note that if p and v are both in the field of the circle,
        then so is q.

        EXAMPLES::

            sage: from flatsurf.geometry.circle import Circle
            sage: c=Circle(vector((0,0)), 25, base_ring=QQ)
            sage: c.other_intersection(vector((3,4)),vector((1,2)))
            (-7/5, -24/5)
        """
        pp = self._V3((p[0], p[1], self._base_ring.one()))
        vv = self._V3((v[0], v[1], self._base_ring.zero()))
        L = pp.cross_product(vv)
        cc = self._V3((self._center[0], self._center[1], self._base_ring.one()))
        vvperp = self._V3((-v[1], v[0], self._base_ring.zero()))
        # line perpendicular to L through center:
        Lperp = cc.cross_product(vvperp)
        # intersection of L and Lperp:
        rr = L.cross_product(Lperp)
        r = self._V2((rr[0] / rr[2], rr[1] / rr[2]))
        return self._V2((2 * r[0] - p[0], 2 * r[1] - p[1]))

    def __rmul__(self, similarity):
        r"""
        Apply a similarity to the circle.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.square_torus()
            sage: c = s.polygon(0).circumscribed_circle()
            sage: c
            Circle((1/2, 1/2), 1/2)
            sage: s.edge_transformation(0,2)
            (x, y) |-> (x, y - 1)
            sage: s.edge_transformation(0,2) * c
            Circle((1/2, -1/2), 1/2)
        """
        from .similarity import SimilarityGroup

        SG = SimilarityGroup(self._base_ring)
        s = SG(similarity)
        return Circle(
            s(self._center), s.det() * self._radius_squared, base_ring=self._base_ring
        )

    def __str__(self):
        return (
            "circle with center "
            + str(self._center)
            + " and radius squared "
            + str(self._radius_squared)
        )

    def __repr__(self):
        return "Circle(" + repr(self._center) + ", " + repr(self._radius_squared) + ")"

    def __hash__(self):
        r"""
        Return a hash value for this circle that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain().relabel()
            sage: hash(S.polygon(0).circumscribed_circle()) == hash(S.polygon(1).circumscribed_circle())
            True

        """
        return hash((self._center, self._radius_squared))

    def __eq__(self, other):
        r"""
        Return whether this circle is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain().relabel()
            sage: S.polygon(0).circumscribed_circle() == S.polygon(1).circumscribed_circle()
            True

        """
        if not isinstance(other, Circle):
            return False

        return (
            self._center == other._center
            and self._radius_squared == other._radius_squared
        )
