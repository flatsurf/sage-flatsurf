r"""
A loose collection of tools for Euclidean geometry in the plane.

.. SEEALSO::

    :mod:`flatsurf.geometry.circle` for everything specific to circles in the plane

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2024 Julian Rüth
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

        sage: from pyexactreal import ExactReals # optional: exactreal  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: R = ExactReals() # optional: exactreal
        sage: is_cosine_sine_of_rational(R.one(), R.zero()) # optional: exactreal
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


def angle(u, v, numerical=False):
    r"""
    Return the angle between the vectors ``u`` and ``v`` divided by `2 \pi`.

    INPUT:

    - ``u``, ``v`` - vectors in the plane

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
    angle = math.acos(cos_uv) / (2 * math.pi)  # rat number between 0 and 1/2

    if numerical:
        return 1.0 - angle if u0 * v1 - u1 * v0 < 0 else angle

    # fast and dirty way using floating point approximation
    # (see below for a slow but exact method)
    from sage.all import RR

    angle_rat = RR(angle).nearby_rational(0.00000001)
    if angle_rat.denominator() > 256:
        raise NotImplementedError(
            "cannot recover a rational angle from these numerical results"
        )
    return 1 - angle_rat if u0 * v1 - u1 * v0 < 0 else angle_rat


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


def line_intersection(l, m):
    r"""
    Return the point of intersection between the lines ``l`` and ``m``. If the
    lines do not have a single point of intersection, returns None.

    INPUT:

    - ``l`` -- a line in the plane given by two points (as vectors) in the plane

    - ``m`` -- a line in the plane given by two points (as vectors) in the plane

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import line_intersection
        sage: line_intersection((vector((-1, 0)), vector((1, 0))), (vector((0, -1)), vector((0, 1))))
        (0, 0)

    Parallel lines have no single point of intersection:: 

        sage: line_intersection((vector((-1, 0)), vector((1, 0))), (vector((-1, 1)), vector((1, 1))))

    Identical lines have no single point of intersection::

        sage: line_intersection((vector((-1, 0)), vector((1, 0))), (vector((-2, 0)), vector((2, 0))))

    """
    Δl = l[1] - l[0]
    Δm = m[1] - m[0]

    # We solve the linear system determining the time when l[0] + t * Δl hits
    # m, i.e., l[0] + t Δl == m[0] + s Δm.
    a, b, c, d = Δl[0], -Δm[0], Δl[1], -Δm[1]
    rhs = m[0] - l[0]

    det = a * d - b * c
    if det == 0:
        # The lines are parallel
        return None

    # Solve the linear system. We only need t (and not s.)
    t = (d * rhs[0] - b * rhs[1]) / det

    return l[0] + t * Δl


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
    intersection = line_intersection((p, p + direction), segment)

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


def is_segment_intersecting(s, t):
    r"""
    Return whether the segments ``s`` and ``t`` intersect.

    INPUT:

    - ``s`` -- a segment given as a pair of endpoints (given as vectors in the plane.)

    - ``t`` -- a segment given as a pair of endpoints (given as vectors in the plane.)

    OUTPUT:

    - ``0`` - do not intersect
    - ``1`` - exactly one endpoint in common
    - ``2`` - non-trivial intersection

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_segment_intersecting
        sage: is_segment_intersecting((vector((0, 0)), vector((1, 0))), (vector((0, 1)), vector((0, 3))))
        0
        sage: is_segment_intersecting((vector((0, 0)), vector((1, 0))), (vector((0, 0)), vector((0, 3))))
        1
        sage: is_segment_intersecting((vector((0, 0)), vector((1, 0))), (vector((0, -1)), vector((0, 3))))
        2
        sage: is_segment_intersecting((vector((-1, -1)), vector((1, 1))), (vector((0, 0)), vector((2, 2))))
        2
        sage: is_segment_intersecting((vector((-1, -1)), vector((1, 1))), (vector((1, 1)), vector((2, 2))))
        1

    """
    turn_from_s = ccw(s[1] - s[0], t[0] - s[0]) * ccw(s[1] - s[0], t[1] - s[0])
    if turn_from_s > 0:
        # Both endpoints of t are on the same side of s
        return 0

    turn_from_t = ccw(t[1] - t[0], s[0] - t[0]) * ccw(t[1] - t[0], s[1] - t[0])
    if turn_from_t > 0:
        # Both endpoints of s are on the same side of t
        return 0

    Δs = s[1] - s[0]
    Δt = t[1] - t[0]
    if ccw(Δs, Δt) == 0:
        # Segments are parallel
        if is_anti_parallel(Δs, Δt):
            # Ensure segments are oriented the same way.
            t = (t[1], t[0])
            Δt = -Δt

        time_to_t0, length = time_on_ray(s[0], Δs, t[0])
        time_to_t1, _ = time_on_ray(s[0], Δs, t[1])

        if time_to_t0 < 0 and time_to_t1 < 0:
            # Both endpoints of t are earlier than s
            return 0

        if time_to_t0 > length and time_to_t1 > length:
            # Both endpoints of t are later than s
            return 0

        if time_to_t0 == length:
            return 1

        if time_to_t1 == 0:
            return 1

        return 2

    if turn_from_t == 0 and turn_from_s == 0:
        return 1

    return 2


def is_between(begin, end, v):
    r"""
    Check whether the vector ``v`` is strictly in the sector formed by the vectors
    ``begin`` and ``end`` (in counter-clockwise order).

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_between
        sage: is_between((1, 0), (1, 1), (2, 1))
        True

        sage: from itertools import product
        sage: vecs = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
        sage: for (i, vi), (j, vj), (k, vk) in product(enumerate(vecs), repeat=3):
        ....:     assert is_between(vi, vj, vk) == ((i == j and i != k) or i < k < j or k < j < i or j < i < k), ((i, vi), (j, vj), (k, vk))
    """
    if begin[0] * end[1] > end[0] * begin[1]:
        # positive determinant
        # [ begin[0] end[0] ]^-1 = [ end[1] -end[0] ]
        # [ begin[1] end[1] ]      [-begin[1]  begin[0] ]
        # v[0] * end[1] - end[0] * v[1] > 0
        # - v[0] * begin[1] + begin[0] * v[1] > 0
        return end[1] * v[0] > end[0] * v[1] and begin[0] * v[1] > begin[1] * v[0]
    elif begin[0] * end[1] == end[0] * begin[1]:
        # aligned vector
        if begin[0] * end[0] >= 0 and begin[1] * end[1] >= 0:
            return v[0] * begin[1] != v[1] * begin[0] or v[0] * begin[0] < 0 or v[1] * begin[1] < 0
        else:
            return begin[0] * v[1] > begin[1] * v[0]
    else:
        # negative determinant
        # [ end[0] begin[0] ]^-1 = [ begin[1] -begin[0] ]
        # [ end[1] begin[1] ]      [-end[1]  end[0] ]
        # v[0] * begin[1] - begin[0] * v[1] > 0
        # - v[0] * end[1] + end[0] * v[1] > 0
        return begin[1] * v[0] < begin[0] * v[1] or end[0] * v[1] < end[1] * v[0]


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
