r"""
Some tools for 2x2 matrices.
"""

from sage.rings.qqbar import AA
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RR

from math import pi as pi_float

from sage.symbolic.constants import pi
from sage.matrix.constructor import matrix, identity_matrix


def number_field_to_AA(a):
    r"""
    It is a mess to convert an element of a number field to the algebraic field
    ``AA``. This is a temporary fix.
    """
    try:
        return AA(a)
    except TypeError:
        return AA.polynomial_root(a.minpoly(), RIF(a))

def is_similarity(m):
    r"""
    Return ``True`` if ``m`` is a similarity and ``False`` otherwise.
    """
    n = m * m.transpose()
    return n[0,1] == 0 and n[1,0] == 0

def homothety_rotation_decomposition(m):
    r"""
    Return a couple composed of the homothety and a rotation matrix.

    The coefficients of the returned pair are either in the ground field of
    ``m`` or in the algebraic field ``AA``.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=1.4142)
        sage: m = matrix([[sqrt2, -sqrt2],[sqrt2,sqrt2]])
        sage: a,rot = homothety_rotation_decomposition(m)
        sage: a
        2.000000000000000?
        sage: rot
        [ 0.7071067811865475? -0.7071067811865475?]
        [ 0.7071067811865475?  0.7071067811865475?]
    """
    if not is_similarity(m):
        raise ValueError("the matrix must be a similarity")

    det = m.det()

    if not det.is_square():
        if not AA.has_coerce_map_from(m.base_ring()):
            l = map(number_field_to_AA,m.list())
            M = MatrixSpace(AA,2)
            m = M(l)
        else:
            m = m.change_ring(AA)

    sqrt_det = det.sqrt()

    return sqrt_det, m / sqrt_det

def similarity_from_vectors(u,v):
    r"""
    Return the unique similarity matrix that maps ``u`` to ``v``.

    EXAMPLES::

        sage: V = VectorSpace(QQ,2)
        sage: u = V((1,0))
        sage: v = V((0,1))
        sage: m = similarity_from_vectors(u,v); m
        [ 0 -1]
        [ 1  0]
        sage: m*u == v
        True

        sage: u = V((2,1))
        sage: v = V((1,-2))
        sage: m = similarity_from_vectors(u,v); m
        [ 0  1]
        [-1  0]
        sage: m * u == v
        True

    An example built from the Pythagorean triple 3^2 + 4^2 = 5^2::

        sage: u2 = V((5,0))
        sage: v2 = V((3,4))
        sage: m = similarity_from_vectors(u2,v2); m
        [ 3/5 -4/5]
        [ 4/5  3/5]
        sage: m * u2 == v2
        True

    Some test over number fields::

        sage: K.<sqrt2> = NumberField(x^2-2, embedding=1.4142)
        sage: V = VectorSpace(K,2)
        sage: u = V((sqrt2,0))
        sage: v = V((1, 1))
        sage: m = similarity_from_vectors(u,v); m
        [ 1/2*sqrt2 -1/2*sqrt2]
        [ 1/2*sqrt2  1/2*sqrt2]
        sage: m*u == v
        True

        sage: m = similarity_from_vectors(u, 2*v); m
        [ sqrt2 -sqrt2]
        [ sqrt2  sqrt2]
        sage: m*u == 2*v
        True
    """
    assert u.parent() is v.parent()

    if u == v:
        return identity_matrix(2,u.base_ring())

    sqnorm_u = u[0]*u[0] + u[1]*u[1]
    sqnorm_v = v[0]*v[0] + v[1]*v[1]

    if sqnorm_u != sqnorm_v:
        r = sqnorm_u / sqnorm_v
        if not r.is_square():
            raise ValueError("there is no similarity in the ground field; consider a suitable field extension.")
        sqrt_r = r.sqrt()
    else:
        sqrt_r = 1

    vv = sqrt_r * v

    cos_uv = (u[0]*vv[0] + u[1]*vv[1]) / sqnorm_u
    sin_uv = (u[0]*vv[1] - u[1]*vv[0]) / sqnorm_u
    return 1/sqrt_r * matrix([[cos_uv, -sin_uv],[sin_uv, cos_uv]])


def rotation_matrix_angle(r):
    r"""
    Return the angle of the rotation matrix ``r`` divided by ``2 pi`.

    EXAMPLES::

        sage: rot_matrix = lambda a: matrix(AA, [[cos(a),-sin(a)],[sin(a),cos(a)]])
        sage: rotation_matrix_angle(rot_matrix(pi/5))
        1/10
        sage: rotation_matrix_angle(rot_matrix(4*pi/7))
        2/7

    .. NOTE:

    The algoritm used is very naive and very slow!
    """
    assert (r * r.transpose()).is_one()

    if not AA.has_coerce_map_from(r.base_ring()):
        r = matrix(AA,2,map(number_field_to_AA,r.list()))

    # first compute the order
    n = 1
    rr = r
    while not rr.is_one():
        rr *= r
        n += 1

    if n == 1:
        return 1

    # then compute which one
    from sage.functions.trig import cos,sin
    c1 = AA(cos(2*pi/n))
    s1 = AA(sin(2*pi/n))
    rr = matrix([[c1,-s1],[s1,c1]])
    m = rr
    for i in xrange(1,n):
        if m == r:
            return Rational((i,n))
        m *= rr

    raise ValueError("an unexpected error occurred!")

def angle(u,v):
    r"""
    Return the angle between the vectors ``u`` and ``v``.

    This will not work if the angle is too small!!

    EXAMPLES:

    As the implementation is dirty, we at least check that it works for all
    denominator up to 100::

        sage: u = vector((AA(1),AA(0)))
        sage: for n in xsrange(1,100):
        ....:     print n
        ....:     for k in xsrange(1,n):
        ....:         v = vector((AA(cos(2*k*pi/n)), AA(sin(2*k*pi/n))))
        ....:         assert angle(u,v) == k/n, "failed for %d/%d get angle(u,v)=%d"%(k,n,angle(u,v))
    """
    # fast and dirty way
    from math import acos,asin,sqrt,pi

    u0 = float(u[0]); u1 = float(u[1])
    v0 = float(v[0]); v1 = float(v[1])

    cos_uv = (u0*v0 + u1*v1) / sqrt((u0*u0 + u1*u1)*(v0*v0 + v1*v1))
    angle = acos(float(cos_uv)) / (2*pi)   # rat number between 0 and 1/2
    angle_rat = RR(angle).nearby_rational(0.000001)
    if u0*v1 - u1*v0 < 0:
        return 1 - angle_rat
    return angle_rat

    # a neater way is provided below by working only with number fields
    #sqnorm_u = u[0]*u[0] + u[1]*u[1]
    #sqnorm_v = v[0]*v[0] + v[1]*v[1]
    #
    #if sqnorm_u != sqnorm_v:
    #    # we need to take a square root in order that u and v have the
    #    # same norm
    #    u = (1 / AA(sqnorm_u)).sqrt() * u.change_ring(AA)
    #    v = (1 / AA(sqnorm_v)).sqrt() * v.change_ring(AA)
    #    sqnorm_u = AA.one()
    #    sqnorm_v = AA.one()
    #
    #cos_uv = (u[0]*v[0] + u[1]*v[1]) / sqnorm_u
    #sin_uv = (u[0]*v[1] - u[1]*v[0]) / sqnorm_u
