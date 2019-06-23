r"""
Some tools for 2x2 matrices and planar geometry.
"""
from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from sage.misc.cachefunc import cached_function

from sage.rings.all import ZZ, QQ, AA, QQbar, RR, CC, RDF, CDF, RIF, CIF
from sage.rings.rational import Rational
from sage.rings.complex_interval_field import ComplexIntervalField


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

    EXAMPLES::

        sage: from flatsurf.geometry.matrix_2x2 import is_similarity

        sage: is_similarity(matrix([[0,1],[1,0]]))
        True
        sage: is_similarity(matrix([[0,-2],[2,0]]))
        True
        sage: is_similarity(matrix([[1,1],[0,1]]))
        False
    """
    n = m * m.transpose()
    return n[0,1].is_zero() and n[1,0].is_zero()

def homothety_rotation_decomposition(m):
    r"""
    Return a couple composed of the homothety and a rotation matrix.

    The coefficients of the returned pair are either in the ground field of
    ``m`` or in the algebraic field ``AA``.

    EXAMPLES::

        sage: from flatsurf.geometry.matrix_2x2 import homothety_rotation_decomposition

        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=1.4142)
        sage: m = matrix([[sqrt2, -sqrt2],[sqrt2,sqrt2]])
        sage: a,rot = homothety_rotation_decomposition(m)
        sage: a
        2
        sage: rot
        [ 1/2*sqrt2 -1/2*sqrt2]
        [ 1/2*sqrt2  1/2*sqrt2]
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

        sage: from flatsurf.geometry.matrix_2x2 import similarity_from_vectors

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
        return identity_matrix(u.base_ring(), n=2)

    sqnorm_u = u[0]*u[0] + u[1]*u[1]
    # Editted by Pat to remove worry about subfield...
    #sqnorm_v = v[0]*v[0] + v[1]*v[1]

    #if sqnorm_u != sqnorm_v:
    #    r = sqnorm_u / sqnorm_v
    #    if not r.is_square():
    #        raise ValueError("there is no similarity in the ground field; consider a suitable field extension. Note: u="+str(u)+" and v="+str(v)+".")
    #    sqrt_r = r.sqrt()
    #else:
    #    sqrt_r = 1

    #vv = sqrt_r * v

    #cos_uv = (u[0]*vv[0] + u[1]*vv[1]) / sqnorm_u
    #sin_uv = (u[0]*vv[1] - u[1]*vv[0]) / sqnorm_u
    # return 1/sqrt_r * matrix([[cos_uv, -sin_uv],[sin_uv, cos_uv]])

    cos_uv = (u[0]*v[0] + u[1]*v[1]) / sqnorm_u
    sin_uv = (u[0]*v[1] - u[1]*v[0]) / sqnorm_u
    return matrix([[cos_uv, -sin_uv],[sin_uv, cos_uv]])


def rotation_matrix_angle(r, check=False):
    r"""
    Return the angle of the rotation matrix ``r`` divided by ``2 pi``.

    EXAMPLES::

        sage: from flatsurf.geometry.matrix_2x2 import rotation_matrix_angle

        sage: def rot_matrix(p, q):
        ....:     z = QQbar.zeta(q) ** p
        ....:     c = z.real()
        ....:     s = z.imag()
        ....:     return matrix(AA, 2, [c,-s,s,c])
        sage: [rotation_matrix_angle(rot_matrix(i, 5)) for i in range(1,5)]
        [1/5, 2/5, 3/5, 4/5]
        sage: [rotation_matrix_angle(rot_matrix(i,7)) for i in range(1,7)]
        [1/7, 2/7, 3/7, 4/7, 5/7, 6/7]

    Some random tests::

        sage: for _ in range(100):
        ....:     r = QQ.random_element(x=0,y=500)
        ....:     r -= r.floor()
        ....:     m = rot_matrix(r.numerator(), r.denominator())
        ....:     assert rotation_matrix_angle(m) == r

    .. NOTE::

        This is using floating point arithmetic and might be wrong.
    """
    e0,e1 = r.change_ring(CDF).eigenvalues()
    m0 = (e0.log() / 2 / CDF.pi()).imag()
    m1 = (e1.log() / 2 / CDF.pi()).imag()
    r0 = RR(m0).nearby_rational(max_denominator=10000)
    r1 = RR(m1).nearby_rational(max_denominator=10000)
    if r0 != -r1:
        raise RuntimeError
    r0 = r0.abs()
    if r[0][1] > 0:
        return QQ.one() - r0
    else:
        return r0

    if check:
        e = r.change_ring(AA).eigenvalues()[0]
        if e.minpoly() != ZZ['x'].cyclotomic_polynomial()(r.denominator()):
            raise RuntimeError
        z = QQbar.zeta(r.denominator())
        if z**r.numerator() != e:
            raise RuntimeError

    return r


def is_cosine_sine_of_rational(c,s):
    r"""
    Check whether the given pair is a cosine and sine of a same rational angle.

    EXAMPLES::

        sage: from flatsurf.geometry.matrix_2x2 import is_cosine_sine_of_rational

        sage: c = s = AA(sqrt(2))/2
        sage: is_cosine_sine_of_rational(c,s)
        True
        sage: c = AA(sqrt(3))/2; s = AA(1/2)
        sage: is_cosine_sine_of_rational(c,s)
        True

        sage: c = AA(sqrt(5)/2); s = (1 - c**2).sqrt()
        sage: c**2 + s**2
        1.000000000000000?
        sage: is_cosine_sine_of_rational(c,s)
        False

        sage: c = (AA(sqrt(5)) + 1)/4; s = (1 - c**2).sqrt()
        sage: is_cosine_sine_of_rational(c,s)
        True

        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: is_cosine_sine_of_rational(K.zero(),-K.one())
        True
    """
    return (QQbar(c) + QQbar.gen() * QQbar(s)).minpoly().is_cyclotomic()

def angle(u, v, assume_rational=False):
    r"""
    Return the angle between the vectors ``u`` and ``v`` divided by `2 \pi`.

    INPUT:

    - ``u``, ``v`` - vectors

    - ``assume_rational`` - whether we assume that the angle is a multiple
      rational of ``pi``. By default it is ``False`` but if it is known in
      advance that the result is rational then setting it to ``True`` might be
      much faster.

    EXAMPLES::

        sage: from flatsurf.geometry.matrix_2x2 import angle

    As the implementation is dirty, we at least check that it works for all
    denominator up to 20::

        sage: u = vector((AA(1),AA(0)))
        sage: for n in xsrange(1,20):       # long time  (10 sec)
        ....:     for k in xsrange(1,n):
        ....:         v = vector((AA(cos(2*k*pi/n)), AA(sin(2*k*pi/n))))
        ....:         assert angle(u,v) == k/n

    And we test up to 50 when setting ``assume_rational`` to ``True``::

        sage: for n in xsrange(1,20):       # long time
        ....:     for k in xsrange(1,n):
        ....:         v = vector((AA(cos(2*k*pi/n)), AA(sin(2*k*pi/n))))
        ....:         assert angle(u,v,assume_rational=True) == k/n

    If the angle is not rational, then the method returns an element in the real
    lazy field::

        sage: v = vector((AA(sqrt(2)), AA(sqrt(3))))
        sage: a = angle(u,v)
        sage: a
        0.1410235542122437?
        sage: exp(2*pi.n()*CC(0,1)*a.n())
        0.632455532033676 + 0.774596669241483*I
        sage: v / v.norm()
        (0.6324555320336758?, 0.774596669241484?)
    """
    if not assume_rational:
        sqnorm_u = u[0] * u[0] + u[1] * u[1]
        sqnorm_v = v[0] * v[0] + v[1] * v[1]

        if sqnorm_u != sqnorm_v:
            uu = u.change_ring(AA)
            vv = (AA(sqnorm_u) / AA(sqnorm_v)).sqrt() * v.change_ring(AA)
        else:
            uu = u
            vv= v

        cos_uv = (uu[0]*vv[0] + uu[1]*vv[1]) / sqnorm_u
        sin_uv = (uu[0]*vv[1] - uu[1]*vv[0]) / sqnorm_u

        is_rational = is_cosine_sine_of_rational(cos_uv, sin_uv)

    else:
        is_rational = True

    if is_rational:
        # fast and dirty way using floating point approximation
        # (see below for a slow but exact method)
        from math import acos,asin,sqrt

        u0 = float(u[0]); u1 = float(u[1])
        v0 = float(v[0]); v1 = float(v[1])

        cos_uv = (u0*v0 + u1*v1) / sqrt((u0*u0 + u1*u1)*(v0*v0 + v1*v1))
        angle = acos(float(cos_uv)) / (2*pi_float)   # rat number between 0 and 1/2
        angle_rat = RR(angle).nearby_rational(0.00000001)
        if angle_rat.denominator() > 100:
            raise NotImplementedError("the numerical method used is not smart enough!")
        if u0*v1 - u1*v0 < 0:
            return 1 - angle_rat
        return angle_rat

    else:
        from sage.functions.trig import acos
        from sage.rings.real_lazy import RLF
        from sage.symbolic.constants import pi

        if sin_uv > 0:
            return acos(RLF(cos_uv)) / RLF(2*pi)
        else:
            return -acos(RLF(cos_uv)) / RLF(2*pi)

    # a neater way is provided below by working only with number fields
    # but this method is slower...
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
