# -*- coding: utf-8 -*-
r"""
Some tools for 2x2 matrices and planar geometry.
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2022 Julian Rüth
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
from __future__ import absolute_import, print_function, division

from sage.rings.all import AA, QQbar, RR

from sage.matrix.constructor import matrix, identity_matrix
from sage.modules.free_module_element import vector

def similarity_from_vectors(u, v, matrix_space=None):
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

    if matrix_space is None:
        from sage.matrix.matrix_space import MatrixSpace
        matrix_space = MatrixSpace(u.base_ring(), 2)

    if u == v:
        return matrix_space.one()

    sqnorm_u = u[0]*u[0] + u[1]*u[1]
    cos_uv = (u[0]*v[0] + u[1]*v[1]) / sqnorm_u
    sin_uv = (u[0]*v[1] - u[1]*v[0]) / sqnorm_u

    m = matrix_space([cos_uv, -sin_uv, sin_uv, cos_uv])
    m.set_immutable()
    return m


def is_cosine_sine_of_rational(c, s):
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

    TESTS::

        sage: from pyexactreal import ExactReals # optional: exactreal
        sage: R = ExactReals() # optional: exactreal
        sage: is_cosine_sine_of_rational(R.one(), R.zero()) # optional: exactreal
        True

    """
    return (QQbar(c) + QQbar.gen() * QQbar(s)).minpoly().is_cyclotomic()


def angle(u, v, numerical=False, assume_rational=False):
    r"""
    Return the angle between the vectors ``u`` and ``v`` divided by `2 \pi`.

    INPUT:

    - ``u``, ``v`` - vectors

    - ``numerical`` - boolean, whether to return floating point numbers

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

    The numerical version (working over floating point numbers)::

        sage: import math
        sage: u = (1, 0)
        sage: for n in xsrange(1,20):
        ....:     for k in xsrange(1,n):
        ....:         a = 2 * k * math.pi / n
        ....:         v = (math.cos(a), math.sin(a))
        ....:         assert abs(angle(u,v,numerical=True) * 2 * math.pi - a) < 1.e-10

    And we test up to 50 when setting ``assume_rational`` to ``True``::

        sage: for n in xsrange(1,20):       # long time
        ....:     for k in xsrange(1,n):
        ....:         v = vector((AA(cos(2*k*pi/n)), AA(sin(2*k*pi/n))))
        ....:         assert angle(u,v,assume_rational=True) == k/n

    If the angle is not rational, then the method returns an element in the real
    lazy field::

        sage: v = vector((AA(sqrt(2)), AA(sqrt(3))))
        sage: a = angle(u, v)
        sage: a    # abs tol 1e-14
        0.14102355421224375
        sage: exp(2*pi.n()*CC(0,1)*a)
        0.632455532033676 + 0.774596669241483*I
        sage: v / v.norm()
        (0.6324555320336758?, 0.774596669241484?)
    """
    if not assume_rational and not numerical:
        sqnorm_u = u[0] * u[0] + u[1] * u[1]
        sqnorm_v = v[0] * v[0] + v[1] * v[1]

        if sqnorm_u != sqnorm_v:
            uu = vector(AA, u)
            vv = (AA(sqnorm_u) / AA(sqnorm_v)).sqrt() * vector(AA, v)
        else:
            uu = u
            vv = v

        cos_uv = (uu[0]*vv[0] + uu[1]*vv[1]) / sqnorm_u
        sin_uv = (uu[0]*vv[1] - uu[1]*vv[0]) / sqnorm_u

        is_rational = is_cosine_sine_of_rational(cos_uv, sin_uv)
    elif assume_rational:
        is_rational = True

    import math

    u0 = float(u[0]); u1 = float(u[1])
    v0 = float(v[0]); v1 = float(v[1])

    cos_uv = (u0*v0 + u1*v1) / math.sqrt((u0*u0 + u1*u1)*(v0*v0 + v1*v1))
    if cos_uv < -1.0:
        assert cos_uv > -1.0000001
        cos_uv = -1.0
    elif cos_uv > 1.0:
        assert cos_uv < 1.0000001
        cos_uv = 1.0
    angle = math.acos(cos_uv) / (2 * math.pi)   # rat number between 0 and 1/2

    if numerical or not is_rational:
        return 1.0 - angle if u0 * v1 - u1*v0 < 0 else angle
    else:
        # fast and dirty way using floating point approximation
        # (see below for a slow but exact method)
        angle_rat = RR(angle).nearby_rational(0.00000001)
        if angle_rat.denominator() > 100:
            raise NotImplementedError("the numerical method used is not smart enough!")
        return 1 - angle_rat if u0*v1 - u1*v0 < 0 else angle_rat

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