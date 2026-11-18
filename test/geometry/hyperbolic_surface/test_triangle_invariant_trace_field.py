r"""
The field `k_0 = \QQ(u_0, u_1, u_2, p)` of a hyperbolic triangle, with
`u_k = 4 \cos^2 \alpha_k` and `p = 8 \cos \alpha_0 \cos \alpha_1 \cos
\alpha_2`, against the fields of ``TRIANGLE_DOUBLES`` and hard-coded degrees.

For the angles `\pi / m_k`, `k_0` is the invariant trace field of the
triangle group of periods `(m_0, m_1, m_2)`, Takeuchi, Arithmetic triangle
groups (1977), Propositions 4 and 5: `\QQ` for `(2, 3, \infty)`. The degrees
were computed independently, with ``qqbar`` for the angles that are not of
the form `\pi / m`.
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

import pytest

from sage.all import AA, QQ, QQbar, NumberField, polygen

from flatsurf.geometry.subfield import (
    _triangle_invariant_trace_field,
    number_field_elements_from_algebraics,
)

from hyperbolic_surface_oracles import TRIANGLE_DOUBLES

y = polygen(QQ, "y")


def cosine(angle):
    r"""
    Return `\cos(2 \pi \theta)` in ``AA`` for the rational ``angle``
    `\theta`.
    """
    angle = QQ(angle)
    return AA((QQbar.zeta(angle.denominator()) ** angle.numerator()).real())


def same_field(K, L):
    r"""
    Return whether ``K`` and ``L``, each ``QQ`` or a number field, are
    isomorphic.
    """
    if K is QQ or L is QQ:
        return K.degree() == 1 and L.degree() == 1
    return K.is_isomorphic(L)


@pytest.mark.parametrize(
    "periods,polynomial",
    [(periods, k0) for periods, _, _, k0 in TRIANGLE_DOUBLES],
)
def test_triangle_doubles(periods, polynomial):
    angles = [QQ(1) / (2 * m) for m in periods]
    K, _, _ = _triangle_invariant_trace_field(angles)
    L = QQ if len(polynomial) == 2 else NumberField(y.parent()(polynomial), "b")
    assert same_field(K, L)


ANGLES = [
    ((QQ(1) / 10, QQ(1) / 14, QQ(1) / 16), 24),
    ((QQ(1) / 4, QQ(1) / 6, QQ(1) / 14), 3),
    ((QQ(1) / 6, QQ(1) / 10, QQ(1) / 14), 6),
    ((QQ(1) / 4, QQ(1) / 8, QQ(1) / 12), 1),
    ((QQ(1) / 6, QQ(1) / 6, QQ(1) / 8), 2),
    ((QQ(0), QQ(1) / 4, QQ(1) / 6), 1),
    ((QQ(2) / 15, QQ(1) / 9, QQ(1) / 10), 12),
    ((QQ(1) / 5, QQ(1) / 10, QQ(1) / 12), 4),
    ((QQ(3) / 14, QQ(1) / 14, QQ(1) / 14), 3),
]


@pytest.mark.parametrize("angles,degree", ANGLES)
def test_degree(angles, degree):
    K, u, p = _triangle_invariant_trace_field(angles)
    assert K.degree() == degree
    assert p**2 == u[0] * u[1] * u[2]


@pytest.mark.parametrize("angles,degree", ANGLES)
def test_values(angles, degree):
    K, u, p = _triangle_invariant_trace_field(angles)
    c = [cosine(angle) for angle in angles]
    for uk, ck in zip(u, c):
        assert AA(uk) == 4 * ck**2
    # With p^2 = u_0 u_1 u_2, p is +-8 c_0 c_1 c_2 exactly, and the sign
    # decides; a direct comparison in AA builds the compositum of the
    # fields of the c_k, of degree up to phi(N) / 2.
    assert p**2 == u[0] * u[1] * u[2]
    q = 8 * c[0] * c[1] * c[2]
    assert AA(p).sign() == q.sign()


def test_cusp():
    K, u, p = _triangle_invariant_trace_field([QQ(0), QQ(1) / 4, QQ(1) / 6])
    assert K is QQ
    assert u == (4, 0, 1)
    assert p == 0


@pytest.mark.parametrize(
    "angles",
    [
        (QQ(2) / 15, QQ(1) / 9, QQ(1) / 10),
        (QQ(1) / 5, QQ(1) / 10, QQ(1) / 12),
        (QQ(3) / 14, QQ(1) / 14, QQ(1) / 14),
    ],
)
def test_qqbar(angles):
    K, _, _ = _triangle_invariant_trace_field(angles)
    c = [cosine(angle) for angle in angles]
    L, _ = number_field_elements_from_algebraics(
        [c[0] ** 2, c[1] ** 2, c[2] ** 2, c[0] * c[1] * c[2]]
    )
    assert same_field(K, L)


@pytest.mark.parametrize("angles,degree", ANGLES)
def test_order(angles, degree):
    _, u, _ = _triangle_invariant_trace_field(angles)
    permuted = (angles[2], angles[0], angles[1])
    _, v, _ = _triangle_invariant_trace_field(permuted)
    assert tuple(AA(x) for x in v) == (AA(u[2]), AA(u[0]), AA(u[1]))
