r"""
Checks of the isomorphisms between hyperbolic planes of the quadratic forms
`x^2 + 8y^2 - 12z^2`, `x^2 + 2y^2 - 3z^2` and `x^2 + y^2 - z^2` against
hard-coded matrices, images, side pairings, angles and areas.
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

import re

import pytest

import sage.misc.sage_unittest

from sage.all import AA, QQ, NumberField, diagonal_matrix, identity_matrix, matrix, polygen

from flatsurf import HyperbolicPlane
from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface


A3 = matrix(
    QQ,
    [
        [QQ(25) / 4, QQ(23) / 2, -QQ(69) / 4],
        [0, 6, -6],
        [-QQ(23) / 4, -QQ(25) / 2, QQ(75) / 4],
    ],
)

KITE_ANGLES = (
    0.3040867239846964,
    0.14102355421224375,
    0.38101825537203543,
    0.14102355421224375,
)

TRIANGLE_ANGLES = (0.25, 0.14875643774533862, 0.09122814004003268)


@pytest.fixture
def H():
    return HyperbolicPlane(QQ)


@pytest.fixture
def H3():
    return HyperbolicPlane(QQ, form=[1, 2, -3])


@pytest.fixture
def H8():
    return HyperbolicPlane(QQ, form=[1, 8, -12])


@pytest.fixture
def kite(H3):
    return H3.polygon(
        vertices=[(0, 0), (QQ(1) / 2, -QQ(1) / 2), (QQ(3) / 4, 0), (QQ(1) / 2, QQ(1) / 2)]
    )


@pytest.fixture
def side_pairing():
    return HyperbolicIsometryGroup(QQ, form=[1, 2, -3])(
        [[-QQ(1) / 3, -QQ(4) / 3, 0], [QQ(2) / 3, -QQ(1) / 3, 0], [0, 0, 1]]
    )


def klein(polygon):
    return [v.coordinates(model="klein") for v in polygon.vertices()]


def sqrt3():
    x = polygen(QQ, "x")
    K = NumberField(x**2 - 3, "a", embedding=AA(3).sqrt())
    return K, K.gen()


def test_rescaling(H8, H3):
    f = H8.isomorphism(H3)
    assert f.matrix() == diagonal_matrix(QQ, [1, 2, 2])
    assert f.matrix().transpose() * H3.quadratic_form() * f.matrix() == H8.quadratic_form()
    assert f.matrix().det() == 4
    p = f(H8.point(QQ(1) / 4, QQ(1) / 8, model="klein"))
    assert p.parent() is H3
    assert p.coordinates(model="klein") == (QQ(1) / 8, QQ(1) / 8)
    assert f.domain() is H8
    assert f.codomain() is H3


def test_rescaling_refused(H3, H):
    with pytest.raises(
        ValueError,
        match=re.escape(
            "isometry is not defined over Rational Field; adjoin a square root of 2"
        ),
    ):
        H3.isomorphism(H)


def test_kite(H3, H, kite):
    f = H3.isomorphism(H, matrix=A3)
    assert A3.transpose() * H.quadratic_form() * A3 == 6 * H3.quadratic_form()
    assert A3.det() == 36
    Q = f(kite)
    assert Q.parent() is H
    assert klein(Q) == [
        (-QQ(23) / 25, -QQ(8) / 25),
        (-QQ(53) / 59, -QQ(24) / 59),
        (-QQ(67) / 77, -QQ(32) / 77),
        (-QQ(67) / 77, -QQ(24) / 77),
    ]
    assert Q._cosh2_sides() == (QQ(4) / 3, QQ(49) / 39, QQ(49) / 39, QQ(4) / 3)
    assert kite._cosh2_sides() == Q._cosh2_sides()
    assert Q.angles(numerical=True) == pytest.approx(KITE_ANGLES)
    assert kite.angles(numerical=True) == pytest.approx(KITE_ANGLES)
    assert Q.area(numerical=True) == pytest.approx(0.032847912218780584)
    assert kite.area(numerical=True) == pytest.approx(0.032847912218780584)


def test_side_pairing(H3, H, side_pairing):
    f = H3.isomorphism(H, matrix=A3)
    g = f(side_pairing)
    assert g.parent() is HyperbolicIsometryGroup(QQ)
    assert g.matrix() == matrix(
        QQ,
        [
            [-QQ(179) / 8, -QQ(71) / 6, -QQ(607) / 24],
            [-QQ(7) / 2, -3, -QQ(9) / 2],
            [QQ(181) / 8, QQ(73) / 6, QQ(617) / 24],
        ],
    )
    assert g.trace() == QQ(1) / 3
    phi = f.isometry_group_morphism()
    assert phi.domain() is HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
    assert phi.codomain() is HyperbolicIsometryGroup(QQ)
    assert phi(side_pairing) == g
    assert (~f)(g) == side_pairing


def test_surface(H3, H, kite, side_pairing):
    S = MutableOrientedHyperbolicSurface(H3)
    S.add_polygon(kite)
    S.glue((0, 0), (0, 3))
    S.glue((0, 1), (0, 2))
    S.set_immutable()

    f = H3.isomorphism(H, matrix=A3)
    T = f(S)
    assert T.hyperbolic_plane() is H
    assert not T.is_mutable()
    assert T.edge_transformation(0, 0) == f(side_pairing)
    assert sorted(T.angles(numerical=True)) == pytest.approx(
        [0.2820471084244875, 0.3040867239846964, 0.38101825537203543]
    )
    sage.misc.sage_unittest.TestSuite(T).run(raise_on_failure=True)

    assert (~f)(T) == S


def test_mutable_surface(H8, H3):
    S = MutableOrientedHyperbolicSurface(H8)
    S.add_polygon(H8.polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]))
    T = H8.isomorphism(H3)(S)
    assert T.is_mutable()
    assert klein(T.polygon(0)) == [(0, 0), (QQ(1) / 4, 0), (0, QQ(1) / 2)]


def test_sqrt3():
    K, a = sqrt3()
    H3K = HyperbolicPlane(K, form=[1, 2, -3])
    HK = HyperbolicPlane(K)
    B = matrix(
        K,
        [[QQ(7) / 12 * a, 0, -QQ(5) / 4], [0, 1, 0], [-QQ(5) / 12 * a, 0, QQ(7) / 4]],
    )
    assert B.transpose() * HK.quadratic_form() * B == H3K.quadratic_form() / 2
    assert B.det() == a / 2
    f = H3K.isomorphism(HK, matrix=B)

    P = H3K.polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)])
    Q = f(P)
    assert klein(Q) == [
        (-QQ(5) / 7, 0),
        (QQ(48) / 563 * a - QQ(385) / 563, 0),
        (-QQ(5) / 7, QQ(2) / 7),
    ]
    assert Q.angles(numerical=True) == pytest.approx(TRIANGLE_ANGLES)
    assert P.angles(numerical=True) == pytest.approx(TRIANGLE_ANGLES)
    assert klein((~f)(Q)) == [(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]
    assert (~f)(Q) == P

    H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
    PQ = H3.polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)])
    assert f(PQ) == Q
    assert f(H3.point(0, 0, model="klein")) == HK.point(-QQ(5) / 7, 0, model="klein")


def test_inverse_and_composition(H8, H3, H):
    f = H8.isomorphism(H3)
    g = H3.isomorphism(H, matrix=A3)
    assert (~f).matrix() == diagonal_matrix(QQ, [1, QQ(1) / 2, QQ(1) / 2])
    assert f.section() == ~f
    assert (~f).domain() is H3
    assert (~f).codomain() is H8

    h = g * f
    assert h.domain() is H8
    assert h.codomain() is H
    assert h.matrix() == A3 * diagonal_matrix(QQ, [1, 2, 2])

    assert (~g) * g == H3.isomorphism(H3)
    assert g * (~g) == H.isomorphism(H)
    assert ((~g) * g).is_identity()
    assert H3.isomorphism(H3).matrix() == identity_matrix(QQ, 3)

    p = H8.point(QQ(1) / 4, QQ(1) / 8, model="klein")
    assert h(p) == g(f(p))


def test_equality_up_to_positive_scalar(H3, H):
    f = H3.isomorphism(H, matrix=A3)
    assert f == H3.isomorphism(H, matrix=2 * A3)
    assert hash(f) == hash(H3.isomorphism(H, matrix=2 * A3))
    assert f == H3.isomorphism(H, matrix=A3 / 3)
    assert f != H3.isomorphism(H, matrix=-A3, check=False)
    assert H3.isomorphism(H3) != H3.isomorphism(H3, matrix=diagonal_matrix(QQ, [-1, -1, 1]))


def test_refusals(H8, H3, H):
    with pytest.raises(ValueError, match="positive multiple of the form of the domain"):
        H8.isomorphism(H3, matrix=identity_matrix(QQ, 3))
    with pytest.raises(ValueError, match="positive last coordinate"):
        H8.isomorphism(H3, matrix=-diagonal_matrix(QQ, [1, 2, 2]))
    with pytest.raises(ValueError, match="positive determinant"):
        H8.isomorphism(H3, matrix=diagonal_matrix(QQ, [1, -2, 2]))
    with pytest.raises(ValueError, match="same base ring"):
        H3.isomorphism(HyperbolicPlane(AA))

    f = H8.isomorphism(H3)
    with pytest.raises(NotImplementedError, match="cannot map Integer"):
        f(QQ(42).numerator())
    with pytest.raises(ValueError, match="does not lie in the domain"):
        f(H.point(0, 0, model="klein"))
    with pytest.raises(ValueError, match="does not lie in the domain"):
        f(HyperbolicIsometryGroup(QQ).one())
