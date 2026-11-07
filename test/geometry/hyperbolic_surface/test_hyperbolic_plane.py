r"""
Checks of the HyperbolicPlane against classical hyperbolic geometry.
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

from sage.all import (
    AA,
    QQ,
    RR,
    I,
    NumberField,
    arctanh,
    cos,
    cosh,
    matrix,
    pi,
    polygen,
    sin,
    sqrt,
)

from flatsurf import HyperbolicPlane


@pytest.fixture(scope="module")
def octagon_field():
    x = polygen(QQ, "x")
    K = NumberField(x**4 + 32 * x**2 - 32, "r", embedding=RR(0.98517143100941))
    return K, K.gen()


@pytest.fixture(scope="module")
def octagon(octagon_field):
    K, r = octagon_field
    s2 = r**2 / 12 + QQ(4) / 3
    H = HyperbolicPlane(K)
    points = [
        (r, 0),
        (r / s2, r / s2),
        (0, r),
        (-r / s2, r / s2),
        (-r, 0),
        (-r / s2, -r / s2),
        (0, -r),
        (r / s2, -r / s2),
    ]
    return H.convex_hull(*[H.point(x, y, model="klein") for (x, y) in points])


def test_rho_klein_rational():
    H = HyperbolicPlane(QQ)
    rho = H.half_circle(0, 1).intersection(H.vertical(-QQ(1) / 2))
    assert rho.coordinates(model="klein") == (-QQ(1) / 2, 0)
    assert rho.apply_isometry(matrix(QQ, [[0, -1], [1, 1]])) == rho


def test_octagon_area(octagon_field, octagon):
    _, r = octagon_field
    s2 = r**2 / 12 + QQ(4) / 3
    assert s2**2 == 2
    assert len(list(octagon.vertices())) == 8
    assert octagon.area(numerical=False) == 2


def test_octagon_angles(octagon):
    E = list(octagon.edges())
    assert len(E) == 8
    for i in range(8):
        turning = E[i].geodesic().angle(E[(i + 1) % 8].geodesic(), numerical=True)
        assert abs(turning - 3 / 8) < 1e-6


def test_octagon_cosine_rule():
    r = RR(0.98517143100941)
    target = 3 + 2 * sqrt(RR(2))
    lhs = (cos(pi / 4) * cos(pi / 8) + cos(pi / 8)) / (sin(pi / 4) * sin(pi / 8))
    assert abs(RR(lhs) - target) < 1e-9
    assert abs(cosh(arctanh(r)) - target) < 1e-9


def test_compact_triangle_gauss_bonnet():
    H = HyperbolicPlane(QQ)
    T = H.convex_hull(H(I), H(2 * I), H(1 + I))
    E = list(T.edges())
    assert len(E) == 3
    turning = [
        E[i].geodesic().angle(E[(i + 1) % 3].geodesic(), numerical=True)
        for i in range(3)
    ]
    assert abs(T.area(numerical=True) - (QQ(1) / 2 - sum(QQ(1) / 2 - t for t in turning))) < 1e-9
