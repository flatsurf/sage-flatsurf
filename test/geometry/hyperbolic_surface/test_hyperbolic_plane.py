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

from sage.all import QQ, RR, I, NumberField, matrix, polygen

from flatsurf import HyperbolicPlane

from hyperbolic_surface_oracles import OCTAGON, RHO


def sqrt2_in_radius(r):
    return sum(c * r**i for i, c in enumerate(OCTAGON["sqrt2_in_radius"]))


@pytest.fixture(scope="module")
def octagon_field():
    x = polygen(QQ, "x")
    f = sum(c * x**i for i, c in enumerate(OCTAGON["klein_radius_polynomial"]))
    K = NumberField(f, "r", embedding=RR(OCTAGON["klein_radius_embedding"]))
    return K, K.gen()


@pytest.fixture(scope="module")
def octagon(octagon_field):
    K, r = octagon_field
    s2 = sqrt2_in_radius(r)
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
    assert rho.coordinates(model="klein") == RHO["klein"]
    assert rho.apply_isometry(matrix(QQ, RHO["matrix"])) == rho


def test_octagon_area(octagon_field, octagon):
    _, r = octagon_field
    s2 = sqrt2_in_radius(r)
    assert s2**2 == 2
    assert len(list(octagon.vertices())) == 8
    assert octagon.area(numerical=False) == OCTAGON["area_over_2pi"]


def test_octagon_angles(octagon):
    E = list(octagon.edges())
    assert len(E) == 8
    for i in range(8):
        turning = E[i].geodesic().angle(E[(i + 1) % 8].geodesic(), numerical=True)
        assert abs(turning - OCTAGON["turning_angle"]) < 1e-6


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
