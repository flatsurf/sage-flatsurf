r"""
Checks of compact hyperbolic polygons given by their vertices against
hard-coded values of classical hyperbolic geometry.
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

from sage.all import QQ, RR, I, NumberField, polygen

from flatsurf import HyperbolicPlane

from hyperbolic_surface_oracles import OCTAGON


# Klein coordinates of a Lambert quadrilateral (three right angles).
QUADRILATERAL = [(0, 0), (QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]

# Klein coordinates of a non-convex quadrilateral, reflex at the origin.
ARROW = [(0, -QQ(1) / 2), (QQ(1) / 2, QQ(1) / 2), (0, 0), (-QQ(1) / 2, QQ(1) / 2)]


@pytest.fixture(scope="module")
def H():
    return HyperbolicPlane(QQ)


@pytest.fixture(scope="module")
def octagon():
    x = polygen(QQ, "x")
    f = sum(c * x**i for i, c in enumerate(OCTAGON["klein_radius_polynomial"]))
    K = NumberField(f, "r", embedding=RR(OCTAGON["klein_radius_embedding"]))
    r = K.gen()
    s2 = sum(c * r**i for i, c in enumerate(OCTAGON["sqrt2_in_radius"]))
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
    return HyperbolicPlane(K).polygon(vertices=points)


def test_quadrilateral(H):
    P = H.polygon(vertices=QUADRILATERAL)
    assert P._cosh2_sides() == (QQ(4) / 3, QQ(3) / 2, QQ(3) / 2, QQ(4) / 3)
    assert [P._angle_key(i)[2] for i in range(4)] == [0, 0, QQ(1) / 9, 0]
    for i in (0, 1, 3):
        assert P.angle(i, numerical=False) == QQ(1) / 4
    expected = [0.25, 0.25, 0.195913276015304, 0.25]
    for angle, value in zip(P.angles(numerical=True), expected):
        assert abs(angle - value) < 1e-9
    assert abs(P.area(numerical=True) - 0.0540867239846964) < 1e-9
    assert abs(H.convex_hull(*P.vertices()).area() - 0.0540867239846964) < 1e-9
    assert P.is_convex()
    assert P.is_simple()
    assert not P.is_degenerate()
    assert not P.is_rational()
    assert not P.is_equilateral()
    assert not P.is_equiangular()
    with pytest.raises(ValueError):
        P.angle(2, numerical=False)


def test_arrow(H):
    P = H.polygon(vertices=ARROW)
    assert P._cosh2_sides() == (QQ(25) / 6, QQ(2), QQ(2), QQ(25) / 6)
    assert P.angle(2, numerical=False) == QQ(3) / 4
    expected = [0.130073469146503, 0.0368407666897357, 0.75, 0.0368407666897357]
    for angle, value in zip(P.angles(numerical=True), expected):
        assert abs(angle - value) < 1e-9
    assert abs(P.area(numerical=True) - 0.0462449974740256) < 1e-9
    assert P.is_simple()
    assert not P.is_convex()
    assert not P.is_rational()
    assert P.describe_polygon() == (
        "a",
        "non-convex quadrilateral",
        "non-convex quadrilaterals",
    )


def test_octagon(octagon):
    assert octagon.is_rational()
    assert octagon.angles() == (OCTAGON["interior_angle"],) * 8
    assert octagon.area() == OCTAGON["area_over_2pi"]
    assert octagon.is_equilateral()
    assert octagon.is_equiangular()
    assert octagon.is_convex(strict=True)
    assert octagon.describe_polygon() == (
        "a",
        "regular octagon",
        "regular octagons",
    )


def test_point_position(H):
    P = H.polygon(vertices=ARROW)
    assert P.get_point_position(H(I)).is_vertex()
    assert (0, -QQ(1) / 4) in P
    assert (0, QQ(1) / 4) not in P


def test_ideal_vertex(H):
    with pytest.raises(ValueError):
        H.polygon(vertices=[(0, 0), (1, 0), (0, QQ(1) / 2)])
    with pytest.raises(ValueError):
        H.polygon(vertices=[H(0), H(I), H(1 + I)])
    with pytest.raises(ValueError):
        H.polygon(polygon=H.convex_hull(H(0), H(I), H(1 + I)))


def test_from_convex_polygon(H):
    hull = H.convex_hull(H(I), H(2 * I), H(1 + I))
    P = H.polygon(polygon=hull)
    assert P.vertices() == tuple(hull.vertices())


def test_inputs(H):
    hull = H.convex_hull(H(I), H(2 * I), H(1 + I))
    with pytest.raises(ValueError):
        H.polygon()
    with pytest.raises(ValueError):
        H.polygon(hull.half_spaces(), vertices=QUADRILATERAL)
    with pytest.raises(ValueError):
        H.polygon(hull.half_spaces(), polygon=hull)
    with pytest.raises(ValueError):
        H.polygon(vertices=QUADRILATERAL, polygon=hull)
