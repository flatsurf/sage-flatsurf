r"""
Checks of the hyperbolic plane with the quadratic form `x^2 + 2y^2 - 3z^2`
against hard-coded points, lines and isometries of its Klein model, the
ellipse `x^2 + 2y^2 < 3`.
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

from sage.all import (
    AA,
    QQ,
    RR,
    I,
    Graphics,
    NumberField,
    diagonal_matrix,
    dumps,
    loads,
    matrix,
    polygen,
)

from flatsurf import HyperbolicPlane


HALF_PLANE = re.escape("the half-plane model needs the form x^2 + y^2 - z^2")


@pytest.fixture
def H3():
    return HyperbolicPlane(QQ, form=[1, 2, -3])


def point(H, x, y, check=True):
    return H.point(QQ(x), QQ(y), model="klein", check=check)


def line(H, a, b, c, oriented=True, check=True):
    return H.geodesic(QQ(a), QQ(b), QQ(c), model="klein", oriented=oriented, check=check)


def test_unique_key(H3):
    assert HyperbolicPlane(QQ, form=[2, 2, -2]) is HyperbolicPlane(QQ)
    assert HyperbolicPlane(QQ, form=[2, 4, -6]) is H3
    assert HyperbolicPlane(QQ, form=diagonal_matrix(QQ, [1, 2, -3])) is H3
    assert H3.quadratic_form() == diagonal_matrix(QQ, [1, 2, -3])
    assert loads(dumps(H3)) is H3
    assert (
        repr(H3) == "Hyperbolic Plane over Rational Field with form x^2 + 2*y^2 - 3*z^2"
    )


def test_invalid_forms():
    with pytest.raises(ValueError, match="positive"):
        HyperbolicPlane(QQ, form=[1, 1, 1])
    with pytest.raises(ValueError, match="diagonal"):
        HyperbolicPlane(QQ, form=matrix(QQ, [[1, 1, 0], [1, 2, 0], [0, 0, -3]]))


def test_change_ring_and_coercion(H3):
    assert H3.change_ring(AA) is HyperbolicPlane(AA, form=[1, 2, -3])
    assert HyperbolicPlane(AA, form=[1, 2, -3]).has_coerce_map_from(H3)
    assert not HyperbolicPlane(AA).has_coerce_map_from(H3)
    assert not H3.has_coerce_map_from(QQ)


def test_point_predicates(H3):
    p = point(H3, QQ(1) / 2, QQ(1) / 2)
    assert p.is_finite() and not p.is_ideal() and not p.is_ultra_ideal()

    for x, y in [(1, 1), (-1, 1)]:
        q = point(H3, x, y)
        assert q.is_ideal() and not q.is_finite() and not q.is_ultra_ideal()

    q = point(H3, 2, 0, check=False)
    assert q.is_ultra_ideal()
    with pytest.raises(ValueError):
        point(H3, 2, 0)

    assert point(HyperbolicPlane(QQ), 1, 1, check=False).is_ultra_ideal()


def test_geodesic_predicates(H3):
    for g in [line(H3, -QQ(1) / 3, 1, 0), line(H3, -1, 0, 1)]:
        assert not g.is_ideal() and not g.is_ultra_ideal()

    tangent = line(H3, -3, 1, 2, check=False)
    assert tangent.is_ideal() and not tangent.is_ultra_ideal()

    outside = line(H3, -2, 1, 0, check=False)
    assert outside.is_ultra_ideal() and not outside.is_ideal()
    with pytest.raises(ValueError):
        line(H3, -2, 1, 0)


def test_pole(H3):
    assert line(H3, -QQ(1) / 3, 1, 0).pole() == point(H3, 9, 0, check=False)
    assert line(H3, -1, 0, 1).pole() == point(H3, 0, QQ(3) / 2, check=False)
    assert line(H3, -3, 1, 2, check=False).pole() == point(H3, 1, 1)


def test_perpendicular(H3):
    diameter = line(H3, 0, 1, -1)
    perpendicular = diameter.perpendicular(point(H3, QQ(1) / 2, QQ(1) / 2))
    assert perpendicular == line(H3, -QQ(3) / 2, 1, 2, oriented=False)

    origin = point(H3, 0, 0)
    assert line(H3, -QQ(1) / 3, 1, 0).perpendicular(origin) == line(
        H3, 0, 0, 1, oriented=False
    )


def test_angles(H3):
    x_third = line(H3, -QQ(1) / 3, 1, 0)
    x_axis = line(H3, 0, 0, 1)
    assert x_third.angle(x_axis, numerical=False) == QQ(1) / 4

    diagonal = line(H3, 0, 1, -1)
    other = line(H3, -QQ(3) / 2, 1, 2)
    assert diagonal.angle(other, numerical=False) == QQ(1) / 4
    assert other.angle(diagonal, numerical=False) == QQ(3) / 4

    import math

    angle = diagonal.angle(line(H3, 0, 1, 0))
    assert abs(math.cos(2 * math.pi * angle) ** 2 - 2 / 3) < 1e-12


def test_an_element(H3):
    assert line(H3, -QQ(3) / 2, 1, 2).an_element() == point(
        H3, QQ(1) / 2, QQ(1) / 2
    )


def test_endpoints(H3):
    g = line(H3, -1, 0, 1)
    assert g.start().coordinates(model="klein") == (-1, 1)
    assert g.end().coordinates(model="klein") == (1, 1)

    g = line(H3, -QQ(1) / 3, 1, 0)
    with pytest.raises(ValueError, match="square root"):
        g.start().coordinates(model="klein")

    x = polygen(QQ)
    K = NumberField(x**2 - 13, "s", embedding=AA(13).sqrt())
    s = K.gen()
    assert g.start().coordinates(model="klein", ring=K) == (QQ(1) / 3, s / 3)
    assert g.end().coordinates(model="klein", ring=K) == (QQ(1) / 3, -s / 3)

    diagonal = H3.geodesic(point(H3, 1, 1), point(H3, -1, -1))
    assert diagonal == line(H3, 0, 1, -1)
    assert diagonal.start() == point(H3, 1, 1)
    assert diagonal.end() == point(H3, -1, -1)


@pytest.mark.parametrize("form,x,y", [([1, 2, -3], 1, 1), ([1, 6, -10], -2, 1)])
def test_half_spaces_of_ideal_point(form, x, y):
    H = HyperbolicPlane(QQ, form=form)
    p = point(H, x, y)
    assert p.is_ideal()
    assert H.intersection(*p.half_spaces()) == p


def test_isometry(H3):
    boost = matrix(QQ, [[2, 0, 3], [0, 1, 0], [1, 0, 2]])
    assert point(H3, -1, 0).apply_isometry(boost, model="klein") == point(H3, 1, 0)
    assert point(H3, 0, 0).apply_isometry(boost, model="klein") == point(
        H3, QQ(3) / 2, 0
    )

    with pytest.raises(ValueError, match="does not define an isometry"):
        point(HyperbolicPlane(QQ), 0, 0).apply_isometry(boost, model="klein")


def test_guards(H3):
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3.infinity()
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3.point(0, 1, model="half_plane")
    with pytest.raises(ValueError, match=HALF_PLANE):
        point(H3, QQ(1) / 2, QQ(1) / 2).coordinates(model="half_plane")
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3(I)
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3(QQ(1))
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3.vertical(0)
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3.half_circle(0, 1)
    with pytest.raises(ValueError, match=HALF_PLANE):
        H3.isometry(point(H3, 0, 0), point(H3, 0, 0))
    with pytest.raises(ValueError, match=HALF_PLANE):
        point(H3, 0, 0).apply_isometry(matrix(QQ, [[1, 1], [0, 1]]))
    with pytest.raises(ValueError, match=HALF_PLANE):
        line(H3, -1, 0, 1).equation(model="half_plane")
    with pytest.raises(ValueError, match=HALF_PLANE):
        line(H3, -1, 0, 1).midpoint()


def test_plot(H3):
    segment = point(H3, 0, 0).segment(point(H3, QQ(1) / 2, QQ(1) / 2))
    geodesic = line(H3, -QQ(1) / 3, 1, 0)
    for subset in [segment, geodesic, geodesic.left_half_space()]:
        for model in ["half_plane", "klein"]:
            assert isinstance(subset.plot(model=model), Graphics)

    carried = H3._plot_point(point(H3, QQ(1) / 2, QQ(1) / 2))
    assert carried.parent() is HyperbolicPlane(RR)
    x, y = carried.coordinates(model="klein")
    assert abs(x - RR(0.288675134594813)) < 1e-12
    assert abs(y - RR(0.408248290463863)) < 1e-12


def test_inexact_geodesic_predicates_are_scaled():
    HR = HyperbolicPlane(RR, form=[1, 2, -3])
    scale = RR(1e-9)

    crossing = HR.geodesic(
        -scale / 3, scale, 0, model="klein", oriented=False
    )
    assert not crossing.is_ideal() and not crossing.is_ultra_ideal()

    tangent = HR.geodesic(-3 * scale, scale, 2 * scale, model="klein", check=False)
    assert tangent.is_ideal() and not tangent.is_ultra_ideal()

    outside = HR.geodesic(-2 * scale, scale, 0, model="klein", check=False)
    assert outside.is_ultra_ideal() and not outside.is_ideal()


def test_klein_ideal_point_small_disk():
    H = HyperbolicPlane(QQ, form=[1, QQ(3) / 4, -1])
    assert H._klein_ideal_point() == (QQ(1) / 2, 1)
    assert len(H.some_subsets()) == 20
    assert HyperbolicPlane(QQ, form=[1, 1, -3])._klein_ideal_point() is None
    assert HyperbolicPlane(QQ, form=[1, 2, -QQ(1) / 2])._klein_ideal_point() == (
        QQ(2) / 3,
        QQ(1) / 6,
    )


def test_conversion_between_rings_without_common_parent():
    x = polygen(QQ)
    K = NumberField(x**2 - 2, "a", embedding=RR(1.414))
    L = NumberField(x**2 - 3, "b", embedding=RR(1.732))
    HK = HyperbolicPlane(K, form=[1, 2, -3])
    HL = HyperbolicPlane(L, form=[1, 2, -3])
    p = HK.point(QQ(1) / 2, QQ(1) / 2, model="klein")

    assert HL(p) == HL.point(QQ(1) / 2, QQ(1) / 2, model="klein")
    assert HL(p).parent() is HL
    assert HK.change_ring(L) is HL

    with pytest.raises(ValueError, match="different forms"):
        HyperbolicPlane(L, form=[1, 2, -5])(p)

    HKa = HyperbolicPlane(K, form=[1, K.gen(), -3])
    with pytest.raises(ValueError, match="different forms"):
        HL(HKa.point(QQ(1) / 2, QQ(1) / 2, model="klein"))
    with pytest.raises(ValueError, match="the form does not lie in"):
        HKa.change_ring(L)


def test_endpoint_maybe_ring(H3):
    g = line(H3, -QQ(1) / 3, 1, 0)
    assert g.start()._coordinates_klein(ring="maybe") is None
    with pytest.raises(ValueError, match="square root of 52/9 not a rational"):
        g.start()._coordinates_klein(ring=None)


def test_form_with_leading_minus_sign_repr():
    x = polygen(QQ)
    K = NumberField(x**2 - 2, "a", embedding=-AA(2).sqrt())
    a = K.gen()
    prefix = "Hyperbolic Plane over Number Field in a with defining polynomial x^2 - 2 with a = -1.414213562373095? with form "
    assert repr(HyperbolicPlane(K, form=[1, -a, -3])) == prefix + "x^2 + (-a)*y^2 - 3*z^2"
    assert (
        repr(HyperbolicPlane(K, form=[1, -a, 2 * a]))
        == prefix + "x^2 + (-a)*y^2 - (-2*a)*z^2"
    )
    assert (
        repr(HyperbolicPlane(K, form=[1, -a / 2, -3]))
        == prefix + "x^2 + (-1/2*a)*y^2 - 3*z^2"
    )
