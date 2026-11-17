r"""
Checks of the group of isometries of the quadratic form `x^2 + 2y^2 - 3z^2`
and of the polygons of the hyperbolic plane of that form against hard-coded
matrices, lengths and angles.
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

from sage.all import AA, QQ, I, diagonal_matrix, dumps, loads, matrix

from flatsurf import HyperbolicPlane
from flatsurf.geometry.hyperbolic_isometry import (
    HyperbolicIsometryGroup,
    segment_isometry,
)


HALF_PLANE = re.escape("the half-plane model needs the form x^2 + y^2 - z^2")
FORM3 = "x^2 + 2*y^2 - 3*z^2"
OTHER_FORM = re.escape(
    "the isometry preserves x^2 + 2*y^2 - 3*z^2 and the plane has the form x^2 + y^2 - z^2"
)


@pytest.fixture
def H3():
    return HyperbolicPlane(QQ, form=[1, 2, -3])


@pytest.fixture
def G3():
    return HyperbolicIsometryGroup(QQ, form=[1, 2, -3])


@pytest.fixture
def B():
    return matrix(QQ, [[2, 0, 3], [0, 1, 0], [1, 0, 2]])


def run_test_suite(parent):
    sage.misc.sage_unittest.TestSuite(parent).run(raise_on_failure=True)


def point(H, x, y):
    return H.point(QQ(x), QQ(y), model="klein")


def test_normalisation_and_printing(G3, H3):
    assert HyperbolicIsometryGroup(QQ, form=[2, 2, -2]) is HyperbolicIsometryGroup(QQ)
    assert HyperbolicIsometryGroup(QQ, form=[2, 4, -6]) is G3
    assert loads(dumps(G3)) is G3
    assert G3.quadratic_form() == diagonal_matrix(QQ, [1, 2, -3])
    assert HyperbolicIsometryGroup(QQ).quadratic_form() == diagonal_matrix(
        QQ, [1, 1, -1]
    )
    assert G3.hyperbolic_plane() is H3
    assert HyperbolicIsometryGroup(QQ).hyperbolic_plane() is HyperbolicPlane(QQ)
    assert (
        repr(G3)
        == "Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form "
        + FORM3
    )
    assert (
        repr(HyperbolicIsometryGroup(QQ))
        == "Group of orientation preserving isometries of the hyperbolic plane over Rational Field"
    )


def test_test_suites(G3):
    run_test_suite(G3)
    run_test_suite(HyperbolicIsometryGroup(QQ, form=[1, 3, -1]))


def test_boost_and_inverse(G3, B):
    g = G3(B)
    assert g.trace() == 5
    assert (~g).matrix() == matrix(QQ, [[2, 0, -3], [0, 1, 0], [-1, 0, 2]])
    assert (~g * g).is_one()


def test_refusals(G3, B):
    G = HyperbolicIsometryGroup(QQ)
    with pytest.raises(ValueError, match=re.escape("matrix is not in SO(2, 1)+")):
        G(B)
    with pytest.raises(ValueError, match=re.escape(f"matrix is not in SO({FORM3})+")):
        G3([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    with pytest.raises(ValueError, match=HALF_PLANE):
        G3(matrix(QQ, [[2, 0], [0, 1]]))
    with pytest.raises(ValueError, match=re.escape(f"matrix is not in SO({FORM3})+")):
        G3(G.an_element())
    with pytest.raises(ValueError, match=re.escape("matrix is not in SO(2, 1)+")):
        G(G3(B))


def test_some_elements(G3):
    assert [g.trace() for g in G3.some_elements()] == [3, -1, QQ(1) / 3, QQ(51) / 13]

    rotation = G3.some_elements()[2]
    assert rotation.matrix() == matrix(
        QQ, [[-QQ(1) / 3, QQ(4) / 3, 0], [-QQ(2) / 3, -QQ(1) / 3, 0], [0, 0, 1]]
    )
    centre = rotation.fixed_point()
    assert centre.parent() is G3.hyperbolic_plane()
    assert centre.coordinates(model="klein") == (0, 0)

    assert G3.an_element().matrix() == matrix(
        QQ,
        [
            [QQ(19) / 13, 0, -QQ(24) / 13],
            [0, 1, 0],
            [-QQ(8) / 13, 0, QQ(19) / 13],
        ],
    )


def test_pushout_and_coercion(G3):
    G3AA = HyperbolicIsometryGroup(AA, form=[1, 2, -3])
    assert (G3.one() * G3AA.one()).parent() is G3AA
    assert G3AA.has_coerce_map_from(G3)
    assert not HyperbolicIsometryGroup(AA).has_coerce_map_from(G3)

    with pytest.raises(TypeError, match="unsupported operand parent"):
        G3.one() * HyperbolicIsometryGroup(QQ).one()

    F, K = G3.construction()
    assert F(K) is G3
    assert F(AA) is G3AA
    assert F != HyperbolicIsometryGroup(QQ).construction()[0]


def test_action(G3, H3, B):
    assert G3(B) * point(H3, 0, 0) == point(H3, QQ(3) / 2, 0)
    with pytest.raises(TypeError, match=OTHER_FORM):
        G3(B) * HyperbolicPlane(QQ)(I)
    with pytest.raises(TypeError, match=OTHER_FORM):
        G3(B)(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]))


def test_lift_to_sl2_refused(G3, B):
    with pytest.raises(ValueError, match=HALF_PLANE):
        G3(B).psl2_matrix()
    with pytest.raises(NotImplementedError, match=re.escape("needs the form x^2 + y^2 - z^2")):
        G3(B)._half_lift()


def test_segment_isometry_on_points(H3):
    p0 = point(H3, 0, 0)
    p1 = point(H3, 1, 0)
    p2 = point(H3, QQ(1) / 3, QQ(2) / 3)

    with pytest.raises(ValueError, match="segments do not have the same length"):
        segment_isometry((p0, point(H3, QQ(1) / 2, 0)), (p0, point(H3, 0, QQ(1) / 2)))

    rotation = matrix(QQ, [[QQ(1) / 3, -QQ(4) / 3, 0], [QQ(2) / 3, QQ(1) / 3, 0], [0, 0, 1]])
    assert segment_isometry((p0, p1), (p0, p2)) == rotation
    assert (
        segment_isometry(
            ((0, 0), (1, 0)), ((0, 0), (QQ(1) / 3, QQ(2) / 3)), form=H3.quadratic_form()
        )
        == rotation
    )

    H = HyperbolicPlane(QQ)
    with pytest.raises(ValueError, match="different forms"):
        segment_isometry((p0, p1), (H(I), H(2 * I)))
    with pytest.raises(ValueError, match="the form differs"):
        segment_isometry((p0, p1), (p0, p2), form=[1, 1, -1])


def test_lambert_quadrilateral(H3):
    L = H3.polygon(vertices=[(0, 0), (1, 0), (1, QQ(1) / 2), (0, QQ(1) / 2)])
    assert [L.edge_cosh_squared(i) for i in range(4)] == [
        QQ(3) / 2,
        QQ(4) / 3,
        QQ(5) / 3,
        QQ(6) / 5,
    ]
    angles = L.angles(numerical=True)
    for angle, expected in zip(angles, (0.25, 0.25, 0.198791808825217, 0.25)):
        assert abs(angle - expected) < 1e-12
    area = L.area(numerical=True)
    assert abs(area - 0.0512081911747834) < 1e-12
    assert abs(H3.convex_hull(*L.vertices()).area() - area) < 1e-12
    assert not L.is_rational()
    run_test_suite(L)


def test_square(H3, G3, B):
    W = H3.polygon(
        vertices=[(0, 0), (QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]
    )
    assert W._cosh2_sides() == (QQ(12) / 11, QQ(11) / 9, QQ(10) / 9, QQ(6) / 5)
    angles = W.angles(numerical=True)
    for angle, expected in zip(angles, (0.25, 0.25, 0.228473982282504, 0.25)):
        assert abs(angle - expected) < 1e-12
    run_test_suite(W)

    assert W.change_ring(AA).parent() is HyperbolicPlane(AA, form=[1, 2, -3])

    image = W.apply_isometry(B)
    assert image.vertices() == tuple(
        point(H3, x, y)
        for (x, y) in [
            (QQ(3) / 2, 0),
            (QQ(8) / 5, 0),
            (QQ(8) / 5, QQ(1) / 5),
            (QQ(3) / 2, QQ(1) / 4),
        ]
    )
    assert image == G3(B)(W) == G3(B) * W
    assert W._transformation_group() is G3

    with pytest.raises(ValueError, match="matrix does not define an isometry"):
        W.apply_isometry(HyperbolicIsometryGroup(QQ).an_element())
    with pytest.raises(ValueError, match=HALF_PLANE):
        W.apply_isometry(matrix(QQ, [[2, 0], [0, 1]]), model="half_plane")


def test_vertex_outside_disk(H3):
    with pytest.raises(
        ValueError,
        match=re.escape("vertex (2, 0) is not a finite point of the hyperbolic plane"),
    ):
        H3.polygon(vertices=[(0, 0), (2, 0), (0, QQ(1) / 2)])


def test_cyclic_quotient():
    H31 = HyperbolicPlane(QQ, form=[1, 3, -1])
    vertices = [(QQ(1) / 2, 0), (-QQ(1) / 4, QQ(1) / 4), (-QQ(1) / 4, -QQ(1) / 4)]
    T = H31.polygon(vertices=vertices)
    assert T._cosh2_sides() == (QQ(9) / 4,) * 3

    Q = T.cyclic_quotient(3)
    assert Q == H31.polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (-QQ(1) / 4, QQ(1) / 4)])
    angles = Q.angles(numerical=True)
    for angle, expected in zip(
        angles, (QQ(1) / 3, 0.0737918088252167, 0.0737918088252167)
    ):
        assert abs(angle - expected) < 1e-12

    with pytest.raises(ValueError, match="not invariant under a rotation of order 3"):
        HyperbolicPlane(QQ).polygon(vertices=vertices).cyclic_quotient(3)


def test_cyclic_quotient_half_turn(H3):
    R = H3.polygon(
        vertices=[
            (-QQ(1) / 2, -QQ(1) / 4),
            (QQ(1) / 2, -QQ(1) / 4),
            (QQ(1) / 2, QQ(1) / 4),
            (-QQ(1) / 2, QQ(1) / 4),
        ]
    )
    assert R.cyclic_quotient(2) == H3.polygon(
        vertices=[(-QQ(1) / 2, -QQ(1) / 4), (QQ(1) / 2, -QQ(1) / 4), (QQ(1) / 2, QQ(1) / 4)]
    )
