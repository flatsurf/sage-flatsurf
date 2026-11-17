r"""
Checks of surfaces glued from polygons of hyperbolic planes of the quadratic
forms `x^2 + 3y^2 - z^2` and `x^2 + 2y^2 - 3z^2` against hard-coded side
pairings, lifts, angles, areas and fields.
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

import itertools
import re

import pytest

import sage.misc.sage_unittest

from sage.all import AA, QQ, NumberField, diagonal_matrix, dumps, loads, matrix, polygen

from flatsurf import HyperbolicPlane, hyperbolic_surfaces
from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface


FORM31 = "x^2 + 3*y^2 - z^2"
FORM3 = "x^2 + 2*y^2 - 3*z^2"
STANDARD = "x^2 + y^2 - z^2"


def run_test_suite(parent):
    sage.misc.sage_unittest.TestSuite(parent).run(raise_on_failure=True)


@pytest.fixture
def H31():
    return HyperbolicPlane(QQ, form=[1, 3, -1])


@pytest.fixture
def H3():
    return HyperbolicPlane(QQ, form=[1, 2, -3])


def triangle(H31):
    return H31.polygon(
        vertices=[(QQ(1) / 2, 0), (-QQ(1) / 4, QQ(1) / 4), (-QQ(1) / 4, -QQ(1) / 4)]
    ).cyclic_quotient(3)


def kite(H3):
    return H3.polygon(
        vertices=[(0, 0), (QQ(1) / 2, -QQ(1) / 2), (QQ(3) / 4, 0), (QQ(1) / 2, QQ(1) / 2)]
    )


@pytest.fixture
def S3(H31):
    S = MutableOrientedHyperbolicSurface(H31)
    S.add_polygon(triangle(H31))
    S.glue((0, 0), (0, 2))
    S.glue((0, 1), (0, 1))
    S.set_immutable()
    return S


@pytest.fixture
def K3(H3):
    S = MutableOrientedHyperbolicSurface(H3)
    S.add_polygon(kite(H3))
    S.glue((0, 0), (0, 3))
    S.glue((0, 1), (0, 2))
    S.set_immutable()
    return S


def octagon_on_H3():
    r"""
    Return the mutable surface of the regular octagon of angle π/4, stretched
    into the disk `x^2 + 2y^2 < 3`, with opposite sides glued, over ``AA``.
    """
    x = polygen(QQ, "x")
    K = NumberField(
        x**4 + 32 * x**2 - 32, "r", embedding=AA(12 * AA(2).sqrt() - 16).sqrt()
    )
    r = K.gen()
    s2 = r**2 / 12 + QQ(4) / 3
    cos = [1, s2 / 2, 0, -s2 / 2, -1, -s2 / 2, 0, s2 / 2]
    sin = [0, s2 / 2, 1, s2 / 2, 0, -s2 / 2, -1, -s2 / 2]
    H = HyperbolicPlane(AA, form=[1, 2, -3])
    O3 = H.polygon(
        vertices=[
            (AA(3).sqrt() * AA(r * c), AA(QQ(3) / 2).sqrt() * AA(r * s))
            for c, s in zip(cos, sin)
        ]
    )
    S = MutableOrientedHyperbolicSurface(H)
    S.add_polygon(O3)
    for i in range(4):
        S.glue((0, i), (0, i + 4))
    return S


def test_sphere_on_H31(S3, H31):
    assert S3.hyperbolic_plane() is H31
    assert repr(S3.hyperbolic_plane()).endswith(FORM31)
    assert S3.polygon(0).vertices() == tuple(
        H31.point(QQ(x), QQ(y), model="klein")
        for x, y in [(0, 0), (QQ(1) / 2, 0), (-QQ(1) / 4, QQ(1) / 4)]
    )

    G31 = HyperbolicIsometryGroup(QQ, form=[1, 3, -1])
    rotation = matrix(QQ, [[-QQ(1) / 2, -QQ(3) / 2, 0], [QQ(1) / 2, -QQ(1) / 2, 0], [0, 0, 1]])
    half_turn = matrix(
        QQ,
        [
            [-QQ(31) / 30, -QQ(1) / 10, QQ(4) / 15],
            [-QQ(1) / 30, -QQ(11) / 10, QQ(4) / 15],
            [-QQ(4) / 15, -QQ(4) / 5, QQ(17) / 15],
        ],
    )
    assert S3.edge_transformation(0, 0) == G31(rotation)
    assert S3.edge_transformation(0, 2) == G31(rotation.inverse())
    assert S3.edge_transformation(0, 1) == G31(half_turn)
    assert S3.edge_transformation(0, 1).trace() == -1
    assert (S3.edge_transformation(0, 0) ** 3).is_one()

    phi = S3.holonomy_representation()
    assert phi.codomain() is G31
    assert repr(phi.codomain()).endswith(FORM31)


def test_sphere_on_H31_angles_and_area(S3):
    assert len(S3.vertices()) == 2
    assert S3(0, 0).angle() == QQ(1) / 3
    assert S3(0, 1) == S3(0, 2)
    assert S3(0, 1).angle(numerical=True) == pytest.approx(0.1475836176504333)
    with pytest.raises(ValueError, match="not a rational multiple"):
        S3(0, 1).angle()
    with pytest.raises(NotImplementedError, match="not a rational multiple"):
        S3.signature()
    with pytest.warns(UserWarning, match="non-rational polygons"):
        assert S3.area() == pytest.approx(0.019083049016233367)


def test_sphere_on_H31_fields(S3):
    assert S3.holonomy_representation().invariant_trace_field() is QQ
    run_test_suite(S3)


def test_kite_on_H3(K3, H3):
    G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
    g = K3.edge_transformation(0, 0)
    assert g.parent() is G3
    assert g.matrix() == matrix(
        QQ, [[-QQ(1) / 3, -QQ(4) / 3, 0], [QQ(2) / 3, -QQ(1) / 3, 0], [0, 0, 1]]
    )
    assert g.trace() == QQ(1) / 3
    h = K3.edge_transformation(0, 1)
    assert h.matrix() == matrix(
        QQ,
        [
            [-QQ(17) / 15, QQ(16) / 15, QQ(8) / 5],
            [-QQ(8) / 15, -QQ(11) / 15, QQ(2) / 5],
            [-QQ(8) / 15, QQ(4) / 15, QQ(7) / 5],
        ],
    )
    assert h.trace() == -QQ(7) / 15
    assert K3.polygon(0)._cosh2_sides() == (
        QQ(4) / 3,
        QQ(49) / 39,
        QQ(49) / 39,
        QQ(4) / 3,
    )
    assert K3.hyperbolic_plane() is H3
    with pytest.warns(UserWarning, match="non-rational polygons"):
        assert K3.area() == pytest.approx(0.032847912218780584)
    assert sorted(K3.angles(numerical=True)) == pytest.approx(
        [0.2820471084244875, 0.3040867239846964, 0.38101825537203543]
    )
    assert K3.holonomy_representation().invariant_trace_field() is QQ
    run_test_suite(K3)


def test_octagon_on_H3():
    S = octagon_on_H3()
    assert all(
        S.edge_transformation(0, i).trace() == 8 * AA(2).sqrt() + 11 for i in range(4)
    )
    S.set_immutable()
    L = S.holonomy_representation().invariant_trace_field()
    assert L.polynomial() == polygen(QQ, "y") ** 2 - 2
    assert S.angles() == [1]
    assert S.area() == 2


def test_octagon_on_H3_minimal_field():
    S = octagon_on_H3()
    S.set_immutable(minimal_field=True)
    K = S.base_ring()
    assert K.polynomial() == polygen(QQ, "y") ** 4 - 18
    a = K.gen()
    assert S.hyperbolic_plane() is HyperbolicPlane(K, form=[1, 2, -3])
    assert S.polygon(0).parent() is S.hyperbolic_plane()
    assert all(
        S.edge_transformation(0, i).trace() == 8 * a**2 / 3 + 11 for i in range(4)
    )
    L = S.holonomy_representation().invariant_trace_field()
    assert L.polynomial() == polygen(QQ, "x") ** 2 - 2
    run_test_suite(S)


def test_polygon_of_other_form_refused(H3):
    P = kite(H3)
    S = MutableOrientedHyperbolicSurface(QQ)
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"polygon lies in a hyperbolic plane with form {FORM3}, the surface in a hyperbolic plane with form {STANDARD}"
        ),
    ):
        S.add_polygon(P)

    T = MutableOrientedHyperbolicSurface(H3)
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"polygon lies in a hyperbolic plane with form {STANDARD}, the surface in a hyperbolic plane with form {FORM3}"
        ),
    ):
        T.add_polygon(HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]))
    with pytest.raises(ValueError, match="polygon lies in a hyperbolic plane"):
        T.replace_polygon(T.add_polygon(P), HyperbolicPlane(QQ).polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]))


def test_plane_given_by_ring_or_plane(H3):
    S = MutableOrientedHyperbolicSurface(QQ)
    assert S == MutableOrientedHyperbolicSurface(HyperbolicPlane(QQ))
    assert S != MutableOrientedHyperbolicSurface(H3)
    assert S.hyperbolic_plane() is HyperbolicPlane(QQ)
    assert MutableOrientedHyperbolicSurface(H3).base_ring() is QQ


def test_equality_and_pickling(S3, K3, H31):
    assert loads(dumps(S3)) == S3
    assert loads(dumps(S3)).hyperbolic_plane() is H31
    assert hash(loads(dumps(S3))) == hash(S3)
    assert S3 != K3


def test_change_ring(S3, K3):
    T = S3.change_ring(AA)
    assert T.hyperbolic_plane() is HyperbolicPlane(AA, form=[1, 3, -1])
    assert T.edge_transformation(0, 0).parent() is HyperbolicIsometryGroup(
        AA, form=[1, 3, -1]
    )
    assert T.change_ring(QQ) == S3
    assert K3.change_ring(AA).hyperbolic_plane().quadratic_form() == diagonal_matrix(
        AA, [1, 2, -3]
    )


def test_lifts_on_G3():
    G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
    one, half_turn, rotation, hyperbolic = G3.some_elements()
    assert half_turn.matrix() == diagonal_matrix(QQ, [-1, -1, 1])

    x, s = one._half_lift()
    assert (x, s) == (4, 16)
    i, j, k = x.parent().gens()
    assert half_turn._half_lift() == (-3 * i, 18)
    assert rotation._half_lift() == (QQ(4) / 3 + QQ(4) / 3 * i, QQ(16) / 3)
    assert hyperbolic._half_lift() == (QQ(64) / 13 + QQ(8) / 13 * k, QQ(256) / 13)


def test_lifts_of_side_pairings(S3, K3):
    x, s = S3.edge_transformation(0, 0)._half_lift()
    i, j, k = x.parent().gens()
    assert (x, s) == (1 - i, 4)
    assert S3.edge_transformation(0, 2)._half_lift() == (1 + i, 4)
    assert S3.edge_transformation(0, 1)._half_lift() == (-8 * i + j - k, 180)

    x, s = K3.edge_transformation(0, 0)._half_lift()
    i, j, k = x.parent().gens()
    assert (x, s) == (QQ(4) / 3 - QQ(4) / 3 * i, QQ(16) / 3)
    assert K3.edge_transformation(0, 1)._half_lift() == (
        QQ(8) / 15 + QQ(16) / 15 * i - QQ(4) / 15 * j,
        QQ(32) / 15,
    )


def test_pillowcase_on_H31(H31):
    C = triangle(H31)
    S = hyperbolic_surfaces.pillowcase(C, mirror_edge=1)
    assert S.hyperbolic_plane() is H31
    assert S.polygon(1).vertices() == tuple(
        H31.point(x, y, model="klein")
        for x, y in [(QQ(4) / 17, QQ(4) / 17), (-QQ(1) / 4, QQ(1) / 4), (QQ(1) / 2, 0)]
    )
    assert S.polygon(1)._cosh2_sides() == (QQ(4) / 3, QQ(9) / 4, QQ(4) / 3)
    assert C._cosh2_sides() == S.polygon(1)._cosh2_sides()
    assert sorted(S.angles(numerical=True)) == pytest.approx(
        [0.1475836176504333, 0.1475836176504333, 0.6666666666666667]
    )
    with pytest.raises(NotImplementedError, match="not a rational multiple"):
        S.signature()
    assert S.holonomy_representation().invariant_trace_field() is QQ
    run_test_suite(S)

    # The default mirror, the edge 0 on the x-axis, is the reflection
    # (x, y) -> (x, -y) whatever the form.
    T = hyperbolic_surfaces.pillowcase(C)
    assert T.polygon(1).vertices() == tuple(
        H31.point(x, y, model="klein")
        for x, y in [(QQ(1) / 2, 0), (0, 0), (-QQ(1) / 4, -QQ(1) / 4)]
    )


@pytest.mark.parametrize("b, c", [(2, 3), (3, 1), (1, 1), (6, 10)])
def test_lift_is_a_homomorphism(b, c):
    G = HyperbolicIsometryGroup(QQ, form=[1, b, -c])
    elements = G.some_elements()
    assert len(elements) <= 5

    def trace(x):
        return x.trace() if G._form is None else x.reduced_trace()

    def is_scalar(x):
        if G._form is None:
            return x.is_scalar(x[0, 0])
        return all(not e for e in x.coefficient_tuple()[1:])

    for size in (1, 2, 3):
        for word in itertools.product(elements, repeat=size):
            g = G.one()
            lift = 1
            for h in word:
                g = g * h
                lift = lift * h._half_lift()[0]
            x, s = g._half_lift()
            assert trace(x) ** 2 / s == g.trace() + 1
            assert is_scalar(x * ~lift)
