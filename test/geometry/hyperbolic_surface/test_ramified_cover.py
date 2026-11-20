r"""
Ramified covers of flat and hyperbolic surfaces: the area, the cone orders,
the invariant trace field and the plane of the covers, the self-glued and
unglued edges and the refused data, against hard-coded values.

The area of a cover of degree `d` is `d` times the area of the surface. The
cone orders and the gluings below were computed with sage-flatsurf.
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

from sage.all import QQ, SymmetricGroup

from flatsurf import (
    HyperbolicPlane,
    MutableOrientedSimilaritySurface,
    euclidean_polygons,
    hyperbolic_polygons,
    hyperbolic_surfaces,
    polygons,
    similarity_surfaces,
)
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface


def hyperbolic_pillowcase():
    return hyperbolic_surfaces.pillowcase(
        hyperbolic_polygons.triangle(QQ(1) / 8, QQ(1) / 8, QQ(1) / 8)
    )


def kite():
    H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
    S = MutableOrientedHyperbolicSurface(H3)
    S.add_polygon(
        H3.polygon(
            vertices=[
                (0, 0),
                (QQ(1) / 2, -QQ(1) / 2),
                (QQ(3) / 4, 0),
                (QQ(1) / 2, QQ(1) / 2),
            ]
        )
    )
    S.glue((0, 0), (0, 3))
    S.glue((0, 1), (0, 2))
    S.set_immutable()
    return S


def square(gluings):
    S = MutableOrientedSimilaritySurface(QQ)
    S.add_polygon(polygons.square())
    for x, y in gluings:
        S.glue(x, y)
    S.set_immutable()
    return S


def test_flat_pillowcase():
    S = similarity_surfaces.pillowcase(euclidean_polygons.triangle(1, 1, 1))
    C = S.ramified_cover(2, {(0, 0): "(1,2)"})
    assert C.is_connected()
    assert len(C.polygons()) == 4
    assert C.signature().cone_orders() == (
        -QQ(2) / 3,
        -QQ(2) / 3,
        -QQ(1) / 3,
        -QQ(1) / 3,
    )
    assert C.area() == 2 * S.area()


@pytest.mark.parametrize(
    "d,data,genus,area,cone_orders,orbifold",
    [
        (
            2,
            {(0, 0): "(1,2)"},
            0,
            QQ(1) / 2,
            [-QQ(3) / 4, -QQ(3) / 4, -QQ(1) / 2, -QQ(1) / 2],
            True,
        ),
        (
            3,
            {(0, 0): "(1,2,3)", (0, 1): "(1,3,2)"},
            1,
            QQ(3) / 4,
            [-QQ(1) / 4] * 3,
            False,
        ),
        (
            4,
            {(0, 0): "(1,2)(3,4)", (0, 1): "(1,3)(2,4)"},
            0,
            QQ(1),
            [-QQ(1) / 2] * 6,
            True,
        ),
    ],
)
def test_hyperbolic_pillowcase(d, data, genus, area, cone_orders, orbifold):
    H = hyperbolic_pillowcase()
    C = H.ramified_cover(d, data)
    assert not C.is_mutable()
    assert len(C.polygons()) == 2 * d
    assert C.area() == area
    signature = C.signature()
    assert signature.genus() == genus
    assert sorted(signature.cone_orders()) == sorted(cone_orders)
    assert C.is_orbifold() == orbifold
    assert C.hyperbolic_plane() == H.hyperbolic_plane()
    field = C.holonomy_representation().invariant_trace_field()
    assert field.is_isomorphic(H.holonomy_representation().invariant_trace_field())
    assert field.degree() == 2


def test_quadratic_form_kept():
    # The angles of the kite are not rational multiples of 2π, so it has
    # no exact area or signature. The kite is a sphere with three vertices,
    # and its double cover is branched over two of them: by Riemann-Hurwitz
    # the cover is again a sphere, with four vertices.
    S = kite()
    C = S.ramified_cover(2, {(0, 0): "(1,2)"})
    assert C.hyperbolic_plane() == S.hyperbolic_plane()
    assert C.hyperbolic_plane() != HyperbolicPlane(QQ)
    assert len(C.polygons()) == 2
    assert len(S.vertices()) == 3
    assert len(C.vertices()) == 4
    assert C.euler_characteristic() == 2
    assert S.ramified_cover(1, {}).hyperbolic_plane() == S.hyperbolic_plane()


def test_self_glued_edge():
    S = square([((0, 0), (0, 2)), ((0, 1), (0, 1)), ((0, 3), (0, 3))])
    C = S.ramified_cover(2, {(0, 0): "(1,2)"})
    assert [C.opposite_edge(label, e) for label in C.labels() for e in range(4)] == [
        (1, 2),
        (0, 1),
        (1, 0),
        (0, 3),
        (0, 2),
        (1, 1),
        (0, 0),
        (1, 3),
    ]
    C = S.ramified_cover(2, {(0, 1): SymmetricGroup(2)("(1,2)")})
    assert [C.opposite_edge(label, 1) for label in C.labels()] == [(1, 1), (0, 1)]


def test_non_involution():
    S = square([((0, 0), (0, 2)), ((0, 1), (0, 1)), ((0, 3), (0, 3))])
    with pytest.raises(ValueError, match="not an involution"):
        S.ramified_cover(3, {(0, 1): "(1,2,3)"})


def test_unglued_edge():
    S = square([((0, 0), (0, 2))])
    C = S.ramified_cover(2, {(0, 0): "(1,2)"})
    assert C.boundary_components() == [((0, 1), (1, 1)), ((0, 3), (1, 3))]


def test_unknown_key():
    S = square([((0, 0), (0, 2)), ((0, 1), (0, 3))])
    with pytest.raises(ValueError, match="not an edge of this surface"):
        S.ramified_cover(2, {(1, 0): "(1,2)"})


def test_inconsistent_data():
    S = square([((0, 0), (0, 2)), ((0, 1), (0, 3))])
    with pytest.raises(ValueError, match="inconsistent covering data"):
        S.ramified_cover(3, {(0, 0): "(1,2)", (0, 2): "(1,3)"})


def test_data_unchanged():
    H = hyperbolic_pillowcase()
    data = {(0, 0): "(1,2,3)", (0, 1): "(1,3,2)"}
    H.ramified_cover(3, data)
    assert data == {(0, 0): "(1,2,3)", (0, 1): "(1,3,2)"}
    assert all(isinstance(value, str) for value in data.values())


@pytest.mark.parametrize(
    "construct",
    [
        hyperbolic_pillowcase,
        lambda: similarity_surfaces.pillowcase(euclidean_polygons.triangle(1, 1, 1)),
    ],
)
def test_degree_one(construct):
    S = construct()
    C = S.ramified_cover(1, {})
    assert C.area() == S.area()
    assert C.signature() == S.signature()
