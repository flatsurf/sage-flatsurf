r"""
Pillowcases of flat and hyperbolic polygons: the gluings, the boundary, the
base ring and the refused arguments, against hard-coded values.

The gluings follow from the reflection in the line of the mirror edge `m`
of a polygon with `n` vertices: the vertex `j` of the mirror image is the
reflection of the vertex `2m + 1 - j mod n`, and the edge `e` is glued to
the edge `2m - e mod n` of the mirror image. The boundary components, base
rings and Klein coordinates below were computed with sage-flatsurf.
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

from sage.all import QQ, NumberField, polygen

from flatsurf import (
    HyperbolicPlane,
    Polygon,
    hyperbolic_polygons,
    hyperbolic_surfaces,
    polygons,
    similarity_surfaces,
)
from flatsurf.geometry.similarity_surface_generators import _pillowcase_gluings


def hyperbolic_triangle():
    return HyperbolicPlane(QQ).polygon(
        vertices=[(0, 0), (QQ(1) / 2, 0), (QQ(1) / 3, QQ(1) / 2)]
    )


def hyperbolic_hexagon():
    return hyperbolic_polygons.right_angled_hexagon(
        QQ(5) / 4, QQ(5) / 4, QQ(145) / 24
    )


def flat_sqrt2_triangle():
    x = polygen(QQ)
    K = NumberField(x**2 - 2, "a", embedding=QQ(14) / 10)
    return Polygon(vertices=[(0, 0), (1, K.gen()), (-1, 1)])


L_SHAPE = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]


@pytest.mark.parametrize(
    "n,glued_edges,mirror_edge,expected",
    [
        (
            5,
            None,
            2,
            (
                2,
                [0, 4, 3, 2, 1],
                [
                    ((0, 0), (1, 4)),
                    ((0, 1), (1, 3)),
                    ((0, 2), (1, 2)),
                    ((0, 3), (1, 1)),
                    ((0, 4), (1, 0)),
                ],
            ),
        ),
        (5, [3, 1], None, (3, [2, 1, 0, 4, 3], [((0, 3), (1, 3)), ((0, 1), (1, 0))])),
        (
            4,
            None,
            None,
            (
                0,
                [1, 0, 3, 2],
                [
                    ((0, 0), (1, 0)),
                    ((0, 1), (1, 3)),
                    ((0, 2), (1, 2)),
                    ((0, 3), (1, 1)),
                ],
            ),
        ),
    ],
)
def test_gluings(n, glued_edges, mirror_edge, expected):
    assert _pillowcase_gluings(n, glued_edges, mirror_edge) == expected


@pytest.mark.parametrize(
    "glued_edges,mirror_edge,message",
    [
        ([], None, "must not be empty"),
        ([0, 2], 1, "mirror_edge 1 is not one of the glued edges"),
        ([0, 2, 0], None, "edge 0 is repeated"),
        ([6], None, "6 is not an edge of a polygon with 6 vertices"),
        (None, 6, "6 is not an edge of a polygon with 6 vertices"),
        (None, -1, "-1 is not an edge of a polygon with 6 vertices"),
    ],
)
def test_refused_arguments(glued_edges, mirror_edge, message):
    with pytest.raises(ValueError, match=message):
        hyperbolic_surfaces.pillowcase(
            hyperbolic_hexagon(), glued_edges=glued_edges, mirror_edge=mirror_edge
        )
    with pytest.raises(ValueError, match=message):
        similarity_surfaces.pillowcase(
            Polygon(vertices=L_SHAPE), glued_edges=glued_edges, mirror_edge=mirror_edge
        )


def test_refused_geometry():
    with pytest.raises(TypeError, match="Euclidean polygon"):
        similarity_surfaces.pillowcase(hyperbolic_triangle())
    with pytest.raises(TypeError, match="hyperbolic polygon"):
        hyperbolic_surfaces.pillowcase(polygons.square())


def test_hyperbolic_triangle():
    S = hyperbolic_surfaces.pillowcase(hyperbolic_triangle(), mirror_edge=1)
    assert S.base_ring() is QQ
    assert not S.is_mutable()
    assert [v.coordinates(model="klein") for v in S.polygon(1).vertices()] == [
        (QQ(36) / 49, QQ(12) / 49),
        (QQ(1) / 3, QQ(1) / 2),
        (QQ(1) / 2, 0),
    ]
    assert [S.opposite_edge(0, e) for e in range(3)] == [(1, 2), (1, 1), (1, 0)]
    assert S.boundary_components() == []


def test_hyperbolic_hexagon():
    S = hyperbolic_surfaces.pillowcase(hyperbolic_hexagon(), glued_edges=[0, 2, 4])
    assert S.base_ring() is QQ
    assert S.boundary_components() == [
        ((0, 1), (1, 5)),
        ((0, 3), (1, 3)),
        ((0, 5), (1, 1)),
    ]


def test_hyperbolic_non_convex():
    half = QQ(1) / 2
    P = HyperbolicPlane(QQ).polygon(
        vertices=[(0, -half), (half, half), (0, 0), (-half, half)]
    )
    S = hyperbolic_surfaces.pillowcase(P)
    assert S.base_ring() is QQ
    assert S(0, 2).angle() == QQ(3) / 2


def test_flat_sqrt2_triangle():
    P = flat_sqrt2_triangle()
    S = similarity_surfaces.pillowcase(P)
    assert S.base_ring() is P.base_ring()
    a = P.base_ring().gen()
    assert tuple(tuple(v) for v in S.polygon(1).vertices()) == (
        (1, a),
        (0, 0),
        (2 * a / 3 + QQ(1) / 3, -2 * a / 3 + QQ(1) / 3),
    )
    assert S.boundary_components() == []


def test_flat_square_with_boundary():
    S = similarity_surfaces.pillowcase(polygons.square(), glued_edges=[0, 1])
    assert S.base_ring() is QQ
    assert S.boundary_components() == [((0, 2), (0, 3), (1, 1), (1, 2))]


def test_flat_non_convex():
    S = similarity_surfaces.pillowcase(Polygon(vertices=L_SHAPE))
    assert list(S.labels()) == [0, 1]
    assert S.boundary_components() == []


@pytest.mark.parametrize(
    "construct,mirror_edge",
    [
        (
            lambda: hyperbolic_surfaces.pillowcase(hyperbolic_triangle(), mirror_edge=1),
            1,
        ),
        (
            lambda: hyperbolic_surfaces.pillowcase(
                hyperbolic_hexagon(), glued_edges=[4, 0, 2]
            ),
            4,
        ),
        (lambda: similarity_surfaces.pillowcase(flat_sqrt2_triangle()), 0),
        (
            lambda: similarity_surfaces.pillowcase(
                Polygon(vertices=L_SHAPE), glued_edges=[1, 3], mirror_edge=3
            ),
            3,
        ),
    ],
)
def test_mirror_edge_identity(construct, mirror_edge):
    S = construct()
    assert S.opposite_edge(0, mirror_edge) == (1, mirror_edge)
    assert S.edge_transformation(0, mirror_edge).is_one()
