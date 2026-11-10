r"""
Checks of surfaces glued from compact hyperbolic polygons against hard-coded
values of classical hyperbolic geometry.
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

from sage.all import QQ, RR, NumberField, matrix, polygen, vector

from flatsurf import HyperbolicPlane, Polygon
from flatsurf.geometry.categories import ConeSurfaces, HyperbolicConeSurfaces
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface

from hyperbolic_surface_oracles import OCTAGON, OCTAGON_GLUINGS, OCTAGON_SIDE_COSH2


# Klein coordinates of a Lambert quadrilateral (three right angles).
QUADRILATERAL = [(0, 0), (QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]


@pytest.fixture(scope="module")
def H():
    return HyperbolicPlane(QQ)


def quadrilateral_surface(H, gluings):
    S = MutableOrientedHyperbolicSurface(QQ)
    S.add_polygon(H.polygon(vertices=QUADRILATERAL))
    for a, b in gluings:
        S.glue((0, a), (0, b))
    return S


@pytest.fixture(scope="module")
def sphere(H):
    S = quadrilateral_surface(H, [(0, 3), (1, 2)])
    S.set_immutable()
    return S


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


def in_radius(r, coefficients):
    return sum(c * r**i for i, c in enumerate(coefficients))


def test_sphere(sphere):
    assert sphere.category().is_subcategory(HyperbolicConeSurfaces().Oriented())
    assert sphere in ConeSurfaces().FiniteType()
    assert sphere.genus() == 0
    assert sphere.euler_characteristic() == 2
    assert len(sphere.vertices()) == 3
    assert not sphere.is_mutable()


def test_sphere_angles(sphere):
    assert sphere(0, 0).angle() == QQ(1) / 4
    assert sphere(0, 1).angle() == QQ(1) / 2
    assert sphere(0, 3) == sphere(0, 1)
    assert abs(sphere(0, 2).angle(numerical=True) - 0.195913276015304) < 1e-9
    with pytest.raises(ValueError):
        sphere(0, 2).angle()
    expected = [0.195913276015304, 0.25, 0.5]
    for angle, value in zip(sorted(sphere.angles(numerical=True)), expected):
        assert abs(angle - value) < 1e-9


def test_sphere_pairings(sphere):
    assert sphere.edge_matrix(0, 0) == matrix(QQ, [[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert sphere.edge_matrix(0, 1) == matrix(
        QQ,
        [
            [0, 1, 0],
            [-QQ(5) / 3, 0, QQ(4) / 3],
            [-QQ(4) / 3, 0, QQ(5) / 3],
        ],
    )
    g = sphere.edge_transformation(0, 1)
    assert g * vector(QQ, (QQ(1) / 2, 0)) == vector(QQ, (0, QQ(1) / 2))
    assert g * vector(QQ, (QQ(1) / 2, QQ(1) / 2)) == vector(
        QQ, (QQ(1) / 2, QQ(1) / 2)
    )
    for label, edge in sphere.edges():
        cross = sphere.opposite_edge(label, edge)
        g = sphere.edge_transformation(label, edge)
        assert g * sphere.edge_transformation(*cross) == g.parent().one()


def test_mismatched_lengths(H):
    with pytest.raises(ValueError, match="do not have the same length"):
        quadrilateral_surface(H, [(0, 1)])
    S = quadrilateral_surface(H, [(0, 3), (1, 2)])
    S.set_immutable()
    assert not S.is_mutable()


def test_self_glued_edge(H):
    with pytest.raises(ValueError, match="adjoin a square root of 4/3"):
        quadrilateral_surface(H, [(0, 0)])


def test_refused_glue_keeps_gluings(H):
    S = quadrilateral_surface(H, [(0, 3)])
    with pytest.raises(ValueError):
        S.glue((0, 0), (0, 1))
    assert S.opposite_edge(0, 0) == (0, 3)
    assert S.opposite_edge(0, 3) == (0, 0)
    assert S.opposite_edge(0, 1) is None


def test_glue_unknown_edge(H):
    S = quadrilateral_surface(H, [])
    with pytest.raises(ValueError, match="there is no polygon 1"):
        S.glue((0, 0), (1, 0))
    with pytest.raises(ValueError, match="polygon 0 has no edge 7"):
        S.glue((0, 7), (0, 0))


def test_replace_polygon_refused(H):
    S = quadrilateral_surface(H, [(0, 3), (1, 2)])
    P = S.polygon(0)
    wider = H.polygon(
        vertices=[(0, 0), (QQ(1) / 3, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]
    )
    with pytest.raises(ValueError, match="cannot replace polygon 0"):
        S.replace_polygon(0, wider)
    assert S.polygon(0) == P

    smaller = H.polygon(
        vertices=[(0, 0), (QQ(1) / 3, 0), (QQ(1) / 3, QQ(1) / 3), (0, QQ(1) / 3)]
    )
    S.replace_polygon(0, smaller)
    assert S.polygon(0) == smaller

    S.set_immutable()
    with pytest.raises(Exception, match="immutable"):
        S.replace_polygon(0, P)


def test_polygon_over_larger_ring(H):
    from sage.all import AA

    S = MutableOrientedHyperbolicSurface(QQ)
    P = HyperbolicPlane(AA).polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)])
    with pytest.raises(ValueError, match="does not coerce into Rational Field"):
        S.add_polygon(P)


def test_euclidean_polygon():
    S = MutableOrientedHyperbolicSurface(QQ)
    with pytest.raises(TypeError):
        S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (0, 1)]))


def test_point_coordinates(H, sphere):
    for p in [
        (QQ(1) / 4, QQ(1) / 4),
        (QQ(1) / 8, QQ(3) / 8),
        (QQ(1) / 4, 0),
        (QQ(1) / 2, QQ(1) / 4),
        (0, 0),
        (QQ(1) / 2, QQ(1) / 2),
    ]:
        point = sphere(0, p)
        assert H.point(p[0], p[1], model="klein") in point.coordinates(0)
        assert sphere(0, H.point(p[0], p[1], model="klein")) == point


@pytest.mark.parametrize("gluing", sorted(OCTAGON_GLUINGS))
def test_octagon(octagon, gluing):
    pairs, trace = OCTAGON_GLUINGS[gluing]
    K = octagon.base_ring()
    r = K.gen()

    assert octagon._cosh2_sides() == (in_radius(r, OCTAGON_SIDE_COSH2),) * 8

    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(octagon)
    for a, b in pairs:
        S.glue((0, a), (0, b))
    S.set_immutable()

    assert S.genus() == 2
    assert len(S.vertices()) == 1
    assert S.angles() == [1]
    assert S.area() == OCTAGON["area_over_2pi"]

    for edge in range(8):
        g = S.edge_transformation(0, edge)
        assert g.parent().base_ring() is K
        assert g.trace() == in_radius(r, trace)


def test_octagon_point(octagon):
    K = octagon.base_ring()
    r = K.gen()
    pairs, _ = OCTAGON_GLUINGS["opposite"]

    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(octagon)
    for a, b in pairs:
        S.glue((0, a), (0, b))
    S.set_immutable()

    a, b = (vector(octagon.vertex(i).coordinates(model="klein")) for i in (0, 1))
    m = (a + b) / 2
    p = S(0, m)
    image = S.edge_transformation(0, 0) * m
    assert octagon.get_point_position(image).get_edge() == 4
    assert set(p.coordinates(0)) == {
        S.hyperbolic_plane().point(m[0], m[1], model="klein"),
        S.hyperbolic_plane().point(image[0], image[1], model="klein"),
    }
    assert S(0, (r / 2, 0)).angle() == 1
