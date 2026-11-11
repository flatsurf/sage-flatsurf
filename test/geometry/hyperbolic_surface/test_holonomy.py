r"""
Holonomy representations of flat orbifolds and of hyperbolic cone surfaces,
checked against hard-coded values of their uniformisation groups.
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

from sage.all import QQ, RR, NumberField, matrix, polygen

from flatsurf import (
    HyperbolicPlane,
    MutableOrientedSimilaritySurface,
    Polygon,
    polygons,
    similarity_surfaces,
    translation_surfaces,
)
from flatsurf.geometry.holonomy import _holonomy_angle
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface

from hyperbolic_surface_oracles import (
    FLAT_ORBIFOLDS,
    OCTAGON_GLUINGS,
    TRIANGLE_DOUBLES,
    octagon_vertices,
    triangle_double_vertices,
)

CONSTRUCTORS = {
    "square_torus": lambda: translation_surfaces.square_torus(),
    "pillowcase": lambda: similarity_surfaces.pillowcase(polygons.square()),
    "equilateral_double": lambda: similarity_surfaces.pillowcase(
        Polygon(angles=(1, 1, 1))
    ),
    "isosceles_right_double": lambda: similarity_surfaces.pillowcase(
        Polygon(angles=(1, 1, 2))
    ),
    "triangle_236_double": lambda: similarity_surfaces.pillowcase(
        Polygon(angles=(1, 2, 3))
    ),
}

# Klein coordinates of a Lambert quadrilateral (three right angles).
QUADRILATERAL = [(0, 0), (QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]


@pytest.fixture(scope="module")
def H():
    return HyperbolicPlane(QQ)


@pytest.fixture(scope="module")
def octagon():
    K, points = octagon_vertices()
    return HyperbolicPlane(K).polygon(vertices=points)


def in_radius(r, coefficients):
    return sum(c * r**i for i, c in enumerate(coefficients))


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_relations(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    phi = S.holonomy_representation()
    loops = phi.domain().vertex_loops()
    assert sorted(v.angle() for v in loops) == angles

    # A cone point of angle 1/m is an orbifold point: the m-th power of its
    # loop maps to the identity.
    for vertex, loop in loops.items():
        a = vertex.angle()
        g = phi(loop)
        assert (g ** int(1 / a)).is_one()
        if a != 1:
            assert not g.is_one()

    phi._test_vertex_loops()


def test_flat_orbifold_images():
    G = CONSTRUCTORS["pillowcase"]().holonomy_representation().image()
    assert G.base_ring() is QQ
    assert G.ngens() == 3

    S = CONSTRUCTORS["triangle_236_double"]()
    G = S.holonomy_representation().image()
    assert G.base_ring() is S.base_ring()
    assert G.base_ring().degree() == 2
    assert G.base_ring().discriminant() == 12
    assert G.ngens() == 2


def test_square_rotations():
    # The square glued by rotations of angle pi/2: its vertices have angles
    # 1/4, 1/4, 1/2 and their loops map to rotations about the developed
    # vertices (1, 0), (0, 1), (1, 1).
    S = MutableOrientedSimilaritySurface(QQ)
    S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), label=0)
    S.glue((0, 0), (0, 1))
    S.glue((0, 2), (0, 3))
    S.set_immutable()

    phi = S.holonomy_representation()
    loops = phi.domain().vertex_loops()
    assert phi(loops[S(0, 1)]).matrix() == matrix(
        QQ, [[0, -1, 1], [1, 0, -1], [0, 0, 1]]
    )
    assert phi(loops[S(0, 3)]).matrix() == matrix(
        QQ, [[0, -1, 1], [1, 0, 1], [0, 0, 1]]
    )
    assert phi(loops[S(0, 0)]).matrix() == matrix(
        QQ, [[-1, 0, 2], [0, -1, 2], [0, 0, 1]]
    )

    with pytest.raises(ValueError):
        S.holonomy_representation(punctured_vertices=False)


def test_sphere_traces(H):
    S = MutableOrientedHyperbolicSurface(QQ)
    S.add_polygon(H.polygon(vertices=QUADRILATERAL))
    S.glue((0, 0), (0, 3))
    S.glue((0, 1), (0, 2))
    S.set_immutable()

    phi = S.holonomy_representation()
    loops = phi.domain().vertex_loops()
    # Vertices of angles 1/4, 1/2 and 0.19591...: traces 1 + 2 cos(2 pi a).
    assert [phi(loops[S(0, i)]).trace() for i in range(3)] == [
        1,
        -1,
        QQ(5) / 3,
    ]
    phi._test_vertex_loops()


@pytest.mark.parametrize("gluing", sorted(OCTAGON_GLUINGS))
def test_octagon(octagon, gluing):
    pairs, trace = OCTAGON_GLUINGS[gluing]
    K = octagon.base_ring()
    r = K.gen()

    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(octagon)
    for a, b in pairs:
        S.glue((0, a), (0, b))
    S.set_immutable()

    phi = S.holonomy_representation()
    F = phi.domain()
    loops = F.vertex_loops()
    assert len(loops) == 1
    for loop in loops.values():
        assert phi(loop).is_one()

    G = phi.image()
    assert G.ngens() == 4
    assert G.base_ring() is K

    # One polygon: the spanning tree is empty and each generator crosses a
    # single edge, so its image is a side pairing.
    for g in F.gens():
        assert len(g._walk) == 1
        assert phi(g) == S.edge_transformation(*g._walk[0])
        assert phi(g).trace() == in_radius(r, trace)

    phi._test_vertex_loops()

    # The vertex is regular, so it can be filled.
    psi = S.holonomy_representation(punctured_vertices=False)
    assert psi.domain() is S.fundamental_group()
    assert psi.image().gens() == G.gens()


@pytest.mark.parametrize("row", TRIANGLE_DOUBLES, ids=lambda row: str(row[0]))
def test_triangle_double_invariant_trace_field(row):
    # The invariant trace field depends on the periods only, not on the vertex
    # at the origin nor on whether a generator is a half-turn.
    K, front, mirror = triangle_double_vertices(row)
    H = HyperbolicPlane(K)
    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(H.polygon(vertices=front))
    S.add_polygon(H.polygon(vertices=mirror))
    for i in range(3):
        S.glue((0, i), (1, -i % 3))
    S.set_immutable()

    field = S.holonomy_representation().invariant_trace_field()
    expected = NumberField(polygen(QQ).parent()(row[3]), "b")
    assert field is not QQ
    assert field.is_isomorphic(expected)


def test_homomorphism():
    S = translation_surfaces.veech_double_n_gon(5)
    phi = S.holonomy_representation()
    gens = phi.domain().gens()
    words = [
        gens[0] * gens[1],
        gens[1] * ~gens[2] * gens[0],
        ~gens[3] * gens[2] * gens[2],
    ]
    for u in words:
        assert phi(~u) == ~phi(u)
        for v in words:
            assert phi(u * v) == phi(u) * phi(v)

    S = CONSTRUCTORS["triangle_236_double"]()
    phi = S.holonomy_representation()
    x, y = phi.domain().gens()
    for u, v in [(x, y), (x * y, ~x), (~y * x, y * y)]:
        assert phi(u * v) == phi(u) * phi(v)
        assert phi(~u) == ~phi(u)


def test_hyperbolic_angle_from_holonomy(H):
    T = MutableOrientedHyperbolicSurface(QQ)
    T.add_polygon(H.polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]))
    T.add_polygon(
        H.polygon(vertices=[(QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)])
    )
    T.glue((0, 1), (1, 2))
    T.glue((0, 0), (0, 2))
    T.glue((1, 0), (1, 1))
    T.set_immutable()

    # The corners at T(1, 0) have irrational angles summing to 1/2, the trace
    # of the loop is -1.
    assert T(1, 0).angle() == QQ(1) / 2
    assert T(0, 0).angle() == QQ(1) / 4
    with pytest.raises(ValueError):
        T(1, 1).angle()


def test_hyperbolic_angle_of_octagon_sector():
    x = polygen(QQ)
    K = NumberField(x**4 + 32 * x**2 - 32, "r", embedding=RR(0.98))
    r = K.gen()
    s2 = r**2 / 12 + QQ(4) / 3
    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(
        HyperbolicPlane(K).polygon(vertices=[(0, 0), (r, 0), (r * s2 / 2, r * s2 / 2)])
    )
    S.glue((0, 0), (0, 2))
    S.set_immutable()

    assert _holonomy_angle(S, S(0, 0)) == QQ(1) / 8
