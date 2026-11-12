r"""
Checks of the boundary components, the boundary loops and the signature of
flat and hyperbolic cone surfaces against hard-coded values.
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

from sage.all import AA, QQ, polygen

from flatsurf import (
    HyperbolicPlane,
    MutableOrientedSimilaritySurface,
    Polygon,
    similarity_surfaces,
    translation_surfaces,
)
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
from flatsurf.geometry.signature import SurfaceSignature

from hyperbolic_surface_oracles import (
    COMPACT_SIGNATURES,
    CUBE_DIVISOR,
    CUSPED_SIGNATURES,
    FLAT_SIGNATURES,
    HEXAGON_SURFACES,
    OCTAGON_GLUINGS,
    hexagon_vertices,
    octagon_vertices,
)

# Klein coordinates of a Lambert quadrilateral (three right angles).
QUADRILATERAL = [(0, 0), (QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]

SQUARE = [(0, 0), (1, 0), (1, 1), (0, 1)]


def flat_surface(vertices, gluings):
    S = MutableOrientedSimilaritySurface(QQ)
    S.add_polygon(Polygon(vertices=vertices))
    for a, b in gluings:
        S.glue((0, a), (0, b))
    S.set_immutable()
    return S


def hexagon_surface(name):
    gluings, _, _ = HEXAGON_SURFACES[name]
    P = HyperbolicPlane(AA).polygon(vertices=hexagon_vertices())
    S = MutableOrientedHyperbolicSurface(AA)
    S.add_polygon(P)
    S.add_polygon(P)
    for a, b in gluings:
        S.glue(a, b)
    S.set_immutable(minimal_field=True)
    return S


def octagon_surface(gluing):
    K, vertices = octagon_vertices()
    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(HyperbolicPlane(K).polygon(vertices=vertices))
    for a, b in OCTAGON_GLUINGS[gluing][0]:
        S.glue((0, a), (0, b))
    S.set_immutable()
    return S


def test_flat_cylinder_boundary():
    S = flat_surface(SQUARE, [(0, 2)])
    assert S.euler_characteristic() == 0
    assert sorted(S.angles()) == [QQ(1) / 2, QQ(1) / 2]
    assert S.boundary_components() == [((0, 1),), ((0, 3),)]

    F = S.fundamental_group(punctured_vertices=True)
    (g,) = F.gens()
    a, b = F.boundary_loops()
    assert a in (g, ~g)
    assert b in (g, ~g)
    assert a == ~b


def test_flat_disk_boundary():
    S = flat_surface([(0, 0), (1, 0), (1, 1)], [])
    assert S.boundary_components() == [((0, 0), (0, 1), (0, 2))]
    for punctured in (False, True):
        (loop,) = S.fundamental_group(punctured_vertices=punctured).boundary_loops()
        assert loop.is_one()


def test_flat_torus_boundary():
    S = translation_surfaces.square_torus()
    assert S.boundary_components() == []
    assert S.fundamental_group().boundary_loops() == []


def test_flat_cylinder_signature():
    S = flat_surface(SQUARE, [(0, 2)])
    s = S.signature()
    assert s == SurfaceSignature(0, boundary_traces=(None, None))
    assert s.genus() == 0
    assert s.cone_orders() == ()
    assert s.boundary_traces() == (None, None)
    assert s.euler_characteristic() == S.euler_characteristic() == 0
    assert s.divisor_euler_characteristic() == 0


def test_flat_corner():
    S = flat_surface([(0, 0), (1, 0), (1, 1)], [])
    assert sorted(S.angles()) == [QQ(1) / 8, QQ(1) / 8, QQ(1) / 4]
    with pytest.raises(ValueError, match="not geodesic"):
        S.signature()


def test_regular_octagon_signature():
    # Stratum H_2(2): a single zero of order 2.
    S = translation_surfaces.regular_octagon()
    assert S.euler_characteristic() == -2
    assert S.angles() == [3]
    s = S.signature()
    assert s.genus() == 2
    assert s.cone_orders() == (2,)
    assert s.euler_characteristic() == S.euler_characteristic() == -2
    assert s.divisor_euler_characteristic() == 0
    assert s.geometry() == "euclidean"


def test_self_glued_polygon_signature():
    S = similarity_surfaces.self_glued_polygon(
        Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
    )
    assert S.euler_characteristic() == 2
    assert S.angles() == [1]
    s = S.signature()
    assert s.genus() == 0
    assert s.cone_orders() == (-QQ(1) / 2,) * 4
    assert s.euler_characteristic() == S.euler_characteristic() == 2
    assert s.divisor_euler_characteristic() == 0
    assert s.orbifold_orders() == (0, (2, 2, 2, 2), 0)


@pytest.mark.parametrize("name", sorted(HEXAGON_SURFACES))
def test_hexagon_surface(name):
    gluings, genus, traces = HEXAGON_SURFACES[name]
    S = hexagon_surface(name)

    K = S.base_ring()
    assert K.polynomial() == polygen(QQ, "y") ** 4 - 4 * polygen(QQ, "y") ** 2 + 1
    assert S.euler_characteristic() == -1
    assert S.area() == 1
    assert S.angles() == [QQ(1) / 2] * 6
    assert S.is_with_boundary()
    assert len(S.boundary_components()) == len(traces)

    s = S.signature()
    assert s == SurfaceSignature(genus, boundary_traces=traces)
    assert s.genus() == genus
    assert s.cone_orders() == ()
    assert s.cusps() == 0
    assert s.boundary_components() == len(traces)
    assert s.boundary_traces() == traces
    assert s.euler_characteristic() == S.euler_characteristic() == -1
    assert s.divisor_euler_characteristic() == -1
    assert s.area() == S.area()
    with pytest.raises(ValueError):
        s.orbifold_orders()

    S._test_gauss_bonnet()


def test_pants_holonomy():
    S = hexagon_surface("pants")
    phi = S.holonomy_representation()
    F = phi.domain()
    x = polygen(QQ, "x")
    for walk in ([(1, 0), (0, 2)], [(1, 0), (0, 4)]):
        M = phi(F(walk)).matrix()
        assert M.trace() == 15
        assert M.charpoly() == x**3 - 15 * x**2 + 15 * x - 1


def test_torus_holonomy():
    S = hexagon_surface("torus")
    phi = S.holonomy_representation()
    F = phi.domain()
    a = F([(1, 0), (0, 2)])
    b = F([(1, 0), (0, 4)])
    assert phi(a).matrix().trace() == 24
    assert phi(b).matrix().trace() == 24
    assert phi(a * b * ~a * ~b).matrix().trace() == 2703


@pytest.mark.parametrize("gluing", sorted(OCTAGON_GLUINGS))
def test_octagon_signature(gluing):
    S = octagon_surface(gluing)
    s = S.signature()
    assert s.genus() == 2
    assert s.cone_orders() == ()
    assert s.area() == 2
    assert S.area() == 2
    assert s.euler_characteristic() == S.euler_characteristic() == -2
    assert s.divisor_euler_characteristic() == -2
    assert s.orbifold_orders() == (2, (), 0)


def test_288_signature():
    K, vertices = octagon_vertices()
    P = HyperbolicPlane(K).polygon(vertices=vertices).cyclic_quotient(8)
    S = MutableOrientedHyperbolicSurface(K)
    S.add_polygon(P)
    S.glue((0, 0), (0, 2))
    S.glue((0, 1), (0, 1))
    S.set_immutable()

    s = S.signature()
    assert s.cone_orders() == (-QQ(7) / 8, -QQ(7) / 8, -QQ(1) / 2)
    assert s.euler_characteristic() == S.euler_characteristic() == 2
    assert s.divisor_euler_characteristic() == -QQ(1) / 4
    assert s.orbifold_orders() == (0, (2, 8, 8), 0)
    assert S.area() == s.area()


def test_lambert_sphere_signature():
    S = MutableOrientedHyperbolicSurface(QQ)
    S.add_polygon(HyperbolicPlane(QQ).polygon(vertices=QUADRILATERAL))
    S.glue((0, 0), (0, 3))
    S.glue((0, 1), (0, 2))
    S.set_immutable()
    with pytest.raises(NotImplementedError):
        S.signature()
    S._test_gauss_bonnet()


@pytest.mark.parametrize("g,periods,area", COMPACT_SIGNATURES)
def test_compact_signatures(g, periods, area):
    s = SurfaceSignature.from_orbifold_orders(g, periods)
    assert s.area() == area / 2
    assert s.geometry() == "hyperbolic"
    assert s.euler_characteristic() == 2 - 2 * g
    assert s.orbifold_orders() == (g, tuple(sorted(periods)), 0)


@pytest.mark.parametrize("g,periods,cusps,area", CUSPED_SIGNATURES)
def test_cusped_signatures(g, periods, cusps, area):
    s = SurfaceSignature.from_orbifold_orders(g, periods, cusps=cusps)
    assert s.cusps() == cusps
    assert s.area() == area / 2
    assert s.euler_characteristic() == 2 - 2 * g - cusps
    assert s.orbifold_orders() == (g, tuple(sorted(periods)), cusps)


@pytest.mark.parametrize("g,periods", FLAT_SIGNATURES)
def test_flat_signatures(g, periods):
    s = SurfaceSignature.from_orbifold_orders(g, periods)
    assert s.divisor_euler_characteristic() == 0
    assert s.geometry() == "euclidean"
    assert s.orbifold_orders() == (g, tuple(sorted(periods)), 0)


def test_cube_divisor():
    s = SurfaceSignature(CUBE_DIVISOR["genus"], CUBE_DIVISOR["orders"])
    assert (
        s.divisor_euler_characteristic()
        == CUBE_DIVISOR["divisor_euler_characteristic"]
    )
    with pytest.raises(ValueError):
        s.orbifold_orders()


def test_signature_constructor():
    with pytest.raises(ValueError):
        SurfaceSignature(-1)
    with pytest.raises(ValueError):
        SurfaceSignature(0, cusps=-1)
    with pytest.raises(ValueError):
        SurfaceSignature(0, [-1])
    with pytest.raises(ValueError):
        SurfaceSignature(0, [-QQ(3) / 2])
    with pytest.raises(ValueError):
        SurfaceSignature(0, [0])
    with pytest.raises(ValueError, match="all None or none of them None"):
        SurfaceSignature(0, boundary_traces=(None, 15))
    with pytest.raises(ValueError, match="all None or none of them None"):
        SurfaceSignature(0, boundary_traces=(15, None))
    assert SurfaceSignature(0, boundary_traces=(None, None)).boundary_traces() == (
        None,
        None,
    )

    s = SurfaceSignature(0, [-QQ(1) / 2, -QQ(7) / 8, -QQ(7) / 8])
    assert repr(s) == "SurfaceSignature(g=0, cone_orders=(-7/8, -7/8, -1/2))"
    assert hash(s) == hash(SurfaceSignature.from_orbifold_orders(0, (8, 2, 8)))
    assert s != SurfaceSignature(0)

    assert repr(SurfaceSignature(2)) == "SurfaceSignature(g=2)"
    t = SurfaceSignature(1, cusps=2, boundary_traces=(15, 3))
    assert repr(t) == "SurfaceSignature(g=1, cusps=2, boundary_traces=(3, 15))"
    assert eval(repr(t)) == t  # noqa: S307
    assert t.euler_characteristic() == -4
    assert t.divisor_euler_characteristic() == -4

    with pytest.raises(ValueError, match="no orbifold orders"):
        t.orbifold_orders()
    with pytest.raises(ValueError, match="period 1 is not"):
        SurfaceSignature.from_orbifold_orders(0, (1, 2))
