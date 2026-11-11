r"""
Pairs of pants and the regular octagon of genus 2: the signature round trip,
the base ring and the refused arguments, against hard-coded values.

A pair of pants with cuffs of traces `t_i` and cone points of orders `k_j`
has the signature of genus 0 with these traces and orders, and the area
`-\chi = len(t) - 2 - \sum k_j`. The base rings and the Klein coordinates
below were computed with sage-flatsurf; the octagon gluings are those of
``OCTAGON_GLUINGS`` in ``hyperbolic_surface_oracles``.
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

from flatsurf import hyperbolic_surfaces
from flatsurf.geometry.signature import SurfaceSignature

y = polygen(QQ, "y")

T = QQ(21) / 4


@pytest.mark.parametrize(
    "boundary_traces,cone_orders,optimize_placement,polynomial,area",
    [
        ((T, T, T), (), True, None, 1),
        ((T, T), (-QQ(2) / 3,), True, None, QQ(2) / 3),
        ((T, T), (-QQ(1) / 2,), True, None, QQ(1) / 2),
        ((T,), (-QQ(2) / 3, -QQ(2) / 3), True, None, QQ(1) / 3),
        ((T,), (-QQ(2) / 3, -QQ(1) / 2), True, None, QQ(1) / 6),
        ((T,), (-QQ(1) / 2, -QQ(2) / 3), True, None, QQ(1) / 6),
        ((T, T, T), (), False, y**2 - 6, 1),
        ((T, T), (-QQ(2) / 3,), False, y**2 - 7, QQ(2) / 3),
        ((T, T), (-QQ(1) / 2,), False, y**2 - 34, QQ(1) / 2),
        ((T,), (-QQ(2) / 3, -QQ(2) / 3), False, y**2 - 3, QQ(1) / 3),
        ((T,), (-QQ(2) / 3, -QQ(1) / 2), False, y**2 - y - 3, QQ(1) / 6),
        ((T, T, QQ(20881) / 144), (), True, None, 1),
        ((T, T, QQ(20881) / 144), (), False, None, 1),
        ((4, 4, 4), (), True, y**4 - 2 * y**2 - 4, 1),
        ((T, T), (-QQ(1) / 4,), True, y**4 - 56 * y**2 - 72 * y + 414, QQ(1) / 4),
        ((), (-QQ(1) / 2, -QQ(3) / 4, -QQ(7) / 8), True, y**4 - 2, QQ(1) / 8),
    ],
)
def test_signature_round_trip(
    boundary_traces, cone_orders, optimize_placement, polynomial, area
):
    S = hyperbolic_surfaces.pants(
        boundary_traces=boundary_traces,
        cone_orders=cone_orders,
        optimize_placement=optimize_placement,
    )
    assert not S.is_mutable()
    assert list(S.labels()) == [0, 1]
    if polynomial is None:
        assert S.base_ring() is QQ
    else:
        assert S.base_ring().polynomial() == polynomial
    assert S.signature() == SurfaceSignature(
        0, cone_orders=cone_orders, boundary_traces=boundary_traces
    )
    assert S.area() == area


def klein_vertices(S, label):
    return [v.coordinates(model="klein") for v in S.polygon(label).vertices()]


def test_moved_hexagon():
    S = hyperbolic_surfaces.pants(boundary_traces=(T, T, T))
    assert klein_vertices(S, 0) == [
        (-QQ(1) / 5, 0),
        (-QQ(1) / 5, -QQ(24) / 25),
        (-QQ(5) / 61, -QQ(60) / 61),
        (QQ(23) / 65, -QQ(12) / 13),
        (QQ(5) / 11, -QQ(48) / 55),
        (QQ(5) / 11, 0),
    ]
    assert S.boundary_components() == [
        ((0, 1), (1, 5)),
        ((0, 3), (1, 3)),
        ((0, 5), (1, 1)),
    ]
    assert S.opposite_edge(0, 0) == (1, 0)
    assert S.edge_transformation(0, 0).is_one()


def test_moved_quadrilateral():
    S = hyperbolic_surfaces.pants(
        boundary_traces=(T,), cone_orders=(-QQ(1) / 2, -QQ(2) / 3)
    )
    assert klein_vertices(S, 0)[0] == (-QQ(6) / 7, 0)
    assert klein_vertices(S, 0)[2:] == [(-QQ(9) / 17, -QQ(10) / 17), (-QQ(9) / 17, 0)]
    assert S.polygon(0).angles() == (QQ(1) / 4, QQ(1) / 4, QQ(1) / 6, QQ(1) / 4)


def test_ring():
    assert hyperbolic_surfaces.pants(boundary_traces=(T, T, T), ring=QQ).base_ring() is QQ
    assert hyperbolic_surfaces.pants(boundary_traces=(T, T, T), ring=AA).base_ring() is AA
    with pytest.raises(ValueError, match="do not lie in Rational Field"):
        hyperbolic_surfaces.pants(
            boundary_traces=(T, T, T), ring=QQ, optimize_placement=False
        )


@pytest.mark.parametrize(
    "boundary_traces,cone_orders,message",
    [
        ((T, T), (), "three cuffs and cone points in total, not 2"),
        ((T, T, T, T), (), "three cuffs and cone points in total, not 4"),
        ((T, T, 3), (), "the boundary trace 3 must be larger than 3"),
        ((T, T), (-1,), "the cone order -1 is not larger than -1"),
        ((T, T), (0,), "a cone order must be nonzero"),
        ((T,), (-QQ(1) / 2, -QQ(1) / 2), "Euler characteristic .* not 0"),
        ((), (-QQ(1) / 2, -QQ(2) / 3, -QQ(5) / 6), "Euler characteristic .* not 0"),
        (
            (QQ(20881) / 144, T),
            (-QQ(1) / 4,),
            r"c2 \+ c1 cos\(π\(1 \+ k0\)\) must be positive",
        ),
        (
            (T, QQ(20881) / 144),
            (-QQ(1) / 4,),
            r"c1 \+ c2 cos\(π\(1 \+ k0\)\) must be positive",
        ),
        (
            (QQ(20881) / 144,),
            (-QQ(1) / 4, -QQ(7) / 8),
            r"cos\(π\(1 \+ k1\)\) \+ c2 cos\(π\(1 \+ k0\)\) must be positive",
        ),
        (
            (QQ(20881) / 144,),
            (-QQ(7) / 8, -QQ(1) / 4),
            r"cos\(π\(1 \+ k0\)\) \+ c2 cos\(π\(1 \+ k1\)\) must be positive",
        ),
    ],
)
def test_refused_arguments(boundary_traces, cone_orders, message):
    with pytest.raises(ValueError, match=message):
        hyperbolic_surfaces.pants(
            boundary_traces=boundary_traces, cone_orders=cone_orders
        )


def test_positive_cone_order():
    with pytest.raises(NotImplementedError, match="larger than 2π"):
        hyperbolic_surfaces.pants(boundary_traces=(T, T), cone_orders=(QQ(1) / 2,))


@pytest.mark.parametrize(
    "gluing,pairs",
    [
        ("opposite", [(0, 4), (1, 5), (2, 6), (3, 7)]),
        ("katok", [(0, 2), (1, 3), (4, 6), (5, 7)]),
    ],
)
def test_regular_octagon(gluing, pairs):
    S = hyperbolic_surfaces.regular_octagon(gluing=gluing)
    assert not S.is_mutable()
    assert S.base_ring().polynomial() == y**4 - 2 * y**2 - 1
    for e, f in pairs:
        assert S.opposite_edge(0, e) == (0, f)
        assert S.opposite_edge(0, f) == (0, e)
    assert S.signature() == SurfaceSignature(2)
    assert S.area() == 2


def test_regular_octagon_refused():
    with pytest.raises(ValueError, match='"opposite" or "katok"'):
        hyperbolic_surfaces.regular_octagon(gluing="adjacent")
    with pytest.raises(ValueError, match="do not lie in Rational Field"):
        hyperbolic_surfaces.regular_octagon(ring=QQ)
