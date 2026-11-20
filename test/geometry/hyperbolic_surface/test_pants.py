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

from sage.all import AA, QQ, RDF, RR, NumberField, polygen

from flatsurf import HyperbolicPlane, hyperbolic_polygons, hyperbolic_surfaces
from flatsurf.geometry.signature import SurfaceSignature

y = polygen(QQ, "y")

T = QQ(21) / 4


# The base ring is the field Q(u_0, u_1, u_2, p) of the slots, u_k = t_k + 1
# for a cuff of trace t_k and 4 cos^2 of the angle of a cone point.
@pytest.mark.parametrize(
    "boundary_traces,cone_orders,polynomial,area",
    [
        ((T, T, T), (), None, 1),
        ((T, T), (-QQ(2) / 3,), None, QQ(2) / 3),
        ((T, T), (-QQ(1) / 2,), None, QQ(1) / 2),
        ((T,), (-QQ(2) / 3, -QQ(2) / 3), None, QQ(1) / 3),
        ((T,), (-QQ(2) / 3, -QQ(1) / 2), None, QQ(1) / 6),
        ((T,), (-QQ(1) / 2, -QQ(2) / 3), None, QQ(1) / 6),
        ((T, T, QQ(20881) / 144), (), None, 1),
        ((4, 4, 4), (), y**2 - y - 1, 1),
        ((5, 5, 5), (), y**2 - 6, 1),
        ((4, 5, T), (), y**2 - 30, 1),
        ((T, 4, 5), (), y**2 - 30, 1),
        ((8, 8, 8), (), None, 1),
        ((15, 15, 15), (), None, 1),
        ((4, 5), (-QQ(2) / 3,), y**2 - 30, QQ(2) / 3),
        ((4,), (-QQ(2) / 3, -QQ(4) / 5), y**2 - y - 1, QQ(7) / 15),
        ((T, T), (-QQ(1) / 4,), y**2 - 2, QQ(1) / 4),
        ((), (-QQ(1) / 2, -QQ(3) / 4, -QQ(7) / 8), y**2 - 2, QQ(1) / 8),
    ],
)
def test_signature_round_trip(boundary_traces, cone_orders, polynomial, area):
    S = hyperbolic_surfaces.pants(
        boundary_traces=boundary_traces, cone_orders=cone_orders
    )
    assert not S.is_mutable()
    assert list(S.labels()) == [0, 1]
    if polynomial is None:
        assert S.base_ring() is QQ
    else:
        assert S.base_ring().is_isomorphic(NumberField(polynomial, "b"))
    assert S.signature() == SurfaceSignature(
        0, cone_orders=cone_orders, boundary_traces=boundary_traces
    )
    assert S.area() == area
    field = S.holonomy_representation().invariant_trace_field()
    assert field is S.base_ring() or field.is_isomorphic(S.base_ring())


# The plane of x^2 + y^2 - z^2: over the base ring when c or c/b is a square
# (21/4 three times, two cuffs with a cone point of angle pi, and 4, 5, 21/4,
# where c/b is a square) or when the line of a half-cuff has its ends over it
# (the cuff of trace 21/4 in the first slot), otherwise over a quadratic
# extension of degree twice its degree.
@pytest.mark.parametrize(
    "boundary_traces,cone_orders,degree",
    [
        ((T, T, T), (), 1),
        ((T, T), (-QQ(1) / 2,), 1),
        ((4, 5, T), (), 2),
        ((T, 4, 5), (), 2),
        ((4,), (-QQ(2) / 3, -QQ(4) / 5), 2),
        ((5, 5, 5), (), 4),
        ((4, 5), (-QQ(2) / 3,), 4),
        ((15, 15, 15), (), 2),
    ],
)
def test_standard(boundary_traces, cone_orders, degree):
    S = hyperbolic_surfaces.pants(
        boundary_traces=boundary_traces, cone_orders=cone_orders, standard=True
    )
    assert S.hyperbolic_plane() is HyperbolicPlane(S.base_ring())
    assert (1 if S.base_ring() is QQ else S.base_ring().degree()) == degree
    assert S.signature() == SurfaceSignature(
        0, cone_orders=cone_orders, boundary_traces=boundary_traces
    )


def test_triangle():
    # Three cone points: the pillowcase of the triangle of angles 1/16, 1/4,
    # 1/8 at its vertices 0, 1, 2, in the plane of its form.
    cone_orders = (-QQ(1) / 2, -QQ(3) / 4, -QQ(7) / 8)
    T = hyperbolic_polygons.triangle(QQ(1) / 16, QQ(1) / 4, QQ(1) / 8)
    S = hyperbolic_surfaces.pants(cone_orders=cone_orders)
    assert S.hyperbolic_plane() is T.parent()
    assert S.polygon(0) == T
    S = hyperbolic_surfaces.pants(cone_orders=cone_orders, base_ring=AA)
    assert S.hyperbolic_plane() is T.parent().change_ring(AA)
    assert S.signature() == SurfaceSignature(0, cone_orders=cone_orders)
    with pytest.raises(ValueError, match="do not lie in Rational Field"):
        hyperbolic_surfaces.pants(cone_orders=cone_orders, base_ring=QQ)


def klein_vertices(S, label):
    return [v.coordinates(model="klein") for v in S.polygon(label).vertices()]


def test_seam_plane_hexagon():
    S = hyperbolic_surfaces.pants(boundary_traces=(T, T, T))
    assert S.hyperbolic_plane() is HyperbolicPlane(QQ, form=[1, 6, -6])
    assert klein_vertices(S, 0) == [
        (0, 0),
        (QQ(12) / 5, 0),
        (QQ(12) / 5, QQ(3) / 25),
        (QQ(60) / 29, QQ(15) / 29),
        (QQ(48) / 25, QQ(3) / 5),
        (0, QQ(3) / 5),
    ]
    assert S.polygon(0).angles() == (QQ(1) / 4,) * 6
    # cosh^2 of the half-cuffs: (t + 1) / 4 = 25/16
    assert S.polygon(0)._cosh2_sides()[1::2] == (QQ(25) / 16,) * 3
    assert S.boundary_components() == [
        ((0, 1), (1, 5)),
        ((0, 3), (1, 3)),
        ((0, 5), (1, 1)),
    ]
    assert S.opposite_edge(0, 0) == (1, 0)
    assert S.edge_transformation(0, 0).is_one()


def test_seam_plane_quadrilateral():
    S = hyperbolic_surfaces.pants(
        boundary_traces=(T,), cone_orders=(-QQ(1) / 2, -QQ(2) / 3)
    )
    assert S.hyperbolic_plane() is HyperbolicPlane(QQ, form=[1, 13, -13])
    assert klein_vertices(S, 0) == [(0, 0), (2, 0), (2, QQ(3) / 5), (0, QQ(3) / 5)]
    assert S.polygon(0).angles() == (QQ(1) / 4, QQ(1) / 4, QQ(1) / 6, QQ(1) / 4)


def test_ring():
    assert (
        hyperbolic_surfaces.pants(boundary_traces=(T, T, T), base_ring=QQ).base_ring()
        is QQ
    )
    S = hyperbolic_surfaces.pants(boundary_traces=(T, T, T), base_ring=AA)
    assert S.hyperbolic_plane() is HyperbolicPlane(AA, form=[1, 6, -6])
    S = hyperbolic_surfaces.pants(boundary_traces=(5, 5, 5), base_ring=AA, standard=True)
    assert S.hyperbolic_plane() is HyperbolicPlane(AA)
    with pytest.raises(ValueError, match="do not lie in Rational Field"):
        hyperbolic_surfaces.pants(boundary_traces=(5, 5, 5), base_ring=QQ)


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
        hyperbolic_surfaces.regular_octagon(base_ring=QQ)
    with pytest.raises(ValueError, match="the ring must be exact, not Real Field"):
        hyperbolic_surfaces.regular_octagon(base_ring=RR)


def test_pants_inexact_ring_refused():
    with pytest.raises(
        ValueError,
        match="the ring must be exact, not Real Field with 53 bits of precision",
    ):
        hyperbolic_surfaces.pants(boundary_traces=(T, T, T), base_ring=RR)
    with pytest.raises(ValueError, match="the ring must be exact, not Real Double"):
        hyperbolic_surfaces.pants(boundary_traces=(T, T, T), base_ring=RDF)
