r"""
Surfaces from a signature: the signature round trip, the area, the base
ring, the number of polygons and the refused signatures, against hard-coded
values.

A surface of signature `(g; k_1, \dots, k_c; b)` with `n` legs is glued from
`2g - 2 + n` pairs of pants, two polygons each, a single one for three cone
points; its area is `-\chi(S, \beta)`. The base rings and invariant trace
fields below were computed with sage-flatsurf. The orbifolds are those of
``COMPACT_SIGNATURES`` in ``hyperbolic_surface_oracles``, whose areas are in
units of `\pi`, here halved to units of `2\pi`.
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

from sage.all import AA, QQ, RR, NumberField, TestSuite, polygen

from flatsurf import HyperbolicPlane, hyperbolic_polygons, hyperbolic_surfaces
from flatsurf.geometry.signature import SurfaceSignature

from hyperbolic_surface_oracles import COMPACT_SIGNATURES

y = polygen(QQ, "y")

T = QQ(21) / 4


def orbifold(g, periods):
    return SurfaceSignature.from_orbifold_orders(g, periods)


def check(S, signature, area, polynomial, polygons):
    assert not S.is_mutable()
    assert S.signature() == signature
    assert S.area() == area
    if polynomial is None:
        assert S.base_ring() is QQ
    else:
        assert S.base_ring().is_isomorphic(NumberField(polynomial, "b"))
    assert len(S.labels()) == polygons


# The polynomials of the base fields are regression values of the current
# construction: they depend on the placement of each pair of pants and on the
# order of the pairs of pants and of the slots in the chain, not only on the
# signature. A triangle, from three cone points, is over the field k_0 of
# its invariant form.
@pytest.mark.parametrize(
    "signature,area,polynomial,polygons",
    [
        (SurfaceSignature(2), 2, None, 4),
        (SurfaceSignature(3), 4, None, 8),
        (orbifold(1, (2,)), QQ(1) / 2, None, 2),
        (orbifold(0, (2, 2, 2, 2, 2)), QQ(1) / 2, None, 2),
        (orbifold(0, (2, 2, 2, 3)), QQ(1) / 6, None, 2),
        (orbifold(0, (3, 3, 3, 3)), QQ(2) / 3, None, 4),
        (orbifold(2, (2,)), QQ(5) / 2, y**2 - 51, 6),
        (SurfaceSignature(1, boundary_traces=(T,)), 1, None, 2),
        (
            SurfaceSignature(0, cone_orders=(-QQ(1) / 2,) * 2, boundary_traces=(T, T)),
            1,
            None,
            2,
        ),
        (orbifold(0, (3, 3, 4)), QQ(1) / 12, y**2 - 2, 2),
        (orbifold(0, (2, 4, 5)), QQ(1) / 20, y**2 - y - 1, 2),
        (orbifold(0, (3, 4, 4)), QQ(1) / 6, None, 2),
        (orbifold(0, (2, 3, 7)), QQ(1) / 42, y**3 - 7 * y - 7, 2),
        (
            SurfaceSignature(0, cone_orders=(-QQ(1) / 3, -QQ(5) / 6, -QQ(7) / 8)),
            QQ(1) / 24,
            y**4 - 16 * y**2 - 24 * y - 2,
            2,
        ),
    ],
)
def test_signature_round_trip(signature, area, polynomial, polygons):
    S = hyperbolic_surfaces.from_signature(signature)
    check(S, signature, area, polynomial, polygons)


def test_triangle_in_form_plane():
    # The pillowcase of the triangle with the cone points of periods 3, 4, 3
    # at its vertices 0, 1, 2, in the plane of its form.
    S = hyperbolic_surfaces.from_signature(orbifold(0, (3, 3, 4)))
    T = hyperbolic_polygons.triangle(QQ(1) / 6, QQ(1) / 8, QQ(1) / 6)
    assert S.hyperbolic_plane() is T.parent()
    assert S.polygon((0, 0)).angles() == T.angles()

    S = hyperbolic_surfaces.from_signature(orbifold(0, (3, 3, 4)), base_ring=AA)
    assert S.base_ring() is AA
    assert S.hyperbolic_plane() is T.parent().change_ring(AA)
    assert S.signature() == orbifold(0, (3, 3, 4))


def test_labels():
    S = hyperbolic_surfaces.from_signature(orbifold(2, (2,)))
    assert list(S.labels()) == [(i, side) for i in range(3) for side in (0, 1)]


def test_default_traces_genus_2_with_cone_point():
    S = hyperbolic_surfaces.from_signature(orbifold(2, (2,)))
    field = S.holonomy_representation().invariant_trace_field()
    assert field.polynomial() == y**2 - 51


def test_cuff_traces_99():
    signature = orbifold(2, (2,))
    S = hyperbolic_surfaces.from_signature(signature, cuff_traces=(99, 99, 99, 99))
    check(S, signature, QQ(5) / 2, None, 6)
    assert S.holonomy_representation().invariant_trace_field() is QQ


@pytest.mark.parametrize(
    "g,periods,area",
    COMPACT_SIGNATURES,
)
def test_compact_signatures(g, periods, area):
    signature = orbifold(g, periods)
    S = hyperbolic_surfaces.from_signature(signature)
    assert S.signature() == signature
    assert S.area() == area / 2


# Genus 2 from the pants (l, l, c) and (c, l', l'), the second developed
# across the cuff c, the edge 5 of the polygon (0, 0) and the edge 1 of the
# polygon (1, 0), into the seam plane of the first: the two front polygons
# share the two ends of that edge, glued by the identity, and the base ring
# is the invariant trace field.
@pytest.mark.parametrize(
    "cuff_traces,polynomial",
    [
        ((T, T, T), None),
        ((4, 4, 4), y**2 - y - 1),
        ((4, 5, T), y**4 - 552 * y**2 - 1200 * y + 54696),
    ],
)
def test_developed_genus_2(cuff_traces, polynomial):
    S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=cuff_traces)
    check(S, SurfaceSignature(2), 2, polynomial, 4)
    front, developed = S.polygon((0, 0)), S.polygon((1, 0))
    assert len([v for v in developed.vertices() if v in front.vertices()]) == 2
    assert S.edge_transformation((0, 0), 5).matrix().is_one()
    field = S.holonomy_representation().invariant_trace_field()
    assert S.base_ring().degree() == field.degree()


# Genus 3 from the pants (l, l, c0), (c0, d, e), (d, e, c1), (c1, l', l'):
# the side pairings of the loops and of the second cuff e of the handle join
# two polygons already in the seam plane of the first pants.
@pytest.mark.parametrize(
    "cuff_traces,polynomial",
    [
        (
            (T, T, 4, 5, T, T),
            y**8
            - 42 * y**6
            - 36 * y**5
            - 331 * y**4
            - 2184 * y**3
            + 1836 * y**2
            + 5076 * y
            + 68881,
        ),
        ((4, T, T, T, T, 5), y**4 - 2 * y**3 - 665 * y**2 + 666 * y + 106599),
    ],
)
def test_genus_3_cuff_traces(cuff_traces, polynomial):
    S = hyperbolic_surfaces.from_signature(SurfaceSignature(3), cuff_traces=cuff_traces)
    check(S, SurfaceSignature(3), 4, polynomial, 8)
    assert S.base_ring().polynomial() == polynomial
    assert S.edge_transformation((0, 0), 5).matrix().is_one()
    TestSuite(S).run()


def test_standard():
    S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), standard=True)
    assert S.hyperbolic_plane() is HyperbolicPlane(QQ)
    assert S.signature() == SurfaceSignature(2)
    S = hyperbolic_surfaces.from_signature(
        SurfaceSignature(2), base_ring=AA, standard=True
    )
    assert S.hyperbolic_plane() is HyperbolicPlane(AA)


def test_ring():
    S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), base_ring=AA)
    assert S.base_ring() is AA
    assert S.signature() == SurfaceSignature(2)
    assert (
        hyperbolic_surfaces.from_signature(
            SurfaceSignature(2), base_ring=QQ
        ).base_ring()
        is QQ
    )
    with pytest.raises(ValueError, match="do not lie in Rational Field"):
        hyperbolic_surfaces.from_signature(orbifold(2, (2,)), base_ring=QQ)
    with pytest.raises(ValueError, match="the ring must be exact, not Real Field"):
        hyperbolic_surfaces.from_signature(SurfaceSignature(2), base_ring=RR)


@pytest.mark.parametrize(
    "signature,error,message",
    [
        (
            SurfaceSignature.from_orbifold_orders(0, (2, 3), cusps=1),
            NotImplementedError,
            "cusps are not implemented",
        ),
        (
            SurfaceSignature(2, cone_orders=(QQ(1) / 2,)),
            NotImplementedError,
            "the cone order 1/2 is positive",
        ),
        (
            SurfaceSignature(1, boundary_traces=(3,)),
            ValueError,
            "the boundary trace 3 must be larger than 3",
        ),
        (
            SurfaceSignature(1, boundary_traces=(None,)),
            ValueError,
            "must be given, not None",
        ),
        (SurfaceSignature(1), ValueError, "is euclidean, not hyperbolic"),
        (orbifold(0, (2, 3, 5)), ValueError, "is spherical, not hyperbolic"),
        (orbifold(0, (2, 2, 2, 2)), ValueError, "is euclidean, not hyperbolic"),
        (
            SurfaceSignature(
                0, cone_orders=(-QQ(1) / 4,), boundary_traces=(T, QQ(20881) / 144)
            ),
            NotImplementedError,
            "may exist without such a decomposition",
        ),
    ],
)
def test_refused_signatures(signature, error, message):
    with pytest.raises(error, match=message):
        hyperbolic_surfaces.from_signature(signature)


def test_refused_cuff_traces():
    with pytest.raises(ValueError, match="must have 3 traces, one per interior cuff, not 1"):
        hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(T,))
    with pytest.raises(ValueError, match="must have 0 traces, one per interior cuff, not 1"):
        hyperbolic_surfaces.from_signature(orbifold(0, (3, 4, 4)), cuff_traces=(T,))
    with pytest.raises(ValueError, match="the cuff trace 2 must be larger than 3"):
        hyperbolic_surfaces.from_signature(SurfaceSignature(2), cuff_traces=(T, T, 2))


def test_not_a_signature():
    with pytest.raises(TypeError, match="must be a SurfaceSignature"):
        hyperbolic_surfaces.from_signature((2,))
