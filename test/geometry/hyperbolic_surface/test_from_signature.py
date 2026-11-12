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

from sage.all import AA, QQ, polygen

from flatsurf import hyperbolic_surfaces
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
        assert S.base_ring().polynomial() == polynomial
    assert len(S.labels()) == polygons


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
        (orbifold(0, (3, 3, 4)), QQ(1) / 12, y**4 - 2, 2),
        (orbifold(0, (2, 4, 5)), QQ(1) / 20, y**4 - y**2 - 1, 2),
        (orbifold(0, (3, 4, 4)), QQ(1) / 6, y**2 - 6, 2),
    ],
)
def test_signature_round_trip(signature, area, polynomial, polygons):
    S = hyperbolic_surfaces.from_signature(signature)
    check(S, signature, area, polynomial, polygons)


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
    [row for row in COMPACT_SIGNATURES if row[:2] != (0, (2, 3, 7))],
)
def test_compact_signatures(g, periods, area):
    signature = orbifold(g, periods)
    S = hyperbolic_surfaces.from_signature(signature)
    assert S.signature() == signature
    assert S.area() == area / 2


def test_ring():
    S = hyperbolic_surfaces.from_signature(SurfaceSignature(2), ring=AA)
    assert S.base_ring() is AA
    assert S.signature() == SurfaceSignature(2)
    assert hyperbolic_surfaces.from_signature(SurfaceSignature(2), ring=QQ).base_ring() is QQ
    with pytest.raises(ValueError, match="isometry is not defined over Rational Field"):
        hyperbolic_surfaces.from_signature(orbifold(2, (2,)), ring=QQ)


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
