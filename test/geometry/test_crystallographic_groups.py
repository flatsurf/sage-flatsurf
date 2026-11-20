r"""
Tests for the Euclidean triangle groups.
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

from sage.all import AA, QQ, RR, matrix

from flatsurf import crystallographic_groups
from flatsurf.geometry.similarity import SimilarityGroup

PERIODS = [
    ((2, 4, 4), True),
    ((4, 2, 4), True),
    ((4, 4, 2), True),
    ((3, 3, 3), False),
    ((2, 3, 6), False),
    ((6, 3, 2), False),
    ((3, 6, 2), False),
    ((2, 6, 3), False),
]


@pytest.mark.parametrize("periods,rational", PERIODS)
def test_orders_and_product(periods, rational):
    G = crystallographic_groups.triangle(*periods)
    assert G.ngens() == 3
    assert G.is_orientable()
    gens = G.gens()
    for g, m in zip(gens, periods):
        assert (g**m).is_one()
        assert not any((g**k).is_one() for k in range(1, m))
    assert (gens[0] * gens[1] * gens[2]).is_one()
    if rational:
        assert G.base_ring() is QQ
    else:
        assert G.base_ring().degree() == 2
        assert G.base_ring().polynomial().list() == [-3, 0, 1]


@pytest.mark.parametrize("periods,rational", PERIODS)
def test_reflections(periods, rational):
    W = crystallographic_groups.triangle(*periods, orientation_preserving=False)
    assert W.ngens() == 3
    assert all(g.is_reflection() for g in W.gens())
    assert not W.is_orientable()
    assert W.orientation_preserving_subgroup().ngens() == 2


def test_generators_244():
    F = SimilarityGroup(QQ)
    G = crystallographic_groups.triangle(2, 4, 4)
    expected = [
        matrix(QQ, 3, [-1, 0, 0, 0, -1, 0, 0, 0, 1]),
        matrix(QQ, 3, [0, -1, 1, 1, 0, 1, 0, 0, 1]),
        matrix(QQ, 3, [0, -1, -1, 1, 0, 1, 0, 0, 1]),
    ]
    assert list(G.gens()) == [F(m) for m in expected]


def test_generators_424():
    F = SimilarityGroup(QQ)
    G = crystallographic_groups.triangle(4, 2, 4)
    expected = [
        matrix(QQ, 3, [0, -1, 0, 1, 0, 0, 0, 0, 1]),
        matrix(QQ, 3, [-1, 0, 2, 0, -1, 2, 0, 0, 1]),
        matrix(QQ, 3, [0, -1, 2, 1, 0, 2, 0, 0, 1]),
    ]
    assert list(G.gens()) == [F(m) for m in expected]


def test_algebraic_real_field():
    G = crystallographic_groups.triangle(2, 3, 6, base_ring=AA)
    assert G.base_ring() is AA
    g1, g2, g3 = G.gens()
    assert (g1 * g2 * g3).is_one()


def test_refusals():
    with pytest.raises(ValueError, match="the period 1 is not an integer at least 2"):
        crystallographic_groups.triangle(1, 4, 4)
    with pytest.raises(ValueError, match=r"\(2, 3, 7\) is hyperbolic, not euclidean"):
        crystallographic_groups.triangle(2, 3, 7)
    with pytest.raises(ValueError, match=r"\(2, 3, 5\) is spherical, not euclidean"):
        crystallographic_groups.triangle(2, 3, 5)
    with pytest.raises(ValueError, match="do not lie in Rational Field"):
        crystallographic_groups.triangle(2, 3, 6, base_ring=QQ)
    with pytest.raises(ValueError, match="the ring must be exact"):
        crystallographic_groups.triangle(2, 4, 4, base_ring=RR)


def test_repr_and_import():
    import flatsurf

    assert repr(crystallographic_groups) == "crystallographic groups"
    assert flatsurf.crystallographic_groups is crystallographic_groups
