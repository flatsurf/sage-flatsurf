r"""
Internal consistency of the reference values of hyperbolic surfaces.
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

from fractions import Fraction

import pytest

from hyperbolic_surface_oracles import (
    BUSER_DIMENSIONS,
    COMPACT_SIGNATURES,
    CUBE_DIVISOR,
    CUSPED_SIGNATURES,
    FLAT_SIGNATURES,
)


@pytest.mark.parametrize("g,periods,area_over_pi", COMPACT_SIGNATURES)
def test_compact_area(g, periods, area_over_pi):
    assert 2 * (2 * g - 2 + sum(1 - Fraction(1, m) for m in periods)) == area_over_pi
    assert -2 * (2 - 2 * g + sum(Fraction(1, m) - 1 for m in periods)) == area_over_pi


@pytest.mark.parametrize("g,periods", FLAT_SIGNATURES)
def test_flat_signature(g, periods):
    assert 2 - 2 * g + sum(Fraction(1, m) - 1 for m in periods) == 0


@pytest.mark.parametrize("g,periods,cusps,area_over_pi", CUSPED_SIGNATURES)
def test_cusped_area(g, periods, cusps, area_over_pi):
    assert (
        2 * (2 * g - 2 + cusps + sum(1 - Fraction(1, m) for m in periods))
        == area_over_pi
    )
    # A cusp is a cone point of angle 0.
    chi = 2 - 2 * g + sum(Fraction(1, m) - 1 for m in periods) - cusps
    assert -2 * chi == area_over_pi


def test_cube_divisor():
    chi = CUBE_DIVISOR["euler_characteristic"]
    assert chi == 2
    assert (
        chi + CUBE_DIVISOR["cone_points"] * (CUBE_DIVISOR["angle"] - 1)
        == CUBE_DIVISOR["euler_characteristic_with_divisor"]
        == 0
    )


@pytest.mark.parametrize("g,p,q,y_pieces,pairs,dimension", BUSER_DIMENSIONS)
def test_buser_dimensions(g, p, q, y_pieces, pairs, dimension):
    assert y_pieces == 2 * g - 2 + p + q
    assert pairs == 3 * g - 3 + p + q
    assert dimension == 6 * g - 6 + 3 * p + 2 * q
    assert dimension == 2 * pairs + p
