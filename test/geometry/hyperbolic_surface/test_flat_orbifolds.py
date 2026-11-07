r"""
Cone angles and Euler characteristic of the flat orbifolds built by the
surface generators.
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

from flatsurf import Polygon, polygons, similarity_surfaces, translation_surfaces

from hyperbolic_surface_oracles import FLAT_ORBIFOLDS

CONSTRUCTORS = {
    "square_torus": lambda: translation_surfaces.square_torus(),
    "pillowcase": lambda: similarity_surfaces.polygon_double(polygons.square()),
    "equilateral_double": lambda: similarity_surfaces.polygon_double(
        Polygon(angles=(1, 1, 1))
    ),
    "isosceles_right_double": lambda: similarity_surfaces.polygon_double(
        Polygon(angles=(1, 1, 2))
    ),
    "triangle_236_double": lambda: similarity_surfaces.polygon_double(
        Polygon(angles=(1, 2, 3))
    ),
}


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_angles(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    observed = sorted(p.angle() for p in S.vertices())
    assert observed == angles

    # A cone point of angle 1/m has period m; the torus has none.
    orders = sorted(int(1 / a) for a in observed if a != 1)
    assert orders == sorted(signature[1])


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_chi(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    assert S.euler_characteristic() == chi
