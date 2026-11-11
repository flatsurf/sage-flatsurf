r"""
Cone angles, orbifold points and Euler characteristic of the flat orbifolds
built by the surface generators.
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


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_angles(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    observed = sorted(p.angle() for p in S.vertices())
    assert observed == angles

    # A cone point of angle 1/m has period m; the torus has none.
    orders = sorted(int(1 / a) for a in observed if a != 1)
    assert orders == sorted(signature[1])


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_points(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    assert S.is_orbifold()
    # None of these surfaces has a self-glued edge: every orbifold point is a
    # vertex.
    points = S.orbifold_points()
    assert all(p.is_vertex() for p, m in points)
    assert sorted(m for p, m in points) == sorted(signature[1])


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_chi(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    assert S.euler_characteristic() == chi


@pytest.mark.parametrize("name,signature,angles,chi", FLAT_ORBIFOLDS)
def test_flat_orbifold_signature(name, signature, angles, chi):
    S = CONSTRUCTORS[name]()
    g, periods = signature
    assert S.signature().katok() == (g, tuple(sorted(periods)), 0)
    assert S.signature().euler_characteristic() == 0
