r"""
Test the names exported by flatsurf.
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


@pytest.mark.parametrize(
    "name, module",
    [
        ("EuclideanPolygons", "flatsurf.geometry.categories"),
        ("HyperbolicPolygons", "flatsurf.geometry.categories"),
        ("EuclideanPolygon", "flatsurf.geometry.polygon"),
        ("HyperbolicPolygon", "flatsurf.geometry.hyperbolic_polygon"),
    ],
)
def test_exported_object(name, module):
    import importlib

    import flatsurf

    assert getattr(flatsurf, name) is getattr(importlib.import_module(module), name)


def test_polygons_is_euclidean_polygons():
    import warnings

    import flatsurf

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Polygons = getattr(flatsurf, "Polygons")

    assert Polygons is flatsurf.EuclideanPolygons


def test_euclidean_membership():
    from sage.all import QQ

    from flatsurf import EuclideanPolygon, EuclideanPolygons, HyperbolicPolygons
    from flatsurf import polygons

    P = polygons.square()
    assert P in EuclideanPolygons(QQ)
    assert isinstance(P, EuclideanPolygon)
    assert P not in HyperbolicPolygons(QQ)


def test_hyperbolic_membership():
    from sage.all import QQ

    from flatsurf import HyperbolicPolygon, HyperbolicPolygons, hyperbolic_polygons

    P = hyperbolic_polygons.regular_ngon(8, QQ(1) / 8)
    assert P in HyperbolicPolygons(P.base_ring())
    assert isinstance(P, HyperbolicPolygon)
