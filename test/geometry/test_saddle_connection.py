r"""
Tests that saddle connections are enumerated correctly.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2024 Julian Rüth
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


def test_L():
    r"""
    Test that saddle connections in an L are counted the same no matter how the
    L is represented.
    """
    from flatsurf import translation_surfaces

    L_from_rectangles = translation_surfaces.mcmullen_L(1, 2, 3, 4)
    connections_from_rectangles = L_from_rectangles.saddle_connections(128)
    assert len(connections_from_rectangles) == 164

    L_from_triangles = L_from_rectangles.triangulate().codomain()
    connections_from_triangles = L_from_triangles.saddle_connections(128)
    assert len(connections_from_triangles) == 164

    from sage.all import QQ
    from flatsurf import MutableOrientedSimilaritySurface, Polygon
    L = MutableOrientedSimilaritySurface(QQ)
    L.add_polygon(Polygon(vertices=[(0, 0), (3, 0), (7, 0), (7, 2), (3, 2), (3, 3), (0, 3), (0, 2)]))
    L.glue((0, 0), (0, 5))
    L.glue((0, 1), (0, 3))
    L.glue((0, 2), (0, 7))
    L.glue((0, 4), (0, 6))
    L.set_immutable()

    connections = L.saddle_connections(128)

    assert len(connections) == 164
