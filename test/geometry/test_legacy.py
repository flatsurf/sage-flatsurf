r"""
Tests that some legacy constructions still work.

We test this here instead of in doctests because these will throw lots of
deprecation warnings that are somewhat annoying to doctest for.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
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


def test_apisa_wright():
    r"""
    Test that defining a surface as exemplified in apisa_wright notebook still
    works after #220.
    """
    from sage.all import Sequence, vector, QuadraticField

    K = QuadraticField(2)
    c = K.gen()
    h24, h3, l15, l6, l7, l8 = 1, 1 + c, 1, c, 1 + c, 2 * c - 1

    from flatsurf import ConvexPolygons, HalfTranslationSurface, Surface_list

    params = [h24, h3, l15, l6, l7, l8]
    K = Sequence(params).universe().fraction_field()
    v24 = vector(K, (0, h24))
    v3 = vector(K, (0, h3))
    v15 = vector(K, (l15, 0))
    v6 = vector(K, (l6, 0))
    v7 = vector(K, (l7, 0))
    v8 = vector(K, (l8, 0))
    C = ConvexPolygons(K)
    P0 = C(edges=[v15, v15, v24, -2 * v15, -v24])
    P1 = C(edges=[2 * v15, v8, v7, v6, v3, -v8, -v7, -v6, -v15, -v15, -v3])
    P2 = C(edges=[v15, v24, -v15, -v24])
    S = Surface_list(base_ring=C.base_ring())
    S.rename("ApisaWrightSurface({})".format(", ".join(map(str, params))))
    S.add_polygons([P0, P1, P2])
    # set_edge_pairing(poly_num1, edge_num1, poly_num2, edge_num2)
    S.set_edge_pairing(0, 0, 0, 1)
    S.set_edge_pairing(0, 2, 0, 4)
    S.set_edge_pairing(0, 3, 1, 0)
    S.set_edge_pairing(1, 1, 1, 5)
    S.set_edge_pairing(1, 2, 1, 6)
    S.set_edge_pairing(1, 3, 1, 7)
    S.set_edge_pairing(1, 4, 1, 10)
    S.set_edge_pairing(1, 8, 2, 2)
    S.set_edge_pairing(1, 9, 2, 0)
    S.set_edge_pairing(2, 1, 2, 3)
    S.set_immutable()
    HalfTranslationSurface(S)
