r"""
Discriminant loci in H(1,1)
"""
# ****************************************************************************
# This file is part of sage-flatsurf.
#
#       Copyright (C) 2020 Vincent Delecroix
#
# sage-flatsurf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# sage-flatsurf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************

import pytest

pytest.importorskip("pyflatsurf")  # noqa

from flatsurf import translation_surfaces
from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf


def test_origami1():
    from sage.all import SymmetricGroup

    G = SymmetricGroup(2)
    r = u = G("(1,2)")
    origami = translation_surfaces.origami(r, u)
    S = to_pyflatsurf(origami)
    assert (
        str(S)
        == "FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, -5, 4)(-2, -6, 5, -4, 6, 3), faces = (1, 2, 3)(-1, 4, 5)(-2, -3, 6)(-4, -5, -6)) with vectors {1: (-1, -1), 2: (1, 0), 3: (0, 1), 4: (-1, 0), 5: (0, -1), 6: (1, 1)}"
    )


def test_origami2():
    from sage.all import SymmetricGroup

    G = SymmetricGroup(3)
    r = G("(1,2,3)")
    u = G("(1,2)")
    origami = translation_surfaces.origami(r, u)
    S = to_pyflatsurf(origami)
    assert (
        str(S)
        == "FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, -5, -9, 8, 5, -4, 6, 3, -2, -6, -7, 9, -8, 7, 4), faces = (1, 2, 3)(-1, 4, 5)(-2, -3, 6)(-4, 7, -6)(-5, 8, 9)(-7, -8, -9)) "
        "with vectors {1: (-1, -1), 2: (1, 0), 3: (0, 1), 4: (-1, 0), 5: (0, -1), 6: (1, 1), 7: (0, 1), 8: (-1, -1), 9: (1, 0)}"
    )


@pytest.mark.parametrize("n", [3, 5, 7])
def test_regular_n_gons(n):
    S = translation_surfaces.veech_double_n_gon(n)
    to_pyflatsurf(S)


@pytest.mark.parametrize("g", [3, 4])
def test_arnoux_yoccoz(g):
    A = translation_surfaces.arnoux_yoccoz(g)
    to_pyflatsurf(A)


def test_ward3():
    W3 = translation_surfaces.ward(3)
    to_pyflatsurf(W3)

    W17 = translation_surfaces.ward(17)
    to_pyflatsurf(W17)
