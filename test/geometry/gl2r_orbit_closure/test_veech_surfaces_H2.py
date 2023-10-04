r"""
Calta-McMullen Veech surfaces in H(2)
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

from flatsurf import translation_surfaces, GL2ROrbitClosure


@pytest.mark.parametrize(
    "w,h,t,e",
    [
        (2, 1, 0, 0),
        (3, 1, 0, 0),
        (3, 1, 0, 1),
        (4, 1, 0, 1),
        (4, 1, 0, 2),
        (3, 2, 0, 0),
        (4, 2, 1, 0),
        (4, 2, 0, 1),
        (4, 2, 1, 1),
    ],
)
def test_H2(w, h, t, e):
    S = translation_surfaces.mcmullen_genus2_prototype(w, h, t, e)
    orbit_closure = GL2ROrbitClosure(S)
    for d in orbit_closure.decompositions(5, 50):
        assert d.parabolic()
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 2
