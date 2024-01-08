r"""
Veech 2n-gons and double n-gons

These two families of examples are from the article of Veech "Teichmüller
curves in moduli space, Eisenstein series and an application to triangular
billiards" (1989). These are examples of what now is called a Veech surface.
These surface satisfies the so-called Veech dichotomy: in any direction of a
saddle connection, the flow is completely periodic and the cylinders have
commensurable moduli ("parabolic direction").
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


@pytest.mark.parametrize("n,bound", [(4, 4), (5, 4), (6, 4), (7, 4), (8, 4)])
def test_veech_2n_gon(n, bound):
    S = translation_surfaces.veech_2n_gon(n)
    orbit_closure = GL2ROrbitClosure(S)
    for d in orbit_closure.decompositions(bound):
        assert d.parabolic()
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 2


@pytest.mark.parametrize("n,bound", [(3, 4), (5, 4), (7, 4), (9, 4)])
def test_veech_double_n_gon(n, bound):
    S = translation_surfaces.veech_double_n_gon(n)
    orbit_closure = GL2ROrbitClosure(S)
    for d in orbit_closure.decompositions(bound):
        assert d.parabolic()
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 2
