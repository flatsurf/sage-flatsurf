r"""
The gothic locus

From the article McMullen-Mukamel-Wright (2017).
"""
######################################################################
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
######################################################################

import pytest

pytest.importorskip("pyflatsurf")  # noqa

from sage.all import QQ, AA, NumberField, polygen
from flatsurf import translation_surfaces, GL2ROrbitClosure
from surface_dynamics import Stratum


def test_gothic_generic():
    x = polygen(QQ)
    K = NumberField(x**3 - 2, "a", embedding=AA(2) ** QQ((1, 3)))
    a = K.gen()
    S = translation_surfaces.cathedral(a, a**2)
    orbit_closure = GL2ROrbitClosure(S)
    assert orbit_closure.ambient_stratum() == Stratum([2, 2, 2], 1)
    for d in orbit_closure.decompositions(4, 50):
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 4
    assert orbit_closure.field_of_definition() == QQ


def test_gothic_exact_reals():
    pytest.importorskip("pyexactreal")

    from pyexactreal import ExactReals

    x = polygen(QQ)
    K = NumberField(x**3 - 2, "a", embedding=AA(2) ** QQ((1, 3)))
    R = ExactReals(K)
    S = translation_surfaces.cathedral(K.gen(), R.random_element([0.1, 0.2]))
    orbit_closure = GL2ROrbitClosure(S)
    assert orbit_closure.ambient_stratum() == Stratum([2, 2, 2], 1)
    for d in orbit_closure.decompositions(4, 50):
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 4


def test_gothic_veech():
    x = polygen(QQ)
    K = NumberField(x**2 - 2, "sqrt2", embedding=AA(2) ** QQ((1, 2)))
    sqrt2 = K.gen()
    x = QQ((1, 2))
    y = 1
    a = x + y * sqrt2
    b = -3 * x - QQ((3, 2)) + 3 * y * sqrt2
    S = translation_surfaces.cathedral(a, b)
    orbit_closure = GL2ROrbitClosure(S)
    assert orbit_closure.ambient_stratum() == Stratum([2, 2, 2], 1)
    for d in orbit_closure.decompositions(4, 50):
        assert d.parabolic()
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 2
    assert (
        orbit_closure.field_of_definition()
        == orbit_closure.V2._isomorphic_vector_space.base_ring()
    )
