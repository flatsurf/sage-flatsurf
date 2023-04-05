r"""
Discriminant loci in H(1,1)
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

from sage.all import polygen, NumberField, AA, QQ
from flatsurf import translation_surfaces, GL2ROrbitClosure


def test_D9_number_field():
    x = polygen(QQ)
    K = NumberField(x**3 - 2, "a", embedding=AA(2) ** QQ((1, 3)))
    a = K.gen()
    S = translation_surfaces.mcmullen_genus2_prototype(2, 1, 0, -1, a / 4)
    orbit_closure = GL2ROrbitClosure(S)
    for d in orbit_closure.decompositions(5):
        ncyl = len(d.cylinders())
        nmin = len(d.minimalComponents())
        nund = len(d.undeterminedComponents())
        assert nund == 0
        assert (nmin == 0) or (ncyl == 0 and 1 <= nmin <= 2)
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == 3
    assert orbit_closure.field_of_definition() == QQ


def test_D9_exact_real():
    pytest.importorskip("pyexactreal")

    from pyexactreal import ExactReals

    R = ExactReals(QQ)
    S = translation_surfaces.mcmullen_genus2_prototype(
        2, 1, 0, -1, R.random_element([0.1, 0.2])
    )
    orbit_closure = GL2ROrbitClosure(S)
    for d in orbit_closure.decompositions(5):
        ncyl = len(d.cylinders())
        nmin = len(d.minimalComponents())
        nund = len(d.undeterminedComponents())
        assert nund == 0
        assert (nmin == 0) or (ncyl == 0 and 1 <= nmin <= 2)
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == 3


def test_D33():
    S = translation_surfaces.mcmullen_genus2_prototype(4, 2, 1, 1, QQ((1, 4)))
    orbit_closure = GL2ROrbitClosure(S)
    for d in orbit_closure.decompositions(5, 100):
        ncyl = len(d.cylinders())
        nmin = len(d.minimalComponents())
        nund = len(d.undeterminedComponents())
        assert nund == 0
        assert (nmin == 0) or (ncyl == 0 and nmin == 2)
        orbit_closure.update_tangent_space_from_flow_decomposition(d)
    assert orbit_closure.dimension() == 3
    K = orbit_closure.field_of_definition()
    assert K.degree() == 2
    assert K.discriminant() == 33
