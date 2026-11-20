r"""
Pickles of surfaces made by release 0.8.0 of sage-flatsurf still load.

The files in ``released_pickles/0.8.0`` are pickles written by that release, by
a script outside the tests. They are kept as bytes, not rebuilt in the tests:
pickles stored by users are loaded with plain ``pickle.loads``, which finds a
class by the name of its module and of its class as they were in the release.
"""

# ####################################################################
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
# ####################################################################

import os
import pickle

import pytest

from sage.all import AA, QQ, NumberField, polygen

from flatsurf import (
    MutableOrientedSimilaritySurface,
    Polygon,
    half_translation_surfaces,
    polygons,
    similarity_surfaces,
)
from flatsurf.geometry.categories import FlatConeSurfaces

PICKLES = os.path.join(os.path.dirname(__file__), "released_pickles", "0.8.0")

DEPRECATED = "cone_surfaces.ConeSurfaces is deprecated"


def load(name, warns=True):
    with open(os.path.join(PICKLES, name + ".pickle"), "rb") as f:
        data = f.read()

    if not warns:
        return pickle.loads(data)

    with pytest.warns(UserWarning, match=DEPRECATED):
        return pickle.loads(data)


def torus(immutable):
    S = MutableOrientedSimilaritySurface(QQ)
    S.add_polygon(polygons.square())
    S.glue((0, 0), (0, 2))
    S.glue((0, 1), (0, 3))
    if immutable:
        S.set_immutable()
    return S


def irrational_triangle():
    x = polygen(QQ)
    K = NumberField(x**3 - 2, "a", embedding=AA(2) ** (QQ(1) / 3))
    return Polygon(vertices=[(0, 0), (1, 0), (0, K.gen())])


def over_ring(surface, ring):
    r"""
    Return ``surface`` with its polygons converted into ``ring``, an immutable
    surface with the same category.

    A number field pickled by another version of SageMath unpickles to a new
    field: its factory caches fields by the version of SageMath and by their
    defining data, and number fields compare by identity. A surface loaded from
    such a pickle is then over a field that differs from the one of a surface
    built fresh, though both fields have the same defining polynomial and
    embedding, and coerce into each other.
    """
    if surface.base_ring() is ring:
        return surface

    assert ring.has_coerce_map_from(surface.base_ring())
    assert surface.base_ring().has_coerce_map_from(ring)

    copy = MutableOrientedSimilaritySurface.from_surface(
        surface.change_ring(ring), category=surface.category()
    )
    copy.set_immutable()
    return copy


# The name of the pickle and the surface it was made from. Each of these
# pickles finds the deprecated name of the category of flat cone surfaces.
CASES = {
    "square_torus": lambda: torus(immutable=True),
    "billiard_irrational": lambda: similarity_surfaces.billiard(
        irrational_triangle()
    ),
    "billiard_rational": lambda: similarity_surfaces.billiard(
        polygons.triangle(2, 3, 4)
    ),
    "half_translation": lambda: half_translation_surfaces.step_billiard(
        [1, 1, 1], [QQ(1) / 2, 1, 2]
    ),
}


@pytest.mark.parametrize("name", CASES)
def test_immutable_pickle(name):
    S = load(name)
    T = over_ring(CASES[name](), S.base_ring())

    assert S == T
    assert hash(S) == hash(T)
    assert S in FlatConeSurfaces()

    if name != "half_translation":
        assert S.category() == T.category()
    # The category of a pickled half-translation surface is a join with the
    # category of half-translation surfaces, whereas a new surface has the
    # category itself; this is so independently of the category of cone
    # surfaces.


def test_mutable_pickle():
    S = load("square_torus_mutable", warns=False)

    assert S == torus(immutable=False)
    assert S.category() == torus(immutable=False).category()
    assert S.is_mutable()

    S.add_polygon(polygons.square(), label=1)
    S.glue((1, 0), (1, 2))
    S.glue((1, 1), (1, 3))
    assert S.polygon(1) == polygons.square()
    assert S.opposite_edge(1, 0) == (1, 2)
