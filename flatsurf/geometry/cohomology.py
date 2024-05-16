r"""
Absolute (simplicial) cohomology of surfaces.

EXAMPLES:

The absolute cohomology of the regular octagon::

    sage: from flatsurf import translation_surfaces, SimplicialCohomology
    sage: S = translation_surfaces.regular_octagon()
    sage: H = SimplicialCohomology(S)

A basis of cohomology::

    sage: H.gens()
    [{B[(0, 1)]: 1}, {B[(0, 2)]: 1}, {B[(0, 3)]: 1}, {B[(0, 0)]: 1}]

The absolute cohomology of the unfolding of the (3, 4, 13) triangle::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(angles=[3, 4, 13])
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")
    sage: H = SimplicialCohomology(S)
    sage: len(H.gens())
    16

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2024 Julian Rüth
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
######################################################################

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method


class SimplicialCohomologyClass(Element):
    def __init__(self, parent, values):
        super().__init__(parent)

        self._values = values

    def _repr_(self):
        return repr(self._values)

    def __call__(self, homology):
        r"""
        Evaluate this class at an element of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: H = SimplicialCohomology(T)

            sage: γ = H.homology().gens()[0]
            sage: f = H({γ: 1.337})
            sage: f(γ)
            1.33700000000000

        """
        return sum(
            self._values.get(gen, 0) * homology.coefficient(gen)
            for gen in self.parent().homology().gens()
        )

    def _add_(self, other):
        r"""
        Return the pointwise sum of two homology classes.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: H = SimplicialCohomology(T)

            sage: γ = H.homology().gens()[0]
            sage: f = H({γ: 1.337})
            sage: (f + f)(γ) == 2*f(γ)
            True

        """
        values = {}
        for gen in self.parent().homology().gens():
            value = self(gen) + other(gen)
            if value:
                values[gen] = value

        return self.parent()(values)

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if self is other:
                return True

            if self.parent() != other.parent():
                return False

            return self._values == other._values

        return super()._richcmp_(other, op)


class SimplicialCohomologyGroup(Parent):
    r"""
    Absolute simplicial cohomology of the ``surface`` with ``coefficients``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))

        sage: SimplicialCohomology(T)
        H¹(Translation Surface in H_1(0) built from a square; Real Field with 53 bits of precision)

    """
    Element = SimplicialCohomologyClass

    def __init__(self, surface, k, coefficients, implementation, category):
        Parent.__init__(self, category=category)

        if surface.is_mutable():
            raise TypeError("surface most be immutable")

        from sage.all import ZZ

        if k not in ZZ:
            raise TypeError("k must be an integer")

        from sage.categories.all import Rings

        if coefficients not in Rings():
            raise TypeError("coefficients must be a ring")

        if implementation != "dual":
            raise NotImplementedError("unknown implementation for cohomology group")

        self._surface = surface
        self._k = k
        self._coefficients = coefficients

    def _repr_(self):
        k = self._k
        if k == 0:
            k = "⁰"
        elif k == 1:
            k = "¹"
        elif k == 2:
            k = "²"
        else:
            k = f"^{k}"

        return f"H{k}({self._surface}; {self._coefficients})"

    def surface(self):
        return self._surface

    def _element_constructor_(self, x):
        if not x:
            x = {}

        if isinstance(x, dict):
            assert all(k in self.homology().gens() for k in x.keys())
            x = {gen: value for (gen, value) in x.items() if value}
            return self.element_class(self, x)

        raise ValueError("cannot create a cohomology class from this data")

    @cached_method
    def homology(self):
        r"""
        Return the homology of the underlying space (with integer
        coefficients.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: H = SimplicialCohomology(T)
            sage: H.homology()
            H₁(Translation Surface in H_1(0) built from a square)

        """
        from flatsurf.geometry.homology import SimplicialHomology

        return SimplicialHomology(self._surface, self._k)

    def gens(self):
        return [self({gen: 1}) for gen in self.homology().gens()]


def SimplicialCohomology(
    surface, k=1, coefficients=None, implementation="dual", category=None
):
    r"""
    TESTS:

    Cohomology is unique and cached::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: SimplicialCohomology(T) is SimplicialCohomology(T)
        True

    """
    return surface.cohomology(k, coefficients, implementation, category)
