r"""
TODO: Document this module.
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
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


class SimplicialCohomology(UniqueRepresentation, Parent):
    r"""
    Absolute simplicial cohomology of the ``surface`` with ``coefficients``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
        sage: T.set_immutable()

    Currently, surfaces must be Delaunay triangulated to compute their homology::

        sage: SimplicialCohomology(T)
        H¹(Translation Surface in H_1(0) built from 2 isosceles triangles; Real Field with 53 bits of precision)

    """
    Element = SimplicialCohomologyClass

    @staticmethod
    def __classcall__(cls, surface, coefficients=None, homology=None, category=None):
        r"""
        Normalize parameters used to construct cohomology.

        TESTS:

        Cohomology is unique and cached::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: SimplicialCohomology(T) is SimplicialCohomology(T)
            True

        """
        from sage.all import RR
        from flatsurf.geometry.homology import SimplicialHomology
        return super().__classcall__(cls, surface, coefficients or RR, homology or SimplicialHomology(surface), category or SetsWithPartialMaps())

    def __init__(self, surface, coefficients, homology, category):
        # TODO: Not checking this anymore. Do we need it?
        # if surface != surface.delaunay_triangulation():
        #     # TODO: This is a silly limitation in here.
        #     raise NotImplementedError("Surface must be Delaunay triangulated")

        Parent.__init__(self, category=category)

        self._surface = surface
        self._homology = homology
        self._coefficients = coefficients

    def _repr_(self):
        return f"H¹({self._surface}; {self._coefficients})"

    def _element_constructor_(self, x):
        if not x:
            x = {}

        if isinstance(x, dict):
            assert all(k in self.homology().gens() for k in x.keys())
            return self.element_class(self, x)

        raise NotImplementedError

    @cached_method
    def homology(self):
        r"""
        Return the homology of the underlying space (with integer
        coefficients.)

        EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
        sage: T.set_immutable()
        sage: H = SimplicialCohomology(T)
        sage: H.homology()
        H₁(Translation Surface in H_1(0) built from 2 isosceles triangles; Integer Ring)

        """
        return self._homology
