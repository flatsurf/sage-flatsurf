r"""
The Veech group of a surface and related concepts.

EXAMPLES:

We write an element of the Veech group as an action on homology::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.square_torus()
    sage: A = S.affine_automorphism_group()
    sage: M = matrix([[1, 2], [0, 1]])
    sage: f = A.derivative().section()(M, check=False)

    sage: from flatsurf import SimplicialHomology
    sage: H = SimplicialHomology(S)
    sage: H.hom(f)
    Induced endomorphism of H₁(Translation Surface in H_1(0) built from a square)
      Defn: Induced by Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 2]
                    [0 1]

"""
######################################################################
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
######################################################################

from sage.groups.matrix_gps.matrix_group import MatrixGroup_generic
from sage.categories.morphism import Morphism
from sage.categories.homset import Homset
from sage.categories.morphism import Morphism


class VeechGroup_base(MatrixGroup_generic):
    def affine_automorphism(self):
        # TODO: Turn the matrix into an affine automorphism
        raise NotImplementedError


class AffineAutomorphismGroup_generic(Homset):
    def __init__(self, surface, category=None):
        super().__init__(surface, surface, category=category, base=surface.base_ring())

    def __repr__(self):
        return f"AffineAutomorphismGroup({self.domain()!r})"

    def derivative(self):
        from sage.all import GL
        return DerivativeMap(self, GL(2, self.base_ring()))

    def _from_matrix(self, matrix):
        return self.__make_element_class__(AffineAutomorphismElement_matrix)(self, matrix)


class AffineAutomorphismElement_matrix(Morphism):
    # TODO: Inherit from the generic morphism machinery.
    def __init__(self, parent, matrix):
        super().__init__(parent)

        if matrix.parent().base_ring() is not parent.base_ring():
            raise TypeError("matrix must be over the surface base ring")

        if matrix.is_mutable():
            raise TypeError("matrix must be immutable")

        self._matrix = matrix

    def _repr_type(self):
        return "Affine"

    def _repr_defn(self):
        return f"Lift of linear action given by\n{self._matrix!r}"


class DerivativeMap(Morphism):
    # TODO: From AffineAutomorphismGroup to GL2R where R is the base ring of the surface.
    def section(self):
        return SectionDerivativeMap(self)


class SectionDerivativeMap(Morphism):
    def __init__(self, derivative):
        super().__init__(derivative.parent().reversed())

        self._derivative = derivative

    def section(self):
        return self._derivative

    def _call_(self, matrix):
        return self(matrix, check=True)

    def _call_with_args(self, matrix, args, kwargs):
        if args:
            raise ValueError("unsupported positional argument")

        check = kwargs.pop("check", True)

        if kwargs:
            raise ValueError("unsupported keyword argument")

        if check:
            raise NotImplementedError("cannot check whether a matrix is in the Veech group yet, see https://github.com/flatsurf/sage-flatsurf/issues/157")

        return self.codomain()._from_matrix(matrix.matrix())
        

def AffineAutomorphismGroup(surface):
    return surface.affine_automorphism_group()


def VeechGroup(surface):
    return surface.veech_group()
