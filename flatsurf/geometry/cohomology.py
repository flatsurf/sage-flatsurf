r"""
Absolute and relative (simplicial) cohomology of surfaces.

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

The relative cohomology, relative to the vertices::

    sage: S = S.erase_marked_points()  # optional: pyflatsurf  # random output due to deprecation warnings
    sage: H = SimplicialCohomology(S, relative=S.vertices())  # optional: pyflatsurf
    sage: len(H.gens())  # optional: pyflatsurf
    17

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
from sage.misc.cachefunc import cached_method


class SimplicialCohomologyClass(Element):
    r"""
    A cohomology class.

    INPUT:

    - ``parent`` -- the cohomology group

    - ``values`` -- a dict; the value at each generator of
      :meth:`SimplicialCohomologyGroup.homology`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: S = translation_surfaces.regular_octagon()
        sage: H = SimplicialCohomology(S)

        sage: f, _, _, _ = H.gens()

        sage: from flatsurf.geometry.cohomology import SimplicialCohomologyClass
        sage: isinstance(f, SimplicialCohomologyClass)
        True

    """

    def __init__(self, parent, values):
        super().__init__(parent)

        self._values = values

    def _repr_(self):
        r"""
        Return a printable representation of this cohomology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: S = translation_surfaces.regular_octagon()
            sage: H = SimplicialCohomology(S)

            sage: f, _, _, _ = H.gens()

            sage: f
            {B[(0, 1)]: 1}

        """
        return repr(self._values)

    def __call__(self, homology):
        r"""
        Evaluate this class at an element of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()
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
        Return the pointwise sum of two cohomology classes.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialCohomology(T)

            sage: γ = H.homology().gens()[0]
            sage: f = H({γ: 1.337})
            sage: (f + f)(γ) == 2*f(γ)
            True

        """
        other = self.parent()(other)

        values = {}
        for gen in self.parent().homology().gens():
            value = self(gen) + other(gen)
            if value:
                values[gen] = value

        return self.parent()(values)

    def _richcmp_(self, other, op):
        r"""
        Return how this cohomoly class compares to other with respect to the
        binary relation ``op``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: S = translation_surfaces.regular_octagon()
            sage: H = SimplicialCohomology(S)

            sage: f, g, _, _ = H.gens()

            sage: f == g
            False

        """
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
    The ``k``-th simplicial cohomology group of the ``surface`` with
    ``coefficients``.

    INPUT:

    - ``surface`` -- a finite type surface without boundary

    - ``k`` -- an integer

    - ``coefficients`` -- a ring

    - ``relative`` -- a subset of points of the ``surface``

    - ``implementation`` -- a string; the algorithm used to compute the
      cohomology, only ``"dual`` is supported at the moment.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: T = translation_surfaces.square_torus()

        sage: SimplicialCohomology(T)
        H¹(Translation Surface in H_1(0) built from a square)

    """

    Element = SimplicialCohomologyClass

    def __init__(self, surface, k, coefficients, relative, implementation, category):
        Parent.__init__(self, category=category)

        if surface.is_mutable():
            raise TypeError("surface most be immutable")

        from sage.all import ZZ

        if k not in ZZ:
            raise TypeError("k must be an integer")

        from sage.categories.all import Rings

        if coefficients not in Rings():
            raise TypeError("coefficients must be a ring")

        if implementation == "dual":
            if not surface.is_compact():
                raise NotImplementedError(
                    "dual implementation can only handle cohomology of compact surfaces"
                )

            if surface.is_with_boundary():
                raise NotImplementedError(
                    "dual implementation can only handle cohomology of surfaces without boundary"
                )

            if coefficients.characteristic() > 0:
                raise NotImplementedError(
                    "dual implementation can only handle cohomology with coefficients of characteristic zero"
                )

            if not coefficients.is_integral_domain():
                raise NotImplementedError(
                    "dual implementation can only handle cohomology with flat coefficient rings"
                )

        else:
            raise NotImplementedError("unknown implementation for cohomology group")

        self._surface = surface
        self._k = k
        self._coefficients = coefficients
        self._relative = relative

    def is_absolute(self):
        r"""
        Return whether this is an absolute cohomology (and not a relative one).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()

            sage: H = SimplicialCohomology(T)
            sage: H.is_absolute()
            True

        """
        return not self._relative

    def _repr_(self):
        r"""
        Return a printable representation of this cohomology group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()

            sage: SimplicialCohomology(T)
            H¹(Translation Surface in H_1(0) built from a square)


            sage: SimplicialCohomology(T, relative=T.vertices())
            H¹(Translation Surface in H_1(0) built from a square, {Vertex 0 of polygon 0})


            sage: SimplicialCohomology(T, coefficients=CC, relative=T.vertices())
            H¹(Translation Surface in H_1(0) built from a square, {Vertex 0 of polygon 0}; Complex Field with 53 bits of precision)


        """
        if self._k == 0:
            k = "⁰"
        elif self._k == 1:
            k = "¹"
        elif self._k == 2:
            k = "²"
        else:
            k = f"^{self._k}"

        Hk = f"H{k}"

        X = repr(self.surface())
        if not self.is_absolute():
            X = f"{X}, {set(self._relative)}"

        from sage.all import RR

        if self._coefficients is not RR:
            sep = ";"
            X = f"{X}{sep} {self._coefficients}"

        return f"{Hk}({X})"

    def surface(self):
        r"""
        Return the surface over which this cohomology is defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()

            sage: H = SimplicialCohomology(T)
            sage: H.surface() is T
            True

        """
        return self._surface

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as a class in this cohomology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()

            sage: H = SimplicialCohomology(T)

            sage: H(0)
            {}
            sage: H({gen: 1 for gen in H.homology().gens()})
            {B[(0, 1)]: 1, B[(0, 0)]: 1}

        """
        if not x:
            x = {}

        if isinstance(x, dict):
            x = {self.homology()(gen): value for (gen, value) in x.items() if value}
            return self.element_class(self, x)

        raise NotImplementedError("cannot create a cohomology class from this data")

    @cached_method
    def homology(self):
        r"""
        Return the homology of the underlying space (with integer
        coefficients).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialCohomology(T)
            sage: H.homology()
            H₁(Translation Surface in H_1(0) built from a square)

        """
        from flatsurf.geometry.homology import SimplicialHomology

        return SimplicialHomology(self._surface, self._k, relative=self._relative)

    def gens(self):
        r"""
        Return generators of this cohomology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialCohomology(T)
            sage: H.gens()
            [{B[(0, 1)]: 1}, {B[(0, 0)]: 1}]

        """
        return [self({gen: 1}) for gen in self.homology().gens()]


def SimplicialCohomology(
    surface, k=1, coefficients=None, relative=None, implementation="dual", category=None
):
    r"""
    Return the ``k``-th simplicial cohomology group of ``surface``.

    INPUT:

    - ``surface`` -- a surface

    - ``k`` -- an integer (default: ``1``)

    - ``coefficients`` -- a ring (default: the reals); consider cohomology with
      coefficients in this ring

    - ``relative`` -- a set (default: the empty set); if non-empty, then
      relative cohomology with respect to this set is constructed.

    - ``implementation`` -- a string (default: ``"dual"``); the algorithm used
      to compute the cohomology groups. Currently only ``"dual"`` is supported,
      i.e., the groups are computed as duals of the generic homology groups
      from SageMath.

    - ``category`` -- a category; if not specified, a category for the
      cohomology group is chosen automatically depending on ``coefficients``.

    TESTS:

    Cohomology is unique and cached::

        sage: from flatsurf import translation_surfaces, SimplicialCohomology
        sage: T = translation_surfaces.square_torus()
        sage: SimplicialCohomology(T) is SimplicialCohomology(T)
        True

    """
    return surface.cohomology(k, coefficients, relative, implementation, category)
