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
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.morphism import (
    SurfaceMorphism,
    SurfaceMorphism_factorization,
    MorphismSpace,
)


class VeechGroup_generic(MatrixGroup_generic):
    r"""
    A generic (currently essentially empty) implementation of the Veech group.

    You should not create such a group directly. Use :func:`VeechGroup` or
    :meth:`flatsurf.geometry.categories.dilation_surfaces.DilationSurfaces.ParentMethods.veech_group`.

    .. NOTE:

        See https://github.com/flatsurf/sage-flatsurf/issues/157 for a related
        issue to actually add some functionality to this object.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: V = S.veech_group()

    TESTS::

        sage: from flatsurf.geometry.veech_group import VeechGroup_generic
        sage: isinstance(V, VeechGroup_generic)
        True

        sage: TestSuite(V).run()

    """

    def __init__(self, surface, category=None):
        self._surface = surface

        from sage.all import ZZ

        super().__init__(ZZ(2), surface.base_ring(), category=category)

    def surface(self):
        r"""
        The surface for which this is the Veech group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.veech_group().surface() is S
            True

        """
        return self._surface

    def _repr_(self):
        r"""
        Return a printable representation of this group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: V = S.veech_group()
            sage: V
            VeechGroup(Translation Surface in H_1(0) built from a square)

        """
        return f"VeechGroup({self._surface!r})"

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of the Veech group from ``x``.

        INPUT:

        - ``x`` -- a 2×2 matrix

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: V = S.veech_group()
            sage: V(matrix([[1, 0], [0, 1]]))
            [1 0]
            [0 1]

        """
        from sage.all import GL

        x = GL(2, self.base_ring())(x)
        x = super()._element_constructor_(x)

        if check:
            if x not in self:
                raise ValueError("this is not an element of the Veech group")

        return x

    def _an_element_(self):
        r"""
        Return an element of the Veech group (for testing).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: V = S.veech_group()
            sage: V.an_element()
            [1 0]
            [0 1]

        """
        from sage.all import identity_matrix

        return self(identity_matrix(2))

    def gens(self):
        r"""
        Return generators of this group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: V = S.veech_group()
            sage: V.gens()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot compute generators of the Veech group yet

        """
        raise NotImplementedError("cannot compute generators of the Veech group yet")

    def __contains__(self, x):
        r"""
        Return whether ``x`` is contained in the Veech group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: V = S.veech_group()
            sage: m = matrix([[1, 0], [0, 1]])

            sage: m in V
            True

        """
        # Unfortunately, MatrixGroup_generic (and therefore GL) does not
        # implement _richcmp_ nicely. Therefore, we cannot rely on the builtin
        # implementation of __contains__, i.e., the problem is that
        # GL(2, QQ)(M) != M

        from sage.all import GL

        x = GL(2, self.base_ring())(x)

        if x.is_one():
            return True

        raise NotImplementedError(
            "cannot decide whether a matrix is in the Veech group yet"
        )

    def __eq__(self, other):
        r"""
        Return whether this group is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.veech_group() == S.veech_group()
            True

        """
        if not isinstance(other, VeechGroup_generic):
            return False

        return self._surface == other._surface

    def __hash__(self):
        r"""
        Return a hash value for this group that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: hash(S.veech_group()) == hash(S.veech_group())
            True

        """
        return hash(self._surface)


class AffineAutomorphismGroup_generic(MorphismSpace):
    r"""
    A generic implementation of the group of affine automorphisms of a surface.

    You should not create such a group directly. Use
    :func:`AffineAutomorphismGroup` or
    :meth:`flatsurf.geometry.categories.dilation_surfaces.DilationSurfaces.ParentMethods.affine_automorphism_group`.

    .. NOTE::

        Currently, this group is mostly empty and just a container to hold on
        to the :meth:`derivative` method.

        We do not provide a ``category`` parameter here because it makes
        serialization overly complicated. (And we do not think that anybody is
        going to use it.)

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()

    TESTS::

        sage: from flatsurf.geometry.veech_group import AffineAutomorphismGroup_generic
        sage: isinstance(A, AffineAutomorphismGroup_generic)
        True

        sage: import cppyy  # optional: pyflatsurf  # random output due to cppyy deprecation warnings
        sage: TestSuite(A).run()  # optional: pyflatsurf

    """

    def __init__(self, surface):
        self._surface = surface

        super().__init__(surface, surface, base=surface.base_ring())

    def surface(self):
        r"""
        The surface for which this is the affine automorphism group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.affine_automorphism_group().surface() is S
            True

        """
        return self._surface

    def _an_element_(self):
        r"""
        Return a morphism in this group (for testing).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.affine_automorphism_group().an_element()
            Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 0]
                    [0 1]

        """
        from sage.all import identity_matrix

        return self.derivative().section()(identity_matrix(2))

    def __repr__(self):
        r"""
        Return a printable repreentation of this group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A
            AffineAutomorphismGroup(Translation Surface in H_1(0) built from a square)

        """
        return f"AffineAutomorphismGroup({self.domain()!r})"

    def derivative(self):
        r"""
        Return the derivative of this group, i.e., the map to
        `GL_2(\mathbb{R})` whose image is the :func:`VeechGroup`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative()
            Derivative morphism:
              From: AffineAutomorphismGroup(Translation Surface in H_1(0) built from a square)
              To:   General Linear Group of degree 2 over Rational Field

        """
        from sage.all import GL, Hom

        codomain = GL(2, self.base_ring())
        parent = Hom(self, codomain)
        return parent.__make_element_class__(DerivativeMap)(self, codomain)

    def _from_matrix(self, matrix):
        r"""
        Return an affine automorphism whose :meth:`derivative` is ``matrix``.

        This is a helper method for :meth:`SectionDerivativeMap._call_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative().section()(matrix([[1, 0], [0, 1]]), check=False)
            Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 0]
                    [0 1]

        """
        return self.__make_element_class__(AffineAutomorphism_matrix)(self, matrix)

    def _test_one(self, **options):
        # Disabled, because we have no proper identity morphism yet
        pass

    def _test_prod(self, **options):
        # Disabled, because we have no proper identity morphism yet
        pass

    def __reduce__(self):
        r"""
        Return a serializable representation of this space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: loads(dumps(A)) == A
            True

        """
        return AffineAutomorphismGroup, (self._surface,)


class AffineAutomorphism_matrix(SurfaceMorphism_factorization):
    r"""
    An affine automorphism of a surface that is an (arbitrary) lift of an
    element of `GL_2(\mathbb{R})`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: f = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
        sage: f
        Affine endomorphism of Translation Surface in H_1(0) built from a square
          Defn: Lift of linear action given by
                [1 3]
                [0 1]

    TESTS::

        sage: from flatsurf.geometry.veech_group import AffineAutomorphism_matrix
        sage: isinstance(f, AffineAutomorphism_matrix)
        True

        sage: TestSuite(f).run()  # optional: pyflatsurf

    """

    def __init__(self, parent, matrix):
        super().__init__(parent)

        if matrix.parent().base_ring() is not parent.base_ring():
            raise TypeError("matrix must be over the surface base ring")

        if matrix.is_mutable():
            raise TypeError("matrix must be immutable")

        self._matrix = matrix

    def _repr_type(self):
        r"""
        Helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 3]
                    [0 1]

        """
        return "Affine"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 3]
                    [0 1]

        """
        return f"Lift of linear action given by\n{self._matrix!r}"

    @cached_method
    def _factorization(self):
        r"""
        Implements :meth:`SurfaceMorphism_factorization._factorization`, the
        actual morphism that does all the heavy lifting for us.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: f = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            sage: f._factorization()  # optional: pyflatsurf
            Composite endomorphism of Translation Surface in H_1(0) built from a square
              Defn:   Linear morphism:
                      From: Translation Surface in H_1(0) built from a square
                      To:   Translation Surface in H_1(0) built from a quadrilateral
                      Defn: [1 3]
                            [0 1]
                    then
                      Delaunay decomposition morphism:
                      From: Translation Surface in H_1(0) built from a quadrilateral
                      To:   Delaunay cell decomposition of Translation Surface in H_1(0) built from a square
                    then
                      Section morphism:
                      From: Delaunay cell decomposition of Translation Surface in H_1(0) built from a square
                      To:   Translation Surface in H_1(0) built from a square
                      Defn: Section of Delaunay decomposition morphism:
                              From: Translation Surface in H_1(0) built from a square
                              To:   Delaunay cell decomposition of Translation Surface in H_1(0) built from a square

        """
        deformation = self.domain().apply_matrix(self._matrix, in_place=False)
        codomain_normalization = self.codomain().delaunay_decompose()
        normalization = deformation.codomain().delaunay_decompose(
            codomain=codomain_normalization.codomain()
        )

        return codomain_normalization.section() * normalization * deformation

    def _test_section_point(self, **options):
        # Disabled, because pyflatsurf backed surfaces cannot have points yet
        pass

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: f = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            sage: g = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            sage: f == g
            True

        """
        if not isinstance(other, AffineAutomorphism_matrix):
            return False

        return self.parent() == other.parent() and self._matrix == other._matrix

    def __hash__(self):
        r"""
        Return a hash value for this automorphism that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: f = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            sage: g = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            sage: hash(f) == hash(g)
            True

        """
        return hash((self.parent(), self._matrix))

    def __reduce__(self):
        r"""
        Return a serializable representation of this automorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: f = A.derivative().section()(matrix([[1, 3], [0, 1]]), check=False)
            sage: loads(dumps(f)) == f
            True

        """
        return AffineAutomorphism_matrix, (self.parent(), self._matrix)


class DerivativeMap(Morphism):
    r"""
    The derivative from the :class:`AffineAutomorphismGroup` to
    `GL_2(\mathbb{R})`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: d = A.derivative()

    TESTS::

        sage: from flatsurf.geometry.veech_group import DerivativeMap
        sage: isinstance(d, DerivativeMap)
        True

        sage: TestSuite(d).run()

    """

    def section(self):
        r"""
        Return a section of this map, i.e., a map that lifts a matrix in
        `GL_2(\mathbb{R})` to an affine automorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: d = A.derivative()
            sage: s = d.section()

        """
        return (
            self.parent().reversed().__make_element_class__(SectionDerivativeMap)(self)
        )

    def _repr_type(self):
        r"""
        Helper for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative()
            Derivative morphism:
              From: ...
              To:   ...

        """
        return "Derivative"

    def __eq__(self, other):
        r"""
        Return whether this derivative is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative() == A.derivative()
            True

        """
        if not isinstance(other, DerivativeMap):
            return False

        return self.parent() == other.parent()

    def __hash__(self):
        r"""
        Return a hash value for this derivative that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: hash(A.derivative()) == hash(A.derivative())
            True

        """
        return hash(self.parent())


class SectionDerivativeMap(Morphism):
    r"""
    A section of a :class:`DerivativeMap`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: d = A.derivative()
        sage: s = d.section()

    TESTS::

        sage: from flatsurf.geometry.veech_group import SectionDerivativeMap
        sage: isinstance(s, SectionDerivativeMap)
        True

        sage: TestSuite(s).run()

    """

    def __init__(self, derivative):
        super().__init__(derivative.parent().reversed())

        self._derivative = derivative

    def section(self):
        r"""
        Return a section of this section, i.e., the original map.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: d = A.derivative()
            sage: s = d.section()

            sage: s.section() is d
            True

        """
        return self._derivative

    def _repr_type(self):
        r"""
        Helper for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: d = A.derivative()
            sage: s = d.section()
            sage: s
            Section morphism:
              From: General Linear Group of degree 2 over Rational Field
              To:   AffineAutomorphismGroup(Translation Surface in H_1(0) built from a square)

        """
        return "Section"

    def _call_(self, matrix):
        r"""
        Return an affine automorphism whose derivative is ``matrix``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: s = A.derivative().section()

        Currently, we cannot check whether an element is in the
        :class:`VeechGroup` so this fails::

            sage: s(matrix([[1, 2], [0, 1]]))
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot decide whether a matrix is in the Veech group yet

        We need to pass ``check=False`` to ignore this check::

            sage: s(matrix([[1, 2], [0, 1]]), check=False)
            Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 2]
                    [0 1]

        """
        return self(matrix, check=True)

    def _call_with_args(self, matrix, args, kwargs):
        r"""
        Helper method for :meth:`_call_` to accept a ``check`` keyword.
        """
        if args:
            raise ValueError("unsupported positional argument")

        check = kwargs.pop("check", True)

        if kwargs:
            raise ValueError("unsupported keyword argument")

        if check:
            if matrix not in VeechGroup(self.codomain().surface()):
                raise ValueError("matrix must be in the Veech group")

        return self.codomain()._from_matrix(matrix.matrix())

    def __eq__(self, other):
        r"""
        Return whether this section is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: A.derivative().section() == A.derivative().section()
            True

        """
        if not isinstance(other, SectionDerivativeMap):
            return False

        return self._derivative == other._derivative

    def __hash__(self):
        r"""
        Return a hash value for this section that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: hash(A.derivative().section()) == hash(A.derivative().section())
            True

        """
        return hash(self._derivative)


def AffineAutomorphismGroup(surface):
    r"""
    Return the group of affine automorphisms of this surface, i.e., the
    group of homeomorphisms that can be locally expressed as affine
    transformations.

    EXAMPLES::

        sage: from flatsurf import dilation_surfaces, AffineAutomorphismGroup
        sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
        sage: AffineAutomorphismGroup(S)
        AffineAutomorphismGroup(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

    .. SEEALSO::

        :meth:`flatsurf.geometry.categories.dilation_surfaces.DilationSurfaces.ParentMethods.affine_automorphism_group`

    """
    return surface.affine_automorphism_group()


def VeechGroup(surface):
    r"""
    Return the Veech group of ``surface``, i.e., the group of matrices that fix
    its vertices.

    EXAMPLES::

        sage: from flatsurf import dilation_surfaces, VeechGroup
        sage: S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
        sage: VeechGroup(S)
        VeechGroup(Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon)

    .. SEEALSO::

        :meth:`flatsurf.geometry.categories.dilation_surfaces.DilationSurfaces.ParentMethods.veech_group`

    """
    return surface.veech_group()
