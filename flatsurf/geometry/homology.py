r"""
Absolute and relative (simplicial) homology of surfaces.

EXAMPLES:

The absolute homology of the regular octagon::

    sage: from flatsurf import translation_surfaces, SimplicialHomology
    sage: S = translation_surfaces.regular_octagon()
    sage: H = SimplicialHomology(S)

A basis of homology, with generators written as (sums of) oriented edges::

    sage: H.gens()
    (B[(0, 1)], B[(0, 2)], B[(0, 3)], B[(0, 0)])

The absolute homology of the unfolding of the (3, 4, 13) triangle::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(angles=[3, 4, 13])
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")
    sage: S.genus()
    8
    sage: H = SimplicialHomology(S)
    sage: len(H.gens())
    16

Relative homology, relative to the singularities of the surface::

    sage: S = S.erase_marked_points()  # optional: pyflatsurf  # random output due to deprecation warnings
    sage: H1 = SimplicialHomology(S, relative=S.vertices())  # optional: pyflatsurf
    sage: len(H1.gens())  # optional: pyflatsurf
    17

We can also form relative `H_0` and `H_2`, though they are not overly
interesting of course::

    sage: H0 = SimplicialHomology(S, relative=S.vertices(), k=0)
    sage: len(H0.gens())
    0

    sage: H2 = SimplicialHomology(S, relative=S.vertices(), k=2)
    sage: len(H2.gens())
    1

We create the homology class corresponding to the core curve of a cylinder (the
interface here is terrible at the moment, see
https://github.com/flatsurf/sage-flatsurf/issues/166)::

    sage: from flatsurf import Polygon, similarity_surfaces, SimplicialHomology, GL2ROrbitClosure

    sage: P = Polygon(angles=[3, 4, 13])
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation").triangulate().codomain()

    sage: from flatsurf.geometry.pyflatsurf.conversion import FlatTriangulationConversion  # optional: pyflatsurf
    sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)  # optional: pyflatsurf
    sage: T = conversion.codomain()  # optional: pyflatsurf
    sage: O = GL2ROrbitClosure(T)  # optional: pyflatsurf

    sage: D = O.decomposition((13, 37))  # optional: pyflatsurf
    sage: cylinder = D.cylinders()[0]  # optional: pyflatsurf

    sage: H = SimplicialHomology(S)
    sage: core = sum(int(str(chain[edge])) * H(conversion.section(edge.positive())) for segment in cylinder.right() for chain in [segment.saddleConnection().chain()] for edge in T.edges())  # optional: pyflatsurf
    sage: core  # optional: pyflatsurf  # random output, the chosen generators vary between operating systems
    972725347814111665129717*B[((0, -1/2*c0, -1/2*c0^2 + 3/2), 2)] + 587352809047576581321682*B[((0, -1/2*c0^2 + 1, -1/2*c0^3 + 3/2*c0), 2)] + 60771110563809382932401*B[((0, -1/2*c0^2 + 1, 1/2*c0^3 - 3/2*c0), 2)] ...

"""

######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2024 Julian Rüth
#                           2023 Julien Boulanger
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
from typing import List, Tuple

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.morphism import Morphism
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.morphism import MorphismSpace


class SimplicialHomologyClass(Element):
    r"""
    An element of a homology group.

    INPUT:

    - ``parent`` -- a :class:`SimplicialHomology`

    - ``chain`` -- an element of the :meth:`SimplicialHomologyGroup.chain_module`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: S = translation_surfaces.regular_octagon()
        sage: H0 = SimplicialHomology(S, k=0)
        sage: g0 = H0.gens()[0]
        sage: g0
        B[Vertex 0 of polygon 0]

        sage: H1 = SimplicialHomology(S, k=1)
        sage: g1 = H1.gens()[0]
        sage: g1
        B[(0, 1)]

        sage: H2 = SimplicialHomology(S, k=2)
        sage: g2 = H2.gens()[0]
        sage: g2
        B[0]

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyClass
        sage: isinstance(g0, SimplicialHomologyClass)
        True

        sage: isinstance(g1, SimplicialHomologyClass)
        True

        sage: isinstance(g2, SimplicialHomologyClass)
        True

    """

    def __init__(self, parent, chain):
        super().__init__(parent)

        self._chain = chain

    def algebraic_intersection(self, other):
        r"""
        Return the algebraic intersection of this class of a closed curve with
        ``other``.

        INPUT:

        - ``other`` - a :class:`SimplicialHomologyClass` in the same homology

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()
            sage: H = SimplicialHomology(S)

            sage: H((0, 0)).algebraic_intersection(H((0, 1)))
            1

            sage: a = H((0, 1))
            sage: b = 3 * H((0, 0)) + 5 * H((0, 2)) - 2 * H((0, 4))
            sage: a.algebraic_intersection(b)
            0

            sage: a = 2 * H((0, 0)) + H((0, 1)) + 3 * H((0, 2)) + H((0, 3)) + H((0, 4)) + H((0, 5)) + H((0, 7))
            sage: b = H((0, 0)) + 2 * H((0, 1)) + H((0, 2)) + H((0, 3)) + 2 * H((0, 4)) + 3 * H((0, 5)) + 4 * H((0, 6)) + 3 * H((0, 7))
            sage: a.algebraic_intersection(b)
            -6

            sage: S = translation_surfaces.cathedral(1, 4)
            sage: H = SimplicialHomology(S)
            sage: a = H((0, 3))
            sage: b = H((2, 1))
            sage: a.algebraic_intersection(b)
            0

            sage: a = H((0, 3))
            sage: b = H((3, 4)) + 3 * H((0, 3)) + 2 * H((0, 0)) - H((1, 7)) + 7 * H((2, 1)) - 2 * H((2, 2))
            sage: a.algebraic_intersection(b)
            2

        """
        if not self.parent().is_absolute():
            raise NotImplementedError(
                "algebraic intersection only available for absolute homology classes"
            )

        other = self.parent()(other)

        if self.parent().degree() != 1:
            raise NotImplementedError(
                "algebraic intersections only available for homology in degree 1"
            )

        intersection = 0

        multiplicities = dict(self._chain)
        other_multiplicities = dict(other._chain)

        for vertex in self.parent().surface().vertices():
            counter = 0
            other_counter = 0

            for edge in [edge for (edge, _) in vertex.edges_ccw()[::2]]:
                opposite_edge = self.parent().surface().opposite_edge(*edge)

                counter += multiplicities.get(edge, 0)
                intersection += counter * other_multiplicities.get(edge, 0)
                intersection -= counter * other_multiplicities.get(opposite_edge, 0)

                counter -= multiplicities.get(opposite_edge, 0)
                other_counter += other_multiplicities.get(edge, 0)
                other_counter -= other_multiplicities.get(opposite_edge, 0)

            if counter:
                raise TypeError("homology class does not correspond to a closed curve")
            if other_counter:
                raise ValueError("homology class does not correspond to a closed curve")

        return intersection

    def _acted_upon_(self, c, self_on_left=None):
        r"""
        Return this homology class scaled by ``c``.

        INPUT:

        - ``c`` -- an element of the base ring of scalars

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: 3 * H.gens()[0]
            3*B[(0, 1)]
            sage: H.gens()[0] * 0
            0

        """
        del self_on_left  # parameter intentionally ignored, the side does not matter
        return self.parent()(c * self._chain)

    @cached_method
    def coefficients(self):
        r"""
        Return the coefficients of this element in terms of the generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0].coefficients()
            (1, 0)
            sage: H.gens()[1].coefficients()
            (0, 1)

        """
        _, _, to_homology = self.parent()._homology()
        return tuple(to_homology(self._chain))

    def _richcmp_(self, other, op):
        r"""
        Return how this class compares to ``other`` with respect to the binary
        relation ``op``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0] == H.gens()[0]
            True
            sage: H.gens()[0] == H.gens()[1]
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

            return self.coefficients() == other.coefficients()

        return super()._richcmp_(other, op)

    def __hash__(self):
        r"""
        Return a hash value of this class that is compatible with
        :meth:`_richcmp_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: hash(H.gens()[0]) == hash(H.gens()[0])
            True

        """
        return hash(self.coefficients())

    def _repr_(self):
        r"""
        Return a printable representation of this class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0]
            B[(0, 1)]

        """
        return repr(self._chain)

    def chain(self):
        r"""
        Return a lift of this element to the
        :meth:`SimplicialHomologyGroup.chain_module`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a.chain()
            B[(0, 1)]

        We can use the chain representation to write a homology class as
        simplices, i.e., edges, with multiplicities::

            sage: coeffs = (a - b).chain().monomial_coefficients()
            sage: coeffs  # random output due to random ordering of edges
            {(0, 1): 1, (0, 0): -1}

        From this representation, we determine the holonomy vector that a chain
        encodes on a translation surface::

            sage: sum(c * T.polygon(label).edge(edge) for ((label, edge), c) in coeffs.items())
            (-1, 1)

        """
        return self._chain

    def coefficient(self, gen):
        r"""
        Return the multiplicity of this class at a generator of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a.coefficient(a)
            1
            sage: a.coefficient(b)
            0

        TESTS::

            sage: a.coefficient(a + b)
            Traceback (most recent call last):
            ...
            ValueError: gen must be a generator not B[(0, 0)] + B[(0, 1)]

        """
        coefficients = gen.coefficients()
        indexes = [i for (i, c) in enumerate(coefficients) if c]

        if len(indexes) != 1 or coefficients[indexes[0]] != 1:
            raise ValueError(f"gen must be a generator not {gen}")

        index = indexes[0]

        return self.coefficients()[index]

    def _add_(self, other):
        r"""
        Return the formal sum of this class and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a + b
            B[(0, 0)] + B[(0, 1)]

        """
        return self.parent()(self._chain + other._chain)

    def _sub_(self, other):
        r"""
        Return the formal difference of this class and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a - b
            -B[(0, 0)] + B[(0, 1)]

        """
        return self.parent()(self._chain - other._chain)

    def _neg_(self):
        r"""
        Return the negative of this homology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a + b
            B[(0, 0)] + B[(0, 1)]
            sage: -(a + b)
            -B[(0, 0)] - B[(0, 1)]

        """
        return self.parent()(-self._chain)

    def surface(self):
        r"""
        Return the surface on which this homology class in defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: h = H.gens()[0]
            sage: h.surface()
            Translation Surface in H_1(0) built from a square

        """
        return self.parent().surface()

    def __bool__(self):
        r"""
        Return whether this class is non-trivial.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: h = H.gens()[0]
            sage: bool(h)
            True
            sage: bool(h-h)
            False

        """
        return bool(self._chain)


class SimplicialHomologyGroup(Parent):
    r"""
    The ``k``-th simplicial homology group of the ``surface`` with
    ``coefficients``.

    .. NOTE:

        This method should not be called directly since it leads to problems
        with pickling and uniqueness. Instead use :meth:`SimplicialHomology` or
        :meth:`homology` on a surface.

    INPUT:

    - ``surface`` -- a finite type surface without boundary

    - ``k`` -- an integer

    - ``coefficients`` -- a ring

    - ``relative`` -- a subset of points of the ``surface``

    - ``implementation`` -- a string; the algorithm used to compute the
      homology, only ``"generic"`` is supported at the moment which uses the
      generic homology machinery of SageMath.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology, MutableOrientedSimilaritySurface
        sage: T = translation_surfaces.square_torus()

    Surfaces must be immutable to compute their homology::

        sage: T = MutableOrientedSimilaritySurface.from_surface(T)
        sage: SimplicialHomology(T)
        Traceback (most recent call last):
        ...
        ValueError: surface must be immutable to compute homology

    ::

        sage: T.set_immutable()
        sage: SimplicialHomology(T)
        H₁(Translation Surface in H_1(0) built from a square)

    TESTS::

        sage: T = translation_surfaces.square_torus()
        sage: H = SimplicialHomology(T, implementation="generic")
        sage: TestSuite(H).run()

    Verify that #310 has been resolved::

        sage: S = translation_surfaces.mcmullen_L(1,1,1,1)
        sage: H = S.homology(coefficients=QQ)
        sage: H.gens()
        (B[(0, 1)], B[(0, 0)], B[(1, 1)], B[(2, 0)])

    """

    Element = SimplicialHomologyClass

    def __init__(self, surface, k, coefficients, relative, implementation, category):
        Parent.__init__(self, base=coefficients, category=category)

        if surface.is_mutable():
            raise TypeError("surface must be immutable")

        from sage.all import ZZ

        if k not in ZZ:
            raise TypeError("k must be an integer")

        from sage.categories.all import Rings

        if coefficients not in Rings():  # pyright: ignore[reportCallIssue]
            raise TypeError("coefficients must be a ring")

        if relative:
            for point in relative:
                if point not in surface.vertices():
                    raise NotImplementedError(
                        "can only compute homology relative to a subset of the vertices"
                    )

        if implementation == "generic":
            if not surface.is_finite_type():
                raise NotImplementedError(
                    "homology only implemented for surfaces with finitely many polygons"
                )

            if surface.is_with_boundary():
                raise NotImplementedError(
                    "homology only implemented for surfaces without boundary"
                )
        else:
            raise NotImplementedError(
                "cannot compute homology with this implementation yet"
            )

        self._surface = surface
        self._k = k
        self._coefficients = coefficients
        self._relative = relative
        self._implementation = implementation

    def is_absolute(self):
        r"""
        Return whether this is absolute homology (and not relative to some set
        of points).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.is_absolute()
            True

        """
        return not self._relative

    def some_elements(self):
        r"""
        Return some typical homology classes (for testing).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.some_elements()
            [0, B[(0, 1)], B[(0, 0)]]

        """
        return [self.zero()] + list(self.gens())

    def surface(self):
        r"""
        Return the surface of which this is the homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.surface() == T
            True

        """
        return self._surface

    @cached_method
    def chain_module(self):
        r"""
        Return the free module of simplicial chains of the surface i.e., formal
        sums of simplicies, e.g., formal sums of edges of the surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.chain_module()
            Free module generated by {(0, 1), (0, 0)} over Integer Ring

        """
        from sage.all import FreeModule

        return FreeModule(self._coefficients, self.simplices())

    @cached_method
    def simplices(self):
        r"""
        Return the simplices that form the generators of :meth:`chain_module`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        In dimension 1, this is the set of edges::

            sage: H = SimplicialHomology(T)
            sage: H.simplices()
            ((0, 1), (0, 0))

        In dimension 0, this is the set of vertices::

            sage: H = SimplicialHomology(T, k=0)
            sage: H.simplices()
            (Vertex 0 of polygon 0,)

        In dimension 2, this is the set of polygons::

            sage: H = SimplicialHomology(T, k=2)
            sage: H.simplices()
            (0,)

        In all other dimensions, there are no simplices::

            sage: H = SimplicialHomology(T, k=12)
            sage: H.simplices()
            ()

        """
        if self._k == 0:
            return tuple(
                vertex
                for vertex in self._surface.vertices()
                if vertex not in self._relative
            )
        if self._k == 1:
            simplices = set()
            for edge in self._surface.edges():
                if self._surface.opposite_edge(*edge) not in simplices:
                    simplices.add(edge)
            return tuple(simplices)
        if self._k == 2:
            return tuple(self._surface.labels())

        return ()

    @cached_method
    def change(self, k=None):
        r"""
        Return this homology but in dimension ``k``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H
            H₁(Translation Surface in H_1(0) built from a square)

            sage: H.change(k=0)
            H₀(Translation Surface in H_1(0) built from a square)

        """
        return SimplicialHomology(
            surface=self._surface,
            k=k if k is not None else self._k,
            coefficients=self._coefficients,
            relative=self._relative,
            implementation=self._implementation,
            category=self.category(),
        )

    def boundary(self, chain):
        r"""
        Return the boundary of ``chain`` as an element of the
        :meth:`chain_module` in lower dimension.

        INPUT:

        - ``chain`` -- an element of :meth:`chain_module`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        ::

            sage: H = SimplicialHomology(T, k=0)
            sage: c = H.chain_module().an_element(); c
            2*B[Vertex 0 of polygon 0]
            sage: H.boundary(c)
            0

        ::

            sage: H = SimplicialHomology(T, k=1)
            sage: c = H.chain_module().an_element(); c
            2*B[(0, 0)] + 2*B[(0, 1)]
            sage: H.boundary(c)
            0

        ::

            sage: H = SimplicialHomology(T, k=2)
            sage: c = H.chain_module().an_element(); c
            2*B[0]
            sage: H.boundary(c)
            0

        """
        chain = self.chain_module()(chain)

        if self._k == 1:
            C0 = self.change(k=0).chain_module()

            def to_C0(point):
                if point in self._relative:
                    return C0.zero()
                return C0(point)

            boundary = C0.zero()
            for edge, coefficient in chain:
                boundary += coefficient * to_C0(
                    self._surface.point(*self._surface.opposite_edge(*edge))
                )
                boundary -= coefficient * to_C0(self._surface.point(*edge))
            return boundary

        if self._k == 2:
            C1 = self.change(k=1).chain_module()
            boundary = C1.zero()
            for face, coefficient in chain:
                for edge in range(len(self._surface.polygon(face).edges())):
                    if (face, edge) in C1.indices():
                        boundary += coefficient * C1((face, edge))
                    else:
                        boundary -= coefficient * C1(
                            self._surface.opposite_edge(face, edge)
                        )
            return boundary

        return self.change(k=self._k - 1).chain_module().zero()

    @cached_method
    def _chain_complex(self):
        r"""
        Return the chain complex of vector spaces that is implementing this
        homology (if the ``"generic"`` implementation has been selected).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._chain_complex()
            Chain complex with at most 3 nonzero terms over Integer Ring

        ::

            sage: H = SimplicialHomology(T, relative=T.vertices())
            sage: H._chain_complex()
            Chain complex with at most 2 nonzero terms over Integer Ring

        """

        def boundary(dimension, chain):
            boundary = self.change(k=dimension).boundary(chain)
            coefficients = boundary.dense_coefficient_list(
                self.change(k=dimension - 1).chain_module().indices()
            )
            return coefficients

        from sage.all import ChainComplex, matrix

        return ChainComplex(
            {
                dimension: matrix(
                    [
                        boundary(dimension, simplex)
                        for simplex in self.change(k=dimension).chain_module().basis()
                    ]
                ).transpose()
                for dimension in range(3)
            },
            base_ring=self._coefficients,
            degree=-1,
        )

    def zero(self) -> SimplicialHomologyClass:
        r"""
        Return the zero homology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.zero()
            0

        """
        return self(self.chain_module().zero())

    @cached_method
    def _homology(self):
        r"""
        Return the free module isomorphic to homology, a lift from that
        module to the chain module, and an inverse (modulo boundaries).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._homology()
            (Finitely generated module V/W over Integer Ring with invariants (0, 0),
             Generic morphism:
               From: Finitely generated module V/W over Integer Ring with invariants (0, 0)
               To:   Free module generated by {(0, 1), (0, 0)} over Integer Ring,
             Generic morphism:
               From: Free module generated by {(0, 1), (0, 0)} over Integer Ring
               To:   Finitely generated module V/W over Integer Ring with invariants (0, 0))

        """
        if self._implementation == "generic":
            C = self._chain_complex()

            cycles = C.differential(self._k).transpose().kernel()
            boundaries = C.differential(self._k + 1).transpose().image()
            homology = cycles.quotient(boundaries)

            F = self.chain_module()

            def lift(x):
                if hasattr(homology, "lift"):
                    # Available for quotients of vector spaces
                    return homology.lift(x)

                from sage.all import vector

                # Available on quotients of other modules
                return vector(x.lift().lift())

            from_homology = homology.module_morphism(
                function=lambda x: F.from_vector(lift(x)),
                codomain=F,
            )

            def _to_homology(x):
                multiplicities = x.dense_coefficient_list(order=F.get_order())
                try:
                    cycle = cycles(multiplicities)
                except TypeError:
                    if multiplicities not in cycles:
                        raise ValueError(
                            "chain is not a cycle so it has no representation in homology"
                        )
                    raise

                return homology(cycle)

            to_homology = F.module_morphism(function=_to_homology, codomain=homology)

            for gen in homology.gens():
                assert to_homology(from_homology(gen)) == gen

            return homology, from_homology, to_homology

        raise NotImplementedError(
            "cannot compute homology with this implementation yet"
        )

    def _test_homology(self, **options):
        r"""
        Test that :meth:`_homology` computes homology correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._test_homology()

        """
        tester = self._tester(**options)

        homology, from_homology, to_homology = self._homology()
        chains = self.chain_module()

        tester.assertEqual(homology, to_homology.codomain())
        tester.assertEqual(homology, from_homology.domain())
        tester.assertEqual(chains, to_homology.domain())
        tester.assertEqual(chains, from_homology.codomain())

        for gen in homology.gens():
            tester.assertEqual(to_homology(from_homology(gen)), gen)

    def _repr_(self):
        r"""
        Return a printable representation of this homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H
            H₁(Translation Surface in H_1(0) built from a square)

        """
        k = self._k
        if k == 0:
            k = "₀"
        elif k == 1:
            k = "₁"
        elif k == 2:
            k = "₂"
        else:
            k = f"_{k}"

        H_k = f"H{k}"

        X = repr(self.surface())
        if not self.is_absolute():
            X = f"{X}, {set(self._relative)}"

        from sage.all import ZZ

        if self._coefficients is not ZZ:
            sep = ";"
            X = f"{X}{sep} {self._coefficients}"

        return f"{H_k}({X})"

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as an element of this homology.

        TESTS::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)

        ::

            sage: H(0)
            0
            sage: H(None)
            0

        ::

            sage: H((0, 0))
            B[(0, 0)]
            sage: H((0, 2))
            -B[(0, 0)]

        ::

            sage: H = SimplicialHomology(T, 0)
            sage: H(H.chain_module().gens()[0])
            B[Vertex 0 of polygon 0]

            sage: H = SimplicialHomology(T, 1)
            sage: H(H.chain_module().gens()[0])
            B[(0, 1)]

            sage: H = SimplicialHomology(T, 2)
            sage: H(H.chain_module().gens()[0])
            B[0]

        """
        if x == 0 or x is None:
            return self.element_class(self, self.chain_module().zero())

        if self._k == 1 and isinstance(x, tuple) and len(x) == 2:
            sgn = 1
            if x not in self.simplices():
                x = self.surface().opposite_edge(*x)
                sgn = -1
            assert x in self.simplices()
            return sgn * self.element_class(self, self.chain_module()(x))

        if x.parent() is self.chain_module():
            return self.element_class(self, x)

        try:
            hom_method = x._homology_
        except AttributeError:
            pass
        else:
            return hom_method(self)

        raise NotImplementedError("cannot convert this element to a homology class yet")

    @cached_method
    def gens(self) -> Tuple[SimplicialHomologyClass, ...]:
        r"""
        Return generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        ::

            sage: H = SimplicialHomology(T)
            sage: H.gens()
            (B[(0, 1)], B[(0, 0)])

        ::

            sage: H = SimplicialHomology(T, 0)
            sage: H.gens()
            (B[Vertex 0 of polygon 0],)

        ::

            sage: H = SimplicialHomology(T, 2)
            sage: H.gens()
            (B[0],)

        """
        if self._k < 0 or self._k > 2:
            return ()

        homology, from_homology, _ = self._homology()
        return tuple(self(from_homology(g)) for g in homology.gens())

    def ngens(self):
        r"""
        Return the Betti number of this homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

            sage: H = T.homology()
            sage: H.ngens()
            2

        """
        return len(self.gens())

    def degree(self):
        r"""
        Return the degree `k` for this homology `H_k`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

            sage: H = SimplicialHomology(T)
            sage: H.degree()
            1

        """
        return self._k

    def symplectic_basis(self) -> List[SimplicialHomologyClass]:
        r"""
        Return a symplectic basis of generators of this homology group.

        A symplectic basis is a basis of the form `(e_1, \ldots, e_n, f_1,
        \ldots, f_n)` such that for the
        :meth:`SimplicialHomologyClass.algebraic_intersection`, `e_i
        \cdot e_j = f_i \cdot f_j = 0` and `e_i \cdot f_j = \delta_{ij}` for
        all `i` and `j`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        ::

            sage: H = SimplicialHomology(T)
            sage: H.symplectic_basis()
            [B[(0, 0)], B[(0, 1)]]

        """
        from sage.all import matrix

        E = matrix(
            self.base_ring(),
            [[g.algebraic_intersection(h) for h in self.gens()] for g in self.gens()],
        )

        from sage.categories.all import Fields

        if self.base_ring() in Fields:
            from sage.matrix.symplectic_basis import symplectic_basis_over_field

            F, C = symplectic_basis_over_field(E)
        else:
            from sage.matrix.symplectic_basis import symplectic_basis_over_ZZ

            F, C = symplectic_basis_over_ZZ(E)

        if any(entry not in [-1, 0, 1] for row in F for entry in row):
            raise NotImplementedError(
                "cannot determine symplectic basis for this homology group over this ring yet"
            )

        return [
            sum(c * g for (c, g) in zip(row, self.gens())) for row in C
        ]  # pyright: ignore

    def _test_symplectic_basis(self, **options):
        r"""
        Verify that :meth:`symplectic_basis` has been implemented correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._test_symplectic_basis()

        """
        tester = self._tester(**options)

        basis = self.symplectic_basis()
        n = len(basis)

        tester.assertEqual(len(self.gens()), n)
        tester.assertEqual(n % 2, 0)

        A = basis[: n // 2]
        B = basis[n // 2 :]

        for i, a in enumerate(A):
            for j, b in enumerate(B):
                tester.assertEqual(a.algebraic_intersection(b), i == j)

        for a in A:
            for aa in A:
                tester.assertEqual(a.algebraic_intersection(aa), 0)

        for b in B:
            for bb in B:
                tester.assertEqual(b.algebraic_intersection(bb), 0)

    def hom(self, f, codomain=None):
        r"""
        Return the homomorphism of homology induced by ``f``.

        INPUT:

        - ``f`` -- a morphism of surfaces or a matrix

        - ``codomain`` -- the simplicial homology this morphism maps into

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: g.matrix()  # optional: pyflatsurf
            [1 0]
            [2 1]

            sage: H.gens()
            (B[(0, 1)], B[(0, 0)])
            sage: [g(h) for h in H.gens()]  # optional: pyflatsurf
            [2*B[(0, 0)] + B[(0, 1)], B[(0, 0)]]

        """
        from flatsurf.geometry.veech_group import SurfaceMorphism
        from sage.matrix.matrix0 import Matrix
        from sage.all import Hom

        if isinstance(f, SurfaceMorphism):
            if codomain is None:
                codomain = f.codomain().homology()

            if codomain.surface() is not f.codomain():
                raise ValueError("codomain must be codomain of morphism or None")

            if f.domain() is self.surface():
                parent = Hom(self, codomain)
                return parent.__make_element_class__(
                    SimplicialHomologyMorphism_induced
                )(parent, f)
        elif isinstance(f, Matrix):
            if codomain is None:
                if f.is_square():
                    codomain = self

            if codomain is None:
                raise NotImplementedError("cannot deduce codomain from this matrix")

            if f.ncols() != self.ngens():
                raise ValueError(
                    "matrix must have one column for each generator of homology"
                )

            if f.nrows() != codomain.ngens():
                raise ValueError(
                    "matrix must have one row for each generator of the codomain"
                )

            parent = Hom(self, codomain)
            return parent.__make_element_class__(SimplicialHomologyMorphism_matrix)(
                parent, f
            )

        raise NotImplementedError(
            "cannot create a morphism in homology from this data yet"
        )

    def _Hom_(self, Y, category=None):
        r"""
                Return the space of morphisms from this homology to ``Y``.
        k        EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: H = S.homology()

                    sage: End(H)
                    Endomorphisms of H₁(Translation Surface in H_1(0) built from a square)

        """
        if isinstance(Y, SimplicialHomologyGroup):
            return SimplicialHomologyMorphismSpace(self, Y, category=category)

        return super()._Hom_(Y, category=category)

    def __eq__(self, other):
        r"""
        Return whether this homology is indistinguishable from ``other``.

        .. NOTE::

            We cannot rely on the builtin `==` by ``id`` since we need to
            detect homologies over equal but distinct surfaces to be equal. See
            :meth:`homology` for ideas on how to fix this.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)

            sage: T = translation_surfaces.square_torus()
            sage: HH = SimplicialHomology(T)

            sage: H == HH
            True

        """
        if not isinstance(other, SimplicialHomologyGroup):
            return False

        return (
            self._surface == other._surface
            and self._coefficients == other._coefficients
            and self._relative == other._relative
            and self._implementation == other._implementation
            and self.category() == other.category()
        )

    def __hash__(self):
        r"""
        Return a hash value for this homology that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)

            sage: T = translation_surfaces.square_torus()
            sage: HH = SimplicialHomology(T)

            sage: hash(H) == hash(HH)
            True

        """
        return hash(
            (
                self._surface,
                self._coefficients,
                self._relative,
                self._implementation,
                self.category(),
            )
        )


class SimplicialHomologyMorphismSpace(MorphismSpace):
    r"""
    The space of homomorphisms from the homology ``domain`` to ``codomain``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: M = matrix([[1, 0], [0, 1]])
        sage: f = A.derivative().section()(M, check=False)

        sage: from flatsurf import SimplicialHomology
        sage: H = SimplicialHomology(S)
        sage: g = H.hom(f)
        sage: G = g.parent()

    Since these are homomorphisms in homology, they preserve the linear
    structure of the homology::

        sage: g.category()
        Category of endsets of modules over Integer Ring

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphismSpace
        sage: isinstance(G, SimplicialHomologyMorphismSpace)
        True

        sage: TestSuite(G).run()

    """

    def __init__(self, domain, codomain, category=None):
        from sage.all import Hom

        super().__init__(
            domain,
            codomain,
            category=category or Hom(domain, codomain).homset_category(),
        )

    def an_element(self):
        r"""
        Return some homomorphism in homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S.homology()).an_element()
            Generic endomorphism of H₁(Translation Surface in H_1(0) built from a square)
              Defn: [1 0]
                    [0 1]

            sage: Hom(S.homology(), T.homology()).an_element()
            Generic morphism:
              From: H₁(Translation Surface in H_1(0) built from a square)
              To:   H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [0 0]
                    [0 0]
                    [0 0]
                    [0 0]

        """
        if self.domain() is self.codomain():
            return self.identity()
        return self.zero()

    def identity(self):
        r"""
        Return the identity homomorphism in this space (if it exists).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S.homology()).identity()
            Generic endomorphism of H₁(Translation Surface in H_1(0) built from a square)
              Defn: [1 0]
                    [0 1]

        """
        if self.is_endomorphism_set():
            from sage.all import identity_matrix

            matrix = identity_matrix(
                self.codomain().base_ring(), self.domain().ngens(), sparse=True
            )
            return self.__make_element_class__(SimplicialHomologyMorphism_matrix)(
                self, matrix
            )
        return super().identity()

    def zero(self):
        r"""
        Return the zero homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: End(S.homology()).zero()
            Generic endomorphism of H₁(Translation Surface in H_1(0) built from a square)
              Defn: [0 0]
                    [0 0]

        """
        from sage.all import zero_matrix

        return self.domain().hom(
            zero_matrix(
                self.codomain().base_ring(),
                nrows=self.codomain().ngens(),
                ncols=self.domain().ngens(),
                sparse=True,
            ),
            codomain=self.codomain(),
        )

    def base_ring(self):
        r"""
        Return the ring over which these homomorphisms are defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: End(S.homology()).base_ring()
            Integer Ring

        """
        if self.domain().base_ring() is self.codomain().base_ring():
            return self.domain().base_ring()

        return super().base_ring()

    def __reduce__(self):
        r"""
        Return a picklable version of this space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: H = End(S.homology()).base_ring()
            sage: loads(dumps(H)) == H
            True

        """
        return SimplicialHomologyMorphismSpace, (
            self.domain(),
            self.codomain(),
            self.homset_category(),
        )

    def __repr__(self):
        r"""
        Return a printable representation of this space of homomorphisms.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: H = End(S.homology())
            sage: H
            Endomorphisms of H₁(Translation Surface in H_1(0) built from a square)

        """
        if self.domain() == self.codomain():
            return f"Endomorphisms of {self.domain()!r}"
        return f"Homomorphisms from {self.domain()!r} to {self.codomain()!r}"


class SimplicialHomologyMorphism_base(Morphism):
    r"""
    Base class for all homomorphisms in homology.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: M = matrix([[1, 0], [0, 1]])
        sage: f = A.derivative().section()(M, check=False)

        sage: from flatsurf import SimplicialHomology
        sage: H = SimplicialHomology(S)
        sage: g = H.hom(f)

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphism_base
        sage: isinstance(g, SimplicialHomologyMorphism_base)
        True

        sage: TestSuite(g).run()  # optional: pyflatsurf

    """

    @cached_method
    def matrix(self):
        r"""
        Return the matrix describing this homomorphism on the generators of
        homology (as a multiplication from the left).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: H.gens()
            (B[(0, 1)], B[(0, 0)], B[(1, 1)], B[(2, 0)])
            sage: g = H.hom(f)

            sage: g.matrix()  # optional: pyflatsurf
            [1 0 0 0]
            [1 1 2 0]
            [0 0 1 0]
            [1 0 0 1]

        """
        from sage.all import matrix

        return matrix(
            [list(self(gen).coefficients()) for gen in self.domain().gens()]
        ).transpose()

    def _add_(self, other):
        r"""
        Return the pointwise sum of this morphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = End(H).one()

            sage: g + h  # optional: pyflatsurf
            Generic endomorphism of H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [2 0 0 0]
                    [1 2 2 0]
                    [0 0 2 0]
                    [1 0 0 2]

        """
        return self.domain().hom(
            self.matrix() + other.matrix(), codomain=self.codomain()
        )

    def _acted_upon_(self, x, self_on_left):
        r"""
        Return the morphism given by pointwise multiplying with ``x``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: (2**1234567 * g).matrix().trace() == 2**1234569  # optional: pyflatsurf
            True

        """
        return self.domain().hom(x * self.matrix())

    def _neg_(self):
        r"""
        Return the pointwise negative of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: -g  # optional: pyflatsurf
            Generic endomorphism of H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [-1  0  0  0]
                    [-1 -1 -2  0]
                    [ 0  0 -1  0]
                    [-1  0  0 -1]

        """
        return self.domain().hom(-self.matrix(), codomain=self.codomain())

    def _composition(self, other):
        r"""
        Return the composition of this homomorphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 2)
            sage: U = translation_surfaces.mcmullen_L(1, 1, 1, 3)

            sage: f = S.homology().hom(2 * identity_matrix(4), codomain=T.homology())
            sage: g = T.homology().hom(3 * identity_matrix(4), codomain=U.homology())

            sage: g * f
            Generic morphism:
              From: H₁(Translation Surface in H_2(2) built from 3 squares)
              To:   H₁(Translation Surface in H_2(2) built from 2 squares and a rectangle)
              Defn: [6 0 0 0]
                    [0 6 0 0]
                    [0 0 6 0]
                    [0 0 0 6]

        """
        return other.domain().hom(
            self.matrix() * other.matrix(), codomain=self.codomain()
        )

    def __bool__(self):
        r"""
        Return whether this is not the homommorphism that is zero everywhere.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: bool(g)  # optional: pyflatsurf
            True

            sage: bool(g.parent().zero())
            False

        """
        return bool(self.matrix())


class SimplicialHomologyMorphism_matrix(SimplicialHomologyMorphism_base):
    r"""
    A homomorphism of homology that is given by a matrix that describes the
    homomorphism on the generators.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
        sage: T = translation_surfaces.square_torus()

        sage: f = S.homology().hom(matrix([[1, 2, 3, 4], [5, 6, 7, 8]]), codomain=T.homology())
        sage: f
        Generic morphism:
          From: H₁(Translation Surface in H_2(2) built from 3 squares)
          To:   H₁(Translation Surface in H_1(0) built from a square)
          Defn: [1 2 3 4]
                [5 6 7 8]

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphism_matrix
        sage: isinstance(f, SimplicialHomologyMorphism_matrix)
        True

        sage: TestSuite(f).run()

    """

    def __init__(self, parent, matrix):
        super().__init__(parent)

        if matrix.is_mutable():
            from sage.all import matrix as copy

            matrix = copy(matrix)
            matrix.set_immutable()

        self._matrix = matrix

    def _call_(self, g):
        r"""
        Return the image of ``g`` under this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = translation_surfaces.square_torus()

            sage: f = S.homology().hom(matrix([[1, 2, 3, 4], [5, 6, 7, 8]]), codomain=T.homology())
            sage: [f(gen) for gen in S.homology().gens()]
            [5*B[(0, 0)] + B[(0, 1)],
             6*B[(0, 0)] + 2*B[(0, 1)],
             7*B[(0, 0)] + 3*B[(0, 1)],
             8*B[(0, 0)] + 4*B[(0, 1)]]

        """
        from sage.all import vector

        image = self._matrix * vector(g.coefficients())

        homology, to_chain, to_homology = self.codomain()._homology()

        image = sum(
            coefficient * gen for (coefficient, gen) in zip(image, homology.gens())
        )

        return self.codomain()(to_chain(image))

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        .. NOTE::

            We cannot override ``_richcmp_`` since our non-uniqueness of
            surfaces breaks the coercion framework in SageMath.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = H.hom(g.matrix())  # optional: pyflatsurf
            sage: h == h  # optional: pyflatsurf
            True

        Note that this determines whether two morphisms are indistinguishable,
        not whether they are pointwise the same::

            sage: h == g  # optional: pyflatsurf
            False

        """
        if not isinstance(other, SimplicialHomologyMorphism_matrix):
            return False

        return self.parent() == other.parent() and self._matrix == other._matrix

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = End(S.homology()).one()
            sage: g = End(S.homology()).one()

            sage: hash(f) == hash(g)
            True

        """
        return hash((self.parent(), self._matrix))

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = End(S.homology()).one()
            sage: f
            Generic endomorphism of H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 1]

        """
        return repr(self._matrix)


class SimplicialHomologyMorphism_induced(SimplicialHomologyMorphism_base):
    r"""
    A homomorphism of homology induced by a morphism of surfaces.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: M = matrix([[1, 2], [0, 1]])
        sage: f = A.derivative().section()(M, check=False)

        sage: from flatsurf import SimplicialHomology
        sage: H = SimplicialHomology(S)
        sage: g = H.hom(f)

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphism_induced
        sage: isinstance(g, SimplicialHomologyMorphism_induced)
        True

        sage: TestSuite(g).run()  # optional: pyflatsurf

    """

    def __init__(self, parent, morphism):
        super().__init__(parent)

        self._morphism = morphism

    def _call_(self, x):
        r"""
        Return the image of the homology class ``x`` under this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: H.gens()
            (B[(0, 1)], B[(0, 0)])
            sage: [g(h) for h in H.gens()]  # optional: pyflatsurf
            [2*B[(0, 0)] + B[(0, 1)], B[(0, 0)]]

        """
        return self._morphism._image_homology(x, codomain=self.codomain())

    def _repr_type(self):
        r"""
        Helper method for :meth:`_repr_` to produce a printable representation
        of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: g._repr_type()
            'Induced'

        """
        return "Induced"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_` to produce a printable representation
        of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: print(g._repr_defn())
            Induced by Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 2]
                    [0 1]

        """
        return f"Induced by {self._morphism!r}"

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        .. NOTE::

            We cannot override ``_richcmp_`` since our non-uniqueness of
            surfaces breaks the coercion framework in SageMath.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = H.hom(f)

            sage: g == h
            True

        Note that this does not compare homomorphisms pointwise::

            sage: h = H.hom(g.matrix())  # optional: pyflatsurf
            sage: g == h  # optional: pyflatsurf
            False

        """
        if not isinstance(other, SimplicialHomologyMorphism_induced):
            return False

        return self.parent() == other.parent() and self._morphism == other._morphism

    def __hash__(self):
        r"""
        Return a hash value for this homomorphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = H.hom(f)
            sage: hash(g) == hash(h)
            True

        """
        return hash((self.parent(), self._morphism))


def SimplicialHomology(
    surface,
    k=1,
    coefficients=None,
    relative=None,
    implementation="generic",
    category=None,
):
    r"""
    Return the ``k``-th simplicial homology group of ``surface``.

    INPUT:

    - ``surface`` -- a surface

    - ``k`` -- an integer (default: ``1``)

    - ``coefficients`` -- a ring (default: the integer ring);
      consider the homology with coefficients in this ring

    - ``relative`` -- a set (default: the empty set); if non-empty,
      then relative homology with respect to this set is
      constructed.

    - ``implementation`` -- a string (default: ``"generic"``); the
      algorithm used to compute the homology groups. Currently only
      ``"generic"`` is supported, i.e., the groups are computed
      with the generic homology machinery from SageMath.

    - ``category`` -- a category; if not specified, a category for
      the homology group is chosen automatically depending on
      ``coefficients``.

    TESTS:

    Homology is unique and cached::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: T = translation_surfaces.square_torus()
        sage: SimplicialHomology(T) is SimplicialHomology(T)
        True

    """
    return surface.homology(k, coefficients, relative, implementation, category)
