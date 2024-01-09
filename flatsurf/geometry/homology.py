r"""
Absolute and relative (simplicial) homology of surfaces.

EXAMPLES:

The absolute homology of the regular octagon::

    sage: from flatsurf import translation_surfaces, SimplicialHomology
    sage: S = translation_surfaces.regular_octagon()
    sage: H = SimplicialHomology(S)

A basis of homology, with generators written as oriented edges::

    sage: H.gens()
    (B[(0, 1)], B[(0, 2)], B[(0, 3)], B[(0, 0)])

The absolute homology of the unfolding of the (3, 4, 13) triangle::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(angles=[3, 4, 13])
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")
    sage: H = SimplicialHomology(S)
    sage: len(H.gens())
    16

Relative homology, relative to the singularities of the surface::

    sage: S = S.erase_marked_points().codomain()
    sage: H1 = SimplicialHomology(S, relative=S.vertices())
    sage: len(H1.gens())
    17
    sage: H0 = SimplicialHomology(S, relative=S.vertices(), k=0)
    sage: len(H0.gens())
    0
    sage: H2 = SimplicialHomology(S, relative=S.vertices(), k=2)
    sage: len(H2.gens())
    1

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

from sage.structure.parent import Parent
from sage.structure.element import Element

from sage.misc.cachefunc import cached_method


class SimplicialHomologyClass(Element):
    def __init__(self, parent, chain):
        super().__init__(parent)

        self._chain = chain

    def algebraic_intersection(self, other):
        r"""
        Return the algebraic intersection of this class of a closed curve with
        ``other``.

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
        intersection = 0

        multiplicities = dict(self._chain)
        other_multiplicities = dict(other._chain)

        for _, adjacent_edges in self.parent().surface().angles(return_adjacent_edges=True):
            counter = 0
            other_counter = 0
            for edge in adjacent_edges:
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
        Return the coefficients of this element in terms of the generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: 3 * H.gens()[0]
            3*B[(0, 1)]

        """
        return self.parent()(c * self._chain)

    @cached_method
    def _homology(self):
        r"""
        Return the coefficients of this element in terms of the generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0]._homology()
            (1, 0)
            sage: H.gens()[1]._homology()
            (0, 1)

        """
        _, _, to_homology = self.parent()._homology()
        return tuple(to_homology(self._chain))

    def _richcmp_(self, other, op):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
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

            return self._homology() == other._homology()

        return super()._richcmp_(other, op)

    def __hash__(self):
        return hash(self._homology())

    def _repr_(self):
        return repr(self._chain)

    def coefficient(self, gen):
        r"""
        Return the multiplicity of this class at a generator of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: a,b = H.gens()
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
        coefficients = gen._homology()
        indexes = [i for (i, c) in enumerate(coefficients) if c]

        if len(indexes) != 1 or coefficients[indexes[0]] != 1:
            raise ValueError(f"gen must be a generator not {gen}")

        index = indexes[0]

        return self._homology()[index]

    def _add_(self, other):
        r"""
        Return the formal sum of homology classes.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a + b
            B[(0, 0)] + B[(0, 1)]

        """
        return self.parent()(self._chain + other._chain)

    def _sub_(self, other):
        return self.parent()(self._chain - other._chain)

    def _neg_(self):
        r"""
        Return the negative of this homology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a + b
            B[(0, 0)] + B[(0, 1)]
            sage: -(a + b)
            -B[(0, 0)] - B[(0, 1)]

        """
        return self.parent()(-self._chain)

    def surface(self):
        return self.parent().surface()

    def __bool__(self):
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

    - ``surface`` -- a finite :class:`flatsurf.geometry.surface.Surface`
      without boundary

    - ``k`` -- an integer

    - ``coefficients`` -- a ring (default: the integers)

    - ``generators`` -- one of ``edge``, ``interior``, ``midpoint`` (default:
      ``edge``) how generators are represented

    - ``relative`` -- a subset of points of the ``surface`` (default: the empty
      set)

    - ``implementation`` -- one of ``"spanning_tree"`` or ``"generic"``
      (default: ``"generic"``); whether the homology is computed with a (very
      efficient) spanning tree algorithm or the generic homology machinery
      provided by SageMath.


    .. TODO::

        Implement the ``"spanning_tree"`` algorithm froom
        ``GL2ROrbitClosure._spanning_tree``.

    .. TODO::

        Use https://github.com/flatsurf/sage-flatsurf/pull/114/files to force
        the representatives to live in particular subgraph of the dual graph.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology, MutableOrientedSimilaritySurface
        sage: T = translation_surfaces.torus((1, 0), (0, 1))

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

        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: H = SimplicialHomology(T, implementation="spanning_tree")  # not tested, spanning_tree not implemented yet
        sage: TestSuite(H).run()  # not tested, spanning_tree not implemented yet

    ::

        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: H = SimplicialHomology(T, implementation="generic")
        sage: TestSuite(H).run()

    """
    Element = SimplicialHomologyClass

    def __init__(self, surface, k, coefficients, generators, relative, implementation, category):
        Parent.__init__(self, base=coefficients, category=category)

        if surface.is_mutable():
            raise TypeError("surface must be immutable")

        from sage.all import ZZ
        if k not in ZZ:
            raise TypeError("k must be an integer")

        from sage.categories.all import Rings
        if coefficients not in Rings():
            raise TypeError("coefficients must be a ring")

        if generators not in ["edge",]:
            raise NotImplementedError("cannot represented homology with these generators yet")

        if relative:
            for point in relative:
                if point not in surface.vertices():
                    raise NotImplementedError("can only compute homology relative to a subset of the vertices")

        if implementation not in ["generic"]:
            raise NotImplementedError("cannot compute homology with this implementation yet")

        self._surface = surface
        self._k = k
        self._coefficients = coefficients
        self._generators = generators
        self._relative = relative
        self._implementation = implementation

    def surface(self):
        r"""
        Return the surface of which this is the homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.surface() == T
            True

        """
        return self._surface

    @cached_method
    def chain_module(self):
        r"""
        Return the free module of simplicial chains of the
        triangulation, i.e., formal sums of simplicies, e.g., formal sums of
        edges of the triangulation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.simplices()
            ((0, 1), (0, 0))

        """
        if self._k == 0:
            return tuple(vertex for vertex in self._surface.vertices() if vertex not in self._relative)
        if self._k == 1:
            simplices = set()
            for edge in self._surface.edges():
                if self._surface.opposite_edge(*edge) not in simplices:
                    simplices.add(edge)
            return tuple(simplices)
        if self._k == 2:
            return tuple(self._surface.labels())

        return tuple()

    @cached_method
    def change(self, k=None):
        return SimplicialHomology(surface=self._surface, k=k if k is not None else self._k, coefficients=self._coefficients, generators=self._generators, relative=self._relative, implementation=self._implementation, category=self.category())

    def boundary(self, chain):
        r"""
        Return the boundary of ``chain`` as an element of the :meth:`chain_module`.

        INPUT:

        - ``chain`` -- an element of :meth:`chain_module`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

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
                boundary += coefficient * to_C0(self._surface.point(*self._surface.opposite_edge(*edge)))
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
                        boundary -= coefficient * C1(self._surface.opposite_edge(face, edge))
            return boundary

        return self.change(k=self._k - 1).chain_module().zero()

    @cached_method
    def _chain_complex(self):
        r"""
        Return the chain complex of vector spaces that is implementing this
        homology (if the ``"generic"`` implementation has been selected.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H._chain_complex()
            Chain complex with at most 3 nonzero terms over Integer Ring

        """
        def boundary(dimension, chain):
            boundary = self.change(k=dimension).boundary(chain)
            coefficients = boundary.dense_coefficient_list(self.change(k=dimension - 1).chain_module().indices())
            return coefficients

        from sage.all import ChainComplex, matrix
        return ChainComplex({
            dimension: matrix([boundary(dimension, simplex) for simplex in self.change(k=dimension).chain_module().basis()]).transpose()
            for dimension in range(3)
        }, base_ring=self._coefficients, degree=-1)

    def zero(self):
        r"""
        Return the zero element of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.zero()
            0

        """
        return self(self.chain_module().zero())

    @cached_method
    def _homology(self):
        r"""
        Return the free module isomorphic to homology, a lift from that
        module to the chain module, and an inverse (modulo boundaries.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
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

            from sage.all import vector
            from_homology = homology.module_morphism(function=lambda x: F.from_vector(vector(list(x.lift().lift()))), codomain=F)

            def _to_homology(x):
                multiplicities = x.dense_coefficient_list(order=F.get_order())
                try:
                    cycle = cycles(multiplicities)
                except TypeError:
                    if multiplicities not in cycles:
                        raise ValueError("chain is not a cycle so it has no representation in homology")
                    raise

                return homology(cycle)


            to_homology = F.module_morphism(function=_to_homology, codomain=homology)

            for gen in homology.gens():
                assert to_homology(from_homology(gen)) == gen

            return homology, from_homology, to_homology

        raise NotImplementedError("cannot compute homology with this implementation yet")

    def _test_homology(self, **options):
        r"""
        Test that :meth:`_homology` compute homology correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
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

        from sage.all import ZZ
        if self._coefficients is not ZZ:
            return f"H{k}({self._surface}; {self._coefficients})"

        return f"H{k}({self._surface})"

    def _element_constructor_(self, x):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)

            sage: H(0)
            0
            sage: H(None)
            0

            sage: H((0, 0))
            B[(0, 0)]
            sage: H((0, 2))
            -B[(0, 0)]

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

        raise NotImplementedError

    @cached_method
    def gens(self):
        r"""
        Return generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

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

        homology, from_homology, to_homology = self._homology()
        return tuple(self(from_homology(g)) for g in homology.gens())

    def __eq__(self, other):
        if not isinstance(other, SimplicialHomologyGroup):
            return False

        return self._surface == other._surface and self._coefficients == other._coefficients and self._generators == other._generators and self._relative == other._relative and self._implementation == other._implementation and self.category() == other.category()

    def __hash__(self):
        return hash((self._surface, self._coefficients, self._generators, self._relative, self._implementation, self.category()))


def SimplicialHomology(surface, k=1, coefficients=None, generators="edge", relative=None, implementation="generic", category=None):
    r"""
    TESTS:

    Homology is unique and cached::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: T.set_immutable()
        sage: SimplicialHomology(T) is SimplicialHomology(T)
        True

    """
    return surface.homology(k, coefficients, generators, relative, implementation, category)
