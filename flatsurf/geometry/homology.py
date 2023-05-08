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

We can also write the generatorls as paths that cross over the edges and
connect points on the interior of neighboring polygons::

    sage: H = SimplicialHomology(S, generators="interior")
    sage: H.gens()

We can also use generators that connect points on the interior of edges of a
polygon which can be advantageous when integrating along such paths while
avoiding to integrate close to the vertices::

    sage: H = SimplicialHomology(S, generators="midpoint")
    sage: H.gens()

Relative homology on the unfolding of the (3, 4, 13) triangle; homology
relative to the subset of vertices::


    sage: from flatsurf import EquiangularPolygons, similarity_surfaces
    sage: P = flatsurf.EquiangularPolygons(3, 4, 13).an_element()
    sage: S = flatsurf.similarity_surfaces.billiard(P, rational=True).minimal_cover(cover_type="translation")
    sage: H = SimplicialHomology(relative=S.singularities())
    sage: H.gens()

TODO: Add examples.
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2023 Julian Rüth
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
from sage.structure.unique_representation import UniqueRepresentation

from sage.misc.cachefunc import cached_method


# TODO: We should have absolute and relative (to a subset of vertices) homology.
# TODO: We implement everything in terms of generators given by the edges of
# the polygons. However, we should have views that pretend that the generators
# are paths between adjacent polygons, and a view with generators given by
# paths between midpoints of edges of polygons.

class SimplicialHomologyClass(Element):
    # TODO: Use the algorithm from GL2ROrbitClosure._spanning_tree to compute a
    # basis of homology and a projection map. Or better, have an algorithm
    # keyword to use the generic implementation and the spanning tree
    # implementation.
    # TODO: Use https://github.com/flatsurf/sage-flatsurf/pull/114/files to
    # force the representatives to live in particular subgraph of the dual
    # graph.
    def __init__(self, parent, chain):
        super().__init__(parent)

        self._chain = chain

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
        Compare homology classes with respect to ``op``.

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
            return not (self == other)

        if op == op_EQ:
            return self._homology() == other._homology()

        raise NotImplementedError

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

        """
        # TODO: Check that gen is a generator.

        for g, c in zip(gen._homology(), self._homology()):
            if g:
                assert g == 1
                return c

        raise ValueError("gen must be a generator of homology")

    # def voronoi_path(self):
    #     r"""
    #     Return this class as a vector over a basis of homology formed by paths
    #     inside Voronoi cells.

    #     EXAMPLES::

    #         sage: from flatsurf import translation_surfaces, SimplicialHomology
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()
    #         sage: H = SimplicialHomology(T)
    #         sage: a, b = H.gens()

    #     The cycle ``a`` is the vertical in the square torus.

    #         sage: a.voronoi_path()
    #         B[((0, 0), (1, 2), (0, 1), (1, 0))]
    #         sage: b.voronoi_path()
    #         B[((0, 1), (1, 0), (0, 2), (1, 1))]

    #     """
    #     # TODO: Don't expose this here but in a separate view that uses different generators.
    #     from sage.all import FreeModule, vector

    #     homology, _, _ = self.parent()._homology()
    #     M = FreeModule(self.parent()._coefficients, self.parent()._paths(voronoi=True))
    #     return M.from_vector(vector(self._homology()))

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

    # TODO: This probably makes no sense.
    # def _path(self, voronoi=False):
    #     r"""
    #     Return this generator as a path.

    #     If ``voronoi``, the path is formed by walking from midpoints of edges while staying inside Voronoi cells.

    #     EXAMPLES::

    #         sage: from flatsurf import translation_surfaces, SimplicialHomology
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()
    #         sage: H = SimplicialHomology(T)
    #         sage: a, b = H.gens()

    #     The chosen generators of homology, correspond to the edge (0, 0), i.e.,
    #     the diagonal with vector (1, 1), and the top horizontal edge (0, 1)
    #     with vector (-1, 0)::

    #         sage: a._path()
    #         ((0, 0),)
    #         sage: b._path()
    #         ((0, 1),)

    #     Lifting the former to a path in Voronoi cells, we consider the midpoint
    #     of (0, 0) as the midpoint of (1, 0) and walk across the triangle 1 to
    #     get to the midpoint of (1, 2). We consider the midpoint of (1, 2) as
    #     the midpoint of (0, 2) and walk across the triangle 0 to the midpoint
    #     of (0, 1). Finally, we consider that midpoint to be the midpoint of (1,
    #     1) and walk across 1 to the midpoint of (1, 0) which is where we
    #     started::

    #         sage: a._path(voronoi=True)
    #         ((0, 0), (1, 2), (0, 1), (1, 0))

    #     Similarly, to write the other path as a path inside Voronoi cells, we
    #     start from the midpoint of (0, 1) and consider it as the midpoint of
    #     (1, 1). We walk across triangle 1 to get to the midpoint of (1, 0). We
    #     consider that to be the midpoint of (0, 0) and walk across 0 to the
    #     midpoint of (0, 2). We consider that to be the midpoint of (1, 2) and
    #     walk across triangle 1 to get to the midpoint of (1, 1) which closes
    #     the loop::

    #         sage: b._path(voronoi=True)
    #         ((0, 1), (1, 0), (0, 2), (1, 1))

    #     TODO: Note from the above discussion that we can have consecutive steps
    #     in the same triangle; in the above the last and the first. We should
    #     collapse these.

    #     ::

    #         sage: from flatsurf import EquiangularPolygons, similarity_surfaces
    #         sage: E = EquiangularPolygons(3, 4, 13)
    #         sage: P = E.an_element()
    #         sage: T = similarity_surfaces.billiard(P, rational=True).minimal_cover(cover_type="translation").erase_marked_points()

    #         sage: H = SimplicialHomology(T)
    #         sage: a = H.gens()[0]; a
    #         B[(0, 0)] - B[(30, 1)]
    #         sage: a._path()
    #         ((31, 2), (0, 0))
    #         sage: a._path(voronoi=True)
    #         ((31, 2), (30, 0), (25, 2), (5, 0), (13, 1), (18, 0), (7, 0), (12, 2), (4, 0), (24, 1), (28, 0), (21, 2), (14, 0), (9, 2), (3, 2),
    #          (0, 0), (3, 1), (23, 1), (26, 1), (27, 1), (16, 0), (19, 2), (18, 1), (13, 0), (14, 1), (21, 1), (6, 2), (2, 0), (20, 1), (11, 0), (17, 0), (30, 1))

    #     """
    #     edges = []

    #     for edge in self._chain.support():
    #         coefficient = self._chain[edge]
    #         if coefficient == 0:
    #             continue
    #         if coefficient < 0:
    #             coefficient *= -1
    #             edge = self.parent()._surface.opposite_edge(*edge)

    #         if coefficient != 1:
    #             raise NotImplementedError

    #         edges.append(edge)

    #     path = [edges.pop()]

    #     surface = self.surface()

    #     while edges:
    #         # Lengthen the path by searching for an edge that starts at the
    #         # vertex at which the constructed path ends.
    #         previous = path[-1]
    #         vertex = surface.singularity(*surface.opposite_edge(*previous))
    #         for edge in edges:
    #             if surface.singularity(*edge) == vertex:
    #                 path.append(edge)
    #                 edges.remove(edge)
    #                 break
    #         else:
    #             raise NotImplementedError

    #     if not voronoi:
    #         return tuple(path)

    #     # Instead of walking the edges of the triangulation, we walk
    #     # across faces between the midpoints of the edges to construct a path
    #     # inside Voronoi cells.
    #     voronoi_path = [path[0]]

    #     for previous, edge in zip(path, path[1:] + path[:1]):
    #         previous = self.surface().opposite_edge(*previous)
    #         while previous != edge:
    #             previous = previous[0], (previous[1] + 2) % 3
    #             voronoi_path.append(previous)
    #             previous = self.surface().opposite_edge(*previous)
    #         voronoi_path.append(edge)

    #     voronoi_path.pop()

    #     return tuple(voronoi_path)

    def surface(self):
        return self.parent().surface()


class SimplicialHomology(UniqueRepresentation, Parent):
    r"""
    Absolute and relative simplicial homology of the ``surface`` with
    ``coefficients``.

    INPUT:

    - ``surface`` -- a finite :class:`flatsurf.geometry.surface.Surface`
      without boundary

    - ``coefficients`` -- a ring (default: the integers)

    - ``generators`` -- one of ``edge``, ``interior``, ``midpoint`` (default:
      ``edge``) how generators are represented

    - ``relative`` -- a subset of points of the ``surface`` (default: the empty
      set)

    - ``implementation`` -- one of ``"spanning_tree"`` or ``"generic"`` (default:
      ``"spanning_tree"``); whether the homology is computed with a (very
      efficient) spanning tree algorithm or the generic homology machinery
      provided by SageMath.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))

    Surfaces must be immutable to compute their homology::

        sage: T = T.copy(mutable=True)
        sage: SimplicialHomology(T)
        Traceback (most recent call last):
        ...
        ValueError: surface must be immutable to compute homology

    ::

        sage: T.set_immutable()
        sage: SimplicialHomology(T)
        H₁(TranslationSurface built from 1 polygon; Integer Ring)

    TESTS::

        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: H = SimplicialHomology(implementation="spanning_tree")
        sage: TestSuite(H).run()

    ::

        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: H = SimplicialHomology(implementation="generic")
        sage: TestSuite(H).run()

    """
    Element = SimplicialHomologyClass

    @staticmethod
    # TODO: implementation should default to spanning_tree
    def __classcall__(cls, surface, coefficients=None, generators="edge", subset=None, implementation="generic", category=None):
        r"""
        Normalize parameters used to construct homology.

        TESTS:

        Homology is unique and cached::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: SimplicialHomology(T) is SimplicialHomology(T)
            True

        """
        if surface.is_mutable():
            raise ValueError("surface must be immutable to compute homology")

        from sage.all import ZZ
        coefficients = coefficients or ZZ

        from sage.all import SetsWithPartialMaps
        category = category or SetsWithPartialMaps()
        subset = frozenset(subset or {})

        return super().__classcall__(cls, surface, coefficients, generators, subset, implementation, category)

    def __init__(self, surface, coefficients, generators, subset, implementation, category):
        Parent.__init__(self, category=category)

        from sage.categories.all import Rings
        if coefficients not in Rings():
            raise TypeError("coefficients must be a ring")

        if generators not in ["edge",]:
            raise NotImplementedError("cannot represented homology with these generators yet")

        if subset:
            raise NotImplementedError("cannot compute relative homology yet")

        if implementation not in ["generic", "spanning_tree"]:
            raise NotImplementedError("cannot compute homology with this implementation yet")

        self._surface = surface
        self._coefficients = coefficients
        self._generators = generators
        self._subset = subset
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
    def chain_module(self, dimension=1):
        r"""
        Return the free module of simplicial ``dimension``-chains of the
        triangulation, i.e., formal sums of ``dimension``-simplicies, e.g., formal
        sums of edges of the triangulation.

        INPUT:

        - ``dimension`` -- an integer (default: ``1``)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.chain_module()
            Free module generated by {(0, 1), (0, 0)} over Integer Ring

        """
        from sage.all import FreeModule
        return FreeModule(self._coefficients, self.simplices(dimension))

    @cached_method
    def simplices(self, dimension=1):
        r"""
        Return the ``dimension``-simplices that form the basis of
        :meth:`chain_module`.

        INPUT:

        - ``dimension`` -- an integer (default: ``1``)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.simplices()
            ((0, 1), (0, 0))

        """
        if dimension == 0:
            return tuple(set(self._surface.singularity(*edge) for edge in self._surface.edge_iterator()))
        if dimension == 1:
            simplices = set()
            for edge in self._surface.edge_iterator():
                if self._surface.opposite_edge(*edge) not in simplices:
                    simplices.add(edge)
            return tuple(simplices)
        if dimension == 2:
            return tuple(self._surface.label_iterator())

        return tuple()

    def boundary(self, chain):
        r"""
        Return the boundary of ``chain`` as an element of the :meth:`chain_module`.

        INPUT:

        - ``chain`` -- an element of :meth:`chain_module`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)

        ::

            sage: c = H.chain_module(dimension=0).an_element(); c
            2*B[singularity with vertex equivalence class frozenset({(0, 1), (0, 2), (0, 3), (0, 0)})]
            sage: H.boundary(c)
            0

        ::

            sage: c = H.chain_module(dimension=1).an_element(); c
            2*B[(0, 0)] + 2*B[(0, 1)]
            sage: H.boundary(c)
            0

        ::

            sage: c = H.chain_module(dimension=2).an_element(); c
            2*B[0]
            sage: H.boundary(c)
            0

        """
        if chain.parent() == self.chain_module(dimension=1):
            C0 = self.chain_module(dimension=0)
            boundary = C0.zero()
            for edge, coefficient in chain:
                boundary += coefficient * C0(self._surface.singularity(*self._surface.opposite_edge(*edge)))
                boundary -= coefficient * C0(self._surface.singularity(*edge))
            return boundary

        if chain.parent() == self.chain_module(dimension=2):
            C1 = self.chain_module(dimension=1)
            boundary = C1.zero()
            for face, coefficient in chain:
                for edge in range(self._surface.polygon(face).num_edges()):
                    if (face, edge) in C1.indices():
                        boundary += coefficient * C1((face, edge))
                    else:
                        boundary -= coefficient * C1(self._surface.opposite_edge(face, edge))
            return boundary

        if chain.parent() == self.chain_module(dimension=0):
            return self.chain_module(dimension=-1).zero()

        # The boundary for all other chains is trivial but we have no way to
        # tell whether the boundary is a 2-dimensional chain (which lives in a
        # non-trivial module) or a chain living in a trivial module.
        raise NotImplementedError("cannot compute boundary of this chain yet")

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
        def boundary(dimension, simplex):
            C = self.chain_module(dimension)
            chain = C.basis()[simplex]
            boundary = self.boundary(chain)
            coefficients = boundary.dense_coefficient_list(self.chain_module(dimension-1).indices())
            return coefficients

        from sage.all import ChainComplex, matrix
        return ChainComplex({
            dimension: matrix([boundary(dimension, simplex) for simplex in self.simplices(dimension)]).transpose()
            for dimension in range(3)
        }, base_ring=self._coefficients, degree=-1)

    @cached_method
    def _homology(self, dimension=1):
        r"""
        Return the a free module isomorphic to homology, a lift from that
        module to the chain module, and an inverse (modulo boundaries.)

        INPUT:

        - ``dimension`` -- an integer (default: ``1``)

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

            cycles = C.differential(dimension).transpose().kernel()
            boundaries = C.differential(dimension + 1).transpose().image()
            homology = cycles.quotient(boundaries)

            F = self.chain_module(dimension)

            from sage.all import vector
            from_homology = homology.module_morphism(function=lambda x: F.from_vector(vector(list(x.lift().lift()))), codomain=F)
            to_homology = F.module_morphism(function=lambda x: homology(cycles(x.dense_coefficient_list(order=F.get_order()))), codomain=homology)

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

        for dimension in [0, 1, 2]:
            homology, from_homology, to_homology = self._homology(dimension=dimension)
            chains = self.chain_module(dimension)

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
            H₁(TranslationSurface built from 1 polygon; Integer Ring)

        """
        return f"H₁({self._surface}; {self._coefficients})"

    def _element_constructor_(self, x):
        if x == 0 or x is None:
            return self.element_class(self, self.chain_module(1).zero())

        if x.parent() in (self.chain_module(0), self.chain_module(1), self.chain_module(2)):
            return self.element_class(self, x)

        raise NotImplementedError

    @cached_method
    def gens(self, dimension=1):
        r"""
        Return generators of homology in ``dimension``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)

        ::

            sage: H.gens(dimension=0)
            (B[singularity with vertex equivalence class frozenset({(0, 1), (0, 2), (0, 3), (0, 0)})],)

        ::

            sage: H.gens(dimension=1)
            (B[(0, 1)], B[(0, 0)])

        ::

            sage: H.gens(dimension=2)
            (B[0],)

        """
        if dimension < 0 or dimension > 2:
            return ()

        homology, from_homology, to_homology = self._homology(dimension)
        return tuple(self(from_homology(g)) for g in homology.gens())

    # @cached_method
    # def _paths(self, voronoi=False):
    #     r"""
    #     Return the generators of homology in dimension 1 as paths.

    #     If ``voronoi``, the paths are inside Voronoi cells without touching the vertices of the triangulation.

    #     EXAMPLES::

    #         sage: from flatsurf import translation_surfaces, SimplicialHomology
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()
    #         sage: H = SimplicialHomology(T)
    #         sage: H._paths()
    #         [((0, 0),), ((0, 1),)]

    #     """
    #     TODO: This does not really make sense probably.
    #     return [gen._path(voronoi=voronoi) for gen in self.gens()]
