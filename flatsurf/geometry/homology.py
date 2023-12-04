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

We can also write the generators as paths that cross over the edges and
connect points on the interior of neighboring polygons::

TODO: We actually don't do that currently anymore. We rather take the actual
edges. But using that form might be useful sometimes.

    sage: H = SimplicialHomology(S, generators="voronoi")
    sage: H.gens()
    (B[(0, 1)], B[(0, 2)], B[(0, 3)], B[(0, 0)])

We can also use generators that connect points on the interior of edges of a
polygon which can be advantageous when integrating along such paths while
avoiding to integrate close to the vertices::

    sage: H = SimplicialHomology(S, generators="midpoint")  # not tested; TODO
    sage: H.gens()  # not tested; TODO

Relative homology on the unfolding of the (3, 4, 13) triangle; homology
relative to the subset of vertices::

    sage: from flatsurf import EuclideanPolygonsWithAngles, similarity_surfaces
    sage: P = EuclideanPolygonsWithAngles(3, 4, 13).an_element()
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")
    sage: H = SimplicialHomology(relative=S.singularities())  # not tested; TODO
    sage: H.gens()  # not tested; TODO

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

    def surface(self):
        return self.parent().surface()


class SimplicialHomologyClass_edge(SimplicialHomologyClass):
    pass


class SimplicialHomologyClass_voronoi(SimplicialHomologyClass):
    r"""
    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: T.set_immutable()
        sage: H = SimplicialHomology(T, generators="voronoi")
        sage: a, b = H.gens()

    The cycle ``a`` is horizontal in the square torus since it
    crosses the vertical edge 1::

        sage: a
        B[(0, 1)]
        sage: b
        B[(0, 0)]

    """
    # def safety(self):
    #     return min(self._safety(label0, edge0) for (label0, edge0) in self._surface.underlying_surface().edges())

    # def _safety(self, label, edge):
    #     polygon = self._surface.polygon(label)
    #     cross_label, cross_edge = self._surface.opposite_edge(label, edge)
    #     cross_polygon = self._surface.polygon(cross_label)
    #     cross_polygon = cross_polygon.translate(polygon.vertex(edge) - cross_polygon.vertex((cross_edge + 1) % cross_polygon.num_edges()))
    #     distance = (polygon.circumscribing_circle().center() - cross_polygon.circumscribing_circle().center()).change_ring(RR).norm()
    #     if distance == 0:
    #         from sage.all import infinity
    #         return infinity
    #     return self._convergence(label) / distance

    # def _convergence(self, label):
    #     r"""
    #     Return the radius of convergence at the point at which we develop the
    #     series for the polygon ``label``.
    #     """
    #     polygon = self._surface.polygon(label)
    #     center = polygon.circumscribing_circle().center()
    #     singularities = self._surface.singularities()

    #     from sage.all import infinity
    #     radius = infinity

    #     if not singularities:
    #         return radius

    #     closest_edges = []
    #     polygons = set([polygon])

    #     def distance(polygon, edge):
    #         # TODO: This is wrong. This is the distance to the endpoints.
    #         return min((center - polygon.vertex(edge)).change_ring(RR).norm(), (center - polygon.vertex((edge + 1) % polygon.num_edges())).change_ring(RR).norm())

    #     from heapq import heappush, heappop
    #     for e in range(polygon.num_edges()):
    #         heappush(closest_edges, (distance(polygon, e), label, str(polygon), polygon, e))

    #     while True:
    #         d, label, _, polygon, edge = heappop(closest_edges)
    #         if d >= radius:
    #             return radius

    #         vertex = self._surface.point(label, self._surface.polygon(label).vertex(edge))
    #         if vertex in singularities:
    #             radius = min(radius, (polygon.vertex(edge) - center).change_ring(RR).norm())

    #         cross_label, cross_edge = self._surface.opposite_edge(label, edge)
    #         cross_polygon = self._surface.polygon(cross_label)
    #         cross_polygon = cross_polygon.translate(polygon.vertex(edge) - cross_polygon.vertex((cross_edge + 1) % cross_polygon.num_edges()))

    #         if cross_polygon in polygons:
    #             continue

    #         polygons.add(cross_polygon)

    #         for e in range(cross_polygon.num_edges()):
    #             if e == cross_label:
    #                 continue

    #             heappush(closest_edges, (distance(cross_polygon, e), cross_label, str(cross_polygon), cross_polygon, e))


class SimplicialHomology(UniqueRepresentation, Parent):
    r"""
    Absolute and relative simplicial homology of the ``surface`` with
    ``coefficients``.

    INPUT:

    - ``surface`` -- a finite :class:`flatsurf.geometry.surface.Surface`
      without boundary

    - ``coefficients`` -- a ring (default: the integers)

    - ``generators`` -- one of ``edge``, ``voronoi`` (default: ``edge``) how
      generators are represented

    - ``relative`` -- a subset of points of the ``surface`` (default: the empty
      set)

    - ``implementation`` -- one of ``"spanning_tree"`` or ``"generic"`` (default:
      ``"spanning_tree"``); whether the homology is computed with a (very
      efficient) spanning tree algorithm or the generic homology machinery
      provided by SageMath.

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
        H₁(Translation Surface in H_1(0) built from a square; Integer Ring)

    TESTS::

        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: H = SimplicialHomology(T, implementation="spanning_tree")  # not tested; TODO
        sage: TestSuite(H).run()  # not tested; TODO

    ::

        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: H = SimplicialHomology(T, implementation="generic")
        sage: TestSuite(H).run()

    """
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

        from sage.all import Sets
        category = category or Sets()
        subset = frozenset(subset or {})

        return super().__classcall__(cls, surface, coefficients, generators, subset, implementation, category)

    def __init__(self, surface, coefficients, generators, subset, implementation, category):
        Parent.__init__(self, category=category)

        from sage.categories.all import Rings
        if coefficients not in Rings():
            raise TypeError("coefficients must be a ring")

        if generators == "edge":
            self.Element = SimplicialHomologyClass_edge
        elif generators == "voronoi":
            self.Element = SimplicialHomologyClass_voronoi
        else:
            raise NotImplementedError("cannot represent homology with these generators yet")

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
            return self._simplices_points()

        if dimension == 1:
            return self._simplices_segments()

        if dimension == 2:
            return self._simplices_polygons()

        return tuple()

    def _simplices_points(self):
        if self._generators == "edge":
            return tuple(set(self._surface.point(label, self._surface.polygon(label).vertex(edge)) for (label, edge) in self._surface.edges()))

        if self._generators == "voronoi":
            for label in self._surface.labels():
                polygon = self._surface.polygon(label)
                # TODO: This fails if the center of the circumscribing circle
                # is not in the polygon; we should be more explicit here and
                # check this condition properly earlier.
                # self._surface.surface_point(label, polygon.circumscribing_circle().center())

            return tuple(self._surface.labels())

        raise NotImplementedError

    def _simplices_segments(self):
        if self._generators in ["edge", "voronoi"]:
            # When "edge", then the edges are the generators.
            # When "voronoi", then the paths crossing the edges are the generators.
            simplices = set()
            for edge in self._surface.edges():
                if self._surface.opposite_edge(*edge) not in simplices:
                    simplices.add(edge)
            return tuple(simplices)

        raise NotImplementedError

    def _simplices_polygons(self):
        if self._generators == "edge":
            return tuple(self._surface.labels())

        if self._generators == "voronoi":
            return tuple(self._surface.edges())

        raise NotImplementedError

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
            2*B[Vertex 0 of polygon 0]
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
            for gen, coefficient in chain:
                boundary += coefficient * self._boundary_segment(gen)

            return boundary

        if chain.parent() == self.chain_module(dimension=2):
            C1 = self.chain_module(dimension=1)
            boundary = C1.zero()
            for gen, coefficient in chain:
                boundary += coefficient * self._boundary_polygon(gen)

            return boundary

        if chain.parent() == self.chain_module(dimension=0):
            return self.chain_module(dimension=-1).zero()

        # The boundary for all other chains is trivial but we have no way to
        # tell whether the boundary is a 2-dimensional chain (which lives in a
        # non-trivial module) or a chain living in a trivial module.
        raise NotImplementedError("cannot compute boundary of this chain yet")

    def _boundary_segment(self, gen):
        if self._generators == "edge":
            C0 = self.chain_module(dimension=0)
            label, edge = gen
            opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)
            return C0(self._surface.point(opposite_label, self._surface.polygon(opposite_label).vertex(opposite_edge))) - C0(self._surface.point(label, self._surface.polygon(label).vertex(edge)))

        if self._generators == "voronoi":
            C0 = self.chain_module(dimension=0)
            label, edge = gen
            opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)

            return C0(opposite_label) - C0(label)

        raise NotImplementedError

    def _boundary_polygon(self, gen):
        if self._generators == "edge":
            C1 = self.chain_module(dimension=1)
            boundary = C1.zero()
            face = gen
            for edge in range(len(self._surface.polygon(face).vertices())):
                if (face, edge) in C1.indices():
                    boundary += C1((face, edge))
                else:
                    boundary -= C1(self._surface.opposite_edge(face, edge))
            return boundary

        if self._generators == "voronoi":
            C1 = self.chain_module(dimension=1)
            boundary = C1.zero()
            label, vertex = gen
            # The counterclockwise walk around "vertex" is a boundary.
            while True:
                edge = (vertex - 1) % len(self._surface.polygon(label).vertices())
                opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)
                if (label, edge) in C1.indices():
                    boundary += C1((label, edge))
                else:
                    boundary -= C1((opposite_label, opposite_edge))

                if (opposite_label, opposite_edge) == gen:
                    return boundary

                label, vertex = opposite_label, opposite_edge

        raise NotImplementedError

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
            H₁(Translation Surface in H_1(0) built from a square; Integer Ring)

        """
        return f"H₁({self._surface}; {self._coefficients})"

    def _element_constructor_(self, x):
        if x == 0 or x is None:
            return self.element_class(self, self.chain_module(1).zero())

        if x.parent() in (self.chain_module(0), self.chain_module(1), self.chain_module(2)):
            return self.element_class(self, x)

        raise NotImplementedError

    def gen(self, n, dimension=1):
        r"""
        Return the ``n``-th generator of homology in ``dimension``, i.e., the
        ``n``-th element of :meth:`gens`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)

        ::

            sage: H.gen(0)
            B[(0, 1)]

        """
        return self.gens()[n]

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
            (B[Vertex 0 of polygon 0],)

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
