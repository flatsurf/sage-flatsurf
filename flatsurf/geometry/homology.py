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
from sage.structure.unique_representation import UniqueRepresentation

from sage.misc.cachefunc import cached_method


class SimplicialHomologyClass(Element):
    def __init__(self, parent, chain):
        super().__init__(parent)

        self._chain = chain

    @cached_method
    def _homology(self):
        r"""
        Return the coefficients of this element in terms of the generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
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
                return c

        raise ValueError("gen must be a generator of homology")

    def voronoi_path(self):
        r"""
        Return this class as a vector over a basis of homology formed by paths
        inside Voronoi cells.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: a,b = H.gens()
            sage: a.voronoi_path()
            B[((0, 0),)]
            sage: b.voronoi_path()
            B[((0, 1),)]

        """
        from sage.all import FreeModule, vector

        homology, _, _ = self.parent()._homology()
        M = FreeModule(self.parent()._coefficients, self.parent()._voronoi_paths())
        return M.from_vector(vector(self._homology()))

    def _voronoi_path(self):
        r"""
        Return this generator as a path inside the Voronoi cells.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: a,b = H.gens()
            sage: a._voronoi_path()
            ((0, 0),)
            sage: b._voronoi_path()
            ((0, 1),)

        """
        edges = []

        for edge in self._chain.support():
            coefficient = self._chain[edge]
            if coefficient == 0:
                continue
            if coefficient < 0:
                coefficient *= -1
                edge = self.parent()._surface.opposite_edge(*edge)

            if coefficient != 1:
                raise NotImplementedError

            edges.append(edge)

        path = [edges.pop()]

        while edges:
            previous = self._surface.singularity(*self._surface.opposite_edge(*path[-1]))
            for edge in edges:
                if self._surface.singularity(*edge) == previous:
                    path.append(edge)
                    edges.remove(edge)
                    break
            else:
                raise NotImplementedError

        return tuple(path)


class SimplicialHomology(UniqueRepresentation, Parent):
    r"""
    Absolute simplicial homology of the ``surface`` with ``coefficients``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))

    Currently, surfaces must be Delaunay triangulated to compute their homology::

        sage: SimplicialHomology(T)
        Traceback (most recent call last):
        ...
        NotImplementedError: surface must be Delaunay triangulated

    Surfaces must be immutable to compute their homology::

        sage: T = T.delaunay_triangulation()
        sage: SimplicialHomology(T)
        Traceback (most recent call last):
        ...
        ValueError: surface must be immutable to compute homology

    ::

        sage: T.set_immutable()
        sage: SimplicialHomology(T)
        H₁(TranslationSurface built from 2 polygons; Integer Ring)

    """
    Element = SimplicialHomologyClass

    @staticmethod
    def __classcall__(cls, surface, coefficients=None, category=None):
        r"""
        Normalize parameters used to construct homology.

        TESTS:

        Homology is unique and cached::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: SimplicialHomology(T) is SimplicialHomology(T)
            True

        """
        if surface.is_mutable():
            raise ValueError("surface must be immutable to compute homology")

        from sage.all import ZZ, SetsWithPartialMaps
        return super().__classcall__(cls, surface, coefficients or ZZ, category or SetsWithPartialMaps())

    def __init__(self, surface, coefficients, category):
        if surface != surface.delaunay_triangulation():
            # TODO: This is a silly limitation in here.
            raise NotImplementedError("surface must be Delaunay triangulated")

        Parent.__init__(self, category=category)

        self._surface = surface
        self._coefficients = coefficients

    @cached_method
    def chain_module(self, dimension=1):
        r"""
        Return the free module of simplicial ``dimension``-chains of the
        triangulation, i.e., formal sums of ``dimension``-simplicies, e.g., formal
        sums of edges of the triangulation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.chain_module()
            Free module generated by {(0, 0), (0, 1), (0, 2)} over Integer Ring

        """
        from sage.all import FreeModule
        return FreeModule(self._coefficients, self.simplices(dimension))

    @cached_method
    def simplices(self, dimension=1):
        r"""
        Return the ``dimension``-simplices that form the basis of
        :meth:`chain_module`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H.simplices()
            ((0, 0), (0, 1), (0, 2))

        """
        if dimension == 0:
            return tuple(set(self._surface.singularity(*edge) for edge in self._surface.edge_iterator()))
        if dimension == 1:
            return tuple(edge for edge in self._surface.edge_iterator() if edge[0] < self._surface.opposite_edge(*edge)[0])
        if dimension == 2:
            return tuple(self._surface.label_iterator())

        return tuple()

    def boundary(self, chain):
        r"""
        Return the boundary of ``chain`` as an element of the :meth:`chain_module`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)

        ::

            sage: c = H.chain_module(dimension=0).an_element(); c
            2*B[singularity with vertex equivalence class frozenset({(0, 1), (0, 2), (1, 2), (0, 0), (1, 0), (1, 1)})]
            sage: H.boundary(c)
            0

        ::

            sage: c = H.chain_module(dimension=1).an_element(); c
            2*B[(0, 0)] + 2*B[(0, 1)] + 3*B[(0, 2)]
            sage: H.boundary(c)
            0

        ::

            sage: c = H.chain_module(dimension=2).an_element(); c
            2*B[0] + 2*B[1]
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
                for edge in range(3):
                    if (face, edge) in C1.indices():
                        boundary += coefficient * C1((face, edge))
                    else:
                        boundary -= coefficient * C1(self._surface.opposite_edge((face, edge)))
            return boundary

        return self.chain_module(dimension=-1).zero()

    @cached_method
    def _chain_complex(self):
        r"""
        Return the chain complex of vector spaces that is implementing this homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
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
        Return the a free module isomorphic to homology, a map from that module
        to the chain module, and an inverse (modulo boundaries.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H._homology()
            (Finitely generated module V/W over Integer Ring with invariants (0, 0),
             Generic morphism:
               From: Finitely generated module V/W over Integer Ring with invariants (0, 0)
               To:   Free module generated by {(0, 0), (0, 1), (0, 2)} over Integer Ring,
             Generic morphism:
               From: Free module generated by {(0, 0), (0, 1), (0, 2)} over Integer Ring
               To:   Finitely generated module V/W over Integer Ring with invariants (0, 0))

        """
        C = self._chain_complex()

        cycles = C.differential(dimension).transpose().kernel()
        boundaries = C.differential(dimension + 1).transpose().image()
        homology = cycles.quotient(boundaries)

        F = self.chain_module(dimension)

        from sage.all import vector
        from_homology = homology.module_morphism(function=lambda x: F.from_vector(vector(list(x.lift().lift()))), codomain=F)
        to_homology = F.module_morphism(function=lambda x: homology(x.dense_coefficient_list()), codomain=homology)

        return homology, from_homology, to_homology

    def _repr_(self):
        return f"H₁({self._surface}; {self._coefficients})"

    def _element_constructor_(self, x):
        if x.parent() in (self.chain_module(0), self.chain_module(1), self.chain_module(2)):
            return self.element_class(self, x)

        raise NotImplementedError

    @cached_method
    def gens(self, dimension=1):
        r"""
        Return generators of homology in ``dimension``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)

        ::

            sage: H.gens(dimension=0)
            (B[singularity with vertex equivalence class frozenset({(0, 1), (0, 2), (1, 2), (0, 0), (1, 0), (1, 1)})],)

        ::

            sage: H.gens(dimension=1)
            (B[(0, 0)], B[(0, 1)])

        ::

            sage: H.gens(dimension=2)
            (B[0] + B[1],)

        """
        if dimension < 0 or dimension > 2:
            return ()

        homology, from_homology, to_homology = self._homology(dimension)
        return tuple(self(from_homology(g)) for g in homology.gens())

    @cached_method
    def _voronoi_paths(self):
        r"""
        Return the generators of homology in dimension 1 as paths in the
        Voronoi cells.

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: H = SimplicialHomology(T)
            sage: H._voronoi_paths()
            [((0, 0),), ((0, 1),)]

        """
        return [gen._voronoi_path() for gen in self.gens()]
