from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation


class SimplicialHomologyClass(Element):
    def __init__(self, parent):
        super().__init__(parent)

    def _repr_(self):
        return "[?]"


class SimplicialHomology(UniqueRepresentation, Parent):
    Element = SimplicialHomologyClass

    @staticmethod
    def __classcall__(cls, surface, coefficients=None, category=None):
        from sage.all import ZZ, SetsWithPartialMaps
        return super().__classcall__(cls, surface, coefficients or ZZ, category or SetsWithPartialMaps())

    def __init__(self, surface, coefficients, category):
        if surface != surface.delaunay_triangulation():
            raise NotImplementedError("Surface must be Delaunay triangulated")

        Parent.__init__(self, category=category)

        self._surface = surface
        self._coefficients = coefficients

        self._bases = {
            0: list(set(self._surface.singularity(*edge) for edge in self._surface.edge_iterator())),
            1: list(edge for edge in self._surface.edge_iterator() if edge[0] < self._surface.opposite_edge(*edge)[0]),
            2: list(self._surface.label_iterator()),
        }

        from sage.all import FreeModule
        self._chains = {
            dimension: FreeModule(coefficients, len(self._bases[dimension]))
            for dimension in self._bases.keys()
        }

        from sage.all import ChainComplex, matrix
        self._chain_complex = ChainComplex({
            0: matrix(0, len(self._bases[0])),
            1: matrix(len(self._bases[1]), len(self._bases[0]), [
                self._boundary(1, edge) for edge in self._bases[1]
                ]).transpose(),
            2: matrix(len(self._bases[2]), len(self._bases[1]), [
                self._boundary(2, edge) for edge in self._bases[2]
                ]).transpose(),
        }, base_ring=coefficients, degree=-1)
        self._homology = self._chain_complex.homology(generators=True)

    def _repr_(self):
        return f"H₁({self._surface}; {self._coefficients})"

    def _element_constructor_(self, x):
        if x != 0:
            raise NotImplementedError
        return self.element_class(self)

    def basis(self):
        for _, generator in self._chain_complex.homology(deg=1, generators=True):
            yield self._lift_to_cycle(generator.vector(degree=1))

    def _lift_to_cycle(self, vector):
        edges = []

        for (i, c) in enumerate(vector):
            if c == 0:
                continue
            if c == -1:
                edges.append(self._surface.opposite_edge(*self._bases[1][i]))
                continue
            if c == 1:
                edges.append(self._bases[1][i])
                continue

            raise NotImplementedError

        cycle = [edges.pop()]

        while edges:
            previous = self._surface.singularity(*self._surface.opposite_edge(*cycle[-1]))
            for edge in edges:
                if self._surface.singularity(*edge) == previous:
                    cycle.append(edge)
                    edges.remove(edge)
                    break
            else:
                raise NotImplementedError

        return cycle

    def _boundary(self, degree, element):
        r"""
        Return the boundary of a basis element of degree as an element
        of the corresponding vector space of the chain complex.
        """
        if degree == 0:
            return []

        if degree == 1:
            return self._vector(0, self._surface.singularity(*self._surface.opposite_edge(*element))) - self._vector(0, self._surface.singularity(*element))

        if degree == 2:
            return self._vector(1, (element, 0)) + self._vector(1, (element, 1)) + self._vector(1, (element, 2))

    def _vector(self, degree, element):
        if degree == 0:
            return self._chains[degree].gen(self._bases[degree].index(element))

        if degree == 1:
            if element in self._bases[degree]:
                return self._chains[degree].gen(self._bases[degree].index(element))
            return -self._chains[degree].gen(self._bases[degree].index(self._surface.opposite_edge(element)))

        if degree == 2:
            return self._chains[degree].gen(self._bases[degree].index(element))
