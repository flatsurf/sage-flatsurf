from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_generic
from sage.all import MatrixGroup

class VeechGroup(FinitelyGeneratedMatrixGroup_generic):
    def __init__(self, idt, category=None):
        self._idt = idt

        gens = list(idt.gens())
        matrix_group = MatrixGroup(gens)
        # TODO: make compatible with gap
        assert isinstance(matrix_group, FinitelyGeneratedMatrixGroup_generic)
        degree = matrix_group.degree()
        base_ring = matrix_group.base_ring()
        category = category or matrix_group.category()
        super().__init__(degree, base_ring, gens, category)

    def genus(self):
        return self._idt.genus()

    def orbifold_points(self):
        return self._idt.orbifold_points()

    def surface(self):
        return self._idt._surface_original

    def cusps(self):
        return self._idt.cusps()

    def ncusps(self):
        return len(self.cusps())

    def orbifold_euler_characteristic(self):
        return self._idt.orbifold_euler_characteristic()