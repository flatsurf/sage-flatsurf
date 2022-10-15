from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation


class SimplicialCohomologyClass(Element):
    def __init__(self, parent):
        super().__init__(parent)

    def _repr_(self):
        return "[?]"

    def vector(self):
        from sage.all import ZZ
        # TODO
        return [ZZ(1) for b in self.parent().homology().basis()]


class SimplicialCohomology(UniqueRepresentation, Parent):
    Element = SimplicialCohomologyClass

    @staticmethod
    def __classcall__(cls, surface, coefficients=None, category=None):
        from sage.all import QQ
        return super().__classcall__(cls, surface, coefficients or QQ, category or SetsWithPartialMaps())

    def __init__(self, surface, coefficients, category):
        if surface != surface.delaunay_triangulation():
            raise NotImplementedError("Surface must be Delaunay triangulated")

        Parent.__init__(self, category=category)

        self._surface = surface
        self._coefficients = coefficients

    def _repr_(self):
        return f"H¹({self._surface}; {self._coefficients})"

    def _element_constructor_(self, x):
        if x != 0:
            raise NotImplementedError
        return self.element_class(self)

    def homology(self):
        from flatsurf.geometry.homology import SimplicialHomology
        return SimplicialHomology(self._surface)
