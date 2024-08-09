from sage.categories.morphism import Morphism


class SurfaceMorphism(Morphism):
    pass


class SurfaceMorphism_factorization(SurfaceMorphism):
    def _image_homology(self, x, codomain=None):
        return self._factorization()._image_homology(x, codomain=codomain)
