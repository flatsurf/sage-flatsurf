from sage.categories.category import Category

class TangentBundles(Category):
    def __init__(self, surfaces):
        self._surfaces = surfaces
        super().__init__()

    def super_categories(self):
        from flatsurf.geometry.categories.surface_category import SurfaceCategory, SurfaceCategoryWithAxiom
        induced = [TangentBundles(surfaces) for surfaces in self._surfaces.super_categories() if isinstance(surfaces, (SurfaceCategory, SurfaceCategoryWithAxiom))]

        if induced:
            return induced

        # We cannot use the category of vector bundles because that requires a
        # fixed base ring. However, our surface categories do not have a fixed
        # base ring.
        from sage.categories.all import Sets
        return [Sets().Topological()]

    @property
    def ParentMethods(self):
        return getattr(self._surfaces, "TangentBundleMethods", None)

    @property
    def ElementMethods(self):
        return getattr(self._surfaces, "TangentVectorMethods", None)
