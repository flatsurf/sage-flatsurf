from flatsurf.geometry.categories.surface_category import SurfaceCategory


class HyperbolicIsometrySurfaces(SurfaceCategory):
    def super_categories(self):
        from flatsurf.geometry.categories.hyperbolic_polygonal_surfaces import HyperbolicPolygonalSurfaces

        return [HyperbolicPolygonalSurfaces()]

    class ParentMethods:
        def cusps(self):
            return set(vertex for vertex in self.vertices() if next(iter(vertex.representatives()))[1].is_ideal())
