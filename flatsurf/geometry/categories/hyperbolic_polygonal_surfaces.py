from flatsurf.geometry.categories.surface_category import SurfaceCategory


class HyperbolicPolygonalSurfaces(SurfaceCategory):
    def super_categories(self):
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces

        return [PolygonalSurfaces()]

    class ParentMethods:
        def plot(self):
            return sum(polygon.plot() for polygon in self.polygons())

        def _an_element_(self):
            label = next(iter(self.labels()))
            polygon = self.polygon(label)

            return self(label, polygon.vertices()[0])
