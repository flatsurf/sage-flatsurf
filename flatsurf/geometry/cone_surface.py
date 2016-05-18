from flatsurf.geometry.similarity_surface import (
    SimilaritySurface_generic, 
    SimilaritySurface_polygons_and_gluings,
    SimilaritySurface_wrapper)

from flatsurf.geometry.surface import SurfaceType, convert_to_type

class ConeSurface_generic(SimilaritySurface_generic):
    r"""
    A Euclidean cone surface.
    """
    def surface_type(self):
        return SurfaceType.CONE

    def angles(self):
        r"""
        Return the set of angles around the vertices of the surface.

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: sfg.translation_surfaces.regular_octagon().angles()
            [3]
        """
        if not self.is_finite():
            raise NotImplementedError("the set of edges is infinite!")

        edges = [(p,e) for p in self.polygon_labels() for e in range(self.polygon(p).num_edges())]
        edges = set(edges)
        angles = []
        while edges:
            p,e = edges.pop()
            angle = self.polygon(p).angle(e)
            pp,ee = self.opposite_edge(p,(e-1)%self.polygon(p).num_edges())
            while pp != p or ee != e:
                edges.remove((pp,ee))
                angle += self.polygon(pp).angle(ee)
                pp,ee = self.opposite_edge(pp,(ee-1)%self.polygon(pp).num_edges())
            angles.append(angle)
        return angles

class ConeSurface_polygons_and_gluings(
        SimilaritySurface_polygons_and_gluings,
        ConeSurface_generic):
    pass

class ConeSurface_wrapper(
        SimilaritySurface_wrapper,
        ConeSurface_generic):
    pass

def convert_to_cone_surface(surface):
    r"""
    Returns a cone surface version of the provided surface.
    """
    if surface.is_finite():
        return ConeSurface_polygons_and_gluings(surface)
    else:
        return ConeSurface_wrapper(surface)

