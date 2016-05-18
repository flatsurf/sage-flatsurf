from flatsurf.geometry.similarity_surface import (
    SimilaritySurface_generic, 
    SimilaritySurface_polygons_and_gluings,
    SimilaritySurface_wrapper)

from flatsurf.geometry.surface import SurfaceType, convert_to_type

class HalfDilationSurface_generic(SimilaritySurface_generic):
    r"""
    A generic half translation surface
    """
    def surface_type(self):
        return SurfaceType.HALF_DILATION
    
    def GL2R_mapping(self, matrix):
        r"""
        Apply a 2x2 matrix to the polygons making up this surface. 
        Returns the flatsurf.geometry.SurfaceMapping from this surface to its image.
        """
        from flatsurf.geometry.mappings import GL2RMapping, IdentityMapping, SurfaceMappingComposition
        mapping = GL2RMapping(self, matrix)
        codomain = mapping.codomain()
        fixed_codomain = convert_to_type(codomain, self.surface_type())
        identity = IdentityMapping(codomain,fixed_codomain)
        return SurfaceMappingComposition(mapping,identity)
        
    def __rmul__(self,matrix):
        r"""
        EXAMPLES::
            sage: from flatsurf.geometry.similarity_surface_generators import InfiniteStaircase
            sage: s=InfiniteStaircase()
            sage: print s
            The infinite staircase
            sage: m=Matrix([[1,2],[0,1]])
            sage: s2=m*s
            sage: print s2.polygon(0)
            Polygon: (0, 0), (1, 0), (3, 1), (2, 1)
        """
        from sage.matrix.matrix import Matrix
        if not isinstance(matrix,Matrix):
            raise NotImplementedError("Only implemented for matrices.")
        if not matrix.dimensions!=(2,2):
            raise NotImplementedError("Only implemented for 2x2 matrices.")
        return self.GL2R_mapping(matrix).codomain()

class HalfDilationSurface_polygons_and_gluings(
        SimilaritySurface_polygons_and_gluings,
        HalfDilationSurface_generic):
    pass
    
class HalfDilationSurface_wrapper(
        SimilaritySurface_wrapper,
        HalfDilationSurface_generic):
    pass

def convert_to_half_dilation_surface(surface):
    r"""
    Returns a dilation surface version of the provided surface.
    """
    if surface.is_finite():
        return HalfDilationSurface_polygons_and_gluings(surface)
    else:
        return HalfDilationSurface_wrapper(surface)

