from __future__ import print_function, absolute_import
from six.moves import range, filter, map

from .surface import Surface
from .half_dilation_surface import HalfDilationSurface
from .rational_cone_surface import RationalConeSurface

class HalfTranslationSurface(HalfDilationSurface, RationalConeSurface):
    r"""
    A half translation surface has gluings between polygons whose monodromy is +I or -I.
    """
    
    def _test_edge_matrix(self, **options):
        r"""
        Check the compatibility condition
        """
        tester = self._tester(**options)
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if self.is_finite():
            it = self.label_iterator()
        else:
            from itertools import islice
            it = islice(self.label_iterator(), 30)

        for lab in it:
            p = self.polygon(lab)
            for e in range(p.num_edges()):
                # Warning: check the matrices computed from the edges,
                # rather the ones overriden by TranslationSurface.
                tester.assertTrue(SimilaritySurface.edge_matrix(self,lab,e).is_one() or \
                    (-SimilaritySurface.edge_matrix(self,lab,e)).is_one(), \
                    "edge_matrix of edge "+str((lab,e))+" is not a translation or rotation by pi.")

    
# This was all implemented in HalfDilationSurface now.
#    
#    def GL2R_mapping(self, matrix):
#        r"""
#        Apply a 2x2 matrix to the polygons making up this surface. 
#        Returns the flatsurf.geometry.SurfaceMapping from this surface to its image.
#        """
#        from flatsurf.geometry.mappings import GL2RMapping, IdentityMapping, SurfaceMappingComposition
#        mapping = GL2RMapping(self, matrix)
#        codomain = mapping.codomain()
#        fixed_codomain = convert_to_type(codomain, self.surface_type())
#        identity = IdentityMapping(codomain,fixed_codomain)
#        return SurfaceMappingComposition(mapping,identity)
#        
#    def __rmul__(self,matrix):
#        r"""
#        EXAMPLES::
#            sage: from flatsurf.geometry.similarity_surface_generators import InfiniteStaircase
#            sage: s=InfiniteStaircase()
#            sage: print s
#            The infinite staircase
#            sage: m=Matrix([[1,2],[0,1]])
#            sage: s2=m*s
#            sage: print s2.polygon(0)
#            Polygon: (0, 0), (1, 0), (3, 1), (2, 1)
#        """
#        from sage.matrix.matrix import Matrix
#        if not isinstance(matrix,Matrix):
#            raise NotImplementedError("Only implemented for matrices.")
#        if not matrix.dimensions!=(2,2):
#            raise NotImplementedError("Only implemented for 2x2 matrices.")
#        return self.GL2R_mapping(matrix).codomain()

