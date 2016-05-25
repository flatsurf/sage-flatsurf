from flatsurf.geometry.surface import Surface
from flatsurf.geometry.similarity_surface import SimilaritySurface
from flatsurf.geometry.mappings import SurfaceMapping, IdentityMapping, SurfaceMappingComposition
from flatsurf.geometry.polygon import Polygons

class HalfDilationSurface(SimilaritySurface):
    r"""
    A generic half translation surface
    """
    def GL2R_mapping(self, matrix):
        r"""
        Apply a 2x2 matrix to the polygons making up this surface. 
        Returns the flatsurf.geometry.SurfaceMapping from this surface to its image.
        """
        return GL2RMapping(self, matrix)
        
    def __rmul__(self,matrix):
        r"""
        EXAMPLES::
            sage: from flatsurf.geometry.similarity_surface_generators import infinite_staircase
            sage: s=infinite_staircase()
            sage: print s.underlying_surface()
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

class GL2RImageSurface(Surface):
    def __init__(self, surface, m, ring=None):
        self._s=surface
        if m.determinant()<=0:
            raise ValueError("Currently only works with matrices of positive determinant.""")
        self._m=m
        if ring is None:
            if m.base_ring() == self._s.base_ring():
                self._base_ring = self._s.base_ring()
            else:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                self._base_ring = cm.common_parent(m.base_ring(), self._s.base_ring())
        else:
            self._base_ring=ring
        self._P=Polygons(self._base_ring)

    def base_ring(self):
        return self._base_ring

    def base_label(self):
        return self._s.base_label()

    def polygon(self, lab):
        p = self._s.polygon(lab)
        edges = [ self._m * p.edge(e) for e in xrange(p.num_edges())]
        return self._P(edges)

    def opposite_edge(self, p, e):
        return self._s.opposite_edge(p,e)

    def is_finite(self):
        return self._s.is_finite()
        

class GL2RMapping(SurfaceMapping):
    r"""
    This class pushes a surface forward under a matrix. 
    
    EXAMPLE::
        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: from flatsurf.geometry.polygon import Polygons
        sage: p = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
        sage: gluings=[((0,i),(0,i+4)) for i in range(4)]
        sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: s=TranslationSurface(Surface_polygons_and_gluings([p], gluings))
        sage: from flatsurf.geometry.half_dilation_surface import GL2RMapping
        sage: mat=Matrix([[2,1],[1,1]])
        sage: m=GL2RMapping(s,mat)
    """
    def __init__(self, s, m, ring=None):
        r"""
        Hit the surface s with the 2x2 matrix m which should have positive determinant.
        """
        codomain = s.__class__(GL2RImageSurface(s,m,ring = ring))
        self._m=m
        self._im=~m
        SurfaceMapping.__init__(self, s, codomain)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        return self.codomain().tangent_vector(
                tangent_vector.polygon_label(), \
                self._m*tangent_vector.point(), \
                self._m*tangent_vector.vector())

    def pull_vector_back(self,tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        return self.domain().tangent_vector(
                tangent_vector.polygon_label(), \
                self._im*tangent_vector.point(), \
                self._im*tangent_vector.vector())
