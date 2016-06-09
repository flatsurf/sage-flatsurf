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

            sage: from flatsurf import *
            sage: s=translation_surfaces.infinite_staircase()
            sage: print s.underlying_surface()
            The infinite staircase
            sage: m=Matrix([[1,2],[0,1]])
            sage: s2=m*s
            sage: TestSuite(s2).run(skip='_test_pickling')
            sage: print s2.polygon(0)
            Polygon: (0, 0), (1, 0), (3, 1), (2, 1)
        """
        from sage.matrix.matrix import Matrix
        if not isinstance(matrix,Matrix):
            raise NotImplementedError("Only implemented for matrices.")
        if not matrix.dimensions!=(2,2):
            raise NotImplementedError("Only implemented for 2x2 matrices.")
        return self.GL2R_mapping(matrix).codomain()

    def _edge_needs_flip_Linfinity(self, p1, e1):
        r"""
        Check whether the provided edge which bouds two triangles should be flipped
        to get closer to the L-infinity Delaunay decomposition.
        """
        p2,e2 = self.opposite_edge(p1,e1)
        poly1 = self.polygon(p1)
        poly2 = self.polygon(p2)
        if poly1.num_edges()!=3 or poly2.num_edges()!=3:
            raise ValueError("Edge must be adjacent to two triangles.")

        sim = self.edge_transformation(l2,e2)
        m = sim.derivative()
        # m is the matrix carrying polygon l2 to polygon l1 along the edge.
        
        # convexity check of the quadrilateral
        if wedge_product(m*poly2.edge(e2-1), poly1.edge(e1+1)) <= 0 or \
           wedge_product(poly1.edge(e1-1), m*poly2.edge(e2+1)) <=0:
            return False

        # compare the norms
        edge1 = poly1.edge(e1)
        edge = m*poly2.edge(e2-1) + poly1.edge(e1+1)
        n1 = max(abs(edge1[0]), abs(edge1[1]))
        n = max(abs(edge[0]), abs(edge[1]))
        return n < n1

    def l_infinity_delaunay_triangulation(self, triangulated=False, in_place=False, limit=None, direction=None):
        r"""
        Returns a L-infinity Delaunay triangulation of a surface, or make some
        triangle flips to get closer to the Delaunay decomposition.
        
        Parameters
        ----------
        triangulated : boolean
            If true, the algorithm assumes the surface is already triangulated. It
            does this without verification.
        in_place : boolean
            If true, the triangulating and the triangle flips are done in place.
            Otherwise, a mutable copy of the surface is made.
        limit : None or Integer
            If None, this will return a Delaunay triangulation. If limit
            is an integer 1 or larger, then at most limit many diagonal flips 
            will be done.
        direction : None or Vector with two entries in the base field
            Used to determine labels when a pair of triangles is flipped. Each triangle
            has a unique separatrix which points in the provided direction or its 
            negation. As such a vector determines a sign for each triangle.
            A pair of adjacent triangles have opposite signs. Labels are chosen
            so that this sign is preserved (as a function of labels).

        EXAMPLES::

            sage: from flatsurf import *
            sage: s0=translation_surfaces.veech_double_n_gon(3)
            sage: field=s0.base_ring()
            sage: a=field.gen()
            sage: from sage.matrix.constructor import Matrix
            sage: m=Matrix([[1,2/a],[0,1]])
            sage: s=(m**5)*s0
            sage: s=s.l_infinity_delaunay_triangulation()
            sage: TestSuite(s).run()
        """
        if not self.is_finite() and limit is None:
            raise NotImplementedError("Not implemented for infinite surfaces unless limit is set")
        if triangulated:
            if in_place:
                s=self
            else:
                from flatsurf.geometry.surface import Surface_dict
                s=self.__class__(Surface_dict(surface=self,mutable=True))
        else:
            from flatsurf.geometry.surface import Surface_list
            s=self.__class__(Surface_list(surface=self.triangulate(in_place=in_place),mutable=True))
        loop=True
        if direction is None:
            base_ring = self.base_ring()
            direction = self.vector_space()( (base_ring.zero(), base_ring.one()) )
        else:
            assert not direction.is_zero()
        count=0
        while loop:
            loop=False
            for (l1,e1),(l2,e2) in s.edge_iterator(gluings=True):
                if (l1<l2 or (l1==l2 and e1<=e2)) and s._edge_needs_flip(l1,e1):
                    s.triangle_flip(l1, e1, in_place=True, direction=direction)
                    count += 1
                    if not limit is None and count>=limit:
                        return s
                    loop=True
                    break
        return s

class GL2RImageSurface(Surface):
    def __init__(self, surface, m, ring=None):

        if surface.is_mutable():
            if surface.is_finite():
                self._s=surface.copy()
            else:
                raise ValueError("Can not apply matrix to mutable infinite surface.")
        else:
            self._s=surface

        if m.determinant()<=0:
            raise ValueError("Currently only works with matrices of positive determinant.")
        self._m=m

        if ring is None:
            if m.base_ring() == self._s.base_ring():
                base_ring = self._s.base_ring()
            else:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                base_ring = cm.common_parent(m.base_ring(), self._s.base_ring())
        else:
            base_ring=ring

        self._P=Polygons(base_ring)

        Surface.__init__(self, base_ring, self._s.base_label(), finite=self._s.is_finite())

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

        sage: from flatsurf import *
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.half_dilation_surface import GL2RMapping
        sage: mat=Matrix([[2,1],[1,1]])
        sage: m=GL2RMapping(s,mat)
        sage: TestSuite(m.codomain()).run()
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
