r"""Mappings between translation surfaces."""

from flatsurf.geometry.polygon import Polygons, wedge_product
from flatsurf.geometry.surface import Surface, Surface_polygons_and_gluings, ExtraLabel
from flatsurf.geometry.similarity_surface import SimilaritySurface
from flatsurf.geometry.translation import TranslationGroup

from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

class SurfaceMapping:
    r"""Abstract class for any mapping between surfaces."""
    
    def __init__(self, domain, codomain):
        self._domain=domain
        self._codomain=codomain
    
    def domain(self):
        r"""
        Return the domain of the mapping.
        """
        return self._domain

    def codomain(self):
        r"""
        Return the range of the mapping.
        """
        return self._codomain

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        raise NotImplementedError

    def pull_vector_back(self,tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        raise NotImplementedError
        
    def __mul__(self,other):
        # Compose SurfaceMappings
        return SurfaceMappingComposition(other,self)
    
    def __rmul__(self,other):
        return SurfaceMappingComposition(self,other)


class FinitelyPerturbedSurface(Surface):
    def __init__(self, surface, polygon_dictionary=None, glue_dictionary=None, base_label=None, ring=None):
        r"""
        
        Warning: No checks are made to make sure the surface is reasonable.
        
        PARAMETERS::
            surface: The surface this is based on.
            polygon_dictionary: A dictionary mapping labels to polygons which will override whatever was on the original surface.
            glue_dictionary: A dictionary mapping edges to edges, which will override whatever was on the original surface. It will automatically be made symmetric.
            base_label: A label representing the base_label on the new surface.
            ring: A new ring containing the vertices.
        """
        self._s=surface
        if polygon_dictionary is None:
            self._pdict={}
        else:
            self._pdict=dict(polygon_dictionary)
        self._gdict={}
        if not glue_dictionary is None:
            for edge1,edge2 in glue_dictionary.iteritems():
                self._gdict[edge1]=edge2
                self._gdict[edge2]=edge1
        if base_label is None:
            self._base_label = surface.base_label()
        else:
            self._base_label = base_label
        if ring is None:
            self._ring = surface.base_ring()
        else:
            self._ring = ring
        self._is_finite = surface.is_finite()
        Surface.__init__(self)

    def base_ring(self):
        return self._ring

    def base_label(self):
        return self._base_label

    def polygon(self, lab):
        p = self._pdict.get(lab)
        if p is None:
            return self._s.polygon(lab)
        return p

    def opposite_edge(self, p, e):
        edge = self._gdict.get((p,e))
        if edge is None:
            return self._s.opposite_edge(p,e)
        return edge

    def is_finite(self):
        return self._is_finite


class BaseLabelChangedSurface(Surface):
    def __init__(self, surface, base_label):
        r"""
        Construct a copy of the provided surface with a new base label.
        """
        self._s=surface
        self._base_label = base_label
        Surface.__init__(self)

    def base_ring(self):
        return self._s.base_ring()

    def base_label(self):
        return self._base_label

    def polygon(self, lab):
        return self._s.polygon(lab)

    def opposite_edge(self, p, e):
        return self._s.opposite_edge(p,e)

    def is_finite(self):
        return self._s.is_finite()

class SurfaceMappingComposition(SurfaceMapping):
    r"""
    Compose two mappings.
    """
    
    def __init__(self, mapping1, mapping2):
        r"""
        Represent the mapping of mapping1 followed by mapping2.
        """
        if mapping1.codomain()!=mapping2.domain():
            raise ValueError("Codomain of mapping1 must be equal to the domain of mapping2")
        self._m1=mapping1
        self._m2=mapping2
        SurfaceMapping.__init__(self, self._m1.domain(), self._m2.codomain())

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        return self._m2.push_vector_forward(self._m1.push_vector_forward(tangent_vector))

    def pull_vector_back(self,tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        return self._m1.pull_vector_back(self._m2.pull_vector_back(tangent_vector))

class IdentityMapping(SurfaceMapping):
    r"""
    Construct an identity map between two `equal' surfaces.
    """
    def __init__(self,domain,codomain):
        SurfaceMapping.__init__(self, domain, codomain)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self._codomain.tangent_vector( \
            tangent_vector.polygon_label(), \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)


    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self._domain.tangent_vector( \
            tangent_vector.polygon_label(), \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)

class MatrixListDeformedSurface(Surface):
    r"""
    Apply a different matrix to each polygon in the surface. 
    Here matrix_function is a python function mapping labels to 2x2 matrices with positive determinant.
    """
    def __init__(self, surface, matrix_function, ring=None):
        self._s=surface
        self._m=matrix_function
        if ring is None:
            self._base_ring = self._s.base_ring()
        else:
            self._base_ring=ring
        self._P=Polygons(self._base_ring)
        Surface.__init__(self)

    def base_ring(self):
        return self._base_ring

    def base_label(self):
        return self._s.base_label()

    def polygon(self, lab):
        p = self._s.polygon(lab)
        edges = [ self._m(lab) * p.edge(e) for e in xrange(p.num_edges())]
        return self._P(edges)

    def opposite_edge(self, p, e):
        return self._s.opposite_edge(p,e)

    def is_finite(self):
        return self._s.is_finite()

class MatrixListDeformedSurfaceMapping(SurfaceMapping):
    r"""
    This mapping applies a possibly different linear matrix to each polygon.
    The matrix is determined by the matrix_function which should be a python
    object.
    """
    def __init__(self, s, matrix_function, ring=None):
        codomain = MatrixListDeformedSurface(s,matrix_function,ring = ring)
        self._m=matrix_function
        SurfaceMapping.__init__(self, s, codomain)

    def push_vector_forward(self,tangent_vector):
        label = tangent_vector.polygon_label()
        m = self._m(label)
        return self.codomain().tangent_vector(
                label, \
                m*tangent_vector.point(), \
                m*tangent_vector.vector())

    def pull_vector_back(self,tangent_vector):
        label = tangent_vector.polygon_label()
        im = ~self._m(label)
        return self.domain().tangent_vector(
                label, \
                im*tangent_vector.point(), \
                im*tangent_vector.vector())

class SimilarityJoinPolygonsMapping(SurfaceMapping):
    r"""
    Return a SurfaceMapping joining two polygons together along the edge provided to the constructor.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import Polygons
        sage: P=Polygons(QQ)
        sage: tri0=P([(1,0),(0,1),(-1,-1)])
        sage: tri1=P([(-1,0),(0,-1),(1,1)])
        sage: gluings=[((0,0),(1,0)),((0,1),(1,1)),((0,2),(1,2))]
        sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: s=TranslationSurface(Surface_polygons_and_gluings([tri0,tri1], gluings))
        sage: from flatsurf.geometry.mappings import *
        sage: m=SimilarityJoinPolygonsMapping(s,0,2)
        sage: s2=m.codomain()
        sage: for label,polygon in s2.label_iterator(polygons=True):
        ...       print "Polygon "+str(label)+" is "+str(polygon)+"."
        Polygon 0 is Polygon: (0, 0), (1, 0), (1, 1), (0, 1).
        sage: for label,edge in s2.edge_iterator():
        ...       print str((label,edge))+" is glued to "+str(s2.opposite_edge(label,edge))+"."
        (0, 0) is glued to (0, 2).
        (0, 1) is glued to (0, 3).
        (0, 2) is glued to (0, 0).
        (0, 3) is glued to (0, 1).
    """
    def __init__(self, s, p1, e1):
        r"""
        Join polygon with label p1 of s to polygon sharing edge e1.
        """
        poly1=s.polygon(p1)
        p2,e2 = s.opposite_edge(p1,e1)
        poly2=s.polygon(p2)
        t=s.edge_transformation(p2, e2)
        dt=t.derivative()
        vs = []
        edge_map={} # Store the pairs for the old edges.
        for i in range(e1):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))
        ne=poly2.num_edges()
        for i in range(1,ne):
            ee=(e2+i)%ne
            edge_map[len(vs)]=(p2,ee)
            vs.append(dt * poly2.edge( ee ))
        for i in range(e1+1, poly1.num_edges()):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))

        inv_edge_map={}
        for key, value in edge_map.iteritems():
            inv_edge_map[value]=(p1,key)

        base_label=s.base_label()
        if base_label==p2:
             base_label=p1
        
        glue_dictionary={}
        for i in range(len(vs)):
            p3,e3 = edge_map[i]
            p4,e4 = s.opposite_edge(p3,e3)
            if p4 == p1 or p4 == p2: 
                glue_dictionary[(p1,i)] = inv_edge_map[(p4,e4)]
            else:
                glue_dictionary[(p1,i)] = (p4,e4)

        s2 = s.__class__(FinitelyPerturbedSurface(
            s, 
            polygon_dictionary={p1: Polygons(s.base_ring())(vs)}, 
            glue_dictionary=glue_dictionary, 
            base_label=base_label, 
            ring = s.base_ring()))

        self._saved_label=p1
        self._removed_label=p2
        self._remove_map = t
        self._remove_map_derivative = dt
        self._glued_edge=e1
        SurfaceMapping.__init__(self, s, s2)

    def removed_label(self):
        r"""
        Return the label that was removed in the joining process.
        """
        return self._removed_label

    def glued_vertices(self):
        r"""
        Return the vertices of the newly glued polygon which bound the diagonal formed by the glue.
        """
        return (self._glued_edge,self._glued_edge+self._domain.polygon(self._removed_label).num_edges())

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._removed_label:
            return self._codomain.tangent_vector( \
                self._saved_label, \
                self._remove_map(tangent_vector.point()), \
                self._remove_map_derivative*tangent_vector.vector(), \
                ring = ring)
        else:
            return self._codomain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

    def pull_vector_back(self,tangent_vector):
        r"""
        Applies the inverse of the mapping to the provided vector.
        """
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._saved_label:
            p=tangent_vector.point()
            v=self._domain.polygon(self._saved_label).vertex(self._glued_edge)
            e=self._domain.polygon(self._saved_label).edge(self._glued_edge)
            from flatsurf.geometry.polygon import wedge_product
            wp = wedge_product(p-v,e)
            if wp > 0:
                # in polygon with the removed label
                return self.domain().tangent_vector( \
                    self._removed_label, \
                    (~ self._remove_map)(tangent_vector.point()), \
                    (~ self._remove_map_derivative)*tangent_vector.vector(), \
                    ring = ring)
            if wp < 0:
                # in polygon with the removed label
                return self.domain().tangent_vector( \
                    self._saved_label, \
                    tangent_vector.point(), \
                    tangent_vector.vector(), \
                    ring = ring)
            # Otherwise wp==0
            w = tangent_vector.vector()
            wp = wedge_product(w,e)
            if wp > 0:
                # in polygon with the removed label
                return self.domain().tangent_vector( \
                    self._removed_label, \
                    (~ self._remove_map)(tangent_vector.point()), \
                    (~ self._remove_map_derivative)*tangent_vector.vector(), \
                    ring = ring)
            return self.domain().tangent_vector( \
                self._saved_label, \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)
        else:
            return self._domain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

class SplitPolygonsMapping(SurfaceMapping):
    r"""
    Class for cutting a polygon along a diagonal.
    
    EXAMPLES::

        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: from flatsurf.geometry.polygon import Polygons
        sage: p = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
        sage: gluings=[((0,i),(0,i+4)) for i in range(4)]
        sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: s=TranslationSurface(Surface_polygons_and_gluings([p], gluings))
        sage: from flatsurf.geometry.mappings import SplitPolygonsMapping
        sage: m = SplitPolygonsMapping(s,0,0,2)
        sage: s2=m.codomain()
        sage: for pair in s2.label_iterator(polygons=True):
        ...       print pair
        (0, Polygon: (0, 0), (1/2*sqrt2 + 1, 1/2*sqrt2), (1/2*sqrt2 + 1, 1/2*sqrt2 + 1), (1, sqrt2 + 1), (0, sqrt2 + 1), (-1/2*sqrt2, 1/2*sqrt2 + 1), (-1/2*sqrt2, 1/2*sqrt2))
        (ExtraLabel(0), Polygon: (0, 0), (-1/2*sqrt2 - 1, -1/2*sqrt2), (-1/2*sqrt2, -1/2*sqrt2))
        sage: for glue in s2.edge_iterator(gluings=True):
        ...       print glue
        ((0, 0), (ExtraLabel(0), 0))
        ((0, 1), (0, 5))
        ((0, 2), (0, 6))
        ((0, 3), (ExtraLabel(0), 1))
        ((0, 4), (ExtraLabel(0), 2))
        ((0, 5), (0, 1))
        ((0, 6), (0, 2))
        ((ExtraLabel(0), 0), (0, 0))
        ((ExtraLabel(0), 1), (0, 3))
        ((ExtraLabel(0), 2), (0, 4))
    """
    
    def __init__(self, s, p, v1, v2, new_label = None):
        r"""
        Split the polygon with label p of surface s along the diagonal joining vertex v1 to vertex v2.
        
        Warning: We do not ensure that new_label is not already in the list of labels unless it is None (as by default).
        """
        if new_label is None:
            new_label = ExtraLabel()
        
        poly=s.polygon(p)
        ne=poly.num_edges()
        if v1<0 or v2<0 or v1>=ne or v2>=ne:
            raise ValueError('Provided vertices out of bounds.')
        if abs(v1-v2)<=1 or abs(v1-v2)>=ne-1:
            raise ValueError('Provided diagonal is not a diagonal.')
        if v2<v1:
            temp=v1
            v1=v2
            v2=temp
            
        newvertices1=[poly.vertex(v2)-poly.vertex(v1)]
        for i in range(v2, v1+ne):
            newvertices1.append(poly.edge(i))
        newpoly1 = Polygons(s.base_ring())(newvertices1)
        
        newvertices2=[poly.vertex(v1)-poly.vertex(v2)]
        for i in range(v1,v2):
            newvertices2.append(poly.edge(i))
        newpoly2 = Polygons(s.base_ring())(newvertices2)
            
        old_to_new_labels={}
        for i in range(ne):
            if i<v1:
                old_to_new_labels[i]=(p,i+ne-v2+1)
            elif i<v2:
                old_to_new_labels[i]=(new_label,i-v1+1)
            else: # i>=v2
                old_to_new_labels[i]=(p,i-v2+1)
        new_to_old_labels={}
        for i,pair in old_to_new_labels.iteritems():
            new_to_old_labels[pair]=i

        glue_dictionary = {(p,0):(new_label,0)}
        for e in range(ne):
            ll,ee = old_to_new_labels[e]
            lll,eee = s.opposite_edge(p,e)
            if lll == p:
                glue_dictionary[(ll,ee)]=old_to_new_labels[eee]
            else:
                glue_dictionary[(ll,ee)]=(lll,eee)
        
        s2 = s.__class__(FinitelyPerturbedSurface(
            s, 
            polygon_dictionary={p: newpoly1, new_label: newpoly2}, 
            glue_dictionary=glue_dictionary, 
            base_label = s.base_label(), 
            ring = s.base_ring()))
        
        self._p=p
        self._v1=v1
        self._v2=v2
        self._new_label=new_label
        from flatsurf.geometry.translation import TranslationGroup
        TG=TranslationGroup(s.base_ring())
        self._tp = TG(-s.polygon(p).vertex(v1))
        self._tnew_label = TG(-s.polygon(p).vertex(v2))
        SurfaceMapping.__init__(self, s, s2)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._p:
            point=tangent_vector.point()
            vertex1=self._domain.polygon(self._p).vertex(self._v1)
            vertex2=self._domain.polygon(self._p).vertex(self._v2)

            wp = wedge_product(vertex2-vertex1,point-vertex1)

            if wp > 0:
                # in new polygon 1
                return self.codomain().tangent_vector( \
                    self._p, \
                    self._tp(tangent_vector.point()), \
                    tangent_vector.vector(), \
                    ring = ring)
            if wp < 0:
                # in new polygon 2
                return self.codomain().tangent_vector( \
                    self._new_label, \
                    self._tnew_label(tangent_vector.point()), \
                    tangent_vector.vector(), \
                    ring = ring)

            # Otherwise wp==0
            w = tangent_vector.vector()
            wp = wedge_product(vertex2-vertex1,w)
            if wp > 0:
                # in new polygon 1
                return self.codomain().tangent_vector( \
                    self._p, \
                    self._tp(tangent_vector.point()), \
                    tangent_vector.vector(), \
                    ring = ring)
            # in new polygon 2
            return self.codomain().tangent_vector( \
                self._new_label, \
                self._tnew_label(tangent_vector.point()), \
                tangent_vector.vector(), \
                ring = ring)
        else:
            # Not in a polygon that was changed. Just copy the data.
            return self._codomain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)


    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._p:
            return self._domain.tangent_vector( \
                self._p, \
                (~ self._tp)(tangent_vector.point()), \
                tangent_vector.vector(), \
                ring = ring)
        elif tangent_vector.polygon_label() == self._new_label:
            return self._domain.tangent_vector( \
                self._p, \
                (~ self._tnew_label)(tangent_vector.point()), \
                tangent_vector.vector(), \
                ring = ring)
        else:
            # Not in a polygon that was changed. Just copy the data.
            return self._domain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

def subdivide_a_polygon(s):
    r"""
    Return a SurfaceMapping which cuts one polygon along a diagonal or None if the surface is triangulated.
    """
    from flatsurf.geometry.polygon import wedge_product
    for l,poly in s.label_iterator(polygons=True):
        n = poly.num_edges() 
        if n>3:
            for i in xrange(n):
                e1=poly.edge(i)
                e2=poly.edge((i+1)%n)
                if wedge_product(e1,e2) != 0:
                    return SplitPolygonsMapping(s,l,i, (i+2)%n)
            raise ValueError("Unable to triangulate polygon with label "+str(l)+\
                ": "+str(poly))
    return None


def triangulation_mapping(s):
    r"""Return a  SurfaceMapping triangulating the provided surface.
    
    EXAMPLES::
        
        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: from flatsurf.geometry.polygon import Polygons
        sage: p = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
        sage: gluings=[((0,i),(0,i+4)) for i in range(4)]
        sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: s=TranslationSurface(Surface_polygons_and_gluings([p], gluings))
        sage: from flatsurf.geometry.mappings import *
        sage: m=triangulation_mapping(s)
        sage: s2=m.codomain()
        sage: for label,polygon in s2.label_iterator(polygons=True):
        ...       print str(polygon)
        Polygon: (0, 0), (-1/2*sqrt2, 1/2*sqrt2 + 1), (-1/2*sqrt2, 1/2*sqrt2)
        Polygon: (0, 0), (1/2*sqrt2, -1/2*sqrt2 - 1), (1/2*sqrt2, 1/2*sqrt2)
        Polygon: (0, 0), (-1/2*sqrt2 - 1, -1/2*sqrt2 - 1), (0, -1)
        Polygon: (0, 0), (-1, -sqrt2 - 1), (1/2*sqrt2, -1/2*sqrt2)
        Polygon: (0, 0), (0, -sqrt2 - 1), (1, 0)
        Polygon: (0, 0), (-1/2*sqrt2 - 1, -1/2*sqrt2), (-1/2*sqrt2, -1/2*sqrt2)
    """
    assert(s.is_finite())
    m=subdivide_a_polygon(s)
    if m is None:
        return None
    s1=m.codomain()
    while True:
        m2=subdivide_a_polygon(s1)
        if m2 is None:
            return m
        s1=m2.codomain()
        m=SurfaceMappingComposition(m,m2)
    
def edge_needs_flip_Linfinity(s, p1, e1):
    r"""
    Check whether the provided edge which bouds two triangles should be flipped
    to get closer to the L-infinity Delaunay decomposition.

    EXAMPLES::

        sage: from flatsurf import *
        sage: t1 = polygons(vertices=[(0,0), (1,0), (1,1)])
        sage: t2 = polygons(vertices=[(0,0), (1,1), (0,1)])
        sage: m = matrix(2, [2,1,1,1])
        sage: t1 = m*t1
        sage: t2 = m*t2
        sage: s = similarity_surfaces([t1,t2], {(0,0):(1,1), (0,1):(1,2), (0,2):(1,0)})

        sage: from flatsurf.geometry.mappings import edge_needs_flip_Linfinity
        sage: edge_needs_flip_Linfinity(s, 0, 0)
        False
        sage: edge_needs_flip_Linfinity(s, 0, 1)
        False
        sage: edge_needs_flip_Linfinity(s, 0, 2)
        True
        sage: edge_needs_flip_Linfinity(s, 1, 0)
        True
        sage: edge_needs_flip_Linfinity(s, 1, 1)
        False
        sage: edge_needs_flip_Linfinity(s, 1, 2)
        False
    """
    p2,e2 = s.opposite_edge(p1,e1)
    poly1 = s.polygon(p1)
    poly2 = s.polygon(p2)
    assert poly1.num_edges() == 3
    assert poly2.num_edges() == 3

    # convexity check of the quadrilateral
    if wedge_product(poly2.edge(e2-1), poly1.edge(e1+1)) <= 0 or \
       wedge_product(poly1.edge(e1-1), poly2.edge(e2+1)) <=0:
        return False

    # compare the norms
    edge1 = poly1.edge(e1)
    edge = poly2.edge(e2-1) + poly1.edge(e1+1)
    n1 = max(abs(edge1[0]), abs(edge1[1]))
    n = max(abs(edge[0]), abs(edge[1]))
    return n < n1

def flip_edge_mapping(s,p1,e1):
    r"""
    Return a mapping whose domain is s which flips the provided edge.
    """
    m1=SimilarityJoinPolygonsMapping(s,p1,e1)
    v1,v2=m1.glued_vertices()
    removed_label = m1.removed_label()
    m2=SplitPolygonsMapping(m1.codomain(), p1, (v1+1)%4, (v1+3)%4, new_label = removed_label)
    return SurfaceMappingComposition(m1,m2)

def one_delaunay_flip_mapping(s):
    r"""
    Returns one delaunay flip, or none if no flips are needed.
    """
    for p,poly in s.label_iterator(polygons=True):
        for e in range(poly.num_edges()):
            if s._edge_needs_flip(p,e):
                return flip_edge_mapping(s,p,e)
    return None

#def edge_needs_join(s,p1,e1):
#    r"""
#    Return if the provided edge which bounds two triangles should be flipped
#    to get closer to the Delaunay decomposition
#    """
#    p2,e2=s.opposite_edge(p1,e1)
#    poly1=s.polygon(p1)
#    poly2=s.polygon(p2)
#    assert poly1.num_edges()==3
#    assert poly2.num_edges()==3
#    from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
#    sim1=similarity_from_vectors(poly1.edge(e1+2),-poly1.edge(e1+1))
#    sim2=similarity_from_vectors(poly2.edge(e2+2),-poly2.edge(e2+1))
#    sim=sim1*sim2
#    return sim[1][0] == 0

def delaunay_triangulation_mapping(s):
    r"""
    Returns a mapping to a Delaunay triangulation or None if the surface already is Delaunay triangulated.
    """
    assert(s.is_finite())
    m=triangulation_mapping(s)
    if m is None:
        s1=s
    else: 
        s1=m.codomain()
    m1=one_delaunay_flip_mapping(s1)
    if m1 is None:
        return m
    if m is None:
        m=m1
    else:
        m=SurfaceMappingComposition(m,m1)
    s1=m1.codomain()
    while True:
        m1=one_delaunay_flip_mapping(s1)
        if m1 is None:
            return m
        s1=m1.codomain()
        m=SurfaceMappingComposition(m,m1)

def delaunay_decomposition_mapping(s):
    r"""
    Returns a mapping to a Delaunay decomposition or possibly None if the surface already is Delaunay.
    """
    m=delaunay_triangulation_mapping(s)
    if m is None:
        s1=s
    else:
        s1=m.codomain()
    edge_vectors=[]
    for p,poly in s1.label_iterator(polygons=True):
        for e in range(poly.num_edges()):
            pp,ee=s1.opposite_edge(p,e)
            if (p<pp or (p==pp and e<ee)) and s1._edge_needs_join(p,e):
                edge_vectors.append( s1.tangent_vector(p,poly.vertex(e),poly.edge(e)) )
    if len(edge_vectors)>0:
        ev=edge_vectors.pop()
        p,e=ev.edge_pointing_along()
        m1=SimilarityJoinPolygonsMapping(s1,p,e)
        s2=m1.codomain()
        while len(edge_vectors)>0:
            ev=edge_vectors.pop()
            ev2=m1.push_vector_forward(ev)
            p,e=ev2.edge_pointing_along()
            mtemp=SimilarityJoinPolygonsMapping(s2,p,e)
            m1=SurfaceMappingComposition(m1,mtemp)
            s2=m1.codomain()
        if m is None:
            return m1
        else:
            return SurfaceMappingComposition(m,m1)
    return m
    
def canonical_first_vertex(polygon):
    r"""
    Return the index of the vertex with smallest y-coordinate.
    If two vertices have the same y-coordinate, then the one with least x-coordinate is returned.
    """
    best=0
    best_pt=polygon.vertex(best)
    for v in range(1,polygon.num_edges()):
        pt=polygon.vertex(v)
        if pt[1]<best_pt[1]:
            best=v
            best_pt=pt
    if best==0:
        if pt[1]==best_pt[1]:
            return v
    return best
   
class CanonicalizePolygonsMapping(SurfaceMapping):
    r"""
    This is a mapping to a surface with the polygon vertices canonically determined.
    A canonical labeling is when the canonocal_first_vertex is the zero vertex.
    """
    def __init__(self, s):
        r"""
        Split the polygon with label p of surface s along the diagonal joining vertex v1 to vertex v2.
        """
        if not s.is_finite():
            raise ValueError("Currently only works with finite surfaces.""")
        ring=s.base_ring()
        T=TranslationGroup(ring)
        P=Polygons(ring)
        cv = {} # dictionary for canonical vertices
        newpolys={} # Polygons for new surfaces
        translations={} # translations bringing the canonical vertex to the origin.
        for l,polygon in s.label_iterator(polygons=True):
            cv[l]=cvcur=canonical_first_vertex(polygon)
            newedges=[]
            for i in range(polygon.num_edges()):
                newedges.append(polygon.edge( (i+cvcur) % polygon.num_edges() ))
            newpolys[l]=P(newedges)
            translations[l]=T( -polygon.vertex(cvcur) )
        newgluing=[]
        for l1,polygon in s.label_iterator(polygons=True):
            for e1 in range(polygon.num_edges()):
                l2,e2=s.opposite_edge(l1,e1)
                ee1= (e1-cv[l1]+polygon.num_edges())%polygon.num_edges()
                polygon2=s.polygon(l2)
                ee2= (e2-cv[l2]+polygon2.num_edges())%polygon2.num_edges()
                newgluing.append( ( (l1,ee1),(l2,ee2) ) )

        s2=s.__class__(Surface_polygons_and_gluings(newpolys,newgluing))
        
        self._cv=cv
        self._translations=translations

        SurfaceMapping.__init__(self, s, s2)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        l=tangent_vector.polygon_label()
        return self.codomain().tangent_vector(l, \
            self._translations[l](tangent_vector.point()), \
            tangent_vector.vector(), \
            ring = ring)

    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        l=tangent_vector.polygon_label()
        return self.domain().tangent_vector(l, \
            (~self._translations[l])(tangent_vector.point()), \
            tangent_vector.vector(), \
            ring = ring)

class ReindexMapping(SurfaceMapping):
    r"""
    Apply a dictionary to relabel the polygons.
    """
    def __init__(self,s,relabler,new_base_label=None):
        r"""
        The parameters should be a surface and a dictionary which takes as input a label and produces a new label.
        """
        if not s.is_finite():
            raise ValueError("Currently only works with finite surfaces.""")
        f = {} # map for labels going forward.
        b = {} # map for labels going backward.
        for l in s.label_iterator():
            if relabler.has_key(l):
                l2=relabler[l]
                f[l]=l2
                if b.has_key(l2):
                    raise ValueError("Provided dictionary has two keys mapping to the same value. Or you are mapping to a label you didn't change.")
                b[l2]=l
            else:
                # If no key then don't change the label
                f[l]=l
                if b.has_key(l):
                    raise ValueError("Provided dictionary has two keys mapping to the same value. Or you are mapping to a label you didn't change.")
                b[l]=l

        self._f=f
        self._b=b
        
        
        SurfaceMapping.__init__(self, s, s.__class__(ReindexMapping.ReindexedSurface(s,self,new_base_label)))

    class ReindexedSurface(Surface):
        def __init__(self, s, reindexmapping,new_base_label=None):
            r"""
            Represents a reindexed similarity surface.
            """
            self._s=s
            self._r=reindexmapping
            if new_base_label is None:
                self._base_label=self._r._f[self._s.base_label()]
            else:
                self._base_label=new_base_label
            Surface.__init__(self)
        
        def base_ring(self):
            return self._s.base_ring()
        
        def base_label(self):
            return self._base_label
        
        def polygon(self, lab):
            return self._s.polygon(self._r._b[lab])
        
        def opposite_edge(self, p, e):
            p_back = self._r._b[p]
            pp_back,ee = self._s.opposite_edge(p_back,e)
            pp = self._r._f[pp_back]
            return (pp,ee)
        
        def is_finite(self):
            return self._s.is_finite()
    
    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        # There is no change- we just move it to the new surface.
        ring = tangent_vector.bundle().base_ring()
        return self.codomain().tangent_vector( \
            tangent_vector.polygon_label(), \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)

    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self.domain().tangent_vector( \
            tangent_vector.polygon_label(), \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)

def my_sgn(val):
    if val>0:
        return 1
    elif val<0:
        return -1
    else:
        return 0

def polygon_compare(poly1,poly2):
    r"""
    Compare two polygons first by area, then by number of sides,
    then by lexigraphical ording on edge vectors."""
    # This should not be used is broken!!
    #from sage.functions.generalized import sgn
    res = my_sgn(-poly1.area()+poly2.area())
    if res!=0:
        return res
    res = my_sgn(poly1.num_edges()-poly2.num_edges())
    if res!=0:
        return res
    ne=poly1.num_edges()
    for i in range(0,ne-1):
        edge_diff = poly1.edge(i) - poly2.edge(i)
        res = my_sgn(edge_diff[0])
        if res!=0:
            return res
        res = my_sgn(edge_diff[1])
        if res!=0:
            return res
    return 0
    
def translation_surface_cmp(s1, s2):
    r"""
    Compare two finite surfaces. 
    The surfaces will be considered equal if and only if there is a translation automorphism
    respecting the polygons and the base_labels.
    """
    if not s1.is_finite() or not s2.is_finite():
        raise NotImplementedError
    lw1=s1.walker()
    lw2=s2.walker()
    from itertools import izip_longest
    for p1,p2 in izip_longest(lw1.polygon_iterator(), lw2.polygon_iterator()):
        if p1 is None:
            # s2 has more polygons
            return -1
        if p2 is None:
            # s1 has more polygons
            return 1
        ret = polygon_compare(p1,p2)
        if ret != 0:
            return ret
    # Polygons are identical. Compare edge gluings.
    for pair1,pair2 in izip_longest(lw1.edge_iterator(), lw2.edge_iterator()):
        l1,e1 = s1.opposite_edge(pair1)
        l2,e2 = s2.opposite_edge(pair2)
        num1 = lw1.label_to_number(l1)
        num2 = lw2.label_to_number(l2)
        ret = cmp(num1,num2)
        if ret!=0:
            return ret
        ret = cmp(e1,e2)
        if ret!=0:
            return ret
    return 0

def canonicalize_translation_surface_mapping(s):
    r"""
    Return the translation surface in a canonical form.
    
    EXAMPLES::

        sage: from flatsurf.geometry.polygon import Polygons
        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
        sage: square1 = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
        sage: square2 = Polygons(K)([(sqrt2/2, sqrt2/2),(-sqrt2/2, sqrt2/2),(-sqrt2/2, -sqrt2/2),(sqrt2/2, -sqrt2/2)])
        sage: gluings=[((1,i),(0, (2*i+4)%8 )) for i in range(4)]
        sage: for i in range(4):
        ...       gluings.append( ((2,i), (0, (2*i+1+4)%8 )) )
        sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: s=TranslationSurface(Surface_polygons_and_gluings([octagon,square1,square2], gluings))
        sage: from flatsurf.geometry.mappings import *
        sage: mat=Matrix([[1,2+sqrt2],[0,1]])
        sage: from flatsurf.geometry.half_dilation_surface import GL2RMapping
        sage: m1=GL2RMapping(s,mat)
        sage: m2=canonicalize_translation_surface_mapping(m1.codomain())
        sage: m=m2*m1
        sage: translation_surface_cmp(m.domain(),m.codomain())==0
        True
        sage: s=m.domain()
        sage: v=s.tangent_vector(0,(0,0),(1,1))
        sage: w=m.push_vector_forward(v)
        sage: print(w)
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (sqrt2 + 3, 1)
    """
    from flatsurf.geometry.translation_surface import TranslationSurface
    if not s.is_finite():
        raise NotImplementedError
    if not isinstance(s,TranslationSurface):
        raise ValueError("Only defined for TranslationSurfaces")
    m1=delaunay_decomposition_mapping(s)
    s2=m1.codomain()
    m2=CanonicalizePolygonsMapping(s2)
    m=SurfaceMappingComposition(m1,m2)
    s2=m.codomain()
    it = s2.label_iterator()
    min_label = it.next()
    smin=TranslationSurface(BaseLabelChangedSurface(s2,min_label))
    for test_label in it:
        stest = TranslationSurface(BaseLabelChangedSurface(s2,test_label))
        c = translation_surface_cmp(smin,stest)
        if c>0:
            min_label = test_label
            smin = stest
    lw=smin.walker()
    lw.find_all_labels()
    m3=ReindexMapping(s2,lw.label_dictionary(),0)
    return SurfaceMappingComposition(m,m3)
    
