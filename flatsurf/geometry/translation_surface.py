r"""
Translation Surfaces.
"""
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.surface import Surface
from flatsurf.geometry.half_translation_surface import HalfTranslationSurface 
from flatsurf.geometry.dilation_surface import DilationSurface

from sage.matrix.constructor import matrix, identity_matrix

class TranslationSurface(HalfTranslationSurface, DilationSurface):
    r"""
    A surface with a flat metric and conical singularities whose cone angles are a multiple of pi.
    """
    
    def minimal_translation_cover(self):
        return self

    def _check_edge_matrix(self):
        r"""
        Check the compatibility condition
        """
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if self.is_finite():
            for lab in self.label_iterator():
                p = self.polygon(lab)
                for e in xrange(p.num_edges()):
                    if not SimilaritySurface.edge_matrix(self,lab,e).is_one():
                        raise ValueError("gluings of (%s,%s) is not through translation"%(lab,e))
        else:
            count = 0
            for lab in self.label_iterator():
                p = self.polygon(lab)
                for e in xrange(p.num_edges()):
                    if notSimilaritySurface.edge_matrix(self,lab,e).is_one():
                        raise ValueError("gluings of (%s,%s) is not through translation"%(lab,e))
                count  = count+1
                if count >= 10:
                    break
    
    def edge_matrix(self, p, e=None):
        if e is None:
            p,e = p
        if e < 0 or e >= self.polygon(p).num_edges():
            raise ValueError
        return identity_matrix(self.base_ring(),2)

    def stratum(self):
        r"""
        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: sfg.translation_surfaces.octagon_and_squares().stratum()
            H(4)
        """
        from sage.dynamics.flat_surfaces.all import AbelianStratum
        from sage.rings.integer_ring import ZZ
        return AbelianStratum([ZZ(a-1) for a in self.angles()])

    #angles = ConeSurface_generic.angles

    def canonicalize_mapping(self):
        r"""
        Return a SurfaceMapping canonicalizing this translation surface.
        """
        from flatsurf.geometry.mappings import canonicalize_translation_surface_mapping, IdentityMapping
        return canonicalize_translation_surface_mapping(self)
        
    def canonicalize(self):
        r"""
        Return a canonical version of this translation surface.
        
        EXAMPLES:

        We will check if an element lies in the Veech group::

            sage: from flatsurf.geometry.polygon import Polygons
            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
            sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
            sage: square1 = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
            sage: square2 = Polygons(K)([(sqrt2/2, sqrt2/2),(-sqrt2/2, sqrt2/2),(-sqrt2/2, -sqrt2/2),(sqrt2/2, -sqrt2/2)])
            sage: gluings=[((1,i),(0, (2*i+4)%8 )) for i in range(4)]
            sage: for i in range(4):
            ...       gluings.append( ((2,i), (0, (2*i+5)%8 )) )
            sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
            sage: from flatsurf.geometry.translation_surface import TranslationSurface
            sage: s=TranslationSurface(Surface_polygons_and_gluings([octagon,square1,square2],gluings))
            sage: print s
            TranslationSurface built from 3 polygons
            sage: mat=Matrix([[1,2+sqrt2],[0,1]])
            sage: s.canonicalize()==(mat*s).canonicalize()
            True
        """
        return self.canonicalize_mapping().codomain()

class MinimalTranslationCover(Surface):
    r"""
    We label copy by cartesian product (polygon from bot, matrix).
    """
    def __init__(self, similarity_surface):
        self._ss = similarity_surface

    def is_finite(self):
        if not self._ss.is_finite():
            return False
        from flatsurf.geometry.rational_cone_surface import RationalConeSurface
        return isinstance(self._ss, RationalConeSurface)

    def base_ring(self):
        return self._ss.base_ring()

    def base_label(self):
        from sage.matrix.constructor import identity_matrix
        I = identity_matrix(self.base_ring(),2)
        I.set_immutable()
        return (self._ss.base_label(), I)

    def polygon(self, lab):
        return lab[1] * self._ss.polygon(lab[0])

    def opposite_edge(self, p, e):
        pp,m = p  # this is the polygon m * ss.polygon(p)
        p2,e2 = self._ss.opposite_edge(pp,e)
        me = self._ss.edge_matrix(pp,e)
        mm = ~me * m
        mm.set_immutable()
        return ((p2,mm),e2)

class AbstractOrigami(Surface):
    r'''Abstract base class for origamis.
    Realization needs just to define a _domain and four cardinal directions.
    '''

    def up(self, label):
        raise NotImplementedError

    def down(self, label):
        raise NotImplementedError

    def right(self, label):
        raise NotImplementedError

    def left(self, label):
        raise NotImplementedError

    def _repr_(self):
        return "Some AbstractOrigami"

    def is_finite(self):
        return self._domain.is_finite()

    def num_polygons(self):
        r"""
        Returns the number of polygons.
        """
        return self._domain.cardinality()

    def base_label(self):
        return self._domain.an_element()

    def polygon_labels(self):
        return self._domain

    def polygon(self, lab):
        if lab not in self._domain:
            #Updated to print a possibly useful error message
            raise ValueError("Label "+str(lab)+" is not in the domain")
        from flatsurf.geometry.polygon import polygons
        return polygons.square()

    @cached_method
    def base_ring(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: o.base_ring()
            Rational Field
        """
        from sage.rings.rational_field import QQ
        return QQ

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e==0:
            return self.down(p),2
        if e==1:
            return self.right(p),3
        if e==2:
            return self.up(p),0
        if e==3:
            return self.left(p),1
        raise ValueError
        
        return self._perms[e](p), (e+2)%4


class Origami(AbstractOrigami):
    def __init__(self, r, u, rr=None, uu=None, domain=None):
        if domain is None:
            self._domain = r.parent().domain()
        else:
            self._domain = domain

        self._r = r
        self._u = u
        if rr is None:
            rr = ~r
        else:
            for a in self._domain.some_elements():
                if r(rr(a)) != a:
                    raise ValueError("r o rr is not identity on %s"%a)
                if rr(r(a)) != a:
                    raise ValueError("rr o r is not identity on %s"%a)
        if uu is None:
            uu = ~u
        else:
            for a in self._domain.some_elements():
                if u(uu(a)) != a:
                    raise ValueError("u o uu is not identity on %s"%a)
                if uu(u(a)) != a:
                    raise ValueError("uu o u is not identity on %s"%a)

        self._perms = [uu,r,u,rr] # down,right,up,left

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e < 0 or e > 3:
            raise ValueError
        return self._perms[e](p), (e+2)%4

    def up(self, label):
        return self.opposite_edge(label,2)[0]

    def down(self, label):
        return self.opposite_edge(label,0)[0]

    def right(self, label):
        return self.opposite_edge(label,1)[0]

    def left(self, label):
        return self.opposite_edge(label,3)[0]

    def _repr_(self):
        return "Origami defined by r=%s and u=%s"%(self._r,self._u)


