from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

ZZ_1 = ZZ(1)
ZZ_2 = ZZ(2)

from flatsurf.geometry.surface import Surface, Surface_list
from flatsurf.geometry.translation_surface import TranslationSurface

def flipper_nf_to_sage(K, name='a'):
    r"""
    Convert a flipper number field into a Sage number field

    .. NOTE::

        Currently, the code is not careful at all with root isolation.

    EXAMPLES::

        sage: import flipper  # optional - flipper
        sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_to_sage
        sage: p = flipper.kernel.Polynomial([-2r] + [0r]*5 + [1r]) # optional - flipper
        sage: r1,r2 = p.real_roots()                               # optional - flipper
        sage: K = flipper.kernel.NumberField(r1)                   # optional - flipper
        sage: K_sage = flipper_nf_to_sage(K)                       # optional - flipper
        sage: K_sage                                               # optional - flipper
        Number Field in a with defining polynomial x^6 - 2
        sage: AA(K_sage.gen())                                     # optional - flipper
        -1.122462048309373?
    """
    from sage.rings.number_field.number_field import NumberField
    from sage.rings.all import QQ,RIF,AA

    r = K.lmbda.interval_approximation()
    l = r.lower * ZZ(10)**(-r.precision)
    u = r.upper * ZZ(10)**(-r.precision)

    p = QQ['x'](K.polynomial.coefficients)
    s = AA.polynomial_root(p, RIF(l,u))
    return NumberField(p, name, embedding=s)

def flipper_nf_element_to_sage(x):
    r"""
    Convert a flipper number field element into Sage

    EXAMPLES::

        sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_element_to_sage
        sage: import flipper                               # optional - flipper
        sage: T = flipper.load('SB_6')                     # optional - flipper
        sage: h = T.mapping_class('s_0S_1S_2s_3s_4s_3S_5') # optional - flipper
        sage: flipper_nf_element_to_sage(h.dilatation())   # optional - flipper
        a
        sage: AA(_)                                        # optional - flipper
        6.45052513748511?
    """
    return flipper_nf_to_sage(x.number_field)(x.linear_combination)

class InfiniteStaircase(Surface):
    r"""
    The infinite staircase.

     ...
     +--+--+
     |  |  |
     +--+--+--+
        |  |  |
        +--+--+--+
           |  |  |
           +--+--+--+
              |  |  |
              +--+--+
                  ...
    """
    def __init__(self):
        Surface.__init__(self)
    
    def _repr_(self):
        r"""
        String representation.
        """
        return "The infinite staircase"

    def base_ring(self):
        r"""
        Return the rational field.
        """
        from sage.rings.rational_field import QQ
        return QQ

    def base_label(self):
        return ZZ.zero()

    def polygon(self, lab):
        r"""
        Return the polygon labeled by ``lab``.
        """
        if lab not in self.polygon_labels():
            raise ValueError("lab (=%s) not a valid label"%lab)
        from flatsurf.geometry.polygon import polygons
        return polygons.square()

    def polygon_labels(self):
        r"""
        The set of labels used for the polygons.
        """
        return ZZ

    def opposite_edge(self, p, e):
        r"""
        Return the pair ``(pp,ee)`` to which the edge ``(p,e)`` is glued to.
        """
        if (p+e) % 2:
            return p+1,(e+2)%4
        else:
            return p-1,(e+2)%4

    def is_finite(self):
        return False

def infinite_staircase():
    return TranslationSurface(InfiniteStaircase())


class EInfinitySurface(Surface):
    r"""
    The surface based on the $E_\infinity$ graph.

     The biparite graph is shown below, with edges numbered:

      0   1   2  -2   3  -3   4  -4 
    *---o---*---o---*---o---*---o---*...
            |
            |-1
            o

    Here, black vertices are colored *, and white o. 
    Black nodes represent vertical cylinders and white nodes
    represent horizontal cylinders.
    """
    def __init__(self,lambda_squared=None, field=None):
        TranslationSurface_generic.__init__(self)
        if lambda_squared==None:
            from sage.rings.number_field.number_field import NumberField
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R=PolynomialRing(ZZ,'x')
            x = R.gen()
            from sage.rings.qqbar import AA
            self._field=NumberField(x**3-ZZ(5)*x**2+ZZ(4)*x-ZZ(1), 'r', embedding=AA(ZZ(4)))
            self._l=self._field.gen()
        else:
            if field is None:
                self._l=lambda_squared
                self._field=lambda_squared.parent()
            else:
                self._field=field
                self._l=field(lambda_squared)
        Surface.__init__(self)

    def _repr_(self):
        r"""
        String representation.
        """
        return "The E-infinity surface"

    def base_ring(self):
        r"""
        Return the rational field.
        """
        return self._field

    def base_label(self):
        return ZZ.zero()

    @cached_method
    def get_white(self,n):
        r"""Get the weight of the white endpoint of edge n."""
        l=self._l
        if n==0 or n==1:
            return l
        if n==-1:
            return l-1
        if n==2:
            return 1-3*l+l**2
        if n>2:
            x=self.get_white(n-1)
            y=self.get_black(n)
            return l*y-x
        return self.get_white(-n)

    @cached_method
    def get_black(self,n):
        r"""Get the weight of the black endpoint of edge n."""
        l=self._l
        if n==0:
            return self._field(1)
        if n==1 or n==-1 or n==2:
            return l-1
        if n>2:
            x=self.get_black(n-1)
            y=self.get_white(n-1)
            return y-x
        return self.get_black(1-n)

    def polygon(self, lab):
        r"""
        Return the polygon labeled by ``lab``.
        """
        if lab not in self.polygon_labels():
            raise ValueError("lab (=%s) not a valid label"%lab)
        from flatsurf.geometry.polygon import rectangle
        return rectangle(2*self.get_black(lab),self.get_white(lab))

    def polygon_labels(self):
        r"""
        The set of labels used for the polygons.
        """
        return ZZ

    def opposite_edge(self, p, e):
        r"""
        Return the pair ``(pp,ee)`` to which the edge ``(p,e)`` is glued to.
        """
        if p==0:
            if e==0:
                return (0,2)
            if e==1:
                return (1,3)
            if e==2:
                return (0,0)
            if e==3:
                return (1,1)
        if p==1:
            if e==0:
                return (-1,2)
            if e==1:
                return (0,3)
            if e==2:
                return (2,0)
            if e==3:
                return (0,1)
        if p==-1:
            if e==0:
                return (2,2)
            if e==1:
                return (-1,3)
            if e==2:
                return (1,0)
            if e==3:
                return (-1,1)
        if p==2:
            if e==0:
                return (1,2)
            if e==1:
                return (-2,3)
            if e==2:
                return (-1,0)
            if e==3:
                return (-2,1)
        if p>2:
            if e%2:
                return -p,(e+2)%4
            else:
                return 1-p,(e+2)%4
        else:
            if e%2:
                return -p,(e+2)%4
            else:
                return 1-p,(e+2)%4

    def is_finite(self):
        return False
        
def e_infinity_surface(lambda_squared=None, field=None):
    r"""
    The translation surface based on the $E_\infinity$ graph.

    The biparite graph is shown below, with edges numbered:

      0   1   2  -2   3  -3   4  -4 
    *---o---*---o---*---o---*---o---*...
            |
            |-1
            o

    Here, black vertices are colored *, and white o. 
    Black nodes represent vertical cylinders and white nodes
    represent horizontal cylinders.
    """
    return TranslationSurface(EInfinitySurface(lambda_squared, field))

class TFractalSurface(Surface):
    r"""
    The TFractal surface.

    The TFractal surface is a translation surface of finite area built from
    infinitely many polygons. The basic building block is the following polygon

     w/r    w     w/r
    +---+------+---+
    | 1 |   2  | 3 | h2
    +---+------+---+
        |   0  | h1
        +------+
            w

    where ``w``, ``h1``, ``h2``, ``r`` are some positive numbers. Default values
    are ``w=h1=h2=1`` and ``r=2``.

    .. TODO::

        In that surface, the linear flow can be computed more efficiently using
        only one affine interval exchange transformation with 5 intervals. But
        the underlying geometric construction is not a covering.

        Warning: we can not play at the same time with tuples and element of a
        cartesian product (see Sage trac ticket #19555)
    """
    def __init__(self, w=ZZ_1, r=ZZ_2, h1=ZZ_1, h2=ZZ_1):
        from sage.structure.sequence import Sequence
        from sage.combinat.words.words import Words

        self._field = Sequence([w,r,h1,h2]).universe()
        if not self._field.is_field():
            self._field = self._field.fraction_field()
        self._w = self._field(w)
        self._r = self._field(r)
        self._h1 = self._field(h1)
        self._h2 = self._field(h2)
        self._words = Words('LR', finite=True, infinite=False)
        self._wL = self._words('L')
        self._wR = self._words('R')
        Surface.__init__(self)

    def _repr_(self):
        return "The T-fractal surface with parameters w=%s, r=%s, h1=%s, h2=%s"%(
                self._w, self._r, self._h1, self._h2)

    def base_label(self):
        return ZZ.zero()

    def base_ring(self):
        return self._field

    @cached_method
    def polygon_labels(self):
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product([self._words, FiniteEnumeratedSet([0,1,2,3])])

    def opposite_edge(self, p, e):
        r"""

        Labeling of polygons

         wl,0             wr,0
        +-----+---------+------+
        |     |         |      |
        | w,1 |   w,2   |  w,3 |
        |     |         |      |
        +-----+---------+------+
              |         |
              |   w,0   |
              |         |
              +---------+
                   w

        and we always have: bot->0, right->1, top->2, left->3

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: T = sfg.tfractal_surface()
            sage: W = T.underlying_surface()._words
            sage: w = W('LLRLRL')
            sage: print T.opposite_edge((w,0),0)
            ((word: LLRLR, 1), 2)
            sage: print T.opposite_edge((w,0),1)
            ((word: LLRLRL, 0), 3)
            sage: print T.opposite_edge((w,0),2)
            ((word: LLRLRL, 2), 0)
            sage: print T.opposite_edge((w,0),3)
            ((word: LLRLRL, 0), 1)
        """
        w,i = p
        w = self._words(w)
        i = int(i)
        e = int(e)

        if e==0: f=2
        elif e==1: f=3
        elif e==2: f=0
        elif e==3: f=1
        else:
            raise ValueError("e (={!r}) must be either 0,1,2 or 3".format(e))

        if i == 0:
            if e == 0:
                if w.is_empty():   lab=(w,2)
                elif w[-1] == 'L': lab=(w[:-1],1)
                elif w[-1] == 'R': lab=(w[:-1],3)
            if e == 1: lab=(w,0)
            if e == 2: lab=(w,2)
            if e == 3: lab=(w,0)
        elif i == 1:
            if e == 0: lab=(w + self._wL, 2)
            if e == 1: lab=(w,2)
            if e == 2: lab=(w + self._wL, 0)
            if e == 3: lab=(w,3)
        elif i == 2:
            if e == 0: lab=(w,0)
            if e == 1: lab=(w,3)
            if e == 2:
                if w.is_empty():   lab=(w,0)
                elif w[-1] == 'L': lab=(w[:-1],1)
                elif w[-1] == 'R': lab=(w[:-1],3)
            if e == 3: lab=(w,1)
        elif i == 3:
            if e == 0: lab=(w + self._wR, 2)
            if e == 1: lab=(w,1)
            if e == 2: lab=(w + self._wR, 0)
            if e == 3: lab=(w,2)
        else:
            raise ValueError("i (={!r}) must be either 0,1,2 or 3".format(i))

        # the fastest label constructor
        lab = self.polygon_labels()._cartesian_product_of_elements(lab)
        return lab,f

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
         w/r         w/r
        +---+------+---+
        | 1 |  2   | 3 |
        |   |      |   |  h2
        +---+------+---+
            |  0   | h1
            +------+
            w

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: T = sfg.tfractal_surface()
            sage: T.polygon(('L',0))
            Polygon: (0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)
            sage: T.polygon(('LRL',0))
            Polygon: (0, 0), (1/8, 0), (1/8, 1/8), (0, 1/8)
        """
        w = self._words(lab[0])
        return (1 / self._r ** w.length()) * self._base_polygon(lab[1])

    @cached_method
    def _base_polygon(self, i):
        from flatsurf.geometry.polygon import Polygons
        if i == 0:
            w = self._w
            h = self._h1
        if i == 1 or i == 3:
            w = self._w / self._r
            h = self._h2
        if i == 2:
            w = self._w
            h = self._h2
        return Polygons(self.base_ring())([(w,0),(0,h),(-w,0),(0,-h)])

    def base_label(self):
        return self.polygon_labels()._cartesian_product_of_elements((self._words(''), 0))

    def is_finite(self):
        return False

def tfractal_surface(w=ZZ_1, r=ZZ_2, h1=ZZ_1, h2=ZZ_1):
    return TranslationSurface(TFractalSurface(w,r,h1,h2))

class SimilaritySurfaceGenerators:
    r"""
    Examples of similarity surfaces.
    """
    @staticmethod
    def example():
        r"""
        Construct a SimilaritySurface from a pair of triangles.

        TESTS::

            sage: from flatsurf import *
            sage: ex = similarity_surfaces.example()
            sage: ex
            SimilaritySurface built from 2 polygons
            sage: TestSuite(ex).run()
        """
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.polygon import polygons
        from sage.rings.rational_field import QQ
        s=Surface_list(base_ring=QQ)
        s.add_polygon(polygons(vertices=[(0,0), (2,-2), (2,0)],ring=QQ)) # gets label 0
        s.add_polygon(polygons(vertices=[(0,0), (2,0), (1,3)],ring=QQ)) # gets label 1
        s.change_polygon_gluings(0, [(1,1), (1,2), (1,0)])
        s.make_immutable()
        return SimilaritySurface(s)

    @staticmethod
    def self_glued_polygon(P):
        r"""
        Return the HalfTranslationSurface formed by gluing all edges of P to themselves.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: p = polygons((2,0),(-1,3),(-1,-3))
            sage: s = similarity_surfaces.self_glued_polygon(p)
            sage: TestSuite(s).run()
        """
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.half_translation_surface import HalfTranslationSurface
        s = Surface_list(base_ring=P.base_ring(), mutable=True)
        s.add_polygon(P,[(0,i) for i in xrange(P.num_edges())])
        s.make_immutable()
        return HalfTranslationSurface(s)

    @staticmethod
    def billiard(P):
        r"""
        Return the ConeSurface associated to the billiard in the polygon ``P``.

        EXAMPLES::

            sage: from flatsurf import *
            sage: P = polygons(vertices=[(0,0), (1,0), (0,1)])
            sage: from flatsurf.geometry.rational_cone_surface import RationalConeSurface
            sage: Q = RationalConeSurface(similarity_surfaces.billiard(P))
            sage: Q
            RationalConeSurface built from 2 polygons
            sage: Q.underlying_surface()
            <class 'flatsurf.geometry.surface.Surface_list'>
            sage: M = Q.minimal_translation_cover()
            sage: M
            TranslationSurface built from 8 polygons
            sage: TestSuite(M).run()
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.cone_surface import ConeSurface
        from sage.matrix.constructor import matrix

        n = P.num_edges()
        r = matrix(2, [-1,0,0,1])
        Q = polygons(edges=[r*v for v in reversed(P.edges())])

        surface = Surface_list(base_ring = P.base_ring())
        surface.add_polygon(P) # gets label 0)
        surface.add_polygon(Q) # gets label 1
        surface.change_polygon_gluings(0,[(1,n-i-1) for i in xrange(n)])
        
        surface.make_immutable()
        s=ConeSurface(surface)
        gs=s.graphical_surface()
        gs.process_options(edge_labels=None,polygon_labels=False)
        gs.make_adjacent(0,0,reverse=True)
        return s

    @staticmethod
    def polygon_double(P):
        r"""
        Return the ConeSurface associated to the billiard in the polygon ``P``.
        Differs from billiard(P) only in the graphical display. Here, we display
        the polygons separately.
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.cone_surface import ConeSurface
        from sage.matrix.constructor import matrix

        n = P.num_edges()
        r = matrix(2, [-1,0,0,1])
        Q = polygons(edges=[r*v for v in reversed(P.edges())])

        surface = Surface_list(base_ring = P.base_ring())
        surface.add_polygon(P) # gets label 0)
        surface.add_polygon(Q) # gets label 1
        surface.change_polygon_gluings(0,[(1,n-i-1) for i in xrange(n)])
        surface.make_immutable()
        return ConeSurface(surface)

    @staticmethod
    def right_angle_triangle(w,h):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: R = similarity_surfaces.right_angle_triangle(2, 3)
            sage: R
            ConeSurface built from 2 polygons
            sage: TestSuite(R).run()
        """
        from sage.structure.sequence import Sequence
        from flatsurf.geometry.polygon import Polygons
        from sage.modules.free_module import VectorSpace
        from sage.modules.free_module_element import vector
        from flatsurf.geometry.cone_surface import ConeSurface

        F = Sequence([w,h]).universe()
        
        if not F.is_field():
            F = F.fraction_field()
        V = VectorSpace(F,2)
        P = Polygons(F)
        s=Surface_list(base_ring=F)
        s.add_polygon(P([V((w,0)),V((-w,h)),V((0,-h))])) # gets label 0
        s.add_polygon(P([V((0,h)),V((-w,-h)),V((w,0))])) # gets label 1
        s.change_polygon_gluings(0,[(1,2),(1,1),(1,0)])
        s.make_immutable()
        return ConeSurface(s)

    # Removed because Surface_polygons_and_gluings is gone.
    #
    #def __call__(self, *args, **kwds):
    #    from flatsurf.geometry.surface import Surface_polygons_and_gluings
    #    from flatsurf.geometry.similarity_surface import SimilaritySurface
    #    return SimilaritySurface(Surface_polygons_and_gluings(*args, **kwds))

similarity_surfaces = SimilaritySurfaceGenerators()

class TranslationSurfaceGenerators:
    r"""
    Common and less common translation surfaces.
    """
    @staticmethod
    def square_torus():
        r"""
        Return flat torus obtained by identification of the opposite sides of a
        square.

        EXAMPLES::

            sage: from flatsurf import *
            sage: T = translation_surfaces.square_torus()
            sage: T
            TranslationSurface built from 1 polygon
            sage: TestSuite(T).run()
            
        Rational directions are completely periodic::

            sage: v = T.tangent_vector(0, (1/33, 1/257), (13,17))
            sage: L = v.straight_line_trajectory()
            sage: L.flow(13+17)
            sage: L.is_closed()
            True

        TESTS::

            sage: TestSuite(T).run()
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.translation_surface import TranslationSurface
        from sage.rings.rational_field import QQ
        s=Surface_list(base_ring=QQ)
        s.add_polygon(polygons.square(),[(0,2),(0,3),(0,0),(0,1)])
        s.make_immutable()
        return TranslationSurface(s)

    @staticmethod
    def veech_2n_gon(n):
        r"""
        The regular 2n-gon with opposite sides identified.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.veech_2n_gon(5)
            sage: s.polygon(0)
            Polygon: (0, 0), (1, 0), (-1/2*a^2 + 5/2, 1/2*a), (-a^2 + 7/2, -1/2*a^3 + 2*a), (-1/2*a^2 + 5/2, -a^3 + 7/2*a), (1, -a^3 + 4*a), (0, -a^3 + 4*a), (1/2*a^2 - 3/2, -a^3 + 7/2*a), (a^2 - 5/2, -1/2*a^3 + 2*a), (1/2*a^2 - 3/2, 1/2*a)
            sage: TestSuite(s).run()
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.translation_surface import TranslationSurface
        p = polygons.regular_ngon(2*n)
        s = Surface_list(base_ring=p.base_ring())
        s.add_polygon(p,[ ( 0, (i+n)%(2*n) ) for i in xrange(2*n)] )
        s.make_immutable()
        return TranslationSurface(s)

    @staticmethod
    def veech_double_n_gon(n):
        r"""
        A pair of regular n-gons with each edge of one identified to an edge of the other to make a translation surface.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.veech_double_n_gon(5)
            sage: TestSuite(s).run()
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.translation_surface import TranslationSurface
        from sage.matrix.constructor import Matrix
        p = polygons.regular_ngon(n)
        s = Surface_list(base_ring=p.base_ring())
        m = Matrix([[-1,0],[0,-1]])
        s.add_polygon(p) # label=0
        s.add_polygon(m*p, [(0,i) for i in xrange(n)])
        s.make_immutable()
        return TranslationSurface(s)


    @staticmethod
    def regular_octagon():
        r"""
        Return the translation surface built from the regular octagon by
        identifying opposite sides.

        EXAMPLES::

            sage: from flatsurf import *
            sage: T = translation_surfaces.regular_octagon()
            sage: T
            TranslationSurface built from 1 polygon
            sage: T.stratum()
            H(2)
            sage: TestSuite(T).run()
        """
        return translation_surfaces.veech_2n_gon(4)

    @staticmethod
    def mcmullen_L(l1,l2,l3,l4):
        r"""
        Return McMullen's L shaped surface with parameters l1, l2, l3, l4.

        Polygon labels and lengths are marked below.
        
        +-----+
        |     |
        |  1  |l1
        |     |
        |     |    l4
        +-----+---------+
        |               |
        |       0       |l2
        |               |
        +-----+---------+
          l3        

        Note that this surface may not work correctly yet due to a non-strictly 
        convex polygon in the representation.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.mcmullen_L(1,1,1,1)
            sage: TestSuite(s).run()
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.translation_surface import TranslationSurface
        from sage.structure.sequence import Sequence
        
        field = Sequence([l1,l2,l3,l4]).universe()
        if not field.is_field():
            field = field.fraction_field()

        s=Surface_list(base_ring=field)
        s.add_polygon(polygons((l3,0),(l4,0),(0,l2),(-l4,0),(-l3,0),(0,-l2), ring=field))
        s.add_polygon(polygons((l3,0),(0,l1),(-l3,0),(0,-l1), ring=field))
        s.change_edge_gluing(0,0,1,2)
        s.change_edge_gluing(0,1,0,3)
        s.change_edge_gluing(0,2,0,5)
        s.change_edge_gluing(0,4,1,0)
        s.change_edge_gluing(1,1,1,3)
        s.make_immutable()
        return TranslationSurface(s)

    @staticmethod
    def ward(n):
        r"""
        Return the surface formed by gluing a regular 2n-gon to two regular n-gons.
        These surfaces have Veech's lattice property due to work of Ward.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.ward(3)
            sage: TestSuite(s).run()
            sage: s=translation_surfaces.ward(7)
            sage: TestSuite(s).run()
        """
        assert n>=3
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_list
        from flatsurf.geometry.translation_surface import TranslationSurface
        o = ZZ_2*polygons.regular_ngon(2*n)
        p1 = polygons(*[o.edge((2*i+n)%(2*n)) for i in xrange(n)])
        p2 = polygons(*[o.edge((2*i+n+1)%(2*n)) for i in xrange(n)])
        s = Surface_list(base_ring=o.parent().field())
        s.add_polygon(o)
        s.add_polygon(p1)
        s.add_polygon(p2)
        s.change_polygon_gluings(1, [(0,2*i) for i in xrange(n)])
        s.change_polygon_gluings(2, [(0,2*i+1) for i in xrange(n)])
        s.make_immutable()
        return TranslationSurface(s)

    @staticmethod
    def octagon_and_squares():
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: os = translation_surfaces.octagon_and_squares()
            sage: os
            TranslationSurface built from 3 polygons
            sage: TestSuite(os).run()
        """
        from flatsurf.geometry.polygon import polygons
        from sage.matrix.matrix_space import MatrixSpace
        from flatsurf.geometry.translation_surface import TranslationSurface
        return translation_surfaces.ward(4)

    @staticmethod
    def from_flipper(h):
        r"""
        Build a translation or half-translation surface from a flipper
        pseudo-Anosov

        EXAMPLES::

            sage: from flatsurf import *
            sage: import flipper                             # optional - flipper

        A torus example::

            sage: t1 = (0r,1r,2r)                            # optional - flipper
            sage: t2 = (~0r,~1r,~2r)                         # optional - flipper
            sage: T = flipper.create_triangulation([t1,t2])  # optional - flipper
            sage: L1 = T.lamination([1r,0r,1r])              # optional - flipper
            sage: L2 = T.lamination([0r,1r,1r])              # optional - flipper
            sage: h1 = L1.encode_twist()                     # optional - flipper
            sage: h2 = L2.encode_twist()                     # optional - flipper
            sage: h = h1*h2^(-1r)                            # optional - flipper
            sage: f = h.flat_structure()                     # optional - flipper
            sage: translation_surfaces.from_flipper(h)       # optional - flipper
            HalfTranslationSurface built from 2 polygons     # optional - flipper

        A non-orientable example::

            sage: T = flipper.load('SB_4')                   # optional - flipper
            sage: h = T.mapping_class('s_0S_1s_2S_3s_1S_2')  # optional - flipper
            sage: h.is_pseudo_anosov()                       # optional - flipper
            True
            sage: S = translation_surfaces.from_flipper(h)   # optional - flipper
            sage: S.num_polygons()                           # optional - flipper
            4
            sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_element_to_sage
            sage: a = flipper_nf_element_to_sage(h.dilatation())  # optional - flipper
        """
        from sage.modules.free_module import VectorSpace
        from flatsurf.geometry.polygon import ConvexPolygons
        from flatsurf.geometry.surface import surface_list_from_polygons_and_gluings
        from flatsurf.geometry.half_translation_surface import HalfTranslationSurface

        f = h.flat_structure()

        x = next(f.edge_vectors.itervalues()).x
        K = flipper_nf_to_sage(x.number_field)
        V = VectorSpace(K, 2)
        edge_vectors = {i: V((K(e.x.linear_combination), K(e.y.linear_combination)))
                for i,e in f.edge_vectors.iteritems()}

        to_polygon_number = {k:(i,2-j) for i,t in enumerate(f.triangulation) for j,k in enumerate(t)}

        C = ConvexPolygons(K)

        polys = []
        adjacencies = {}
        for i,t in enumerate(f.triangulation):
            for j,k in enumerate(t):
                adjacencies[(i,2-j)] = to_polygon_number[~k]
            try:
                poly = C([edge_vectors[i] for i in reversed(tuple(t))])
            except ValueError:
                raise ValueError("t = {}, edges = {}".format(
                    t, [edge_vectors[i].n(digits=6) for i in t]))
            polys.append(poly)

        return HalfTranslationSurface(surface_list_from_polygons_and_gluings(polys, adjacencies))

    @staticmethod
    def origami(r,u,rr=None,uu=None,domain=None):
        r"""
        Return the origami defined by the permutations ``r`` and ``u``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: o.underlying_surface()
            Origami defined by r=(1,2) and u=(1,3)
            sage: o.stratum()
            H(2)
            sage: TestSuite(o).run()
        """
        from flatsurf.geometry.translation_surface import Origami
        return TranslationSurface(Origami(r,u,rr,uu,domain))

    @staticmethod
    def infinite_staircase1():
        r"""
        Return the infinite staircase

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase1()
            sage: S
            TranslationSurface built from infinitely many polygons
            sage: TestSuite(S).run(skip='_test_pickling')
        """
        return infinite_staircase()

    @staticmethod
    def infinite_staircase2():
        r"""
        Return the infinite staircase built as an origami

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: S = translation_surfaces.infinite_staircase2()
            sage: S.underlying_surface()
            Origami defined by r=<function <lambda> at ...> and
            u=<function <lambda> at ...>
            sage: TestSuite(S).run(skip='_test_pickling')
        """
        from flatsurf.geometry.translation_surface import Origami, TranslationSurface
        return TranslationSurface(Origami(
                lambda x: x+1 if x%2 else x-1,  # r  (edge 1)
                lambda x: x-1 if x%2 else x+1,  # u  (edge 2)
                lambda x: x+1 if x%2 else x-1,  # rr (edge 3)
                lambda x: x-1 if x%2 else x+1,  # uu (edge 0)
                domain = ZZ))

    @staticmethod
    def t_fractal(w=ZZ_1, r=ZZ_2, h1=ZZ_1, h2=ZZ_1):
        r"""
        Return the T-fractal with parameters ``w``, ``r``, ``h1``, ``h2``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: tf = translation_surfaces.t_fractal().underlying_surface()
            sage: tf
            The T-fractal surface with parameters w=1, r=2, h1=1, h2=1
            sage: TestSuite(tf).run(skip='_test_pickling')
        """
        return tfractal_surface(w,r,h1,h2)

    @staticmethod
    def chamanara(alpha):
        r"""
        Return the Chamanara surface with parameter ``alpha``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C
            TranslationSurface built from infinitely many polygons
            sage: TestSuite(C).run(skip='_test_pickling')
        """
        from flatsurf.geometry.chamanara import chamanara_surface
        return chamanara_surface(alpha)

translation_surfaces = TranslationSurfaceGenerators()
