from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method

ZZ_1 = ZZ(1)
ZZ_2 = ZZ(2)

from flatsurf.geometry.surface import Surface
from flatsurf.geometry.translation_surface import TranslationSurface

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
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.polygon import polygons
        p0 = polygons(vertices=[(0,0), (2,-2), (2,0)])
        p1 = polygons(vertices=[(0,0), (2,0), (1,3)])
        ps = (p0,p1)
        glue={ (0,2):(1,0), (0,0):(1,1), (0,1):(1,2), (1,0):(0,2), (1,1):(0,0), (1,2):(0,1) }
        return SimilaritySurface(Surface_polygons_and_gluings(ps,glue))


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
            <class 'flatsurf.geometry.surface.Surface_polygons_and_gluings'>
            sage: M = Q.minimal_translation_cover()
            sage: M
            TranslationSurface built from 8 polygons
            sage: TestSuite(M).run()
        """
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.cone_surface import ConeSurface
        from sage.matrix.constructor import matrix

        n = P.num_edges()
        r = matrix(2, [-1,0,0,1])
        Q = polygons(*[r*v for v in reversed(P.edges())])

        glue = {(0,i): (1,n-i-1) for i in range(n)}
        glue.update({(1,i): (0,n-i-1) for i in range(P.num_edges())})
        return ConeSurface(Surface_polygons_and_gluings((P,P), glue))

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
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.cone_surface import ConeSurface

        F = Sequence([w,h]).universe()
        if not F.is_field():
            F = F.fraction_field()
        V = VectorSpace(F,2)
        P = Polygons(F)
        p1 = P([V((w,0)),V((-w,h)),V((0,-h))])
        p2 = P([V((0,h)),V((-w,-h)),V((w,0))])
        ps = (p1,p2)
        glue = {(0,0):(1,2),(0,1):(1,1),(0,2):(1,0)}
        return ConeSurface(Surface_polygons_and_gluings(ps,glue))

    def __call__(self, *args, **kwds):
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        return SimilaritySurface(Surface_polygons_and_gluings(*args, **kwds))

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
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.translation_surface import TranslationSurface
        polygons = [polygons.square()]
        identifications = {(0,0):(0,2), (0,1):(0,3)}
        return TranslationSurface(Surface_polygons_and_gluings(polygons, identifications))

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
        from flatsurf.geometry.polygon import polygons
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.translation_surface import TranslationSurface
        polygons = [polygons.regular_ngon(8)]
        identifications = {}
        identifications.update(dict(((0,i),(0,i+4)) for i in xrange(4)))
        return TranslationSurface(Surface_polygons_and_gluings(polygons, identifications))

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
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.translation_surface import TranslationSurface
        polygons = {
            0:polygons((l3,0),(l4,0),(0,l2),(-l4,0),(-l3,0),(0,-l2)), 
            1:polygons((l3,0),(0,l1),(-l3,0),(0,-l1)) }
        identifications = {
            (0,0):(1,2),
            (0,1):(0,3),
            (0,2):(0,5),
            (0,4):(1,0),
            (1,1):(1,3)}
        return TranslationSurface(Surface_polygons_and_gluings(polygons, identifications))

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
        from flatsurf.geometry.surface import Surface_polygons_and_gluings
        from flatsurf.geometry.translation_surface import TranslationSurface

        o = polygons.regular_ngon(8)
        K = o.parent().field()
        sqrt2 = K.gen()

        rot = MatrixSpace(K,2)([[sqrt2/ZZ_2,-sqrt2/ZZ_2],[sqrt2/ZZ_2,sqrt2/ZZ_2]])

        s = ZZ_2 * polygons.square(field=K)

        polygons = [o, s, rot * s]
        identifications = {
            (0,0): (1,3),
            (0,1): (2,3),
            (0,2): (1,0),
            (0,3): (2,0),
            (0,4): (1,1),
            (0,5): (2,1),
            (0,6): (1,2),
            (0,7): (2,2),
            }
        return TranslationSurface(Surface_polygons_and_gluings(polygons, identifications))

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
