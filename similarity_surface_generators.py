from sage.rings.integer_ring import ZZ

ZZ_1 = ZZ(1)
ZZ_2 = ZZ(2)
ZZ_3 = ZZ(3)

from similarity_surface import TranslationSurface_generic
class InfiniteStaircase(TranslationSurface_generic):
    def _repr_(self):
        return "The infinite staircase"
    def base_ring(self):
        from sage.rings.rational_field import QQ
        return QQ

    def polygon(self, lab):
        if lab not in self.polygon_labels():
            raise ValueError("lab (=%s) not a valid label"%lab)
        from polygon import square
        return square()

    def polygon_labels(self):
        return ZZ

    def opposite_edge(self, p, e):
        if ((p%2) + (e%2)) % 2:
            return p+1,(e+2)%4
        else:
            return p-1,(e+2)%4

class TFractal(TranslationSurface_generic):
    r"""
     w/r         w/r
    +---+------+---+
    | 1 |  2   | 3 |
    |   |      |   |  h2
    +---+------+---+
        |  0   | h1
        +------+
            w

    where r is some ratio
    """
    def __init__(self, w=ZZ_1, r=ZZ_3, h1=ZZ_1, h2=ZZ_1):
        from sage.structure.sequence import Sequence
        from sage.combinat.words.words import Words

        self._field = Sequence([w,r,h1,h2]).universe()
        if not self._field.is_field():
            self._field = self._field.fraction_field()
        self._w = self._field(w)
        self._r = self._field(r)
        self._h1 = self._field(h1)
        self._h2 = self._field(h2)
        self._words = Words('LR')

    def _repr_(self):
        return "The T-fractal surface with parameters w=%s, r=%s, h1=%s, h2=%s"%(
                self._w, self._r, self._h1, self._h2)

    def base_ring(self):
        return self._field

    def polygon_labels(self):
        from sage.combinat.words.words import Words
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        from cartesian_product import CartesianProduct # custom cartesian product

        return CartesianProduct([self._words, FiniteEnumeratedSet([0,1,2,3])])

    def opposite_edge(self, p, e):
        r"""
         w/r         w/r
        +---+------+---+
        | 1 |  2   | 3 |
        |   |      |   |  h2
        +---+------+---+
            |  0   | h1
            +------+
            w
        """
        w,i = p
        if i == 0:
            if e == 0:
                if w.is_empty():   return (w,2),2
                elif w[-1] == 'L': return (w[:-1],1),2
                elif w[-1] == 'R': return (w[:-1],3),2
            if e == 1: return (w,0),3
            if e == 2: return (w,2),0
            if e == 3: return (w,0),1
        if i == 1:
            if e == 0: return (w + self._words('L'), 2), 2
            if e == 1: return (w,2),3
            if e == 2: return (w + self._words('L'), 0), 0
            if e == 3: return (w,3), 1
        if i == 2:
            if e == 0: return (w,0),2
            if e == 1: return (w,3),3
            if e == 2:
                if w.is_empty():   return (w,0),0
                elif w[-1] == 'L': return (w[:-1],1),0
                elif w[-1] == 'R': return (w[:-1],3),0
            if e == 3: return (w,1),1
        if i == 3:
            if e == 0: return (w + self._words('R'), 2), 2
            if e == 1: return (w,1),3
            if e == 2: return (w + self._words('R'), 0), 0
            if e == 3: return (w,2),1

    def polygon(self, lab):
        r"""
         w/r         w/r
        +---+------+---+
        | 1 |  2   | 3 |
        |   |      |   |  h2
        +---+------+---+
            |  0   | h1
            +------+
            w
        """
        from polygon import Polygons
        w,i = lab
        n = w.length()
        if i == 0:
            w = self._w / self._r**n
            h = self._h1 / self._r**n
        if i == 1 or i == 3:
            w = self._w / self._r**(n+1)
            h = self._h2 / self._r**n
        if i == 2:
            w = self._w / self._r**n
            h = self._h2 / self._r**n
        return Polygons(self.base_ring())([(w,0),(0,h),(-w,0),(0,-h)])

    def base_label(self):
        return (self._words(''), 0)


class SimilaritySurfaceGenerators:
    r"""
    Examples of similarity surfaces.
    """
    @staticmethod
    def example():
        r"""
        Construct a SimilaritySurface from a pair of triangles.
        """
        from similarity_surface import SimilaritySurface_polygons_and_gluings
        from polygon import PolygonCreator
        pc=PolygonCreator()
        pc.add_vertex((0,0))
        pc.add_vertex((2,-2))
        pc.add_vertex((2,0))
        p0=pc.get_polygon()
        pc=PolygonCreator()
        pc.add_vertex((0,0))
        pc.add_vertex((2,0))
        pc.add_vertex((1,3))
        p1=pc.get_polygon()
        ps=(p0,p1)
        glue={ (0,2):(1,0), (0,0):(1,1), (0,1):(1,2), (1,0):(0,2), (1,1):(0,0), (1,2):(0,1) }
        return SimilaritySurface_polygons_and_gluings(ps,glue)


class TranslationSurfaceGenerators:
    r"""
    Common and less common translation surfaces.
    """
    @staticmethod
    def regular_octagon():
        r"""
        Return the translation surface built from the regular octagon by
        identifying opposite sides.

        EXAMPLES::

            sage: T = translation_surfaces.regular_octagon()
            sage: T
            Translation surface built from the regular octagon
            sage: T.stratum()
            H_2(2)
        """
        from polygon import regular_octagon
        from similarity_surface import TranslationSurface_polygons_and_gluings
        polygons = [regular_octagon()]
        identifications = {}
        identifications.update(dict(((0,i),(0,i+4)) for i in xrange(4)))
        return TranslationSurface_polygons_and_gluings(polygons=polygons, identifications=identifications)

    @staticmethod
    def octagon_and_squares():
        from polygon import square, regular_octagon
        from sage.matrix.matrix_space import MatrixSpace
        from similarity_surface import TranslationSurface_polygons_and_gluings

        o = regular_octagon()
        K = o.parent().field()
        sqrt2 = K.gen()
        rot = MatrixSpace(K,2)([[sqrt2/ZZ_2,-sqrt2/ZZ_2],[sqrt2/ZZ_2,sqrt2/ZZ_2]])
        polygons = [regular_octagon(), ZZ_2*square(K), ZZ_2*rot*square(K)]
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
        return TranslationSurface_polygons_and_gluings(polygons=polygons, identifications=identifications)

    @staticmethod
    def origami(r,u,rr=None,uu=None,domain=None):
        r"""
        Return the origami defined by the permutations ``r`` and ``u``.

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: o
            Origami defined by r=(1,2) and u=(1,3)
            sage: o.stratum()
            H_2(2)
        """
        from similarity_surface import Origami
        return Origami(r,u,rr,uu,domain)

    @staticmethod
    def infinite_origami_example():
        from similarity_surface import Origami
        return Origami(
            lambda x: x+1,
            lambda x: x-1,
            lambda x: x-1,
            lambda x: x+1,
            ZZ)

    @staticmethod
    def infinite_staircase():
        return InfiniteStaircase()

    @staticmethod
    def t_fractal(w=ZZ_1, r=ZZ_3, h1=ZZ_1, h2=ZZ_1):
        return TFractal(w,r,h1,h2)
