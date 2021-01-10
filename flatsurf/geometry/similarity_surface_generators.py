# -*- coding: utf-8 -*-
#*********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020      Julian Rüth
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
#*********************************************************************
from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems, itervalues

from sage.rings.all import ZZ, QQ, RIF, AA, NumberField, polygen
from sage.modules.all import VectorSpace, vector
from sage.structure.coerce import py_scalar_parent
from sage.structure.element import get_coercion_model, parent
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

ZZ_1 = ZZ(1)
ZZ_2 = ZZ(2)

from .polygon import polygons, ConvexPolygons, Polygon, ConvexPolygon, build_faces

from .surface import Surface, Surface_list
from .translation_surface import TranslationSurface
from .dilation_surface import DilationSurface
from .similarity_surface import SimilaritySurface
from .half_translation_surface import HalfTranslationSurface
from .cone_surface import ConeSurface
from .rational_cone_surface import RationalConeSurface


def flipper_nf_to_sage(K, name='a'):
    r"""
    Convert a flipper number field into a Sage number field

    .. NOTE::

        Currently, the code is not careful at all with root isolation.

    EXAMPLES::

        sage: import flipper  # optional - flipper
        sage: import realalg  # optional - flipper
        sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_to_sage
        sage: K = realalg.RealNumberField([-2r] + [0r]*5 + [1r])   # optional - flipper
        sage: K_sage = flipper_nf_to_sage(K)                       # optional - flipper
        sage: K_sage                                               # optional - flipper
        Number Field in a with defining polynomial x^6 - 2 with a = 1.122462048309373?
        sage: AA(K_sage.gen())                                     # optional - flipper
        1.122462048309373?
    """
    r = K.lmbda.interval()
    l = r.lower * ZZ(10)**(-r.precision)
    u = r.upper * ZZ(10)**(-r.precision)

    p = QQ['x'](K.coefficients)
    s = AA.polynomial_root(p, RIF(l,u))
    return NumberField(p, name, embedding=s)

def flipper_nf_element_to_sage(x, K=None):
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
    if K is None:
        K = flipper_nf_to_sage(x.field)
    coeffs = list(map(QQ, x.coefficients))
    coeffs.extend([0] * (K.degree() - len(coeffs)))
    return K(coeffs)

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
        if lambda_squared==None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R=PolynomialRing(ZZ,'x')
            x = R.gen()
            field=NumberField(x**3-ZZ(5)*x**2+ZZ(4)*x-ZZ(1), 'r', embedding=AA(ZZ(4)))
            self._l=field.gen()
        else:
            if field is None:
                self._l=lambda_squared
                field=lambda_squared.parent()
            else:
                self._l=field(lambda_squared)
        Surface.__init__(self,field, ZZ.zero(), finite=False)

    def _repr_(self):
        r"""
        String representation.
        """
        return "The E-infinity surface"

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
            return self.base_ring().one()
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
        return polygons.rectangle(2*self.get_black(lab),self.get_white(lab))

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
        from sage.combinat.words.words import Words

        field = Sequence([w,r,h1,h2]).universe()
        if not field.is_field():
            field = field.fraction_field()
        self._w = field(w)
        self._r = field(r)
        self._h1 = field(h1)
        self._h2 = field(h2)
        self._words = Words('LR', finite=True, infinite=False)
        self._wL = self._words('L')
        self._wR = self._words('R')

        base_label=self.polygon_labels()._cartesian_product_of_elements((self._words(''), 0))

        Surface.__init__(self, field, base_label, finite=False)

    def _repr_(self):
        return "The T-fractal surface with parameters w=%s, r=%s, h1=%s, h2=%s"%(
                self._w, self._r, self._h1, self._h2)

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
            sage: T.opposite_edge((w,0),0)
            ((word: LLRLR, 1), 2)
            sage: T.opposite_edge((w,0),1)
            ((word: LLRLRL, 0), 3)
            sage: T.opposite_edge((w,0),2)
            ((word: LLRLRL, 2), 0)
            sage: T.opposite_edge((w,0),3)
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
        if i == 0:
            w = self._w
            h = self._h1
        if i == 1 or i == 3:
            w = self._w / self._r
            h = self._h2
        if i == 2:
            w = self._w
            h = self._h2
        return ConvexPolygons(self.base_ring())([(w,0),(0,h),(-w,0),(0,-h)])

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
        s = Surface_list(base_ring=QQ)
        s.add_polygon(polygons(vertices=[(0,0), (2,-2), (2,0)],ring=QQ)) # gets label 0
        s.add_polygon(polygons(vertices=[(0,0), (2,0), (1,3)],ring=QQ)) # gets label 1
        s.change_polygon_gluings(0, [(1,1), (1,2), (1,0)])
        s.set_immutable()
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
        s = Surface_list(base_ring=P.base_ring(), mutable=True)
        s.add_polygon(P,[(0,i) for i in range(P.num_edges())])
        s.set_immutable()
        return HalfTranslationSurface(s)

    @staticmethod
    def billiard(P, rational=False):
        r"""
        Return the ConeSurface associated to the billiard in the polygon ``P``.

        INPUT:

        - ``P`` - a polygon

        - ``rational`` (boolean, default ``False``) - whether to assume that the polygon
          has all its angle rational multiple of pi.

        EXAMPLES::

            sage: from flatsurf import *

            sage: P = polygons(vertices=[(0,0), (1,0), (0,1)])
            sage: from flatsurf.geometry.rational_cone_surface import RationalConeSurface
            sage: Q = similarity_surfaces.billiard(P, rational=True)
            sage: Q
            RationalConeSurface built from 2 polygons
            sage: M = Q.minimal_cover(cover_type="translation")
            sage: M
            TranslationSurface built from 8 polygons
            sage: TestSuite(M).run()

        A non-convex examples (L-shape polygon)::

            sage: P = polygons(vertices=[(0,0), (2,0), (2,1), (1,1), (1,2), (0,2)], convex=False)
            sage: Q = similarity_surfaces.billiard(P, rational=True)
            sage: TestSuite(Q).run()
            sage: M = Q.minimal_cover(cover_type="translation")
            sage: TestSuite(M).run()
            sage: M.stratum()
            H_2(2, 0^5)

        A quadrilateral from Eskin-McMullen-Mukamel-Wright::

            sage: E = EquiangularPolygons(1, 1, 1, 7)
            sage: P = E.an_element()
            sage: S = similarity_surfaces.billiard(P)
            sage: TestSuite(S).run()
            sage: S = S.minimal_cover(cover_type="translation")
            sage: TestSuite(S).run()
            sage: S = S.erase_marked_points() # optional: pyflatsurf
            sage: TestSuite(S).run()
            sage: S, _ = S.normalized_coordinates()
            sage: TestSuite(S).run()

        Unfolding a triangle with non-algebraic lengths::

            sage: E = EquiangularPolygons(3, 3, 5)
            sage: from pyexactreal import ExactReals # optional: exactreal
            sage: R = ExactReals(E.base_ring()) # optional: exactreal
            sage: P = E(R.random_element()) # optional: exactreal
            sage: S = similarity_surfaces.billiard(P); S # optional: exactreal
            ConeSurface built from 2 polygons
            sage: TestSuite(S).run() # optional: exactreal

        """
        if not isinstance(P, Polygon):
            raise TypeError("invalid input")

        V = P.module()

        if not isinstance(P, ConvexPolygon):
            # triangulate non-convex ones
            base_ring = P.base_ring()
            C = ConvexPolygons(base_ring)
            comb_edges = P.triangulation()
            vertices = P.vertices()
            comb_triangles = build_faces(len(vertices), comb_edges)
            triangles = []
            internal_edges = []  # list (p1, e1, p2, e2)
            external_edges = []  # list (p1, e1)
            edge_to_lab = {}
            for num,(i,j,k) in enumerate(comb_triangles):
                triangles.append(C(vertices=[vertices[i], vertices[j], vertices[k]]))
                edge_to_lab[(i,j)] = (num, 0)
                edge_to_lab[(j,k)] = (num, 1)
                edge_to_lab[(k,i)] = (num, 2)
            for num,(i,j,k) in enumerate(comb_triangles):
                if (j, i) in edge_to_lab:
                    num2, e2 = edge_to_lab[j, i]
                    internal_edges.append((num, 0, num2, e2))
                else:
                    external_edges.append((num, 0))
                if (k, j) in edge_to_lab:
                    num2, e2 = edge_to_lab[k, j]
                    internal_edges.append((num, 1, num2, e2))
                else:
                    external_edges.append((num, 1))
                if (i, k) in edge_to_lab:
                    num2, e2 = edge_to_lab[i, k]
                    internal_edges.append((num, 2, num2, e2))
                else:
                    external_edges.append((num, 1))
            P = triangles
        else:
            internal_edges = []
            external_edges = [(0, i) for i in range(P.num_edges())]
            base_ring = P.base_ring()
            P = [P]

        m = len(P)
        Q = []
        surface = Surface_list(base_ring=base_ring)
        for p in P:
            surface.add_polygon(p)
        for p in P:
            surface.add_polygon(polygons(edges=[V((-x,y)) for x,y in reversed(p.edges())]))
        for (p1,e1,p2,e2) in internal_edges:
            surface.set_edge_pairing(p1, e1, p2, e2)
            ne1 = surface.polygon(p1).num_edges()
            ne2 = surface.polygon(p2).num_edges()
            surface.set_edge_pairing(m + p1, ne1-e1-1, m + p2, ne2-e2-1)
        for p, e in external_edges:
            ne = surface.polygon(p).num_edges()
            surface.set_edge_pairing(p, e, m+p, ne-e-1)

        surface.set_immutable()
        if rational:
            s = RationalConeSurface(surface)
        else:
            s = ConeSurface(surface)
        return s

    @staticmethod
    def polygon_double(P):
        r"""
        Return the ConeSurface associated to the billiard in the polygon ``P``.
        Differs from billiard(P) only in the graphical display. Here, we display
        the polygons separately.
        """
        from sage.matrix.constructor import matrix

        n = P.num_edges()
        r = matrix(2, [-1,0,0,1])
        Q = polygons(edges=[r*v for v in reversed(P.edges())])

        surface = Surface_list(base_ring = P.base_ring())
        surface.add_polygon(P) # gets label 0)
        surface.add_polygon(Q) # gets label 1
        surface.change_polygon_gluings(0,[(1,n-i-1) for i in range(n)])
        surface.set_immutable()
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
        from sage.modules.free_module_element import vector

        F = Sequence([w,h]).universe()

        if not F.is_field():
            F = F.fraction_field()
        V = VectorSpace(F,2)
        P = ConvexPolygons(F)
        s = Surface_list(base_ring=F)
        s.add_polygon(P([V((w,0)),V((-w,h)),V((0,-h))])) # gets label 0
        s.add_polygon(P([V((0,h)),V((-w,-h)),V((w,0))])) # gets label 1
        s.change_polygon_gluings(0,[(1,2),(1,1),(1,0)])
        s.set_immutable()
        return ConeSurface(s)

    # Removed because Surface_polygons_and_gluings is gone.
    #
    #def __call__(self, *args, **kwds):
    #    from flatsurf.geometry.surface import Surface_polygons_and_gluings
    #    from flatsurf.geometry.similarity_surface import SimilaritySurface
    #    return SimilaritySurface(Surface_polygons_and_gluings(*args, **kwds))

similarity_surfaces = SimilaritySurfaceGenerators()

class DilationSurfaceGenerators:
    @staticmethod
    def basic_dilation_torus(a):
        r"""
        Return a dilation torus built from a $1 \times 1$ square and a $a \times 1$ rectangle.
        Each edge of the square is glued to the opposite edge of the rectangle. This results
        in horizontal edges glued by a dilation with a scaling factor of a, and vertical
        edges being glued by translation.

            b       a
          +----+---------+
          | 0  | 1       |
        c |    |         | c
          +----+---------+
            a       b

        EXAMPLES::

            sage: from flatsurf import *
            sage: ds = dilation_surfaces.basic_dilation_torus(AA(sqrt(2)))
            sage: ds
            DilationSurface built from 2 polygons
            sage: TestSuite(ds).run()
        """
        s = Surface_list(base_ring=a.parent().fraction_field())
        CP = ConvexPolygons(s.base_ring())
        s.add_polygon(CP(edges=[(0,1),(-1,0),(0,-1),(1,0)])) # label 0
        s.add_polygon(CP(edges=[(0,1),(-a,0),(0,-1),(a,0)])) # label 1
        s.change_edge_gluing(0, 0, 1, 2)
        s.change_edge_gluing(0, 1, 1, 3)
        s.change_edge_gluing(0, 2, 1, 0)
        s.change_edge_gluing(0, 3, 1, 1)
        s.change_base_label(0)
        s.set_immutable()
        return DilationSurface(s)

    @staticmethod
    def genus_two_square(a, b, c, d):
        r"""A genus two dilation surface is returned.

        The unit square is made into an octagon by marking a point on
        each of its edges. Then opposite sides of this octagon are
        glued together by translation. (Since we currently require strictly
        convex polygons, we subdivide the square into a hexagon and two
        triangles as depicted below.) The parameters $a$, $b$, $c$,
        and $d$ should be real numbers strictly between zero and one.
        These represent the lengths of an edge of the resulting octagon,
        as below.
                 c
           +--+-------+
         d |2/        |
           |/         |
           +    0     +
           |         /|
           |        /1| b
           +-------+--+
              a
        The other edges will have length $1-a$, $1-b$, $1-c$, and $1-d$.
        Dilations used to glue edges will be by factors $c/a$, $d/b$,
        $(1-c)/(1-a)$ and $(1-d)/(1-b)$.

        EXAMPLES::

            sage: from flatsurf import *
            sage: ds = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
            sage: ds
            DilationSurface built from 3 polygons
            sage: TestSuite(ds).run()
        """
        field = Sequence([a, b, c, d]).universe().fraction_field()
        s = Surface_list(base_ring=QQ)
        CP = ConvexPolygons(field)
        hexagon = CP(edges=[(a,0), (1-a,b), (0,1-b), (-c,0), (c-1,-d), (0,d-1)])
        s.add_polygon(hexagon) # polygon 0
        s.change_base_label(0)
        triangle1 = CP(edges=[(1-a,0), (0,b), (a-1,-b)])
        s.add_polygon(triangle1) # polygon 1
        triangle2 = CP(edges=[(1-c,d), (c-1,0), (0,-d)])
        s.add_polygon(triangle2) # polygon 2
        s.change_edge_gluing(0, 0, 0, 3)
        s.change_edge_gluing(0, 2, 0, 5)
        s.change_edge_gluing(0, 1, 1, 2)
        s.change_edge_gluing(0, 4, 2, 0)
        s.change_edge_gluing(1, 0, 2, 1)
        s.change_edge_gluing(1, 1, 2, 2)
        s.set_immutable()
        return DilationSurface(s)

dilation_surfaces = DilationSurfaceGenerators()

class HalfTranslationSurfaceGenerators:
    # TODO: ideally, we should be able to construct a non-convex polygon and make the construction
    # below as a special case of billiard unfolding.
    @staticmethod
    def step_billiard(w, h):
        r"""
        Return a (finite) step billiard associated to the given widths ``w`` and heights ``h``.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1,1,1,1], [1,1/2,1/3,1/5])
            sage: S
            HalfTranslationSurface built from 8 polygons
            sage: TestSuite(S).run()
        """
        n = len(h)
        assert len(w) == n
        if n < 2:
            raise ValueError("w and h must have length at least 2")
        H = sum(h)
        W = sum(w)

        R = Sequence(w + h).universe()
        C = ConvexPolygons(R.fraction_field())

        P = []
        Prev = []
        x = 0
        y = H
        for i in range(n - 1):
            P.append(C(vertices=[(x, 0), (x + w[i], 0), (x + w[i], y - h[i]), (x + w[i], y), (x, y)]))
            x += w[i]
            y -= h[i]
        assert x == W - w[-1]
        assert y == h[-1]
        P.append(C(vertices=[(x, 0), (x + w[-1], 0), (x + w[-1], y), (x, y)]))

        Prev = [C(vertices=[(x, -y) for x,y in reversed(p.vertices())]) for p in P]

        S = Surface_list(base_ring = C.base_ring())
        S.rename("StepBilliard(w=[%s], h=[%s])" % (', '.join(map(str, w)), ', '.join(map(str, h))))
        S.add_polygons(P)    # get labels 0, ..., n-1
        S.add_polygons(Prev) # get labels n, n+1, ..., 2n-1

        # reflection gluings
        # (gluings between the polygon and its reflection)
        S.set_edge_pairing(0, 4, n, 4)
        S.set_edge_pairing(n-1, 0, 2*n-1, 2)
        S.set_edge_pairing(n-1, 1, 2*n-1, 1)
        S.set_edge_pairing(n-1, 2, 2*n-1, 0)
        for i in range(n-1):
            # set_edge_pairing(polygon1, edge1, polygon2, edge2)
            S.set_edge_pairing(i, 0, n+i, 3)
            S.set_edge_pairing(i, 2, n+i, 1)
            S.set_edge_pairing(i, 3, n+i, 0)

        # translation gluings
        S.set_edge_pairing(n-2, 1, n-1, 3)
        S.set_edge_pairing(2*n-2, 2, 2*n-1, 3)
        for i in range(n-2):
            S.set_edge_pairing(i, 1, i+1, 4)
            S.set_edge_pairing(n+i, 2, n+i+1, 4)

        S.set_immutable()
        return HalfTranslationSurface(S)

half_translation_surfaces = HalfTranslationSurfaceGenerators()

class TranslationSurfaceGenerators:
    r"""
    Common and less common translation surfaces.
    """
    @staticmethod
    def square_torus(a=1):
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
        return TranslationSurfaceGenerators.torus((a,0), (0,a))

    @staticmethod
    def torus(u, v):
        r"""
        Return the flat torus obtained as the quotient of the plane by the pair
        of vectors ``u`` and ``v``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, AA(2).sqrt()), (AA(3).sqrt(), 3))
            sage: T
            TranslationSurface built from 1 polygon
            sage: T.polygon(0)
            Polygon: (0, 0), (1, 1.414213562373095?), (2.732050807568878?, 4.414213562373095?), (1.732050807568878?, 3)
        """
        u = vector(u)
        v = vector(v)
        field = Sequence([u,v]).universe().base_ring()
        if isinstance(field, type):
            field = py_scalar_parent(field)
        if not field.is_field():
            field = field.fraction_field()
        s = Surface_list(base_ring=field)
        p = polygons(vertices=[(0,0), u, u+v, v], base_ring=field)
        s.add_polygon(p, [(0,2),(0,3),(0,0),(0,1)])
        s.set_immutable()
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
        p = polygons.regular_ngon(2*n)
        s = Surface_list(base_ring=p.base_ring())
        s.add_polygon(p,[ ( 0, (i+n)%(2*n) ) for i in range(2*n)] )
        s.set_immutable()
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
        from sage.matrix.constructor import Matrix
        p = polygons.regular_ngon(n)
        s = Surface_list(base_ring=p.base_ring())
        m = Matrix([[-1,0],[0,-1]])
        s.add_polygon(p) # label=0
        s.add_polygon(m*p, [(0,i) for i in range(n)])
        s.set_immutable()
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
            H_2(2)
            sage: TestSuite(T).run()
        """
        return translation_surfaces.veech_2n_gon(4)

    @staticmethod
    def mcmullen_genus2_prototype(w, h, t, e, rel=0, base_ring=None):
        r"""
        McMullen prototypes in the stratum H(2).

        These prototype appear at least in McMullen "Teichmüller curves in genus
        two: Discriminant and spin" (2004). The notation from that paper are
        quadruple ``(a, b, c, e)`` which translates in our notation as
        ``w = b``, ``h = c``, ``t = a`` (and ``e = e``).

        The associated discriminant is `D = e^2 + 4 wh`.

        If ``rel`` is a positive parameter (less than w-lambda) the surface belongs
        to the eigenform locus in H(1,1).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from surface_dynamics import AbelianStratum

            sage: prototypes = {
            ....:      5: [(1,1,0,-1)],
            ....:      8: [(1,1,0,-2), (2,1,0,0)],
            ....:      9: [(2,1,0,-1)],
            ....:     12: [(1,2,0,-2), (2,1,0,-2), (3,1,0,0)],
            ....:     13: [(1,1,0,-3), (3,1,0,-1), (3,1,0,1)],
            ....:     16: [(3,1,0,-2), (4,1,0,0)],
            ....:     17: [(1,2,0,-3), (2,1,0,-3), (2,2,0,-1), (2,2,1,-1), (4,1,0,-1), (4,1,0,1)],
            ....:     20: [(1,1,0,-4), (2,2,1,-2), (4,1,0,-2), (4,1,0,2)],
            ....:     21: [(1,3,0,-3), (3,1,0,-3)],
            ....:     24: [(1,2,0,-4), (2,1,0,-4), (3,2,0,0)],
            ....:     25: [(2,2,0,-3), (2,2,1,-3), (3,2,0,-1), (4,1,0,-3)]}

            sage: for D in sorted(prototypes):
            ....:     for w,h,t,e in prototypes[D]:
            ....:          T = translation_surfaces.mcmullen_genus2_prototype(w,h,t,e)
            ....:          assert T.stratum() == AbelianStratum(2)
            ....:          assert (D.is_square() and T.base_ring() is QQ) or (T.base_ring().polynomial().discriminant() == D)

        An example with some relative homology::

            sage: U8 = translation_surfaces.mcmullen_genus2_prototype(2,1,0,0,1/4)    # discriminant 8
            sage: U12 = translation_surfaces.mcmullen_genus2_prototype(3,1,0,0,3/10)   # discriminant 12

            sage: U8.stratum()
            H_2(1^2)
            sage: U8.base_ring().polynomial().discriminant()
            8
            sage: U8.j_invariant()
            (
                      [4 0]
            (0), (0), [0 2]
            )

            sage: U12.stratum()
            H_2(1^2)
            sage: U12.base_ring().polynomial().discriminant()
            12
            sage: U12.j_invariant()
            (
                      [6 0]
            (0), (0), [0 2]
            )
        """
        w = ZZ(w)
        h = ZZ(h)
        t = ZZ(t)
        e = ZZ(e)
        g = w.gcd(h)
        gg = g.gcd(t).gcd(e)
        if w <= 0 or h <= 0 or t < 0 or t >= g or not g.gcd(t).gcd(e).is_one() or e+h >= w:
            raise ValueError("invalid parameters")

        x = polygen(QQ)
        poly = x**2 - e * x - w*h
        if poly.is_irreducible():
            if base_ring is None:
                emb = AA.polynomial_root(poly, RIF(0,w))
                K = NumberField(poly, 'l', embedding=emb)
                l = K.gen()
            else:
                K = base_ring
                roots = poly.roots(K, multiplicities=False)
                if len(roots) != 2:
                    raise ValueError("invalid base ring")
                roots.sort(key=lambda x: x.numerical_approx())
                assert roots[0] < 0 and roots[0] > 0
                l = roots[1]
        else:
            if base_ring is None:
                K = QQ
            else:
                K = base_ring
            D = e**2 + 4 * w*h
            d = D.sqrt()
            l = (e + d) / 2

        try:
            rel = K(rel)
        except TypeError:
            K = get_coercion_model().common_parent(K, parent(rel))
            l = K(l)
            rel = K(rel)

        # (lambda,lambda) square on top
        # twisted (w,0), (t,h)
        s = Surface_list(base_ring=K)
        if rel:
            if rel < 0 or rel > w - l:
                raise ValueError("invalid rel argument")
            s.add_polygon(polygons(vertices=[(0,0),(l,0),(l+rel,l),(rel,l)], ring=K))
            s.add_polygon(polygons(vertices=[(0,0),(rel,0),(rel+l,0),(w,0),(w+t,h),(l+rel+t,h),(t+l,h),(t,h)], ring=K))
            s.set_edge_pairing(0, 1, 0, 3)
            s.set_edge_pairing(0, 0, 1, 6)
            s.set_edge_pairing(0, 2, 1, 1)
            s.set_edge_pairing(1, 2, 1, 4)
            s.set_edge_pairing(1, 3, 1, 7)
            s.set_edge_pairing(1, 0, 1, 5)
        else:
            s.add_polygon(polygons(vertices=[(0,0),(l,0),(l,l),(0,l)], ring=K))
            s.add_polygon(polygons(vertices=[(0,0),(l,0),(w,0),(w+t,h),(l+t,h),(t,h)], ring=K))
            s.set_edge_pairing(0, 1, 0, 3)
            s.set_edge_pairing(0, 0, 1, 4)
            s.set_edge_pairing(0, 2, 1, 0)
            s.set_edge_pairing(1, 1, 1, 3)
            s.set_edge_pairing(1, 2, 1, 5)
        s.set_immutable()
        return TranslationSurface(s)

    @staticmethod
    def mcmullen_L(l1, l2, l3, l4):
        r"""
        Return McMullen's L shaped surface with parameters l1, l2, l3, l4.

        Polygon labels and lengths are marked below::

            +-----+
            |     |
            |  1  |l1
            |     |
            |     |    l4
            +-----+---------+
            |     |         |
            |  0  |    2    |l2
            |     |         |
            +-----+---------+
              l3

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.mcmullen_L(1,1,1,1)
            sage: TestSuite(s).run()

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: translation_surfaces.mcmullen_L(1r, 1r, 1r, 1r)
            TranslationSurface built from 3 polygons
        """
        field = Sequence([l1,l2,l3,l4]).universe()
        if isinstance(field, type):
            field = py_scalar_parent(field)
        if not field.is_field():
            field = field.fraction_field()

        s = Surface_list(base_ring=field)
        s.add_polygon(polygons((l3,0),(0,l2),(-l3,0),(0,-l2), ring=field))
        s.add_polygon(polygons((l3,0),(0,l1),(-l3,0),(0,-l1), ring=field))
        s.add_polygon(polygons((l4,0),(0,l2),(-l4,0),(0,-l2), ring=field))
        s.change_edge_gluing(0,0,1,2)
        s.change_edge_gluing(0,1,2,3)
        s.change_edge_gluing(0,2,1,0)
        s.change_edge_gluing(0,3,2,1)
        s.change_edge_gluing(1,1,1,3)
        s.change_edge_gluing(2,0,2,2)
        s.set_immutable()
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
        o = ZZ_2*polygons.regular_ngon(2*n)
        p1 = polygons(*[o.edge((2*i+n)%(2*n)) for i in range(n)])
        p2 = polygons(*[o.edge((2*i+n+1)%(2*n)) for i in range(n)])
        s = Surface_list(base_ring=o.parent().field())
        s.add_polygon(o)
        s.add_polygon(p1)
        s.add_polygon(p2)
        s.change_polygon_gluings(1, [(0,2*i) for i in range(n)])
        s.change_polygon_gluings(2, [(0,2*i+1) for i in range(n)])
        s.set_immutable()
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
        return translation_surfaces.ward(4)

    @staticmethod
    def cathedral(a, b):
        r"""
        Return the cathedral surface with parameters ``a`` and ``b``.

        For any parameter ``a`` and ``b``, the cathedral surface belongs to the
        so-called Gothic locus described in McMullen, Mukamel, Wright "Cubic
        curves and totally geodesic subvarieties of moduli space" (2017).

                     1
                   <--->

                    /\           2a
                   /  \      +------+
               a  b|   | a  /        \
             +----+    +---+          +
             |    |    |   |          |
            1| P0 |P1  |P2 |  P3      |
             |    |    |   |          |
             +----+    +---+          +
                 b|    |    \        /
                   \  /      +------+
                    \/

        If a and b satisfies

        .. MATH::

            a = x + y \sqrt(d) \qquad b = -3x -3/2 + 3y \sqrt(d)

        for some rational x,y and d >= 0 then it is a Teichmüller curve.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.cathedral(1,2)
            sage: C.stratum()
            H_4(2^3)
            sage: TestSuite(C).run()

            sage: from pyexactreal import ExactReals # optional: exactreal
            sage: K = QuadraticField(5, embedding=AA(5).sqrt())
            sage: R = ExactReals(K) # optional: exactreal
            sage: C = translation_surfaces.cathedral(K.gen(), R.random_element([0.1, 0.2])) # optional: exactreal
            sage: C.stratum() # optional: exactreal
            H_4(2^3)
            sage: TestSuite(C).run() # optional: exactreal
        """
        ring = Sequence([a,b]).universe()
        if isinstance(ring, type):
            ring = py_scalar_parent(ring)
        if not ring.has_coerce_map_from(QQ):
            ring = ring.fraction_field()
        a = ring(a)
        b = ring(b)
        P = ConvexPolygons(ring)
        s = Surface_list(base_ring=ring)
        half = QQ((1,2))
        p0 = P(vertices=[(0,0),(a,0),(a,1),(0,1)])
        p1 = P(vertices=[(a,0),(a,-b),(a+half,-b-half),(a+1,-b),(a+1,0),(a+1,1),(a+1,b+1),(a+half,b+1+half),(a,b+1),(a,1)])
        p2 = P(vertices=[(a+1,0),(2*a+1,0),(2*a+1,1),(a+1,1)])
        p3 = P(vertices=[(2*a+1,0), (2*a+1+half,-half),(4*a+1+half,-half),(4*a+2,0),(4*a+2,1),(4*a+1+half,1+half),(2*a+1+half,1+half),(2*a+1,1)])
        s.add_polygon(p0)
        s.add_polygon(p1)
        s.add_polygon(p2)
        s.add_polygon(p3)
        s.set_edge_pairing(0, 0, 0, 2)
        s.set_edge_pairing(0, 1, 1, 9)
        s.set_edge_pairing(0, 3, 3, 3)
        s.set_edge_pairing(1, 0, 1, 3)
        s.set_edge_pairing(1, 1, 3, 4)
        s.set_edge_pairing(1, 2, 3, 6)
        s.set_edge_pairing(1, 4, 2, 3)
        s.set_edge_pairing(1, 5, 1, 8)
        s.set_edge_pairing(1, 6, 3, 0)
        s.set_edge_pairing(1, 7, 3, 2)
        s.set_edge_pairing(2, 0, 2, 2)
        s.set_edge_pairing(2, 1, 3, 7)
        s.set_edge_pairing(3, 1, 3, 5)
        s.set_immutable()
        return TranslationSurface(s)

    @staticmethod
    def arnoux_yoccoz(genus):
        r"""
        Construct the Arnoux-Yoccoz surface of genus 3 or greater.

        This presentation of the surface follows Section 2.3 of
        Joshua P. Bowman's paper "The Complete Family of Arnoux-Yoccoz
        Surfaces."

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.arnoux_yoccoz(4)
            sage: TestSuite(s).run()
            sage: s.is_delaunay_decomposed()
            True
            sage: s = s.canonicalize()
            sage: field=s.base_ring()
            sage: a = field.gen()
            sage: from sage.matrix.constructor import Matrix
            sage: m = Matrix([[a,0],[0,~a]])
            sage: ss = m*s
            sage: ss = ss.canonicalize()
            sage: s.cmp(ss) == 0
            True

        The Arnoux-Yoccoz pseudo-Anosov are known to have (minimal) invariant
        foliations with SAF=0::

            sage: S3 = translation_surfaces.arnoux_yoccoz(3)
            sage: Jxx, Jyy, Jxy = S3.j_invariant()
            sage: Jxx.is_zero() and Jyy.is_zero()
            True
            sage: Jxy
            [ 0  2  0]
            [ 2 -2  0]
            [ 0  0  2]

            sage: S4 = translation_surfaces.arnoux_yoccoz(4)
            sage: Jxx, Jyy, Jxy = S4.j_invariant()
            sage: Jxx.is_zero() and Jyy.is_zero()
            True
            sage: Jxy
            [ 0  2  0  0]
            [ 2 -2  0  0]
            [ 0  0  2  2]
            [ 0  0  2  0]
        """
        g=ZZ(genus)
        assert g>=3
        x = polygen(AA)
        p=sum([x**i for i in range(1,g+1)])-1
        cp = AA.common_polynomial(p)
        alpha_AA = AA.polynomial_root(cp, RIF(1/2, 1))
        field=NumberField(alpha_AA.minpoly(),'alpha',embedding=alpha_AA)
        a=field.gen()
        V = VectorSpace(field,2)
        p=[None for i in range(g+1)]
        q=[None for i in range(g+1)]
        p[0]=V(( (1-a**g)/2, a**2/(1-a) ))
        q[0]=V(( -a**g/2, a ))
        p[1]=V(( -(a**(g-1)+a**g)/2, (a-a**2+a**3)/(1-a) ))
        p[g]=V(( 1+(a-a**g)/2, (3*a-1-a**2)/(1-a) ))
        for i in range(2,g):
            p[i]=V(( (a-a**i)/(1-a) , a/(1-a) ))
        for i in range(1,g+1):
            q[i]=V(( (2*a-a**i-a**(i+1))/(2*(1-a)), (a-a**(g-i+2))/(1-a) ))
        P = ConvexPolygons(field)
        s = Surface_list(field)
        T = [None] * (2*g+1)
        Tp = [None] * (2*g+1)
        from sage.matrix.constructor import Matrix
        m=Matrix([[1,0],[0,-1]])
        for i in range(1,g+1):
            # T_i is (P_0,Q_i,Q_{i-1})
            T[i]=s.add_polygon(P(edges=[ q[i]-p[0], q[i-1]-q[i], p[0]-q[i-1] ]))
            # T_{g+i} is (P_i,Q_{i-1},Q_{i})
            T[g+i]=s.add_polygon(P(edges=[ q[i-1]-p[i], q[i]-q[i-1], p[i]-q[i] ]))
            # T'_i is (P'_0,Q'_{i-1},Q'_i)
            Tp[i]=s.add_polygon(m*s.polygon(T[i]))
            # T'_{g+i} is (P'_i,Q'_i, Q'_{i-1})
            Tp[g+i]=s.add_polygon(m*s.polygon(T[g+i]))
        for i in range(1,g):
            s.change_edge_gluing(T[i],0,T[i+1],2)
            s.change_edge_gluing(Tp[i],2,Tp[i+1],0)
        for i in range(1,g+1):
            s.change_edge_gluing(T[i],1,T[g+i],1)
            s.change_edge_gluing(Tp[i],1,Tp[g+i],1)
        #P 0 Q 0 is paired with P' 0 Q' 0, ...
        s.change_edge_gluing(T[1],2,Tp[g],2)
        s.change_edge_gluing(Tp[1],0,T[g],0)
        # P1Q1 is paired with P'_g Q_{g-1}
        s.change_edge_gluing(T[g+1],2,Tp[2*g],2)
        s.change_edge_gluing(Tp[g+1],0,T[2*g],0)
        # P1Q0 is paired with P_{g-1} Q_{g-1}
        s.change_edge_gluing(T[g+1],0,T[2*g-1],2)
        s.change_edge_gluing(Tp[g+1],2,Tp[2*g-1],0)
        # PgQg is paired with Q1P2
        s.change_edge_gluing(T[2*g],2,T[g+2],0)
        s.change_edge_gluing(Tp[2*g],0,Tp[g+2],2)
        for i in range(2,g-1):
            # PiQi is paired with Q'_i P'_{i+1}
            s.change_edge_gluing(T[g+i],2,Tp[g+i+1],2)
            s.change_edge_gluing(Tp[g+i],0,T[g+i+1],0)
        s.set_immutable()
        return TranslationSurface(s)

    @staticmethod
    def from_flipper(h):
        r"""
        Build a (half-)translation surface from a flipper pseudo-Anosov.

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
            sage: ts = translation_surfaces.from_flipper(h)  # optional - flipper
            sage: ts                                         # optional - flipper
            HalfTranslationSurface built from 2 polygons
            sage: TestSuite(ts).run()                        # optional - flipper

        A non-orientable example::

            sage: T = flipper.load('SB_4')                   # optional - flipper
            sage: h = T.mapping_class('s_0S_1s_2S_3s_1S_2')  # optional - flipper
            sage: h.is_pseudo_anosov()                       # optional - flipper
            True
            sage: S = translation_surfaces.from_flipper(h)   # optional - flipper
            sage: TestSuite(S).run()                         # optional - flipper
            sage: S.num_polygons()                           # optional - flipper
            4
            sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_element_to_sage
            sage: a = flipper_nf_element_to_sage(h.dilatation())  # optional - flipper
        """
        from .surface import surface_list_from_polygons_and_gluings

        f = h.flat_structure()

        x = next(itervalues(f.edge_vectors)).x
        K = flipper_nf_to_sage(x.field)
        V = VectorSpace(K, 2)
        edge_vectors = {i: V((flipper_nf_element_to_sage(e.x, K),
                              flipper_nf_element_to_sage(e.y, K)))
                for i,e in iteritems(f.edge_vectors)}

        to_polygon_number = {k:(i,j) for i,t in enumerate(f.triangulation) for j,k in enumerate(t)}

        C = ConvexPolygons(K)

        polys = []
        adjacencies = {}
        for i,t in enumerate(f.triangulation):
            for j,k in enumerate(t):
                adjacencies[(i,j)] = to_polygon_number[~k]
            try:
                poly = C([edge_vectors[i] for i in tuple(t)])
            except ValueError:
                raise ValueError("t = {}, edges = {}".format(
                    t, [edge_vectors[i].n(digits=6) for i in t]))
            polys.append(poly)

        return HalfTranslationSurface(surface_list_from_polygons_and_gluings(polys, adjacencies))

    @staticmethod
    def origami(r, u, rr=None, uu=None, domain=None):
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
            H_2(2)
            sage: TestSuite(o).run()
        """
        from .translation_surface import Origami
        return TranslationSurface(Origami(r,u,rr,uu,domain))

    @staticmethod
    def infinite_staircase():
        r"""
        Return the infinite staircase built as an origami.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: S = translation_surfaces.infinite_staircase()
            sage: S.underlying_surface()
            The infinite staircase
            sage: TestSuite(S).run(skip='_test_pickling')
        """
        from .translation_surface import Origami
        o = Origami(
                lambda x: x+1 if x%2 else x-1,  # r  (edge 1)
                lambda x: x-1 if x%2 else x+1,  # u  (edge 2)
                lambda x: x+1 if x%2 else x-1,  # rr (edge 3)
                lambda x: x-1 if x%2 else x+1,  # uu (edge 0)
                domain = ZZ,
                base_label=ZZ(0))
        o.rename("The infinite staircase")
        s = TranslationSurface(o)
        from flatsurf.geometry.similarity import SimilarityGroup
        SG = SimilarityGroup(QQ)
        def pos(n):
            if n%2 == 0:
                return SG((n/2,n/2))
            else:
                return SG((n//2,n//2+1))
        gs=s.graphical_surface(default_position_function = pos)
        gs.make_all_visible(limit = 10)
        return s

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

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.e_infinity_surface()
            sage: TestSuite(s).run(skip='_test_pickling')
        """
        return TranslationSurface(EInfinitySurface(lambda_squared, field))


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
        from .chamanara import chamanara_surface
        return chamanara_surface(alpha)

translation_surfaces = TranslationSurfaceGenerators()
