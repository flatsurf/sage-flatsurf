r"""
This file implements polygons with
 - action of matrices in GL^+(2,R)
 - conversion between ground fields

EXAMPLES::

    sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
    sage: p = Polygons(K)([(1,0),(-sqrt2,1+sqrt2),(sqrt2-1,-1-sqrt2)])
    sage: p
    Polygon: (0, 0), (1, 0), (sqrt2, sqrt2)

    sage: M = MatrixSpace(K,2)
    sage: m = M([[1,1+sqrt2],[0,1]])
    sage: m * p
    Polygon: (0, 0), (1, 0), (2*sqrt2 + 2, sqrt2)

    sage: s = Polygons(QQ)([(0,0),(1,0),(2,2),(-2,3)])
"""

from sage.categories.sets_cat import Sets

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.rings.qqbar import AA

from matrix_2x2 import (angle,
                        is_similarity,
                        homothety_rotation_decomposition,
                        similarity_from_vectors,
                        rotation_matrix_angle)

# we implement action of GL(2,K) on polygons
from sage.categories.action import Action
import operator

def wedge_product(v,w):
    return v[0]*w[1]-v[1]*w[0]

class ActionOnPolygons(Action):
    def __init__(self, polygons):
        from sage.matrix.matrix_space import MatrixSpace
        K = polygons.field()
        Action.__init__(self, MatrixSpace(K,2), polygons, True, operator.mul)

    def _call_(self, g, x):
        if g.det() <= 0:
            # Maybe we can allow an action, which also reverses the edge ordering? -Pat
            raise ValueError("can not act with matrix with negative determinant")
        return x.parent()([g*e for e in x.edges()])

from sage.structure.element import Element
class Polygon(Element):
    r"""
    A polygon in RR^2 defined up to translation.

    .. NOTE::

    Should we precompute the angles ?
    """
    def __init__(self, parent, edges):
        r"""
        INPUT:

        - ``parent`` -- a parent

        - ``edges`` -- a list of vectors or couples whose sum is zero
        """
        Element.__init__(self, parent)
        field = parent.field()
        self._x = [field(e[0]) for e in edges]
        self._y = [field(e[1]) for e in edges]

        if sum(self._x) or sum(self._y):
            raise ValueError("the sum over the edges do not sum up to 0")

        # the following is very long because of the angle buisness!!!
        self._check()

    def _check(self):
        # check convexity
        edges = self.edges()
        for i in xrange(len(edges)-1):
            if wedge_product(edges[i],edges[i+1]) <= 0:
                raise ValueError("not convex!")
        if wedge_product(edges[-1],edges[0]) <= 0:
            raise ValueError("not convex!")

        # check angles (relatively long time)
        sum(self.angle(i) for i in xrange(self.num_edges())) == self.num_edges() - 3

    def base_ring(self):
        return self.parent().base_ring()

    field=base_ring

    def num_edges(self):
        return len(self._x)

    def _repr_(self):
        r"""
        String representation.
        """
        return "Polygon: " + ", ".join(map(str,self.vertices()))

    def vertices(self, translation=None):
        r"""
        Return the set of vertices as vectors.
        """
        V = self.parent().vector_space()
        if translation is None:
            zero = V.zero()
        else:
            try:
                translation + V.zero()
            except StandardError:
                from sage.modules.free_module_element import vector
                translation = vector(translation)
                try:
                    translation + V.zero()
                except StandardError:
                    raise ValueError("can not convert translation to a vector in R^2")
            zero = translation
        res = [zero]
        # This code returns a last vertex equal to the first vertex:
        #for x,y in zip(self._x,self._y):
        #    res.append(res[-1] + V((x,y)))
        for i in range(self.num_edges()-1):
            res.append(res[-1] + V((self._x[i],self._y[i])))
        return res

    def __iter__(self):
        return iter(self.vertices())

    def edges(self):
        r"""
        Return the set of edges as vectors.
        """
        V = self.parent().vector_space()
        return map(V, zip(self._x,self._y))

    def edge(self, i):
        r"""
        Return the ``i``-th edge of that polygon.
        """
        V = self.parent().vector_space()
        return V((self._x[i], self._y[i]))

    def plot(self, translation=None):
        r"""
        Plot the polygon with the origine at ``translation``.
        """
        from sage.plot.point import point2d
        from sage.plot.line import line2d
        from sage.plot.polygon import polygon2d
        from sage.modules.free_module import VectorSpace
        V = VectorSpace(RR,2)
        P = self.vertices(translation)
        return point2d(P, color='red') + line2d(P + [P[0]], color='orange') + polygon2d(P, alpha=0.3)

    def angle(self, e):
        r"""
        Return the angle at the begining of the start point of the edge ``e``.

        EXAMPLES::

            sage: square().angle(0)
            1/4
            sage: regular_octagon().angle(0)
            3/8
        """
        return angle(self.edge(e), - self.edge((e-1)%self.num_edges()))

    def area(self):
        r"""
        Return the area of self.

        EXAMPLES::

            sage: regular_octagon().area()
            8*sqrt2 + 8
            sage: square().area()
            1
            sage: (2*square()).area()
            4
        """
        x = self._x
        y = self._y
        zero = self.field().zero()
        x0 = self._x[0]
        y0 = self._y[0]
        x1 = x0 + self._x[1]
        y1 = y0 + self._y[1]
        a = 0
        for i in xrange(2,len(self._x)):
            a += (x0*y1 - x1*y0)/2
            x0 = x1; y0 = y1
            x1 += x[i]
            y1 += y[i]
        return a

from sage.structure.parent import Parent
class Polygons(Parent):
    Element = Polygon
    def __init__(self, field):
        Parent.__init__(self, category=Sets())

        if not AA.has_coerce_map_from(field):
            raise ValueError("the field must have a coercion to AA")
        self._field = field

        self.register_action(ActionOnPolygons(self))

    def has_coerce_map_from(self, other):
        return isinstance(other, Polygons) and self.field().has_coerce_map_from(other.field())

    def _an_element_(self):
        return self([(1,0),(0,1),(-1,0),(0,-1)])

    def base_ring(self):
        return self._field

    field = base_ring

    def _repr_(self):
        return "polygons with coordinates in %s"%self.base_ring()

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), 2)

    def _element_constructor_(self, data):
        return self.element_class(self, data)


def square(field=None):
    if field is None:
        field = QQ
    return Polygons(field)([(1,0),(0,1),(-1,0),(0,-1)])

def regular_octagon(field=None):
    r"""
    Return a regular octagon with sides of length 2.

    TODO: implement regular_ngons.
    """
    if field is None:
        from sage.rings.number_field.number_field import NumberField
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        R = PolynomialRing(ZZ,'x')
        x = R.gen()
        field = NumberField(x**2 - 2, 'sqrt2', embedding=AA(2).sqrt())
        sqrt2 = field.gen()

        # a faster way would be
        # field, sqrt2, _ = AA(2).sqrt().as_number_field_element()
        # but proceeding that way we get a field with no embedding in AA!!!!!
    else:
        sqrt = field.gen()

    edges = [(0,2), (-sqrt2,sqrt2), (-2,0), (-sqrt2,-sqrt2), (0,-2), (sqrt2,-sqrt2),
            (2,0), (sqrt2,sqrt2)]

    return Polygons(field)(edges)

class PolygonCreator():
    r"""
    Class for iteratively constructing a polygon over the field.
    """
    def __init__(self, field = QQ):
        self._v=[]
        self._w=[]
        self._field=field

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self._field, 2)

    def add_vertex(self, new_vertex):
        r"""
        Add a vertex to the polygon.
        Returns 1 if successful and 0 if not, in which case the resulting
        polygon would not have been convex.
        """
        V=self.vector_space()
        newv=V(new_vertex)
        if (len(self._v)==0):
            self._v.append(newv)
            self._w.append(V.zero())
            return 1
        if (len(self._v)==1):
            if (self._v[0]==newv):
                return 0
            else:
                self._w[-1]=newv-self._v[-1]
                self._w.append(self._v[0]-newv)
                self._v.append(newv)
                return 1
        if (len(self._v)>=2):
            neww1=newv-self._v[-1]
            if wedge_product(self._w[-2],neww1) <= 0:
                return 0
            neww2=self._v[0]-newv
            if wedge_product(neww1,neww2)<= 0:
                return 0
            if wedge_product(neww2,self._w[0])<= 0:
                return 0
            self._w[-1]=newv-self._v[-1]
            self._w.append(self._v[0]-newv)
            self._v.append(newv)
            return 1

    def get_polygon(self):
        r"""
        Return the polygon.
        Raises a ValueError if less than three vertices have been accepted.
        """
        if len(self._v)<2:
            raise ValueError("Not enough vertices!")
        return Polygons(self._field)(self._w)


