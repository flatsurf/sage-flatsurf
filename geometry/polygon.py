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

from geometry.matrix_2x2 import (angle,
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

class PolygonPosition:
    r"""
    Class for describing the position of a point within or outside of a polygon.
    """
    # Position Types:
    OUTSIDE = 0
    INTERIOR = 1
    EDGE_INTERIOR = 2
    VERTEX = 3

    def __init__(self, position_type, edge = None, vertex = None):
        self._position_type=position_type
        if self.is_vertex():
            if vertex is None:
                raise ValueError("Constructed vertex position with no specified vertex.")
            self._vertex=vertex
        if self.is_in_edge_interior():
            if edge is None:
                raise ValueError("Constructed edge position with no specified edge.")
            self._edge=edge
    
    def __repr__(self):
        if self.is_outside():
            return "point positioned outside polygon"
        if self.is_in_interior():
            return "point positioned in interior of polygon"
        if self.is_in_edge_interior():
            return "point positioned on interior of edge "+str(self._edge)+" of polygon"
        return "point positioned on vertex "+str(self._vertex)+" of polygon"

    def is_outside(self):
        return self._position_type == PolygonPosition.OUTSIDE

    def is_inside(self):
        r""" 
        Return true if the position is not outside the closure of the polygon
        """
        return bool(self._position_type)

    def is_in_interior(self):
        return self._position_type == PolygonPosition.INTERIOR

    def is_in_edge_interior(self):
        return self._position_type == PolygonPosition.EDGE_INTERIOR

    def is_vertex(self):
        return self._position_type == PolygonPosition.VERTEX

    def get_position_type(self):
        return self._position_type

    def get_edge(self):
        return self._edge
    
    def get_vertex(self):
        return self._vertex

from sage.structure.element import Element
class Polygon(Element):
    r"""
    A polygon in RR^2 defined up to translation.

    .. NOTE::

    Should we precompute the angles ?
    No, in my oppinion. -Pat
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

        # Linear Time Sanity Checks
        if sum(self._x) or sum(self._y):
            raise ValueError("the sum over the edges do not sum up to 0")
        self._convexity_check()

    def _convexity_check(self):
        # Updated convexity check
        edges = self.edges()
        pc=PolygonCreator(field=self.base_ring())
        for v in self.vertices():
            if not pc.add_vertex(v):
                raise ValueError("not convex!")

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

    def vector_space(self):
        return self.parent().vector_space()

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

    def contains_point(self,point,translation=None):
        r"""
        Return true if the point is within the polygon (after the polygon is possibly translated)
        """
        return self.get_point_position(point,translation=translation).is_inside()

    def get_point_position(self,point,translation=None):
        r"""
        Get a combinatorial position of a points position compared to the polygon

        INPUT:

        - ``point`` -- a point in the plane (vector over the underlying base_ring())

        - ``translation`` -- optional translation to applied to the polygon (vector over the underlying base_ring())

        OUTPUT:

        - a PolygonPosition object

        EXAMPLES::

            sage: s=square()
            sage: V=s.parent().vector_space()
            sage: s.get_point_position(V((1/2,1/2)))
            point positioned in interior of polygon
            sage: s.get_point_position(V((1,0)))
            point positioned on vertex 1 of polygon
            sage: s.get_point_position(V((1,1/2)))
            point positioned on interior of edge 1 of polygon
            sage: s.get_point_position(V((1,3/2)))
            point positioned outside polygon
        """
        V = self.parent().vector_space()
        if translation is None:
            v1=V.zero()
        else:
            v1=translation
        for i in range(self.num_edges()):
            v0=v1
            e=V((self._x[i],self._y[i]))
            v1=v0+e
            w=wedge_product(e,point-v0)
            if w < 0:
                return PolygonPosition(PolygonPosition.OUTSIDE)
            if w == 0:
                # Lies on the line through edge i!
                n=self.num_edges()
                # index and edge after v1 
                ip1=(i+1)%n
                e=V((self._x[ip1],self._y[ip1]))
                w=wedge_product(e,point-v1)
                if w<0:
                    return PolygonPosition(PolygonPosition.OUTSIDE)
                if w==0:
                    # Found vertex ip1!
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=ip1)
                # index, edge and vertex prior to v0
                im1=(i+n-1)%n
                e=V((self._x[im1],self._y[im1]))
                vm1=v0-e
                w=wedge_product(e,point-vm1)
                if w<0:
                    return PolygonPosition(PolygonPosition.OUTSIDE)
                if w==0:
                    # Found vertex i!
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=i)
                # Otherwise we found the interior of edge i:
                return PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)
        # Loop terminated (on positive side of each edge)
        return PolygonPosition(PolygonPosition.INTERIOR)

    def flow(self,point,holonomy,translation=None):
        r"""
        Flow a point in the direction of holonomy for the length of the holonomy, or until the point leaves the polygon.
        Note that ValueErrors may be thrown if the point is not in the polygon, or if it is on the boundary and the 
        holonomy does not point into the polygon.

        INPUT:

        - ``point`` -- a point in the closure of the polygon (vector over the underlying base_ring())

        - ``holonomy`` -- direction and magnitude of motion (vector over the underlying base_ring())

        - ``translation`` -- optional translation to applied to the polygon (vector over the underlying base_ring())

        OUTPUT:

        - The point within the polygon where the motion stops (or leaves the polygon)

        - The amount of holonomy left to flow

        - a PolygonPosition object representing the combinatorial position of the stopping point

        EXAMPLES::

            sage: s=square()
            sage: V=s.parent().vector_space()
            sage: p=V((1/2,1/2))
            sage: w=V((2,0))
            sage: s.flow(p,w)
            ((1, 1/2), (3/2, 0), point positioned on interior of edge 1 of polygon)
        """
        V = self.parent().vector_space()
        if holonomy == V.zero():
            # not flowing at all!
            return point, V.zero(), self.get_point_position(point,translation=translation)
        from sage.matrix.constructor import matrix
        if translation is None:
            v0=V.zero()
        else:
            v0=translation
        w=holonomy
        for i in range(self.num_edges()):
            e=self.edge(i)
            #print "i="+str(i)+" e="+str(e)
            m=matrix([[e[0], -holonomy[0]],[e[1], -holonomy[1]]])
            #print "m="+str(m)
            #print "diff="+str(point-v0)
            try:
                ret=m.inverse()*(point-v0)
                s=ret[0]
                t=ret[1]
                #print "s="+str(s)+" and t="+str(t)
                # What if the matrix is non-invertible?

                # s is location it intersects on edge, t is the portion of the holonomy to reach this intersection
                if t>0 and 0<=s and s<=1:
                    # The ray passes through edge i.
                    if t>1:
                        # the segment from point with the given holonomy stays within the polygon
                        return point+holonomy, V.zero(), PolygonPosition(PolygonPosition.INTERIOR)
                    if s==1:
                        # exits through vertex i+1
                        v0=v0+e
                        return v0, point+holonomy-v0, PolygonPosition(PolygonPosition.VERTEX, vertex= (i+1)%self.num_edges())
                    if s==0:
                        # exits through vertex i
                        return v0, point+holonomy-v0, PolygonPosition(PolygonPosition.VERTEX, vertex= i)
                        # exits through vertex i
                    # exits through interior of edge i
                    prod=t*holonomy
                    return point+prod, holonomy-prod, PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)
            except ZeroDivisionError:
                # can safely ignore this error. It means that the edge and the holonomy are parallel.
                pass
            v0=v0+e
        # Our loop has terminated. This can mean one of several errors...
        pos = self.get_point_position(point,translation=translation)
        if pos.is_outside():
            raise ValueError("Started with point outside polygon")
        raise ValueError("Point on boundary of polygon and holonomy not pointed into the polygon.")

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
        Plot the polygon with the origin at ``translation``.
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

        #if not AA.has_coerce_map_from(field):
        #    raise ValueError("the field must have a coercion to AA")
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

def rectangle(width,height):
    F=width.parent()
    return Polygons(F)([(width,F(0)),(F(0),height),(-width,F(0)),(F(0),-height)])

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


