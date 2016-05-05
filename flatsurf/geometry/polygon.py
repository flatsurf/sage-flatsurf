r"""
Convex polygons in the plane (R^2)

This file implements convex polygons with

 - action of matrices in GL^+(2,R)
 - conversion between ground fields

EXAMPLES::

    sage: from flatsurf.geometry.polygon import polygons

    sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
    sage: p = polygons((1,0), (-sqrt2,1+sqrt2), (sqrt2-1,-1-sqrt2))
    sage: p
    Polygon: (0, 0), (1, 0), (-sqrt2 + 1, sqrt2 + 1)

    sage: M = MatrixSpace(K,2)
    sage: m = M([[1,1+sqrt2],[0,1]])
    sage: m * p
    Polygon: (0, 0), (1, 0), (sqrt2 + 4, sqrt2 + 1)
"""


import operator

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.rings.qqbar import AA
from sage.modules.free_module import VectorSpace
from sage.categories.action import Action

from flatsurf.geometry.matrix_2x2 import angle

# we implement action of GL(2,K) on polygons

ZZ_0=ZZ.zero()
ZZ_2=ZZ(2)

def wedge_product(v,w):
    return v[0]*w[1]-v[1]*w[0]

def is_same_direction(v,w):
    if wedge_product(v,w)!=ZZ_0:
        return False
    return v[0]*w[0]>0 or v[1]*w[1]>0

def is_opposite_direction(v,w):
    if wedge_product(v,w)!=ZZ_0:
        return False
    return v[0]*w[0]<0 or v[1]*w[1]<0

#def is_same_direction(v,w):
#    return v and w and \
#           not wedge_product(v,w) and \
#           (v[0]*w[0] > 0 or v[1]*w[1] > 0)
#
#def is_opposite_direction(v,w):
#    return v and w and \
#           not wedge_product(v,w) and \
#           v[0]*w[0] < 0 or v[1]*w[1] < 0

class MatrixActionOnPolygons(Action):
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

    def is_in_boundary(self):
        r"""
        Return true if the position is in the boundary of the polygon
        (either the interior of an edge or a vertex).
        """
        return self._position_type == PolygonPosition.EDGE_INTERIOR or \
            self._position_type == PolygonPosition.VERTEX

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

class ConvexPolygon(Element):
    r"""
    A convex polygon in the plane RR^2 defined up to translation.
    """
    def __init__(self, parent, edges=None, vertices=None, translation=None):
        r"""
        To construct the polygon you should either use a list of edge vectors 
        or a list of vertices. Using both will result in a ValueError. The polygon
        needs to be convex with postively oriented boundary.
        INPUT:

        - ``parent`` -- a parent
        - ``edges`` -- a list of vectors or couples whose sum is zero
        - ``vertices`` -- a list of vertices of the polygon
        - ``translation`` -- a vector representing an addition translation to be applied to the resulting polygon
        """
        Element.__init__(self, parent)
        V = parent.vector_space()
        if translation is None:
            t = V.zero()
        else:
            t = V(translation)
        if edges is not None:
            if vertices is not None:
                raise ValueError("both edges and vertices are defined which is not allowed")
            self._v = [V.zero()]
            total=t
            for i in range(len(edges)-1):
                total += V(edges[i])
                self._v.append(total)
            # Linear Time Sanity Checks
            total += V( edges[len(edges)-1] )
            if total != V.zero():
                raise ValueError("the sum over the edges do not sum up to 0")
        elif vertices is not None:
            self._v = [V(x)+t for x in vertices]
        else:
            raise ValueError("edges and vertices can't both be None")
        # Make the polgon immutable:
        self._v = tuple(self._v)
        self._convexity_check()
        for vv in self._v:
            vv.set_immutable()

    def __hash__(self):
        return hash(self._v)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from flatsurf.geometry.polygon import polygons
            sage: p1 = polygons.square()
            sage: p2 = polygons((1,0),(0,1),(-1,0),(0,-1), ring=QQbar)
            sage: p1 == p2
            True

            sage: p3 = polygons((2,0),(-1,1),(-1,-1))
            sage: p1 == p3
            False
        """
        if not isinstance(other, ConvexPolygon):
            raise TypeError
        return self._v == other._v

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from flatsurf.geometry.polygon import polygons
            sage: p1 = polygons.square()
            sage: p2 = polygons((1,0),(0,1),(-1,0),(0,-1), ring=QQbar)
            sage: p1 != p2
            False

            sage: p3 = polygons((2,0),(-1,1),(-1,-1))
            sage: p1 != p3
            True
        """
        if not isinstance(other, ConvexPolygon):
            raise TypeError
        return self._v != other._v

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
        return len(self._v)

    def _repr_(self):
        r"""
        String representation.
        """
        return "Polygon: " + ", ".join(map(str,self.vertices()))

    def vector_space(self):
        r"""Return the vector space containing the vertices."""
        return self.parent().vector_space()

    def vertices(self, translation=None):
        r"""
        Return the set of vertices as vectors.
        """
        if translation is not None:
            raise RuntimeError("the 'translation' argument is ignored and should not be used")
        return self._v

    def vertex(self,index):
        return self._v[index % self.num_edges()]


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

            sage: from flatsurf.geometry.polygon import polygons
            sage: s = polygons.square()
            sage: V = s.parent().vector_space()
            sage: s.get_point_position(V((1/2,1/2)))
            point positioned in interior of polygon
            sage: s.get_point_position(V((1,0)))
            point positioned on vertex 1 of polygon
            sage: s.get_point_position(V((1,1/2)))
            point positioned on interior of edge 1 of polygon
            sage: s.get_point_position(V((1,3/2)))
            point positioned outside polygon
        """
        V = self.vector_space()
        if translation is None:
            # Since we allow the initial vertex to be non-zero, this changed:
            v1=self.vertex(0)
        else:
            # Since we allow the initial vertex to be non-zero, this changed:
            v1=translation+self.vertex(0)
        # Below, we only make use of edge vectors:
        for i in range(self.num_edges()):
            v0=v1
            e=self.edge(i)
            v1=v0+e
            w=wedge_product(e,point-v0)
            if w < 0:
                return PolygonPosition(PolygonPosition.OUTSIDE)
            if w == 0:
                # Lies on the line through edge i!
                n=self.num_edges()
                # index and edge after v1 
                ip1=(i+1)%n
                e=self.edge(ip1)
                w=wedge_product(e,point-v1)
                if w<0:
                    return PolygonPosition(PolygonPosition.OUTSIDE)
                if w==0:
                    # Found vertex ip1!
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=ip1)
                # index, edge and vertex prior to v0
                im1=(i+n-1)%n
                e=self.edge(im1)
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

    def flow_to_exit(self,point,direction):
        r"""
        Flow a point in the direction of holonomy until the point leaves the
        polygon.  Note that ValueErrors may be thrown if the point is not in the
        polygon, or if it is on the boundary and the holonomy does not point
        into the polygon.

        INPUT:

        - ``point`` -- a point in the closure of the polygon (as a vector)

        - ``holonomy`` -- direction of motion (a vector of non-zero length)

        OUTPUT:

        - The point in the boundary of the polygon where the trajectory exits

        - a PolygonPosition object representing the combinatorial position of the stopping point
        """
        V = self.parent().vector_space()
        if direction == V.zero():
            raise ValueError("Zero vector provided as direction.")
        v0=self.vertex(0)
        w=direction
        from sage.matrix.constructor import matrix
        for i in range(self.num_edges()):
            e=self.edge(i)
            #print "i="+str(i)+" e="+str(e)+" direction="+str(direction)
            m=matrix([[e[0], -direction[0]],[e[1], -direction[1]]])
            try:
                ret=m.inverse()*(point-v0)
                s=ret[0]
                t=ret[1]
                #print "s="+str(s)+" and t="+str(t)
                # What if the matrix is non-invertible?
                
                # Answer: You'll get a ZeroDivisionError which means that the edge is parallel
                # to the direction.

                # s is location it intersects on edge, t is the portion of the direction to reach this intersection
                if t>0 and 0<=s and s<=1:
                    # The ray passes through edge i.
                    if s==1:
                        # exits through vertex i+1
                        v0=v0+e
                        return v0, PolygonPosition(PolygonPosition.VERTEX, vertex= (i+1)%self.num_edges())
                    if s==0:
                        # exits through vertex i
                        return v0, PolygonPosition(PolygonPosition.VERTEX, vertex= i)
                        # exits through vertex i
                    # exits through interior of edge i
                    prod=t*direction
                    return point+prod, PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)
            except ZeroDivisionError:
                # Here we know the edge and the direction are parallel
                if wedge_product(e,point-v0)==0:
                    # In this case point lies on the edge. 
                    # We need to work out which direction to move in.
                    if is_same_direction(e,point-v0):
                        # exits through vertex i+1
                        v0=v0+e
                        return v0, PolygonPosition(PolygonPosition.VERTEX, vertex= (i+1)%self.num_edges())
                    else:
                        # exits through vertex i
                        return v0, PolygonPosition(PolygonPosition.VERTEX, vertex= i)
                pass
            v0=v0+e
        # Our loop has terminated. This can mean one of several errors...
        pos = self.get_point_position(point)
        if pos.is_outside():
            raise ValueError("Started with point outside polygon")
        raise ValueError("Point on boundary of polygon and direction not pointed into the polygon.")


    def flow(self,point,holonomy,translation=None):
        r"""
        Flow a point in the direction of holonomy for the length of the
        holonomy, or until the point leaves the polygon.  Note that ValueErrors
        may be thrown if the point is not in the polygon, or if it is on the
        boundary and the holonomy does not point into the polygon.

        INPUT:

        - ``point`` -- a point in the closure of the polygon (vector over the underlying base_ring())

        - ``holonomy`` -- direction and magnitude of motion (vector over the underlying base_ring())

        - ``translation`` -- optional translation to applied to the polygon (vector over the underlying base_ring())

        OUTPUT:

        - The point within the polygon where the motion stops (or leaves the polygon)

        - The amount of holonomy left to flow

        - a PolygonPosition object representing the combinatorial position of the stopping point

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: s = polygons.square()
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
            v0=self.vertex(0)
        else:
            v0=self.vertex(0)+translation
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
        Return the list of edges as vectors.
        """
        return [self.edge(i) for i in range(self.num_edges())]

    def edge(self, i):
        r"""
        Return a vector representing the ``i``-th edge of the polygon.
        """
        return self.vertex(i+1)-self.vertex(i)

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
        return point2d(P, color='red') + line2d(P + (P[0],), color='orange') + polygon2d(P, alpha=0.3)

    def angle(self, e):
        r"""
        Return the angle at the begining of the start point of the edge ``e``.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: polygons.square().angle(0)
            1/4
            sage: polygons.regular_ngon(8).angle(0)
            3/8
        """
        return angle(self.edge(e), - self.edge((e-1)%self.num_edges()))

    def area(self):
        r"""
        Return the area of this polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: polygons.regular_ngon(8).area()
            2*a + 2
            sage: _ == 2*AA(2).sqrt() + 2
            True

            sage: AA(polygons.regular_ngon(11).area())
            9.36563990694544?

            sage: polygons.square().area()
            1
            sage: (2*polygons.square()).area()
            4
        """
        # Will use an area formula obtainable from Green's theorem. See for instance:
        # http://math.blogoverflow.com/2014/06/04/greens-theorem-and-area-of-polygons/
        total = self.field().zero()
        for i in range(self.num_edges()):
            total += (self.vertex(i)[0]+self.vertex(i+1)[0])*self.edge(i)[1]
        return total/ZZ_2


class ConvexPolygons(Parent):
    Element = ConvexPolygon
    def __init__(self, field):
        Parent.__init__(self, category=Sets())
        self._field = field
        self.register_action(MatrixActionOnPolygons(self))

    def has_coerce_map_from(self, other):
        return isinstance(other, Polygons) and self.field().has_coerce_map_from(other.field())

    def _an_element_(self):
        return self([(1,0),(0,1),(-1,0),(0,-1)])

    def base_ring(self):
        return self._field

    field = base_ring

    def _repr_(self):
        return "polygons with coordinates in %s"%self.base_ring()

    @cached_method
    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), 2)

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

Polygons = ConvexPolygons

def number_field_elements_from_algebraics(elts, name='a'):
    r"""
    The native Sage function ``number_field_elements_from_algebraics`` currently
    returns number field *without* embedding. This function return field with
    embedding!

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import number_field_elements_from_algebraics
        sage: z = QQbar.zeta(5)
        sage: c = z.real()
        sage: s = z.imag()
        sage: number_field_elements_from_algebraics((c,s))
        (Number Field in a with defining polynomial y^4 - 5*y^2 + 5,
         [1/2*a^2 - 3/2, 1/2*a])
    """
    from sage.rings.qqbar import number_field_elements_from_algebraics
    from sage.rings.number_field.number_field import NumberField
    field,elts,phi = number_field_elements_from_algebraics(elts, minimal=True)

    polys = [x.polynomial() for x in elts]
    K = NumberField(field.polynomial(), name, embedding=AA(phi(field.gen())))
    gen = K.gen()

    return K, [x.polynomial()(gen) for x in elts]

class PolygonsConstructor:
    def square(self, side=1, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: polygons.square()
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: polygons.square(field=QQbar).parent()
            polygons with coordinates in Algebraic Field
        """
        return self.rectangle(side,side,**kwds)

    def rectangle(self, width, height, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: polygons.rectangle(1,2)
            Polygon: (0, 0), (1, 0), (1, 2), (0, 2)

            sage: K.<sqrt2> = QuadraticField(2)
            sage: polygons.rectangle(1,sqrt2)
            Polygon: (0, 0), (1, 0), (1, sqrt2), (0, sqrt2)
            sage: _.parent()
            polygons with coordinates in Number Field in sqrt2 with defining
            polynomial x^2 - 2
        """
        return self((width,0),(0,height),(-width,0),(0,-height), **kwds)

    @staticmethod
    def regular_ngon(n):
        r"""
        Return a regular n-gon.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: p = polygons.regular_ngon(17)
            sage: p
            Polygon: (0, 0), (1, 0), ..., (-1/2*a^14 + 15/2*a^12 - 45*a^10 + 275/2*a^8 - 225*a^6 + 189*a^4 - 70*a^2 + 15/2, 1/2*a)
        """
        from sage.rings.qqbar import QQbar

        c = QQbar.zeta(n).real()
        s = QQbar.zeta(n).imag()

        field, (c,s) = number_field_elements_from_algebraics((c,s))

        cn = field.one()
        sn = field.zero()
        edges = [(cn,sn)]
        for _ in range(n-1):
            cn,sn = c*cn - s*sn, c*sn + s*cn
            edges.append((cn,sn))

        return Polygons(field)(edges)

    def __call__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: polygons((1,0),(0,1),(-1,0),(0,-1))
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: polygons((1,0),(0,1),(-1,0),(0,-1), ring=QQbar)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: _.parent()
            polygons with coordinates in Algebraic Field
        """
        base_ring = None
        if 'ring' in kwds:
            base_ring = kwds.pop('ring')
        if 'base_ring' in kwds:
            base_ring = kwds.pop('base_ring')
        if 'field' in kwds:
            base_ring = kwds.pop('field')

        if base_ring is None:
            from sage.structure.sequence import Sequence
            from sage.modules.free_module_element import vector

            s = Sequence(map(vector, args))
            V = s.universe()
            base_ring = V.base_ring()
        else:
            from sage.modules.free_module import VectorSpace
            V = VectorSpace(base_ring,2)
            s = map(V, args)

        if base_ring not in Fields():
            base_ring = base_ring.fraction_field()

        return ConvexPolygons(base_ring)(s, **kwds)

polygons = PolygonsConstructor()

def regular_octagon(field=None):
    from sage.misc.superseded import deprecation
    deprecation(33, "Do not use this function anymore but regular_ngon")
    return polygons.regular_ngon(8)

class PolygonCreator():
    r"""
    Class for iteratively constructing a polygon over the field.
    """
    def __init__(self, field = QQ):
        r"""Create a polygon in the provided field."""
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


