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

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

import operator

from sage.misc.cachefunc import cached_method

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields
from sage.categories.action import Action

from sage.rings.all import ZZ, QQ, AA, RR

from sage.modules.free_module_element import vector
from sage.modules.free_module_element import free_module_element
from sage.matrix.constructor import matrix
from sage.rings.polynomial.polynomial_ring import polygen

from .matrix_2x2 import angle

# we implement action of GL(2,K) on polygons

ZZ_0 = ZZ.zero()
ZZ_2 = ZZ(2)

def dot_product(v,w):
    return v[0]*w[0]+v[1]*w[1]

def wedge_product(v,w):
    return v[0]*w[1]-v[1]*w[0]

def wedge(u, v):
    r"""
    General wedge product of two vectors.
    """
    d = len(u)
    R = u.base_ring()
    assert len(u) == len(v) and v.base_ring() == R
    return free_module_element(R, d*(d-1)//2, [(u[i]*v[j] - u[j]*v[i]) for i in range(d-1) for j in range(i+1,d)])

def tensor(u, v):
    r"""
    General tensor product of two vectors.
    """
    d = len(u)
    R = u.base_ring()
    assert len(u) == len(v) and v.base_ring() == R
    return matrix(R, d, [u[i]*v[j] for j in range(d) for i in range(d)])

def line_intersection(p1,p2,q1,q2):
    r"""
    Return the point of intersection between the line joining p1 to p2
    and the line joining q1 to q2. If the lines are parallel we return
    None. Here p1, p2, q1 and q2 should be vectors in the plane.
    """
    if wedge_product(p2-p1,q2-q1) == 0:
        return None
    # Since the wedge product is non-zero, the following is invertible:
    m=matrix([[p2[0]-p1[0], q1[0]-q2[0]],[p2[1]-p1[1], q1[1]-q2[1]]])
    return p1+(m.inverse()*(q1-p1))[0] * (p2-p1)

def is_same_direction(v,w,zero=None):
    r"""
    EXAMPLES::

        sage: from flatsurf.geometry.polygon import is_same_direction
        sage: V = QQ**2

        sage: is_same_direction(V((0,1)), V((0,2)))
        True
        sage: is_same_direction(V((1,-1)), V((2,-2)))
        True
        sage: is_same_direction(V((4,-2)), V((2,-1)))
        True
        sage: is_same_direction(V((1,2)), V((2,4)))
        True
        sage: is_same_direction(V((0,2)), V((0,1)))
        True

        sage: is_same_direction(V((1,1)), V((1,2)))
        False
        sage: is_same_direction(V((1,2)), V((2,1)))
        False
        sage: is_same_direction(V((1,2)), V((1,-2)))
        False
        sage: is_same_direction(V((1,2)), V((-1,-2)))
        False
        sage: is_same_direction(V((2,-1)), V((-2,1)))
        False

        sage: is_same_direction(V((1,0)), V.zero())
        Traceback (most recent call last):
        ...
        TypeError: zero vector has no direction

        sage: for _ in range(100):
        ....:    v = V.random_element()
        ....:    if not v: continue
        ....:    assert is_same_direction(v, 2*v)
        ....:    assert not is_same_direction(v, -v)
    """
    if not v or not w:
        raise TypeError("zero vector has no direction")
    return not wedge_product(v,w) and (v[0]*w[0] > 0 or v[1]*w[1] > 0)

def is_opposite_direction(v,w):
    r"""
    EXAMPLES::

        sage: from flatsurf.geometry.polygon import is_opposite_direction
        sage: V = QQ**2

        sage: is_opposite_direction(V((0,1)), V((0,-2)))
        True
        sage: is_opposite_direction(V((1,-1)), V((-2,2)))
        True
        sage: is_opposite_direction(V((4,-2)), V((-2,1)))
        True
        sage: is_opposite_direction(V((-1,-2)), V((2,4)))
        True

        sage: is_opposite_direction(V((1,1)), V((1,2)))
        False
        sage: is_opposite_direction(V((1,2)), V((2,1)))
        False
        sage: is_opposite_direction(V((0,2)), V((0,1)))
        False
        sage: is_opposite_direction(V((1,2)), V((1,-2)))
        False
        sage: is_opposite_direction(V((1,2)), V((-1,2)))
        False
        sage: is_opposite_direction(V((2,-1)), V((-2,-1)))
        False

        sage: is_opposite_direction(V((1,0)), V.zero())
        Traceback (most recent call last):
        ...
        TypeError: zero vector has no direction

        sage: for _ in range(100):
        ....:    v = V.random_element()
        ....:    if not v: continue
        ....:    assert not is_opposite_direction(v, v)
        ....:    assert not is_opposite_direction(v,2*v)
        ....:    assert is_opposite_direction(v, -v)
    """
    if not v or not w:
        raise TypeError("zero vector has no direction")
    return not wedge_product(v,w) and (v[0]*w[0] < 0 or v[1]*w[1] < 0)

def solve(x,u,y,v):
    r"""
    Return (a,b) so that: x + au = y + bv

    INPUT:

    - ``x``, ``u``, ``y``, ``v`` -- two dimensional vectors

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import solve
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: V = VectorSpace(K,2)
        sage: x = V((1,-sqrt2))
        sage: y = V((1,1))
        sage: a = V((0,1))
        sage: b = V((-sqrt2, sqrt2+1))
        sage: u = V((0,1))
        sage: v = V((-sqrt2, sqrt2+1))
        sage: a, b = solve(x,u,y,v)
        sage: x + a*u == y + b*v
        True

        sage: u = V((1,1))
        sage: v = V((1,sqrt2))
        sage: a, b = solve(x,u,y,v)
        sage: x + a*u == y + b*v
        True
    """
    d = -u[0] * v[1] + u[1] * v[0]
    if d.is_zero():
        raise ValueError("parallel vectors")
    a = v[1] * (x[0]-y[0]) + v[0] * (y[1] - x[1])
    b = u[1] * (x[0]-y[0]) + u[0] * (y[1] - x[1])
    return (a/d, b/d)

class MatrixActionOnPolygons(Action):
    def __init__(self, polygons):
        from sage.matrix.matrix_space import MatrixSpace
        K = polygons.field()
        Action.__init__(self, MatrixSpace(K,2), polygons, True, operator.mul)

    def _act_(self, g, x):
        r"""
        Apply the 2x2 matrix g to the polygon.
        The matrix must have non-zero determinant. If the determinant is negative, then the vertices and edges
        are relabeled according to the involutions $v \mapsto (n-v)%n$ and  $e \mapsto n-1-e$ respectively.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import Polygons, ConvexPolygon
            sage: P=Polygons(QQ)
            sage: p=ConvexPolygon(P,[(1,0),(0,1),(-1,-1)])
            sage: print(p)
            Polygon: (1, 0), (0, 1), (-1, -1)
            sage: r=matrix(ZZ,[[0,1],[1,0]])
            sage: print(r*p)
            Polygon: (0, 1), (-1, -1), (1, 0)
        """
        det = g.det()
        if det > 0:
            return x.parent()(vertices=[g*v for v in x.vertices()])
        if det < 0:
            # Note that in this case we reverse the order
            vertices=[g*x.vertex(0)]
            for i in range(x.num_edges()-1,0,-1):
                vertices.append(g*x.vertex(i))
            return x.parent()(vertices=vertices)
        raise ValueError("Can not act on a polygon with matrix with zero determinant")


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
        if not self.is_in_edge_interior():
            raise ValueError("Asked for edge when not in edge interior.")
        return self._edge

    def get_vertex(self):
        if not self.is_vertex():
            raise ValueError("Asked for vertex when not a vertex.")
        return self._vertex

class ConvexPolygon(Element):
    r"""
    A convex polygon in the plane RR^2
    """
    def __init__(self, parent, vertices):
        r"""
        To construct the polygon you should either use a list of edge vectors
        or a list of vertices. Using both will result in a ValueError. The polygon
        needs to be convex with postively oriented boundary.

        INPUT:

        - ``parent`` -- a parent

        - ``vertices`` -- a list of vertices of the polygon
        """
        Element.__init__(self, parent)

        V = parent.vector_space()
        self._v = tuple(map(V, vertices))
        for vv in self._v: vv.set_immutable()
        self._convexity_check()

    def __hash__(self):
        # Apparently tuples do not cache their hash!
        if not hasattr(self,"_hash"):
            self._hash = hash(self._v)
        return self._hash

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
        return isinstance(other, ConvexPolygon) and self._v == other._v

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
        return not isinstance(other, ConvexPolygon) or self._v != other._v

    def __cmp__(self, other):
        if not isinstance(other,ConvexPolygon):
            raise ValueError("__cmp__ only implemented for ConvexPolygons")
        if not self.parent().base_ring()==other.parent().base_ring():
            raise ValueError("__cmp__ only implemented for ConvexPolygons defined over the same base_ring")
        sign = self.num_edges() - other.num_edges()
        if sign>0:
            return 1
        if sign<0:
            return -1
        sign = self.area() - other.area()
        if sign>self.base_ring().zero():
            return 1
        if sign<self.base_ring().zero():
            return -1
        for v in range(1,self.num_edges()):
            p=self.vertex(v)
            q=other.vertex(v)
            sign = p[0]-q[0]
            if sign>self.base_ring().zero():
                return 1
            if sign<self.base_ring().zero():
                return -1
            sign = p[1]-q[1]
            if sign>self.base_ring().zero():
                return 1
            if sign<self.base_ring().zero():
                return -1
        return 0

    def is_strictly_convex(self):
        r"""
        Check whether the polygon is strictly convex

        EXAMPLES::

            sage: from flatsurf import *
            sage: polygons(vertices=[(0,0), (1,0), (1,1)]).is_strictly_convex()
            True
            sage: polygons(vertices=[(0,0), (1,0), (2,0), (1,1)]).is_strictly_convex()
            False
        """
        for i in range(self.num_edges()):
            if wedge_product(self.edge(i), self.edge(i+1)).is_zero():
                return False
        return True

    def _convexity_check(self):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: polygons(vertices=[(0,0),(1,0)])
            Traceback (most recent call last):
            ...
            ValueError: a polygon should have more than two edges!
            sage: polygons(vertices=[(0,0),(1,2),(0,1),(-1,2)])
            Traceback (most recent call last):
            ...
            ValueError: not convex
            sage: polygons(vertices=[(0,0),(1,0),(2,0)])
            Traceback (most recent call last):
            ...
            ValueError: degenerate polygon
        """
        if self.num_edges() <= 2:
            raise ValueError("a polygon should have more than two edges!")

        if not sum(self.edges()).is_zero():
            raise ValueError("the sum over the edges do not sum up to 0")

        for i in range(self.num_edges()):
            if wedge_product(self.edge(i), self.edge(i+1)) < 0:
                raise ValueError("not convex")
            if is_opposite_direction(self.edge(i), self.edge(i+1)):
                raise ValueError("degenerate polygon")

    def find_separatrix(self, direction=None, start_vertex=0):
        r"""
        Returns a pair (v,same) where v is a vertex and dir is a boolean.
        The provided parameter "direction" should be a non-zero vector with
        two entries, or by default direction=(0,1).

        A separatrix is a ray leaving a vertex and entering the polygon.

        The vertex v will have a separatrix leaving it which is parallel to
        direction. The returned value "same" answers the question if this separatrix
        points in the same direction as "direction". There is a boundary case:
        we allow the separatrix to be an edge if and only if traveling along
        the sepatrix from the vertex would travel in a counter-clockwise
        direction about the polygon.

        The vertex returned is uniquely defined from the above if the polygon
        is a triangle. Otherwise, we return the first vertex with this property
        obtained by inspecting starting at start_vertex (defaults to 0) and
        then moving in the counter-clockwise direction.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: p=polygons.square()
            sage: print(p.find_separatrix())
            (1, True)
            sage: print(p.find_separatrix(start_vertex=2))
            (3, False)
        """
        if direction is None:
            direction = self.vector_space()((self.base_ring().zero(), self.base_ring().one()))
        else:
            assert not direction.is_zero()
        v=start_vertex
        n=self.num_edges()
        zero=self.base_ring().zero()
        for i in range(self.num_edges()):
            if wedge_product(self.edge(v),direction) >= zero and \
                wedge_product(self.edge(v+n-1),direction) > zero:
                return v,True
            if wedge_product(self.edge(v),direction) <= zero and \
                wedge_product(self.edge(v+n-1),direction) < zero:
                return v,False
            v=v+1%n
        raise RuntimeError("Failed to find a separatrix")

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
        r"""
        Return the vector space containing the vertices.
        """
        return self.parent().vector_space()

    def vertices(self, translation=None):
        r"""
        Return the set of vertices as vectors.
        """
        if translation is None:
            return self._v

        translation = self.parent().vector_space()(translation)
        return [t+v for v in self.vertices()]

    def vertex(self, i):
        r"""
        Return the ``i``-th vertex as a vector
        """
        return self._v[i % len(self._v)]

    def __iter__(self):
        return iter(self.vertices())

    def contains_point(self, point, translation=None):
        r"""
        Return true if the point is within the polygon (after the polygon is possibly translated)
        """
        return self.get_point_position(point,translation=translation).is_inside()

    def get_point_position(self, point, translation=None):
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

            sage: p=polygons(edges=[(1,0),(1,0),(1,0),(0,1),(-3,0),(0,-1)])
            sage: V=p.vector_space()
            sage: p.get_point_position(V([10,0]))
            point positioned outside polygon
            sage: p.get_point_position(V([1/2,0]))
            point positioned on interior of edge 0 of polygon
            sage: p.get_point_position(V([3/2,0]))
            point positioned on interior of edge 1 of polygon
            sage: p.get_point_position(V([2,0]))
            point positioned on vertex 2 of polygon
            sage: p.get_point_position(V([5/2,0]))
            point positioned on interior of edge 2 of polygon
            sage: p.get_point_position(V([5/2,1/4]))
            point positioned in interior of polygon
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
                dp1 = dot_product(e,point-v0)
                if dp1 == 0:
                    return PolygonPosition(PolygonPosition.VERTEX, vertex=i)
                dp2 = dot_product(e,e)
                if 0 < dp1 and dp1 < dp2:
                    return PolygonPosition(PolygonPosition.EDGE_INTERIOR, edge=i)
        # Loop terminated (on inside of each edge)
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
        for i in range(self.num_edges()):
            e=self.edge(i)
            m=matrix([[e[0], -direction[0]],[e[1], -direction[1]]])
            try:
                ret=m.inverse()*(point-v0)
                s=ret[0]
                t=ret[1]
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
                    if (point-v0).is_zero() or is_same_direction(e,point-v0):
                        # exits through vertex i+1
                        return self.vertex(i+1), PolygonPosition(PolygonPosition.VERTEX, vertex= (i+1)%self.num_edges())
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

    def flow_map(self, direction):
        r"""
        Return a polygonal map associated to the flow in ``direction`` in this
        polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: S = polygons(vertices=[(0,0),(2,0),(2,2),(1,2),(0,2),(0,1)])
            sage: S.flow_map((0,1))
             Flow polygon map:
              3 2
              0
             top lengths: [1, 1]
             bot lengths: [2]
            sage: S.flow_map((1,1))
            Flow polygon map:
             3 2 1
             4 5 0
            top lengths: [1, 1, 2]
            bot lengths: [1, 1, 2]
            sage: S.flow_map((-1,-1))
            Flow polygon map:
             0 5 4
             1 2 3
            top lengths: [2, 1, 1]
            bot lengths: [2, 1, 1]

            sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
            sage: S.flow_map((sqrt2,1))
            Flow polygon map:
             3 2 1
             4 5 0
            top lengths: [1, 1, 2*sqrt2]
            bot lengths: [sqrt2, sqrt2, 2]
        """
        direction = vector(direction)
        DP = direction.parent()
        P = self.vector_space()
        if DP != P:
            from sage.structure.element import get_coercion_model
            cm = get_coercion_model()
            P = cm.common_parent(DP,P)
            ring = P.base_ring()
            direction = direction.change_ring(ring)
        else:
            ring = P.base_ring()

        # first compute the transversal length of each edge
        t = P([direction[1], -direction[0]])
        lengths = [t.dot_product(e) for e in self.edges()]
        n = len(lengths)
        for i in range(n):
            j = (i+1)%len(lengths)
            l0 = lengths[i]
            l1 = lengths[j]
            if l0 >= 0 and l1 <  0: rt = j
            if l0 >  0 and l1 <= 0: rb = j
            if l0 <= 0 and l1 >  0: lb = j
            if l0 <  0 and l1 >= 0: lt = j

        if rt < lt:
            top_lengths = lengths[rt:lt]
            top_labels = list(range(rt,lt))
        else:
            top_lengths = lengths[rt:] + lengths[:lt]
            top_labels = list(range(rt,n)) + list(range(lt))
        top_lengths = [-x for x in reversed(top_lengths)]
        top_labels.reverse()

        if lb < rb:
            bot_lengths = lengths[lb:rb]
            bot_labels = list(range(lb,rb))
        else:
            bot_lengths = lengths[lb:] + lengths[:rb]
            bot_labels = list(range(lb,n)) + list(range(rb))

        from .interval_exchange_transformation import FlowPolygonMap
        return FlowPolygonMap(ring, bot_labels, bot_lengths,
                                                    top_labels, top_lengths)

    def flow(self, point, holonomy, translation=None):
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
            sage: V = s.parent().vector_space()
            sage: p = V((1/2,1/2))
            sage: w = V((2,0))
            sage: s.flow(p,w)
            ((1, 1/2), (3/2, 0), point positioned on interior of edge 1 of polygon)
        """
        V = self.parent().vector_space()
        if holonomy == V.zero():
            # not flowing at all!
            return point, V.zero(), self.get_point_position(point,translation=translation)
        if translation is None:
            v0=self.vertex(0)
        else:
            v0=self.vertex(0)+translation
        w=holonomy
        for i in range(self.num_edges()):
            e=self.edge(i)
            m=matrix([[e[0], -holonomy[0]],[e[1], -holonomy[1]]])
            try:
                ret=m.inverse()*(point-v0)
                s=ret[0]
                t=ret[1]
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
        Return an iterator overt the edges
        """
        return [self.edge(i) for i in range(self.num_edges())]

    def edge(self, i):
        r"""
        Return a vector representing the ``i``-th edge of the polygon.
        """
        return self.vertex(i+1) - self.vertex(i)

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

    def circumscribing_circle(self):
        r"""
        Returns the circle which circumscribes this polygon.
        Raises a ValueError if the polygon is not circumscribed by a circle.

        EXAMPLES::

            sage: from flatsurf import *
            sage: P = polygons(vertices=[(0,0),(1,0),(2,1),(-1,1)])
            sage: P.circumscribing_circle()
            Circle((1/2, 3/2), 5/2)
        """
        from .circle import circle_from_three_points
        circle = circle_from_three_points(self.vertex(0), self.vertex(1), self.vertex(2), self.base_ring())
        for i in range(3,self.num_edges()):
            if not circle.point_position(self.vertex(i))==0:
                raise ValueError("Vertex "+str(i)+" is not on the circle.")
        return circle

    def j_invariant(self):
        r"""
        Return the Kenyon-Smille J-invariant of this polygon.

        The base ring of the polygon must be a number field.

        The output is a triple ``(Jxx, Jyy, Jxy)`` that corresponds
        respectively to the Sah-Arnoux-Fathi invariant of the vertical flow,
        the Sah-Arnoux-Fathi invariant of the horizontal flow and the `xy`-matrix.

        EXAMPLES::

            sage: from flatsurf import *

            sage: polygons.right_triangle(1/3,1).j_invariant()
            (
                      [0 0]
            (0), (0), [1 0]
            )

        The regular 8-gon::

            sage: polygons.regular_ngon(8).j_invariant()
            (
                      [2 2]
            (0), (0), [2 1]
            )

            (
                           [  0 3/2]
            (1/2), (-1/2), [3/2   0]
            )

        Some extra debugging::

            sage: from flatsurf.geometry.polygon import wedge
            sage: K.<a> = NumberField(x^3 - 2, embedding=AA(2)**(1/3))
            sage: ux = 1 + a + a**2
            sage: uy = -2/3 + a
            sage: vx = 1/5 - a**2
            sage: vy = a + 7/13*a**2
            sage: p = polygons((ux, uy), (vx,vy), (-ux-vx,-uy-vy), ring=K)
            sage: Jxx, Jyy, Jxy = p.j_invariant()
            sage: wedge(ux.vector(), vx.vector()) == Jxx
            True
            sage: wedge(uy.vector(), vy.vector()) == Jyy
            True
        """
        if self.base_ring() is QQ:
            raise NotImplementedError

        K = self.base_ring()
        try:
            V, from_V, to_V = K.vector_space()
        except (AttributeError, ValueError):
            raise ValueError("the surface needs to be define over a number field")

        dim = K.degree()
        Jxx = Jyy = free_module_element(K, dim*(dim-1)//2)
        Jxy = matrix(K, dim)
        vertices = list(self.vertices())
        vertices.append(vertices[0])
        for i in range(len(vertices) - 1):
            a = to_V(vertices[i][0])
            b = to_V(vertices[i][1])
            c = to_V(vertices[i+1][0])
            d = to_V(vertices[i+1][1])
            Jxx += wedge(a, c)
            Jyy += wedge(b, d)
            Jxy += tensor(a, d)
            Jxy -= tensor(c, b)

        return (Jxx, Jyy, Jxy)

class ConvexPolygons(UniqueRepresentation, Parent):
    r"""
    The set of convex polygons with a fixed base field.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import ConvexPolygons
        sage: C = ConvexPolygons(QQ)
        sage: C(vertices=[(0,0), (2,0), (1,1)])
        Polygon: (0, 0), (2, 0), (1, 1)
        sage: C(edges=[(1,0), (0,1), (-1,0), (0,-1)])
        Polygon: (0, 0), (1, 0), (1, 1), (0, 1)

    TESTS::

        sage: ConvexPolygons(QQ) is ConvexPolygons(QQ)
        True
        sage: TestSuite(ConvexPolygons(QQ)).run()
        sage: TestSuite(ConvexPolygons(QQbar)).run()
    """
    Element = ConvexPolygon

    def __init__(self, field):
        Parent.__init__(self, category=Sets())
        if not field in Fields():
            raise ValueError("'field' must be a field")
        self._field = field
        self.register_action(MatrixActionOnPolygons(self))

    def has_coerce_map_from(self, other):
        r"""
        TESTS::

            sage: from flatsurf.geometry.polygon import ConvexPolygons
            sage: C1 = ConvexPolygons(QQ)
            sage: C2 = ConvexPolygons(AA)
            sage: C2.has_coerce_map_from(C1)
            True
            sage: C1.has_coerce_map_from(C2)
            False
        """
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
        Return the vector space in which this polygon embeds.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import ConvexPolygons
            sage: C = ConvexPolygons(QQ)
            sage: C.vector_space()
            Vector space of dimension 2 over Rational Field
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), 2)

    def _element_constructor_(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf.geometry.polygon import ConvexPolygons

            sage: C = ConvexPolygons(QQ)
            sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
            sage: p
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: C(p) is p
            True

            sage: D = ConvexPolygons(QQbar)
            sage: D(p)
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: D(vertices=p.vertices())
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: D(edges=p.edges())
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
        """
        if len(args) == 1 and isinstance(args[0], ConvexPolygon):
            a = args[0]
            if a.parent() is self:
                return a
            vertices = map(self.vector_space(), a.vertices())
            args = None

        else:
            vertices = kwds.get('vertices')
            edges = kwds.get('edges')

            if vertices is None:
                if edges is None:
                    if not args:
                        raise ValueError("need something!")
                    if len(args) == 1:
                        edges = args[0]
                    else:
                        edges = args
                if edges is not None:
                    v = self.vector_space().zero()
                    vertices = []
                    for e in map(self.vector_space(), edges):
                        vertices.append(v)
                        v += e
                else:
                    raise ValueError("either vertices or edges should be provided")

            if vertices is None and edges is None:
                raise ValueError("exactly one of 'vertices' or 'edges' should be provided")

        return self.element_class(self, vertices)

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
        (Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.902113032590308?,
         [1/2*a^2 - 3/2, 1/2*a])
    """
    # case when all elements are rationals
    if all(x in QQ for x in elts):
        return QQ, [QQ(x) for x in elts]

    # general case
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
            polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
        """
        return self((width,0),(0,height),(-width,0),(0,-height), **kwds)

    def triangle(self, a, b, c):
        """
        Return the triangle with angles a*pi/N,b*pi/N,c*pi/N where N=a+b+c.
        
        INPUT:
        
        - ``a``, ``b``, ``c`` -- integers
        
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: T = polygons.triangle(3,4,5)
            sage: T
            Polygon: (0, 0), (1, 0), (1/2*a + 3/2, 1/2*a + 3/2)
            sage: T.base_ring()
            Number Field in a with defining polynomial y^2 - 3 with a = -1.732050807568878?
        """
        from sage.rings.qqbar import QQbar
        from sage.rings.qqbar import number_field_elements_from_algebraics
        zN = QQbar.zeta(2 * (a + b + c))
        L = polygen(QQbar, 'L')
        z = zN**a
        rotated = (1 - L * z) * zN**b
        rotated_conjugate = rotated.map_coefficients(lambda z: z.conjugate())
        real = rotated - rotated_conjugate
        p = -real[0].imag() / real[1].imag() * z
        r = p.real()
        i = p.imag()
        field, (r, i), phi = number_field_elements_from_algebraics((r, i), embedded=True)
        return polygons(vertices=[(0, 0), (1, 0), (r, i)], ring=field)

    @staticmethod
    def regular_ngon(n, field=None):
        r"""
        Return a regular n-gon with unit length edges, first edge horizontal, and other vertices lying above this edge.

        Assuming field is None (by default) the polygon is defined over a NumberField (the minimal number field determined by n).
        Otherwise you can set field equal to AA to define the polygon over the Algebraic Reals. Other values for the field
        parameter will result in a ValueError.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: p = polygons.regular_ngon(17)
            sage: p
            Polygon: (0, 0), (1, 0), ..., (-1/2*a^14 + 15/2*a^12 - 45*a^10 + 275/2*a^8 - 225*a^6 + 189*a^4 - 70*a^2 + 15/2, 1/2*a)

            sage: polygons.regular_ngon(3,field=AA)
            Polygon: (0, 0), (1, 0), (1/2, 0.866025403784439?)
        """
        # The code below crashes for n=4!
        if n==4:
            return polygons.square(QQ(1), field=field)

        from sage.rings.qqbar import QQbar

        c = QQbar.zeta(n).real()
        s = QQbar.zeta(n).imag()

        if field is None:
            field, (c,s) = number_field_elements_from_algebraics((c,s))
        cn = field.one()
        sn = field.zero()
        edges = [(cn,sn)]
        for _ in range(n-1):
            cn,sn = c*cn - s*sn, c*sn + s*cn
            edges.append((cn,sn))

        return Polygons(field)(edges=edges)

    @staticmethod
    def right_triangle(angle, leg0=None, leg1=None, hypotenuse=None):
        r"""
        Return a right triangle in a number field with an angle of pi*angle.

        You can specify the length of the first leg (``leg0``), the second leg (``leg1``),
        or the ``hypotenuse``.

        EXAMPLES::

            sage: from flatsurf import *

            sage: P = polygons.right_triangle(1/3, 1)
            sage: P
            Polygon: (0, 0), (1, 0), (1, a)
            sage: P.base_ring()
            Number Field in a with defining polynomial y^2 - 3 with a = 1.732050807568878?

            sage: polygons.right_triangle(1/4,1)
            Polygon: (0, 0), (1, 0), (1, 1)
            sage: polygons.right_triangle(1/4,1).field()
            Rational Field
        """
        from sage.rings.qqbar import QQbar

        angle = QQ(angle)
        if angle <= 0 or angle > QQ((1,2)):
            raise ValueError('angle must be in ]0,1/2]')

        z = QQbar.zeta(2*angle.denom())**angle.numer()
        c = z.real()
        s = z.imag()

        nargs = (leg0 is not None) + (leg1 is not None) + (hypotenuse is not None)

        if nargs == 0:
            leg0 = 1
        elif nargs > 1:
            raise ValueError('only one length can be specified')

        if leg0 is not None:
            c,s = leg0*c/c, leg0*s/c
        elif leg1 is not None:
            c,s = leg1*c/s, leg1*s/s
        elif hypotenuse is not None:
            c,s = hypotenuse*c, hypotenuse*s

        field, (c,s) = number_field_elements_from_algebraics((c,s))
        return Polygons(field)(edges=[(c,field.zero()),(field.zero(),s),(-c,-s)])

    def __call__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf import *

            sage: polygons((1,0),(0,1),(-1,0),(0,-1))
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: polygons((1,0),(0,1),(-1,0),(0,-1), ring=QQbar)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: _.parent()
            polygons with coordinates in Algebraic Field

            sage: polygons(vertices=[(0,0), (1,0), (0,1)])
            Polygon: (0, 0), (1, 0), (0, 1)
        """
        from sage.modules.free_module_element import vector
        from sage.modules.free_module import VectorSpace
        from sage.structure.sequence import Sequence

        base_ring = None
        if 'ring' in kwds:
            base_ring = kwds.pop('ring')
        if 'base_ring' in kwds:
            base_ring = kwds.pop('base_ring')
        if 'field' in kwds:
            base_ring = kwds.pop('field')

        vertices = edges = None
        if 'edges' in kwds:
            edges = kwds.pop('edges')
        elif 'vertices' in kwds:
            vertices = kwds.pop('vertices')
        elif args:
            edges = args
        else:
            raise ValueError

        if vertices is not None:
            vertices = list(map(vector, vertices))
            if base_ring is None:
                base_ring = Sequence(vertices).universe().base_ring()
        if edges is not None:
            edges = list(map(vector, edges))
            if base_ring is None:
                base_ring = Sequence(edges).universe().base_ring()

        if base_ring not in Fields():
            base_ring = base_ring.fraction_field()

        return ConvexPolygons(base_ring)(vertices=vertices, edges=edges)

polygons = PolygonsConstructor()

def regular_octagon(field=None):
    from sage.misc.superseded import deprecation
    deprecation(33, "Do not use this function anymore but regular_ngon(8)")
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


