r"""
Polygons embedded in the plane R^2.

This file implements polygons with

 - action of matrices in GL^+(2,R)
 - conversion between ground fields

The emphasis is mostly on convex polygons but there is some limited support
for non-convex polygons.

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

from sage.all import cached_method, Parent, UniqueRepresentation, Sets,\
                     Fields, ZZ, QQ, AA, RR, QQbar, matrix, polygen, vector,\
                     free_module_element
from sage.structure.element import get_coercion_model, Vector
from sage.structure.coerce import py_scalar_parent
cm = get_coercion_model()
from sage.structure.element import Element
from sage.categories.action import Action
from sage.arith.functions import lcm
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.structure.sequence import Sequence


from .matrix_2x2 import angle
from .subfield import number_field_elements_from_algebraics

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

def segment_intersect(e1, e2, base_ring=None):
    r"""
    Return whether the segments ``e1`` and ``e2`` intersect.

    OUTPUT:

    - ``0`` - do not intersect
    - ``1`` - one endpoint in common
    - ``2`` - non-trivial intersection

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import segment_intersect
        sage: for (u,v,ans) in [(((0,0),(1,0)),((0,1),(0,3)),0),
        ....:         (((0,0),(1,0)),((0,0),(0,3)),1), (((0,0),(1,0)),((0,-1),(0,3)),2),
        ....:         (((-1,-1),(1,1)),((0,0),(2,2)),2), (((-1,-1),(1,1)),((1,1),(2,2)),1)]:
        ....:     assert segment_intersect(u,v) == ans

        sage: for _ in range(4000):
        ....:     us = (randint(-4, 4), randint(-4, 4))
        ....:     ut = (randint(-4, 4), randint(-4, 4))
        ....:     vs = (randint(-4, 4), randint(-4, 4))
        ....:     vt = (randint(-4, 4), randint(-4, 4))
        ....:     if us == ut or vs == vt:
        ....:         continue
        ....:     ans1 = segment_intersect((us,ut),(vs,vt))
        ....:     ans2 = segment_intersect((ut,us),(vs,vt))
        ....:     ans3 = segment_intersect((us,ut),(vt,vs))
        ....:     ans4 = segment_intersect((ut,us),(vt,vs))
        ....:     ans5 = segment_intersect((vs,vt),(us,ut))
        ....:     ans6 = segment_intersect((vt,vs),(us,ut))
        ....:     ans7 = segment_intersect((vs,vt),(ut,us))
        ....:     ans8 = segment_intersect((vt,vs),(ut,us))
        ....:     assert (ans1 == ans2 == ans3 == ans4 == ans5 == ans6 == ans7 == ans8), (us, ut, vs, vt, ans1, ans2, ans3, ans4, ans5, ans6, ans7, ans8)
    """
    if e1[0] == e1[1] or e2[0] == e2[1]:
        raise ValueError("degenerate segments")

    if base_ring is None:
        elts = [e[i][j] for e in (e1,e2) for i in (0,1) for j in (0,1)]
        base_ring = cm.common_parent(*elts)
        if isinstance(base_ring, type):
            base_ring = py_scalar_parent(base_ring)
    m = matrix(base_ring, 3)
    xs1, ys1 = map(base_ring, e1[0])
    xt1, yt1 = map(base_ring, e1[1])
    xs2, ys2 = map(base_ring, e2[0])
    xt2, yt2 = map(base_ring, e2[1])

    m[0] = [xs1, ys1, 1]
    m[1] = [xt1, yt1, 1]
    m[2] = [xs2, ys2, 1]
    s0 = m.det()
    m[2] = [xt2, yt2, 1]
    s1 = m.det()
    if (s0 > 0 and s1 > 0) or (s0 < 0 and s1 < 0):
        # e2 stands on one side of the line generated by e1
        return 0

    m[0] = [xs2, ys2, 1]
    m[1] = [xt2, yt2, 1]
    m[2] = [xs1, ys1, 1]
    s2 = m.det()
    m[2] = [xt1, yt1, 1]
    s3 = m.det()
    if (s2 > 0 and s3 > 0) or (s2 < 0 and s3 < 0):
        # e1 stands on one side of the line generated by e2
        return 0

    if s0 == 0 and s1 == 0:
        assert s2 == 0 and s3 == 0
        if xt1 < xs1 or (xt1 == xs1 and yt1 < ys1):
            xs1,xt1 = xt1,xs1
            ys1,yt1 = yt1,ys1
        if xt2 < xs2 or (xt2 == xs2 and yt2 < ys2):
            xs2,xt2 = xt2,xs2
            ys2,yt2 = yt2,ys2

        if xs1 == xt1 == xs2 == xt2:
            xs1, xt1, xs2, xt2 = ys1, yt1, ys2, yt2

        assert xs1 < xt1 and xs2 < xt2, (xs1, xt1, xs2, xt2)

        if (xs2 > xt1) or (xt2 < xs1):
            return 0 # no intersection
        elif (xs2 == xt1) or (xt2 == xs1):
            return 1 # one endpoint in common
        else:
            assert xs1 <= xs2 < xt1 or xs1 < xt2 <= xt1 or \
                   (xs2 < xs1 and xt2 > xt1) or \
                   (xs2 > xs1 and xt2 < xt1), (xs1, xt1, xs2, xt2)
            return 2 # one dimensional

    elif s0 == 0 or s1 == 0:
        # treat alignment here
        if s2 == 0 or s3 == 0:
            return 1 # one endpoint in common
        else:
            return 2 # intersection in the middle

    return 2 # middle intersection

def is_between(e0, e1, f):
    r"""
    Check whether the vector ``f`` is strictly in the sector formed by the vectors
    ``e0`` and ``e1`` (in counter-clockwise order).

    TESTS::

        sage: from flatsurf.geometry.polygon import is_between
        sage: V = ZZ^2
        sage: vecs = [V((1,0)), V((2,1)), V((1,1)), V((1,2)),
        ....:  V((0,1)), V((-1,2)), V((-1,1)), V((-2,1)),
        ....:  V((-1,0)), V((-2,-1)), V((-1,-1)), V((-1,-2)),
        ....:  V((0,-1)), V((1,-2)), V((1,-1)), V((2,-1))]
        sage: for i in range(len(vecs)):
        ....:     for j in range(len(vecs)):
        ....:         if i == j: continue
        ....:         for k in range(len(vecs)):
        ....:             if k == i or k == j: continue
        ....:             ans = is_between(vecs[i], vecs[j], vecs[k])
        ....:             expected = i < k < j or k < j < i or j < i < k
        ....:             if ans != expected:
        ....:                 print((i,j,k),expected,ans)
    """
    if e0[0] * e1[1] > e1[0] * e0[1]:
        # positive determinant
        # [ e0[0] e1[0] ]^-1 = [ e1[1] -e1[0] ]
        # [ e0[1] e1[1] ]      [-e0[1]  e0[0] ]
        # f[0] * e1[1] - e1[0] * f[1] > 0
        # - f[0] * e0[1] + e0[0] * f[1] > 0
        return e1[1] * f[0] > e1[0] * f[1] and e0[0] * f[1] > e0[1] * f[0]
    elif e0[0] * e1[1] == e1[0] * e0[1]:
        # aligned vector
        return e0[0] * f[1] > e0[1] * f[0]
    else:
        # negative determinant
        # [ e1[0] e0[0] ]^-1 = [ e0[1] -e0[0] ]
        # [ e1[1] e0[1] ]      [-e1[1]  e1[0] ]
        # f[0] * e0[1] - e0[0] * f[1] > 0
        # - f[0] * e1[1] + e1[0] * f[1] > 0
        return e0[1] * f[0] <= e0[0] * f[1] or e1[0] * f[1] <= e1[1] * f[0]

def triangulate(vertices):
    r"""
    Return a triangulation of the list of vectors ``vertices``.

    This function assumes that ``vertices`` form the vertices of a polygon
    enumerated in counter-clockwise order.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import triangulate
        sage: V = ZZ**2
        sage: verts = list(map(V, [(0,0), (1,0), (1,1), (0,1)]))
        sage: triangulate(verts)
        [(0, 2)]

        sage: quad = [(0,0), (1,-1), (0,1), (-1,-1)]
        sage: quad = list(map(V, quad))
        sage: for i in range(4):
        ....:     print(triangulate(quad[i:] + quad[:i]))
        [(0, 2)]
        [(1, 3)]
        [(0, 2)]
        [(1, 3)]

        sage: poly = [(0,0),(1,1),(2,0),(3,1),(4,0),(4,2),
        ....:     (-4,2),(-4,0),(-3,1),(-2,0),(-1,1)]
        sage: poly = list(map(V, poly))
        sage: triangulate(poly)
        [(1, 3), (3, 5), (5, 8), (6, 8), (8, 10), (10, 1), (1, 5), (5, 10)]
        sage: for i in range(len(poly)):
        ....:     _ = triangulate(poly[i:] + poly[:i])

        sage: poly = [(0,0), (1,0), (2,0), (2,1), (2,2), (1,2), (0,2), (0,1)]
        sage: poly = list(map(V, poly))
        sage: edges = triangulate(poly)
        sage: edges
        [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
        sage: for i in range(len(poly)):
        ....:     _ = triangulate(poly[i:] + poly[:i])

        sage: poly = [(0,0), (1,2), (3,3), (1,4), (0,6), (-1,4), (-3,-3), (-1,2)]
        sage: poly = list(map(V, poly))
        sage: triangulate(poly)
        [(0, 3), (1, 3), (3, 5), (5, 7), (7, 3)]
        sage: for i in range(len(poly)):
        ....:     _ = triangulate(poly[i:] + poly[:i])


        sage: x = polygen(QQ)
        sage: p = x^4 - 5*x^2 + 5
        sage: r = AA.polynomial_root(p, RIF(1.17,1.18))
        sage: K.<a> = NumberField(p, embedding=r)
        sage: V = K**2
        sage: poly = [(1/2*a^2 - 3/2, 1/2*a),
        ....:  (-a^3 + 2*a^2 + 2*a - 4, 0),
        ....:  (1/2*a^2 - 3/2, -1/2*a),
        ....:   (1/2*a^3 - a^2 - 1/2*a + 1, 1/2*a^2 - a),
        ....:   (-1/2*a^2 + 1, 1/2*a^3 - 3/2*a),
        ....:   (-1/2*a + 1, a^3 - 3/2*a^2 - 2*a + 5/2),
        ....:   (1, 0),
        ....:   (-1/2*a + 1, -a^3 + 3/2*a^2 + 2*a - 5/2),
        ....:   (-1/2*a^2 + 1, -1/2*a^3 + 3/2*a),
        ....:   (1/2*a^3 - a^2 - 1/2*a + 1, -1/2*a^2 + a)]
        sage: poly = list(map(V, poly))
        sage: triangulate(poly)
        [(0, 3), (1, 3), (3, 5), (5, 7), (7, 9), (9, 3), (3, 7)]
    """
    n = len(vertices)
    if n < 3:
        raise ValueError
    if n == 3:
        return []
    for i in range(n - 1):
        eiright = vertices[(i+1)%n] - vertices[i]
        eileft = vertices[(i-1)%n] - vertices[i]
        for j in range(i + 2, (n if i else n-1)):
            ejright = vertices[(j+1)%n] - vertices[j]
            ejleft = vertices[(j-1)%n] - vertices[j]
            chord = vertices[j] - vertices[i]

            # check angles with neighbouring edges
            if not (is_between(eiright, eileft, chord) and \
                    is_between(ejright, ejleft, -chord)):
                continue

            # check intersection with other edges
            e = (vertices[i], vertices[j])
            good = True
            for k in range(n):
                f = (vertices[k], vertices[(k+1)%n])
                if k == (i - 1) % n or k == i or \
                   k == (j - 1) % n or k == j:
                       assert segment_intersect(e, f) == 1
                       continue
                res = segment_intersect(e, f)
                if res:
                    assert res == 2
                    good = False
                    break
            if good:
                part0 = [(s+i, t+i) for s,t in triangulate(vertices[i:j+1])]
                part1 = []
                for (s,t) in triangulate(vertices[j:] + vertices[:i+1]):
                    if s < n-j:
                        s += j
                    else:
                        s -= n - j
                    if t < n-j:
                        t += j
                    else:
                        t -= n - j
                    part1.append((s, t))
                return [(i, j)] + part0 + part1
    raise RuntimeError("input {} must be wrong".format(vertices))

def build_faces(n, edges):
    r"""
    Given a combinatorial list of pairs ``edges`` forming a cell-decomposition
    of a polygon (with vertices labeled from ``0`` to ``n-1``) return the list
    of cells.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import build_faces
        sage: build_faces(4, [(0,2)])
        [[0, 1, 2], [2, 3, 0]]
        sage: build_faces(4, [(1,3)])
        [[1, 2, 3], [3, 0, 1]]
        sage: build_faces(5, [(0,2), (0,3)])
        [[0, 1, 2], [3, 4, 0], [0, 2, 3]]
        sage: build_faces(5, [(0,2)])
        [[0, 1, 2], [2, 3, 4, 0]]
        sage: build_faces(5, [(1,4)])
        [[1, 2, 3, 4], [4, 0, 1]]
        sage: build_faces(5, [(1,3),(3,0)])
        [[1, 2, 3], [3, 4, 0], [0, 1, 3]]
    """
    polygons = [list(range(n))]
    for u,v in edges:
        j = None
        for i,p in enumerate(polygons):
            if u in p and v in p:
                if j is not None:
                    raise RuntimeError
                j = i
        if j is None:
            raise RuntimeError
        p = polygons[j]
        i0 = p.index(u)
        i1 = p.index(v)
        if i0 > i1:
            i0, i1 = i1, i0
        polygons[j] = p[i0:i1+1]
        polygons.append(p[i1:] + p[:i0+1])
    return polygons

class MatrixActionOnPolygons(Action):
    def __init__(self, polygons):
        from sage.matrix.matrix_space import MatrixSpace
        K = polygons.field()
        Action.__init__(self, MatrixSpace(K,2), polygons, True, operator.mul)

    def _act_(self, g, x):
        r"""
        Apply the 2x2 matrix `g` to the polygon `x`.

        The matrix must have non-zero determinant. If the determinant is
        negative, then the vertices and edges are relabeled according to the
        involutions `v \mapsto (n-v)%n` and  `e \mapsto n-1-e` respectively.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: p = polygons(vertices = [(1,0),(0,1),(-1,-1)])
            sage: print(p)
            Polygon: (1, 0), (0, 1), (-1, -1)
            sage: r = matrix(ZZ,[[0,1], [1,0]])
            sage: print(r*p)
            Polygon: (0, 1), (-1, -1), (1, 0)
        """
        det = g.det()
        if det > 0:
            return x.parent()(vertices=[g*v for v in x.vertices()])
        if det < 0:
            # Note that in this case we reverse the order
            vertices = [g*x.vertex(0)]
            for i in range(x.num_edges() - 1, 0, -1):
                vertices.append(g * x.vertex(i))
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


class Polygon(Element):
    def __init__(self, parent, vertices, check=True):
        Element.__init__(self, parent)
        V = parent.vector_space()
        self._v = tuple(map(V, vertices))
        for vv in self._v: vv.set_immutable()
        if check:
            self._inside_outside_check()
            self._non_intersection_check()

    def _inside_outside_check(self):
        r"""
        TESTS::

            sage: from flatsurf import Polygons
            sage: P = Polygons(QQ)
            sage: P(vertices=[(0,0),(-1,-1),(0,1),(1,-1)])
            Traceback (most recent call last):
            ...
            ValueError: the vertices are in clockwise order
        """
        # NOTE: should we do something more efficient?
        if self.area() < 0:
            raise ValueError("the vertices are in clockwise order")

    def _non_intersection_check(self):
        r"""
        TESTS::

            sage: from flatsurf import Polygons
            sage: P = Polygons(QQ)
            sage: P(vertices=[(0,0),(2,0),(1,1),(1,-1)])
            Traceback (most recent call last):
            ...
            ValueError: edge 0 (= ((0, 0), (2, 0))) and edge 2 (= ((1, 1), (1, -1))) intersects
        """
        n = len(self._v)
        for i in range(n-1):
            ei = (self._v[i], self._v[i+1])
            for j in range(i + 1, n):
                ej = (self._v[j], self._v[(j+1)%n])
                res = segment_intersect(ei, ej)
                if j == i+1 or (i == 0 and j == n-1):
                    if res > 1:
                        raise ValueError("edge %d (= %s) and edge %d (= %s) backtrack" % (i, ei, j, ej))
                elif res > 0:
                    raise ValueError("edge %d (= %s) and edge %d (= %s) intersects" % (i, ei, j, ej))

    def __hash__(self):
        # Apparently tuples do not cache their hash!
        try:
            return self._hash
        except AttributeError:
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
        if not isinstance(self, Polygon) or not isinstance(other, Polygon):
            return NotImplemented
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
        if not isinstance(self, Polygon) or not isinstance(other, Polygon):
            return NotImplemented
        return self._v != other._v

    def __lt__(self, other): raise TypeError
    __le__ = __lt__
    __gt__ = __lt__
    __ge__ = __lt__

    def cmp(self, other):
        r"""
        Implement a total order on polygons
        """
        if not isinstance(other, Polygon):
            raise TypeError("__cmp__ only implemented for ConvexPolygons")
        if not self.parent().base_ring() == other.parent().base_ring():
            raise ValueError("__cmp__ only implemented for ConvexPolygons defined over the same base_ring")
        sign = self.num_edges() - other.num_edges()
        if sign > 0:
            return 1
        if sign < 0:
            return -1
        sign = self.area() - other.area()
        if sign > self.base_ring().zero():
            return 1
        if sign < self.base_ring().zero():
            return -1
        for v in range(1,self.num_edges()):
            p = self.vertex(v)
            q = other.vertex(v)
            sign = p[0]-q[0]
            if sign > self.base_ring().zero():
                return 1
            if sign < self.base_ring().zero():
                return -1
            sign = p[1]-q[1]
            if sign > self.base_ring().zero():
                return 1
            if sign < self.base_ring().zero():
                return -1
        return 0

    def triangulation(self):
        r"""
        Return a list of pairs of indices of vertices that together with the boundary
        form a triangulation.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1), (0,2), (-1,2), (-1,1), (-2,1),
            ....:                    (-2,0), (-1,0), (-1,-1), (0,-1)], convex=False)
            sage: P.triangulation()
            [(0, 2), (2, 8), (3, 5), (6, 8), (8, 3), (3, 6), (9, 11), (0, 9), (2, 9)]
        """
        if len(self._v) == 3:
            return []
        return triangulate(self._v)

    def translate(self, u):
        r"""
        TESTS::

            sage: from flatsurf import polygons
            sage: polygons(vertices=[(0,0), (2,0), (1,1)]).translate((3,-2))
            Polygon: (3, -2), (5, -2), (4, -1)
        """
        P = self.parent()
        u = P.vector_space()(u)
        return P.element_class(P, [u+v for v in self._v], check=False)

    def change_ring(self, R):
        r"""
        Return an equal polygon over the base ring ``R``.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons.square()
            sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2)**(1/2))
            sage: S.change_ring(K)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: S.change_ring(K).base_ring()
            Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.4142...
        """
        if R is self.base_ring():
            return self
        return ConvexPolygons(R)(vertices=self._v, check=False)

    def is_convex(self):
        for i in range(self.num_edges()):
            if wedge_product(self.edge(i), self.edge(i+1)) < 0:
                return False
        return True

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
        V = VectorSpace(RR,2)
        P = self.vertices(translation)
        return point2d(P, color='red') + line2d(P + (P[0],), color='orange') + polygon2d(P, alpha=0.3)

    def angle(self, e, numerical=False, assume_rational=False):
        r"""
        Return the angle at the begining of the start point of the edge ``e``.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: polygons.square().angle(0)
            1/4
            sage: polygons.regular_ngon(8).angle(0)
            3/8

            sage: T = polygons(vertices=[(0,0), (3,1), (1,5)])
            sage: [T.angle(i, numerical=True) for i in range(3)]
            [0.16737532973071603, 0.22741638234956674, 0.10520828791971722]
            sage: sum(T.angle(i, numerical=True) for i in range(3))   # abs tol 1e-13
            0.5
        """
        return angle(self.edge(e), - self.edge((e-1)%self.num_edges()), numerical=numerical, assume_rational=assume_rational)

    def angles(self, numerical=False, assume_rational=False):
        r"""
        Return the list of angles of this polygon (divided by `2 \pi`).

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: T = polygons(angles=[1,2,3])
            sage: [T.angle(i) for i in range(3)]
            [1/12, 1/6, 1/4]
            sage: sum(T.angle(i) for i in range(3))
            1/2
        """
        return [self.angle(i) for i in range(self.num_edges())]

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

    def is_isometric(self, other, certificate=False):
        r"""
        Return whether ``self`` and ``other`` are isometric convex polygons via an orientation
        preserving isometry.

        If ``certificate`` is set to ``True`` return also a pair ``(index, rotation)``
        of an integer ``index`` and a matrix ``rotation`` such that the given rotation
        matrix identifies this polygon with the other and the edges 0 in this polygon
        is mapped to the edge ``index`` in the other.

        .. TODO::

            Implement ``is_linearly_equivalent`` and ``is_similar``.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons.square()
            sage: S.is_isometric(S)
            True
            sage: U = matrix(2,[0,-1,1,0]) * S
            sage: U.is_isometric(S)
            True

            sage: x = polygen(QQ)
            sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2)**(1/2))
            sage: S = S.change_ring(K)
            sage: U = matrix(2, [sqrt2/2, -sqrt2/2, sqrt2/2, sqrt2/2]) * S
            sage: U.is_isometric(S)
            True

            sage: U2 = polygons((1,0), (sqrt2/2, sqrt2/2), (-1,0), (-sqrt2/2, -sqrt2/2))
            sage: U2.is_isometric(U)
            False
            sage: U2.is_isometric(U, certificate=True)
            (False, None)

            sage: S = polygons((1,0), (sqrt2/2, 3), (-2,3), (-sqrt2/2+1, -6))
            sage: T = polygons((sqrt2/2,3), (-2,3), (-sqrt2/2+1, -6), (1,0))
            sage: ans, certif = S.is_isometric(T, certificate=True)
            sage: assert ans
            sage: shift, rot = certif
            sage: polygons(edges=[rot * S.edge((k + shift) % 4) for k in range(4)], base_point=T.vertex(0)) == T
            True


            sage: T = (matrix(2, [sqrt2/2, -sqrt2/2, sqrt2/2, sqrt2/2]) * S).translate((3,2))
            sage: ans, certif = S.is_isometric(T, certificate=True)
            sage: assert ans
            sage: shift, rot = certif
            sage: polygons(edges=[rot * S.edge(k + shift) for k in range(4)], base_point=T.vertex(0)) == T
            True
        """
        if type(self) is not type(other):
            raise TypeError

        n = self.num_edges()
        if other.num_edges() != n:
            return False
        sedges = self.edges()
        oedges = other.edges()

        slengths = [x**2 + y**2 for x,y in sedges]
        olengths = [x**2 + y**2 for x,y in oedges]
        for i in range(n):
            if slengths == olengths:
                # we have a match of lengths after a shift by i
                xs,ys = sedges[0]
                xo,yo = oedges[0]
                ms = matrix(2, [xs, -ys, ys, xs])
                mo = matrix(2, [xo, -yo, yo, xo])
                rot = mo * ~ms
                assert rot.det() == 1 and (rot * rot.transpose()).is_one()
                assert oedges[0] == rot * sedges[0]
                if all(oedges[i] == rot * sedges[i] for i in range(1,n)):
                    return (True, (0 if i == 0 else n-i, rot)) if certificate else True
            olengths.append(olengths.pop(0))
            oedges.append(oedges.pop(0))
        return (False, None) if certificate else False

    def is_translate(self, other, certificate=False):
        r"""
        Return whether ``other`` is a translate of ``self``.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons(vertices=[(0,0), (3,0), (1,1)])
            sage: T1 = S.translate((2,3))
            sage: S.is_translate(T1)
            True
            sage: T2 = polygons(vertices=[(-1,1), (1,0), (2,1)])
            sage: S.is_translate(T2)
            False
            sage: T3 = polygons(vertices=[(0,0), (3,0), (2,1)])
            sage: S.is_translate(T3)
            False

            sage: S.is_translate(T1, certificate=True)
            (True, (0, 1))
            sage: S.is_translate(T2, certificate=True)
            (False, None)
            sage: S.is_translate(T3, certificate=True)
            (False, None)
        """
        if type(self) is not type(other):
            raise TypeError

        n = self.num_edges()
        if other.num_edges() != n:
            return False
        sedges = self.edges()
        oedges = other.edges()
        for i in range(n):
            if sedges == oedges:
                return (True, (i, 1)) if certificate else True
            oedges.append(oedges.pop(0))
        return (False, None) if certificate else False

    def is_half_translate(self, other, certificate=False):
        r"""
        Return whether ``other`` is a translate or half-translate of ``self``.

        If ``certificate`` is set to ``True`` then return also a pair ``(orientation, index)``.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons(vertices=[(0,0), (3,0), (1,1)])
            sage: T1 = S.translate((2,3))
            sage: S.is_half_translate(T1)
            True
            sage: T2 = polygons(vertices=[(-1,1), (1,0), (2,1)])
            sage: S.is_half_translate(T2)
            True
            sage: T3 = polygons(vertices=[(0,0), (3,0), (2,1)])
            sage: S.is_half_translate(T3)
            False

            sage: S.is_half_translate(T1, certificate=True)
            (True, (0, 1))
            sage: ans, certif = S.is_half_translate(T2, certificate=True)
            sage: assert ans
            sage: shift, rot = certif
            sage: polygons(edges=[rot * S.edge(k + shift) for k in range(3)], base_point=T2.vertex(0)) == T2
            True
            sage: S.is_half_translate(T3, certificate=True)
            (False, None)
        """
        if type(self) is not type(other):
            raise TypeError

        n = self.num_edges()
        if other.num_edges() != n:
            return False

        sedges = self.edges()
        oedges = other.edges()
        for i in range(n):
            if sedges == oedges:
                return (True, (i, 1)) if certificate else True
            oedges.append(oedges.pop(0))

        assert oedges == other.edges()
        oedges = [-e for e in oedges]
        for i in range(n):
            if sedges == oedges:
                return (True, (0 if i == 0 else n-i, -1)) if certificate else True
            oedges.append(oedges.pop(0))

        return (False, None) if certificate else False

class ConvexPolygon(Polygon):
    r"""
    A convex polygon in the plane RR^2
    """
    def __init__(self, parent, vertices, check=True):
        r"""
        To construct the polygon you should either use a list of edge vectors
        or a list of vertices. Using both will result in a ValueError. The polygon
        needs to be convex with postively oriented boundary.

        INPUT:

        - ``parent`` -- a parent

        - ``vertices`` -- a list of vertices of the polygon
        """
        Polygon.__init__(self, parent, vertices, check=False)
        if check:
            self._convexity_check()

    def is_convex(self):
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

class Polygons(UniqueRepresentation, Parent):
    Element = Polygon

    def __init__(self, field):
        Parent.__init__(self, category=Sets())
        if not field in Fields():
            raise ValueError("'field' must be a field")
        self._field = field
        self.register_action(MatrixActionOnPolygons(self))

    def base_ring(self):
        return self._field

    field = base_ring

    @cached_method
    def vector_space(self):
        r"""
        Return the vector space in which this polygon embeds.

        EXAMPLES::

            sage: from flatsurf import Polygons, ConvexPolygons
            sage: P = Polygons(QQ)
            sage: P.vector_space()
            Vector space of dimension 2 over Rational Field
            sage: C = ConvexPolygons(QQ)
            sage: C.vector_space()
            Vector space of dimension 2 over Rational Field
        """
        return VectorSpace(self.base_ring(), 2)

    def _repr_(self):
        return "Polygons(%s)"%self.base_ring()

    def _element_constructor_(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf import Polygons, ConvexPolygons

            sage: C = Polygons(QQ)
            sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
            sage: p
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: C(p) is p
            True
            sage: C((1,0), (0,1), (-1, 1))
            Traceback (most recent call last):
            ...
            ValueError: the polygon does not close up

            sage: D = ConvexPolygons(QQbar)
            sage: D(p)
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: D(vertices=p.vertices())
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: D(edges=p.edges())
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
        """
        check = kwds.pop('check', True)

        if len(args) == 1 and isinstance(args[0], Polygon):
            vertices = map(self.vector_space(), args[0].vertices())
            args = ()

        else:
            vertices = kwds.pop('vertices', None)
            edges = kwds.pop('edges', None)
            base_point = kwds.pop('base_point', (0,0))

            if (vertices is None) and (edges is None):
                if len(args) == 1:
                    edges = args[0]
                elif args:
                    edges = args
                else:
                    raise ValueError("exactly one of 'vertices' or 'edges' must be provided")
            if kwds:
                raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

            if edges is not None:
                v = self.vector_space()(base_point)
                vertices = []
                for e in map(self.vector_space(), edges):
                    vertices.append(v)
                    v += e
                if v != vertices[0]:
                    raise ValueError("the polygon does not close up")

        return self.element_class(self, vertices, check)

class ConvexPolygons(Polygons):
    r"""
    The set of convex polygons with a fixed base field.

    EXAMPLES::

        sage: from flatsurf import ConvexPolygons
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

            sage: from flatsurf import ConvexPolygons
            sage: C1 = ConvexPolygons(QQ)
            sage: C2 = ConvexPolygons(AA)
            sage: C2.has_coerce_map_from(C1)
            True
            sage: C1.has_coerce_map_from(C2)
            False
        """
        return isinstance(other, ConvexPolygons) and self.field().has_coerce_map_from(other.field())

    def _an_element_(self):
        return self([(1,0),(0,1),(-1,0),(0,-1)])

    def _repr_(self):
        return "ConvexPolygons(%s)"%self.base_ring()

    def _element_constructor_(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf import ConvexPolygons

            sage: C = ConvexPolygons(QQ)
            sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
            sage: p
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: C(p) is p
            True
            sage: C((1,0), (0,1), (-1, 1))
            Traceback (most recent call last):
            ...
            ValueError: the polygon does not close up

            sage: D = ConvexPolygons(QQbar)
            sage: D(p)
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: D(vertices=p.vertices())
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)
            sage: D(edges=p.edges())
            Polygon: (0, 0), (1, 0), (2, 0), (1, 1)

        """
        check = kwds.pop('check', True)

        if len(args) == 1 and isinstance(args[0], Polygon):
            vertices = map(self.vector_space(), args[0].vertices())
            args = ()

        else:
            vertices = kwds.pop('vertices', None)
            edges = kwds.pop('edges', None)
            base_point = kwds.pop('base_point', (0,0))

            if (vertices is None) and (edges is None):
                if len(args) == 1:
                    edges = args[0]
                elif args:
                    edges = args
                else:
                    raise ValueError("exactly one of 'vertices' or 'edges' must be provided")
            if kwds:
                raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

            if edges is not None:
                v = self.vector_space()(base_point)
                vertices = []
                for e in map(self.vector_space(), edges):
                    vertices.append(v)
                    v += e
                if v != vertices[0]:
                    raise ValueError("the polygon does not close up")

        return self.element_class(self, vertices, check)

class EquiangularPolygons:
    r"""
    Polygons with fixed (rational) angles.

    EXAMPLES::

        sage: from flatsurf import EquiangularPolygons

        sage: P = EquiangularPolygons(1,2,5)
        sage: P
        EquiangularPolygons(1, 2, 5) over Algebraic Real Field
        sage: L = P.lengths_polytope()    # polytope of admissible lengths for edges
        sage: L
        A 1-dimensional polyhedron in AA^3 defined as the convex hull of 1 vertex and 1 ray
        sage: lengths = L.rays()[0].vector()
        sage: lengths
        (1, 0.4142135623730951?, 0.7653668647301795?)
        sage: p = P(*lengths)    # build one polygon with the given lengths
        sage: p
        Polygon: (0, 0), (1, 0), (0.7071067811865475?, 0.2928932188134525?)
        sage: p.angles()
        [1/16, 1/8, 5/16]
        sage: P.angles(integral=False)
        [1/16, 1/8, 5/16]

        sage: P = EquiangularPolygons(1, 2, 1, 2, 2, 1)
        sage: L = P.lengths_polytope()
        sage: L
        A 4-dimensional polyhedron in AA^6 defined as the convex hull of 1 vertex and 6 rays
        sage: rays = [r.vector() for r in L.rays()]
        sage: rays
        [(1, 0, 0, 0, 1, 0.3472963553338607?),
         (0, 1, 0, 0, 0.8793852415718168?, 0.6527036446661393?),
         (0.6527036446661393?, 0, 1, 0, 0, 0.8793852415718168?),
         (0.8793852415718168?, 0, 0, 1, 0, 0.6527036446661393?),
         (0, 0.6527036446661393?, 0.8793852415718168?, 0, 0, 1),
         (0, 0.8793852415718168?, 0, 0.8793852415718168?, 0, 0.8793852415718168?)]
        sage: lengths = 3*rays[0] + rays[2] + 2*rays[3] + rays[4]
        sage: p = P(*lengths)
        sage: p
        Polygon: (0, 0), (5.411474127809773?, 0), (6.024814926262611?, 0.2232377940978994?), (5.085122305476703?, 1.850833156796647?), (3.553033419238747?, 3.136408376169726?), (0.7339555568810219?, 4.162468806146732?)
        sage: p.angles()
        [2/9, 4/9, 2/9, 4/9, 4/9, 2/9]
    """
    def __init__(self, *angles, **kwds):
        number_field = kwds.pop('number_field', False)
        if kwds:
            raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))
        if len(angles) == 1 and isinstance(angles[0], (tuple, list)):
            angles = angles[0]

        n = len(angles)
        if n < 3:
            raise ValueError("'angles' should be a list of at least 3 numbers")
        angles = [QQ.coerce(a) for a in angles]  # total sum of angle should be
        if any(angle <= 0 for angle in angles):
            raise ValueError("'angles' must be positive rational numbers")

        # normalize the sum so that it is (n-2)/2 (ie in multiple of 2pi)
        s = ZZ(n - 2) / sum(angles) / ZZ(2)
        if s != 1:
            angles = [s * a for a in angles]
        if any(angle <= 0 or angle >= 1 for angle in angles):
            raise ValueError("each angle must be > 0 and < 2 pi")
        self._angles = angles

        zetas = [QQbar.zeta(a.denominator()) ** a.numerator() for a in angles]
        cosines = [z.real() for z in zetas]
        sines = [z.imag() for z in zetas]
        if number_field:
            base_ring, elts = number_field_elements_from_algebraics(cosines + sines, name='a')
            cosines = elts[:n]
            sines = elts[n:]
        else:
            base_ring = AA
        self._cosines = cosines
        self._sines = sines
        self._base_ring = base_ring
        assert len(cosines) == n and len(sines) == n

    def convexity(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons
            sage: EquiangularPolygons(1, 2, 5).convexity()
            True
            sage: EquiangularPolygons(2, 2, 3, 13).convexity()
            False
        """
        return all(2 * a <= 1 for a in self._angles)

    def strict_convexity(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons
            sage: E = EquiangularPolygons([1, 1, 1, 1, 2])
            sage: E.angles()
            [1/4, 1/4, 1/4, 1/4, 1/2]
            sage: E.convexity()
            True
            sage: E.strict_convexity()
            False
        """
        return all(2 * a < 1 for a in self._angles)

    def angles(self, integral=False):
        angles = self._angles
        if integral:
            D = lcm([a.denominator() for a in self._angles])
            angles = [(D * a).numerator() for a in self._angles]
        return angles

    def __repr__(self):
        r"""
        TESTS::

            sage: from flatsurf import EquiangularPolygons
            sage: EquiangularPolygons(1, 2, 3)
            EquiangularPolygons(1, 2, 3) over Algebraic Real Field
            sage: EquiangularPolygons(1, 2, 3, number_field=True)
            EquiangularPolygons(1, 2, 3) over Number Field in a with defining polynomial y^2 - 3 with a = 1.732050807568878?
        """
        return "EquiangularPolygons({}) over {}".format(", ".join(map(str,self.angles(True))), self._base_ring)

    def vector_space(self):
        return VectorSpace(self._base_ring, 2)

    def slopes(self, e0=(1,0), cosines=None, sines=None):
        r"""
        List of slopes of the edges as a list of unit vectors.

        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons
            sage: EquiangularPolygons(1, 2, 1, 2).slopes()
            [(1, 0),
             (0.500000000000000?, 0.866025403784439?),
             (-1.000000000000000?, 0.?e-18),
             (-0.500000000000000?, -0.866025403784439?)]
        """
        V = self.vector_space()
        if cosines is None:
            cosines = self._cosines
        if sines is None:
            sines = self._sines
        n = len(cosines)
        v = V.zero()
        e = V(e0)
        edges = [e]
        for i in range(n-1):
            e = V((-cosines[i+1] * e[0] - sines[i+1] * e[1], sines[i+1] * e[0] - cosines[i+1] * e[1]))
            edges.append(e)
        return edges

    def lengths_polytope(self):
        r"""
        Return the polytope parametrizing the admissible length data.

        Be careful that even though the lengths are admissible, they may not
        define a polygon without intersection.

        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons
            sage: EquiangularPolygons(1, 2, 1, 2).lengths_polytope()
            A 2-dimensional polyhedron in AA^4 defined as the convex hull of 1 vertex and 2 rays
        """
        n = len(self._angles)
        slopes = self.slopes()
        eqns = [[0] + [s[0] for s in slopes], [0] + [s[1] for s in slopes]]
        ieqs = []
        for i in range(n):
            ieq = [0] * (n + 1)
            ieq[i+1] = 1
            ieqs.append(ieq)

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(eqns=eqns, ieqs=ieqs, base_ring=self._base_ring)

    def an_element(self, *args, **kwds):
        r"""
        Return a polygon in this family.

        Note that this might fail due to intersection.

        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons
            sage: EquiangularPolygons(4, 3, 4, 4, 3, 4).an_element()
            Polygon: (0, 0), (2, 0), (2.284629676546571?, 1.979642883761866?), (0.7720837808745704?, 3.725213899997057?), (-1.227916219125430?, 3.725213899997057?), (-1.512545895672000?, 1.745571016235191?)
        """
        return self(sum(r.vector() for r in self.lengths_polytope().rays()), *args, **kwds)

    def __call__(self, *lengths, **kwds):
        r"""
        TESTS::

            sage: from flatsurf import EquiangularPolygons
            sage: P = EquiangularPolygons(1, 2, 1, 2)
            sage: L = P.lengths_polytope()
            sage: r0, r1 = [r.vector() for r in L.rays()]
            sage: lengths = r0 + r1
            sage: P(*lengths)
            Polygon: (0, 0), (1, 0), (3/2, 0.866025403784439?), (1/2, 0.866025403784439?)
            sage: P(*lengths[:-2])
            Polygon: (0, 0), (1, 0), (3/2, 0.866025403784439?), (1/2, 0.866025403784439?)

            sage: P = EquiangularPolygons(2, 2, 3, 13)
            sage: P([1, 1])
            Polygon: (0, 0), (1, 0), (0.1909830056250526?, 0.5877852522924731?), (0.1909830056250526?, 0.1387572757128878?)
        """
        number_field = kwds.pop('number_field', False)
        if kwds:
            raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))
        if len(lengths) == 1 and isinstance(lengths[0], (tuple, list, Vector)):
            lengths = lengths[0]

        n = len(self._angles)
        if len(lengths) != n-2 and len(lengths) != n:
            raise ValueError("invalid 'lengths' argument")

        cosines = self._cosines
        sines = self._sines
        base_ring = self._base_ring
        if number_field:
            if base_ring is AA:
                base_ring, elts = number_field_elements_from_algebraics(cosines + sines + list(lengths), name='a')
                cosines = elts[:n]
                sines = elts[n:2*n]
                lengths = elts[2*n:]
            elif any(l not in self._base_ring for l in lengths):
                raise NotImplementedError

        V = self.vector_space()
        slopes = self.slopes(cosines=cosines, sines=sines)
        v = V((0,0))
        vertices = [v]

        if len(lengths) == n - 2:
            for i in range(n-2):
                v += lengths[i] * slopes[i]
                vertices.append(v)
            s,t = matrix([slopes[-1],slopes[n-2]]).solve_left(vertices[0] - vertices[n-2])
            assert vertices[0] - s*slopes[-1] == vertices[n-2] + t*slopes[n-2]
            if s <= 0 or t <= 0:
                raise ValueError("the provided lengths do not give rise to a polygon")
            vertices.append(vertices[0] - s*slopes[-1])

        elif len(lengths) == n:
            for i in range(n):
                v += lengths[i] * slopes[i]
                vertices.append(v)
            if not vertices[-1].is_zero():
                raise ValueError("the provided lengths do not give rise to a polygon")
            vertices.pop(-1)

        if self.convexity():
            return ConvexPolygons(base_ring)(vertices=vertices)
        else:
            return Polygons(base_ring)(vertices=vertices)

    def billiard_unfolding_angles(self, cover_type="translation"):
        r"""
        Return the angles of the unfolding rational, half-translation or translation surface.

        INPUT:

        - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
          ``"half-translation"`` or ``"translation"``

        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons

            sage: E = EquiangularPolygons(1, 2, 5)
            sage: E.billiard_unfolding_angles(cover_type="rational")
            {1/8: 1, 1/4: 1, 5/8: 1}
            sage: (1/8 - 1) + (1/4 - 1) + (5/8 - 1)  # Euler characteristic (of the sphere)
            -2
            sage: E.billiard_unfolding_angles(cover_type="half-translation")
            {1/2: 3, 5/2: 1}
            sage: E.billiard_unfolding_angles(cover_type="translation")
            {1: 3, 5: 1}

            sage: E = EquiangularPolygons(1, 3, 1, 7)
            sage: E.billiard_unfolding_angles(cover_type="rational")
            {1/6: 2, 1/2: 1, 7/6: 1}
            sage: 2 * (1/6 - 1) + (1/2 - 1) + (7/6 - 1) # Euler characteristic
            -2
            sage: E.billiard_unfolding_angles(cover_type="half-translation")
            {1/2: 5, 7/2: 1}
            sage: E.billiard_unfolding_angles(cover_type="translation")
            {1: 5, 7: 1}

            sage: E = EquiangularPolygons(1, 3, 5, 7)
            sage: E.billiard_unfolding_angles(cover_type="rational")
            {1/8: 1, 3/8: 1, 5/8: 1, 7/8: 1}
            sage: (1/8 - 1) + (3/8 - 1) + (5/8 - 1) + (7/8 - 1) # Euler characteristic
            -2
            sage: E.billiard_unfolding_angles(cover_type="half-translation")
            {1/2: 1, 3/2: 1, 5/2: 1, 7/2: 1}
            sage: E.billiard_unfolding_angles(cover_type="translation")
            {1: 1, 3: 1, 5: 1, 7: 1}

            sage: E = EquiangularPolygons(1, 2, 8)
            sage: E.billiard_unfolding_angles(cover_type="rational")
            {1/11: 1, 2/11: 1, 8/11: 1}
            sage: (1/11 - 1) + (2/11 - 1) + (8/11 - 1) # Euler characteristic
            -2
            sage: E.billiard_unfolding_angles(cover_type="half-translation")
            {1: 1, 2: 1, 8: 1}
            sage: E.billiard_unfolding_angles(cover_type="translation")
            {1: 1, 2: 1, 8: 1}
        """
        rat_angles = {}
        for a in self.angles():
            if 2*a in rat_angles:
                rat_angles[2*a] += 1
            else:
                rat_angles[2*a] = 1
        if cover_type == "rational":
            return rat_angles

        N = lcm([x.denominator() for x in rat_angles])
        if N % 2:
            N *= 2

        cov_angles = {}
        for x, mult in rat_angles.items():
            y = x.numerator()
            d = x.denominator()
            if d%2:
                d *= 2
            else:
                y = y/2
            assert N % d == 0
            if y in cov_angles:
                cov_angles[y] += mult * N//d
            else:
                cov_angles[y] = mult * N//d

        if cover_type == "translation" and any(y.denominator() == 2 for y in cov_angles):
            covcov_angles = {}
            for y,mult in cov_angles.items():
                yy = y.numerator()
                if yy not in covcov_angles:
                    covcov_angles[yy] = 0
                covcov_angles[yy] += 2//y.denominator() * mult
            return covcov_angles
        elif cover_type == "half-translation" or cover_type == "translation":
            return cov_angles
        else:
            raise ValueError("unknown 'cover_type' {!r}".format(cover_type))

    def billiard_unfolding_stratum(self, cover_type="translation", marked_points=False):
        r"""
        Return the stratum of quadratic or Abelian differential obtained by
        unfolding a billiard in a polygon of this equiangular family.

        INPUT:

        - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
          ``"half-translation"`` or ``"translation"``

        - ``marked_poins`` (optional, default ``False``) - whether the stratum should
          have regular marked points

        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons, similarity_surfaces

            sage: E = EquiangularPolygons(1, 2, 5)
            sage: E.billiard_unfolding_stratum("half-translation")
            Q_1(3, -1^3)
            sage: E.billiard_unfolding_stratum("translation")
            H_3(4)
            sage: E.billiard_unfolding_stratum("half-translation", True)
            Q_1(3, -1^3)
            sage: E.billiard_unfolding_stratum("translation", True)
            H_3(4, 0^3)

            sage: E = EquiangularPolygons(1, 3, 1, 7)
            sage: E.billiard_unfolding_stratum("half-translation")
            Q_1(5, -1^5)
            sage: E.billiard_unfolding_stratum("translation")
            H_4(6)
            sage: E.billiard_unfolding_stratum("half-translation", True)
            Q_1(5, -1^5)
            sage: E.billiard_unfolding_stratum("translation", True)
            H_4(6, 0^5)

            sage: P = E.an_element()
            sage: S = similarity_surfaces.billiard(P)
            sage: S.minimal_cover("half-translation").stratum()
            Q_1(5, -1^5)
            sage: S.minimal_cover("translation").stratum()
            H_4(6, 0^5)

            sage: E = EquiangularPolygons(1, 3, 5, 7)
            sage: E.billiard_unfolding_stratum("half-translation")
            Q_3(5, 3, 1, -1)
            sage: E.billiard_unfolding_stratum("translation")
            H_7(6, 4, 2)

            sage: P = E.an_element()
            sage: S = similarity_surfaces.billiard(P)
            sage: S.minimal_cover("half-translation").stratum()
            Q_3(5, 3, 1, -1)
            sage: S.minimal_cover("translation").stratum()
            H_7(6, 4, 2, 0)

            sage: E = EquiangularPolygons(1, 2, 8)
            sage: E.billiard_unfolding_stratum("half-translation")
            H_5(7, 1)
            sage: E.billiard_unfolding_stratum("translation")
            H_5(7, 1)

            sage: E.billiard_unfolding_stratum("half-translation", True)
            H_5(7, 1, 0)
            sage: E.billiard_unfolding_stratum("translation", True)
            H_5(7, 1, 0)

            sage: E = EquiangularPolygons(9, 6, 3, 2)
            sage: p = E.an_element()
            sage: B = similarity_surfaces.billiard(p)
            sage: B.minimal_cover("half-translation").stratum()
            Q_4(7, 4, 1, 0)
            sage: E.billiard_unfolding_stratum("half-translation", True)
            Q_4(7, 4, 1, 0)
            sage: B.minimal_cover("translation").stratum()
            H_8(8, 2^3, 0^2)
            sage: E.billiard_unfolding_stratum("translation", True)
            H_8(8, 2^3, 0^2)
        """
        angles = self.billiard_unfolding_angles(cover_type)
        if all(a.is_integer() for a in angles):
            from surface_dynamics import AbelianStratum
            if not marked_points and len(angles) == 1 and 1 in angles:
                return AbelianStratum([0])
            else:
                return AbelianStratum({ZZ(a-1): mult for a,mult in angles.items() if marked_points or a != 1})
        else:
            from surface_dynamics import QuadraticStratum
            return QuadraticStratum({ZZ(2*(a-1)): mult for a,mult in angles.items() if marked_points or a != 1})

    def billiard_unfolding_stratum_dimension(self, cover_type="translation", marked_points=False):
        r"""
        Return the dimension of the stratum of quadratic or Abelian differential
        obtained by unfolding a billiard in a polygon of this equiangular family.

        INPUT:

        - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
          ``"half-translation"`` or ``"translation"``

        - ``marked_poins`` (optional, default ``False``) - whether the stratum should
          have marked regular points

        EXAMPLES::

            sage: from flatsurf import EquiangularPolygons

            sage: E = EquiangularPolygons(1, 2, 5)
            sage: E.billiard_unfolding_stratum_dimension("half-translation")
            4
            sage: E.billiard_unfolding_stratum("half-translation").dimension()
            4
            sage: E.billiard_unfolding_stratum_dimension(cover_type="translation")
            6
            sage: E.billiard_unfolding_stratum("translation").dimension()
            6
            sage: E.billiard_unfolding_stratum_dimension("translation", True)
            9
            sage: E.billiard_unfolding_stratum("translation", True).dimension()
            9

            sage: E = EquiangularPolygons(1, 3, 5)
            sage: E.billiard_unfolding_stratum_dimension("half-translation")
            6
            sage: E.billiard_unfolding_stratum("half-translation").dimension()
            6
            sage: E.billiard_unfolding_stratum_dimension("translation")
            6
            sage: E.billiard_unfolding_stratum("translation").dimension()
            6

            sage: E = EquiangularPolygons(1, 3, 1, 7)
            sage: E.billiard_unfolding_stratum_dimension("half-translation")
            6

            sage: E = EquiangularPolygons(1, 3, 5, 7)
            sage: E.billiard_unfolding_stratum_dimension("half-translation")
            8

            sage: E = EquiangularPolygons(1, 2, 8)
            sage: E.billiard_unfolding_stratum_dimension()
            11
            sage: E.billiard_unfolding_stratum().dimension()
            11
            sage: E.billiard_unfolding_stratum_dimension(marked_points=True)
            12
            sage: E.billiard_unfolding_stratum(marked_points=True).dimension()
            12
        """
        if cover_type == "rational":
            raise NotImplementedError
        if cover_type != "translation" and cover_type != "half-translation":
            raise ValueError

        angles = self.billiard_unfolding_angles(cover_type)
        if not marked_points:
            if 1 in angles:
                del angles[1]
            if not angles:
                angles[1] = 1

        abelian = all(a.is_integer() for a in angles)
        s = sum(angles.values())
        chi = sum(mult * (a - 1) for a, mult in angles.items())
        assert chi.denominator() == 1
        chi = ZZ(chi)
        assert chi%2 == 0
        g = chi//2 + 1
        return (2*g + s - 1 if abelian else 2*g + s - 2)

class PolygonsConstructor:
    def square(self, side=1, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: polygons.square()
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: polygons.square(field=QQbar).parent()
            ConvexPolygons(Algebraic Field)
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
            ConvexPolygons(Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?)
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

            sage: polygons.triangle(1,2,3).angles()
            [1/12, 1/6, 1/4]
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

        return ConvexPolygons(field)(edges=edges)

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
        return ConvexPolygons(field)(edges=[(c,field.zero()),(field.zero(),s),(-c,-s)])

    def __call__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf import *

            sage: polygons((1,0),(0,1),(-1,0),(0,-1))
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: polygons((1,0),(0,1),(-1,0),(0,-1), ring=QQbar)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: _.parent()
            ConvexPolygons(Algebraic Field)

            sage: polygons(vertices=[(0,0), (1,0), (0,1)])
            Polygon: (0, 0), (1, 0), (0, 1)

            sage: polygons(edges=[(2,0),(-1,1),(-1,-1)], base_point=(3,3))
            Polygon: (3, 3), (5, 3), (4, 4)
            sage: polygons(vertices=[(0,0),(2,0),(1,1)], base_point=(3,3))
            Traceback (most recent call last):
            ...
            ValueError: invalid keyword 'base_point'


            sage: polygons(angles=[1,1,1,2], length=1)
            Polygon: (0, 0), (1, 0), (-1/2*a^2 + 5/2, 1/2*a), (-1/2*a^2 + 2, 1/2*a^3 - 3/2*a)
            sage: polygons(angles=[1,1,1,2], length=2)
            Polygon: (0, 0), (2, 0), (-a^2 + 5, a), (-a^2 + 4, a^3 - 3*a)
            sage: polygons(angles=[1,1,1,2], length=AA(2)**(1/2))
            Polygon: (0, 0), (a^5 - 5*a^3 + 5*a, 0), (1/2*a^7 - 5/2*a^5 + 3/2*a^3 + 3*a, -1/2*a^7 + 7/2*a^5 - 15/2*a^3 + 5*a), (1/2*a^7 - 3*a^5 + 4*a^3 + 1/2*a, -1/2*a^7 + 4*a^5 - 9*a^3 + 7/2*a)

            sage: polygons(angles=[1]*5).angles()
            [3/10, 3/10, 3/10, 3/10, 3/10]
            sage: polygons(angles=[1]*8).angles()
            [3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 3/8]

            sage: P = polygons(angles=[1,1,3,3], lengths=[3,1])
            sage: P.angles()
            [1/8, 1/8, 3/8, 3/8]
            sage: e0 = P.edge(0); assert e0[0]**2 + e0[1]**2 == 3**2
            sage: e1 = P.edge(1); assert e1[0]**2 + e1[1]**2 == 1

            sage: polygons(angles=[1,1,1,2])
            Polygon: (0, 0), (1, 0), (-1/2*a^2 + 5/2, 1/2*a), (-1/2*a^2 + 2, 1/2*a^3 - 3/2*a)

            sage: polygons(angles=[1,1,1,8])
            Traceback (most recent call last):
            ...
            ValueError: 'angles' do not determine convex polygon; you might want to set the option 'convex=False'
            sage: polygons(angles=[1,1,1,8], convex=False)
            Traceback (most recent call last):
            ...
            ValueError: non-convex equiangular polygon; lengths must be provided
            sage: polygons(angles=[1,1,1,8], lengths=[1,1], convex=False)
            Polygon: (0, 0), (1, 0), (1/2*a^8 - 9/2*a^6 + 27/2*a^4 - 15*a^2 + 11/2, 1/2*a), (1/2*a^8 - 9/2*a^6 + 13*a^4 - 13*a^2 + 4, -1/2*a^9 + 9/2*a^7 - 13*a^5 + 13*a^3 - 7/2*a)

        TESTS::

            sage: from itertools import product
            sage: for a,b,c in product(range(1,5), repeat=3):
            ....:     if gcd([a,b,c]) != 1:
            ....:         continue
            ....:     T = polygons(angles=[a,b,c])
            ....:     D = 2*(a+b+c)
            ....:     assert T.angles() == [a/D, b/D, c/D]
            ....:     assert T.edge(0) == T.vector_space()((1,0))
        """
        base_ring = None
        if 'ring' in kwds:
            base_ring = kwds.pop('ring')
        elif 'base_ring' in kwds:
            base_ring = kwds.pop('base_ring')
        elif 'field' in kwds:
            base_ring = kwds.pop('field')

        convex = kwds.pop('convex', True)

        vertices = edges = angles = base_point = None
        if 'edges' in kwds:
            edges = kwds.pop('edges')
            base_point = kwds.pop('base_point', (0,0))
        elif 'vertices' in kwds:
            vertices = kwds.pop('vertices')
        elif 'angles' in kwds:
            angles = kwds.pop('angles')
            lengths = kwds.pop('lengths', None)
            length = kwds.pop('length', None)
            base_point = kwds.pop('base_point', (0,0))
            number_field = kwds.pop('number_field', True)
        elif args:
            edges = args
            args = ()
            base_point = kwds.pop('base_point', (0,0))

        if (vertices is not None) + (edges is not None) + (angles is not None) != 1:
            raise ValueError("exactly one of 'vertices', 'edges' or 'angles' should be provided")

        if vertices is None and edges is None and angles is None and lengths is None:
            raise ValueError("either vertices, edges or angles should be provided")
        if args:
            raise ValueError("invalid argument {!r}".format(args))
        if kwds:
            raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

        if vertices is not None:
            vertices = list(map(vector, vertices))
            if base_ring is None:
                base_ring = Sequence([x for x,_ in vertices] + [y for _,y in vertices]).universe()
                if isinstance(base_ring, type):
                    base_ring = py_scalar_parent(base_ring)

        elif edges is not None:
            edges = list(map(vector, edges))
            if base_ring is None:
                base_ring = Sequence([x for x,_ in edges] + [y for _,y in edges] + list(base_point)).universe()
                if isinstance(base_ring, type):
                    base_ring = py_scalar_parent(base_ring)

        elif angles is not None:
            E = EquiangularPolygons(*angles, number_field=False)
            if convex and not E.convexity():
                raise ValueError("'angles' do not determine convex polygon; you might want to set the option 'convex=False'")
            n = len(angles)
            if length is None and lengths is None:
                if not E.convexity():
                    raise ValueError("non-convex equiangular polygon; lengths must be provided")
                lengths = [AA(1)] * (n-2)
            elif lengths is None:
                if not E.convexity():
                    raise ValueError("non-convex equiangular polygon; lengths must be provided")
                lengths = [AA.coerce(length)] * (n-2)
            elif length is None:
                if len(lengths) != n-2:
                    raise ValueError("'lengths' must be a list of n-2 numbers (one less than 'angles')")
                lengths = [AA.coerce(length) for length in lengths]
            else:
                raise ValueError("only one of 'length' or 'lengths' can be set together with 'angles'")

            return E(lengths, number_field=True)

        if base_ring not in Fields():
            base_ring = base_ring.fraction_field()

        if convex:
            return ConvexPolygons(base_ring)(vertices=vertices, edges=edges, base_point=base_point)
        else:
            return Polygons(base_ring)(vertices=vertices, edges=edges, base_point=base_point)

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
        return ConvexPolygons(self._field)(self._w)


