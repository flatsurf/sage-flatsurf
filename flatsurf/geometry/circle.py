r"""
This class contains methods useful for working with circles. 

This will be used to build a LazyDelaunayTriangulation class which will compute the
Delaunay decomposition for infinite surfaces.
"""

from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector

def circle_from_three_points(p,q,r,base_ring=None):
    r"""
    Construct a circle from three points on the circle.
    """
    if base_ring is None:
        base_ring=p.base_ring()
    V2 = VectorSpace(base_ring,2)
    V3 = VectorSpace(base_ring,3)
    
    v1=V3((p[0]+q[0],p[1]+q[1],2))
    v2=V3((p[1]-q[1],q[0]-p[0],0))
    line1=v1.cross_product(v2)
    v1=V3((p[0]+r[0],p[1]+r[1],2))
    v2=V3((p[1]-r[1],r[0]-p[0],0))
    line2=v1.cross_product(v2)
    center_3 = line1.cross_product(line2)
    if center_3[2].is_zero():
        raise ValueError("The three points lie on a line.")
    center = V2( (center_3[0]/center_3[2], center_3[1]/center_3[2]) )
    return Circle(center, (p[0]-center[0])**2+(p[1]-center[1])**2 )

class Circle:
    def __init__(self, center, radius_squared, base_ring=None):
        r"""
        Construct a circle from a Vector representing the center, and the
        radius squared.
        """
        if base_ring is None:
            self._base_ring = radius_squared.parent()
        else:
            self._base_ring = base_ring

        # for calculations:
        self._V2 = VectorSpace(self._base_ring,2)
        self._V3 = VectorSpace(self._base_ring,3)

        self._center = self._V2(center)
        self._radius_squared = self._base_ring(radius_squared)
        
        
    def center(self):
        r"""
        Return the center of the circle as a vector.
        """
        return self._center
        
    def radius_squared(self):
        r"""
        Return the square of the radius of the circle.
        """
        return self._radius_squared
    
    def point_position(self, point):
        r"""
        Return 1 if point lies in the circle, 0 if the point lies on the circle,
        and -1 if the point lies outide the circle.
        """
        value = (point[0]-self._center[0])**2 + (point[1]-self._center[1])**2 - \
            self._radius_squared
        if value > self._base_ring.zero():
            return -1
        if value < self._base_ring.zero():
            return 1
        return 0
    
    def tangent_vector(self, point):
        r"""
        Return a vector based at the provided point (which must lie on the circle)
        which is tangent to the circle and points in the counter-clockwise
        direction.

        EXAMPLES::

            sage: from flatsurf.geometry.circle import Circle
            sage: c=Circle(vector((0,0)), 2, base_ring=QQ)
            sage: c.tangent_vector(vector((1,1)))
            (-1, 1)
        """
        if not self.point_position(point) == 0:
            raise ValueError("point not on circle.")
        return vector((self._center[1]-point[1], point[0]-self._center[0]))
    
    def other_intersection(self, p, v):
        r"""
        Consider a point p on the circle and a vector v. Let L be the line
        through p in direction v. Then L intersects the circle at another
        point q. This method returns q.
        
        Note that if p and v are both in the field of the circle,
        then so is q.
        
        EXAMPLES::

            sage: from flatsurf.geometry.circle import Circle
            sage: c=Circle(vector((0,0)), 25, base_ring=QQ)
            sage: c.other_intersection(vector((3,4)),vector((1,2)))
            (-7/5, -24/5)
        """
        pp=self._V3((p[0],p[1],self._base_ring.one()))
        vv=self._V3((v[0],v[1],self._base_ring.zero()))
        L = pp.cross_product(vv)
        cc=self._V3((self._center[0],self._center[1],self._base_ring.one()))
        vvperp=self._V3((-v[1],v[0],self._base_ring.zero()))
        # line perpendicular to L through center:
        Lperp = cc.cross_product(vvperp)
        # intersection of L and Lperp:
        rr = L.cross_product(Lperp)
        r = self._V2((rr[0]/rr[2],rr[1]/rr[2]))
        return self._V2((2*r[0]-p[0], 2*r[1]-p[1]))
    def __str__(self):
        return "circle with center "+str(self._center)+" and radius squared " + \
            str(self._radius_squared)
    
    def __repr__(self):
        return "Circle("+repr(self._center)+", "+repr(self._radius_squared)+")"
    
