r"""
This class contains methods useful for working with circles. 

This will be used to build a LazyDelaunayTriangulation class which will compute the
Delaunay decomposition for infinite surfaces.
"""
#*****************************************************************************
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector

def circle_from_three_points(p,q,r,base_ring=None):
    r"""
    Construct a circle from three points on the circle.
    """
    if base_ring is None:
        base_ring=p.base_ring()
    V2 = VectorSpace(base_ring.fraction_field(), 2)
    V3 = VectorSpace(base_ring.fraction_field(), 3)
    
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
    
    def closest_point_on_line(self, point, direction_vector):
        r"""
        Consider the line through the provided point in the given direction.
        Return the closest point on this line to the center of the circle.
        """
        cc = self._V3((self._center[0],self._center[1], self._base_ring.one()))
        # point at infinite orthogonal to direction_vector:
        dd = self._V3((direction_vector[1],-direction_vector[0], self._base_ring.zero()))
        l1 = cc.cross_product(dd)
        
        pp = self._V3((point[0],point[1], self._base_ring.one()))
        # direction_vector pushed to infinity
        ee = self._V3((direction_vector[0],direction_vector[1], self._base_ring.zero()))
        l2 = pp.cross_product(ee)
        
        # This is the point we want to return
        rr = l1.cross_product(l2)
        try:
            return self._V2((rr[0]/rr[2], rr[1]/rr[2]))
        except ZeroDivisionError:
            raise ValueError("Division by zero error. Perhaps direction is zero. "+\
                "point="+str(point)+" direction="+str(direction_vector)+" circle="+\
                str(self))


    def line_position(self, point, direction_vector):
        r"""
        Consider the line through the provided point in the given direction.
        We return 1 if the line passes through the circle, 0 if it is tangent
        to the circle and -1 if the line does not intersect the circle.
        """
        return self.point_position(self.closest_point_on_line(point,direction_vector))

    def line_segment_position(self, p, q):
        r"""
        Consider the open line segment pq.We return 1 if the line segment
        enters the interior of the circle, zero if it touches the circle 
        tangentially (at a point in the interior of the segment) and
        and -1 if it does not touch the circle or its interior.
        """
        if self.point_position(p)==1:
            return 1
        if self.point_position(q)==1:
            return 1
        r=self.closest_point_on_line(p,q-p)
        pos = self.point_position(r)
        if pos ==-1:
            return -1
        # This checks if r lies in the interior of pq
        if p[0]==q[0]:
            if (p[1]<r[1] and r[1]<q[1]) or (p[1]>r[1] and r[1]>q[1]):
                return pos
        elif (p[0]<r[0] and r[0]<q[0]) or (p[0]>r[0] and r[0]>q[0]):
            return pos
        # It does not lie in the interior.
        return -1

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

    def __rmul__(self, similarity):
        r"""
        Apply a similarity to the circle.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.circle import *
            sage: s = translation_surfaces.square_torus()
            sage: c = s.polygon(0).circumscribing_circle()
            sage: c
            Circle((1/2, 1/2), 1/2)
            sage: s.edge_transformation(0,2)
            (x, y) |-> (x, y - 1)
            sage: s.edge_transformation(0,2) * c
            Circle((1/2, -1/2), 1/2)
        """
        from .similarity import SimilarityGroup
        SG = SimilarityGroup(self._base_ring)
        s = SG(similarity)
        return Circle(s(self._center), \
            s.det()*self._radius_squared, \
            base_ring=self._base_ring)

    def __str__(self):
        return "circle with center "+str(self._center)+" and radius squared " + \
            str(self._radius_squared)
    
    def __repr__(self):
        return "Circle("+repr(self._center)+", "+repr(self._radius_squared)+")"
    
