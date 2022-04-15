r"""
Two dimensional hyperbolic geometry.

EXAMPLES::

    sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

    sage: H2 = HyperbolicPlane(QQ)
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
#                      2022 Sam Freedman
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
######################################################################

from sage.structure.parent import Parent
from sage.all import QQ

class HyperbolicPlane(Parent):
    r"""
    The hyperbolic plane over a base ring.
   
    We do not use a fixed representation of the hyperbolic plane internally.
    However, we mostly think of this as the upper half plane in the Poincaré
    model.
    
    All objects in the plane must be specified over the given base ring. Note
    that, in some representations, objects might appear to live in a larger
    ring. E.g., when specifying a line by giving a midpoint and the square of
    its radius in the half plane model, then the ideal endpoints of this line
    might have coordinates in the ring after adjoining a square root.
    
    The implemented elements of the plane are convex subsets such as (finite
    and infinite) points, geodesics, closed half planes, and closed convex
    polygons.
    """
    def __init__(self, base_ring=QQ, category=None):
        raise NotImplementedError
        
    def infinity(self):
        r"""
        Return the point at infinity in the Poincaré half plane model.
        """
        return self.projective(0, 1)
        
    def real(self, r):
        r"""
        Return the ideal point ``r`` on the real axis in the Poincaré half
        plane model.
        """
        return self.projective(r, 1)
    
    def projective(self, p, q):
        r"""
        Return the ideal point with projective coordinates ``[p: q]`` in the
        Poincaré half plane model.
        """
        raise NotImplementedError 
        
    def half_circle(self, center, radius_squared):
        r"""
        Return the geodesic centered around the real ``center`` and with
        ``radius_squared`` in the Poincaré half plane model. The geodesic is
        oriented such that the point at infinity is to its left.
        
        Use the ``-`` operator to pass to the geodesic with opposite
        orientation.
        """
        raise NotImplementedError 
        
    def vertical(self, real):
        r"""
        Return the vertical geodesic at the ``real`` ideal point in the
        Poincaré half plane model. The geodesic is oriented such that it goes
        towards from ``real`` to the point at infinity.
        
        Use the ``-`` operator to pass to the geodesic with opposite
        orientation.
        """
        raise NotImplementedError 
        
    def chord(self, a, b):
        r"""
        Return the geodesic from the point on the unit circle whose argument is
        `a` to the point whose argument is `b` in the Klein model.
        
        Both `a` and `b` are understood as rational multiples of 2π, i.e., they
        are taken mod 1.
        """
        raise NotImplementedError
        
    def half_plane(self, geodesic):
        r"""
        Return the closed half plane that is on the left of ``geodesic``.
        
        Use the ``-`` operator to pass to the half plane on the right.
        """
        raise NotImplementedError
        
    def intersection(self, other):
        return other

        
class HyperbolicHalfPlane:
    r"""
    A closed half plane of the hyperbolic plane.
    
    Use :meth:`HyperbolicPlane.half_plane` to create a half plane.
    """
    def __init__(self, geodesic):
        self._geodesic = geodesic
        
    def __neg__(self):
        raise NotImplementedError
        
    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this half plane.
        """
        raise NotImplementedError
        
    def is_subset(self, half_plane):
        r"""
        Return whether this half plane is contained in ``half_plane``.
        """
        raise NotImplementedError
        
    def intersection(self, other):
        raise NotImplementedError

    
class HyperbolicPoint:
    r"""
    A (possibly infinite) point in the :class:`HyperbolicPlane`, namely the
    unique intersection point of the geodesics ``g`` and ``h``.
    """
    def __init__(self, parent, g, h):
        if g == h or g == -h:
            raise ValueError("geodesics must have a unique point of intersection")
            
        raise NotImplementedError
        
    def is_finite(self):
        r"""
        Return whether this is a finite point.
        """
        raise NotImplementedError
        
    def coordinates(self, model="half_plane", ring=None):
        r"""
        Return coordinates of this point in ``ring``.
        
        If ``model`` is ``"half_plane"``, return projective coordinates in the
        Poincaré half plane model.
        
        If ``model`` is ``"klein"``, return polar coordinates in the Klein model.
        
        If no ``ring`` has been specified, an appropriate extension of the base
        ring of the :class:`HyperbolicPlane` is chosen where these coordinates
        live.
        """
        raise NotImplementedError
        
    def intersection(self, other):
        r"""
        Return the intersection of this point and the convex object ``other``.
        """
        raise NotImplementedError

    
class HyperbolicConvexPolygon:
    r"""
    A (possibly unbounded) closed polygon in the :class:`HyperbolicPlane`,
    i.e., the intersection of a finite number of :class:`HyperbolicHalfPlane`s.
    """
    def __init__(self, parent, half_planes, assume_normalized=False):
        raise NotImplementedError
        
    def _normalize(self):
        r"""
        Normalize the internal list of half planes so that they describe the
        :meth:`boundary`.
        """
        raise NotImplementedError
        
    def intersection(self, other):
        r"""
        Return the intersection of this polygon with ``other`` where ``other``
        is another polygon or a half plane.
        """
        raise NotImplementedError
        
    def boundary(self):
        r"""
        Return the boundary of this polygon as a list of (oriented) geodesics.
        
        The output is minimal and sorted by angle in the Klein model.
        """
        raise NotImplementedError
        
    def vertices(self):
        r"""
        Return the vertices of this polygon, i.e., the points of intersection
        of the :meth:`boundary` geodesics.
        """
        raise NotImplementedError
        
    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this closed polygon.
        """
        raise NotImplementedError
        
    def intersection(self, other):
        r"""
        Return the intersection of this polygon and the convex object ``other``.
        """
        raise NotImplementedError
        

class HyperbolicGeodesic:
    r"""
    An oriented geodesic in the hyperbolic plane.
    
    We internally represent geodesics as the solutions to the equation `a(x^2 +
    y^2) + bx + c = 0` for `(x, y)` in the upper half plane.
    """
    def __init__(self, parent, a, b, c):
        raise NotImplementedError
        
    def intersection(self, other):
        r"""
        Return the (possibly infinite) point of intersection of this geodesic and ``other``.
        
        Returns the geodesic itself if they are equal as unoriented geodesics.
        """
        if self == other or self == -other:
            return self
        
        raise NotImplementedError
        
    def start(self):
        r"""
        Return the ideal starting point of this geodesic.
        
        Note that this is only possible if the radius of this geodesic is a
        square in the base ring of the :class:`HyperbolicPlane`.
        """
        raise NotImplementedError
        
    def end(self):
        r"""
        Return the ideal end point of this geodesic.
        
        Note that this is only possible if the radius of this geodesic is a
        square in the base ring of the :class:`HyperbolicPlane`.
        """
        return (-self).start()
   
    def _richcmp_(self, other, op):
        r"""
        Return how this geodesic compares to ``other``.
        
        Geodesics are partially ordered by their angle in [0, 1] in the Klein
        model, see :meth:`HyperbolicPlane.chord`.
        """
        
    def _neg_(self):
        raise NotImplementedError
        
    def intersection(self, other):
        r"""
        Return the intersection of this geodesic and the convex object ``other``.
        """
        raise NotImplementedError
        

class HyperbolicEdge:
    r"""
    An oriented (possibly infinite) segment in the hyperbolic plane such as a
    boundary edge of a :class:`HyperbolicConvexPolygon`.
    """
    def __init__(self, geodesic, start=None, end=None):
        raise NotImplementedError
        
    def intersection(self, other):
        r"""
        Return the intersection of this edge and the convex object ``other``.
        """
        raise NotImplementedError
