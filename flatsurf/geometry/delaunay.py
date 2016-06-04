r"""
This file contains methods useful for working with the Delaunay decomposition.
"""

from flatsurf.geometry.surface import Surface, Surface_list

class LazyTriangulatedSurface(Surface):
    def __init__(self, similarity_surface):
        self._ss = similarity_surface
        self._s = similarity_surface.underlying_surface()
        if self._s.is_mutable():
            raise ValueError("Surface must be immutable.")

        # This surface will converge to the Delaunay (but the user
        # will not know we don't have the whole thing already!)
        self._sl = Surface_list(base_ring = self._s.base_ring(), \
            surface=self._ss, copy=False, mutable=True)
        # Wrap surface._sl so we can use higher-level methods.
        self._ssl = self._ss.__class__(self._sl)
        
        Surface.__init__(self, mutable=False)
        
    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        return self._s.base_ring()

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        polygon = self._sl.polygon(lab)
        if polygon.num_edges() > 3:
            self._ssl.triangulate(in_place=True, label = lab)
            return self._sl.polygon(lab)
        else:
            return polygon

    def base_label(self):
        r"""
        Always returns the same label.
        """
        return 0

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        pp,ee = self._sl.opposite_edge(p,e)
        polygon = self._sl.polygon(pp)
        if polygon.num_edges() > 3:
            self.polygon(pp)
            return self._sl.opposite_edge(p,e)
        else:
            return (pp,ee)

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        return self._s.is_finite()

#def LazyDelaunayTriangulatedSurface(LazyTriangulatedSurface):
#    def __init__(self, similarity_surface, direction=None):
#        r"""
#        
#        Direction defaults to the upward unit vector.
#        """
#        self._ss = similarity_surface
#        self._s = similarity_surface.underlying_surface()
#        if self._s.is_mutable():
#            raise ValueError("Surface must be immutable.")

#        # This surface will converge to the Delaunay (but the user
#        # will not know we don't have the whole thing already!)
#        self._sl = Surface_list(base_ring = self._s.base_ring(), \
#            surface=self._ss, copy=False, mutable=True)
#        # Wrap surface._sl so we can use higher-level methods.
#        self._ssl = self._ss.__class__(self._sl)

#        # Our Delaunay will respect the provided direction.
#        if direction is None:
#            direction = self._s.vector_space()( \
#                (self._s.base_ring().zero(), self._s.base_ring().one() ) )
#        else:
#            self._direction = self._s.vector_space()(direction)

#        # We will now try to get the base_polygon.
#        self._ssl.triangulate(in_place=True, label=0)
#        polygon0 = self._sl.polygon(0)
#        circle = polygon0.circumscribing_circle()
#        loop = True
#        while loop:
#            loop = False
#            for e in xrange(polygon0.num_edges()):
#                pp,ee = self._sl.opposite_edge(p,e)
#                

#        Surface.__init__(self, mutable=False)

#    def __init__(self, similarity_surface):
#        
#        
#    def base_ring(self):
#        r"""
#        The field on which the coordinates of ``self`` live.

#        This method must be overriden in subclasses!
#        """
#        return self._s.base_ring()

#    def polygon(self, lab):
#        r"""
#        Return the polygon with label ``lab``.
#        """
#        polygon = self._sl.polygon(lab)
#        if polygon.num_edges() > 3:
#            self._ssl.triangulate(in_place=True, label = lab)
#            return self._sl.polygon(lab)
#        else:
#            return polygon

#    def base_label(self):
#        r"""
#        Always returns the same label.
#        """
#        return 0

#    def opposite_edge(self, p, e):
#        r"""
#        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
#        returns the pair (``pp``, ``ee``) to which this edge is glued.
#        """
#        pp,ee = self._sl.opposite_edge(p,e)
#        polygon = self._sl.polygon(pp)
#        if polygon.num_edges() > 3:
#            self.polygon(pp)
#            return self._sl.opposite_edge(p,e)
#        else:
#            return (pp,ee)

#    def is_finite(self):
#        r"""
#        Return whether or not the surface is finite.
#        """
#        return self._s.is_finite()

