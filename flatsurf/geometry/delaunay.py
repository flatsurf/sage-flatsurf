r"""
This file contains classes implementing Surface which are used useful for
triangulating, Delaunay triangulating, and Delaunay decomposing infinite
surfaces.
"""

from flatsurf.geometry.surface import Surface, Surface_list

class LazyTriangulatedSurface(Surface):
    r"""
    Surface class used to triangulate an infinite surface.
    
    EXAMPLES::

    Example with relabel=False::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyTriangulatedSurface(s,relabel=False))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run(skip="_test_pickling")

    Example with relabel=True::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyTriangulatedSurface(s,relabel=True))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run(skip="_test_pickling")
    """
    def __init__(self, similarity_surface, relabel=True):

        if similarity_surface.is_mutable():
            raise ValueError("Surface must be immutable.")
        
        # This surface will converge to the Delaunay Triangulation
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, \
            mutable=True)

        Surface.__init__(self, self._s.base_ring(), self._s.base_label(), \
            mutable=False, finite=self._s.is_finite())

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        polygon = self._s.polygon(lab)
        if polygon.num_edges() > 3:
            self._s.triangulate(in_place=True, label = lab)
            return self._s.polygon(lab)
        else:
            return polygon

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        pp,ee = self._s.opposite_edge(p,e)
        polygon = self._s.polygon(pp)
        if polygon.num_edges() > 3:
            self.polygon(pp)
            return self._s.opposite_edge(p,e)
        else:
            return (pp,ee)

class LazyDelaunayTriangulatedSurface(Surface):
    r"""
    Surface class used to find a Delaunay triangulation of an infinite surface.
    
    EXAMPLES::

    Example with relabel=False::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyDelaunayTriangulatedSurface(s,relabel=False))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run(skip="_test_pickling")

    Example with relabel=True::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.delaunay import *
        sage: s=translation_surfaces.infinite_staircase()
        sage: ss=TranslationSurface(LazyDelaunayTriangulatedSurface(s,relabel=True))
        sage: ss.polygon(0).num_edges()
        3
        sage: TestSuite(ss).run(skip="_test_pickling")
    """

    def _setup_direction(self, direction):
        # Our Delaunay will respect the provided direction.
        if direction is None:
            self._direction = self._s.vector_space()( \
                (self._s.base_ring().zero(), self._s.base_ring().one() ) )
        else:
            self._direction = self._ss.vector_space()(direction)

    def __init__(self, similarity_surface, direction=None, relabel=True):
        r"""
        Construct a lazy Delaunay triangulation of the provided similarity_surface.
        """
        if similarity_surface.underlying_surface().is_mutable():
            raise ValueError("Surface must be immutable.")

        # This surface will converge to the Delaunay Triangulation
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, mutable=True)

        self._setup_direction(direction)

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels=set()

        # We will now try to get the base_polygon.
        base_label=self._s.base_label()
        
        self._s.triangulate(in_place=True, label=base_label)
        failed=True
        while failed:
            failed=False
            polygon = self._s.polygon(base_label)
            circle = polygon.circumscribing_circle()
            #print polygon
            #print circle
            # Reversing order because the previously flipped edge is always zero
            for e in xrange(2,-1,-1):
                #print("e=",e)
                if self._certify_or_improve(circle,base_label,e) != 1:
                    #print("flipping base_label,",e,self._s.opposite_edge(base_label,e))
                    self._s.triangle_flip(base_label,e,in_place=True, \
                        direction=self._direction)
                    failed=True
                    break
        self._certified_labels.add(base_label)
        Surface.__init__(self, self._s.base_ring(), base_label, \
            finite=self._s.is_finite(), mutable=False)

    def polygon(self, label):
        if label in self._certified_labels:
            return self._s.polygon(label)
        else:
            raise ValueError("Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface.")
    
    def opposite_edge(self, label, edge):
        if label in self._certified_labels:
            ll,ee = self._s.opposite_edge(label,edge)
            if ll in self._certified_labels:
                return ll,ee
            #print("Searching for opposite to ",label,edge)
            result=-1
            while result!=1:
                ll,ee = self._s.opposite_edge(label,edge)
                pp=self._s.polygon(ll)
                circle=pp.circumscribing_circle()
                #print("oe ",pp,circle,ee)
                result = self._certify_or_improve(self._s.edge_transformation(ll,ee)*circle,label,edge)
            self._certified_labels.add(ll)
            return ll,ee
        else:
            raise ValueError("Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface.")
    
    def _certify_or_improve(self, circle, l, e):
        r"""
        This method attempts to develop the disk bounded by circle through the 
        edge (l,e). Here we are using the coordinates of the polygon with label l,
        and the interior of the disk intersected with the line continuing the edge
        is contained in the interior of the edge (l,e). 
        
        To continue the development from polygon(l), we pass through edge e
        entering a new polygon, polygon(ll) through edge ee. There are several possible
        outcomes:

        1) We may notice that the vertex opposite edge ee is contained in the interior
        of the disk. This means the disk can not be fully developed into the surface.
        We return -1. We then expect to do a "cascade of flips" and whatever method
        originally asked about the disk will be given a negative answer.
        
        2) We may notice that the edge we are passing through requires a diagonal flip
        to get closer to the Delaunay decomposition. We return 0 and expect the calling
        method to perform the flip. 
        
        3) It could be true that 1 and 2 don't occur but the disk might pass though
        some new edges of polygon(ll). In this case we proceed inductively, calling
        _certify_or_improve on these edges. If any of these methods return 0 we do a
        single flip and retry. If they return -1, then we perform the flip and pass 
        the message back (return -1). If all of these calls return 1, we have
        successfully developed the disk through the edge. We return 1.
        """
        failed=True
        while failed:
            failed=False
            ll,ee = self._s.opposite_edge(l,e)
            if ll in self._certified_labels:
                # Already certified.
                return 1
            pp=self._s.polygon(ll)
            if pp.num_edges()>3:
                # not triangulated!
                self._s.triangulate(in_place=True, label=ll)
                ll,ee = self._s.opposite_edge(l,e)
                pp=self._s.polygon(ll)
            # The polygon pp is now a triangle.
            c = self._s.edge_transformation(l,e)*circle
            #print("c_or_i",pp,c,ee)
            if c.point_position(pp.vertex((ee+2)%3))==1:
                # This causes a cascade of flips since the circle is not embedded.
                return -1
            if not l in self._certified_labels and self._s._edge_needs_flip(l,e):
                # Indicate that the edge needs to be flipped,
                # but no point in the circle yet.
                # Just one flip needed.
                return 0
            lp = c.line_segment_position(pp.vertex((ee+1)%3),pp.vertex((ee+2)%3))
            if lp==1:
                # Crosses through the circle.
                #print("crosses edge ",(ee+1)%3)
                result = self._certify_or_improve(c,ll,(ee+1)%3)
                if result != 1:
                    #print("flip edge ",ll,(ee+1)%3,self._s.opposite_edge(ll,(ee+1)%3))
                    self._s.triangle_flip(ll,(ee+1)%3,in_place=True, \
                        direction=self._direction)
                    if result == -1:
                        # Continue the cascade of flips:
                        return -1
                    # Otherwise just one flip is needed, but we need to restart
                    failed=True
                    continue
            lp = c.line_segment_position(pp.vertex((ee+2)%3),pp.vertex(ee))
            if lp==1:
                # Crosses through the circle.
                #print("crosses edge ",(ee+2)%3)
                result = self._certify_or_improve(c,ll,(ee+2)%3)
                if result != 1:
                    #print("flip edge ",ll,(ee+2)%3,self._s.opposite_edge(ll,(ee+1)%3))
                    self._s.triangle_flip(ll,(ee+2)%3,in_place=True, \
                        direction=self._direction)
                    if result == -1:
                        # Continue the cascade of flips:
                        return -1
                    # Otherwise just one flip is needed, but we need to restart
                    failed=True
                    continue
        # At this point we've successfully verified that this portion of the circle
        # embeds.
        return 1
        
class LazyDelaunaySurface(LazyDelaunayTriangulatedSurface):
    # We just inherit to use some methods.

    r"""
    This is an implementation of Surface. It takes a surface (typically
    infinite) from the constructor. This class respresents the 
    Delaunay decomposition of this surface. We compute this decomposition
    lazily so that it works for infinite surfaces.
    
    EXAMPLES::

        sage: from flatsurf import*
        sage: from flatsurf.geometry.delaunay import*
        sage: s=translation_surfaces.infinite_staircase()
        sage: m=matrix([[2,1],[1,1]])
        sage: ss=TranslationSurface(LazyDelaunaySurface(m*s,relabel=False))
        sage: ss.polygon(ss.base_label())
        Polygon: (0, 0), (0, -1), (1, -1), (1, 0)
        sage: TestSuite(ss).run(skip="_test_pickling")
    """


    def __init__(self, similarity_surface, direction=None, relabel=True):
        r"""
        Construct a lazy Delaunay triangulation of the provided similarity_surface.
        
        """
        if similarity_surface.underlying_surface().is_mutable():
            raise ValueError("Surface must be immutable.")

        # This surface will converge to the Delaunay Decomposition
        self._s = similarity_surface.copy(relabel=relabel, lazy=True, mutable=True)

        self._setup_direction(direction)

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels=set()

        base_label=self._s.base_label()

        # We will now try to get the base_polygon.
        self._s.triangulate(in_place=True, label=base_label)
        failed=True
        while failed:
            failed=False
            polygon = self._s.polygon(base_label)
            circle = polygon.circumscribing_circle()
            #print polygon
            #print circle
            # Reversing order because the previously flipped edge is always zero
            for e in xrange(2,-1,-1):
                #print("e=",e)
                if self._certify_or_improve(circle,base_label,e) != 1:
                    #print("flipping base_label,",e,self._s.opposite_edge(base_label,e))
                    self._s.triangle_flip(base_label,e,in_place=True, \
                        direction=self._direction)
                    failed=True
                    break

        self._join_adjacent_delaunay(base_label)
        self._certified_labels.add(base_label)
        
        Surface.__init__(self, self._s.base_ring(), base_label, \
            finite=self._s.is_finite(), mutable=False)

    def _join_adjacent_delaunay(self, label):
        r"""
        This glues adjacent polygons to this one when the adjacent polygon
        has the same circumscribing circle.
        """
        joined=True
        while joined:
            joined=False
            p = self._s.polygon(label)
            for e in xrange(p.num_edges()):
                if self._s._edge_needs_join(label,e):
                    self._s.join_polygons(label,e,in_place=True)
                    joined=True
                    break

    def polygon(self, label):
        if label in self._certified_labels:
            return self._s.polygon(label)
        else:
            raise ValueError("Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface.")
    
    def opposite_edge(self, label, edge):
        if label in self._certified_labels:
            ll,ee = self._s.opposite_edge(label,edge)
            if ll in self._certified_labels:
                return ll,ee

            #print("Searching for opposite to ",label,edge)
            # This does flips until the polygon opposite (label, edge) is also Delaunay
            result=-1
            while result!=1:
                ll,ee = self._s.opposite_edge(label,edge)
                pp=self._s.polygon(ll)
                circle=pp.circumscribing_circle()
                #print("oe ",pp,circle,ee)
                result = self._certify_or_improve( \
                    self._s.edge_transformation(ll,ee)*circle,label,edge)
                    
            self._join_adjacent_delaunay(ll)
            # Now the opposite edge has changed.
            ll,ee = self._s.opposite_edge(label,edge)
            self._certified_labels.add(ll)
            return ll,ee
        else:
            raise ValueError("Asked for polygon not known to be Delaunay. Make sure you obtain polygon labels by walking through the surface.")

