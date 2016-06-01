r"""
Similarity surfaces.
"""

from sage.misc.cachefunc import cached_method

from sage.structure.sage_object import SageObject

from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.infinity import Infinity

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.qqbar import AA
from sage.rings.real_mpfr import RR
from sage.rings.real_mpfi import RIF

from sage.modules.free_module_element import vector

from sage.matrix.constructor import matrix, identity_matrix

ZZ_1 = Integer(1)
ZZ_2 = Integer(2)

from flatsurf import *

from flatsurf.geometry.matrix_2x2 import (is_similarity,
                    homothety_rotation_decomposition,
                    similarity_from_vectors,
                    rotation_matrix_angle,
                    is_cosine_sine_of_rational)
                    
from flatsurf.geometry.similarity import SimilarityGroup
from flatsurf.geometry.polygon import Polygons, wedge_product
from flatsurf.geometry.surface import Surface, LabelWalker, ExtraLabel

class SimilaritySurface(SageObject):
    r"""
    An oriented surface built from a set of polygons and edges identified with
    similarities (i.e. composition of homothety, rotations and translations).

    Each polygon is identified with a unique key (its label). The choice of the
    label of the polygons is done at startup. If the set is finite then by
    default the labels are the first non-negative integers 0,1,...

    The edge are identified by a couple (polygon label, edge number).

    .. NOTE::

        This class is abstract and should not be called directly. Instead you
        can either use SimilaritySurface_from_polygons_and_identifications or
        inherit from SimilaritySurface_generic and implement the methods:

        - base_ring(self): the base ring in which coordinates lives
        - polygon(self, lab): the polygon associated to the label ``lab``
        - base_label(self): return a first label
        - opposite_edge(self, lab, edge): a couple (``other_label``, ``other_edge``) representing the edge being glued
        - is_finite(self): return true if the surface is built from finitely many labeled polygons
    """
    
    def __init__(self, surface):
        if isinstance(surface,SimilaritySurface):
            self._s=surface.underlying_surface()
        else:
            self._s=surface
        assert isinstance(self._s,Surface)

    def underlying_surface(self):
        r"""
        Return the surface underlying this SimilaritySurface.
        """
        return self._s

    def _check(self):
        r"""
        DEPRECATED

        Just use the standard test suite to implement tests

        EXAMPLES::

            sage: from flatsurf import *
            sage: n = 6
            sage: ps = [polygons.regular_ngon(2*n)]
            sage: gluings = [((0,i),(0,i+n)) for i in range(n)]
            sage: s = TranslationSurface(Surface_polygons_and_gluings(ps,gluings))
            sage: TestSuite(s).run(verbose=True)
            running ._test_category() . . . pass
            running ._test_edge_matrix() . . . pass
            running ._test_gluings() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
        """
        from sage.misc.superseded import deprecation
        deprecation(33, "Just use TestSuite!...")
        from sage.misc.sage_unittest import TestSuite
        TestSuite(self).run()

    def _test_gluings(self, **options):
        # iterate over pairs with pair1 glued to pair2
        tester = self._tester(**options)
        if self.is_finite():
            it = self.edge_iterator(gluings=True)
        else:
            from itertools import islice
            it = islice(self.edge_iterator(gluings=True), 30)

        for pair1,pair2 in it:
            tester.assertEqual(self.opposite_edge(pair2), pair1,
                "edges not glued correctly:\n%s -> %s -> %s"%(pair1,pair2,self.opposite_edge(pair2)))

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
        return self._s.polygon(lab)

    def base_label(self):
        r"""
        Always returns the same label.
        """
        return self._s.base_label()

    def opposite_edge(self, p, e=None):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        if e is None:
            return self._s.opposite_edge(p[0],p[1])
        return self._s.opposite_edge(p,e)

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        return self._s.is_finite()

    def is_mutable(self):
        r"""
        Return if the underlying surface is mutable.
        """
        return self._s.is_mutable()

    # 
    # generic methods
    #
    
    #def compute_surface_type_from_gluings(self,limit=None):
    #    r"""
    #    Compute the surface type by looking at the edge gluings. 
    #    If limit is defined, we try to guess the type by looking at limit many edges.
    #    """
    #    if limit is None:
    #        if not self.is_finite():
    #            raise ValueError("Need a limit when working with an infinite surface.")
    #        it = self.edge_iterator()
    #        label,edge = it.next()
    #        # Use honest matrices!
    #        m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #        surface_type = surface_type_from_matrix(m)
    #        for label,edge in it:
    #            # Use honest matrices!
    #            m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #            surface_type = combine_surface_types(surface_type, surface_type_from_matrix(m))
    #        return surface_type
    #    else:
    #        count=0
    #        it = self.edge_iterator()
    #        label,edge = it.next()
    #        # Use honest matrices!
    #        m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #        surface_type = surface_type_from_matrix(m)
    #        for label,edge in it:
    #            # Use honest matrices!
    #            m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #            surface_type = combine_surface_types(surface_type, surface_type_from_matrix(m))
    #            count=count+1
    #            if count >= limit:
    #                return surface_type
    #        return surface_type

    def walker(self):
        return self._s.walker()

    def label_iterator(self, polygons=False):
        r"""
        Iterator over all polygon labels.
        
        If the keyword polygons is True then we return pairs (label, polygon)
        instead of just labels.
        """
        if polygons:
            return self._s.label_polygon_iterator()
        else:
            return self._s.label_iterator()

    def edge_iterator(self, gluings=False):
        r"""
        Iterate over the edges of polygons, which are pairs (l,e) where l is a polygon label, 0 <= e < N and N is the number of edges of the polygon with label l.
        
        If the keyword gluings is set to true, then we iterate over ordered 
        pairs of edges ((l,e),(ll,ee)) where edge (l,e) is glued to (ll,ee).
        
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import Polygons
            sage: P=Polygons(QQ)
            sage: tri0=P([(1,0),(0,1),(-1,-1)])
            sage: tri1=P([(-1,0),(0,-1),(1,1)])
            sage: gluings=[((0,0),(1,0)),((0,1),(1,1)),((0,2),(1,2))]
            sage: from flatsurf.geometry.surface import Surface_polygons_and_gluings
            sage: from flatsurf.geometry.translation_surface import TranslationSurface
            sage: s=TranslationSurface(Surface_polygons_and_gluings([tri0,tri1], gluings))
            sage: for edge in s.edge_iterator():
            ...       print edge
            (0, 0)
            (0, 1)
            (0, 2)
            (1, 0)
            (1, 1)
            (1, 2)
        """
        if gluings:
            return self._s.edge_gluing_iterator()
        else:
            return self._s.edge_iterator()

    def num_polygons(self):
        r"""
        Return the number of polygons.
        """
        return self._s.num_polygons()

    def num_edges(self):
        r"""
        Return the total number of edges of all polygons used.
        """
        return self._s.num_edges()
        
    def num_singularities(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import *

            sage: translation_surfaces.regular_octagon().num_singularities()
            1

            sage: S = SymmetricGroup(4)
            sage: r = S('(1,2)(3,4)')
            sage: u = S('(2,3)')
            sage: translation_surfaces.origami(r,u).num_singularities()
            2

            sage: S = SymmetricGroup(8)
            sage: r = S('(1,2,3,4,5,6,7,8)')
            sage: u = S('(1,8,5,4)(2,3)(6,7)')
            sage: translation_surfaces.origami(r,u).num_singularities()
            4
        """
        if not self.is_finite():
            raise ValueError("the method only work for finite surfaces")

        # NOTE:
        # the very same code is implemented in the method angles (translation
        # surfaces). we should factor out the code
        edges = set((p,e) for p in self.label_iterator() for e in range(self.polygon(p).num_edges()))

        n = ZZ(0)
        while edges:
            p,e = edges.pop()
            n += 1
            ee = (e-1) % self.polygon(p).num_edges()
            p,e = self.opposite_edge(p,ee)
            while (p,e) in edges:
                edges.remove((p,e))
                ee = (e-1) % self.polygon(p).num_edges()
                p,e = self.opposite_edge(p,ee)
        return n

    def _repr_(self):
        if self.num_polygons() == Infinity:
            num = 'infinitely many'
        else:
            num = str(self.num_polygons())

        if self.num_polygons() == 1:
            end = ""
        else:
            end = "s"

        return "{} built from {} polygon{}".format(self.__class__.__name__, num, end)

    def edge_matrix(self, p, e=None):
        r"""
        Return the edge to which this edge is identified and the matrix to be
        applied.
        """
        if e is None:
            p,e = p
        u = self.polygon(p).edge(e)
        pp,ee = self.opposite_edge(p,e)
        v = self.polygon(pp).edge(ee)

        # be careful, because of the orientation, it is -v and not v
        return similarity_from_vectors(u,-v)

    def edge_transformation(self, p, e):
        r"""
        Return the similarity bringing the provided edge to the opposite edge.

        EXAMPLES::
        
            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: print(s.polygon(0))
            Polygon: (0, 0), (2, -2), (2, 0)
            sage: print(s.polygon(1))
            Polygon: (0, 0), (2, 0), (1, 3)
            sage: print(s.opposite_edge(0,0))
            (1, 1)
            sage: g = s.edge_transformation(0,0)
            sage: g((0,0))
            (1, 3)
            sage: g((2,-2))
            (2, 0)
        """
        q=self.polygon(p)
        G=SimilarityGroup(self.base_ring())
        a=self.polygon(p).vertex(e)
        b=self.polygon(p).vertex(e+1)
        # This is the similarity carrying the origin to a and (1,0) to b:
        g=G(b[0]-a[0],b[1]-a[1],a[0],a[1])

        pp,ee = self.opposite_edge(p,e)
        qq=self.polygon(pp)
        # Be careful here: opposite vertices are identified
        aa=qq.vertex(ee+1)
        bb=qq.vertex(ee)
        # This is the similarity carrying the origin to aa and (1,0) to bb:
        gg=G(bb[0]-aa[0],bb[1]-aa[1],aa[0],aa[1])

        # This is the similarity carrying (a,b) to (aa,bb):
        return gg*(~g)

    def mutable_copy(self, dictionary=True):
        r"""
        Returns a mutable copy of this surface.
        
        If dictionary is false, labels will be changed, but the resulting
        surface will be slightly more efficient (as a list will be used
        for storing polygons rather than a dictionary). See Surface_fast.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: ss=translation_surfaces.ward(3)
            sage: print(ss.is_mutable())
            False
            sage: from flatsurf.geometry.surface import Surface_fast
            sage: s=ss.mutable_copy()
            sage: print(s.is_mutable())
            True
            sage: TestSuite(s).run()
            sage: print(s==ss)
            True
        """
        if self.is_finite():
            from flatsurf.geometry.surface import Surface_fast
            return self.__class__(Surface_fast(surface=self,\
                mutable=True,dictionary=dictionary))
        else:
            raise NotImplementedError("Mutable copy not implemented yet for infinite surfaces.")
    
    def triangle_flip(self, l1, e1, in_place=False, test=False):
        r"""
        Flips the diagonal of the quadrilateral formed by two triangles
        glued together along the provided edge (l1,e1). This can be broken 
        into two steps: join along the edge to form a convex quadilateral,
        then cut along the other diagonal. Raises a ValueError if this 
        quadrilateral would be non-convex.
        
        Parameters
        ----------
        l1
            label of polygon
        e1 : integer
            edge of the polygon
        in_place : boolean
            if True do the flip to the current surface which must be mutable. 
            In this case the updated surface will be returned.
            Otherwise a mutable copy is made and then an edge is flipped, which is then returned.
        test : boolean
            if True we don't actually flip, and we return True or False depending
            on whether or not the flip would be successful.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: s=similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
            sage: print(s.polygon(0))
            Polygon: (0, 0), (1, 0), (0, 1)
            sage: print s.triangle_flip(0, 0, test=True)
            False
            sage: print s.triangle_flip(0, 1, test=True)
            True
            sage: print s.triangle_flip(0, 2, test=True)
            False

            sage: from flatsurf import *
            sage: s=similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
            sage: from flatsurf.geometry.surface import Surface_fast
            sage: s=s.__class__(Surface_fast(s, mutable=True))
            sage: try:
            ....:     s.triangle_flip(0,0,in_place=True)
            ....: except ValueError as e:
            ....:     print e
            Gluing triangles along this edge yields a non-convex quadrilateral.
            sage: print s.triangle_flip(0,1,in_place=True)
            ConeSurface built from 2 polygons
            sage: print s.polygon(0)
            Polygon: (0, 0), (1, 1), (0, 1)
            sage: print s.polygon(1)
            Polygon: (0, 0), (-1, -1), (0, -1)
            sage: for p in s.edge_iterator(gluings=True):
            ....:     print p
            ((0, 0), (1, 0))
            ((0, 1), (0, 2))
            ((0, 2), (0, 1))
            ((1, 0), (0, 0))
            ((1, 1), (1, 2))
            ((1, 2), (1, 1))
            sage: try:
            ....:     s.triangle_flip(0,2,in_place=True)
            ....: except ValueError as e:
            ....:     print e
            ....: 
            Gluing triangles along this edge yields a non-convex quadrilateral.

            sage: from flatsurf import *
            sage: p=polygons((2,0),(-1,3),(-1,-3))
            sage: s=similarity_surfaces.self_glued_polygon(p)
            sage: from flatsurf.geometry.surface import Surface_fast
            sage: s=s.__class__(Surface_fast(s,mutable=True))
            sage: s.triangle_flip(0,1,in_place=True)
            HalfTranslationSurface built from 1 polygon
            sage: for x in s.label_iterator(polygons=True):
            ....:     print x
            (0, Polygon: (0, 0), (3, 3), (1, 3))
            sage: for x in s.edge_iterator(gluings=True):
            ....:     print x
            ((0, 0), (0, 0))
            ((0, 1), (0, 1))
            ((0, 2), (0, 2))
            sage: TestSuite(s).run()
        """
        if test:
            # Just test if the flip would be successful
            p1=self.polygon(l1)
            if not p1.num_edges()==3:
                return false
            l2,e2 = self.opposite_edge(l1,e1)
            p2 = self.polygon(l2)
            if not p2.num_edges()==3:
                return false
            sim = self.edge_transformation(l2,e2)
            hol = sim( p2.vertex( (e2+2)%3 ) - p1.vertex((e1+2)%3) )
            from flatsurf.geometry.polygon import wedge_product
            return wedge_product(p1.edge((e1+2)%3), hol) > 0 and \
                wedge_product(p1.edge((e1+1)%3), hol) > 0
        if in_place:
            s=self.underlying_surface()
        else:
            from flatsurf.geometry.surface import Surface_fast
            s=Surface_fast(surface=self.underlying_surface(), \
                mutable=True, dictionary=True)
        p1=self.polygon(l1)
        if not p1.num_edges()==3:
            raise ValueError("The polygon with the provided label is not a triangle.")
        l2,e2 = self.opposite_edge(l1,e1)
        p2 = self.polygon(l2)
        if not p2.num_edges()==3:
            raise ValueError("The polygon opposite the provided edge is not a triangle.")
        sim = self.edge_transformation(l2,e2)
        m = sim.derivative()
        hol = sim( p2.vertex( (e2+2)%3 ) - p1.vertex((e1+2)%3) )
        from flatsurf import polygons
        # The new polygons
        #print [hol, m*p2.edge((e2+2)%3), p1.edge((e1+1)%3)]
        #print [-hol, p1.edge((e1+2)%3), m*p2.edge((e2+1)%3)]
        try:
            np1 = polygons(edges=[hol, m * p2.edge((e2+2)%3), p1.edge((e1+1)%3)])
            np2 = polygons(edges=[-hol, p1.edge((e1+2)%3), m * p2.edge((e2+1)%3)])
        except (ValueError, TypeError):
            raise ValueError("Gluing triangles along this edge yields a non-convex quadrilateral.")
        # Old gluings:
        pairs = [self.opposite_edge(l2,(e2+2)%3), \
            self.opposite_edge(l1,(e1+1)%3), \
            self.opposite_edge(l1,(e1+2)%3), \
            self.opposite_edge(l2,(e2+1)%3)]
        for i, (l,e) in enumerate(pairs):
            if l==l1:
                if e==(e1+1)%3:
                    pairs[i]=(l1,2)
                elif e==(e1+2)%3:
                    pairs[i]=(l2,1)
                else:
                    raise ValueError("Surfaced passed has errors in polygon gluings.")
            elif l==l2:
                if e==(e2+1)%3:
                    pairs[i]=(l2,2)
                elif e==(e2+2)%3:
                    pairs[i]=(l1,1)
                else:
                    raise ValueError("Surfaced passed has errors in polygon gluings.")
        if l1==l2:
            s.change_polygon(l1,np1)
            s.change_edge_gluing(l1,0,l1,0)
            s.change_edge_gluing(l1,1,pairs[0][0],pairs[0][1])
            s.change_edge_gluing(l1,2,pairs[1][0],pairs[1][1])
        else:
            s.change_polygon(l1,np1)
            s.change_polygon(l2,np2)
            s.change_edge_gluing(l1,0,l2,0)
            s.change_edge_gluing(l1,1,pairs[0][0],pairs[0][1])
            s.change_edge_gluing(l1,2,pairs[1][0],pairs[1][1])
            s.change_edge_gluing(l2,1,pairs[2][0],pairs[2][1])
            s.change_edge_gluing(l2,2,pairs[3][0],pairs[3][1])
        if in_place:
            return self
        else:
            return self.__class__(s)

    def join_polygons(self, p1, e1, test=False, in_place=False):
        r"""
        Join polygons across the provided edge (p1,e1). By default,
        it returns the surface obtained by joining the two polygons 
        together. It raises a ValueError if gluing the two polygons
        together results in a non-convex polygon. This is done to the 
        current surface if in_place is True, and otherwise a mutable 
        copy is made and then modified. 
        
        If test is True then instead of changing the surface, it just
        checks to see if the change would be successful and returns
        True if successful or False if not.
        
        EXAMPLES::

            sage: from flatsurf import *
            sage: ss=translation_surfaces.ward(3)
            sage: from flatsurf.geometry.surface import Surface_fast
            sage: s=ss.mutable_copy(dictionary=False)
            sage: s.join_polygons(0,0, in_place=True)
            TranslationSurface built from 2 polygons
            sage: print(s.polygon(0))
            Polygon: (0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (0, 2*a), (-1, a)
            sage: s.join_polygons(0,4, in_place=True)
            TranslationSurface built from 1 polygon
            sage: print(s.polygon(0))
            Polygon: (0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (1, 3*a), (0, 2*a), (-1, a)
        """
        poly1=self.polygon(p1)
        p2,e2 = self.opposite_edge(p1,e1)
        poly2=self.polygon(p2)
        if p1==p2:
            if test:
                return False
            else:
                raise ValueError("Can't glue polygon to itself.")
        t=self.edge_transformation(p2, e2)
        dt=t.derivative()
        vs = []
        edge_map={} # Store the pairs for the old edges.
        for i in range(e1):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))
        ne=poly2.num_edges()
        for i in range(1,ne):
            ee=(e2+i)%ne
            edge_map[len(vs)]=(p2,ee)
            vs.append(dt * poly2.edge( ee ))
        for i in range(e1+1, poly1.num_edges()):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))

        from flatsurf.geometry.polygon import Polygons
        try:
            new_polygon = Polygons(self.base_ring())(vs)
        except (ValueError, TypeError):
            if test:
                return False
            else:
                raise ValueError("Joining polygons along this edge results in a non-convex polygon.")
        
        if test:
            # Gluing would be successful
            return True
        
        # Now no longer testing. Do the gluing.
        if in_place:
            ss=self
        else:
            ss=self.mutable_copy()
        s=ss.underlying_surface()

        inv_edge_map={}
        for key, value in edge_map.iteritems():
            inv_edge_map[value]=(p1,key)
        
        glue_list=[]
        for i in range(len(vs)):
            p3,e3 = edge_map[i]
            p4,e4 = self.opposite_edge(p3,e3)
            if p4 == p1 or p4 == p2: 
                glue_list.append(inv_edge_map[(p4,e4)])
            else:
                glue_list.append((p4,e4))

        if s.base_label()==p2:
             s.change_base_polygon(p1)
        s.remove_polygon(p2)
        
        s.change_polygon(p1, new_polygon, glue_list)
        
        return ss

    def subdivide_polygon(self, p, v1, v2, new_label = None, test=False):
        r"""
        Cut the polygon with label p along the diagonal joining vertex
        v1 to vertex v2. This cuts p into two polygons, one will keep the same
        label. The other will get a new label, which can be provided
        via new_label. Otherwise a default new label will be provided.
        If test=False, then the surface will be changed (in place). If
        test=True, then it just checks to see if the change would be successful
        
        The change will be done in place.
        """        
        poly=self.polygon(p)
        ne=poly.num_edges()
        if v1<0 or v2<0 or v1>=ne or v2>=ne:
            if test:
                return False
            else:
                raise ValueError('Provided vertices out of bounds.')
        if abs(v1-v2)<=1 or abs(v1-v2)>=ne-1:
            if test:
                return False
            else:
                raise ValueError('Provided diagonal is not a diagonal.')

        if v2<v1:
            temp=v1
            v1=v2
            v2=temp
            
        newvertices1=[poly.vertex(v2)-poly.vertex(v1)]
        for i in range(v2, v1+ne):
            newvertices1.append(poly.edge(i))
        newpoly1 = Polygons(self.base_ring())(newvertices1)
        
        newvertices2=[poly.vertex(v1)-poly.vertex(v2)]
        for i in range(v1,v2):
            newvertices2.append(poly.edge(i))
        newpoly2 = Polygons(self.base_ring())(newvertices2)

        if new_label is None:
            new_label = self.underlying_surface().add_polygon(None)
            
        old_to_new_labels={}
        for i in range(ne):
            if i<v1:
                old_to_new_labels[i]=(p,i+ne-v2+1)
            elif i<v2:
                old_to_new_labels[i]=(new_label,i-v1+1)
            else: # i>=v2
                old_to_new_labels[i]=(p,i-v2+1)
        new_to_old_labels={}
        for i,pair in old_to_new_labels.iteritems():
            new_to_old_labels[pair]=i

        glue_dictionary = {(p,0):(new_label,0)}

        glue1 = [None for i in xrange(newpoly1.num_edges())]
        glue1[0]=(new_label,0)
        glue2 = [None for i in xrange(newpoly2.num_edges())]
        glue2[0]=(p,0)
        for e in range(ne):
            ll,ee = old_to_new_labels[e]
            lll,eee = self.opposite_edge(p,e)
            if lll == p:
                if ll==p:
                    glue1[ee]=old_to_new_labels[eee]
                else:
                    glue2[ee]=old_to_new_labels[eee]
            else:
                if ll==p:
                    glue1[ee]=(lll,eee)
                else:
                    glue2[ee]=(lll,eee)
        
        self.underlying_surface().change_polygon(new_label, newpoly2, glue2)
        self.underlying_surface().change_polygon(p, newpoly1, glue1)

    def minimal_translation_cover(self):
        r"""
        Return the minimal translation cover.

        "Be careful that if the surface is not built from one polygon, this is
        not the smallest translation cover of the surface." - Vincent 
        
        "I disagree with the prior statement. Can you provide an example?" -Pat

        EXAMPLES::

            sage: from flatsurf import *
            sage: S = similarity_surfaces.example()
            sage: T = S.minimal_translation_cover()
            sage: T
            TranslationSurface built from infinitely many polygons
            sage: T.polygon(T.base_label())
            Polygon: (0, 0), (2, -2), (2, 0)
        """
        from flatsurf.geometry.translation_surface import MinimalTranslationCover, TranslationSurface
        return TranslationSurface(MinimalTranslationCover(self))

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), 2)

    def fundamental_group(self, base_label=None):
        r"""
        Return the fundamental group of this surface.
        """
        if not self.is_finite():
            raise ValueError("the method only work for finite surfaces")
        if base_label is None:
            base_label = self.base_label()
        from fundamental_group import FundamentalGroup
        return FundamentalGroup(self, base_label)

    def tangent_bundle(self, ring=None):
        r"""
        Return the tangent bundle

        INPUT:

        - ``ring`` -- an optional field (defaults to the coordinate field of the
          surface)
        """
        if ring is None:
            ring = self.base_ring()

        try:
            return self._tangent_bundle_cache[ring]
        except AttributeError:
            self._tangent_bundle_cache = {}
        except KeyError:
            pass

        from tangent_bundle import SimilaritySurfaceTangentBundle
        self._tangent_bundle_cache[ring] = SimilaritySurfaceTangentBundle(self, ring)
        return self._tangent_bundle_cache[ring]

    def tangent_vector(self, lab, p, v, ring=None):
        r"""
        Return a tangent vector.

        INPUT:

        - ``lab`` -- label of a polygon

        - ``p`` -- coordinates of a point in the polygon

        - ``v`` -- coordinates of a vector in R^2
        
        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S.tangent_vector(S.base_label(), (1/2,1/2), (1,1))
            SimilaritySurfaceTangentVector in polygon (1, [-1  0]
            [ 0 -1]) based at (1/2, -3/2) with vector (1, 1)
            sage: K.<sqrt2> = QuadraticField(2)
            sage: S.tangent_vector(S.base_label(), (1/2,1/2), (1,sqrt2))
            SimilaritySurfaceTangentVector in polygon (1, [-1  0]
            [ 0 -1]) based at (1/2, -3/2) with vector (1, sqrt2)
        """
        p = vector(p)
        v = vector(v)

        if p.parent().dimension() != 2 or v.parent().dimension() != 2:
            raise ValueError("p (={!r}) and v (={!v}) should have two coordinates")

        if ring is None:
            R = p.base_ring()
            if R != v.base_ring():
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                R = cm.common_parent(R, v.base_ring())
                p = p.change_ring(R)
                v = v.change_ring(R)
    
            R2 = self.base_ring()
            if R != R2:
                if R2.has_coerce_map_from(R):
                    p = p.change_ring(R2)
                    v = v.change_ring(R2)
                    R = R2
                elif not R.has_coerce_map_from(R2):
                    raise ValueError("not able to find a common ring for arguments")
            return self.tangent_bundle(R)(lab, p, v)
        else:
            return self.tangent_bundle(ring)(lab, p, v)
    
    
    
    def triangulation_mapping(self):
        r"""
        Return a SurfaceMapping triangulating the suface or None if the surface is already triangulated.
        """
        from flatsurf.geometry.mappings import triangulation_mapping
        return triangulation_mapping(self)
    
    def triangulate(self, in_place=False):
        r"""
        Return a triangulated version of this surface.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.mcmullen_L(1,1,1,1)
            sage: ss=s.triangulate()
            sage: gs=ss.graphical_surface()
            sage: gs.make_all_visible()
            sage: print(gs)
            Graphical version of Similarity Surface TranslationSurface built from 6 polygons
        """
        if in_place:
            s=self
        else:
            s=self.mutable_copy()
        from flatsurf.geometry.polygon import wedge_product
        loop=True
        while loop:
            loop=False
            for l,poly in s.label_iterator(polygons=True):
                n = poly.num_edges() 
                if n>3:
                    for i in xrange(n):
                        e1=poly.edge(i)
                        e2=poly.edge((i+1)%n)
                        if wedge_product(e1,e2) != 0:
                            s.subdivide_polygon(l,i,(i+2)%n)
                            loop=True
                            break
                    if loop:
                        break
                    else:
                        # This should never happen:
                        raise ValueError("Unable to triangulate polygon with label "+ \
                            str(l)+": "+str(poly))
        return s
    
    def _edge_needs_flip(self,p1,e1):
        r"""
        Returns -1 if the the provided edge incident to two triangles which 
        should be flipped to get closer to the Delaunay decomposition. 
        Returns 0 if the quadrilateral formed by the triangles is inscribed 
        in a circle, and returns 1 otherwise.
        
        A ValueError is raised if the edge is not indident to two triangles.
        """
        p2,e2=self.opposite_edge(p1,e1)
        poly1=self.polygon(p1)
        poly2=self.polygon(p2)
        if poly1.num_edges()!=3 or poly2.num_edges()!=3:
            raise ValueError("Edge must be adjacent to two triangles.")
        from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
        sim1=similarity_from_vectors(poly1.edge(e1+2),-poly1.edge(e1+1))
        sim2=similarity_from_vectors(poly2.edge(e2+2),-poly2.edge(e2+1))
        sim=sim1*sim2
        return sim[1][0]<0

    def _edge_needs_join(self,p1,e1):
        r"""
        Returns -1 if the the provided edge incident to two triangles which 
        should be flipped to get closer to the Delaunay decomposition. 
        Returns 0 if the quadrilateral formed by the triangles is inscribed 
        in a circle, and returns 1 otherwise.
        
        A ValueError is raised if the edge is not indident to two triangles.
        """
        p2,e2=self.opposite_edge(p1,e1)
        poly1=self.polygon(p1)
        poly2=self.polygon(p2)
        from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
        sim1=similarity_from_vectors(poly1.vertex(e1) - poly1.vertex(e1+2),\
            -poly1.edge(e1+1))
        sim2=similarity_from_vectors(poly2.vertex(e2) - poly2.vertex(e2+2),\
            -poly2.edge(e2+1))
        sim=sim1*sim2
        from sage.functions.generalized import sgn
        return sim[1][0]==0
    
    def delaunay_triangulation(self, triangulated=False, in_place=False):
        if not self.is_finite():
            raise NotImplementedError("Not implemented for infinite surfaces.")
        if triangulated:
            if in_place:
                s=self
            else:
                from flatsurf.geometry.surface import Surface_fast
                s=self.__class__(Surface_fast(self,mutable=True))
        else:
            from flatsurf.geometry.surface import Surface_fast
            s=self.__class__(Surface_fast(self.triangulate(in_place=in_place),mutable=True))
        loop=True
        while loop:
            loop=False
            for (l1,e1),(l2,e2) in s.edge_iterator(gluings=True):
                if (l1<l2 or (l1==l2 and e1<=e2)) and s._edge_needs_flip(l1,e1):
                    s.triangle_flip(l1, e1, in_place=True)
                    loop=True
                    break
        return s
    
    def delaunay_decomposition(self, triangulated=False, \
            delaunay_triangulated=False, in_place=False):
        r"""
        Return the Delaunay Decomposition of this surface.
    
        EXAMPLES::

            sage: from flatsurf import *
            sage: s0=translation_surfaces.octagon_and_squares()
            sage: a=s0.base_ring().gens()[0]
            sage: m=Matrix([[1,2+a],[0,1]])
            sage: s=m*s0
            sage: s=s.triangulate()
            sage: ss=s.delaunay_decomposition(triangulated=True)
            sage: ss.polygon(0)
            Polygon: (0, 0), (0, -2), (a, -a - 2), (a + 2, -a - 2), (2*a + 2, -2), (2*a + 2, 0), (a + 2, a), (a, a)
            sage: ss.polygon(1)
            Polygon: (0, 0), (0, -2), (2, -2), (2, 0)
            sage: ss.polygon(2)
            Polygon: (0, 0), (-a, a), (-2*a, 0), (-a, -a)
        """
        if not self.is_finite():
            raise NotImplementedError("Not implemented for infinite surfaces.")
        if in_place:
            s=self
        else:
            s=self.mutable_copy()
        if not delaunay_triangulated:
            s.delaunay_triangulation(triangulated=triangulated,in_place=True)
        # Now s is the Delaunay Triangulated
        loop=True
        while loop:
            loop=False
            for (l1,e1),(l2,e2) in s.edge_iterator(gluings=True):
                if (l1<l2 or (l1==l2 and e1<=e2)) and s._edge_needs_join(l1,e1):
                    s.join_polygons(l1, e1, in_place=True)
                    loop=True
                    break
        return s
    
    def graphical_surface(self, *args, **kwds):
        r"""
        Return a GraphicalSurface representing this surface.
        
        By default this returns a cached version of the GraphicalSurface. If
        ``cached=False'' is provided as a keyword option then a new 
        GraphicalSurface is returned. Other keyword options:

        INPUT:

        - ``polygon_labels`` -- a boolean (default ``True``) whether the label
          of polygons are displayed

        - ``edge_labels`` -- option to control the display of edge labels. It
          can be one of

            - ``False`` or ``None`` for no labels

            - ``'gluings'`` -- to put on each side of each non-adjacent edge, the
              name of the polygon to which it is glued

            - ``'number'`` -- to put on each side of each edge the number of the
              edge

            - ``'gluings and numbers'`` -- full information

        EXAMPLES::

            sage: # Test the difference between the cached graphical_surface and the uncached version.
            sage: from flatsurf import *
            sage: s=translation_surfaces.octagon_and_squares()
            sage: print(s.plot())
            Graphics object consisting of 32 graphics primitives
            sage: print(s.graphical_surface(cached=False,adjacencies=[]).plot())
            Graphics object consisting of 18 graphics primitives
        """
        from flatsurf.graphical.surface import GraphicalSurface
        if kwds.has_key("cached"):
            if not kwds["cached"]:
                # cached=False: return a new surface.
                kwds.pop("cached",None)
                return GraphicalSurface(self, *args, **kwds)
            kwds.pop("cached",None)
        if hasattr(self, '_gs'):
            self._gs.process_options(*args, **kwds)
        else:
            self._gs = GraphicalSurface(self, *args, **kwds)
        return self._gs

    def plot(self, *args, **kwds):
        return self.graphical_surface(*args, **kwds).plot()

# I'm not sure we want to support this...
#
#    def minimize_monodromy_mapping(self):
#        r"""
#        Return a mapping from this surface to a similarity surface
#        with a minimal monodromy group. 
#        Note that this may be slow for infinite surfaces.
#        
#        EXAMPLES::
#            sage: from flatsurf.geometry.polygon import Polygons
#            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
#            sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
#            sage: square = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
#            sage: gluings = [((0,i),(1+(i%2),i//2)) for i in range(8)]
#            sage: from flatsurf.geometry.surface import surface_from_polygons_and_gluings
#            sage: s=surface_from_polygons_and_gluings([octagon,square,square],gluings)
#            sage: print s
#            Rational cone surface built from 3 polygons
#            sage: m=s.minimize_monodromy_mapping()
#            sage: s2=m.codomain()
#            sage: print s2
#            Translation surface built from 3 polygons
#            sage: v=s.tangent_vector(2,(0,0),(1,0))
#            sage: print m.push_vector_forward(v)
#            SimilaritySurfaceTangentVector in polygon 2 based at (0, 0) with vector (-1/2*sqrt2, -1/2*sqrt2)
#            sage: w=s2.tangent_vector(2,(0,0),(0,-1))
#            sage: print m.pull_vector_back(w)
#            SimilaritySurfaceTangentVector in polygon 2 based at (0, 0) with vector (1/2*sqrt2, 1/2*sqrt2)
#        """
#        lw = self.walker()
#        class MatrixFunction:
#            def __init__(self, lw):
#                self._lw=lw
#                from sage.matrix.constructor import identity_matrix
#                self._d = {lw.surface().base_label():
#                    identity_matrix(lw.surface().base_ring(), n=2)}
#            def __call__(self, label):
#                try:
#                    return self._d[label]
#                except KeyError:
#                    e = self._lw.edge_back(label)
#                    label2,e2 = self._lw.surface().opposite_edge(label,e)
#                    m=self._lw.surface().edge_matrix(label,e) * self(label2)
#                    self._d[label]=m
#                    return m
#        mf = MatrixFunction(lw)
#        from flatsurf.geometry.mappings import (
#            MatrixListDeformedSurfaceMapping,
#            IdentityMapping)
#        mapping = MatrixListDeformedSurfaceMapping(self, mf)
#        surface_type = mapping.codomain().compute_surface_type_from_gluings(limit=100)
#        new_codomain = convert_to_type(mapping.codomain(),surface_type)
#        identity = IdentityMapping(mapping.codomain(), new_codomain)
#        return identity * mapping
#    
#    def minimal_monodromy_surface(self):
#        r"""
#        Return an equivalent similarity surface with minimal monodromy.
#        Note that this may be slow for infinite surfaces.
#        
#        EXAMPLES::
#            sage: from flatsurf.geometry.polygon import Polygons
#            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
#            sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
#            sage: square = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
#            sage: gluings = [((0,i),(1+(i%2),i//2)) for i in range(8)]
#            sage: from flatsurf.geometry.surface import surface_from_polygons_and_gluings
#            sage: s=surface_from_polygons_and_gluings([octagon,square,square],gluings)
#            sage: print s
#            Rational cone surface built from 3 polygons
#            sage: s2=s.minimal_monodromy_surface()
#            sage: print s2
#            Translation surface built from 3 polygons
#        """
#        return self.minimize_monodromy_mapping().codomain()
    
    def __eq__(self, other):
        r"""
        Implements a naive notion of equality where two finite surfaces are equal if:
        - their base rings are equal,
        - their base labels are equal,
        - their polygons are equal and labeled and glued in the same way.
        For infinite surfaces we use reference equality.
        """
        if not self.is_finite():
            return self is other
        if self is other:
            return True
        if not isinstance(other, SimilaritySurface):
            raise TypeError
        if not other.is_finite():
            raise ValueError("Can not compare infinite surfaces.")
        if not self.is_mutable() and not other.is_mutable():
            hash1 = hash(self)
            hash2 = hash(other)
            if hash1 != hash2:
                return False
        if self.base_ring() != other.base_ring():
            return False
        if self.base_label() != other.base_label():
            return False
        if self.num_polygons() != other.num_polygons():
            return False
        for label,polygon in self.label_iterator(polygons=True):
            try:
                polygon2 = other.polygon(label)
            except ValueError:
                return False
            if polygon != polygon2:
                return False
            for edge in xrange(polygon.num_edges()):
                if self.opposite_edge(label,edge) != other.opposite_edge(label,edge):
                    return False
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        r"""
        Hash compatible with equals.
        """
        if self._s.is_mutable():
            raise ValueError("Attempting to hash with mutable underlying surface.")
        h = 17*hash(self.base_ring())+23*hash(self.base_label())
        for pair in self.label_iterator(polygons=True):
            h = h + 7*hash(pair)
        for edgepair in self.edge_iterator(gluings=True):
            h = h + 3*hash(edgepair)
        return h



