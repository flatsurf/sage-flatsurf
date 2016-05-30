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

from flatsurf.geometry.matrix_2x2 import (is_similarity,
                    homothety_rotation_decomposition,
                    similarity_from_vectors,
                    rotation_matrix_angle,
                    is_cosine_sine_of_rational)
                    
from flatsurf.geometry.similarity import SimilarityGroup

from flatsurf.geometry.surface import Surface, LabelWalker

class SimilaritySurface(Surface):
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
            it = self.edge_gluing_iterator()
        else:
            from itertools import islice
            it = islice(self.edge_gluing_iterator(), 30)

        for pair1,pair2 in it:
            tester.assertEqual(self.opposite_edge_pair(pair2), pair1,
                "edges not glued correctly:\n%s -> %s -> %s"%(pair1,pair2,self.opposite_edge_pair(pair2)))

    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        return self._s.base_ring()

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.

        This method must be overriden in subclasses.
        """
        return self._s.polygon(lab)

    def base_label(self):
        r"""
        Always returns the same label.
        """
        return self._s.base_label()

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.

        This method must be overriden in subclasses.
        """
        return self._s.opposite_edge(p,e)

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        return self._s.is_finite()

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

    def opposite_edge_pair(self,label_edge_pair):
        r"""
        Same as opposite_edge(label, pair) except taking a pair as input.
        """
        return self.opposite_edge(label_edge_pair[0], label_edge_pair[1])

    def label_walker(self):
        try:
            return self._lw
        except AttributeError:
            self._lw = LabelWalker(self)
        return self._lw

    def label_iterator(self):
        r"""
        Iterator over the polygon labels.
        """
        return iter(self.label_walker())

    def polygon_iterator(self):
        r"""
        Iterate over the polygons.
        """
        return self.label_walker().polygon_iterator()

    def label_polygon_iterator(self):
        r"""
        Iterate over pairs (label,polygon).
        """
        return self.label_walker().label_polygon_iterator()

    def edge_iterator(self):
        r"""
        Iterate over the edges of polygons, which are pairs (l,e) where l is a polygon label, 0 <= e < N and N is the number of edges of the polygon with label l.
        
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
        return self.label_walker().edge_iterator()

    def edge_gluing_iterator(self):
        r"""
        Iterate over the ordered pairs of edges being glued.
        """
        for label,edge in self.edge_iterator():
            label2,edge2 = self.opposite_edge(label,edge)
            yield ((label,edge),(label2,edge2))

    def num_polygons(self):
        r"""
        Return the number of polygons.
        """
        if self.is_finite():
            lw=self.label_walker()
            lw.find_all_labels()
            return len(lw)
        else:
            from sage.rings.infinity import Infinity
            return Infinity

    def num_edges(self):
        r"""
        Return the total number of edges of all polygons used.
        """
        if self.is_finite():
            lw=self.label_walker()
            lw.find_all_labels()
            return sum(p.num_edges() for p in self.polygon_iterator())
        else:
            from sage.rings.infinity import Infinity
            return Infinity

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

    def area(self):
        if self.num_polygons.is_finite():
            return sum(p.area() for p in self.polygon_iterator())
        raise NotImplementedError("area is not implemented for surfaces built from an infinite number of polygons")

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

    def edge_dict(self):
        if not self.is_finite():
            raise ValueError("the surface must be finite")
        edges = {}
        for l,p in self.label_polygon_iterator():
            for e in xrange(p.num_edges()):
                ll,ee = self.opposite_edge(l,e)
                edges[(l,e)] = (ll,ee)

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

    def fundamental_group_basis(self):
        r"""
        Return a basis for the fundamental group as a sequence of paths:

        [vertex0, edge0, vertex1, edge1, ...].
        """
        raise NotImplementedError
        if not self.is_finite():
            raise ValueError("the method would dramatically fails for infinite surfaces!!!")

        tree = {}   # goes from leaves to root self.polygon_labels()
        basis = []

        p = self.base_label() # the root of the tree
        tree[p] = (None,None)

        wait = [] # list of triples p1 -- e --> p2
        for e in xrange(self.polygon(p).num_edges()):
            pp,ee = self.opposite_edge(p,e)
            wait.append((pp,ee,p,e))
        while wait:
            p1,e1,p2,e2 = wait.pop()
            if p1 in tree: # new cycle
                if p1 < p2 or (p1 == p2 and e1 < e2):
                    i = p1
                    p1_to_root = [i]
                    while i != None:
                        i,e = tree[i]
                        p1_to_root.append(e)
                        p1_to_root.append(i)
            else:
                tree[p1] = (p2,e)

        return tree,bridges

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
        Return a SurfaceMapping triangulating the suface.
        """
        from flatsurf.geometry.mappings import triangulation_mapping
        return triangulation_mapping(self)
    
    def triangulate(self):
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
        return self.triangulation_mapping().codomain()
        
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
        if "cached" in kwds:
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
#        lw = self.label_walker()
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
        if self.base_ring() != other.base_ring():
            return False
        if self.base_label() != other.base_label():
            return False
        if self.num_polygons() != other.num_polygons():
            return False
        for label,polygon in self.label_polygon_iterator():
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
        h = 17*hash(self.base_ring())+23*hash(self.base_label())
        for pair in self.label_polygon_iterator():
            h = h + 7*hash(pair)
        for edgepair in self.edge_gluing_iterator():
            h = h + 3*hash(edgepair)
        return h

