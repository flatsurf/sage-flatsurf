r"""
Similarity surfaces.
"""

from sage.misc.cachefunc import cached_method

from sage.structure.sage_object import SageObject

from sage.sets.family import Family

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

from flatsurf.geometry.surface import *

class SimilaritySurface_generic(SageObject):
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

        Also if your type is more specific than Similarity Surface then you should override:

        - surface_type(self): return the appropriate SurfaceType (see flatsurf.geometry.surface)
    """
    _plot_options = {}

    def _check(self):
        r"""
        Run all the methods that start with _check
        """
        for name in dir(self):
            if name.startswith('_check') and name != '_check':
                print name, "...",
                getattr(self, name)()
                print "done"

    def _check_gluings(self):
        # iterate over pairs with pair1 glued to pair2
        for pair1,pair2 in self.edge_gluing_iterator():
            if not self.opposite_edge_pair(pair2)==pair1:
                raise ValueError("edges not glued correctly:\n%s -> %s -> %s"%(pair1,pair2,self.opposite_edge_pair(pair2)))

    def _check_type(self):
        claimed_type = self.surface_type()
        computed_type = self.compute_surface_type_from_gluings(limit=100)
        if claimed_type != combine_surface_types(claimed_type,computed_type):
            raise ValueError("Surface computed to be of type %s which is not more specific than the claimed type of %s."%(
                surface_type_to_str(computed_type), surface_type_to_str(claimed_type)))
    
    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        raise NotImplementedError

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.

        This method must be overriden in subclasses.
        """
        raise NotImplementedError

    def base_label(self):
        r"""
        Always returns the same label.
        """
        raise NotImplementedError

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.

        This method must be overriden in subclasses.
        """
        raise NotImplementedError

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        raise NotImplementedError

    def surface_type(self):
        r"""
        Return an integer representing the surface's type. The convention is that a surface has the given type 
        if and only if all its 2x2 edge gluing matrixes lie in a certain matrix group.
        """
        return SurfaceType.SIMILARITY

    # 
    # generic methods
    #
    
    def compute_surface_type_from_gluings(self,limit=None):
        r"""
        Compute the surface type by looking at the edge gluings. 
        If limit is defined, we try to guess the type by looking at limit many edges.
        """
        if limit is None:
            if not self.is_finite():
                raise ValueError("Need a limit when working with an infinite surface.")
            it = self.edge_iterator()
            label,edge = it.next()
            # Use honest matrices!
            m = SimilaritySurface_generic.edge_matrix(self,label,edge)
            surface_type = surface_type_from_matrix(m)
            for label,edge in it:
                # Use honest matrices!
                m = SimilaritySurface_generic.edge_matrix(self,label,edge)
                surface_type = combine_surface_types(surface_type, surface_type_from_matrix(m))
            return surface_type
        else:
            count=0
            it = self.edge_iterator()
            label,edge = it.next()
            # Use honest matrices!
            m = SimilaritySurface_generic.edge_matrix(self,label,edge)
            surface_type = surface_type_from_matrix(m)
            for label,edge in it:
                # Use honest matrices!
                m = SimilaritySurface_generic.edge_matrix(self,label,edge)
                surface_type = combine_surface_types(surface_type, surface_type_from_matrix(m))
                count=count+1
                if count >= limit:
                    return surface_type
            return surface_type

    def is_half_dilation_surface(self):
        r"""Return if all the 2x2 gluing matrices are diagonal matrices."""
        return is_half_dilation_surface_type(self.surface_type())

    def is_dilation_surface(self):
        r"""Return if all the 2x2 gluing matrices are positive diagonal matrices."""
        return is_dilation_surface_type(self.surface_type())

    def is_cone_surface(self):
        r"""Return if all the 2x2 gluing matrices lie in O(2)."""
        return is_cone_surface_type(self.surface_type())

    def is_rational_cone_surface(self):
        r"""Return if all the 2x2 gluing matrices are finite order elements of O(2)."""
        return is_rational_cone_surface_type(self.surface_type())

    def is_half_translation_surface(self):
        r"""Return if all the 2x2 gluing matrices are in {I, -I}."""
        return is_half_translation_surface_type(self.surface_type())

    def is_translation_surface(self):
        r"""Return if all the 2x2 gluing matrices are in {I, -I}."""
        return is_translation_surface_type(self.surface_type())

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
            sage: from flatsurf.geometry.translation_surface import TranslationSurface_polygons_and_gluings
            sage: s=TranslationSurface_polygons_and_gluings([tri0,tri1], gluings)
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

        st = surface_type_to_str(self.surface_type())

        if self.num_polygons() == 1:
            end = ""
        else:
            end = "s"

        return "{} built from {} polygon{}".format(st, num, end)

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
            Translation surface built from infinitely many polygons
            sage: T.polygon(T.base_label())
            Polygon: (0, 0), (2, -2), (2, 0)
        """
        from flatsurf.geometry.translation_surface import MinimalTranslationCover
        return MinimalTranslationCover(self)

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

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.tangent_vector(0, (1/2,1/2), (1,1))
            SimilaritySurfaceTangentVector in polygon 1 based at (-1/2, 3/2) with vector
            (-1, -1)
            sage: K.<sqrt2> = QuadraticField(2)
            sage: S.tangent_vector(0, (1/2,1/2), (1,sqrt2))
            SimilaritySurfaceTangentVector in polygon 1 based at (-1/2, 3/2) with vector
            (-1, -sqrt2)

            sage: S = ChamanaraSurface(sqrt2/2)
            sage: S.tangent_vector(1, (0,1), (1,1))
            SimilaritySurfaceTangentVector in polygon 0 based at (-sqrt2, sqrt2
            + 1) with vector (-1, -1)
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
        
    def graphical_surface(self, *args, **kwds):
        r"""
        Return a GraphicalSurface representing this surface.
        """
        opt = self._plot_options.copy()
        opt.update(kwds)
        from flatsurf.graphical.surface import GraphicalSurface
        return GraphicalSurface(self, *args, **opt)

    surface_plot = graphical_surface

    def plot(self, *args, **kwds):
        return self.surface_plot(*args, **kwds).plot()

    def minimize_monodromy_mapping(self):
        r"""
        Return a mapping from this surface to a similarity surface
        with a minimal monodromy group. 
        Note that this may be slow for infinite surfaces.
        
        EXAMPLES::
            sage: from flatsurf.geometry.polygon import Polygons
            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
            sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
            sage: square = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
            sage: gluings = [((0,i),(1+(i%2),i//2)) for i in range(8)]
            sage: from flatsurf.geometry.surface import surface_from_polygons_and_gluings
            sage: s=surface_from_polygons_and_gluings([octagon,square,square],gluings)
            sage: print s
            Rational cone surface built from 3 polygons
            sage: m=s.minimize_monodromy_mapping()
            sage: s2=m.codomain()
            sage: print s2
            Translation surface built from 3 polygons
            sage: v=s.tangent_vector(2,(0,0),(1,0))
            sage: print m.push_vector_forward(v)
            SimilaritySurfaceTangentVector in polygon 2 based at (0, 0) with vector (-1/2*sqrt2, -1/2*sqrt2)
            sage: w=s2.tangent_vector(2,(0,0),(0,-1))
            sage: print m.pull_vector_back(w)
            SimilaritySurfaceTangentVector in polygon 2 based at (0, 0) with vector (1/2*sqrt2, 1/2*sqrt2)
        """
        lw = self.label_walker()
        class MatrixFunction:
            def __init__(self, lw):
                self._lw=lw
                from sage.matrix.constructor import identity_matrix
                self._d = {lw.surface().base_label():
                    identity_matrix(lw.surface().base_ring(), n=2)}
            def __call__(self, label):
                try:
                    return self._d[label]
                except KeyError:
                    e = self._lw.edge_back(label)
                    label2,e2 = self._lw.surface().opposite_edge(label,e)
                    m=self._lw.surface().edge_matrix(label,e) * self(label2)
                    self._d[label]=m
                    return m
        mf = MatrixFunction(lw)
        from flatsurf.geometry.mappings import (
            MatrixListDeformedSurfaceMapping,
            IdentityMapping)
        mapping = MatrixListDeformedSurfaceMapping(self, mf)
        surface_type = mapping.codomain().compute_surface_type_from_gluings(limit=100)
        new_codomain = convert_to_type(mapping.codomain(),surface_type)
        identity = IdentityMapping(mapping.codomain(), new_codomain)
        return identity * mapping
    
    def minimal_monodromy_surface(self):
        r"""
        Return an equivalent similarity surface with minimal monodromy.
        Note that this may be slow for infinite surfaces.
        
        EXAMPLES::
            sage: from flatsurf.geometry.polygon import Polygons
            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
            sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
            sage: square = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
            sage: gluings = [((0,i),(1+(i%2),i//2)) for i in range(8)]
            sage: from flatsurf.geometry.surface import surface_from_polygons_and_gluings
            sage: s=surface_from_polygons_and_gluings([octagon,square,square],gluings)
            sage: print s
            Rational cone surface built from 3 polygons
            sage: s2=s.minimal_monodromy_surface()
            sage: print s2
            Translation surface built from 3 polygons
        """
        return self.minimize_monodromy_mapping().codomain()
        
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
        if not isinstance(other, SimilaritySurface_generic):
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

class SimilaritySurface_polygons_and_gluings(SimilaritySurface_generic):
    r"""
    Similarity surface build from a list of polygons and gluings.
    """
    def __init__(self, *args):
        r"""
        The constructor either acts as a copy constructor for a finite surface, or you can pass two options:
        polygons and identifications.
        
        INPUT:

        - ``polygons`` - a family of polygons (might be a list, a dictionnary
          label -> polygon or more generally a Family)

        - ``identifications`` - the identification of the edges. A list of pairs
          ((p0,e0),(p1,e1)) or

        EXAMPLES::
            sage: from flatsurf.geometry.polygon import Polygons
            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
            sage: octagon = Polygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
            sage: square1 = Polygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
            sage: square2 = Polygons(K)([(sqrt2/2, sqrt2/2),(-sqrt2/2, sqrt2/2),(-sqrt2/2, -sqrt2/2),(sqrt2/2, -sqrt2/2)])
            sage: gluings=[(('b',i),('a', (2*i+4)%8 )) for i in range(4)]
            sage: for i in range(4):
            ...       gluings.append( (('c',i), ('a', (2*i+5)%8 )) )
            sage: from flatsurf.geometry.similarity_surface import SimilaritySurface_polygons_and_gluings
            sage: s=SimilaritySurface_polygons_and_gluings({'a':octagon,'b':square1,'c':square2}, gluings)
            sage: s2=SimilaritySurface_polygons_and_gluings(s)
            sage: s2==s
            True
            sage: hash(s2)==hash(s)
            True
        """
        if len(args)==2:
            polygons=args[0]
            identifications=args[1]

            self._polygons = Family(polygons)

            if self._polygons.cardinality() == 0:
                raise ValueError("there should be at least one polygon")

            self._field = self._polygons.an_element().parent().field()

            n = 0
            for p in self._polygons:
                if p.parent().field() != self._field:
                    raise ValueError("the field must be the same for all polygons")
                n += 1
                if n > 10:  # the number of polygons may be infinite...
                    break

            if isinstance(identifications, (list,dict)):
                edge_identifications = {}
                if isinstance(identifications, dict):
                    it = identifications.iteritems()
                else:
                    it = iter(identifications)
                for e0,e1 in it:
                    edge_identifications[e0] = e1
                    # Check that e0 makes sense. 
                    assert e0[1]>=0 and e0[1]<self._polygons[e0[0]].num_edges()
                    # Check that e1 makes sense. 
                    assert e1[1]>=0 and e1[1]<self._polygons[e1[0]].num_edges()
                    
                    if e1 in edge_identifications:
                        assert edge_identifications[e1] == e0
                    else:
                        edge_identifications[e1] = e0
            else:
                edge_identifications = identifications

            self._edge_identifications = edge_identifications
        elif len(args)==1:
            # Copy constructor for finite surface.
            s = args[0]
            if not s.is_finite():
                raise ValueError("Can only copy finite surface.")
            polygons = {}
            edge_identifications = {}
            for label,polygon in s.label_polygon_iterator():
                polygons[label]=polygon
                for edge in xrange(polygon.num_edges()):
                    edge_identifications[(label,edge)]=s.opposite_edge(label,edge)
            self._field = s.base_ring()
            self._polygons=Family(polygons)
            self._edge_identifications = edge_identifications
        else:
            raise ValueError("Can only be called with one or two arguments.")
    
    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        from sage.rings.infinity import Infinity
        return self._polygons.cardinality() != Infinity

    def num_polygons(self):
        return self._polygons.cardinality()

    def base_ring(self):
        return self._field

    def base_label(self):
        return self.polygon_labels().an_element()

    def polygon_labels(self):
        from sage.combinat.words.alphabet import build_alphabet
        return build_alphabet(self._polygons.keys())

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        return self._polygons[lab]

    def opposite_edge(self, p, e):
        if (p,e) not in self._edge_identifications:
            e = e % self._polygons[p].num_edges()
            if (p,e) not in self._edge_identifications:
                raise ValueError("The pair"+str((p,e))+" is not a valid edge identifier.")
        return self._edge_identifications[(p,e)]
    
class SimilaritySurface_wrapper(SimilaritySurface_generic):
    r"""
    Wraps a surface in a SimilaritySurface package. This is primarily for changing the type of infinite surfaces.
    
    EXAMPLES::
        sage: from flatsurf.geometry.similarity_surface_generators import InfiniteStaircase
        sage: s=InfiniteStaircase()
        sage: from flatsurf.geometry.similarity_surface import convert_to_similarity_surface
        sage: s2=convert_to_similarity_surface(s)
        sage: print s2
        Similarity surface built from infinitely many polygons
    """
    def __init__(self, surface):
        if isinstance(surface, SimilaritySurface_wrapper):
            # It's a party. We'll share data!
            self._s = surface._s
            self._base_ring = surface._base_ring
            self._base_label = surface._base_label
            self._polygons= surface._polygons
            self._edge_identifications = surface._edge_identifications
            self._is_finite = s._is_finite()
        else:
            self._s = surface
            self._base_ring = surface.base_ring()
            self._base_label = surface.base_label()
            self._polygons={}
            self._edge_identifications={}
            self._is_finite = surface.is_finite()

    def base_ring(self):
        return self._base_ring

    def polygon(self, lab):
        try:
            return self._polygons[lab]
        except KeyError:
            p = self._s.polygon(lab)
            self._polygons[lab]=p
            return p

    def base_label(self):
        return self._base_label

    def opposite_edge(self, p, e):
        try:
            return self._edge_identifications[(p,e)]
        except KeyError:
            pair = self._s.opposite_edge(p,e)
            self._edge_identifications[(p,e)]=pair
            self._edge_identifications[pair]=(p,e)
            return pair

    def is_finite(self):
        return self._is_finite

def convert_to_similarity_surface(surface):
    r"""
    Returns a similarity surface version of the provided surface.
    """
    if surface.is_finite():
        return SimilaritySurface_polygons_and_gluings(surface)
    else:
        return SimilaritySurface_wrapper(surface)

