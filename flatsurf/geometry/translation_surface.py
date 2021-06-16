r"""
Translation Surfaces.
"""

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

from sage.matrix.constructor import identity_matrix

from .surface import Surface
from .half_translation_surface import HalfTranslationSurface
from .dilation_surface import DilationSurface

class TranslationSurface(HalfTranslationSurface, DilationSurface):
    r"""
    A surface with a flat metric and conical singularities whose cone angles are a multiple of pi.
    """

    def minimal_translation_cover(self):
        return self

    def _test_edge_matrix(self, **options):
        r"""
        Check the compatibility condition
        """
        tester = self._tester(**options)

        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if self.is_finite():
            it = self.label_iterator()
        else:
            from itertools import islice
            it = islice(self.label_iterator(), 30)

        for lab in it:
            p = self.polygon(lab)
            for e in range(p.num_edges()):
                # Warning: check the matrices computed from the edges,
                # rather the ones overriden by TranslationSurface.
                m=SimilaritySurface.edge_matrix(self,lab,e)
                tester.assertTrue(m.is_one(), \
                    "edge_matrix of edge "+str((lab,e))+" is not a translation.")

    def edge_matrix(self, p, e=None):
        if e is None:
            p,e = p
        if e < 0 or e >= self.polygon(p).num_edges():
            raise ValueError
        return identity_matrix(self.base_ring(),2)

    def stratum(self):
        r"""
        Return the stratum this surface belongs to.

        This uses the package ``surface_dynamics``
        (see http://www.labri.fr/perso/vdelecro/flatsurf_sage.html)

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: sfg.translation_surfaces.octagon_and_squares().stratum()
            H_3(4)
        """
        from surface_dynamics import AbelianStratum
        from sage.rings.integer_ring import ZZ
        return AbelianStratum([ZZ(a-1) for a in self.angles()])

    def _canonical_first_vertex(polygon):
        r"""
        Return the index of the vertex with smallest y-coordinate.
        If two vertices have the same y-coordinate, then the one with least x-coordinate is returned.
        """
        best=0
        best_pt=polygon.vertex(best)
        for v in range(1,polygon.num_edges()):
            pt=polygon.vertex(v)
            if pt[1]<best_pt[1]:
                best=v
                best_pt=pt
        if best==0:
            if pt[1]==best_pt[1]:
                return v
        return best

    def standardize_polygons(self, in_place=False):
        r"""
        Replaces each polygon with a polygon with a new polygon which differs by translation
        and reindexing. The new polygon will have the property that vertex zero is the origin,
        and all vertices lie either in the upper half plane, or on the x-axis with non-negative
        x-coordinate.

        This is done to the current surface if in_place=True. A mutable copy is created and returned
        if in_place=False (as default).

        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.veech_double_n_gon(4)
            sage: s.polygon(1)
            Polygon: (0, 0), (-1, 0), (-1, -1), (0, -1)
            sage: [s.opposite_edge(0,i) for i in range(4)]
            [(1, 0), (1, 1), (1, 2), (1, 3)]
            sage: ss=s.standardize_polygons()
            sage: ss.polygon(1)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: [ss.opposite_edge(0,i) for i in range(4)]
            [(1, 2), (1, 3), (1, 0), (1, 1)]
            sage: TestSuite(ss).run()

        Make sure first vertex is sent to origin::
            sage: from flatsurf import *
            sage: P = ConvexPolygons(QQ)
            sage: p = P(vertices = ([(1,1),(2,1),(2,2),(1,2)]))
            sage: s = Surface_list(QQ)
            sage: s.add_polygon(p)
            0
            sage: s.change_polygon_gluings(0, [(0,2),(0,3),(0,0),(0,1)])
            sage: s.change_base_label(0)
            sage: ts = TranslationSurface(s)
            sage: ts.standardize_polygons().polygon(0)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
        """
        if self.is_finite():
            if in_place:
                if not self.is_mutable():
                    raise ValueError("An in_place call for standardize_polygons can only be done for a mutable surface.")
                s=self
            else:
                s=self.copy(mutable=True)
            cv = {} # dictionary for non-zero canonical vertices
            translations={} # translations bringing the canonical vertex to the origin.
            for l,polygon in s.label_iterator(polygons=True):
                best=0
                best_pt=polygon.vertex(best)
                for v in range(1,polygon.num_edges()):
                    pt=polygon.vertex(v)
                    if (pt[1]<best_pt[1]) or (pt[1]==best_pt[1] and pt[0]<best_pt[0]):
                        best=v
                        best_pt=pt
                # We replace the polygon if the best vertex is not the zero vertex, or
                # if the coordinates of the best vertex differs from the origin.
                if not (best==0 and best_pt.is_zero()):
                    cv[l]=best
            for l,v in iteritems(cv):
                s.set_vertex_zero(l,v,in_place=True)
            return s
        else:
            assert in_place == False, "In place standardization only available for finite surfaces."
            return TranslationSurface(LazyStandardizedPolygonSurface(self))

    def cmp(self, s2, limit=None):
        r"""
        Compare two surfaces. This is an ordering returning -1, 0, or 1.

        The surfaces will be considered equal if and only if there is a translation automorphism
        respecting the polygons and the base_labels.

        If the two surfaces are infinite, we just examine the first limit polygons.
        """
        if self.is_finite():
            if s2.is_finite():
                assert limit is None, "Limit only enabled for finite surfaces."

                #print("comparing number of polygons")
                sign = self.num_polygons()-s2.num_polygons()
                if sign>0:
                    return 1
                if sign<0:
                    return -1
                #print("comparing polygons")
                lw1=self.walker()
                lw2=s2.walker()
                for p1,p2 in zip(lw1.polygon_iterator(), lw2.polygon_iterator()):
                    # Uses Polygon.cmp:
                    ret = p1.cmp(p2)
                    if ret != 0:
                        return ret
                # Polygons are identical. Compare edge gluings.
                #print("comparing edge gluings")
                for pair1,pair2 in zip(lw1.edge_iterator(), lw2.edge_iterator()):
                    l1,e1 = self.opposite_edge(pair1)
                    l2,e2 = s2.opposite_edge(pair2)
                    num1 = lw1.label_to_number(l1)
                    num2 = lw2.label_to_number(l2)
                    ret = (num1 > num2) - (num1 < num2)
                    if ret:
                        return ret
                    ret = (e1 > e2) - (e1 < e2)
                    if ret:
                        return ret
                return 0
            else:
                # s1 is finite but s2 is infinite.
                return -1
        else:
            if s2.is_finite():
                # s1 is infinite but s2 is finite.
                return 1
            else:
                # both surfaces are infinite.
                lw1=self.walker()
                lw2=s2.walker()
                count = 0
                for (l1,p1),(l2,p2) in zip(lw1.label_polygon_iterator(), lw2.label_polygon_iterator()):
                    # Uses Polygon.cmp:
                    ret = p1.cmp(p2)
                    if ret != 0:
                        print("Polygons differ")
                        return ret
                    # If here the number of edges should be equal.
                    for e in range(p1.num_edges()):
                        ll1,ee1 = self.opposite_edge(l1,e)
                        ll2,ee2 = s2.opposite_edge(l2,e)
                        num1 = lw1.label_to_number(ll1, search=True, limit=limit)
                        num2 = lw2.label_to_number(ll2, search=True, limit=limit)
                        ret = (num1 > num2) - (num1 < num2)
                        if ret:
                            return ret
                        ret = (ee1 > ee2) - (ee1 < ee2)
                        if ret:
                            return ret
                    if count >= limit:
                        break
                    count += 1
                return 0

    # TODO: deprecation
    cmp_translation_surface = cmp

    def canonicalize_mapping(self):
        r"""
        Return a SurfaceMapping canonicalizing this translation surface.
        """
        from flatsurf.geometry.mappings import canonicalize_translation_surface_mapping, IdentityMapping
        return canonicalize_translation_surface_mapping(self)

    def canonicalize(self, in_place=False):
        r"""
        Return a canonical version of this translation surface.

        EXAMPLES:

        We will check if an element lies in the Veech group::

            sage: from flatsurf import *
            sage: s = translation_surfaces.octagon_and_squares()
            sage: s
            TranslationSurface built from 3 polygons
            sage: a = s.base_ring().gen()
            sage: mat = Matrix([[1,2+a],[0,1]])
            sage: s1 = s.canonicalize()
            sage: s1.underlying_surface().set_immutable()
            sage: s2 = (mat*s).canonicalize()
            sage: s2.underlying_surface().set_immutable()
            sage: s1.cmp(s2) == 0
            True
            sage: hash(s1) == hash(s2)
            True
        """
        # Old version
        #return self.canonicalize_mapping().codomain()
        if in_place:
            if not self.is_mutable():
                raise ValueError("canonicalize with in_place=True is only defined for mutable translation surfaces.")
            s=self
        else:
            s=self.copy(mutable=True)
        if not s.is_finite():
            raise ValueError("canonicalize is only defined for finite translation surfaces.")
        ret=s.delaunay_decomposition(in_place=True)
        s.standardize_polygons(in_place=True)
        ss=s.copy(mutable=True)
        labels={label for label in s.label_iterator()}
        labels.remove(s.base_label())
        for label in labels:
            ss.underlying_surface().change_base_label(label)
            if ss.cmp(s)>0:
                s.underlying_surface().change_base_label(label)
        # We now have the base_label correct.
        # We will use the label walker to generate the canonical labeling of polygons.
        w=s.walker()
        w.find_all_labels()
        s.relabel(w.label_dictionary(), in_place=True)
        # Set immutable
        s.underlying_surface().set_immutable()
        return s

    def rel_deformation(self, deformation, local=False, limit=100):
        r"""
        Perform a rel deformation of the surface and return the result.

        This algorithm currently assumes that all polygons affected by this deformation are
        triangles. That should be fixable in the future.

        INPUT:

        - ``deformation`` (dictionary) - A dictionary mapping singularities of
          the surface to deformation vectors (in some 2-dimensional vector
          space). The rel deformation being done will move the singularities
          (relative to each other) linearly to the provided vector for each
          vertex. If a singularity is not included in the dictionary then the
          vector will be treated as zero.

        - ``local`` - (boolean) - If true, the algorithm attempts to deform all
          the triangles making up the surface without destroying any of them.
          So, the area of the triangle must be positive along the full interval
          of time of the deformation.  If false, then the deformation must have
          a particular form: all vectors for the deformation must be paralell.
          In this case we achieve the deformation with the help of the SL(2,R)
          action and Delaunay triangulations.

        - ``limit`` (integer) - Restricts the length of the size of SL(2,R)
          deformations considered. The algorithm should be roughly worst time
          linear in limit.

        TODO:

        - Support arbitrary rel deformations.
        - Remove the requirement that triangles be used.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.arnoux_yoccoz(4)
            sage: field = s.base_ring()
            sage: a = field.gen()
            sage: V = VectorSpace(field,2)
            sage: deformation1 = {s.singularity(0,0):V((1,0))}
            sage: s1 = s.rel_deformation(deformation1).canonicalize()
            sage: deformation2 = {s.singularity(0,0):V((a,0))}
            sage: s2 = s.rel_deformation(deformation2).canonicalize()
            sage: m = Matrix([[a,0],[0,~a]])
            sage: s2.cmp((m*s1).canonicalize())
            0
        """
        s=self
        # Find a common field
        field=s.base_ring()
        for singularity, v in iteritems(deformation):
            if v.parent().base_field() != field:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                field = cm.common_parent(field, v.parent().base_field())
        from sage.modules.free_module import VectorSpace
        vector_space = VectorSpace(field,2)

        from collections import defaultdict
        vertex_deformation=defaultdict(vector_space.zero) # dictionary associating the vertices.
        deformed_labels=set() # list of polygon labels being deformed.

        for singularity, vect in iteritems(deformation):
            # assert s==singularity.surface()
            for label,v in singularity.vertex_set():
                vertex_deformation[(label,v)]=vect
                deformed_labels.add(label)
                assert s.polygon(label).num_edges()==3

        from flatsurf.geometry.polygon import wedge_product, ConvexPolygons

        if local:

            ss=s.copy(mutable=True, new_field=field)
            us=ss.underlying_surface()

            P = ConvexPolygons(field)
            for label in deformed_labels:
                polygon=s.polygon(label)
                a0=vector_space(polygon.vertex(1))
                b0=vector_space(polygon.vertex(2))
                v0=vector_space(vertex_deformation[(label,0)])
                v1=vector_space(vertex_deformation[(label,1)])
                v2=vector_space(vertex_deformation[(label,2)])
                a1=v1-v0
                b1=v2-v0
                # We deform by changing the triangle so that its vertices 1 and 2 have the form
                # a0+t*a1 and b0+t*b1
                # respectively. We are deforming from t=0 to t=1.
                # We worry that the triangle degenerates along the way.
                # The area of the deforming triangle has the form
                # A0 + A1*t + A2*t^2.
                A0 = wedge_product(a0,b0)
                A1 = wedge_product(a0,b1)+wedge_product(a1,b0)
                A2 = wedge_product(a1,b1)
                if A2 != field.zero():
                    # Critical point of area function
                    c = A1/(-2*A2)
                    if field.zero()<c and c<field.one():
                        assert A0+A1*c+A2*c**2 > field.zero(), "Triangle with label %r degenerates at critical point before endpoint" % label
                assert A0+A1+A2 > field.zero(), "Triangle with label %r degenerates at or before endpoint" % label
                # Triangle does not degenerate.
                us.change_polygon(label,P(vertices=[vector_space.zero(),a0+a1,b0+b1]))
            return ss

        else: # Non local deformation
            # We can only do this deformation if all the rel vector are parallel.
            # Check for this.
            nonzero=None
            for singularity, vect in iteritems(deformation):
                vvect=vector_space(vect)
                if vvect!=vector_space.zero():
                    if nonzero is None:
                        nonzero=vvect
                    else:
                        assert wedge_product(nonzero,vvect)==0, \
                            "In non-local deformation all deformation vectos must be parallel"
            assert nonzero is not None, "Deformation appears to be trivial."
            from sage.matrix.constructor import Matrix
            m=Matrix([[nonzero[0],-nonzero[1]],[nonzero[1],nonzero[0]]])
            mi=~m
            g=Matrix([[1,0],[0,2]],ring=field)
            prod=m*g*mi
            ss=None
            k=0
            while True:
                if ss is None:
                    ss=s.copy(mutable=True, new_field=field)
                else:
                    # In place matrix deformation
                    ss.apply_matrix(prod)
                ss.delaunay_triangulation(direction=nonzero,in_place=True)
                deformation2={}
                for singularity, vect in iteritems(deformation):
                    found_start=None
                    for label,v in singularity.vertex_set():
                        if wedge_product(s.polygon(label).edge(v),nonzero) >= 0 and \
                        wedge_product(nonzero,-s.polygon(label).edge((v+2)%3)) > 0:
                            found_start=(label,v)
                            found=None
                            for vv in range(3):
                                if wedge_product(ss.polygon(label).edge(vv),nonzero) >= 0 and \
                                wedge_product(nonzero,-ss.polygon(label).edge((vv+2)%3)) > 0:
                                    found=vv
                                    deformation2[ss.singularity(label,vv)]=vect
                                    break
                            assert found is not None
                            break
                    assert found_start is not None
                try:
                    sss=ss.rel_deformation(deformation2,local=True)
                    sss.apply_matrix(mi*g**(-k)*m)
                    sss.delaunay_triangulation(direction=nonzero,in_place=True)
                    return sss
                except AssertionError as e:
                    pass
                k=k+1
                if limit is not None and k>=limit:
                    assert False, "Exeeded limit iterations"

    def j_invariant(self):
        r"""
        Return the Kenyon-Smillie J-invariant of this translation surface.

        It is assumed that the coordinates are defined over a number field.

        EXAMPLES::

            sage: from flatsurf import *
            sage: O = translation_surfaces.regular_octagon()
            sage: O.j_invariant()
            (
                      [2 2]
            (0), (0), [2 1]
            )
        """
        it = self.label_iterator()
        lab = next(it)
        P = self.polygon(lab)
        Jxx, Jyy, Jxy = P.j_invariant()
        for lab in it:
            xx,yy,xy = self.polygon(lab).j_invariant()
            Jxx += xx
            Jyy += yy
            Jxy += xy
        return (Jxx, Jyy, Jxy)

class MinimalTranslationCover(Surface):
    r"""
    Do not use translation_surface.MinimalTranslationCover. Use 
    minimal_cover.MinimalTranslationCover instead. This class is being
    deprecated.
    """
    def __init__(self, similarity_surface):
        if similarity_surface.underlying_surface().is_mutable():
            if similarity_surface.is_finite():
                self._ss=similarity_surface.copy()
            else:
                raise ValueError("Can not construct MinimalTranslationCover of a surface that is mutable and infinite.")
        else:
            self._ss = similarity_surface

        # We are finite if and only if self._ss is a finite RationalConeSurface.
        if not self._ss.is_finite():
            finite = False
        else:
            try:
                from flatsurf.geometry.rational_cone_surface import RationalConeSurface
                ss_copy = self._ss.reposition_polygons(relabel=True)
                rcs = RationalConeSurface(ss_copy)
                rcs._test_edge_matrix()
                finite=True
            except AssertionError:
                # print("Warning: Could be indicating infinite surface falsely.")
                finite=False

        I = identity_matrix(self._ss.base_ring(),2)
        I.set_immutable()
        base_label=(self._ss.base_label(), I)

        Surface.__init__(self, self._ss.base_ring(), base_label, finite=finite, mutable=False)

    def polygon(self, lab):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C.polygon('a')
            Traceback (most recent call last):
            ...
            ValueError: invalid label 'a'
        """
        if not isinstance(lab, tuple) or len(lab) != 2:
            raise ValueError("invalid label {!r}".format(lab))
        p = self._ss.polygon(lab[0])
        return lab[1] * self._ss.polygon(lab[0])

    def opposite_edge(self, p, e):
        pp,m = p  # this is the polygon m * ss.polygon(p)
        p2,e2 = self._ss.opposite_edge(pp,e)
        me = self._ss.edge_matrix(pp,e)
        mm = ~me * m
        mm.set_immutable()
        return ((p2,mm),e2)

class AbstractOrigami(Surface):
    r'''Abstract base class for origamis.
    Realization needs just to define a _domain and four cardinal directions.
    '''
    def __init__(self, domain, base_label=None):
        self._domain=domain
        if base_label is None:
            base_label = domain.an_element()
        from sage.rings.rational_field import QQ
        Surface.__init__(self,QQ,base_label,finite=domain.is_finite(),mutable=False)

    def up(self, label):
        raise NotImplementedError

    def down(self, label):
        raise NotImplementedError

    def right(self, label):
        raise NotImplementedError

    def left(self, label):
        raise NotImplementedError

    def _repr_(self):
        return "Some AbstractOrigami"

    def num_polygons(self):
        r"""
        Returns the number of polygons.
        """
        return self._domain.cardinality()

    def polygon_labels(self):
        return self._domain

    def polygon(self, lab):
        if lab not in self._domain:
            #Updated to print a possibly useful error message
            raise ValueError("Label "+str(lab)+" is not in the domain")
        from flatsurf.geometry.polygon import polygons
        return polygons.square()

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e==0:
            return self.down(p),2
        if e==1:
            return self.right(p),3
        if e==2:
            return self.up(p),0
        if e==3:
            return self.left(p),1
        raise ValueError


class Origami(AbstractOrigami):
    def __init__(self, r, u, rr=None, uu=None, domain=None, base_label=None):
        if domain is None:
            domain = r.parent().domain()

        self._r = r
        self._u = u
        if rr is None:
            rr = ~r
        else:
            for a in domain.some_elements():
                if r(rr(a)) != a:
                    raise ValueError("r o rr is not identity on %s"%a)
                if rr(r(a)) != a:
                    raise ValueError("rr o r is not identity on %s"%a)
        if uu is None:
            uu = ~u
        else:
            for a in domain.some_elements():
                if u(uu(a)) != a:
                    raise ValueError("u o uu is not identity on %s"%a)
                if uu(u(a)) != a:
                    raise ValueError("uu o u is not identity on %s"%a)

        self._perms = [uu,r,u,rr] # down,right,up,left
        AbstractOrigami.__init__(self,domain,base_label)

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError("Polygon label p="+str(p)+" is not in domain="+str(self._domain))
        if e < 0 or e > 3:
            raise ValueError("Edge value e="+str(e)+" does not satisfy 0<=e<4.")
        return self._perms[e](p), (e+2)%4

    def up(self, label):
        return self.opposite_edge(label,2)[0]

    def down(self, label):
        return self.opposite_edge(label,0)[0]

    def right(self, label):
        return self.opposite_edge(label,1)[0]

    def left(self, label):
        return self.opposite_edge(label,3)[0]

    def _repr_(self):
        return "Origami defined by r=%s and u=%s"%(self._r,self._u)

class LazyStandardizedPolygonSurface(Surface):
    r"""
    This class handles standardizing polygons for infinite translation surfaces.
    See the TranslationSurface.standardize_polygons method.

    This class should not be instantiated directly.
    Instead use TranslationSurface.standardize_polygons.
    """

    def __init__(self, surface, relabel=False):
        self._s = surface.copy(mutable=True, relabel=relabel)
        self._labels = set()
        Surface.__init__(self, self._s.base_ring(), self._s.base_label(), finite=self._s.is_finite(), mutable=False)

    def standardize(self, label):
        best=0
        polygon = self._s.polygon(label)
        best_pt=polygon.vertex(best)
        for v in range(1,polygon.num_edges()):
            pt=polygon.vertex(v)
            if (pt[1]<best_pt[1]) or (pt[1]==best_pt[1] and pt[0]<best_pt[0]):
                best=v
                best_pt=pt
        if best!=0:
            self._s.set_vertex_zero(label,best,in_place=True)
        self._labels.add(label)

    def polygon(self, label):
        r"""
        Return the polygon with the provided label.

        This method must be overriden in subclasses.
        """
        if label in self._labels:
            return self._s.polygon(label)
        else:
            self.standardize(label)
            return self._s.polygon(label)

    def opposite_edge(self, l, e):
        r"""
        Given the label ``l`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``ll``, ``ee``) to which this edge is glued.

        This method must be overriden in subclasses.
        """
        if l not in self._labels:
            self.standardize(l)
        ll,ee = self._s.opposite_edge(l,e)
        if ll in self._labels:
            return (ll,ee)
        else:
            self.standardize(ll)
            return self._s.opposite_edge(l,e)
