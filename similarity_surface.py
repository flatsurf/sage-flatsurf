r"""
Translation surfaces.
"""

from sage.misc.cachefunc import cached_method

from sage.structure.sage_object import SageObject

from sage.sets.family import Family

from sage.rings.integer import Integer
from sage.rings.rational import Rational

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.qqbar import AA
from sage.rings.real_mpfr import RR
from sage.rings.real_mpfi import RIF

from sage.matrix.constructor import matrix, identity_matrix

ZZ_1 = Integer(1)
ZZ_2 = Integer(2)

from polygon import Polygon,Polygons
from matrix_2x2 import (is_similarity,
                    homothety_rotation_decomposition,
                    similarity_from_vectors,
                    rotation_matrix_angle)

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
        inherit from it and implement the methods:

        - num_polygons(self): the number of polygons in that surface
        - base_ring(self): the base ring in which coordinates lives
        - polygon(self, lab): the polygon associated to the label ``lab``

    It might also be good to implement

        - base_polygon(self): a couple (``label``, ``polygon``) that is a
          somewhat fixed polygon
        - opposite_edge(self, lab, edege): a couple (``other_lab``, ``other_edge``)
    """
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
        for lab in self.polygon_labels().some_elements():
            p = self.polygon(lab)
            for e in xrange(p.num_edges()):
                llab,ee = self.opposite_edge(lab,e)
                lllab,eee = self.opposite_edge(llab,ee)
                if (lllab,eee) != (lab,e):
                    raise ValueError("edges not glued correctly:\n(%s,%s) -> (%s,%s) -> (%s,%s)"%(lab,e,llab,ee,lllab,eee))

    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        raise NotImplementedError

    def polygon_labels(self):
        r"""
        The set of labels used by the polygons.

        This method must be overriden in subclasses.
        """
        raise NotImplementedError

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.

        This method must be overriden in subclasses.
        """
        raise NotImplementedError

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.

        This method must be overriden in subclasses.
        """
        raise NotImplementedError


    #
    # generic methods
    #

    def num_polygons(self):
        r"""
        Return the number of polygons.
        """
        return self.polygon_labels().cardinality()

    def is_finite(self):
        from sage.rings.infinity import Infinity
        return self.num_polygons() != Infinity

    def polygon_iterator(self):
        r"""
        Iterator over the polygons.
        """
        from itertools import imap
        return imap(self.polygon, self.polygon_labels())

    def num_edges(self):
        r"""
        Return the number of edges.
        """
        if self.polygon_labels().is_finite():
            return sum(p.num_edges() for p in self.polygon_iterator())
        else:
            from sage.rings.infinity import Infinity
            return Infinity

    def base_label(self):
        return self.polygon_labels().an_element()

    def _repr_(self):
        if self.num_polygons() == 1:
            end = ""
        else:
            end = "s"
        return "Similarity surface built from %s polygon"%self._polygons.cardinality() + end

    def area(self):
        if self.num_polygons.is_finite():
            return sum(p.area() for p in self.polygon_iterator())
        raise NotImplementedError("area is not implemented for surfaces built from an infinite number of polygons")

    def edge_matrix(self, p, e=None):
        r"""
        Return the edge to which this edge is identified and the matrix to be
        applied.

        EXAMPLES::

            sage: from
            sage: SimilaritySurface([square()]
        """
        if e is None:
            # What does the following line do?
            # -Pat
            p,e = p
        u = self.polygon(p).edge(e)
        pp,ee = self.opposite_edge(p,e)
        v = self.polygon(pp).edge(ee)

        # be careful, because of the orientation, it is -v and not v
        res = similarity_from_vectors(u,-v)
        return similarity_from_vectors(u,-v)

    def edge_dict(self):
        if not self.is_finite():
            raise ValueError("the surface must be finite")
        edges = {}
        for p in self.polygon_labels():
            for e in xrange(self.polygon(p).num_edges()):
                pp,ee = self.opposite_edge(p,e)
                edges[(p,e)] = (pp,ee)

    def minimal_translation_cover(self):
        return MinimalTranslationCover(self)

    def get_bundle(self):
        r"""
        Return a pair (sm,sb), where sm is the active SurfaceManipulator, and sb is the surface
        bundle for this surface (which is added if neccessary to the SurfaceManipulator).
        If necessary, we create one or both objects.
        """
        from surface_manipulator import SurfaceManipulator
        sm = SurfaceManipulator.launch()
        sb = sm.find_bundle(self)
        if sb is None:
            from similarity_surface_bundle import SimilaritySurfaceBundle
            sb = SimilaritySurfaceBundle(self, editor=sm)
            sm.add_surface(sb)
        return sm, sb

    def edit(self):
        r"""
        Launch the tk editor to interactively modify ``self``.
        """
        sm,sb = self.get_bundle()
        sm.set_surface(sb)
        #old version below
        #from translation_surface_editor import TranslationSurfaceEditor
        #fse = TranslationSurfaceEditor(self)
        #fse.window.mainloop()

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), 2)

    def edge_iterator(self):
        for p in self.polygon_labels():
            for e in xrange(self.polygon(p).num_edges()):
                yield p,e

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


class SimilaritySurface_polygons_and_gluings(SimilaritySurface_generic):
    r"""
    Similarity surface build from a list of polygons and gluings.
    """
    def __init__(self, polygons=None, identifications=None):
        r"""
        INPUT:

        - ``polygons`` - a family of polygons (might be a list, a dictionnary
          label -> polygon or more generally a Family)

        - ``identifications`` - the identification of the edges. A list of pairs
          ((p0,e0),(p1,e1)) or

        EXAMPLES::

            sage: from polygon import square
        """
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
                if e1 in edge_identifications:
                    assert edge_identifications[e1] == e0
                else:
                    edge_identifications[e1] = e0
        else:
            edge_identifications = identifications

        self._edge_identifications = edge_identifications

    def num_polygons(self):
        return self._polygons.cardinality()

    def base_ring(self):
        return self._field

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
            raise ValueError("not a valid edge identifier")
        return self._edge_identifications[(p,e)]

class ConicSurface(SimilaritySurface_generic):
    r"""
    A conic surface.
    """
    def angles(self):
        r"""
        Return the set of angles around the vertices of the surface.
        """
        if not self.is_finite():
            raise NotImplementedError("the set of edges is infinite!")

        edges = self.edges()
        edges = set(self.edges())
        angles = []
        while edges:
            p,e = edges.pop()
            angle = self.polygon(p).angle(e)
            pp,ee = self.opposite_edge(p,(e-1)%self.polygon(p).num_edges())
            while pp != p or ee != e:
                edges.remove((pp,ee))
                angle += self.polygon(pp).angle(ee)
                pp,ee = self.opposite_edge(pp,(ee-1)%self.polygon(pp).num_edges())
            angles.append(angle)
        return angles

class TranslationSurface_generic(ConicSurface):
    r"""
    A surface with a flat metric and conical singularities (not necessarily
    multiple angle of pi or 2pi).

    - polygon = polygon + vertex (or equivalently, canonical ordering)

    A translation surface is:

    - field embedded in R
    - index set for the (convex) polygons + favorite polygon
    - edges: ((t1,e1),(t2,e2))

    For finite case:

    - canonical labelings of polygons
    - Delaunay triangulation
    """
    def _check_edge_matrix(self):
        r"""
        Check the compatibility condition
        """
        for lab in self.polygon_labels().some_elements():
            p = self.polygon(lab)
            for e in xrange(p.num_edges()):
                if not self.edge_matrix(lab,e).is_one():
                    raise ValueError("gluings of (%s,%s) is not through translation"%(lab,e))

    def _repr_(self):
        if self.num_polygons() == 1:
            end = ""
        else:
            end = "s"
        return "Translation surface built from %s polygon"%(self.num_polygons()) + end

    def plot(self):
        return TranslationSurfacePlot(self).plot()

    def edge_matrix(self, p, e=None):
        if e is None:
            p,e = p
        if p not in self.polygon_labels():
            raise ValueError
        elif e < 0 or e >= self.polygon(p).num_edges():
            raise ValueError
        return identity_matrix(self.base_ring(),2)

    def stratum(self):
        from sage.dynamics.flat_surfaces.all import AbelianStratum
        return AbelianStratum([a-1 for a in self.angles()])

class TranslationSurface_polygons_and_gluings(
        TranslationSurface_generic,
        SimilaritySurface_polygons_and_gluings):
    pass

class MinimalTranslationCover(TranslationSurface_generic):
    r"""
    We label copy by cartesian product (polygon from bot, matrix).
    """
    def __init__(self, similarity_surface):
        self._ss = similarity_surface

        from sage.matrix.matrix_space import MatrixSpace
        from sage.categories.cartesian_product import cartesian_product
        from sage.rings.semirings.non_negative_integer_semiring import NN

    def base_ring(self):
        return self._ss.base_ring()

    def base_label(self):
        from sage.matrix.constructor import identity_matrix
        I = identity_matrix(self.base_ring(),2)
        I.set_immutable()
        return (self._ss.base_label(), I)

    def polygon(self, lab):
        return lab[1] * self._ss.polygon(lab[0])

    def polygon_labels(self):
        r"""
        Return the set of polygons used for the labels.
        """
        from cartesian_product import CartesianProduct
        from finitely_generated_matrix_group import FinitelyGeneratedMatrixSubgroup

        ss = self._ss

        M = [ss.edge_matrix(p,e) for p,e in ss.edge_iterator()]
        for m in M: m.set_immutable()
        M = sorted(set(M))
        for m in M: m.set_immutable()
        G = FinitelyGeneratedMatrixSubgroup(M)
        return CartesianProduct([ss.polygon_labels(), G])

    def opposite_edge(self, p, e):
        pp,m = p  # this is the polygon m * ss.polygon(p)
        p2,e2 = self._ss.opposite_edge(pp,e)
        me = self._ss.edge_matrix(pp,e)
        mm = ~me * m
        mm.set_immutable()
        return ((p2,mm),e2)

class Origami(TranslationSurface_generic):
    def __init__(self, r, u, rr=None, uu=None, domain=None):
        if domain is None:
            self._domain = r.parent().domain()
        else:
            self._domain = domain

        self._r = r
        self._u = u
        if rr is None:
            rr = ~r
        else:
            for a in self._domain.some_elements():
                if r(rr(a)) != a:
                    raise ValueError("r o rr is not identity on %s"%a)
                if rr(r(a)) != a:
                    raise ValueError("rr o r is not identity on %s"%a)
        if uu is None:
            uu = ~u
        else:
            for a in self._domain.some_elements():
                if u(uu(a)) != a:
                    raise ValueError("u o uu is not identity on %s"%a)
                if uu(u(a)) != a:
                    raise ValueError("uu o u is not identity on %s"%a)

        self._perms = [uu,r,u,rr] # down,right,up,left

    def num_polygons(self):
        r"""
        Returns the number of polygons.
        """
        return self._domain.cardinality()

    def polygon_labels(self):
        return self._domain

    def polygon(self, lab):
        if lab not in self._domain:
            raise ValueError
        from polygon import square
        return square()

    def base_ring(self):
        return QQ

    def _repr_(self):
        return "Origami defined by r=%s and u=%s"%(self._r,self._u)

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e < 0 or e > 3:
            raise ValueError
        return self._perms[e](p), (e+2)%4

