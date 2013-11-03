r"""
Translation surface in Sage.
Stupid edit here from Pat.

EXAMPLES:

Right now we can do from sage::

    sage: import translation_surface
    sage: t = TranslationSurface()
    sage: t.edit()

    sage: translation_surfaces.octagon()
    Translation surface built from an octagon

A flat torus::

        sage: S = ConicSurface([square()], [((0,0),(0,2)), ((0,1),(0,3))])
        sage: S
        Similarity surface built from 1 polygon
        sage: S.edge_matrix(0,0)
        [1 0]
        [0 1]
        sage: S.angles()

A surface with conical singularities::

        sage: S = ConicSurface([square()], [((0,0),(0,1)),((0,2),(0,3))])
        sage: S.edge_matrix(0,0)
        [ 0  1]
        [-1  0]
        sage: S.angles()
        [4]
"""

from sage.misc.cachefunc import cached_method

from sage.structure.sage_object import SageObject

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


class SimilaritySurface(SageObject):
    r"""
    An oriented surface built from a set of polygons and edges identified with
    similarities (i.e. composition of homothety, rotations and translations).

    Each polygon is identified with a unique key (its label). The choice of the
    label of the polygons is done at startup. If the set is finite then by
    default the labels are the first non-negative integers 0,1,...

    The edge are identified by a couple (polygon label, edge number).

    Remark:

    If an edge is glued to itself we get a pole in the middle... it seems to be
    not ennoying.
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
        from sage.sets.family import Family
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

    def _repr_(self):
        if self.polygons().cardinality() == 1:
            end = ""
        else:
            end = "s"
        return "Similarity surface built from %s polygon"%self._polygons.cardinality() + end

    def base_ring(self):
        return self._field

    def area(self):
        if self._polygons.is_finite():
            return sum(p.area() for p in self._polygons)
        raise NotImplementedError("area is not implemented for surfaces built from an infinite number of polygons")

    def polygons(self):
        r"""
        The family of polygons that constitute self.
        """
        return self._polygons

    def edges(self):
        r"""
        Return the set of edges.
        """
        return self._edge_identifications.keys()

    def polygon_labels(self):
        r"""
        Return the set of labels used for the polygons.
        """
        return self._polygons.keys()

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        return self._polygons[lab]

    def edge_to_identified(self, p, e):
        if (p,e) not in self._edge_identifications:
            raise ValueError("not a valid edge identifier")
        return self._edge_identifications[(p,e)]

    @cached_method
    def edge_matrix(self, p, e=None):
        r"""
        Return the edge to which this edge is identified and the matrix to be
        applied.

        EXAMPLES::

            sage: from
            sage: SimilaritySurface([square()]
        """
        if e is None:
            p,e = p
        u = self.polygon(p).edge(e)
        pp,ee = self.edge_to_identified(p,e)
        v = self.polygon(pp).edge(ee)

        # be careful, because of the orientation, it is -v and not v
        return similarity_from_vectors(u,-v)

class ConicSurface(SimilaritySurface):
    r"""
    A conic surface.
    """
    def _check(self):
        r"""
        We check that the matrices are isometries.
        """

    def angles(self):
        r"""
        Return the set of angles around the vertices of the surface.
        """
        if not self.polygons().is_finite():
            raise NotImplementedError("the set of edges is infinite!")

        edges = self.edges()
        edges = set(self.edges())
        angles = []
        while edges:
            p,e = edges.pop()
            angle = self.polygon(p).angle(e)
            pp,ee = self.edge_to_identified(p,(e-1)%self.polygon(p).num_edges())
            while pp != p or ee != e:
                edges.remove((pp,ee))
                angle += self.polygon(pp).angle(ee)
                pp,ee = self.edge_to_identified(pp,(ee-1)%self.polygon(pp).num_edges())
            angles.append(angle)
        return angles

class TranslationSurface(ConicSurface):
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
    def _repr_(self):
        return "Translation surface built from %d polygons"%self.num_polygons()

    def to_str(self):
        r"""
        Return a string that can be used to reconstruct self.
        """
        K = self._field
        P = K.defining_polynomial()
        print "var:%s"%P.variable_name()
        print "poly:%s"%P
        print "embedding:%s"%K.gen_embedding()
        print "reliable_comparison:%s"%self._reliable_comparison
        # WHAT ABOUT THE index_set and edge_labels (should be a command)!!!
        print "index_set:XXX"
        print "edge_labels:XXX"
        print "START polygon"
        for i in self.polygons():
            "%s %s"

    def edit(self):
        r"""
        Launch the tk editor to interactively modify ``self``.
        """
        raise NotImplementedError("Pat should be working on it")
        #old version below
        #from translation_surface_editor import TranslationSurfaceEditor
        #fse = TranslationSurfaceEditor(self)
        #fse.window.mainloop()

    def plot(self):
        return TranslationSurfacePlot(self).plot()

    def edge_matrix(self, p, e=None):
        if e is None:
            p,e = p
        if p not in self.polygon_labels():
            raise ValueError
        elif e < 0 or e >= self.polygons()[p].num_edges():
            raise ValueError
        return identity_matrix(self.field())

    def stratum(self):
        from sage.dynamics.flat_surfaces.all import AbelianStratum
        return AbelianStratum([a-1 for a in self.angles()])

class Origami(TranslationSurface):
    def __init__(self, r, u):
        self._domain = r.parent().domain()
        self._r = r
        self._u = u
        self._perms = [~u,r,u,~r] # down,right,up,left

    def _repr_(self):
        return "Origami defined by r=%s and u=%s"%(self._r,self._u)

    def polygons(self):
        from polygon import square
        return Family(self._domain, lambda x: square())

    def polygon(self, lab):
        if lab not in self._domain:
            raise ValueError
        from polygon import square
        return square()

    def edges(self):
        return [(p,j) for p in self._domain for j in xrange(4)]

    def edge_to_identified(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e < 0 or e > 3:
            raise ValueError
        return self._perms[e](p), (e+2)%4

class TranslationSurfaceGenerators:
    r"""
    Common and less common translation surfaces.
    """
    @staticmethod
    def regular_octagon():
        r"""
        Return the translation surface built from the regular octagon by
        identifying opposite sides.

        EXAMPLES::

            sage: T = translation_surfaces.regular_octagon()
            sage: T
            Translation surface built from the regular octagon
            sage: T.stratum()
            H_2(2)
        """
        from polygon import regular_octagon
        polygons = [regular_octagon()]
        identifications = {}
        identifications.update(dict(((0,i),(0,i+4)) for i in xrange(4)))
        return TranslationSurface(polygons=polygons, identifications=identifications)

    @staticmethod
    def octagon_and_squares():
        from polygon import square, regular_octagon
        from sage.matrix.matrix_space import MatrixSpace

        o = regular_octagon()
        K = o.parent().field()
        sqrt2 = K.gen()
        rot = MatrixSpace(K,2)([[sqrt2/ZZ_2,-sqrt2/ZZ_2],[sqrt2/ZZ_2,sqrt2/ZZ_2]])
        polygons = [regular_octagon(), ZZ_2*square(K), ZZ_2*rot*square(K)]
        identifications = {
            (0,0): (1,3),
            (0,1): (2,3),
            (0,2): (1,0),
            (0,3): (2,0),
            (0,4): (1,1),
            (0,5): (2,1),
            (0,6): (1,2),
            (0,7): (2,2),
            }
        return TranslationSurface(polygons=polygons, identifications=identifications)

    @staticmethod
    def origami(r,u):
        r"""
        Return the origami defined by the permutations ``r`` and ``u``.

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: o
            Origami defined by r=(1,2) and u=(1,3)
            sage: o.stratum()
            H_2(2)
        """
        return Origami(r,u)

translation_surfaces = TranslationSurfaceGenerators()
