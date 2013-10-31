r"""
Translation surface in Sage.
Stupid edit here from Pat.

EXAMPLES:

Right now we can do from sage::

    sage: import translation_surface
    sage: t = TranslationSurface()
    sage: t.edit()
"""

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR

ZZ_1 = Integer(1)
ZZ_2 = Integer(2)

from polygon import Polygon,Polygons

class TranslationSurfaceDisplay:
    r"""
    Specify position of polygons in the plane.

    Internally it is a mapping from the index set of polygon to K^2 zhere K is a
    field embedded in RR^2.

    TODO:

        Think about a class for plotting saddle connections, cylinders, circles.
    """

class TranslationSurface(SageObject):
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
    def __init__(self, polygons=None, index_set=None, identifications=None, pos=None, edge_labels=None):
        r"""
        - ``data`` - a mapping from the index set to a set of polygons (it might
          be a list in which case the index set is considered as being the
          integer from ``0`` to ``leng(polygons)-1`` or a dictionnary or a
          function)

        - ``index_set`` - an index set (optional if ``data`` is a list or a
          dictionnary)

        - ``field`` - an optional field (optional if ``data`` is a list or a
          dictionnary)

        - ``edge_labels`` - a dictionnary or a function (polygon,edge_number) -> label

        - ``identifications`` - a dictionnary or a function
          (polygon,edge_number) -> (polygon',edge_number',matrix) that is an involution
          without fixed point reduced to the two first coordinates and such that
          the matrix map the edge to the other.

        - ``pos`` - an optional mapping from the index set to positions
        """
        if isinstance(polygons, list):
            from sage.sets.integer_range import IntegerRange
            self._polygons = polygons
            self._index_set = IntegerRange(len(polygons))
        elif isinstance(polygons, dict):
            from sage.sets.totally_ordered_finite_set import FiniteEnumeratedSets
            self._polygons = polygons
            self._index_set = FiniteEnumeratedSets(sorted(polygons.keys()))
        else:
            raise NotImplementedError("not yet implemented!! wait a day or a week, or a month or a year...")

        # the following is a bit dirty to determined the common field of the
        # polygons
        from sage.structure.sequence import Sequence
        common_field = Sequence([self._polygons[a].parent().field().an_element() for a in self._index_set]).universe()
        self._field = common_field
        for a in self._index_set:
            if self._polygons[a].parent().field() is not common_field:
                self._polygons[a] = Polygons(common_field)(self._polygon[a])

        self._identifications = identifications
        self._pos = pos
        self._edge_labels = None          # optional (for cutting sequences)

    def _repr_(self):
        return "Translation surface built from %d polygons"%self.num_polygons()

    def num_polygons(self):
        return len(self._polygons)

    def field(self):
        return self._field

    def volume(self):
        raise NotImplementedError

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
        from translation_surface_editor import TranslationSurfaceEditor
        fse = TranslationSurfaceEditor(self)
        fse.window.mainloop()

    def _set_edge_labels(self):
        self._edge_labels = {}
        m = 0
        for i in self._index_set:
            for j in xrange(self._polygons[i].num_edges()):
                if (i,j) not in self._edge_labels:
                    ii,jj = self._identifications[i,j]
                    self._edge_labels[i,j] = self._edge_labels[ii,jj] = m
                m += 1

    def _set_positions(self):
        if isinstance(self._polygons, dict):
            self._pos = {}
        elif isinstance(self._polygons, list):
            self._pos = [None] * len(self._polygons)
        for a in self._index_set:
            self._pos[a] = 10 * self.field().random_element()

    def plot(self):
        from sage.plot.graphics import Graphics

        if self._pos is None:
            self._set_positions()
        if self._edge_labels is None:
            self._set_edge_labels()

        G = Graphics()
        for a in self._index_set:
            p = self._polygons[a]
            v = p.vertices(translation=self._pos[a])
            G += p.plot(translation=self._pos[a])
            for i in xrange(p.num_edges()):
                G += text(str(self._edge_labels[a,i]),
                        (v[i]+v[(i+1)%p.num_edges()])/2,
                        color='black')

        # then we possibly plot the lable inside each polygon
        if len(self._index_set) != 1:
            for a in self._index_set:
                p = self._polygons[a]
                m = sum(p.vertices(translation=self._pos[a])) / p.num_edges()
                G += text(str(a), m, color='black')
        return G

class TranslationSurfaces:
    @staticmethod
    def regular_octagon():
        from polygon import regular_octagon
        polygons = [regular_octagon()]
        identifications = {}
        identifications.update(dict(((0,i),(0,i+4)) for i in xrange(4)))
        identifications.update(dict(((0,i+4),(0,i)) for i in xrange(4)))
        edge_labels = {(0,0): 'A', (4,0): 'A', (1,0): 'B', (5,0): 'B',
                     (2,0): 'C', (6,0): 'C', (3,0): 'D', (7,0): 'D'}
        return TranslationSurface(polygons=polygons,
                identifications=identifications,
                edge_labels=edge_labels)

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
        for k in identifications.keys():
            identifications[identifications[k]] = k
        pos = [(0,0),(0,0),(0,2)]
        return TranslationSurface(polygons=polygons,
                identifications=identifications, pos=pos)

translation_surfaces = TranslationSurfaces()
