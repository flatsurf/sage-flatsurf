from sage.structure.sage_object import SageObject

from sage.sets.family import Family

class Surface(SageObject):
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
    
    def __init__(self, mutable=False):
        self._mutable = mutable
    
    # Do we really want to inherit from SageObject?
    
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

    # 
    # generic methods: you might consider overriding these for speed or for enabling mutability.
    #

    def num_polygons(self):
        r""" 
        Return the number of polygons making up the surface, or sage.rings.infinity.Infinity if the surface is infinite.
        
        This is a generic method. On a finite surface it will be linear time in the edges the first time it is run, then constant time (assuming not mutation occurs).
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
            raise RunTimeWarning("This default implementation of Surface.num_edges() is slow when called multiple times.")
            return sum(p.num_edges() for p in self.polygon_iterator())
        else:
            from sage.rings.infinity import Infinity
            return Infinity

    def label_iterator(self):
        r"""
        Iterator over all polygon labels.
        """
        return iter(self.label_walker())

    def label_polygon_iterator(self):
        r"""
        Iterate over pairs (label,polygon).
        """
        for label in self.label_iterator():
            yield label, self.polygon(label)

    def edge_iterator(self):
        r"""
        Iterate over the edges of polygons, which are pairs (l,e) where l is a polygon label, 0 <= e < N and N is the number of edges of the polygon with label l.
        """ 
        for label,polygon in self.label_polygon_iterator():
            for edge in xrange(polygon.num_edges()):
                yield label, edge

    def edge_gluing_iterator(self):
        r"""
        Iterate over the ordered pairs of edges being glued.
        """
        for label_edge_pair in self.edge_iterator():
            yield (label_edge_pair, self.opposite_edge_pair(label_edge_pair))

    def area(self):
        r"""
        Return the area of this surface.
        """
        if self.is_finite():
            return sum(p.area() for p in self.polygon_iterator())
        raise NotImplementedError("area is not implemented for surfaces built from an infinite number of polygons")


    #
    # Methods which you probably do not want to override.
    #

    def opposite_edge_pair(self,label_edge_pair):
        r"""
        Same as opposite_edge(label, pair) except taking a pair as input.
        """
        return self.opposite_edge(label_edge_pair[0], label_edge_pair[1])

    #
    # Methods which should not be overriden
    #

    def is_mutable(self):
        r"""
        Return if this surface is mutable.
        """
        return self._mutable

    def make_immutable(self):
        r"""
        Mark this surface as immutable. 
        """
        self._mutable = False

    def label_walker(self):
        try:
            return self._lw
        except AttributeError:
            self._lw = LabelWalker(self)
        return self._lw

    def _mutate(self):
        r"""
        Call before a mutation occurs.
        """
        assert(self.is_mutable())
        # Remove the label walker because it will be invalid.
        del self._lw


class Surface_polygons_and_gluings(Surface):
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
          ((p0,e0),(p1,e1)).
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
    
    def base_ring(self):
        return self._field

    def polygon_labels(self):
        from sage.combinat.words.alphabet import build_alphabet
        return build_alphabet(self._polygons.keys())

    def base_label(self):
        return self.polygon_labels().an_element()

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

#####
##### LABEL WALKER
#####

from collections import deque

class LabelWalker:
    r"""
    Take a canonical walk around the surface and find the labels of polygons.

    We start at the base_label().
    Then the labels are visited in order involving the combinatorial distance from the base_label(),
    where combinatorial distance measures the minimal number of edges which need to be crossed to reach the
    polygon with a givel label. Ties are broken using lexigraphical order on the numbers associated to edges crossed
    (labels are not used in this lexigraphical ordering).
    """
    
    class LabelWalkerIterator:
        def __init__(self, label_walker):
            self._lw = label_walker
            self._i = 0
        
        def next(self):
            if self._i < len(self._lw):
                label = self._lw.number_to_label(self._i)
                self._i = self._i +1
                return label
            if self._i == len(self._lw):
                label = self._lw.find_a_new_label()
                if label is None:
                    raise StopIteration()
                self._i = self._i+1
                return label
            raise StopIteration()
            
        def __iter__(self):
            return self
    
    def __init__(self, surface):
        self._s=surface
        self._labels=[self._s.base_label()]
        self._label_dict={self._s.base_label():0}
        
        # This will stores an edge to move through to get to a polygon closer to the base_polygon
        self._label_edge_back = {self._s.base_label():None} 
        
        self._walk=deque()
        self._walk.append((self._s.base_label(),0))
    
    def label_dictionary(self):
        r""" 
        Return a dictionary mapping labels to integers which gives a canonical order on labels.
        """
        return self._label_dict
    
    def edge_back(self, label, limit=None):
        r"""
        Return the `canonical' edge to walk through to get closer to the base_label, 
        or None if label already is the base_label.
        
        Remark: This could be slow on infinite surfaces!
        """
        try:
            return self._label_edge_back[label]
        except KeyError:
            if limit is None:
                if not self._s.is_finite():
                    limit=1000
                else:
                    limit=self._s.num_polygons()
            for i in xrange(limit):
                new_label=self.find_a_new_label()
                if label == new_label:
                    return self._label_edge_back[label]
        # Maybe the surface is not connected?
        raise KeyError("Unable to find label %s. Are you sure the surface is connected?"%(label))
    
    def __iter__(self):
        return LabelWalker.LabelWalkerIterator(self)

    def polygon_iterator(self):
        for label in self:
            yield self._s.polygon(label)

    def label_polygon_iterator(self):
        for label in self:
            yield label, self._s.polygon(label)

    def edge_iterator(self):
        for label,polygon in self.label_polygon_iterator():
            for e in xrange(polygon.num_edges()):
                yield label,e

    def __len__(self):
        r"""
        Return the number of labels found.
        """
        return len(self._labels)

    def find_a_new_label(self):
        r"""
        Finds a new label, stores it, and returns it. Returns None if we have already found all labels.
        """
        while len(self._walk)>0:
            label,e = self._walk.popleft()
            opposite_label,opposite_edge=self._s.opposite_edge(label,e)
            e=e+1
            if e < self._s.polygon(label).num_edges():
                self._walk.appendleft((label,e))
            if not opposite_label in self._label_dict:
                n=len(self._labels)
                self._labels.append(opposite_label)
                self._label_dict[opposite_label]=n
                self._walk.append((opposite_label,0))
                self._label_edge_back[opposite_label]=opposite_edge
                return opposite_label
        return None

    def find_new_labels(self,n):
        r"""
        Look for n new labels. Return the list of labels found.
        """
        new_labels = []
        for i in range(n):
            label = self.find_a_new_label()
            if label is None:
                return new_labels
            else:
                new_labels.append(label)
        return new_labels
        
        
    def find_all_labels(self):
        assert(self._s.is_finite())
        label = self.find_a_new_label()
        while not label is None:
            label = self.find_a_new_label()
            
    def number_to_label(self, n):
        r"""
        Return the n-th label where n is less than the length.
        """
        return self._labels[n]
    
    def label_to_number(self, label):
        r"""
        Return the number associated to the provided label (which must have already been found).
        """
        return self._label_dict[label]

    def surface(self):
        return self._s
