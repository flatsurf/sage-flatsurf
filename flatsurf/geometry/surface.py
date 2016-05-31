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
    
    #class CachedData:
    #    r"""
    #    This class is just for storing cached data, e.g. the area, number of polygons, number of edges, etc.
    #    """
    #    pass
    
    def __init__(self, mutable=False):
        self._mutable = mutable
        #self._cache=Surface.CachedData()
        self._cache = {}
        
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
    # If the surface can be changed, implement the following methods:
    #

        
    def _change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Internal method used by change_polygon(). Should not be called directly.
        """
        raise NotImplementedError
    
    def _change_edge_gluing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.
        """
        raise NotImplementedError

    def _add_polygon(self, new_polygon, gluing_list=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.
        """
        raise NotImplementedError

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        raise NotImplementedError

    def _change_base_polygon(self, new_base_label):
        r"""
        Internal method for change_base_label. Should not be called directly.
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
            lw=self.walker()
            lw.find_all_labels()
            return len(lw)
        else:
            from sage.rings.infinity import Infinity
            return Infinity

    def label_iterator(self):
        r"""
        Iterator over all polygon labels.
        """
        return iter(self.walker())

    def label_polygon_iterator(self):
        r"""
        Iterate over pairs (label,polygon).
        """
        for label in self.label_iterator():
            yield label, self.polygon(label)

    #
    # Methods which you probably do not want to override.
    #

    def num_edges(self):
        r"""
        Return the total number of edges of all polygons used.
        
        Not cached. Will be likely be linear in the number of edges.
        """
        if self.is_finite():
            try:
                return self._cache["num_edges"]
            except KeyError:
                num_edges = sum(p.num_edges() for l,p in self.label_polygon_iterator())
                self._cache["num_edges"] = num_edges
                return num_edges
        else:
            from sage.rings.infinity import Infinity
            return Infinity

    def area(self):
        r"""
        Return the area of this surface.
        
        By default this method is not cached.
        """
        if self.is_finite():
            try:
                return self._cache["area"]
            except KeyError:
                area = sum(p.area() for l,p in self.label_polygon_iterator())
                self._cache["area"] = area
                return area
        raise NotImplementedError("area is not implemented for surfaces built from an infinite number of polygons")

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
            yield (label_edge_pair, \
                self.opposite_edge(label_edge_pair[0], label_edge_pair[1]))

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

    def walker(self):
        r"""
        Return a LabelWalker which walks over the surface in a canonical way.
        """
        try:
            return self._cache["lw"]
        except KeyError:
            lw = LabelWalker(self)
            self._cache["lw"] = lw
            return lw

    def __mutate(self):
        r"""
        Called before a mutation occurs. Do not call directly.
        """
        assert(self.is_mutable())
        # Remove the cache which will likely be invalidated.
        #self._cache=CachedData()
        self._cache = {}

    def change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Assuming label is currently in the list of labels, change the 
        poygon assigned to the provided label to new_polygon, and 
        glue the edges according to gluing_list (which must be a list
        of pairs of length equal to number of edges of the polygon). 
        
        Warning: the gluing_list may be incorporated by reference.
        """
        self.__mutate()
        assert gluing_list is None or new_polygon.num_edges() == len(gluing_list)
        self._change_polygon(label, new_polygon, gluing_list)

    def change_edge_gluing(self, label1, edge1, label2, edge2):
        r"""
        Updates the gluing so that (label,edge1) is glued to (label2, edge2) 
        and vice versa.
        """
        self.__mutate()
        self._change_edge_gluing(label1, edge1, label2, edge2)
    
    def add_polygon(self, new_polygon, gluing_list=None):
        r"""
        Adds a the provided polygon to the surface. Utilizes gluing_list
        for the gluing data for edges (which must be a list
        of pairs of length equal to number of edges of the polygon).
        
        Return the label assigned to the new_polygon.
        """
        self.__mutate()
        assert gluing_list is None or new_polygon.num_edges() == len(gluing_list)
        return self._add_polygon(new_polygon, gluing_list)

    def remove_polygon(self, label):
        r"""
        Remove the polygon with the provided label.
        """
        self.__mutate()
        return self._remove_polygon(label)

    def change_base_polygon(self, new_base_label):
        r"""
        Change the base_label to the provided label.
        """
        self.__mutate()
        return self._change_base_polygon(new_base_label)

class Surface_fast(Surface):
    r"""
    A fast mutable implementation of surface that does not support polygon removal.
    """
    def __init__(self, surface = None, base_ring=None, mutable=None):
        r"""
        Needs documentation
        """
        self._base_label = 0
        if surface is None:
            if base_ring is None:
                raise ValueError("Either surface or base_ring must be provided.")
            self._base_ring= base_ring
            if not mutable is None:
                if not mutable:
                    raise ValueError("If no surface is provided, then mutable must be true.")
            self._p = []
            Surface.__init__(self, mutable=True)
        else:
            from flatsurf.geometry.similarity_surface import SimilaritySurface
            if isinstance(surface,SimilaritySurface):
                surface=surface.underlying_surface()
            if not isinstance(surface,Surface):
                raise ValueError("surface must be either a Surface or SimilaritySurface")
            if not base_ring is None:
                raise ValueError("You currently can not provide both a surface and a base_ring.")
            self._base_ring = surface.base_ring()
            if not surface.is_finite():
                raise ValueError("Can not copy an infinite surface.")
            Surface.__init__(self, mutable=True)
            self._p=[]
            label_dict = {}
            p = 0
            for label,polygon in surface.label_polygon_iterator():
                label_dict[label] = p
                p=p+1
                self.add_polygon(polygon)
            for pair1,pair2 in surface.edge_gluing_iterator():
                l1,e1=pair1
                l2,e2=pair2
                ll1=label_dict[l1]
                ll2=label_dict[l2]
                self._p[ll1][1][e1]=(ll2,e2)
            if mutable is None or not mutable:
                self.make_immutable()

    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        return self._base_ring

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        return self._p[lab][0]

    def base_label(self):
        r"""
        Always returns the same label.
        """
        return self._base_label

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        return self._p[p][1][e]

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        return True

    # Methods for changing the surface

    def _change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Internal method used by change_polygon(). Should not be called directly.
        """
        self._p[label][0]=new_polygon
        if gluing_list is None:
            self._p[label][1]=[None for e in xrange(new_polygon.num_edges())]
        else:
            if not isinstance(gluing_list,list):
                raise ValueError("gluing_list must be None or a list.")
            self._p[label][1]=gluing_list
            for e in range(new_polygon.num_sides()):
                pair = gluing_list[e]
                if not pair is None:
                    l2,e2 = pair
                    self._p[l2][1][e2] = (label, e)
    
    def _change_edge_gluing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.
        """
        self._p[label1][1][edge1]=(label2,edge2)
        self._p[label2][1][edge2]=(label1,edge1)
            

    def _add_polygon(self, new_polygon, gluing_list=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.
        """
        new_label = len(self._p)
        if gluing_list is None:
            self._p.append([new_polygon,\
                [None for i in xrange(new_polygon.num_edges())]])
        else:
            if not isinstance(gluing_list,list):
                raise ValueError("gluing_list must be None or a list.")
            self._p.append([new_polygon,gluing_list])
    
    def num_polygons(self):
        r""" 
        Return the number of polygons making up the surface in constant time.
        """
        return len(self._p)
    
    def label_iterator(self):
        r"""
        Iterator over all polygon labels.
        """
        return iter(xrange(self.num_polygons()))

    def _change_base_polygon(self, new_base_label):
        r"""
        Internal method for change_base_label. Should not be called directly.
        """
        self._base_label =  new_base_label


class Surface_polygons_and_gluings(Surface):
    r"""
    Similarity surface build from a list of polygons and gluings.
    """
    def __init__(self, *args, **kwds):
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
        Surface.__init__(self)
    
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
