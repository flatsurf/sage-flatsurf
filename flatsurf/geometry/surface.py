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

    def _add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.
        """
        raise NotImplementedError

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        raise NotImplementedError

    def _change_base_label(self, new_base_label):
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
    
    def change_polygon_gluings(self, label, glue_list):
        r"""
        Updates the list of glued polygons according to the provided list,
        which is a list of pairs (pp,ee) whose position in the list
        describes the edge of the polygon with the provided label.
        """
        self.__mutate()
        for e,(pp,ee) in enumerate(glue_list):
            self._change_edge_gluing(label, e, pp, ee)
    
    def add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Adds a the provided polygon to the surface. Utilizes gluing_list
        for the gluing data for edges (which must be a list
        of pairs of length equal to number of edges of the polygon).

        If the parameter label is provided, the Surface attempts to use
        this as the label for the new_polygon. However, this may fail 
        depending on the implementation.

        Returns the label assigned to the new_polygon (which may differ
        from the label provided).
        """
        self.__mutate()
        assert gluing_list is None or new_polygon.num_edges() == len(gluing_list)
        return self._add_polygon(new_polygon, gluing_list,label)

    def remove_polygon(self, label):
        r"""
        Remove the polygon with the provided label.
        """
        self.__mutate()
        return self._remove_polygon(label)

    def change_base_label(self, new_base_label):
        r"""
        Change the base_label to the provided label.
        """
        self.__mutate()
        return self._change_base_label(new_base_label)

    def __hash__(self):
        r"""
        Hash compatible with equals.
        """
        if self._s.is_mutable():
            raise ValueError("Attempting to hash mutable surface.")
        if not self._s.is_finite:
            raise ValueError("Attempting to hash infinite surface.")
        h = 73+17*hash(self.base_ring())+23*hash(self.base_label())
        for pair in self.label_iterator(polygons=True):
            h = h + 7*hash(pair)
        for edgepair in self.edge_iterator(gluings=True):
            h = h + 3*hash(edgepair)
        return h


    def __eq__(self, other):
        r"""
        Implements a naive notion of equality where two finite surfaces are equal if:
        - their base rings are equal,
        - their base labels are equal,
        - their polygons are equal and labeled and glued in the same way.
        For infinite surfaces we use reference equality.
        """
        if self is other:
            return True
        if not isinstance(other, Surface):
            raise TypeError
        if not self.is_finite():
            if other.is_finite():
                return False
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

####
#### Surface_list
####

class Surface_list(Surface):
    r"""
    A fast mutable implementation of surface.
    """
    def __init__(self, base_ring=None, surface = None, copy=True, mutable=None):
        r"""
        Surface_list is a Surface implementation which uses int for labels.
        (Internally, things are stored in a list.)
        
        Parameters
        ----------
        base_ring : field or None
            Field containing the vertices of the polygons. If left None (as default)
            the base_ring is copied from the surface. If surface is also None
            a ValueError will be raised.
        surface : Surface or SimilaritySurface or None
            A surface to be copied or referenced to obtain the initial state. 
            If surface is None (as default), the surface is initialized to be 
            empty and mutability is forced.
        copy : boolean
            If copy is true and surface is not None, then a copy is made of the 
            surface. This will result in an error if surface is infinite. If
            copy is true, an reference to surface is used instead. Results of 
            querying are cached, therefore trying to store a reference of a
            mutable surface is not allowed. The main purpose of the option
            copy=False is to obtain a "copy" of an infinite surface. In case
            copy=False, the labels for polygons are decided on in a lazy way.
            In particular, you need to explore the resulting surface before
            calling polygon(100). (You can explore the surface using any 
            of the iterators.)
        mutable : boolean or None
            If mutable is true, the resulting surface will be mutable. If mutable
            is false, then the resulting surface will not be mutable. If mutable
            is left at its default value of None, then the surface will be mutable
            if and only if a surface is not provided.

        EXAMPLES::

            sage: from flatsurf import *
            sage: print("We will surgically add a square into an infinite billiard surface")
            We will surgically add a square into an infinite billiard surface
            sage: p = polygons(vertices=[(0,0),(4,0),(0,3)])
            sage: s = similarity_surfaces.billiard(p)
            sage: ts=s.minimal_translation_cover().copy(relabel=True, mutable=True)
            Warning: Could be indicating infinite surface falsely.
            sage: # Explore the surface a bit
            sage: ts.polygon(0)
            Polygon: (0, 0), (4, 0), (0, 3)
            sage: ts.opposite_edge(0,0)
            (1, 2)
            sage: ts.polygon(1)
            Polygon: (0, 0), (0, -3), (4, 0)
            sage: s = ts.underlying_surface()
            sage: l=s.add_polygon(polygons.square(side=4))
            sage: s.change_edge_gluing(0,0,l,2)
            sage: s.change_edge_gluing(1,2,l,0)
            sage: s.change_edge_gluing(l,1,l,3)
            sage: print("Glued in label is "+str(l))
            Glued in label is 2
            sage: count = 0
            sage: for x in ts.edge_iterator(gluings=True):
            ....:     print(x)
            ....:     count=count+1
            ....:     if count>15:
            ....:         break
            ((0, 0), (2, 2))
            ((0, 1), (3, 1))
            ((0, 2), (4, 0))
            ((2, 0), (1, 2))
            ((2, 1), (2, 3))
            ((2, 2), (0, 0))
            ((2, 3), (2, 1))
            ((3, 0), (5, 2))
            ((3, 1), (0, 1))
            ((3, 2), (6, 0))
            ((4, 0), (0, 2))
            ((4, 1), (7, 1))
            ((4, 2), (8, 0))
            ((1, 0), (8, 2))
            ((1, 1), (9, 1))
            ((1, 2), (2, 0))
            sage: count = 0
            sage: for l,p in ts.label_iterator(polygons=True):
            ....:     print(str(l)+" -> "+str(p))
            ....:     count=count+1
            ....:     if count>5:
            ....:         break
            0 -> Polygon: (0, 0), (4, 0), (0, 3)
            2 -> Polygon: (0, 0), (4, 0), (4, 4), (0, 4)
            3 -> Polygon: (0, 0), (-72/25, -21/25), (28/25, -96/25)
            4 -> Polygon: (0, 0), (0, 3), (-4, 0)
            1 -> Polygon: (0, 0), (0, -3), (4, 0)
            5 -> Polygon: (0, 0), (-28/25, 96/25), (-72/25, -21/25)
        """
        self._base_label = 0
        self._p = []
        self._reference_surface = None # Whether or not we store a reference surface
        self._removed_labels = []
        if surface is None:
            if base_ring is None:
                raise ValueError("Either surface or base_ring must be provided.")
            self._base_ring= base_ring
            if not mutable is None and not mutable:
                raise ValueError("If no surface is provided, then mutable must be true.")
            self._num_polygons=0
            Surface.__init__(self, mutable=True)
        else:
            from flatsurf.geometry.similarity_surface import SimilaritySurface
            if isinstance(surface,SimilaritySurface):
                surface=surface.underlying_surface()
            if not isinstance(surface,Surface):
                raise ValueError("surface must be either a Surface or SimilaritySurface")
            if not base_ring is None and base_ring != surface.base_ring():
                raise ValueError("You currently can not provide both a surface and a base_ring.")
            self._base_ring = surface.base_ring()
            if copy==True:
                if not surface.is_finite():
                    raise ValueError("Can not copy an infinite surface.")
                Surface.__init__(self, mutable=True)
                self._num_polygons=0
                label_dict = {}
                p = 0
                for label,polygon in surface.label_polygon_iterator():
                    label_dict[label] = p
                    p=p+1
                    # This automatically adds one to _num_polygons:
                    self.add_polygon(polygon)
                for pair1,pair2 in surface.edge_gluing_iterator():
                    l1,e1=pair1
                    l2,e2=pair2
                    ll1=label_dict[l1]
                    ll2=label_dict[l2]
                    self._p[ll1][1][e1]=(ll2,e2)
                if mutable is None or not mutable:
                    self.make_immutable()
            else:
                if surface.is_mutable():
                    raise ValueError("Surface_list will not store reference to a mutable surface.")
                self._reference_surface = surface
                self._ref_to_int={}
                self._int_to_ref=[]
                self.__get_label(surface.base_label())
 
                # Cache the base polygon
                polygon = surface.polygon(surface.base_label())
                self._p[0]=[polygon, [None for i in xrange(polygon.num_edges())]]

                self._num_polygons = self._reference_surface.num_polygons()

                # Take care of is_finite!
                self.is_finite = self._reference_surface.is_finite
                
                Surface.__init__(self, mutable=not mutable is None and mutable)

                # Because of the reference, we can't do better than exploring the surface:
                #self.label_iterator = super(Surface, self).label_iterator

    def __get_label(self, ref_label):
        r"""
        Returns a corresponding label. Creates a new label if necessary.
        """
        try:
            return self._ref_to_int[ref_label]
        except KeyError:
            polygon = self._reference_surface.polygon(ref_label)
            data = [polygon,[None for i in xrange(polygon.num_edges())]]
            if len(self._removed_labels)>0:
                i = self._removed_labels.pop()
                self._p[i]=data
            else:
                i = len(self._p)
                self._p.append(data)
            self._ref_to_int[ref_label]=i
            self._int_to_ref.append(ref_label)
            return i

    def base_ring(self):
        r"""
        The field in which the coordinates of polygons live.
        """
        return self._base_ring

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        try:
            data = self._p[lab]
        except IndexError:
            raise ValueError("No known polygon with provided label "+str(lab)+". "+\
                "This can be caused by failing to explore your surface. "+\
                "See the documentation in flatsurf.geometry.surface.Surface_list.")
        try:
            return data[0]
        except TypeError:
            raise ValueError("Provided label was removed.")

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
        try:
            data = self._p[p]
        except KeyError:
             raise ValueError("No known polygon with provided label")
        try:
            glue = data[1]
        except TypeError:
            raise ValueError("Provided label was removed.")
        try:
            oe = glue[e]
        except KeyError:
            raise ValueError("Invalid edge")
        if oe is None:
            if self._reference_surface is None:
                # Perhaps the user of this class left an edge unglued?
                return None
            else:
                ref_p = self._int_to_ref[p]
                ref_pp, ref_ee = self._reference_surface.opposite_edge(ref_p, e)
                pp= self.__get_label(ref_pp)
                return_value = (pp, ref_ee)
                glue[e] = return_value
                return return_value
        else:
            # Sucessfully return edge data
            return oe

    def is_finite(self):
        # Implementation Note:
        # If we are storing a reference, this version of is_finite gets overriden.
        # (See the __init__ method.)
        return True

    # Methods for changing the surface

    def _change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Internal method used by change_polygon(). Should not be called directly.
        """
        try:
            data = self._p[label]
        except KeyError:
            raise ValueError("No known polygon with provided label")
        if data is None:
            raise ValueError("Provided label was removed from the surface.")
        data[0]=new_polygon
        if data[1] is None or new_polygon.num_edges() != len(data[1]):
            data[1]=[None for e in xrange(new_polygon.num_edges())]
        if not gluing_list is None:
            self.change_polygon_gluings(label,gluing_list)

    def _change_edge_gluing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.
        """
        try:
            data = self._p[label1]
        except KeyError:
            raise ValueError("No known polygon with provided label1="+str(label1))
        if data is None:
            raise ValueError("Provided label1="+str(label1)+" was removed from the surface.")
        data[1][edge1]=(label2,edge2)
        try:
            data = self._p[label2]
        except KeyError:
            raise ValueError("No known polygon with provided label2="+str(label2))
        if data is None:
            raise ValueError("Provided label2="+str(label2)+" was removed from the surface.")
        data[1][edge2]=(label1,edge1)

    def _add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.
        """
        if new_polygon is None:
            data=[None,None]
        else:
            data=[new_polygon, [None for i in xrange(new_polygon.num_edges())] ]
        if len(self._removed_labels)>0:
            new_label = self._removed_labels.pop()
            self._p[new_label]=data
        else:
            new_label = len(self._p)
            self._p.append(data)
            if not self._reference_surface is None:
                # Need a blank in this list for algorithmic reasons
                self._int_to_ref.append(None)
        if not gluing_list is None:
            self.change_polygon_gluings(new_label,gluing_list)
        self._num_polygons += 1
        return new_label
            
    def num_polygons(self):
        r""" 
        Return the number of polygons making up the surface in constant time.
        """
        return self._num_polygons
    
    def label_iterator(self):
        r"""
        Iterator over all polygon labels.
        """
        if not self._reference_surface is None:
            for i in Surface.label_iterator(self):
                yield i
        elif self._num_polygons == len(self._p):
            for i in xrange(self.num_polygons()):
                yield i
        else:
            # We've removed some labels
            found=0
            i=0
            while found < self._num_polygons:
                if not self._p[i] is None:
                    #print "found",i
                    found += 1
                    yield i
                i += 1

    def _change_base_label(self, new_base_label):
        r"""
        Internal method for change_base_label. Should not be called directly.
        """
        self._base_label =  new_base_label

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        if label == len(self._p)-1:
            last = self._p.pop()
        else:
            self._p[label]=None
            self._removed_labels.append(label)
        self._num_polygons -= 1

####
#### Surface_dict
####

class Surface_dict(Surface):
    r"""
    A fast mutable implementation of surface using a dictionary.

    EXAMPLES::

        from flatsurf import *
        from flatsurf.geometry.surface import Surface_dict
        p=polygons.regular_ngon(10)
        s=Surface_dict(base_ring=p.base_ring())
        s.add_polygon(p,label="A")
        s.change_polygon_gluings("A",[("A",(e+5)%10) for e in xrange(10)])
        s.change_base_label("A")
        TestSuite(s).run()
    """
    ###
    ### Brief summary of internal workings.
    ###
    #
    # The Surface_dict maintains a dictionary self._p for storing polygons together with gluing data.
    # Here self._p[label] is typically a list of two elements [polygon, gluing_list].
    # The gluing_list is then a list of pairs (other_label, other_edge) so that typically
    # self.opposite_edge(label, edge) returns self._p[label][1][edge].
    #
    # If constructed with a surface parameter which is not None and copy=False, then Surface_dict
    # stores a reference to the provided surface as self._reference_surface. 
    # (Otherwise self._reference_surface is None). If we have a reference surface, then to represent a removed 
    # surface we set self._[label]=None.
    #
    def __init__(self, base_ring=None, surface = None, copy=True, mutable=None):
        r"""
        Surface_list is a Surface implementation which uses int for labels.
        (Internally, things are stored in a list.)
        
        Parameters
        ----------
        base_ring : field or None
            Field containing the vertices of the polygons. If left None (as default)
            the base_ring is copied from the surface. If surface is also None
            a ValueError will be raised.
        surface : Surface or SimilaritySurface or None
            A surface to be copied or referenced to obtain the initial state. 
            If surface is None (as default), the surface is initialized to be 
            empty and mutability is forced.
        copy : boolean
            If copy is true and surface is not None, then a copy is made of the 
            surface. This will result in an error if surface is infinite. If
            copy is true, an reference to surface is used instead. Results of 
            querying are cached, therefore trying to store a reference of a
            mutable surface is not allowed. The main purpose of the option
            copy=False is to obtain a "copy" of an infinite surface. In case
            copy=False, the labels for polygons are decided on in a lazy way.
            In particular, you need to explore the resulting surface before
            calling polygon(100). (You can explore the surface using any 
            of the iterators.)
        mutable : boolean or None
            If mutable is true, the resulting surface will be mutable. If mutable
            is fals e, then the resulting surface will not be mutable. If mutable
            is left at its default value of None, then the surface will be mutable
            if and only if a surface is not provided.
        """
        self._p = {}
        self._reference_surface = None
        self._base_label = None
        if surface is None:
            if base_ring is None:
                raise ValueError("Either surface or base_ring must be provided.")
            self._base_ring = base_ring
            if not mutable is None and not mutable:
                raise ValueError("If no surface is provided, then mutable must be true.")
            Surface.__init__(self, mutable=True)
        else:
            from flatsurf.geometry.similarity_surface import SimilaritySurface
            if isinstance(surface,SimilaritySurface):
                surface=surface.underlying_surface()
            if not isinstance(surface,Surface):
                raise ValueError("surface must be either a Surface or SimilaritySurface")
            if not base_ring is None and base_ring != surface.base_ring():
                raise ValueError("You currently can not provide both a surface and a base_ring.")
            self._base_ring = surface.base_ring()
            self._base_label = surface.base_label()
            if copy==True:
                if not surface.is_finite():
                    raise ValueError("Can not copy an infinite surface.")
                for label,polygon in surface.label_polygon_iterator():
                    self._p[label]=[polygon, \
                        [ surface.opposite_edge(label,e) for e in xrange(polygon.num_edges()) ] ]
                # The only way we're mutable is if mutable=True:
                Surface.__init__(self, mutable=not mutable is None and mutable)
            else:
                if surface.is_mutable():
                    raise ValueError("Surface_dict will not store reference to a mutable surface.")
                self._reference_surface = surface 
                # Take care of is_finite!
                self.is_finite = self._reference_surface.is_finite
                # The only way we're mutable is if mutable=True:
                Surface.__init__(self, mutable=not mutable is None and mutable)

    def base_ring(self):
        r"""
        The field in which the coordinates of polygons live.
        """
        return self._base_ring

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        try:
            data = self._p[lab]
        except KeyError:
            # No label in dictionary.
            if self._reference_surface is None:
               raise ValueError("No known polygon with provided label "+str(lab)+".")
            else:
                polygon = self._reference_surface.polygon(lab)
                data = [ self._reference_surface.polygon(lab), \
                    [ self._reference_surface.opposite_edge(lab,e) for e in xrange(polygon.num_edges()) ] ]
                self._p[lab] = data
        if data is None:
            raise ValueError("Label "+str(label)+" was removed from the surface.")
        return data[0]

    def base_label(self):
        r"""
        Returns the base label.
        """
        return self._base_label

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        try:
            data = self._p[p]
        except KeyError:
             self.polygon(p)
             data=self._p[p]
        if data is None:
            raise ValueError("Label "+str(label)+" was removed from the surface.")
        gluing_data=data[1]
        try:
            return gluing_data[e]
        except IndexError:
            raise ValueError("Edge e="+str(e)+" is out of range in polygon with label "+str(p))

    def is_finite(self):
        # Implementation Note:
        # If we are storing a reference, this version of is_finite gets overriden.
        # (See the __init__ method.)
        return True

    # Methods for changing the surface

    def _change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Internal method used by change_polygon(). Should not be called directly.
        """
        try:
            data = self._p[label]
            if data is None:
                raise ValueError("Label "+str(label)+" was removed from the surface.")
            data[0]=new_polygon
        except KeyError:
            # Polygon probably lies in reference surface
            if self._reference_surface is None:
                raise ValueError("No known polygon with provided label")
            else:
                # Ensure the reference surface had a polygon with the provided label:
                old_polygon = self._reference_surface.polygon(label)
                if old_polygon.num_edges() == new_polygon.num_edges():
                    self._p[label]=[new_polygon, \
                        [self._reference_surface.opposite_edge(label,e) for e in xrange(new_polygon.num_edges())] ]
                else:
                    self._p[label]=[new_polygon, [None for e in xrange(new_polygon.num_edges())] ]
        if len(data[1]) != new_polygon.num_edges():
            data[1] = [None for e in xrange(new_polygon.num_edges())]
        if not gluing_list is None:
            self.change_polygon_gluings(label,gluing_list)

    def _change_edge_gluing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.
        """
        try:
            data = self._p[label1]
        except KeyError:
            raise ValueError("No known polygon with provided label1="+str(label1))
        try:
            data[1][edge1]=(label2,edge2)
        except IndexError:
            # break down error
            if data is None:
                raise ValueError("Polygon with label1="+str(label1)+" was removed.")
            try:
                data1=data[1]
            except IndexError:
                raise RuntimeError("This index error should not happen.")
            try:
                data1[edge1]=(label2,edge2)
            except IndexError:
                raise ValueError("Edge "+str(edge1)+" is out of range in polygon with label1="+str(label1))
        try:
            data = self._p[label2]
        except KeyError:
            raise ValueError("No known polygon with provided label2="+str(label2))
        try:
            data[1][edge2]=(label1,edge1)
        except IndexError:
            # break down error
            if data is None:
                raise ValueError("Polygon with label1="+str(label1)+" was removed.")
            try:
                data1=data[1]
            except IndexError:
                raise RuntimeError("This index error should not happen.")
            try:
                data1[edge2]=(label1,edge1)
            except IndexError:
                raise ValueError("Edge "+str(edge2)+" is out of range in polygon with label1="+str(label2))

    def _add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.
        """
        data=[new_polygon, [None for i in xrange(new_polygon.num_edges())] ]
        if label is None:
            new_label = ExtraLabel()
        else:
            if label in self._p and not self._p[label] is None:
                raise ValueError("We already have a polygon with label="+str(label))
            elif not self._reference_surface is None:
                raise ValueError("Can not assign this label to a Surface_dict containing a reference surface. ")
            new_label=label
        self._p[new_label]=data
        if not gluing_list is None:
            self.change_polygon_gluings(new_label,gluing_list)
        return new_label
    
    def num_polygons(self):
        r""" 
        Return the number of polygons making up the surface in constant time.
        """
        if self.is_finite():
            if self._reference_surface is None:
                return len(self._p)
            else:
                # Unfortunately, I don't see a good way to compute this.
                return Surface.num_polygons(self)
        else:
            from sage.rings.infinity import Infinity
            return Infinity

    def label_iterator(self):
        r"""
        Iterator over all polygon labels.
        """
        if self._reference_surface is None:
            for i in self._p:
                yield i
        else:
            for i in Surface.label_iterator(self):
                yield i

    def _change_base_label(self, new_base_label):
        r"""
        Internal method for change_base_label. Should not be called directly.
        """
        self._base_label =  new_base_label

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        if self._reference_surface is None:
            try:
                data = self._p[label]
            except KeyError:
                raise ValueError("Label "+str(label)+" is not in the surface.")
            if data is None:
                raise ValueError("Label "+str(label)+" was already removed from the surface.")
            del self._p[label]
        else:
            try:
                data = self._p[label]
                # Success. 
                if data is None:
                    raise ValueError("Label "+str(label)+" was already removed from the surface.")
                self._p[label] = None
            except KeyError:
                # Assume on faith we are removing a polygon in the base_surface.
                self._p[label] = None


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

######
###### ExtraLabels
######

class ExtraLabel(SageObject):
    r""" 
    Used to spit out new labels.
    """
    _next=int(0)
    
    def __init__(self):
        r"""
        Construct a new label.
        """
        self._label = int(ExtraLabel._next)
        ExtraLabel._next = ExtraLabel._next + 1
    
    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self._label == other._label)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(23*self._label)
        
    def __str__(self):
        return "E"+str(self._label)
    
    def __repr__(self):
        return "ExtraLabel("+str(self._label)+")"

