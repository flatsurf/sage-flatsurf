from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

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
        can either use Surface_list or Surface_dict. Another option is to
        inherit from Surface and implement the two methods:

        - polygon(self, label): the polygon associated to the provided label
        - opposite_edge(self, label, edge): a couple (``other_label``, ``other_edge``)
          representing the edge being glued

        If you want the surface to be mutable, you should override the methods:
        - _change_polygon(self, label, new_polygon, gluing_list=None)
        - _change_edge_gluing(self, label1, edge1, label2, edge2)
        - _add_polygon(self, new_polygon, gluing_list=None, label=None)
        - _remove_polygon(self, label)
        See the documentation of those methods for details.
    """
    def __init__(self, base_ring, base_label, finite=None, mutable=False):
        r"""
        Represents a surface defined using polygons whose vertices lie
        in the provided base_ring.

        Parameters
        ----------
        base_ring : field
            Field containing the vertices of the polygons.
        base_label :
            A preferred label for a polygon in the surface.
        finite : boolean
            The truth value of the statement "The surface is finite."
        mutable : boolean
            If mutable is true, the resulting surface will be mutable.
        """
        if finite is None:
            raise ValueError("finite must be either True or False")
        self._base_ring = base_ring
        self._base_label = base_label
        self._finite=finite
        self._mutable = mutable
        self._cache = {}

    def polygon(self, label):
        r"""
        Return the polygon with the provided label.

        This method must be overriden in subclasses.
        """
        raise NotImplementedError

    def opposite_edge(self, l, e):
        r"""
        Given the label ``l`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``ll``, ``ee``) to which this edge is glued.

        This method must be overriden in subclasses.
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

    #
    # generic methods: override these for speed if possible
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
            for edge in range(polygon.num_edges()):
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

    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        return self._base_ring

    def base_label(self):
        r"""
        Always returns the same label.
        """
        return self._base_label

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        return self._finite

    def is_mutable(self):
        r"""
        Return if this surface is mutable.
        """
        return self._mutable

    def set_immutable(self):
        r"""
        Mark this surface as immutable.
        """
        self._mutable = False


    def make_immutable(self):
        r"""
        Mark this surface as immutable.
        """
        from sage.misc.superseded import deprecation
        deprecation(13109, "Do not use .make_immutable(). Use .set_immutable() instead.")
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

        This method updates both the edges of the polygon with label "label"
        and updates the edges listed in the glue_list.
        """
        self.__mutate()
        p=self.polygon(label)
        if p.num_edges() != len(glue_list):
            raise ValueEror("len(glue_list)="+str(len(glue_list))+\
                " and number of sides of polygon="+str(p.num_edges())+\
                " should be the same.")
        for e,(pp,ee) in enumerate(glue_list):
            self._change_edge_gluing(label, e, pp, ee)

    def add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Adds a the provided polygon to the surface. Utilizes gluing_list
        for the gluing data for edges (which must be a list
        of pairs representing edges of length equal to number of edges
        of the polygon).

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
        Remove the polygon with the provided label. Causes a ValueError
        if the base_label is removed.
        """
        if label==self._base_label:
            raise ValueError("Can not remove the base_label.")
        self.__mutate()
        return self._remove_polygon(label)

    def change_base_label(self, new_base_label):
        r"""
        Change the base_label to the provided label.
        """
        self.__mutate()
        self._base_label=new_base_label

    def __hash__(self):
        r"""
        Hash compatible with equals.
        """
        if hasattr(self,"_hash"):
            return self._hash
        if self.is_mutable():
            raise ValueError("Attempting to hash mutable surface.")
        if not self.is_finite():
            raise ValueError("Attempting to hash infinite surface.")
        h = 73+17*hash(self.base_ring())+23*hash(self.base_label())
        for pair in self.label_polygon_iterator():
            h = h + 7*hash(pair)
        for edgepair in self.edge_gluing_iterator():
            h = h + 3*hash(edgepair)
        self._hash=h
        return h


    def __eq__(self, other):
        r"""
        Implements a naive notion of equality where two finite surfaces are equal if:
        - their base rings are equal,
        - their base labels are equal,
        - their polygons are equal and labeled and glued in the same way.
        For infinite surfaces we use reference equality.
        Raises a value error if the surfaces are defined over different rings.
        """
        if self is other:
            return True
        if not isinstance(other, Surface):
            raise TypeError
        if not self.is_finite():
            if other.is_finite():
                return False
            else:
                raise ValueError("Can not compare infinite surfaces.")
        if self.base_ring() != other.base_ring():
            raise ValueError("Refusing to compare surfaces with different base rings.")
        if not self.is_mutable() and not other.is_mutable():
            hash1 = hash(self)
            hash2 = hash(other)
            if hash1 != hash2:
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
            for edge in range(polygon.num_edges()):
                if self.opposite_edge(label,edge) != other.opposite_edge(label,edge):
                    return False
        return True

    def __ne__(self, other):
        return not self == other

    def _test_base_ring(self, **options):
        # Test that the base_label is associated to a polygon
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        from sage.rings.ring import Field
        tester.assertTrue(isinstance(self.base_ring(), Field), \
            "base_ring="+str(self.base_ring())+" is not a Field.")

    def _test_base_label(self, **options):
        # Test that the base_label is associated to a polygon
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        from .polygon import ConvexPolygon
        tester.assertTrue(isinstance(self.polygon(self.base_label()), ConvexPolygon), \
            "polygon(base_label) does not return a ConvexPolygon. "+\
            "Here base_label="+str(self.base_label()))

    def _test_gluings(self, **options):
        # iterate over pairs with pair1 glued to pair2
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        if self.is_finite():
            it = self.edge_gluing_iterator()
        else:
            from itertools import islice
            it = islice(self.edge_gluing_iterator(), 30)

        for pair1,pair2 in it:
            tester.assertEqual(self.opposite_edge(pair2[0], pair2[1]), pair1,
                "edges not glued correctly:\n%s -> %s -> %s"%(pair1,pair2,self.opposite_edge(pair2[0], pair2[1])))

    def _test_override(self, **options):
        # Test that the required methods have been overridden and that some other methods have not been overridden.

        # Of course, we don't care if the methods are overridden or not we just want to warn the programmer.
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        # Build a naive Surface.
        from sage.rings.rational_field import QQ
        s=Surface(QQ,0,finite=True)

        # Check for override:
        tester.assertNotEqual(self.polygon.__func__,
                              s.polygon.__func__,
            "Method polygon of Surface must be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertNotEqual(self.opposite_edge.__func__, s.opposite_edge.__func__,
            "Method opposite_edge of Surface must be overridden. The Surface is of type "+str(type(self))+".")

        # Check not overridden:
        tester.assertEqual(self.base_ring.__func__, s.base_ring.__func__, \
            "Method base_ring of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.base_label.__func__, s.base_label.__func__, \
            "Method base_label of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.is_finite.__func__, s.is_finite.__func__, \
            "Method is_finite of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.is_mutable.__func__, s.is_mutable.__func__, \
            "Method is_mutable of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.set_immutable.__func__, s.set_immutable.__func__, \
            "Method set_immutable of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.walker.__func__, s.walker.__func__, \
            "Method walker of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.change_polygon.__func__, s.change_polygon.__func__, \
            "Method change_polygon of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.change_polygon_gluings.__func__, s.change_polygon_gluings.__func__, \
            "Method change_polygon_gluings of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.change_edge_gluing.__func__, s.change_edge_gluing.__func__, \
            "Method change_edge_gluing of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.add_polygon.__func__, s.add_polygon.__func__, \
            "Method add_polygon of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.remove_polygon.__func__, s.remove_polygon.__func__, \
            "Method remove_polygon of Surface should not be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertEqual(self.change_base_label.__func__, s.change_base_label.__func__, \
            "Method change_base_label of Surface should not be overridden. The Surface is of type "+str(type(self))+".")

        if self.is_mutable():
            # Check for override:
            tester.assertNotEqual(self._change_polygon.__func__,s._change_polygon.__func__,\
                "Method _change_polygon of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")
            tester.assertNotEqual(self._change_edge_gluing.__func__,s._change_edge_gluing.__func__,\
                "Method _change_edge_gluing of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")
            tester.assertNotEqual(self._add_polygon.__func__,s._add_polygon.__func__,"Method _add_polygon of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")
            tester.assertNotEqual(self._remove_polygon.__func__,s._remove_polygon.__func__,"Method _remove_polygon of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")

    def _test_polygons(self, **options):
        # Test that the base_label is associated to a polygon
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)
        from .polygon import ConvexPolygon
        if self.is_finite():
            it = self.label_iterator()
        else:
            from itertools import islice
            it = islice(self.label_iterator(), 30)
        for label in it:
            tester.assertTrue(isinstance(self.polygon(label), ConvexPolygon), \
                "polygon(label) does not return a ConvexPolygon when label="+str(label))

####
#### Surface_list
####

class Surface_list(Surface):
    r"""
    A fast mutable implementation of surface using a list to store polygons and gluings.

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.surface import Surface_list
        sage: p=polygons.regular_ngon(5)
        sage: s=Surface_list(base_ring=p.base_ring())
        sage: s.add_polygon(p) # gets label 0
        0
        sage: s.add_polygon( (-matrix.identity(2))*p ) # gets label 1
        1
        sage: s.change_polygon_gluings(0,[(1,e) for e in range(5)])
        sage: # base label defaults to zero.
        sage: s.set_immutable()
        sage: TestSuite(s).run()
    """
    ###
    ### Brief summary of internal workings.
    ###
    #
    # The Surface_list maintains a list self._p for storing polygons together
    # with gluing data.
    #
    # Here self._p[label] is typically a list of two elements
    # [polygon, gluing_list]. The gluing_list is then a list of pairs
    # (other_label, other_edge) so that typically
    # self.opposite_edge(label, edge) returns self._p[label][1][edge].
    #
    # If constructed with a surface parameter which is not None and copy=False,
    # then Surface_list stores a reference to the provided surface as
    # self._reference_surface. (Otherwise self._reference_surface is None).
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
            is false, then the resulting surface will not be mutable. If mutable
            is left at its default value of None, then the surface will be mutable
            if and only if a surface is not provided.

        EXAMPLES::

            sage: from flatsurf import *
            sage: print("We will surgically add a square into an infinite billiard surface")
            We will surgically add a square into an infinite billiard surface
            sage: p = polygons(vertices=[(0,0),(4,0),(0,3)])
            sage: s = similarity_surfaces.billiard(p)
            sage: ts=s.minimal_cover(cover_type="translation").copy(relabel=True, mutable=True)
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
        self._p = []
        self._reference_surface = None # Whether or not we store a reference surface
        self._removed_labels = []
        if surface is None:
            if base_ring is None:
                raise ValueError("Either surface or base_ring must be provided.")
            if not mutable is None and not mutable:
                raise ValueError("If no surface is provided, then mutable must be true.")
            self._num_polygons=0
            # default label is zero.
            Surface.__init__(self, base_ring, 0, finite=True, mutable=True)
        else:
            from .similarity_surface import SimilaritySurface
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
                # Temporarily set base_label to none. Update below.
                Surface.__init__(self, surface.base_ring(), None, finite=True, mutable=True)
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
                self.change_base_label(label_dict[surface.base_label()])
                if mutable is None or not mutable:
                    self.set_immutable()
            else:
                if surface.is_mutable():
                    raise ValueError("Surface_list will not store reference to a mutable surface.")
                self._reference_surface = surface
                self._ref_to_int={}
                self._int_to_ref=[]
                self.__get_label(surface.base_label())

                # Cache the base polygon
                polygon = surface.polygon(surface.base_label())
                self._p[0]=[polygon, [None for i in range(polygon.num_edges())]]

                self._num_polygons = self._reference_surface.num_polygons()
                # Set base label to zero.
                Surface.__init__(self, surface.base_ring(), 0, finite=surface.is_finite(), mutable=not mutable is None and mutable)

    def __get_label(self, ref_label):
        r"""
        Returns a corresponding label. Creates a new label if necessary.
        """
        try:
            return self._ref_to_int[ref_label]
        except KeyError:
            polygon = self._reference_surface.polygon(ref_label)
            data = [polygon,[None for i in range(polygon.num_edges())]]
            if len(self._removed_labels)>0:
                i = self._removed_labels.pop()
                self._p[i]=data
                self._ref_to_int[ref_label]=i
                self._int_to_ref[i]=ref_label
            else:
                i = len(self._p)
                if i!=len(self._int_to_ref):
                    raise RuntimeError("length of self._int_to_ref is "+\
                        str(len(self._int_to_ref))+" should be the same as "+\
                        "i="+str(i))
                self._p.append(data)
                self._ref_to_int[ref_label]=i
                self._int_to_ref.append(ref_label)
            return i

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
            # Here data is probably None...
            raise ValueError("Provided label was removed.")

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
            # Here data=None
            raise ValueError("Provided label was removed.")
        try:
            oe = glue[e]
        except KeyError:
            raise ValueError("Edge out of range of polygon.")
        if oe is None:
            if self._reference_surface is None:
                # Perhaps the user of this class left an edge unglued?
                raise RuntimeError("Obtained None as opposite_edge.")
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
            data[1]=[None for e in range(new_polygon.num_edges())]
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

        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.surface import Surface_list
            sage: p=polygons.regular_ngon(5)
            sage: s=Surface_list(base_ring=p.base_ring())
            sage: s.add_polygon(p, label=3)
            3
            sage: s.add_polygon( (-matrix.identity(2))*p, label=30)
            30
            sage: s.change_polygon_gluings(3,[(30,e) for e in range(5)])
            sage: s.change_base_label(30)
            sage: s.num_polygons()
            2
            sage: TestSuite(s).run()
            sage: s.remove_polygon(3)
            sage: s.add_polygon(p, label=6)
            6
            sage: s.change_polygon_gluings(6,[(30,e) for e in range(5)])
            sage: s.num_polygons()
            2
            sage: TestSuite(s).run()
            sage: s.change_base_label(6)
            sage: s.remove_polygon(30)
            sage: label = s.add_polygon((-matrix.identity(2))*p)
            sage: s.change_polygon_gluings(6,[(label,e) for e in range(5)])
            sage: TestSuite(s).run()
        """
        if new_polygon is None:
            data=[None,None]
        else:
            data=[new_polygon, [None for i in range(new_polygon.num_edges())] ]
        if label is None:
            if len(self._removed_labels)>0:
                new_label = self._removed_labels.pop()
                self._p[new_label]=data
            else:
                new_label = len(self._p)
                self._p.append(data)
                if not self._reference_surface is None:
                    # Need a blank in this list for algorithmic reasons
                    self._int_to_ref.append(None)
        else:
            new_label=int(label)
            if new_label<len(self._p):
                if not self._p[new_label] is None:
                    raise ValueError("Trying to add a polygon with label="+str(label)+" which already indexes a polygon.")
                self._p[new_label]=data
            else:
                if new_label-len(self._p)>100:
                    raise ValueError("Adding a polygon with label="+str(label)+" would add more than 100 entries in our list.")
                for i in range(len(self._p),new_label):
                    self._p.append(None)
                    self._removed_labels.append(i)
                    if not self._reference_surface is None:
                        # Need a blank in this list for algorithmic reasons
                        self._int_to_ref.append(None)

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
            for i in range(self.num_polygons()):
                yield i
        else:
            # We've removed some labels
            found=0
            i=0
            while found < self._num_polygons:
                if not self._p[i] is None:
                    found += 1
                    yield i
                i += 1

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        if label == len(self._p)-1:
            self._p.pop()
            if not self._reference_surface is None:
                ref_label = self._int_to_ref.pop()
                assert(len(self._int_to_ref)==label)
                if not ref_label is None:
                    del self._ref_to_int[ref_label]
        else:
            self._p[label]=None
            self._removed_labels.append(label)
            if not self._reference_surface is None:
                ref_label = self._int_to_ref[label]
                if not ref_label is None:
                    self._int_to_ref[label]=None
                    del self._ref_to_int[ref_label]
        self._num_polygons -= 1

def surface_list_from_polygons_and_gluings(polygons, gluings, mutable=False):
    r"""
    Take a list of polygons and gluings (given either as a list of pairs of edges, or as a dictionary),
    and produce a Surface_list from it. The mutable parameter determines the mutability of the resulting
    surface.
    """
    if not (isinstance(polygons,list) or isinstance(polygons,tuple)):
        raise ValueError("polygons must be a list or tuple.")
    field = polygons[0].parent().field()
    s=Surface_list(base_ring=field)
    for p in polygons:
        s.add_polygon(p)
    try:
        # dict case:
        it = iteritems(gluings)
    except AttributeError:
        # list case:
        it = gluings
    for (l1,e1),(l2,e2) in it:
        s.change_edge_gluing(l1,e1,l2,e2)
    if not mutable:
        s.set_immutable()
    return s

####
#### Surface_dict
####

class Surface_dict(Surface):
    r"""
    A mutable implementation of surface using a dictionary to store polygons and gluings. The dictionary
    implementation has the advantage that many label types are supported.

    The example below indicates how to construct a Surface_dict.

    EXAMPLES::

        sage: from flatsurf import *
        sage: p=polygons.regular_ngon(10)
        sage: s=Surface_dict(base_ring=p.base_ring())
        sage: s.add_polygon(p,label="A")
        'A'
        sage: s.change_polygon_gluings("A",[("A",(e+5)%10) for e in range(10)])
        sage: s.change_base_label("A")
        sage: s.set_immutable()
        sage: TestSuite(s).run()
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
    # polygon we set self._p[label]=None.
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
        if surface is None:
            if base_ring is None:
                raise ValueError("Either surface or base_ring must be provided.")
            if not mutable is None and not mutable:
                raise ValueError("If no surface is provided, then mutable must be true.")
            Surface.__init__(self, base_ring, None, finite=True, mutable=True)
        else:
            from .similarity_surface import SimilaritySurface
            if isinstance(surface,SimilaritySurface):
                surface=surface.underlying_surface()
            if not isinstance(surface,Surface):
                raise ValueError("surface must be either a Surface or SimilaritySurface")
            if not base_ring is None and base_ring != surface.base_ring():
                raise ValueError("You currently can not provide both a surface and a base_ring.")
            if copy==True:
                if not surface.is_finite():
                    raise ValueError("Can not copy an infinite surface.")
                for label,polygon in surface.label_polygon_iterator():
                    self._p[label]=[polygon, \
                        [ surface.opposite_edge(label,e) for e in range(polygon.num_edges()) ] ]
                # The only way we're mutable is if mutable=True:
                Surface.__init__(self, surface.base_ring(), surface.base_label(), finite=True, mutable=not mutable is None and mutable)
            else:
                if surface.is_mutable():
                    raise ValueError("Surface_dict will not store reference to a mutable surface.")
                self._reference_surface = surface
                # The only way we're mutable is if mutable=True:
                Surface.__init__(self, surface.base_ring(), surface.base_label(), \
                    finite=surface.is_finite(), mutable=not mutable is None and mutable)

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
                    [ self._reference_surface.opposite_edge(lab,e) for e in range(polygon.num_edges()) ] ]
                self._p[lab] = data
        if data is None:
            raise ValueError("Label "+str(label)+" was removed from the surface.")
        return data[0]

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
                    data=[new_polygon, \
                        [self._reference_surface.opposite_edge(label,e) for e in range(new_polygon.num_edges())] ]
                    self._p[label]=data
                else:
                    data=[new_polygon, [None for e in range(new_polygon.num_edges())] ]
                    self._p[label]=data
        if len(data[1]) != new_polygon.num_edges():
            data[1] = [None for e in range(new_polygon.num_edges())]
        if not gluing_list is None:
            self.change_polygon_gluings(label,gluing_list)

    def _change_edge_gluing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.
        """
        try:
            data = self._p[label1]
        except KeyError:
            if self._reference_surface is None:
                raise ValueError("No known polygon with provided label1="+str(label1))
            else:
                # Failure likely because reference_surface contains the polygon.
                # import the data into this surface.
                polygon = self._reference_surface.polygon(label1)
                data = [polygon, [self._reference_surface.opposite_edge(label1,e) for e in range(polygon.num_edges())]]
                self._p[label1]=data
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
            if self._reference_surface is None:
                raise ValueError("No known polygon with provided label2="+str(label2))
            else:
                # Failure likely because reference_surface contains the polygon.
                # import the data into this surface.
                polygon = self._reference_surface.polygon(label2)
                data = [polygon, [self._reference_surface.opposite_edge(label2,e) for e in range(polygon.num_edges())]]
                self._p[label2]=data
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
        data=[new_polygon, [None for i in range(new_polygon.num_edges())] ]
        if label is None:
            new_label = ExtraLabel()
        else:
            try:
                old_data = self._p[label]
                if old_data is None:
                    # already removed this polygon. That's good, we can add.
                    new_label=label
                else:
                    raise ValueError("We already have a polygon with label="+str(label))
            except KeyError:
                # This seems inconvienient to enforce:
                #
                #if not self._reference_surface is None:
                #    # Can not be sure we are not overwriting a polygon in the reference surface.
                #    raise ValueError("Can not assign this label to a Surface_dict containing a reference surface,"+\
                #        "which may already contain this label.")
                new_label = label
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

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        if self._reference_surface is None:
            try:
                data = self._p[label]
            except KeyError:
                raise ValueError("Label "+str(label)+" is not in the surface.")
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

class BaseRingChangedSurface(Surface):
    r"""
    A surface with a different base_ring.
    """
    def __init__(self, surface, ring):
        self._s=surface
        self._base_ring=ring
        from flatsurf.geometry.polygon import Polygons
        self._P=Polygons(self._base_ring)
        Surface.__init__(self, ring, self._s.base_label(), mutable=False, finite=self._s.is_finite())

    def polygon(self, lab):
        return self._P( self._s.polygon(lab) )

    def opposite_edge(self, p, e):
        return self._s.opposite_edge(p,e)


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

        def __next__(self):
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

        next = __next__ # for Python 2

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
        Return the "canonical" edge to walk through to get closer to the base_label,
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
            for i in range(limit):
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
            for e in range(polygon.num_edges()):
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

    def label_to_number(self, label, search=False, limit=100):
        r"""
        Return the number associated to the provided label.

        Returns an error if the label has not already been found by the walker
        unless search=True in which case we look for the label. We look by
        continuing to look for the label by walking over the surface visiting
        the next limit many polygons.
        """
        if not search:
            return self._label_dict[label]
        else:
            if label in self._label_dict:
                return self._label_dict[label]
            for i in range(limit):
                l = self.find_a_new_label()
                if label == l:
                    return self._label_dict[label]
            raise ValueError("Failed to find label even after searching. limit="+str(limit))

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

    def __init__(self, value=None):
        r"""
        Construct a new label.
        """
        if value is None:
            self._label = int(ExtraLabel._next)
            ExtraLabel._next = ExtraLabel._next + 1
        else:
            self._label = value

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

