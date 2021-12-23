# -*- coding: utf-8 -*-
r"""
Data structures for surfaces built from polygons.

All surfaces in sage-flatsurf are built from polygons whose sides are
identified by similarities. This module provides data structures to describe
such surfaces. Currently, there are two fundamental such data structures,
namely :class:`Surface_list` and `Surface_dict`. The former labels the polygons
that make up a surface by non-negative integers and the latter can use
arbitrary labels. Additionally, there are lots of other surface representations
that are not really implementing data structures but essentially just wrap
these two, e.g., a :class:`MinimalTranslationCover`.

All these surface implementations inherit from :class:`Surface` which describes
the contract that all surfaces must satisfy. As an absolute minimum, they
implement :meth:`Surface.polygon` which maps polygon labels to actual polygons,
and :meth:`Surface.opposite_edge` which describes the gluing of polygons.

EXAMPLES:

We built a torus by gluing the opposite sides of a square::

    sage: from flatsurf import polygons
    sage: from flatsurf.geometry.surface import Surface_list

    sage: S = Surface_list(QQ)
    sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
    0
    sage: S.set_edge_pairing(0, 0, 0, 2)
    sage: S.set_edge_pairing(0, 1, 0, 3)

    sage: S.polygon(0)
    Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
    sage: S.opposite_edge(0, 0)
    (0, 2)

There are two separate hierarchies of surfaces in sage-flatsurf. The underlying
data of a surface described by the subclasses of :class:`Surface` here and the
:class:`SimilaritySurface` and its subclasses which wrap a :class:`Surface`.
While a :class:`Surface` essentially provides the raw data of a surface, a
:class:`SimilaritySurface` then adds mathematical knowledge to that data
structure, e.g., by declaring that the data describes a
:class:`TranslationSurface`::

    sage: from flatsurf import TranslationSurface
    sage: T = TranslationSurface(S)

We can recover the underlying surface again::

    sage: T.underlying_surface() is S
    True

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Pat Hooper
#                      2019-2020 Vincent Delecroix
#                      2020-2021 Julian Rüth
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# ********************************************************************
from __future__ import absolute_import, print_function, division
from six.moves import range
from six import iteritems
from collections import deque

from sage.structure.sage_object import SageObject


class Surface(SageObject):
    r"""
    Abstract base class of all surfaces that are built from a set of polygons
    with edges identified. The identifications are compositions of homothety,
    rotations and translations, i.e., similarities that ensure that the surface
    is oriented.

    Concrete implementations of a surface must implement at least
    :meth:`polygon` and :meth:`opposite_edge`.

    To be able to modify a surface, subclasses should also implement
    :meth:`_change_polygon`, :meth:`_set_edge_pairing`, :meth:`_add_polygon`,
    :meth:`_remove_polygon`.

    For concrete implementations of a Surface, see, e.g., :class:`Surface_list`
    and :class:`Surface_dict`.

    INPUT:

    - ``base_ring`` -- the ring containing the coordinates of the vertices of
      the polygons

    - ``base_label`` -- the label of a chosen special polygon in the surface,
      see :meth:`base_label`

    - ``finite`` -- whether this is a finite surface, see :meth:`is_finite`

    - ``mutable`` -- whether this is a mutable surface; can be changed later
      with :meth:`set_immutable`

    EXAMPLES::

        sage: from flatsurf.geometry.surface import Surface, Surface_list, Surface_dict

        sage: S = Surface_list(QQ)
        sage: isinstance(S, Surface)
        True

        sage: S = Surface_dict(QQ)
        sage: isinstance(S, Surface)
        True

    """

    def __init__(self, base_ring, base_label, finite, mutable):
        if finite not in [False, True]:
            raise ValueError("finite must be either True or False")

        from sage.all import Rings
        if base_ring not in Rings():
            raise ValueError("base_ring must be a ring")

        if mutable not in [False, True]:
            raise ValueError("mutable must be either True or False")

        self._base_ring = base_ring
        self._base_label = base_label
        self._finite = finite
        self._mutable = mutable

        self._cache = {}

    def is_triangulated(self, limit=None):
        r"""
        EXAMPLES::

            sage: import flatsurf
            sage: G = SymmetricGroup(4)
            sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: S.is_triangulated()
            False
            sage: S.triangulate().is_triangulated()
            True
        """
        it = self.label_iterator()
        if not self.is_finite():
            if limit is None:
                raise ValueError("for infinite polygon, 'limit' must be set to a positive integer")
            else:
                from itertools import islice
                it = islice(it, limit)
        for p in it:
            if self.polygon(p).num_edges() != 3:
                return False
        return True

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

    def _change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Internal method used by change_polygon(). Should not be called directly.

        Mutable surfaces should implement this method.
        """
        raise NotImplementedError

    def _set_edge_pairing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.

        Mutable surfaces should implement this method.
        """
        raise NotImplementedError

    def _add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.

        Mutable surfaces should implement this method.
        """
        raise NotImplementedError

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.

        Mutable surfaces should implement this method.
        """
        raise NotImplementedError

    def num_polygons(self):
        r"""
        Return the number of polygons making up the surface, or
        sage.rings.infinity.Infinity if the surface is infinite.

        This is a generic method. On a finite surface it will be linear time in
        the edges the first time it is run, then constant time (assuming no
        mutation occurs).

        Subclasses should consider overriding this method for increased
        performance.
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

        Subclasses should consider overriding this method for increased
        performance.
        """
        return iter(self.walker())

    def label_polygon_iterator(self):
        r"""
        Iterate over pairs (label, polygon).

        Subclasses should consider overriding this method for increased
        performance.
        """
        for label in self.label_iterator():
            yield label, self.polygon(label)

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

    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        return self._base_ring

    def base_label(self):
        r"""
        Return the label of a special chosen polygon in the surface.

        When no specific choice was made with :meth:`set_base_label`, this
        might just be the initial polygon, i.e., the one that was first added
        in the construction of this surface.

        This label is for example used as the starting position when walking
        the surface in a canonical order with :meth:`walker`.

        EXAMPLES::

            sage: import flatsurf
            sage: G = SymmetricGroup(4)
            sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: S.base_label()
            1

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

    def set_edge_pairing(self, label1, edge1, label2, edge2):
        r"""
        Updates the gluing so that (label,edge1) is glued to (label2, edge2).
        """
        self.__mutate()
        self._set_edge_pairing(label1, edge1, label2, edge2)

    # TODO: deprecation alias?
    change_edge_gluing = set_edge_pairing

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
            raise ValueError("len(glue_list)="+str(len(glue_list))+\
                " and number of sides of polygon="+str(p.num_edges())+\
                " should be the same.")
        for e,(pp,ee) in enumerate(glue_list):
            self._set_edge_pairing(label, e, pp, ee)

    def add_polygons(self, polygons):
        return [self.add_polygon(p) for p in polygons]

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
        return self._add_polygon(new_polygon, gluing_list, label)

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
        tester = self._tester(**options)
        from sage.all import Rings
        tester.assertTrue(self.base_ring() in Rings())

    def _test_base_label(self, **options):
        # Test that the base_label is associated to a polygon
        tester = self._tester(**options)
        from .polygon import ConvexPolygon
        tester.assertTrue(isinstance(self.polygon(self.base_label()), ConvexPolygon), \
            "polygon(base_label) does not return a ConvexPolygon. "+\
            "Here base_label="+str(self.base_label()))

    def _test_gluings(self, **options):
        # iterate over pairs with pair1 glued to pair2
        tester = self._tester(**options)

        if self.is_finite():
            it = self.label_iterator()
        else:
            from itertools import islice
            it = islice(self.label_iterator(), 30)

        for lab in it:
            p = self.polygon(lab)
            for k in range(p.num_edges()):
                e = (lab, k)
                f = self.opposite_edge(lab, k)
                tester.assertFalse(f is None,
                        "edge ({}, {}) is not glued".format(lab, k))
                g = self.opposite_edge(f[0], f[1])
                tester.assertEqual(e, g,
                        "edge gluing is not a pairing:\n{} -> {} -> {}".format(e, f, g))

    def _test_override(self, **options):
        # Test that the required methods have been overridden and that some other methods have not been overridden.

        # Of course, we don't care if the methods are overridden or not we just want to warn the programmer.
        if 'tester' in options:
            tester = options['tester']
        else:
            tester = self._tester(**options)

        # Check for override:
        tester.assertNotEqual(self.polygon.__func__,
                              Surface.polygon,
            "Method polygon of Surface must be overridden. The Surface is of type "+str(type(self))+".")
        tester.assertNotEqual(self.opposite_edge.__func__, Surface.opposite_edge,
            "Method opposite_edge of Surface must be overridden. The Surface is of type "+str(type(self))+".")

        if self.is_mutable():
            # Check for override:
            tester.assertNotEqual(self._change_polygon.__func__, Surface._change_polygon,\
                "Method _change_polygon of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")
            tester.assertNotEqual(self._set_edge_pairing.__func__, Surface._set_edge_pairing,\
                "Method _set_edge_pairing of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")
            tester.assertNotEqual(self._add_polygon.__func__, Surface._add_polygon, "Method _add_polygon of Surface must be overridden in a mutable surface. "+\
                "The Surface is of type "+str(type(self))+".")
            tester.assertNotEqual(self._remove_polygon.__func__, Surface._remove_polygon, "Method _remove_polygon of Surface must be overridden in a mutable surface. "+\
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


class Surface_list(Surface):
    r"""
    A fast mutable :class:`Surface` using a list to store polygons and gluings.

    ALGORITHM:

    Internally, we maintain a list ``_p`` for storing polygons together with
    gluing data.

    Each ``_p[label]`` is typically a pair ``(polygon, gluing_list)`` where
    ``gluing_list`` is a list of pairs ``(other_label, other_edge)`` such that
    :meth:`opposite_edge(label, edge)` returns ``_p[label][1][edge]``.

    INPUT:

    - ``base_ring`` -- ring or ``None`` (default: ``None``); the ring
      containing the coordinates of the vertices of the polygons. If ``None``,
      the :meth:`base_ring` will be the one of ``surface``.

    - ``surface`` -- :class:`Surface`, :class:`SimilaritySurface`, or
      ``None`` (default: ``None``); a surface to be copied or referenced (see
      ``copy``). If ``None``, the surface is initially empty.

    - ``copy`` -- boolean or ``None`` (default: ``None``); whether the data
      underlying ``surface`` is copied into this surface or just a reference to
      that surface is kept. If ``None``, a sensible default is chosen, namely
      ``surface.is_mutable()``.

    - ``mutable`` -- boolean or ``None`` (default: ``None``); whether this
      surface is mutable. When ``None``, the surface will be mutable iff
      ``surface`` is ``None``.

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

    We surgically add a square into an infinite billiard surface::

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

    def __init__(self, base_ring=None, surface=None, copy=None, mutable=None):
        self._p = []  # list of pairs (polygon, gluings)
        self._reference_surface = None
        self._removed_labels = []
        self._num_polygons = 0

        # Validate input parameters and fill in defaults
        base_ring, surface, copy, mutable, finite = Surface_list._validate_init_parameters(base_ring=base_ring, surface=surface, copy=copy, mutable=mutable, finite=None)

        Surface.__init__(self, base_ring, base_label=0, finite=finite, mutable=True)

        # Initialize surface from reference surface
        if surface is not None:
            if copy is True:
                reference_label_to_label = {
                    label: self.add_polygon(polygon)
                    for label, polygon in surface.label_polygon_iterator()
                }

                for ((label, edge), (glued_label, glued_edge)) in surface.edge_gluing_iterator():
                    self.set_edge_pairing(reference_label_to_label[label], edge, reference_label_to_label[glued_label], glued_edge)

                self.change_base_label(reference_label_to_label[surface.base_label()])
            else:
                self._reference_surface = surface
                self._ref_to_int = {}
                self._int_to_ref = []

                self._num_polygons = surface.num_polygons()

                # Cache the base polygon
                self.change_base_label(self.__get_label(surface.base_label()))
                assert self.base_label() == 0

        if not mutable:
            self.set_immutable()

    @classmethod
    def _validate_init_parameters(cls, base_ring, surface, copy, mutable, finite):
        r"""
        Helper method for ``__init__`` that validates the parameters and
        returns them in the same order with defaults filled in.
        """
        if surface is None and base_ring is None:
            raise ValueError("Either surface or base_ring must be provided.")

        if surface is None:
            if copy is not None:
                raise ValueError("Cannot copy when surface was provided.")

            if mutable is None:
                mutable = True

            finite = True
        else:
            from .similarity_surface import SimilaritySurface
            if isinstance(surface, SimilaritySurface):
                surface = surface.underlying_surface()

            if not isinstance(surface, Surface):
                raise TypeError("surface must be a Surface or a SimilaritySurface")

            if not surface.is_finite() and surface.is_mutable():
                raise NotImplementedError("Cannot create surface from infinite mutable surface yet.")

            if base_ring is None:
                base_ring = surface.base_ring()

            if base_ring != surface.base_ring():
                raise NotImplementedError("Cannot provide both a surface and a base_ring yet.")

            if mutable is None:
                mutable = True

            if copy is None:
                copy = surface.is_mutable()

            if copy and not surface.is_finite():
                raise ValueError("Cannot copy infinite surface.")

            if surface.is_mutable() and not copy:
                raise ValueError("Cannot reference mutable surface.")

            finite = surface.is_finite()

        return base_ring, surface, copy, mutable, finite

    def __get_label(self, ref_label):
        r"""
        Returns a corresponding label. Creates a new label if necessary.
        """
        try:
            return self._ref_to_int[ref_label]
        except KeyError:
            polygon = self._reference_surface.polygon(ref_label)
            data = [polygon, [None for i in range(polygon.num_edges())]]
            if len(self._removed_labels) > 0:
                i = self._removed_labels.pop()
                self._p[i] = data
                self._ref_to_int[ref_label] = i
                self._int_to_ref[i] = ref_label
            else:
                i = len(self._p)
                if i != len(self._int_to_ref):
                    raise RuntimeError("length of self._int_to_ref is " + str(len(self._int_to_ref))+" should be the same as i=" + str(i))
                self._p.append(data)
                self._ref_to_int[ref_label] = i
                self._int_to_ref.append(ref_label)
            return i

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        try:
            data = self._p[lab]
        except IndexError:
            if self._reference_surface is None:
                raise ValueError(f"No polygon with label {lab}.")

            for label in self.label_iterator():
                if label >= lab:
                    break

            if lab >= len(self._p):
                raise ValueError(f"no polygon with label {lab}")

            data = self._p[lab]

        if data is None:
            raise ValueError("Provided label was removed.")

        return data[0]

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        try:
            data = self._p[p]
        except KeyError:
             raise ValueError("No known polygon with provided label")
        if data is None:
            raise ValueError("Provided label was removed.")
        glue = data[1]
        try:
            oe = glue[e]
        except KeyError:
            raise ValueError("Edge out of range of polygon.")
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

    def _set_edge_pairing(self, label1, edge1, label2, edge2):
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

    # TODO: deprecation alias?
    _change_edge_gluing = _set_edge_pairing

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

    def ramified_cover(self, d, data):
        r"""
        Make a ramified cover of this surface.

        INPUT:

        - ``d`` - integer (the degree of the cover)

        - ``data`` - a dictionary which to a pair ``(label, edge_num)`` associates a permutation
          of {1,...,d}
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        G = SymmetricGroup(d)
        for k in data:
            data[k] = G(data[k])
        cover = Surface_list(base_ring=self.base_ring())
        labels = list(self.label_iterator())
        edges = set(self.edge_iterator())
        cover_labels = {}
        for i in range(1,d+1):
            for lab in self.label_iterator():
                cover_labels[(lab, i)] = cover.add_polygon(self.polygon(lab))
        while edges:
            lab, e = elab = edges.pop()
            llab, ee = eelab = self.opposite_edge(lab, e)
            edges.remove(eelab)
            if elab in data:
                if eelab in data:
                    if not (data[elab] * data[eelab]).is_one():
                        raise ValueError("inconsistent covering data")
                s = data[elab]
            elif eelab in data:
                s = ~data[eelab]
            else:
                s = G.one()

            for i in range(1, d+1):
                p0 = cover_labels[(lab, i)]
                p1 = cover_labels[(lab, s(i))]
                cover.set_edge_pairing(p0, e, p1, ee)
        return cover


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


class Surface_dict(Surface):
    r"""
    A mutable :class:`Surface` using a dict to store polygons and gluings.

    Unlike :class:`Surface_list`, this surface is not limited to integer
    labels. However, :class:`Surface_list` is likely more efficient for most
    applications.

    ALGORITHM:

    Internally, we maintain a dict ``_p`` for storing polygons together with
    gluing data.

    Each ``_p[label]`` is typically a pair ``(polygon, gluing_dict)`` where
    ``gluing_dict`` is maps ``other_label`` to ``other_edge`` such that
    :meth:`opposite_edge(label, edge)` returns ``_p[label][1][edge]``.

    INPUT:

    - ``base_ring`` -- ring or ``None`` (default: ``None``); the ring
      containing the coordinates of the vertices of the polygons. If ``None``,
      the :meth:`base_ring` will be the one of ``surface``.

    - ``surface`` -- :class:`Surface`, :class:`SimilaritySurface`, or
      ``None`` (default: ``None``); a surface to be copied or referenced (see
      ``copy``). If ``None``, the surface is initially empty.

    - ``copy`` -- boolean or ``None`` (default: ``None``); whether the data
      underlying ``surface`` is copied into this surface or just a reference to
      that surface is kept. If ``None``, a sensible default is chosen, namely
      ``surface.is_mutable()``.

    - ``mutable`` -- boolean or ``None`` (default: ``None``); whether this
      surface is mutable. When ``None``, the surface will be mutable iff
      ``surface`` is ``None``.

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

    def __init__(self, base_ring=None, surface=None, copy=None, mutable=None):
        self._p = {}
        self._reference_surface = None

        # Validate input parameters and fill in defaults
        base_ring, surface, copy, mutable, finite = Surface_list._validate_init_parameters(base_ring=base_ring, surface=surface, copy=copy, mutable=mutable, finite=None)

        Surface.__init__(self, base_ring, base_label=None, finite=finite, mutable=True)

        # Initialize surface from reference surface
        if surface is not None:
            if copy is True:
                reference_label_to_label = {
                    label: self.add_polygon(polygon, label=label)
                    for label, polygon in surface.label_polygon_iterator()
                }

                for ((label, edge), (glued_label, glued_edge)) in surface.edge_gluing_iterator():
                    self.set_edge_pairing(reference_label_to_label[label], edge, reference_label_to_label[glued_label], glued_edge)

                self.change_base_label(reference_label_to_label[surface.base_label()])
            else:
                self._reference_surface = surface

                self.change_base_label(surface.base_label())

        if not mutable:
            self.set_immutable()

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        try:
            data = self._p[lab]
        except KeyError:
            if self._reference_surface is None:
                raise ValueError(f"No polygon with label {lab}.")

            polygon = self._reference_surface.polygon(lab)
            data = [self._reference_surface.polygon(lab),
                    [self._reference_surface.opposite_edge(lab, e) for e in range(polygon.num_edges())]]
            self._p[lab] = data

        if data is None:
            raise ValueError("Label "+str(lab)+" was removed from the surface.")
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
            raise ValueError("Label "+str(p)+" was removed from the surface.")
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

    def _set_edge_pairing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by set_edge_pairing(). Should not be called directly.
        """
        try:
            data = self._p[label1]
        except KeyError:
            if self._reference_surface is None:
                raise ValueError("No known polygon with provided label1 = {}".format(label1))
            else:
                # Failure likely because reference_surface contains the polygon.
                # import the data into this surface.
                polygon = self._reference_surface.polygon(label1)
                data = [polygon, [self._reference_surface.opposite_edge(label1, e) for e in range(polygon.num_edges())]]
                self._p[label1] = data
        try:
            data[1][edge1] = (label2,edge2)
        except IndexError:
            # break down error
            if data is None:
                raise ValueError("polygon with label1={} was removed".format(label1))
            data1 = data[1]
            try:
                data1[edge1] = (label2, edge2)
            except IndexError:
                raise ValueError("edge1={} is out of range in polygon with label1={}".format(edge1, label1))
        try:
            data = self._p[label2]
        except KeyError:
            if self._reference_surface is None:
                raise ValueError("no polygon with label2={}".format(label2))
            else:
                # Failure likely because reference_surface contains the polygon.
                # import the data into this surface.
                polygon = self._reference_surface.polygon(label2)
                data = [polygon, [self._reference_surface.opposite_edge(label2,e) for e in range(polygon.num_edges())]]
                self._p[label2]=data
        try:
            data[1][edge2] = (label1,edge1)
        except IndexError:
            # break down error
            if data is None:
                raise ValueError("polygon with label1={} was removed".format(label1))
            data1 = data[1]
            try:
                data1[edge2] = (label1, edge1)
            except IndexError:
                raise ValueError("edge {} is out of range in polygon with label2={}".format(edge2, label2))

    _change_edge_gluing = _set_edge_pairing

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
                    raise ValueError("label={} already used by another polygon".format(label))
            except KeyError:
                # This seems inconvienient to enforce:
                #
                #if not self._reference_surface is None:
                #    # Can not be sure we are not overwriting a polygon in the reference surface.
                #    raise ValueError("Can not assign this label to a Surface_dict containing a reference surface,"+\
                #        "which may already contain this label.")
                new_label = label
        self._p[new_label]=data
        if gluing_list is not None:
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
        from flatsurf.geometry.polygon import ConvexPolygons
        self._P=ConvexPolygons(self._base_ring)
        Surface.__init__(self, ring, self._s.base_label(), mutable=False, finite=self._s.is_finite())

    def polygon(self, lab):
        return self._P( self._s.polygon(lab) )

    def opposite_edge(self, p, e):
        return self._s.opposite_edge(p,e)


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


class LabelComparator(object):
    r"""
    Implements a total ordering on labels, which may be of varying types.
    
    We use hashes, so if hash(label1) < hash(label2) we declare label1 < label2.
    For objects with the same hash, we store an arbitary ordering.
    """
    def __init__(self):
        r"""
        Initialize a label comparator.
        """
        self._hash_collision_resolver = {}

    def _get_resolver_index(self, label_hash, label):
        try:
            lst = self._hash_collision_resolver[label_hash]
        except KeyError:
            lst = []
            self._hash_collision_resolver[label_hash] = lst
        for i, stored_label in enumerate(lst):
            if label == stored_label:
                return i
        # At this point we know label is not in lst
        lst.append(label)
        return len(lst)-1
            
    def lt(self, l1, l2):
        r"""
        Return the truth value of l1 < l2.
        """
        h1 = hash(l1)
        h2 = hash(l2)
        if h1 < h2:
            return True
        if h1 > h2:
            return False
        # Otherwise the hashes are equal.
        if l1 == l2:
            return False
        return self._get_resolver_index(h1, l1) < self._get_resolver_index(h1, l2)
    
    def le(self, l1, l2):
        r"""
        Return the truth value of l1 <= l2.
        """
        return self.lt(l1, l2) or l1 == l2
    
    def gt(self, l1, l2):
        r"""
        Return the truth value of l1 > l2.
        """
        return self.lt(l2, l1)
    
    def ge(self, l1, l2):
        r"""
        Return the truth value of l1 >= l2.
        r"""
        return self.lt(l2, l1) or l1 == l2
