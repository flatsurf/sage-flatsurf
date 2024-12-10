r"""
Legacy data structures for surfaces built from polygons.

All surfaces in sage-flatsurf are built from polygons whose sides are
identified by similarities. This module provides deprecated data structures to
describe such surfaces. Currently, there are two fundamental such data
structures, namely :class:`Surface_list` and `Surface_dict`. The former labels
the polygons that make up a surface by non-negative integers and the latter can
use arbitrary labels.

In principle there is nothing wrong with this approach. However, the
implementation is plagued by feature creep so we tried to clean things up in
2023 and reimplemented these in :mod:`flatsurf.geometry.surface`.

All the surfaces here inherit from :class:`Surface` which describes the
contract that surfaces used to satisfy. As an absolute minimum, they implement
:meth:`Surface.polygon` which maps polygon labels to actual polygons, and
:meth:`Surface.opposite_edge` which describes the gluing of polygons.

EXAMPLES:

We built a torus by gluing the opposite sides of a square::

    sage: from flatsurf import Polygon
    sage: from flatsurf.geometry.surface import Surface_list

    sage: S = Surface_list(QQ)
    doctest:warning
    ...
    UserWarning: Surface_list has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
    sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
    0
    sage: S.set_edge_pairing(0, 0, 0, 2)
    sage: S.set_edge_pairing(0, 1, 0, 3)

    sage: S.polygon(0)
    Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
    sage: S.opposite_edge(0, 0)
    (0, 2)

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 W. Patrick Hooper
#                      2019-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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
import collections.abc

from sage.misc.cachefunc import cached_method

from flatsurf.geometry.surface import OrientedSimilaritySurface


class Surface(OrientedSimilaritySurface):
    r"""
    Abstract base class of all surfaces that are built from a set of polygons
    with edges identified. The identifications are compositions of homothety,
    rotations and translations, i.e., similarities that ensure that the surface
    is oriented.

    Concrete implementations of a surface must implement at least
    :meth:`polygon` and :meth:`opposite_edge`.

    To be able to modify a surface, subclasses should also implement
    ``_change_polygon``, ``_set_edge_pairing``, ``_add_polygon``,
    ``_remove_polygon``.

    For concrete implementations of a Surface, see, e.g., :class:`Surface_list`
    and :class:`Surface_dict`.

    INPUT:

    - ``base_ring`` -- the ring containing the coordinates of the vertices of
      the polygons

    - ``base_label`` -- the label of a chosen special polygon in the surface,
      see :meth:`.base_label`

    - ``finite`` -- whether this is a finite surface, see :meth:`is_finite_type`

    - ``mutable`` -- whether this is a mutable surface; can be changed later
      with :meth:`.set_immutable`

    EXAMPLES::

        sage: from flatsurf.geometry.surface import Surface, Surface_list, Surface_dict

        sage: S = Surface_list(QQ)
        sage: isinstance(S, Surface)
        True

        sage: S = Surface_dict(QQ)
        doctest:warning
        ...
        UserWarning: Surface_dict has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
        sage: isinstance(S, Surface)
        True

    TESTS:

    Users are being warned if they try to define a surface over an inexact ring::

        sage: S = Surface_list(RR)
        ...
        UserWarning: surface defined over an inexact ring; many operations in sage-flatsurf are not going to work correctly over this ring
        sage: isinstance(S, Surface)
        True

    """

    def __init__(
        self,
        base_ring,
        base_label,
        finite,
        mutable,
        category=None,
        deprecation_warning=True,
    ):
        if deprecation_warning:
            import warnings

            warnings.warn(
                "base class Surface has been deprecated and will be removed in a future version of sage-flatsurf; use OrientedSimilaritySurface instead"
            )

        if finite not in [False, True]:
            raise ValueError("finite must be either True or False")

        from sage.all import Rings

        if base_ring not in Rings():
            raise ValueError("base_ring must be a ring")

        if not base_ring.is_exact():
            from warnings import warn

            warn(
                "surface defined over an inexact ring; many operations in sage-flatsurf are not going to work correctly over this ring"
            )

        if mutable not in [False, True]:
            raise ValueError("mutable must be either True or False")

        self._base_label = base_label
        self._finite = finite
        self._mutable = mutable

        self._cache = {}

        from flatsurf.geometry.categories import SimilaritySurfaces

        if category is None:
            category = SimilaritySurfaces()

        # Previously, all surfaces were assumed to be connected and without
        # boundary (even though it was possible to construct non-connected
        # surfaces but only the base_label component was really functional
        # then).
        category &= SimilaritySurfaces().Oriented().WithoutBoundary().Connected()

        if finite:
            category = category.FiniteType()

        OrientedSimilaritySurface.__init__(self, base=base_ring, category=category)

        if not mutable:
            self._refine_category_(self.refined_category())

    def _test_refined_category(self, **options):
        if self.is_mutable():
            return

        super()._test_refined_category(**options)

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf.geometry.surface import Surface, Surface_list, Surface_dict
            sage: S = Surface_list(QQ)
            sage: S
            Surface built from 0 polygons

        """
        if not self.is_finite_type():
            return "Surface built from infinitely many polygons"
        if self.num_polygons() == 1:
            return "Surface built from 1 polygon"

        return "Surface built from {} polygons".format(self.num_polygons())

    def labels(self):
        return LabelsView(self)

    def is_triangulated(self, limit=None):
        r"""
        EXAMPLES::

            sage: import flatsurf
            sage: G = SymmetricGroup(4)
            sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: S.is_triangulated()
            False
            sage: S.triangulate().codomain().is_triangulated()
            True
        """
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
            )

        it = self.label_iterator()
        if not self.is_finite_type():
            if limit is None:
                raise ValueError(
                    "for infinite polygon, 'limit' must be set to a positive integer"
                )
            else:
                from itertools import islice

                it = islice(it, limit)
        for p in it:
            if len(self.polygon(p).vertices()) != 3:
                return False
        return True

    def polygon(self, label):
        r"""
        Return the polygon with the provided label.

        This method must be overridden in subclasses.
        """
        raise NotImplementedError

    def opposite_edge(self, label, e=None):
        r"""
        Given the label ``label`` of a polygon and an edge ``e`` in that
        polygon returns the pair (``ll``, ``ee``) to which this edge is glued.

        This method must be overridden in subclasses.
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
        if self.is_finite_type():
            lw = self.walker()
            lw.find_all_labels()
            return len(lw)
        else:
            from sage.rings.infinity import Infinity

            return Infinity

    def label_iterator(self, polygons=False):
        r"""
        Iterator over all polygon labels.

        Subclasses should consider overriding this method for increased
        performance.
        """
        if polygons:
            return zip(self.labels(), self.polygons())
        return iter(self.walker())

    def roots(self):
        if self._base_label is None:
            return ()
        return (self._base_label,)

    def is_finite_type(self):
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

        self._refine_category_(self.refined_category())

    def walker(self):
        r"""
        Return a LabelWalker which walks over the surface in a canonical way.
        """
        try:
            return self._cache["lw"]
        except KeyError:
            lw = LabelWalker(self, deprecation_warning=False)
            self._cache["lw"] = lw
            return lw

    def __mutate(self):
        r"""
        Called before a mutation occurs. Do not call directly.
        """
        if not self.is_mutable():
            raise Exception("surface must be mutable")

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
        if not (gluing_list is None or len(new_polygon.vertices()) == len(gluing_list)):
            raise ValueError
        self._change_polygon(label, new_polygon, gluing_list)

    def set_edge_pairing(self, label1, edge1, label2, edge2):
        r"""
        Update the gluing so that (``label1``, ``edge1``) is glued to
        (``label2``, ``edge2``).
        """
        self.__mutate()
        self._set_edge_pairing(label1, edge1, label2, edge2)

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
        p = self.polygon(label)
        if len(p.vertices()) != len(glue_list):
            raise ValueError(
                "len(glue_list)="
                + str(len(glue_list))
                + " and number of sides of polygon="
                + str(len(p.vertices()))
                + " should be the same."
            )
        for e, (pp, ee) in enumerate(glue_list):
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
        if not (gluing_list is None or len(new_polygon.vertices()) == len(gluing_list)):
            raise ValueError
        label = self._add_polygon(new_polygon, gluing_list, label)
        if self._base_label is None:
            self.change_base_label(label)
        return label

    def remove_polygon(self, label):
        r"""
        Remove the polygon with the provided label. Causes a ValueError
        if the base_label is removed.
        """
        if label == self._base_label:
            raise ValueError("Can not remove the base_label.")
        self.__mutate()
        return self._remove_polygon(label)

    def change_base_label(self, new_base_label):
        r"""
        Change the base_label to the provided label.
        """
        self.__mutate()
        self._base_label = new_base_label

    @cached_method
    def __hash__(self):
        r"""
        Hash compatible with equals.
        """
        if self.is_mutable():
            raise TypeError("mutable surface is not hashable")

        if not self.is_finite_type():
            raise TypeError("cannot hash this infinite surface")

        return hash(
            (
                self.base_ring(),
                self.root(),
                tuple(zip(self.labels(), self.polygons())),
                tuple(self.gluings()),
            )
        )

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.surface import Surface_dict
            sage: from flatsurf.geometry.polygon import Polygon, ConvexPolygons

            sage: S = Surface_dict(QQ)
            sage: P = ConvexPolygons(QQ)
            doctest:warning
            ...
            UserWarning: ConvexPolygons() has been deprecated and will be removed from a future version of sage-flatsurf; use Polygon() to create polygons.
            If you really need the category of convex polygons over a ring use EuclideanPolygons(ring).Simple().Convex() instead.
            sage: S.add_polygon(P([(1, 0), (0, 1), (-1, -1)]), label=0)
            doctest:warning
            ...
            UserWarning: ConvexPolygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead
            0
            sage: S == S
            True

            sage: T = Surface_dict(QQ)
            sage: S == T
            False

        TESTS::

            sage: S == 42
            False

        """
        if self is other:
            return True

        if not isinstance(other, Surface):
            return False

        if self.is_mutable() != other.is_mutable():
            return False

        if not self._eq_oriented_similarity_surfaces(other):
            return False

        if self.num_polygons() == 0:
            # Only compare base labels when the surfaces are not empty.
            return True

        if self.root() != other.root():
            return False

        return True

    def __ne__(self, other):
        return not self == other

    def _eq_oriented_similarity_surfaces(self, other):
        # Whether this surface equals other in terms of oriented similarity surfaces
        if self is other:
            return True

        if not isinstance(other, OrientedSimilaritySurface):
            return False

        if self.base_ring() != other.base_ring():
            return False

        if self.category() != other.category():
            return False

        if self.is_finite_type() != other.is_finite_type():
            return False

        if self.is_finite_type():
            if len(self.polygons()) == 0:
                return len(other.polygons()) == 0
            if len(other.polygons()) == 0:
                return False

        if self.roots() != other.roots():
            return False

        for label in self.roots():
            if self.polygon(label) != other.polygon(label):
                return False

        if not self.is_finite_type():
            raise NotImplementedError("cannot compare these infinite surfaces yet")

        if len(self.polygons()) != len(other.polygons()):
            return False

        for label, polygon in zip(self.labels(), self.polygons()):
            try:
                polygon2 = other.polygon(label)
            except ValueError:
                return False
            if polygon != polygon2:
                return False
            for edge in range(len(polygon.vertices())):
                if self.opposite_edge(label, edge) != other.opposite_edge(label, edge):
                    return False

        return True

    def _test_base_ring(self, **options):
        # Test that the base_label is associated to a polygon
        tester = self._tester(**options)
        from sage.all import Rings

        tester.assertTrue(self.base_ring() in Rings())

    def _test_base_label(self, **options):
        # Test that the base_label is associated to a polygon
        tester = self._tester(**options)

        tester.assertTrue(
            self.polygon(self.root()).is_convex(),
            "polygon(base_label) does not return a ConvexPolygon. "
            + "Here base_label="
            + str(self.root()),
        )

    def _test_override(self, **options):
        # Test that the required methods have been overridden and that some other methods have not been overridden.

        # Of course, we don't care if the methods are overridden or not we just want to warn the programmer.
        if "tester" in options:
            tester = options["tester"]
        else:
            tester = self._tester(**options)

        # Check for override:
        tester.assertNotEqual(
            self.polygon.__func__,
            Surface.polygon,
            "Method polygon of Surface must be overridden. The Surface is of type "
            + str(type(self))
            + ".",
        )
        tester.assertNotEqual(
            self.opposite_edge.__func__,
            Surface.opposite_edge,
            "Method opposite_edge of Surface must be overridden. The Surface is of type "
            + str(type(self))
            + ".",
        )

        if self.is_mutable():
            # Check for override:
            tester.assertNotEqual(
                self._change_polygon.__func__,
                Surface._change_polygon,
                "Method _change_polygon of Surface must be overridden in a mutable surface. "
                + "The Surface is of type "
                + str(type(self))
                + ".",
            )
            tester.assertNotEqual(
                self._set_edge_pairing.__func__,
                Surface._set_edge_pairing,
                "Method _set_edge_pairing of Surface must be overridden in a mutable surface. "
                + "The Surface is of type "
                + str(type(self))
                + ".",
            )
            tester.assertNotEqual(
                self._add_polygon.__func__,
                Surface._add_polygon,
                "Method _add_polygon of Surface must be overridden in a mutable surface. "
                + "The Surface is of type "
                + str(type(self))
                + ".",
            )
            tester.assertNotEqual(
                self._remove_polygon.__func__,
                Surface._remove_polygon,
                "Method _remove_polygon of Surface must be overridden in a mutable surface. "
                + "The Surface is of type "
                + str(type(self))
                + ".",
            )

    def _test_polygons(self, **options):
        # Test that the base_label is associated to a polygon
        if "tester" in options:
            tester = options["tester"]
        else:
            tester = self._tester(**options)

        if self.is_finite_type():
            it = self.label_iterator()
        else:
            from itertools import islice

            it = islice(self.label_iterator(), 30)
        for label in it:
            tester.assertTrue(
                self.polygon(label).is_convex(),
                "polygon(label) does not return a ConvexPolygon when label="
                + str(label),
            )

    def point(self, label, position, limit=None, ring=None):
        r"""
        Return the :class:`flatsurf.geometry.surface_objects.SurfacePoint` of
        this surface at ``position`` in the polygon ``label``.

        INPUT:

        - ``label`` -- a label of a polygon in this surface, see :meth:`label_iterator`

        - ``position`` -- a vector with coordinates in this surface's base ring

        EXAMPLES::

            sage: from flatsurf import Polygon
            sage: from flatsurf.geometry.surface import Surface_list

            sage: S = Surface_list(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.set_edge_pairing(0, 0, 0, 2)
            sage: S.set_edge_pairing(0, 1, 0, 3)

            sage: S.point(0, (0, 0))
            Vertex 0 of polygon 0

        """
        return self(label, position, limit=limit, ring=ring)

    def _an_element_(self):
        r"""
        Return a point of this surface.

        EXAMPLES::

            sage: from flatsurf import Polygon
            sage: from flatsurf.geometry.surface import Surface_list

            sage: S = Surface_list(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.set_edge_pairing(0, 0, 0, 2)
            sage: S.set_edge_pairing(0, 1, 0, 3)

            sage: S.an_element()
            Point (1/2, 1/2) of polygon 0

        """
        label = next(self.label_iterator())
        polygon = self.polygon(label)

        # We use a point that can be constructed without problems on an
        # infinite surface.
        if polygon.is_convex():
            coordinates = polygon.centroid()
        else:
            # Sometimes, this is not implemented because it requires the edge
            # transformation to be known, so we prefer the centroid.
            coordinates = polygon.edge(0) / 2
        return self.point(label, coordinates)

    def set_default_graphical_surface(self, graphical_surface):
        r"""
        Replace the default graphical surface with the provided GraphicalSurface.
        """
        from flatsurf.graphical.surface import GraphicalSurface

        if not isinstance(graphical_surface, GraphicalSurface):
            raise ValueError("graphical_surface must be a GraphicalSurface")
        if self != graphical_surface.get_surface():
            raise ValueError(
                "The provided graphical_surface renders a different surface!"
            )
        self._gs = graphical_surface

    def graphical_surface(self, *args, **kwds):
        r"""
        Return a GraphicalSurface representing this surface.

        By default this returns a cached version of the GraphicalSurface. If
        ``cached=False`` is provided as a keyword option then a new
        GraphicalSurface is returned.

        All other parameters are passed on to
        :class:`~flatsurf.graphical.surface.GraphicalSurface` or its
        :meth:`~flatsurf.graphical.surface.GraphicalSurface.process_options`.
        Note that this mutates the cached graphical surface for future calls.

        EXAMPLES:

        Test the difference between the cached graphical_surface and the uncached version:

        .. jupyter-execute::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.octagon_and_squares()
            sage: s.plot()
            ...Graphics object consisting of 32 graphics primitives

        .. jupyter-execute::

            sage: s.graphical_surface(adjacencies=[]).plot()
            ...Graphics object consisting of 18 graphics primitives

        """
        from flatsurf.graphical.surface import GraphicalSurface

        if "cached" in kwds:
            if not kwds["cached"]:
                # cached=False: return a new surface.
                kwds.pop("cached", None)
                return GraphicalSurface(self, *args, **kwds)
            kwds.pop("cached", None)
        if hasattr(self, "_gs"):
            self._gs.process_options(*args, **kwds)
        else:
            self._gs = GraphicalSurface(self, *args, **kwds)
        return self._gs

    def plot(self, *args, **kwds):
        r"""
        Return a plot of the surface.

        The parameters are passed on to :meth:`graphical_surface` and
        :meth:`flatsurf.graphical.surface.GraphicalSurface.plot`. Consult their
        documentation for details.

        EXAMPLES:

        .. jupyter-execute::

            sage: import flatsurf
            sage: S = flatsurf.translation_surfaces.veech_double_n_gon(5)
            sage: S.plot()
            ...Graphics object consisting of 21 graphics primitives

        Keywords are passed on to the underlying plotting routines, see
        :meth:`flatsurf.graphical.surface.GraphicalSurface.plot` for details:

        .. jupyter-execute::

            sage: S.plot(fill=None)
            ...Graphics object consisting of 21 graphics primitives

        Note that some keywords mutate the underlying cached graphical surface,
        see :meth:`graphical_surface`:

        .. jupyter-execute::

            sage: S.plot(edge_labels='gluings and number')
            ...Graphics object consisting of 23 graphics primitives

        """
        if len(args) > 1:
            raise ValueError("plot() can take at most one non-keyword argument")

        graphical_surface_keywords = {
            key: kwds.pop(key)
            for key in [
                "cached",
                "adjacencies",
                "polygon_labels",
                "edge_labels",
                "default_position_function",
            ]
            if key in kwds
        }

        if len(args) == 1:
            from flatsurf.graphical.surface import GraphicalSurface

            if not isinstance(args[0], GraphicalSurface):
                raise TypeError("non-keyword argument must be a GraphicalSurface")

            import warnings

            warnings.warn(
                "Passing a GraphicalSurface to plot() is deprecated because it mutates that GraphicalSurface. This functionality will be removed in a future version of sage-flatsurf. "
                "Call process_options() and plot() on the GraphicalSurface explicitly instead."
            )

            gs = args[0]
            gs.process_options(**graphical_surface_keywords)
        else:
            # It's very surprising that plot mutates the underlying cached
            # graphical surface. We should change that and make the graphical
            # surface not cached. See
            # https://github.com/flatsurf/sage-flatsurf/issues/97
            gs = self.graphical_surface(**graphical_surface_keywords)

        return gs.plot(**kwds)


class Surface_list(Surface):
    r"""
    A fast mutable :class:`Surface` using a list to store polygons and gluings.

    ALGORITHM:

    Internally, we maintain a list ``_p`` for storing polygons together with
    gluing data.

    Each ``_p[label]`` is typically a pair ``(polygon, gluing_list)`` where
    ``gluing_list`` is a list of pairs ``(other_label, other_edge)`` such that
    :meth:`opposite_edge(label, edge) <Surface.opposite_edge>` returns
    ``_p[label][1][edge]``.

    INPUT:

    - ``base_ring`` -- ring or ``None`` (default: ``None``); the ring
      containing the coordinates of the vertices of the polygons. If ``None``,
      the base ring will be the one of ``surface``.

    - ``surface`` -- a surface or ``None`` (default: ``None``); a surface to be
      copied or referenced (see ``copy``). If ``None``, the surface is
      initially empty.

    - ``copy`` -- boolean or ``None`` (default: ``None``); whether the data
      underlying ``surface`` is copied into this surface or just a reference to
      that surface is kept. If ``None``, a sensible default is chosen, namely
      ``surface.is_mutable()``.

    - ``mutable`` -- boolean or ``None`` (default: ``None``); whether this
      surface is mutable. When ``None``, the surface will be mutable iff
      ``surface`` is ``None``.

    EXAMPLES::

        sage: from flatsurf import polygons, Surface_list, Polygon, similarity_surfaces
        sage: p=polygons.regular_ngon(5)
        sage: s=Surface_list(base_ring=p.base_ring())
        doctest:warning
        ...
        UserWarning: Surface_list has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
        sage: s.add_polygon(p) # gets label 0
        0
        sage: s.add_polygon( (-matrix.identity(2))*p ) # gets label 1
        1
        sage: s.change_polygon_gluings(0,[(1,e) for e in range(5)])
        sage: # base label defaults to zero.
        sage: s.set_immutable()
        sage: TestSuite(s).run()

    We surgically add a square into an infinite billiard surface::

        sage: p = Polygon(vertices=[(0,0),(4,0),(0,3)])
        sage: s = similarity_surfaces.billiard(p)
        sage: ts=s.minimal_cover(cover_type="translation").copy(relabel=True, mutable=True)
        doctest:warning
        ...
        UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead.
        Use relabel({old: new for (new, old) in enumerate(surface.labels())}) for integer labels. However, there is no immediate replacement for lazy copying of infinite surfaces.
        Have a look at the implementation of flatsurf.geometry.delaunay.LazyMutableSurface and adapt it to your needs.
        sage: # Explore the surface a bit
        sage: ts.polygon(0)
        Polygon(vertices=[(0, 0), (4, 0), (0, 3)])
        sage: ts.opposite_edge(0,0)
        (1, 2)
        sage: ts.polygon(1)
        Polygon(vertices=[(0, 0), (0, -3), (4, 0)])
        sage: s = ts
        sage: l=s.add_polygon(polygons.square(side=4))
        sage: s.change_edge_gluing(0,0,l,2)
        sage: s.change_edge_gluing(1,2,l,0)
        sage: s.change_edge_gluing(l,1,l,3)
        sage: print("Glued in label is "+str(l))
        Glued in label is 2
        sage: count = 0
        sage: for x in ts.gluings():
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
        0 -> Polygon(vertices=[(0, 0), (4, 0), (0, 3)])
        2 -> Polygon(vertices=[(0, 0), (4, 0), (4, 4), (0, 4)])
        3 -> Polygon(vertices=[(0, 0), (-72/25, -21/25), (28/25, -96/25)])
        4 -> Polygon(vertices=[(0, 0), (0, 3), (-4, 0)])
        1 -> Polygon(vertices=[(0, 0), (0, -3), (4, 0)])
        5 -> Polygon(vertices=[(0, 0), (-28/25, 96/25), (-72/25, -21/25)])

    TESTS:

    Verify that the replacement for this class implements all the features that
    used to be around::

        sage: from flatsurf import MutableOrientedSimilaritySurface, Surface_list
        sage: legacy_methods = set(dir(Surface_list(QQ)))
        doctest:warning
        ...
        UserWarning: Surface_list has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
        sage: non_legacy_methods = set(dir(MutableOrientedSimilaritySurface(QQ)))
        sage: [method for method in legacy_methods if method not in non_legacy_methods and not method.startswith("_")]
        []

    """

    def __init__(
        self,
        base_ring=None,
        surface=None,
        copy=None,
        mutable=None,
        category=None,
        deprecation_warning=True,
    ):
        if deprecation_warning:
            import warnings

            warnings.warn(
                "Surface_list has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead"
            )

        self._p = []  # list of pairs (polygon, gluings)
        self._reference_surface = None
        self._removed_labels = []
        self._num_polygons = 0

        # Validate input parameters and fill in defaults
        (
            base_ring,
            surface,
            copy,
            mutable,
            finite,
        ) = Surface_list._validate_init_parameters(
            base_ring=base_ring,
            surface=surface,
            copy=copy,
            mutable=mutable,
            finite=None,
        )

        Surface.__init__(
            self,
            base_ring=base_ring,
            base_label=0,
            finite=finite,
            mutable=True,
            category=category,
            deprecation_warning=False,
        )

        # Initialize surface from reference surface
        if surface is not None:
            if copy is True:
                reference_label_to_label = {
                    label: self.add_polygon(polygon)
                    for label, polygon in zip(surface.labels(), surface.polygons())
                }

                for (
                    (label, edge),
                    (glued_label, glued_edge),
                ) in surface.gluings():
                    self.set_edge_pairing(
                        reference_label_to_label[label],
                        edge,
                        reference_label_to_label[glued_label],
                        glued_edge,
                    )

                self.change_base_label(reference_label_to_label[surface.root()])
            else:
                self._reference_surface = surface
                self._ref_to_int = {}
                self._int_to_ref = []

                if not surface.is_finite_type():
                    from sage.all import infinity

                    self._num_polygons = infinity
                else:
                    self._num_polygons = len(surface.polygons())

                # Cache the base polygon
                self.change_base_label(self.__get_label(surface.root()))
                assert self.root() == 0

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
            from flatsurf.geometry.surface import OrientedSimilaritySurface

            if not isinstance(surface, OrientedSimilaritySurface):
                raise TypeError("surface must be an OrientedSimilaritySurface")

            if not surface.is_finite_type() and surface.is_mutable():
                raise NotImplementedError(
                    "Cannot create surface from infinite mutable surface yet."
                )

            if base_ring is None:
                base_ring = surface.base_ring()

            if base_ring != surface.base_ring():
                raise NotImplementedError(
                    "Cannot provide both a surface and a base_ring yet."
                )

            if mutable is None:
                mutable = True

            if copy is None:
                copy = surface.is_mutable()

            if copy and not surface.is_finite_type():
                raise ValueError("Cannot copy infinite surface.")

            if surface.is_mutable() and not copy:
                raise ValueError("Cannot reference mutable surface.")

            finite = surface.is_finite_type()

        return base_ring, surface, copy, mutable, finite

    def __get_label(self, ref_label):
        r"""
        Returns a corresponding label. Creates a new label if necessary.
        """
        try:
            return self._ref_to_int[ref_label]
        except KeyError:
            polygon = self._reference_surface.polygon(ref_label)
            data = [polygon, [None for i in range(len(polygon.vertices()))]]
            if len(self._removed_labels) > 0:
                i = self._removed_labels.pop()
                self._p[i] = data
                self._ref_to_int[ref_label] = i
                self._int_to_ref[i] = ref_label
            else:
                i = len(self._p)
                if i != len(self._int_to_ref):
                    raise RuntimeError(
                        "length of self._int_to_ref is "
                        + str(len(self._int_to_ref))
                        + " should be the same as i="
                        + str(i)
                    )
                self._p.append(data)
                self._ref_to_int[ref_label] = i
                self._int_to_ref.append(ref_label)
            return i

    def is_compact(self):
        if self._reference_surface is not None:
            # Since we are only modifying a finite number of polygons, this
            # surface is compact iff its reference surface is.
            return self._reference_surface.is_compact()

        return True

    def change_base_label(self, new_base_label):
        r"""
        Change the base_label to the provided label.
        """
        super().change_base_label(int(new_base_label))

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

    def opposite_edge(self, p, e=None):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        if e is None:
            p, e = p
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
                pp = self.__get_label(ref_pp)
                return_value = (pp, ref_ee)
                glue[e] = return_value
                return return_value
        else:
            # Successfully return edge data
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
        data[0] = new_polygon
        if data[1] is None or len(new_polygon.vertices()) != len(data[1]):
            data[1] = [None for e in range(len(new_polygon.vertices()))]
        if gluing_list is not None:
            self.change_polygon_gluings(label, gluing_list)

    def _set_edge_pairing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by change_edge_gluing(). Should not be called directly.
        """
        try:
            data = self._p[label1]
        except KeyError:
            raise ValueError("No known polygon with provided label1=" + str(label1))
        if data is None:
            raise ValueError(
                "Provided label1=" + str(label1) + " was removed from the surface."
            )
        data[1][edge1] = (label2, edge2)
        try:
            data = self._p[label2]
        except KeyError:
            raise ValueError("No known polygon with provided label2=" + str(label2))
        if data is None:
            raise ValueError(
                "Provided label2=" + str(label2) + " was removed from the surface."
            )
        data[1][edge2] = (label1, edge1)

    _change_edge_gluing = _set_edge_pairing

    def _add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.

        EXAMPLES::

            sage: from flatsurf import polygons, Surface_list
            sage: p=polygons.regular_ngon(5)
            sage: s=Surface_list(base_ring=p.base_ring())
            sage: s.add_polygon(p, label=3)
            3
            sage: s.add_polygon( (-matrix.identity(2))*p, label=30)
            30
            sage: s.change_polygon_gluings(3,[(30,e) for e in range(5)])
            sage: s.change_base_label(30)
            sage: len(s.polygons())
            2
            sage: TestSuite(s).run()
            sage: s.remove_polygon(3)
            sage: s.add_polygon(p, label=6)
            6
            sage: s.change_polygon_gluings(6,[(30,e) for e in range(5)])
            sage: len(s.polygons())
            2
            sage: TestSuite(s).run()
            sage: s.change_base_label(6)
            sage: s.remove_polygon(30)
            sage: label = s.add_polygon((-matrix.identity(2))*p)
            sage: s.change_polygon_gluings(6,[(label,e) for e in range(5)])
            sage: TestSuite(s).run()
        """
        if new_polygon is None:
            data = [None, None]
        else:
            data = [new_polygon, [None for i in range(len(new_polygon.vertices()))]]
        if label is None:
            if len(self._removed_labels) > 0:
                new_label = self._removed_labels.pop()
                self._p[new_label] = data
            else:
                new_label = len(self._p)
                self._p.append(data)
                if self._reference_surface is not None:
                    # Need a blank in this list for algorithmic reasons
                    self._int_to_ref.append(None)
        else:
            new_label = int(label)
            if new_label < len(self._p):
                if self._p[new_label] is not None:
                    raise ValueError(
                        "Trying to add a polygon with label="
                        + str(label)
                        + " which already indexes a polygon."
                    )
                self._p[new_label] = data
            else:
                if new_label - len(self._p) > 100:
                    raise ValueError(
                        "Adding a polygon with label="
                        + str(label)
                        + " would add more than 100 entries in our list."
                    )
                for i in range(len(self._p), new_label):
                    self._p.append(None)
                    self._removed_labels.append(i)
                    if self._reference_surface is not None:
                        # Need a blank in this list for algorithmic reasons
                        self._int_to_ref.append(None)

                self._p.append(data)
                if self._reference_surface is not None:
                    # Need a blank in this list for algorithmic reasons
                    self._int_to_ref.append(None)

        if gluing_list is not None:
            self.change_polygon_gluings(new_label, gluing_list)
        self._num_polygons += 1

        return new_label

    def num_polygons(self):
        r"""
        Return the number of polygons making up the surface in constant time.
        """
        return self._num_polygons

    def _test_roots(self, **options):
        # Surface_list does not iterate labels in a canonical order from the
        # roots(). Instead, it iterates by increasing labels.
        pass

    def label_iterator(self, polygons=False):
        r"""
        Iterator over all polygon labels.
        """
        if polygons:
            yield from zip(self.labels(), self.polygons())
            return
        if self._reference_surface is not None:
            yield from Surface.label_iterator(self)
        elif self._num_polygons == len(self._p):
            yield from range(self.num_polygons())
        else:
            # We've removed some labels
            found = 0
            i = 0
            while found < self._num_polygons:
                if self._p[i] is not None:
                    found += 1
                    yield i
                i += 1

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        if label == len(self._p) - 1:
            self._p.pop()
            if self._reference_surface is not None:
                ref_label = self._int_to_ref.pop()
                assert len(self._int_to_ref) == label
                if ref_label is not None:
                    del self._ref_to_int[ref_label]
        else:
            self._p[label] = None
            self._removed_labels.append(label)
            if self._reference_surface is not None:
                ref_label = self._int_to_ref[label]
                if ref_label is not None:
                    self._int_to_ref[label] = None
                    del self._ref_to_int[ref_label]
        self._num_polygons -= 1

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import Surface_list, polygons
            sage: P=polygons.regular_ngon(5)
            sage: S = Surface_list(base_ring=P.base_ring())
            doctest:warning
            ...
            UserWarning: Surface_list has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
            sage: T = Surface_list(base_ring=P.base_ring())

            sage: S == T
            True

            sage: S.add_polygon(P, label=3)
            3

            sage: S == T
            False

        """
        if not isinstance(other, Surface_list):
            return False

        if self._reference_surface is not None:
            equal = self._eq_reference_surface(other)
            if equal is True:
                return True
            if equal is False:
                return False

        return super().__eq__(other)

    def _eq_reference_surface(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other`` by
        comparing their reference surfaces.

        Returns ``None``, when no conclusion could be reached.

        This is a helper method for :meth:`__eq__`.
        """
        if self._reference_surface != other._reference_surface:
            return None

        for label in range(len(self._p)):
            if self._p[label] is None:
                if label >= len(other._p) or other._p[label] is not None:
                    return None
                continue
            try:
                if self.polygon(label) != other.polygon(label):
                    return False
            except ValueError:
                return False

        for label in range(len(other._p)):
            if other._p[label] is None:
                if label >= len(self._p) or self._p[label] is not None:
                    return None
                continue
            try:
                if self.polygon(label) != other.polygon(label):
                    return False
            except ValueError:
                return False

        if self.base_ring() != other.base_ring():
            return False
        if self.is_mutable() != other.is_mutable():
            return False
        if self.base_label() != other.base_label():
            return False

        return True


def surface_list_from_polygons_and_gluings(polygons, gluings, mutable=False):
    r"""
    Take a list of polygons and gluings (given either as a list of pairs of edges, or as a dictionary),
    and produce a Surface_list from it. The mutable parameter determines the mutability of the resulting
    surface.
    """
    import warnings

    warnings.warn(
        "surface_list_from_polygons_and_gluings() has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead and add polygons and gluings explicitly"
    )

    if not (isinstance(polygons, list) or isinstance(polygons, tuple)):
        raise ValueError("polygons must be a list or tuple.")
    field = polygons[0].parent().field()
    s = Surface_list(base_ring=field, deprecation_warning=False)
    for p in polygons:
        s.add_polygon(p)
    try:
        # dict case:
        it = gluings.items()
    except AttributeError:
        # list case:
        it = gluings
    for (l1, e1), (l2, e2) in it:
        s.change_edge_gluing(l1, e1, l2, e2)
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
    :meth:`opposite_edge(label, edge) <Surface.opposite_edge>` returns
    ``_p[label][1][edge]``.

    INPUT:

    - ``base_ring`` -- ring or ``None`` (default: ``None``); the ring
      containing the coordinates of the vertices of the polygons. If ``None``,
      the base ring will be the one of ``surface``.

    - ``surface`` -- a surface or ``None`` (default: ``None``); a surface to be
      copied or referenced (see ``copy``). If ``None``, the surface is
      initially empty.

    - ``copy`` -- boolean or ``None`` (default: ``None``); whether the data
      underlying ``surface`` is copied into this surface or just a reference to
      that surface is kept. If ``None``, a sensible default is chosen, namely
      ``surface.is_mutable()``.

    - ``mutable`` -- boolean or ``None`` (default: ``None``); whether this
      surface is mutable. When ``None``, the surface will be mutable iff
      ``surface`` is ``None``.

    EXAMPLES::

        sage: from flatsurf import polygons, Surface_dict
        sage: p=polygons.regular_ngon(10)
        sage: s=Surface_dict(base_ring=p.base_ring())
        doctest:warning
        ...
        UserWarning: Surface_dict has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
        sage: s.add_polygon(p,label="A")
        'A'
        sage: s.change_polygon_gluings("A",[("A",(e+5)%10) for e in range(10)])
        sage: s.change_base_label("A")
        sage: s.set_immutable()
        sage: TestSuite(s).run()

    TESTS:

    Verify that the replacement for this class implements all the features that
    used to be around::

        sage: from flatsurf import MutableOrientedSimilaritySurface, Surface_dict
        sage: legacy_methods = set(dir(Surface_dict(QQ)))
        doctest:warning
        ...
        UserWarning: Surface_dict has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
        sage: non_legacy_methods = set(dir(MutableOrientedSimilaritySurface(QQ)))
        sage: [method for method in legacy_methods if method not in non_legacy_methods and not method.startswith("_")]
        []

    """

    def __init__(
        self,
        base_ring=None,
        surface=None,
        copy=None,
        mutable=None,
        category=None,
        deprecation_warning=True,
    ):
        if deprecation_warning:
            import warnings

            warnings.warn(
                "Surface_dict has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead"
            )

        self._p = {}
        self._reference_surface = None

        # Validate input parameters and fill in defaults
        (
            base_ring,
            surface,
            copy,
            mutable,
            finite,
        ) = Surface_list._validate_init_parameters(
            base_ring=base_ring,
            surface=surface,
            copy=copy,
            mutable=mutable,
            finite=None,
        )

        Surface.__init__(
            self,
            base_ring=base_ring,
            base_label=None,
            finite=finite,
            mutable=True,
            category=category,
            deprecation_warning=False,
        )

        # Initialize surface from reference surface
        if surface is not None:
            base_label = surface.root()

            if copy is True:
                reference_label_to_label = {
                    label: self.add_polygon(polygon, label=label)
                    for label, polygon in zip(surface.labels(), surface.polygons())
                }

                for (
                    (label, edge),
                    (glued_label, glued_edge),
                ) in surface.gluings():
                    self.set_edge_pairing(
                        reference_label_to_label[label],
                        edge,
                        reference_label_to_label[glued_label],
                        glued_edge,
                    )

                self.change_base_label(reference_label_to_label[base_label])
            else:
                self._reference_surface = surface

                self.change_base_label(base_label)

        if not mutable:
            self.set_immutable()

    def is_compact(self):
        if self._reference_surface is not None:
            return self._reference_surface.is_compact()

        return True

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
            data = [
                self._reference_surface.polygon(lab),
                [
                    self._reference_surface.opposite_edge(lab, e)
                    for e in range(len(polygon.vertices()))
                ],
            ]
            self._p[lab] = data

        if data is None:
            raise ValueError("Label " + str(lab) + " was removed from the surface.")
        return data[0]

    def opposite_edge(self, p, e=None):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        if e is None:
            p, e = p
        try:
            data = self._p[p]
        except KeyError:
            self.polygon(p)
            data = self._p[p]
        if data is None:
            raise ValueError("Label " + str(p) + " was removed from the surface.")
        gluing_data = data[1]
        try:
            return gluing_data[e]
        except IndexError:
            raise ValueError(
                "Edge e=" + str(e) + " is out of range in polygon with label " + str(p)
            )

    # Methods for changing the surface

    def _change_polygon(self, label, new_polygon, gluing_list=None):
        r"""
        Internal method used by change_polygon(). Should not be called directly.
        """
        try:
            data = self._p[label]
            if data is None:
                raise ValueError(
                    "Label " + str(label) + " was removed from the surface."
                )
            data[0] = new_polygon
        except KeyError:
            # Polygon probably lies in reference surface
            if self._reference_surface is None:
                raise ValueError("No known polygon with provided label")
            else:
                # Ensure the reference surface had a polygon with the provided label:
                old_polygon = self._reference_surface.polygon(label)
                if len(old_polygon.vertices()) == len(new_polygon.vertices()):
                    data = [
                        new_polygon,
                        [
                            self._reference_surface.opposite_edge(label, e)
                            for e in range(len(new_polygon.vertices()))
                        ],
                    ]
                    self._p[label] = data
                else:
                    data = [
                        new_polygon,
                        [None for e in range(len(new_polygon.vertices()))],
                    ]
                    self._p[label] = data
        if len(data[1]) != len(new_polygon.vertices()):
            data[1] = [None for e in range(len(new_polygon.vertices()))]
        if gluing_list is not None:
            self.change_polygon_gluings(label, gluing_list)

    def _set_edge_pairing(self, label1, edge1, label2, edge2):
        r"""
        Internal method used by set_edge_pairing(). Should not be called directly.
        """
        try:
            data = self._p[label1]
        except KeyError:
            if self._reference_surface is None:
                raise ValueError(
                    "No known polygon with provided label1 = {}".format(label1)
                )
            else:
                # Failure likely because reference_surface contains the polygon.
                # import the data into this surface.
                polygon = self._reference_surface.polygon(label1)
                data = [
                    polygon,
                    [
                        self._reference_surface.opposite_edge(label1, e)
                        for e in range(len(polygon.vertices()))
                    ],
                ]
                self._p[label1] = data
        try:
            data[1][edge1] = (label2, edge2)
        except IndexError:
            # break down error
            if data is None:
                raise ValueError("polygon with label1={} was removed".format(label1))
            data1 = data[1]
            try:
                data1[edge1] = (label2, edge2)
            except IndexError:
                raise ValueError(
                    "edge1={} is out of range in polygon with label1={}".format(
                        edge1, label1
                    )
                )
        try:
            data = self._p[label2]
        except KeyError:
            if self._reference_surface is None:
                raise ValueError("no polygon with label2={}".format(label2))
            else:
                # Failure likely because reference_surface contains the polygon.
                # import the data into this surface.
                polygon = self._reference_surface.polygon(label2)
                data = [
                    polygon,
                    [
                        self._reference_surface.opposite_edge(label2, e)
                        for e in range(len(polygon.vertices()))
                    ],
                ]
                self._p[label2] = data
        try:
            data[1][edge2] = (label1, edge1)
        except IndexError:
            # break down error
            if data is None:
                raise ValueError("polygon with label1={} was removed".format(label1))
            data1 = data[1]
            try:
                data1[edge2] = (label1, edge1)
            except IndexError:
                raise ValueError(
                    "edge {} is out of range in polygon with label2={}".format(
                        edge2, label2
                    )
                )

    _change_edge_gluing = _set_edge_pairing

    def _add_polygon(self, new_polygon, gluing_list=None, label=None):
        r"""
        Internal method used by add_polygon(). Should not be called directly.
        """
        data = [new_polygon, [None for i in range(len(new_polygon.vertices()))]]
        if label is None:
            new_label = ExtraLabel()
        else:
            try:
                old_data = self._p[label]
                if old_data is None:
                    # already removed this polygon. That's good, we can add.
                    new_label = label
                else:
                    raise ValueError(
                        "label={} already used by another polygon".format(label)
                    )
            except KeyError:
                # This seems inconvenient to enforce:
                #
                # if not self._reference_surface is None:
                #    # Can not be sure we are not overwriting a polygon in the reference surface.
                #    raise ValueError("Can not assign this label to a Surface_dict containing a reference surface,"+\
                #        "which may already contain this label.")
                new_label = label
        self._p[new_label] = data
        if gluing_list is not None:
            self.change_polygon_gluings(new_label, gluing_list)
        return new_label

    def num_polygons(self):
        r"""
        Return the number of polygons making up the surface in constant time.
        """
        if self.is_finite_type():
            if self._reference_surface is None:
                return len(self._p)
            else:
                # Unfortunately, I don't see a good way to compute this.
                return Surface.num_polygons(self)
        else:
            from sage.rings.infinity import Infinity

            return Infinity

    def label_iterator(self, polygons=False):
        r"""
        Iterator over all polygon labels.
        """
        if polygons:
            yield from zip(self.labels(), self.polygons())
            return
        if self._reference_surface is None:
            yield from self._p
        else:
            yield from Surface.label_iterator(self)

    def _remove_polygon(self, label):
        r"""
        Internal method used by remove_polygon(). Should not be called directly.
        """
        if self._reference_surface is None:
            try:
                data = self._p[label]
            except KeyError:
                raise ValueError("Label " + str(label) + " is not in the surface.")
            del self._p[label]
        else:
            try:
                data = self._p[label]
                # Success.
                if data is None:
                    raise ValueError(
                        "Label " + str(label) + " was already removed from the surface."
                    )
                self._p[label] = None
            except KeyError:
                # Assume on faith we are removing a polygon in the base_surface.
                self._p[label] = None

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import Surface_dict, polygons
            sage: P=polygons.regular_ngon(5)
            sage: S = Surface_dict(base_ring=P.base_ring())
            sage: T = Surface_dict(base_ring=P.base_ring())

            sage: S == T
            True

            sage: S.add_polygon(P, label=3)
            3

            sage: S == T
            False

        """
        if not isinstance(other, Surface_dict):
            return False

        if self._reference_surface is not None:
            equal = self._eq_reference_surface(other)
            if equal is True:
                return True
            if equal is False:
                return False

        return super().__eq__(other)

    def _eq_reference_surface(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other`` by
        comparing their reference surfaces.

        Returns ``None``, when no conclusion could be reached.

        This is a helper method for :meth:`__eq__`.
        """
        if self._reference_surface != other._reference_surface:
            return None

        for label, polygon in self._p.items():
            if polygon is None:
                if label not in other._p or other._p[label] is not None:
                    return None
                continue
            try:
                if self.polygon(label) != other.polygon(label):
                    return False
            except ValueError:
                return False

        for label, polygon in other._p.items():
            if polygon is None:
                if label not in self._p or self._p[label] is not None:
                    return None
                continue
            try:
                if self.polygon(label) != other.polygon(label):
                    return False
            except ValueError:
                return False

        if self.base_ring() != other.base_ring():
            return False
        if self.is_mutable() != other.is_mutable():
            return False
        if self.base_label() != other.base_label():
            return False

        return True


class LabelWalker:
    r"""
    Take a canonical walk around the surface and find the labels of polygons.

    We start at the base_label().
    Then the labels are visited in order involving the combinatorial distance from the base_label(),
    where combinatorial distance measures the minimal number of edges which need to be crossed to reach the
    polygon with a givel label. Ties are broken using lexicographical order on the numbers associated to edges crossed
    (labels are not used in this lexicographical ordering).
    """

    class LabelWalkerIterator:
        def __init__(self, label_walker):
            self._lw = label_walker
            self._i = 0

        def __next__(self):
            if self._i < len(self._lw):
                label = self._lw.number_to_label(self._i)
                self._i = self._i + 1
                return label
            if self._i == len(self._lw):
                label = self._lw.find_a_new_label()
                if label is None:
                    raise StopIteration()
                self._i = self._i + 1
                return label
            raise StopIteration()

        next = __next__  # for Python 2

        def __iter__(self):
            return self

    def __init__(self, surface, deprecation_warning=True):
        if deprecation_warning:
            import warnings

            warnings.warn(
                "LabelWalker has been deprecated and will be removed in a future version of sage-flatsurf; use labels() instead"
            )

        self._s = surface
        self._labels = [self._s.root()]
        self._label_dict = {self._labels[0]: 0}

        # This will stores an edge to move through to get to a polygon closer to the base_polygon
        self._label_edge_back = {self._labels[0]: None}

        from collections import deque

        self._walk = deque()
        self._walk.append((self._labels[0], 0))

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
                if not self._s.is_finite_type():
                    limit = 1000
                else:
                    limit = self._s.num_polygons()
            for i in range(limit):
                new_label = self.find_a_new_label()
                if label == new_label:
                    return self._label_edge_back[label]
        # Maybe the surface is not connected?
        raise KeyError(
            f"Unable to find label {label}. Are you sure the surface is connected?"
        )

    def __iter__(self):
        return LabelWalker.LabelWalkerIterator(self)

    def polygon_iterator(self):
        for label in self:
            yield self._s.polygon(label)

    def label_polygon_iterator(self):
        for label in self:
            yield label, self._s.polygon(label)

    def edge_iterator(self, gluings=False):
        if gluings:
            yield from self._s.gluings()
            return
        for label, polygon in self.label_polygon_iterator():
            for e in range(len(polygon.vertices())):
                yield label, e

    def __len__(self):
        r"""
        Return the number of labels found.
        """
        return len(self._labels)

    def find_a_new_label(self):
        r"""
        Finds a new label, stores it, and returns it. Returns None if we have already found all labels.
        """
        while len(self._walk) > 0:
            label, e = self._walk.popleft()
            opposite_label, opposite_edge = self._s.opposite_edge(label, e)
            e = e + 1
            if e < len(self._s.polygon(label).vertices()):
                self._walk.appendleft((label, e))
            if opposite_label not in self._label_dict:
                n = len(self._labels)
                self._labels.append(opposite_label)
                self._label_dict[opposite_label] = n
                self._walk.append((opposite_label, 0))
                self._label_edge_back[opposite_label] = opposite_edge
                return opposite_label
        return None

    def find_new_labels(self, n):
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
        if not self._s.is_finite_type():
            raise NotImplementedError
        label = self.find_a_new_label()
        while label is not None:
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
                new_label = self.find_a_new_label()
                if label == new_label:
                    return self._label_dict[label]
            raise ValueError(
                "Failed to find label even after searching. limit=" + str(limit)
            )

    def surface(self):
        return self._s


class ExtraLabel:
    r"""
    Used to spit out new labels.
    """

    _next = int(0)

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
        return isinstance(other, self.__class__) and self._label == other._label

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(23 * self._label)

    def __str__(self):
        return "E" + str(self._label)

    def __repr__(self):
        return "ExtraLabel(" + str(self._label) + ")"


class LabelsView(collections.abc.Set):
    def __init__(self, surface):
        self._surface = surface

    def __contains__(self, x):
        try:
            self._surface.polygon(x)
        except ValueError:
            return False

        return True

    def __iter__(self):
        return self._surface.label_iterator()

    def __len__(self):
        return self._surface.num_polygons()

    def __repr__(self):
        if self._surface.is_finite_type():
            return repr(tuple(self))

        from itertools import islice

        return f"({', '.join(str(x) for x in islice(self, 16))}, …)"


def SurfaceClass(surface, name, category, *args, **kwargs):
    category = category.Oriented().Connected().WithoutBoundary()

    message = f"{name} has been deprecated and will be removed in a future "
    message += "version of sage-flatsurf; there is no distinction between "
    message += "an (underlying) Surface and the SimilaritySurface types anymore."

    if surface.is_finite_type():
        message += f" Calling set_immutable() on this surface should determine the category of this surface automatically so calling {name} should not be necessary in this case."
    else:
        message += " For this surface of infinite type, you should create a subclass of OrientedSimilaritySurface and set the category in the __init__ method; see flatsurf.geometry.similarity_surface_generators.EInfinitySurface for an example"

    message += f" You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended. We will now refine the category of this surface to make sure that it is in the {category}."

    if args or kwargs:
        message += " Note that we are ignoring any other parameters that were passed to this function."

    import warnings

    warnings.warn(message)

    surface._refine_category_(category)

    return surface


def SimilaritySurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, SimilaritySurface
        sage: S = Surface_list(QQ)
        doctest:warning
        ...
        UserWarning: Surface_list has been deprecated and will be removed in a future version of sage-flatsurf; use MutableOrientedSimilaritySurface instead
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = SimilaritySurface(S)
        doctest:warning
        ...
        UserWarning: SimilaritySurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling SimilaritySurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary oriented similarity surfaces.
        sage: S.category()
        Category of connected without boundary finite type oriented similarity surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import SimilaritySurfaces

    return SurfaceClass(
        surface, "SimilaritySurface()", SimilaritySurfaces(), *args, **kwargs
    )


def HalfDilationSurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, HalfDilationSurface
        sage: S = Surface_list(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = HalfDilationSurface(S)
        doctest:warning
        ...
        UserWarning: HalfDilationSurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling HalfDilationSurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary oriented dilation surfaces.
        sage: S.category()
        Category of connected without boundary finite type oriented dilation surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import DilationSurfaces

    return SurfaceClass(
        surface, "HalfDilationSurface()", DilationSurfaces(), *args, **kwargs
    )


def DilationSurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, DilationSurface
        sage: S = Surface_list(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = DilationSurface(S)
        doctest:warning
        ...
        UserWarning: DilationSurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling DilationSurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary oriented positive dilation surfaces.
        sage: S.category()
        Category of connected without boundary finite type oriented positive dilation surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import DilationSurfaces

    return SurfaceClass(
        surface, "DilationSurface()", DilationSurfaces().Positive(), *args, **kwargs
    )


def ConeSurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, ConeSurface
        sage: S = Surface_list(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = ConeSurface(S)
        doctest:warning
        ...
        UserWarning: ConeSurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling ConeSurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary oriented cone surfaces.
        sage: S.category()
        Category of connected without boundary finite type oriented cone surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import ConeSurfaces

    return SurfaceClass(surface, "ConeSurface()", ConeSurfaces(), *args, **kwargs)


def RationalConeSurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, RationalConeSurface
        sage: S = Surface_list(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = RationalConeSurface(S)
        doctest:warning
        ...
        UserWarning: RationalConeSurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling RationalConeSurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary oriented rational cone surfaces.
        sage: S.category()
        Category of connected without boundary finite type oriented rational cone surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import ConeSurfaces

    return SurfaceClass(
        surface, "RationalConeSurface()", ConeSurfaces().Rational(), *args, **kwargs
    )


def HalfTranslationSurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, HalfTranslationSurface
        sage: S = Surface_list(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = HalfTranslationSurface(S)
        doctest:warning
        ...
        UserWarning: HalfTranslationSurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling HalfTranslationSurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary oriented half translation surfaces.
        sage: S.category()
        Category of connected without boundary finite type oriented half translation surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import HalfTranslationSurfaces

    return SurfaceClass(
        surface, "HalfTranslationSurface()", HalfTranslationSurfaces(), *args, **kwargs
    )


def TranslationSurface(surface, *args, **kwargs):
    r"""
    Refine the category of ``surface``.

    This function is deprecated and should not be used anymore.

    EXAMPLES::

        sage: from flatsurf import Surface_list, polygons, TranslationSurface
        sage: S = Surface_list(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S.set_edge_pairing(0, 0, 0, 2)
        sage: S.set_edge_pairing(0, 1, 0, 3)
        sage: S = TranslationSurface(S)
        doctest:warning
        ...
        UserWarning: TranslationSurface() has been deprecated and will be removed in a future version of sage-flatsurf; there is no distinction between an (underlying) Surface and the SimilaritySurface types anymore.
        Calling set_immutable() on this surface should determine the category of this surface automatically so calling TranslationSurface() should not be necessary in this case.
        You can still explicitly refine the category of a surface with _refine_category_() but this is not recommended.
        We will now refine the category of this surface to make sure that it is in the Category of connected without boundary translation surfaces.
        sage: S.category()
        Category of connected without boundary finite type translation surfaces
        sage: S.set_immutable()
        sage: S.category()
        Category of connected without boundary finite type translation surfaces

    """
    from flatsurf.geometry.categories import TranslationSurfaces

    return SurfaceClass(
        surface, "TranslationSurface()", TranslationSurfaces(), *args, **kwargs
    )
