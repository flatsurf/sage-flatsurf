r"""
Generic mutable and immutable surfaces

This module provides base classes and implementations of surfaces. Most
surfaces in sage-flatsurf inherit from some of the classes in this module.

The most important class in this module is
:class:`MutableOrientedSimilaritySurface` which allows you to create a surface
by gluing polygons with similarities.

EXAMPLES:

We build a translation surface by gluing two hexagons, labeled 0 and 1::

    sage: from flatsurf import MutableOrientedSimilaritySurface, polygons
    sage: S = MutableOrientedSimilaritySurface(QuadraticField(3))

    sage: S.add_polygon(polygons.regular_ngon(6))
    0
    sage: S.add_polygon(polygons.regular_ngon(6))
    1

    sage: S.glue((0, 0), (1, 3))
    sage: S.glue((0, 1), (1, 4))
    sage: S.glue((0, 2), (1, 5))
    sage: S.glue((0, 3), (1, 0))
    sage: S.glue((0, 4), (1, 1))
    sage: S.glue((0, 5), (1, 2))

    sage: S
    Translation Surface built from 2 regular hexagons

We signal that the construction is complete. This refines the category of the
surface and makes more functionality available::

    sage: S.set_immutable()
    sage: S
    Translation Surface in H_2(1^2) built from 2 regular hexagons

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                           2023 Julian Rüth
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

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.surface_objects import SurfacePoint


class Surface_base(Parent):
    r"""
    A base class for all surfaces in sage-flatsurf.

    This class patches bits of the category framework in SageMath that assume
    that all parent structures are immutable.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()

        sage: from flatsurf.geometry.surface import Surface_base
        sage: isinstance(S, Surface_base)
        True

    """

    def _refine_category_(self, category):
        r"""
        Refine the category of this surface to a subcategory ``category``.

        We need to override this method from ``Parent`` since we need to
        disable a hashing check that is otherwise enabled when doctesting.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.category()
            Category of finite type oriented similarity surfaces

            sage: S._refine_category_(S.refined_category())
            sage: S.category()
            Category of connected without boundary finite type translation surfaces

        """
        from sage.structure.debug_options import debug

        old_refine_category_hash_check = debug.refine_category_hash_check
        debug.refine_category_hash_check = False
        try:
            super()._refine_category_(category)
        finally:
            debug.refine_category_hash_check = old_refine_category_hash_check

        # The (cached) an_element is going to fail its test suite because it has the wrong category now.
        # Make sure that we create a new copy of an_element when requested.
        self._cache_an_element = None


class MutablePolygonalSurface(Surface_base):
    r"""
    A base class for mutable surfaces that are built by gluing polygons.

    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface
        sage: S = MutableOrientedSimilaritySurface(QQ)

        sage: from flatsurf.geometry.surface import MutablePolygonalSurface
        sage: isinstance(S, MutablePolygonalSurface)
        True

    """

    def __init__(self, base, category=None):
        from sage.all import ZZ

        self._next_label = ZZ(0)
        self._roots = ()

        self._polygons = {}

        self._mutable = True

        super().__init__(base, category=category)

    def _test_refined_category(self, **options):
        r"""
        Test that this surface has been refined to its best possible
        subcategory (that can be computed cheaply.)

        We override this method here to disable this check for mutable
        surfaces. Mutable surfaces have not been refined yet since changes in
        the surface might require a widening of the category which is not
        possible.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S._test_refined_category()
            sage: S.set_immutable()
            sage: S._test_refined_category()

        """
        if self._mutable:
            return

        super()._test_refined_category(**options)

    def add_polygon(self, polygon, *, label=None):
        r"""
        Add an unglued polygon to this surface and return its label.

        INPUT:

        - ``polygon`` -- a simple Euclidean polygon

        - ``label`` -- a hashable identifier or ``None`` (default: ``None``);
          if ``None`` an integer identifier is automatically selected

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.add_polygon(polygons.square(), label=0)
            Traceback (most recent call last):
            ...
            ValueError: polygon label already present in this surface

            sage: S.add_polygon(polygons.square(), label='X')
            'X'

        """
        if not self._mutable:
            raise Exception("cannot modify an immutable surface")

        if label is None:
            while self._next_label in self._polygons:
                self._next_label += 1
            label = self._next_label

        if label in self._polygons:
            raise ValueError("polygon label already present in this surface")

        self._polygons[label] = polygon.change_ring(self.base_ring())
        return label

    def add_polygons(self, polygons):
        r"""
        Add several polygons with automatically assigned labels at once.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygons([polygons.square(), polygons.square()])
            doctest:warning
            ...
            UserWarning: add_polygons() has been deprecated and will be removed in a future version of sage-flatsurf; use labels = [add_polygon(p) for p in polygons] instead
            [0, 1]

        """
        import warnings

        warnings.warn(
            "add_polygons() has been deprecated and will be removed in a future version of sage-flatsurf; use labels = [add_polygon(p) for p in polygons] instead"
        )

        return [self.add_polygon(p) for p in polygons]

    def set_default_graphical_surface(self, graphical_surface):
        r"""
        EXAMPLES:

        This has been disabled because it tends to break caching::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.set_default_graphical_surface(S.graphical_surface())
            Traceback (most recent call last):
            ...
            NotImplementedError: set_default_graphical_surface() has been removed from this version of sage-flatsurf. If you want to change the default plotting of a surface, create a subclass and override graphical_surface() instead

        """
        raise NotImplementedError(
            "set_default_graphical_surface() has been removed from this version of sage-flatsurf. If you want to change the default plotting of a surface, create a subclass and override graphical_surface() instead"
        )

    def remove_polygon(self, label):
        r"""
        Remove the polygon with label ``label`` from this surface (and all data
        associated to it.)

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.remove_polygon(0)

            sage: S.add_polygon(polygons.square())
            0

        """
        if not self._mutable:
            raise Exception("cannot modify an immutable surface")

        self._polygons.pop(label)
        self._roots = tuple(root for root in self._roots if root != label)

    def roots(self):
        r"""
        Return a label for each connected component on this surface.

        This implements :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.roots()
            ()

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.roots()
            (0,)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            1

            sage: S.roots()
            (0, 1)

            sage: S.glue((0, 0), (1, 0))
            sage: S.roots()
            (0,)

        .. SEEALSO::

            :meth:`components`

        """
        return LabeledView(
            self, RootedComponents_MutablePolygonalSurface(self).keys(), finite=True
        )

    def components(self):
        r"""
        Return the connected components as the sequence of their respective
        polygon labels.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.components()
            ()

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.components()
            ((0,),)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            1

            sage: S.components()
            ((0,), (1,))

            sage: S.glue((0, 0), (1, 0))
            sage: S.components()
            ((0, 1),)

        """
        return LabeledView(
            self, RootedComponents_MutablePolygonalSurface(self).values(), finite=True
        )

    def polygon(self, label):
        r"""
        Return the polygon with label ``label`` in this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.polygon(0)
            Traceback (most recent call last):
            ...
            KeyError: 0

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.polygon(0)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        return self._polygons[label]

    def set_immutable(self):
        r"""
        Make this surface immutable.

        Any mutation attempts from now on will be an error.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))

        Note that declaring a surface immutable refines its category and
        thereby unlocks more methods that are available to such a surface::

            sage: S.category()
            Category of finite type oriented similarity surfaces
            sage: old_methods = set(method for method in dir(S) if not method.startswith('_'))

            sage: S.set_immutable()
            sage: S.category()
            Category of connected without boundary finite type translation surfaces
            sage: new_methods = set(method for method in dir(S) if not method.startswith('_'))
            sage: new_methods - old_methods
            {'angles',
             'apply_matrix',
             'area',
             'canonicalize',
             'canonicalize_mapping',
             'erase_marked_points',
             'holonomy_field',
             'j_invariant',
             'l_infinity_delaunay_triangulation',
             'minimal_translation_cover',
             'normalized_coordinates',
             'rel_deformation',
             'stratum'}

        An immutable surface cannot be mutated anymore::

            sage: S.remove_polygon(0)
            Traceback (most recent call last):
            ...
            Exception: cannot modify an immutable surface

        However, the category of an immutable might be further refined as
        (expensive to determine) features of the surface are deduced.

        """
        if self._mutable:
            self.set_roots(self.roots())
            self._mutable = False

        self._refine_category_(self.refined_category())

    def is_mutable(self):
        r"""
        Return whether this surface can be modified.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.is_mutable()
            True

            sage: S.set_immutable()
            sage: S.is_mutable()
            False

        """
        return self._mutable

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See
        :meth:`~.categories.similarity_surfaces.SimilaritySurfaces.FiniteType.ParentMethods._test_eq_surface`
        for details on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: T = MutableOrientedSimilaritySurface(QQ)

            sage: S == T
            True

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S == T
            False

        TESTS::

            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: T = MutableOrientedSimilaritySurface(QQ)

            sage: S != T
            False

            sage: S.add_polygon(polygons.square())
            0

            sage: S != T
            True

        """
        if not isinstance(other, MutablePolygonalSurface):
            return False

        if self.base() != other.base():
            return False

        if self._polygons != other._polygons:
            return False

        # Note that the order of the root labels matters since it changes the order of iteration in labels()
        if self._roots != other._roots:
            return False

        if self._mutable != other._mutable:
            return False

        return True

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: T = MutableOrientedSimilaritySurface(QQ)
            sage: hash(S) == hash(T)
            Traceback (most recent call last):
            ...
            TypeError: cannot hash a mutable surface

            sage: S.set_immutable()
            sage: T.set_immutable()
            sage: hash(S) == hash(T)
            True

        """
        if self._mutable:
            raise TypeError("cannot hash a mutable surface")

        return hash((tuple(self.labels()), tuple(self.polygons()), self._roots))

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S
            Empty Surface

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S
            Translation Surface with boundary built from a square

        """
        if not self.is_finite_type():
            return "Surface built from infinitely many polygons"

        if len(self.labels()) == 0:
            return "Empty Surface"

        return f"{self._describe_surface()} built from {self._describe_polygons()}"

    def _describe_surface(self):
        r"""
        Return a string describing this kind of surface.

        This is a helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S._describe_surface()
            'Translation Surface'

        """
        return "Surface"

    def _describe_polygons(self):
        r"""
        Return a string describing the nature of the polygons that make up this surface.

        This is a helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S._describe_polygons()
            ''

        """
        polygons = [
            (-len(p.erase_marked_vertices().vertices()), p.describe_polygon())
            for p in self.polygons()
        ]
        polygons.sort()
        polygons = [description for (edges, description) in polygons]

        if not polygons:
            return ""

        collated = []
        while polygons:
            count = 1
            polygon = polygons.pop()
            while polygons and polygons[-1] == polygon:
                count += 1
                polygons.pop()

            if count == 1:
                collated.append(f"{polygon[0]} {polygon[1]}")
            else:
                collated.append(f"{count} {polygon[2]}")

        description = collated.pop()

        if collated:
            description = ", ".join(collated) + " and " + description

        return description

    def set_root(self, root):
        r"""
        Set ``root`` as the label at which exploration of a connected component
        starts.

        This method can be called for connected and disconnected surfaces. In
        either case, it establishes ``root`` as the new label from which
        enumeration of the connected component containing it starts. If another
        label for this component had been set earlier, it is replaced.

        .. NOTE::

            After roots have been declared explicitly, gluing operations come
            at an additional cost since the root labels have to be updated
            sometimes. It is therefore good practice to declare the root labels
            after all the gluings have been established when creating a
            surface.

        INPUT:

        - ``root`` -- a polygon label in this surface

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.set_root(0)
            Traceback (most recent call last):
            ...
            ValueError: root must be a label in the surface

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S.add_polygon(polygons.square())
            1

            sage: S.set_root(0)
            sage: S.set_root(1)

            sage: S.roots()
            (0, 1)

        Note that the roots get updated automatically when merging components::

            sage: S.glue((0, 0), (1, 0))
            sage: S.roots()
            (0,)

        The purpose of changing the root label is to modify the order of
        exploration, e.g., in :meth:`labels`::

            sage: S.labels()
            (0, 1)

            sage: S.set_root(1)
            sage: S.labels()
            (1, 0)

        .. SEEALSO::

            :meth:`set_roots` to replace all the root labels

        """
        if not self._mutable:
            raise Exception("cannot modify an immutable surface")

        root = [label for label in self.labels() if label == root]
        if not root:
            raise ValueError("root must be a label in the surface")
        assert len(root) == 1
        root = root[0]

        component = [component for component in self.components() if root in component]
        assert len(component) == 1
        component = component[0]

        self._roots = tuple(r for r in self._roots if r not in component) + (root,)

    def set_roots(self, roots):
        r"""
        Declare that the labels in ``roots`` are the labels from which their
        corresponding connected components should be enumerated.

        There must be at most one label for each connected component in
        ``roots``. Components that have no label set explicitly will have their
        label chosen automatically.

        INPUT:

        - ``roots`` -- a sequence of polygon labels in this surface

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S.add_polygon(polygons.square())
            1
            sage: S.add_polygon(polygons.square())
            2

            sage: S.glue((0, 0), (1, 0))

            sage: S.set_roots([1])
            sage: S.roots()
            (1, 2)

        Setting the roots of connected components affects their enumeration in :meth:`labels`::

            sage: S.labels()
            (1, 0, 2)

            sage: S.set_roots([0, 2])
            sage: S.labels()
            (0, 1, 2)

        There must be at most one root per component::

            sage: S.set_roots([0, 1, 2])
            Traceback (most recent call last):
            ...
            ValueError: there must be at most one root label for each connected component

        """
        if not self._mutable:
            raise Exception("cannot modify an immutable surface")

        roots = [[label for label in self.labels() if label == root] for root in roots]

        if any(len(root) == 0 for root in roots):
            raise ValueError("roots must be existing labels in the surface")

        assert all(len(root) == 1 for root in roots)

        roots = tuple(root[0] for root in roots)

        for component in self.components():
            if len([root for root in roots if root in component]) > 1:
                raise ValueError(
                    "there must be at most one root label for each connected component"
                )

        self._roots = tuple(roots)

    def change_base_label(self, label):
        r"""
        This is a deprecated alias for :meth:`set_root`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.change_base_label(0)
            doctest:warning
            ...
            UserWarning: change_base_label() has been deprecated and will be removed in a future version of sage-flatsurf; use set_root() instead

        """
        import warnings

        warnings.warn(
            "change_base_label() has been deprecated and will be removed in a future version of sage-flatsurf; use set_root() instead"
        )

        self.set_root(label)

    @cached_method
    def labels(self):
        r"""
        Return the polygon labels in this surface.

        This replaces the generic
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`
        method with a more efficient implementation.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.labels()
            ()

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0

            sage: S.labels()
            (0,)

        .. SEEALSO::

            :meth:`polygon` to translate polygon labels to the corresponding polygons

            :meth:`polygons` for the corresponding sequence of polygons

        """
        return LabelsFromView(self, self._polygons.keys(), finite=True)

    @cached_method
    def polygons(self):
        r"""
        Return the polygons that make up this surface.

        The order the polygons are returned is guaranteed to be compatible with
        the order of the labels in :meth:`~.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        This replaces the generic
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygons`
        with a more efficient implementation.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: S.polygons()
            ()

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S.add_polygon(polygons.square())
            1
            sage: S.add_polygon(polygons.square())
            2

            sage: S.polygons()
            (Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))

        .. SEEALSO::

            :meth:`polygon` to get a single polygon

        """
        return Polygons_MutableOrientedSimilaritySurface(self)


class OrientedSimilaritySurface(Surface_base):
    r"""
    Base class for surfaces built from Euclidean polygons that are glued with
    orientation preserving similarities.

    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface
        sage: S = MutableOrientedSimilaritySurface(QQ)

        sage: from flatsurf.geometry.surface import OrientedSimilaritySurface
        sage: isinstance(S, OrientedSimilaritySurface)
        True

    """
    Element = SurfacePoint

    def __init__(self, base, category=None):
        from sage.categories.all import Rings

        if base not in Rings():
            raise TypeError("base ring must be a ring")

        from flatsurf.geometry.categories import SimilaritySurfaces

        if category is None:
            category = SimilaritySurfaces().Oriented()

        category &= SimilaritySurfaces().Oriented()

        super().__init__(base, category=category)

    def _describe_surface(self):
        if not self.is_finite_type():
            return "Surface built from infinitely many polygons"

        if not self.is_connected():
            # Many checks do not work yet if a surface is not connected, so we stop here.
            return "Disconnected Surface"

        if self.is_translation_surface(positive=True):
            description = "Translation Surface"
        elif self.is_translation_surface(positive=False):
            description = "Half-Translation Surface"
        elif self.is_dilation_surface(positive=True):
            description = "Positive Dilation Surface"
        elif self.is_dilation_surface(positive=False):
            description = "Dilation Surface"
        elif self.is_cone_surface():
            description = "Cone Surface"
            if self.is_rational_surface():
                description = f"Rational {description}"
        else:
            description = "Surface"

        if hasattr(self, "stratum"):
            try:
                description += f" in {self.stratum()}"
            except NotImplementedError:
                # Computation of the stratum might fail due to #227.
                pass
        elif self.genus is not NotImplemented:
            description = f"Genus {self.genus()} {description}"

        if self.is_with_boundary():
            description += " with boundary"

        return description


class MutableOrientedSimilaritySurface_base(OrientedSimilaritySurface):
    r"""
    Base class for surface built from Euclidean polyogns glued by orientation
    preserving similarities.

    This provides the features of :class:`MutableOrientedSimilaritySurface`
    without making a choice about how data is stored internally; it is a
    generic base class for other implementations of mutable surfaces.

    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface
        sage: S = MutableOrientedSimilaritySurface(QQ)

        sage: from flatsurf.geometry.surface import MutableOrientedSimilaritySurface_base
        sage: isinstance(S, MutableOrientedSimilaritySurface_base)
        True

    """

    def triangle_flip(self, l1, e1, in_place=False, test=False, direction=None):
        r"""
        Overrides
        :meth:`.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.triangle_flip`
        to provide in-place flipping of triangles.

        See that method for details.
        """
        if not in_place:
            return super().triangle_flip(
                l1=l1, e1=e1, in_place=in_place, test=test, direction=direction
            )

        s = self

        p1 = s.polygon(l1)
        if not len(p1.vertices()) == 3:
            raise ValueError("The polygon with the provided label is not a triangle.")
        l2, e2 = s.opposite_edge(l1, e1)

        sim = s.edge_transformation(l2, e2)
        p2 = s.polygon(l2)
        if not len(p2.vertices()) == 3:
            raise ValueError(
                "The polygon opposite the provided edge is not a triangle."
            )

        from flatsurf import Polygon

        p2 = Polygon(vertices=[sim(v) for v in p2.vertices()], base_ring=p1.base_ring())

        if direction is None:
            direction = (s.base_ring() ** 2)((0, 1))
        # Get vertices corresponding to separatices in the provided direction.
        v1 = p1.find_separatrix(direction=direction)[0]
        v2 = p2.find_separatrix(direction=direction)[0]
        # Our quadrilateral has vertices labeled:
        # * 0=p1.vertex(e1+1)=p2.vertex(e2)
        # * 1=p1.vertex(e1+2)
        # * 2=p1.vertex(e1)=p2.vertex(e2+1)
        # * 3=p2.vertex(e2+2)
        # Record the corresponding vertices of this quadrilateral.
        q1 = (3 + v1 - e1 - 1) % 3
        q2 = (2 + (3 + v2 - e2 - 1) % 3) % 4

        new_diagonal = p2.vertex((e2 + 2) % 3) - p1.vertex((e1 + 2) % 3)
        # This list will store the new triangles which are being glued in.
        # (Unfortunately, they may not be cyclically labeled in the correct way.)
        new_triangle = []
        try:
            new_triangle.append(
                Polygon(
                    edges=[
                        p1.edge((e1 + 2) % 3),
                        p2.edge((e2 + 1) % 3),
                        -new_diagonal,
                    ],
                    base_ring=p1.base_ring(),
                )
            )
            new_triangle.append(
                Polygon(
                    edges=[
                        p2.edge((e2 + 2) % 3),
                        p1.edge((e1 + 1) % 3),
                        new_diagonal,
                    ],
                    base_ring=p1.base_ring(),
                )
            )
            # The above triangles would be glued along edge 2 to form the diagonal of the quadrilateral being removed.
        except ValueError:
            raise ValueError(
                "Gluing triangles along this edge yields a non-convex quadrilateral."
            )

        # Find the separatrices of the two new triangles, and in particular which way they point.
        new_sep = []
        new_sep.append(new_triangle[0].find_separatrix(direction=direction)[0])
        new_sep.append(new_triangle[1].find_separatrix(direction=direction)[0])
        # The quadrilateral vertices corresponding to these separatrices are
        # new_sep[0]+1 and (new_sep[1]+3)%4 respectively.

        # i=0 if the new_triangle[0] should be labeled l1 and new_triangle[1] should be labeled l2.
        # i=1 indicates the opposite labeling.
        if new_sep[0] + 1 == q1:
            assert (new_sep[1] + 3) % 4 == q2
            i = 0
        else:
            assert (new_sep[1] + 3) % 4 == q1
            assert new_sep[0] + 1 == q2
            i = 1

        # These quantities represent the cyclic relabeling of triangles needed.
        cycle1 = (new_sep[i] - v1 + 3) % 3
        cycle2 = (new_sep[1 - i] - v2 + 3) % 3

        # This will be the new triangle with label l1:
        tri1 = Polygon(
            edges=[
                new_triangle[i].edge(cycle1),
                new_triangle[i].edge((cycle1 + 1) % 3),
                new_triangle[i].edge((cycle1 + 2) % 3),
            ],
            base_ring=p1.base_ring(),
        )
        # This will be the new triangle with label l2:
        tri2 = Polygon(
            edges=[
                new_triangle[1 - i].edge(cycle2),
                new_triangle[1 - i].edge((cycle2 + 1) % 3),
                new_triangle[1 - i].edge((cycle2 + 2) % 3),
            ],
            base_ring=p1.base_ring(),
        )
        # In the above, edge 2-cycle1 of tri1 would be glued to edge 2-cycle2 of tri2
        diagonal_glue_e1 = 2 - cycle1
        diagonal_glue_e2 = 2 - cycle2

        assert p1.find_separatrix(direction=direction) == tri1.find_separatrix(
            direction=direction
        )
        assert p2.find_separatrix(direction=direction) == tri2.find_separatrix(
            direction=direction
        )

        # Two opposite edges will not change their labels (label,edge) under our regluing operation.
        # The other two opposite ones will change and in fact they change labels.
        # The following finds them (there are two cases).
        # At the end of the if statement, the following will be true:
        # * new_glue_e1 and new_glue_e2 will be the edges of the new triangle with label l1 and l2 which need regluing.
        # * old_e1 and old_e2 will be the corresponding edges of the old triangles.
        # (Note that labels are swapped between the pair. The appending 1 or 2 refers to the label used for the triangle.)
        if p1.edge(v1) == tri1.edge(v1):
            # We don't have to worry about changing gluings on edge v1 of the triangles with label l1
            # We do have to worry about the following edge:
            new_glue_e1 = (
                3 - diagonal_glue_e1 - v1
            )  # returns the edge which is neither diagonal_glue_e1 nor v1.
            # This corresponded to the following old edge:
            old_e1 = 3 - e1 - v1  # Again this finds the edge which is neither e1 nor v1
        else:
            temp = (v1 + 2) % 3
            assert p1.edge(temp) == tri1.edge(temp)
            # We don't have to worry about changing gluings on edge (v1+2)%3 of the triangles with label l1
            # We do have to worry about the following edge:
            new_glue_e1 = (
                3 - diagonal_glue_e1 - temp
            )  # returns the edge which is neither diagonal_glue_e1 nor temp.
            # This corresponded to the following old edge:
            old_e1 = (
                3 - e1 - temp
            )  # Again this finds the edge which is neither e1 nor temp
        if p2.edge(v2) == tri2.edge(v2):
            # We don't have to worry about changing gluings on edge v2 of the triangles with label l2
            # We do have to worry about the following edge:
            new_glue_e2 = (
                3 - diagonal_glue_e2 - v2
            )  # returns the edge which is neither diagonal_glue_e2 nor v2.
            # This corresponded to the following old edge:
            old_e2 = 3 - e2 - v2  # Again this finds the edge which is neither e2 nor v2
        else:
            temp = (v2 + 2) % 3
            assert p2.edge(temp) == tri2.edge(temp)
            # We don't have to worry about changing gluings on edge (v2+2)%3 of the triangles with label l2
            # We do have to worry about the following edge:
            new_glue_e2 = (
                3 - diagonal_glue_e2 - temp
            )  # returns the edge which is neither diagonal_glue_e2 nor temp.
            # This corresponded to the following old edge:
            old_e2 = (
                3 - e2 - temp
            )  # Again this finds the edge which is neither e2 nor temp

        # remember the old gluings.
        old_opposite1 = s.opposite_edge(l1, old_e1)
        old_opposite2 = s.opposite_edge(l2, old_e2)

        us = s

        # Replace the triangles.
        us.replace_polygon(l1, tri1)
        us.replace_polygon(l2, tri2)
        # Glue along the new diagonal of the quadrilateral
        us.glue((l1, diagonal_glue_e1), (l2, diagonal_glue_e2))
        # Now we deal with that pair of opposite edges of the quadrilateral that need regluing.
        # There are some special cases:
        if old_opposite1 == (l2, old_e2):
            # These opposite edges were glued to each other.
            # Do the same in the new surface:
            us.glue((l1, new_glue_e1), (l2, new_glue_e2))
        else:
            if old_opposite1 == (l1, old_e1):
                # That edge was "self-glued".
                us.glue((l2, new_glue_e2), (l2, new_glue_e2))
            else:
                # The edge (l1,old_e1) was glued in a standard way.
                # That edge now corresponds to (l2,new_glue_e2):
                us.glue((l2, new_glue_e2), (old_opposite1[0], old_opposite1[1]))
            if old_opposite2 == (l2, old_e2):
                # That edge was "self-glued".
                us.glue((l1, new_glue_e1), (l1, new_glue_e1))
            else:
                # The edge (l2,old_e2) was glued in a standard way.
                # That edge now corresponds to (l1,new_glue_e1):
                us.glue((l1, new_glue_e1), (old_opposite2[0], old_opposite2[1]))
        return s

    def standardize_polygons(self, in_place=False):
        r"""
        Replace each polygon with a new polygon which differs by
        translation and reindexing. The new polygon will have the property
        that vertex zero is the origin, and all vertices lie either in the
        upper half plane, or on the x-axis with non-negative x-coordinate.

        This is done to the current surface if in_place=True. A mutable
        copy is created and returned if in_place=False (as default).

        This overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.FiniteType.Oriented.ParentMethods.standardize_polygons`
        to provide in-place standardizing of surfaces.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
            sage: p = Polygon(vertices = ([(1,1),(2,1),(2,2),(1,2)]))
            sage: s = MutableOrientedSimilaritySurface(QQ)
            sage: s.add_polygon(p)
            0
            sage: s.glue((0, 0), (0, 2))
            sage: s.glue((0, 1), (0, 3))
            sage: s.set_root(0)
            sage: s.set_immutable()

            sage: s.standardize_polygons().polygon(0)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        if not in_place:
            S = MutableOrientedSimilaritySurface.from_surface(self)
            S.standardize_polygons(in_place=True)
            return S

        cv = {}  # dictionary for non-zero canonical vertices
        for label, polygon in zip(self.labels(), self.polygons()):
            best = 0
            best_pt = polygon.vertex(best)
            for v in range(1, len(polygon.vertices())):
                pt = polygon.vertex(v)
                if (pt[1] < best_pt[1]) or (pt[1] == best_pt[1] and pt[0] < best_pt[0]):
                    best = v
                    best_pt = pt
            # We replace the polygon if the best vertex is not the zero vertex, or
            # if the coordinates of the best vertex differs from the origin.
            if not (best == 0 and best_pt.is_zero()):
                cv[label] = best
        for label, v in cv.items():
            self.set_vertex_zero(label, v, in_place=True)

        return self


class MutableOrientedSimilaritySurface(
    MutableOrientedSimilaritySurface_base, MutablePolygonalSurface
):
    r"""
    A surface built from Euclidean polyogns glued by orientation preserving
    similarities.

    This is the main tool to create new surfaces of finite type in
    sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface
        sage: S = MutableOrientedSimilaritySurface(QQ)

        sage: from flatsurf import polygons
        sage: S.add_polygon(polygons.square())
        0

        sage: S.glue((0, 0), (0, 2))
        sage: S.glue((0, 1), (0, 3))

        sage: S.set_immutable()

        sage: S
        Translation Surface in H_1(0) built from a square

    TESTS::

        sage: TestSuite(S).run()

    """

    def __init__(self, base, category=None):
        self._gluings = {}

        from flatsurf.geometry.categories import SimilaritySurfaces

        if category is None:
            category = SimilaritySurfaces().Oriented().FiniteType()

        category &= SimilaritySurfaces().Oriented().FiniteType()

        super().__init__(base, category=category)

    @classmethod
    def from_surface(cls, surface, category=None):
        r"""
        Return a mutable copy of ``surface``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface, polygons

            sage: T = translation_surfaces.square_torus()

        We build a surface that is made from two tori:

            sage: S = MutableOrientedSimilaritySurface.from_surface(T)
            sage: S.add_polygon(polygons.square())
            1
            sage: S.glue((1, 0), (1, 2))
            sage: S.glue((1, 1), (1, 3))

            sage: S.set_immutable()

            sage: S
            Disconnected Surface built from 2 squares

        """
        if not surface.is_finite_type():
            raise TypeError
        self = MutableOrientedSimilaritySurface(surface.base_ring(), category=category)

        for label in surface.labels():
            self.add_polygon(surface.polygon(label), label=label)

        for label in surface.labels():
            for edge in range(len(surface.polygon(label).vertices())):
                cross = surface.opposite_edge(label, edge)
                if cross:
                    self.glue((label, edge), cross)

        if isinstance(surface, MutablePolygonalSurface):
            # Only copy explicitly set roots over
            self._roots = surface._roots
        else:
            self.set_roots(surface.roots())

        return self

    def add_polygon(self, polygon, *, label=None):
        # Overrides add_polygon from MutablePolygonalSurface
        label = super().add_polygon(polygon, label=label)
        assert label not in self._gluings
        self._gluings[label] = [None] * len(polygon.vertices())

        if self._roots:
            self._roots = self._roots + (label,)

        return label

    def remove_polygon(self, label):
        # Overrides remove_polygon from MutablePolygonalSurface
        self._unglue_polygon(label)
        self._gluings.pop(label)

        super().remove_polygon(label)

    def glue(self, x, y):
        r"""
        Glue ``x`` and ``y`` with an (orientation preserving) similarity.

        INPUT:

        - ``x`` -- a pair consisting of a polygon label and an edge index for
          that polygon

        - ``y`` -- a pair consisting of a polygon label and an edge index for
          that polygon

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface, polygons

            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygons.square())
            0

        Glue two opposite sides of the square to each other::

            sage: S.glue((0, 1), (0, 3))

        Glue the other sides of the square to themselves::

            sage: S.glue((0, 0), (0, 0))
            sage: S.glue((0, 2), (0, 2))

        Note that existing gluings are removed when gluing already glued
        sides::

            sage: S.glue((0, 0), (0, 2))
            sage: S.set_immutable()

            sage: S
            Translation Surface in H_1(0) built from a square

        """
        if not self._mutable:
            raise Exception(
                "cannot modify immutable surface; create a copy with MutableOrientedSimilaritySurface.from_surface()"
            )

        if x[0] not in self._polygons:
            raise ValueError

        if y[0] not in self._polygons:
            raise ValueError

        self.unglue(*x)
        self.unglue(*y)

        if self._roots:
            component = set(self.component(x[0]))
            if y[0] not in component:
                # Gluing will join two connected components.
                cross_component = set(self.component(y[0]))
                for root in reversed(self._roots):
                    if root in component or root in cross_component:
                        self._roots = tuple([r for r in self._roots if r != root])
                        break
                else:
                    assert False, "did not find any root to eliminate"

        self._gluings[x[0]][x[1]] = y
        self._gluings[y[0]][y[1]] = x

    def unglue(self, label, edge):
        r"""
        Unglue the side ``edge`` of the polygon ``label`` if it is glued.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface, translation_surfaces

            sage: T = translation_surfaces.square_torus()

            sage: S = MutableOrientedSimilaritySurface.from_surface(T)

            sage: S.unglue(0, 0)

            sage: S.gluings()
            (((0, 1), (0, 3)), ((0, 3), (0, 1)))

            sage: S.set_immutable()
            sage: S
            Translation Surface with boundary built from a square

        """
        if not self._mutable:
            raise Exception(
                "cannot modify immutable surface; create a copy with MutableOrientedSimilaritySurface.from_surface()"
            )

        cross = self._gluings[label][edge]
        if cross is not None:
            self._gluings[cross[0]][cross[1]] = None

        self._gluings[label][edge] = None

        if cross is not None and self._roots:
            component = set(self.component(label))
            if cross[0] not in component:
                # Ungluing created a new connected component.
                cross_component = set(self.component(cross[0]))
                assert label not in cross_component
                for root in self._roots:
                    if root in component:
                        self._roots = self._roots + (
                            LabeledView(
                                surface=self, view=cross_component, finite=True
                            ).min(),
                        )
                        break
                    if root in cross_component:
                        self._roots = self._roots + (
                            LabeledView(
                                surface=self, view=component, finite=True
                            ).min(),
                        )
                        break
                else:
                    assert False, "did not find any root to split"

    def _unglue_polygon(self, label):
        r"""
        Remove all gluigns from polygon ``label``.

        This is a helper method to completely unglue a polygon before removing
        or replacing it.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface, translation_surfaces

            sage: T = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(T)

            sage: S._unglue_polygon(0)
            sage: S.gluings()
            ()

        """
        for edge, cross in enumerate(self._gluings[label]):
            if cross is None:
                continue
            cross_label, cross_edge = cross
            self._gluings[cross_label][cross_edge] = None
        self._gluings[label] = [None] * len(self.polygon(label).vertices())

    def set_edge_pairing(self, label0, edge0, label1, edge1):
        r"""
        TESTS::

            sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygons.square())
            0
            sage: S.set_edge_pairing(0, 0, 0, 2)
            doctest:warning
            ...
            UserWarning: set_edge_pairing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead
            sage: S.set_edge_pairing(0, 1, 0, 3)

            sage: S.gluings()
            (((0, 0), (0, 2)), ((0, 1), (0, 3)), ((0, 2), (0, 0)), ((0, 3), (0, 1)))

        """
        import warnings

        warnings.warn(
            "set_edge_pairing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead"
        )
        return self.glue((label0, edge0), (label1, edge1))

    def change_edge_gluing(self, label0, edge0, label1, edge1):
        r"""
        TESTS::

            sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygons.square())
            0
            sage: S.change_edge_gluing(0, 0, 0, 2)
            doctest:warning
            ...
            UserWarning: change_edge_gluing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead
            sage: S.change_edge_gluing(0, 1, 0, 3)

            sage: S.gluings()
            (((0, 0), (0, 2)), ((0, 1), (0, 3)), ((0, 2), (0, 0)), ((0, 3), (0, 1)))

        """
        import warnings

        warnings.warn(
            "change_edge_gluing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead"
        )
        return self.glue((label0, edge0), (label1, edge1))

    def change_polygon_gluings(self, label, gluings):
        r"""
        TESTS::

            sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygons.square())
            0
            sage: S.change_polygon_gluings(0, [(0, 2), (0, 3), (0, 0), (0, 1)])
            doctest:warning
            ...
            UserWarning: change_polygon_gluings() has been deprecated and will be removed in a future version of sage-flatsurf; use glue() in a loop instead

            sage: S.gluings()
            (((0, 0), (0, 2)), ((0, 1), (0, 3)), ((0, 2), (0, 0)), ((0, 3), (0, 1)))

        """
        import warnings

        warnings.warn(
            "change_polygon_gluings() has been deprecated and will be removed in a future version of sage-flatsurf; use glue() in a loop instead"
        )

        for edge0, cross in enumerate(gluings):
            if cross is None:
                self.unglue(label, edge0)
            else:
                self.glue((label, edge0), cross)

    def change_polygon(self, label, polygon, gluing_list=None):
        r"""
        TESTS::

            sage: from flatsurf import MutableOrientedSimilaritySurface, translation_surfaces

            sage: T = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(T)

            sage: S.change_polygon(0, 2 * S.polygon(0))
            doctest:warning
            ...
            UserWarning: change_polygon() has been deprecated and will be removed in a future version of sage-flatsurf; use replace_polygon() or remove_polygon() and add_polygon() instead

            sage: S.polygon(0)
            Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)])

        """
        import warnings

        warnings.warn(
            "change_polygon() has been deprecated and will be removed in a future version of sage-flatsurf; use replace_polygon() or remove_polygon() and add_polygon() instead"
        )

        if not self._mutable:
            raise Exception(
                "cannot modify immutable surface; create a copy with MutableOrientedSimilaritySurface.from_surface()"
            )

        # Note that this obscure feature. If the number of edges is unchanged, we keep the gluings, otherwise we trash them all.
        if len(polygon.vertices()) != len(self.polygon(label).vertices()):
            self._unglue_polygon(label)
            self._gluings[label] = [None] * len(polygon.vertices())

        self._polygons[label] = polygon

        if gluing_list is not None:
            for i, cross in enumerate(gluing_list):
                self.glue((label, i), cross)

    def replace_polygon(self, label, polygon):
        r"""
        Replace the polygon ``label`` with ``polygon`` while keeping its
        gluings intact.

        INPUT:

        - ``label`` -- an element of :meth:`~.MutablePolygonalSurface.labels`

        - ``polygon`` -- a Euclidean polygon

        EXAMPLES::

            sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))

            sage: S.replace_polygon(0, Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)]))

        The replacement of a polygon must have the same number of sides::

            sage: S.replace_polygon(0, Polygon(vertices=[(0, 0), (2, 0), (2, 2)]))
            Traceback (most recent call last):
            ...
            ValueError: polygon must be a quadrilateral

        To replace the polygon without keeping its glueings, remove the polygon
        first and then add a new one::

            sage: S.remove_polygon(0)
            sage: S.add_polygon(Polygon(vertices=[(0, 0), (2, 0), (2, 2)]), label=0)
            0

        """
        old = self.polygon(label)

        if len(old.vertices()) != len(polygon.vertices()):
            from flatsurf.geometry.categories.polygons import Polygons

            article, singular, plural = Polygons._describe_polygon(len(old.vertices()))
            raise ValueError(f"polygon must be {article} {singular}")

        self._polygons[label] = polygon

    def opposite_edge(self, label, edge=None):
        r"""
        Return the edge that ``edge`` of ``label`` is glued to or ``None`` if this edge is unglued.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        INPUT:

        - ``label`` -- one of the labels included in :meth:`~.MutablePolygonalSurface.labels`

        - ``edge`` -- a non-negative integer to specify an edge (the edges
          of a polygon are numbered starting from zero.)

        EXAMPLES::

            sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0

            sage: S.glue((0, 0), (0, 1))
            sage: S.glue((0, 2), (0, 2))

            sage: S.opposite_edge(0, 0)
            (0, 1)
            sage: S.opposite_edge(0, 1)
            (0, 0)
            sage: S.opposite_edge(0, 2)
            (0, 2)
            sage: S.opposite_edge(0, 3)

            sage: S.opposite_edge((0, 0))
            doctest:warning
            ...
            UserWarning: calling opposite_edge() with a single argument has been deprecated and will be removed in a future version of sage-flatsurf; use opposite_edge(label, edge) instead
            (0, 1)

        """
        if edge is None:
            import warnings

            warnings.warn(
                "calling opposite_edge() with a single argument has been deprecated and will be removed in a future version of sage-flatsurf; use opposite_edge(label, edge) instead"
            )
            label, edge = label
        return self._gluings[label][edge]

    def set_vertex_zero(self, label, v, in_place=False):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.set_vertex_zero`
        to make it possible to set the zero vertex in-place.
        """
        if not in_place:
            return super().set_vertex_zero(label, v, in_place=in_place)

        us = self
        if not us.is_mutable():
            raise ValueError(
                "set_vertex_zero can only be done in_place for a mutable surface."
            )
        p = us.polygon(label)
        n = len(p.vertices())
        if not (0 <= v < n):
            raise ValueError
        glue = []

        from flatsurf import Polygon

        pp = Polygon(
            edges=[p.edge((i + v) % n) for i in range(n)], base_ring=us.base_ring()
        )

        for i in range(n):
            e = (v + i) % n
            ll, ee = us.opposite_edge(label, e)
            if ll == label:
                ee = (ee + n - v) % n
            glue.append((ll, ee))

        us.remove_polygon(label)
        us.add_polygon(pp, label=label)
        for e, cross in enumerate(glue):
            us.glue((label, e), cross)
        return self

    def relabel(self, relabeling_map, in_place=False):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.relabel`
        to allow relabeling in-place.
        """
        if not in_place:
            return super().relabel(relabeling_map=relabeling_map, in_place=in_place)

        us = self
        if not isinstance(relabeling_map, dict):
            raise NotImplementedError(
                "Currently relabeling is only implemented via a dictionary."
            )
        domain = set()
        codomain = set()
        data = {}
        for l1, l2 in relabeling_map.items():
            p = us.polygon(l1)
            glue = []
            for e in range(len(p.vertices())):
                ll, ee = us.opposite_edge(l1, e)
                try:
                    lll = relabeling_map[ll]
                except KeyError:
                    lll = ll
                glue.append((lll, ee))
            data[l2] = (p, glue)
            domain.add(l1)
            codomain.add(l2)
        if len(domain) != len(codomain):
            raise ValueError(
                "The relabeling_map must be injective. Received " + str(relabeling_map)
            )
        changed_labels = domain.intersection(codomain)
        added_labels = codomain.difference(domain)
        removed_labels = domain.difference(codomain)
        # Pass to add_polygons
        roots = list(us.roots())
        relabel_errors = {}
        for l2 in added_labels:
            p, glue = data[l2]
            l3 = us.add_polygon(p, label=l2)
            if not l2 == l3:
                # This means the label l2 could not be added for some reason.
                # Perhaps the implementation does not support this type of label.
                # Or perhaps there is already a polygon with this label.
                relabel_errors[l2] = l3
        # Pass to change polygons
        for l2 in changed_labels:
            p, glue = data[l2]
            us.remove_polygon(l2)
            us.add_polygon(p, label=l2)
            us.replace_polygon(l2, p)
        # Deal with the component roots
        roots = [relabeling_map.get(label, label) for label in roots]
        roots = [relabel_errors.get(label, label) for label in roots]
        # Pass to remove polygons:
        for l1 in removed_labels:
            us.remove_polygon(l1)
        # Pass to update the edge gluings
        if len(relabel_errors) == 0:
            # No problems. Update the gluings.
            for l2 in codomain:
                p, glue = data[l2]
                for e, cross in enumerate(glue):
                    us.glue((l2, e), cross)
        else:
            # Use the gluings provided by relabel_errors when necessary
            for l2 in codomain:
                p, glue = data[l2]
                for e in range(len(p.vertices())):
                    ll, ee = glue[e]
                    try:
                        # First try the error dictionary
                        us.glue((l2, e), (relabel_errors[ll], ee))
                    except KeyError:
                        us.glue((l2, e), (ll, ee))
        us.set_roots(roots)
        return self, len(relabel_errors) == 0

    def join_polygons(self, p1, e1, test=False, in_place=False):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.join_polygons`
        to allow joining in-place.
        """
        if test:
            in_place = False

        if not in_place:
            return super().join_polygon(p1, e1, test=test, in_place=in_place)

        poly1 = self.polygon(p1)
        p2, e2 = self.opposite_edge(p1, e1)
        poly2 = self.polygon(p2)
        if p1 == p2:
            raise ValueError("Can't glue polygon to itself.")
        t = self.edge_transformation(p2, e2)
        dt = t.derivative()
        es = []
        edge_map = {}  # Store the pairs for the old edges.
        for i in range(e1):
            edge_map[len(es)] = (p1, i)
            es.append(poly1.edge(i))
        ne = len(poly2.vertices())
        for i in range(1, ne):
            ee = (e2 + i) % ne
            edge_map[len(es)] = (p2, ee)
            es.append(dt * poly2.edge(ee))
        for i in range(e1 + 1, len(poly1.vertices())):
            edge_map[len(es)] = (p1, i)
            es.append(poly1.edge(i))

        from flatsurf import Polygon

        new_polygon = Polygon(edges=es, base_ring=self.base_ring())

        # Do the gluing.
        ss = self
        s = ss

        inv_edge_map = {}
        for key, value in edge_map.items():
            inv_edge_map[value] = (p1, key)

        glue_list = []
        for i in range(len(es)):
            p3, e3 = edge_map[i]
            p4, e4 = self.opposite_edge(p3, e3)
            if p4 == p1 or p4 == p2:
                glue_list.append(inv_edge_map[(p4, e4)])
            else:
                glue_list.append((p4, e4))

        if p2 in s.roots():
            s.set_roots((p1 if label == p2 else label for label in s.roots()))

        s.remove_polygon(p2)

        s.remove_polygon(p1)
        s.add_polygon(new_polygon, label=p1)
        for e, cross in enumerate(glue_list):
            s.glue((p1, e), cross)

        return s

    def subdivide_polygon(self, p, v1, v2, test=False, new_label=None):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.subdivide_polygon`
        to allow subdividing in-place.
        """
        if test:
            return super().subdivide_polygon(
                p=p, v1=v1, v2=v2, test=test, new_label=new_label
            )

        poly = self.polygon(p)
        ne = len(poly.vertices())
        if v1 < 0 or v2 < 0 or v1 >= ne or v2 >= ne:
            raise ValueError("Provided vertices out of bounds.")
        if abs(v1 - v2) <= 1 or abs(v1 - v2) >= ne - 1:
            raise ValueError("Provided diagonal is not actually a diagonal.")

        if v2 < v1:
            v2 = v2 + ne

        newedges1 = [poly.vertex(v2) - poly.vertex(v1)]
        for i in range(v2, v1 + ne):
            newedges1.append(poly.edge(i))

        from flatsurf import Polygon

        newpoly1 = Polygon(edges=newedges1, base_ring=self.base_ring())

        newedges2 = [poly.vertex(v1) - poly.vertex(v2)]
        for i in range(v1, v2):
            newedges2.append(poly.edge(i))
        newpoly2 = Polygon(edges=newedges2, base_ring=self.base_ring())

        # Store the old gluings
        old_gluings = {(p, i): self.opposite_edge(p, i) for i in range(ne)}

        # Update the polygon with label p, add a new polygon.
        self.remove_polygon(p)
        self.add_polygon(newpoly1, label=p)
        if new_label is None:
            new_label = self.add_polygon(newpoly2)
        else:
            new_label = self.add_polygon(newpoly2, label=new_label)
        # This gluing is the diagonal we used.
        self.glue((p, 0), (new_label, 0))

        # Setup conversion from old to new labels.
        old_to_new_labels = {}
        for i in range(v1, v2):
            old_to_new_labels[(p, i % ne)] = (new_label, i - v1 + 1)
        for i in range(v2, ne + v1):
            old_to_new_labels[(p, i % ne)] = (p, i - v2 + 1)

        for e in range(1, len(newpoly1.vertices())):
            pair = old_gluings[(p, (v2 + e - 1) % ne)]
            if pair in old_to_new_labels:
                pair = old_to_new_labels[pair]
            self.glue((p, e), (pair[0], pair[1]))

        for e in range(1, len(newpoly2.vertices())):
            pair = old_gluings[(p, (v1 + e - 1) % ne)]
            if pair in old_to_new_labels:
                pair = old_to_new_labels[pair]
            self.glue((new_label, e), (pair[0], pair[1]))

    def reposition_polygons(self, in_place=False, relabel=None):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.FiniteType.Oriented.ParentMethods.reposition_polygons`
        to allow normalizing in-place.
        """
        if not in_place:
            return super().reposition_polygons(in_place=in_place, relabel=relabel)

        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from reposition_polygon; use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to reposition_polygons()"
                )

        s = self

        labels = list(s.labels())
        from flatsurf.geometry.similarity import SimilarityGroup

        S = SimilarityGroup(self.base_ring())
        identity = S.one()
        it = iter(labels)
        label = next(it)
        changes = {label: identity}
        for label in it:
            polygon = self.polygon(label)
            adjacencies = {
                edge: self.opposite_edge(label, edge)[0]
                for edge in range(len(polygon.vertices()))
            }
            edge = min(
                adjacencies,
                # pylint: disable-next=cell-var-from-loop
                key=lambda edge: labels.index(adjacencies[edge]),
            )
            label2, edge2 = s.opposite_edge(label, edge)
            changes[label] = changes[label2] * s.edge_transformation(label, edge)
        it = iter(labels)
        # Skip the base label:
        label = next(it)
        for label in it:
            p = s.polygon(label)
            p = changes[label].derivative() * p
            s.replace_polygon(label, p)
        return s

    def triangulate(self, in_place=False, label=None, relabel=None):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.triangulate`
        to allow triangulating in-place.
        """
        if relabel is not None:
            import warnings

            warnings.warn(
                "the relabel keyword argument of triangulate() is ignored, it has been deprecated and will be removed in a future version of sage-flatsurf"
            )

        if not in_place:
            return super().triangulate(in_place=in_place, label=label)

        if label is None:
            # We triangulate the whole surface
            # Store the current labels.
            labels = [label for label in self.labels()]
            s = self
            # Subdivide each polygon in turn.
            for label in labels:
                s = s.triangulate(in_place=True, label=label)
            return s

        poly = self.polygon(label)
        n = len(poly.vertices())
        if n > 3:
            s = self
        else:
            # This polygon is already a triangle.
            return self
        from flatsurf.geometry.euclidean import ccw

        for i in range(n - 3):
            poly = s.polygon(label)
            n = len(poly.vertices())
            for i in range(n):
                e1 = poly.edge(i)
                e2 = poly.edge((i + 1) % n)
                if ccw(e1, e2) != 0:
                    # This is in case the polygon is a triangle with subdivided edge.
                    e3 = poly.edge((i + 2) % n)
                    if ccw(e1 + e2, e3) != 0:
                        s.subdivide_polygon(label, i, (i + 2) % n)
                        break
        return s

    def delaunay_single_flip(self):
        r"""
        Perform a single in place flip of a triangulated mutable surface
        in-place.
        """
        lc = self._label_comparator()
        for (l1, e1), (l2, e2) in self.gluings():
            if (
                lc.lt(l1, l2) or (l1 == l2 and e1 <= e2)
            ) and self._delaunay_edge_needs_flip(l1, e1):
                self.triangle_flip(l1, e1, in_place=True)
                return True
        return False

    def delaunay_triangulation(
        self,
        triangulated=False,
        in_place=False,
        direction=None,
        relabel=None,
    ):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.delaunay_triangulation`
        to allow triangulating in-place.
        """
        if not in_place:
            return super().delaunay_triangulation(
                triangulated=triangulated,
                in_place=in_place,
                direction=direction,
                relabel=relabel,
            )

        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from delaunay_triangulation(); use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_triangulation()"
                )

        if triangulated:
            s = self
        else:
            s = self
            self.triangulate(in_place=True)

        if direction is None:
            direction = (self.base_ring() ** 2)((0, 1))

        if direction.is_zero():
            raise ValueError

        from collections import deque

        unchecked_labels = deque(s.labels())
        checked_labels = set()
        while unchecked_labels:
            label = unchecked_labels.popleft()
            flipped = False
            for edge in range(3):
                if s._delaunay_edge_needs_flip(label, edge):
                    # Record the current opposite edge:
                    label2, edge2 = s.opposite_edge(label, edge)
                    # Perform the flip.
                    s.triangle_flip(label, edge, in_place=True, direction=direction)
                    # Move the opposite polygon to the list of labels we need to check.
                    if label2 != label:
                        try:
                            checked_labels.remove(label2)
                            unchecked_labels.append(label2)
                        except KeyError:
                            # Occurs if label2 is not in checked_labels
                            pass
                    flipped = True
                    break
            if flipped:
                unchecked_labels.append(label)
            else:
                checked_labels.add(label)
        return s

    def delaunay_decomposition(
        self,
        triangulated=False,
        delaunay_triangulated=False,
        in_place=False,
        direction=None,
        relabel=None,
    ):
        r"""
        Overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.delaunay_decomposition`
        to allow normalizing in-place.
        """
        if not in_place:
            return super().delaunay_decomposition(
                triangulated=triangulated,
                delaunay_triangulated=delaunay_triangulated,
                in_place=in_place,
                direction=direction,
                relabel=relabel,
            )

        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from delaunay_decomposition(); use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_decomposition()"
                )

        s = self
        if not delaunay_triangulated:
            s = s.delaunay_triangulation(
                triangulated=triangulated,
                in_place=True,
                direction=direction,
                relabel=relabel,
            )

        while True:
            for (l1, e1), (l2, e2) in s.gluings():
                if s._delaunay_edge_needs_join(l1, e1):
                    s.join_polygons(l1, e1, in_place=True)
                    break
            else:
                return s

    def cmp(self, s2, limit=None):
        r"""
        Compare two surfaces. This is an ordering returning -1, 0, or 1.

        The surfaces will be considered equal if and only if there is a translation automorphism
        respecting the polygons and the root labels.

        If the two surfaces are infinite, we just examine the first limit polygons.
        """
        if self.is_finite_type():
            if s2.is_finite_type():
                if limit is not None:
                    raise ValueError("limit only enabled for finite surfaces")

                sign = len(self.polygons()) - len(s2.polygons())
                if sign > 0:
                    return 1
                if sign < 0:
                    return -1

                lw1 = self.labels()
                labels1 = list(lw1)

                lw2 = s2.labels()
                labels2 = list(lw2)

                for l1, l2 in zip(lw1, lw2):
                    ret = self.polygon(l1).cmp(self.polygon(l2))
                    if ret != 0:
                        return ret

                    for e in range(len(self.polygon(l1).vertices())):
                        ll1, e1 = self.opposite_edge(l1, e)
                        ll2, e2 = s2.opposite_edge(l2, e)
                        num1 = labels1.index(ll1)
                        num2 = labels2.index(ll2)

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
            if s2.is_finite_type():
                # s1 is infinite but s2 is finite.
                return 1
            else:
                if limit is None:
                    raise NotImplementedError

                # both surfaces are infinite.
                from itertools import islice

                lw1 = self.labels()
                labels1 = list(islice(lw1, limit))

                lw2 = s2.labels()
                labels2 = list(islice(lw2, limit))

                count = 0
                for l1, l2 in zip(lw1, lw2):
                    ret = self.polygon(l1).cmp(s2.polygon(l2))
                    if ret != 0:
                        return ret

                    for e in range(len(self.polygon(l1).vertices())):
                        ll1, ee1 = self.opposite_edge(l1, e)
                        ll2, ee2 = s2.opposite_edge(l2, e)
                        num1 = labels1.index(ll1)
                        num2 = labels2.index(ll2)
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

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See
        :meth:`~.categories.similarity_surfaces.SimilaritySurfaces.FiniteType.ParentMethods._test_eq_surface`
        for details on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface

            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: T = MutableOrientedSimilaritySurface(AA)

            sage: S == T
            False

            sage: S == S
            True

        """
        if not isinstance(other, MutableOrientedSimilaritySurface):
            return False

        if not super().__eq__(other):
            return False

        if self._gluings != other._gluings:
            return False

        return True

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface

            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: T = MutableOrientedSimilaritySurface(QQ)

            sage: hash(S) == hash(T)
            Traceback (most recent call last):
            ...
            TypeError: cannot hash a mutable surface

            sage: S.set_immutable()
            sage: T.set_immutable()
            sage: hash(S) == hash(T)
            True

        """
        if self._mutable:
            raise TypeError("cannot hash a mutable surface")

        return hash((super().__hash__(), tuple(self.gluings())))


class BaseRingChangedSurface(OrientedSimilaritySurface):
    r"""
    Changes the ring over which a surface is defined.

    EXAMPLES:

    This class is used in the implementation of
    :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.Oriented.ParentMethods.change_ring`::

        sage: from flatsurf import translation_surfaces
        sage: T = translation_surfaces.square_torus()
        sage: S = T.change_ring(AA)

        sage: from flatsurf.geometry.surface import BaseRingChangedSurface
        sage: isinstance(S, BaseRingChangedSurface)
        True

        sage: TestSuite(S).run()

    """

    def __init__(self, surface, ring, category=None):
        if surface.is_mutable():
            raise NotImplementedError("surface must be immutable")

        self._reference = surface
        super().__init__(ring, category=category or surface.category())

    def is_mutable(self):
        r"""
        Return whether this surface can be modified, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: S = T.change_ring(AA)

            sage: S.is_mutable()
            False

        """
        return False

    def labels(self):
        r"""
        Return the labels of the polygons of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: S = T.change_ring(AA)
            sage: S.labels()
            (0,)

        """
        return self._reference.labels()

    def roots(self):
        r"""
        Return a label for each connected component on this surface.

        This implements :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: S = T.change_ring(AA)

            sage: S.roots()
            (0,)

        """
        return self._reference.roots()

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: S = T.change_ring(AA)

            sage: p = S.polygon(0)
            sage: p.base_ring()
            Algebraic Real Field

        """
        return self._reference.polygon(label).change_ring(self.base_ring())

    def opposite_edge(self, label, edge):
        r"""
        Return the edge that ``edge`` of ``label`` is glued to or ``None`` if this edge is unglued.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: S = T.change_ring(AA)

            sage: S.opposite_edge(0, 0)
            (0, 2)

        """
        return self._reference.opposite_edge(label, edge)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See
        :meth:`~.categories.similarity_surfaces.SimilaritySurfaces.FiniteType.ParentMethods._test_eq_surface`
        for details on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: T.change_ring(AA) == T.change_ring(AA)
            True

        """
        if not isinstance(other, BaseRingChangedSurface):
            return False

        return self._reference == other._reference and self.base() == other.base()

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: hash(T.change_ring(AA)) == hash(T.change_ring(AA))
            True

        """
        return hash((self._reference, self.base()))


class RootedComponents_MutablePolygonalSurface(collections.abc.Mapping):
    r"""
    Connected components of a :class:`MutablePolygonalSurface`.

    The components are represented as a mapping that maps the root labels to
    the labels of the corresponding component.

    This is a helper method for :meth:`MutablePolygonalSurface.components` and
    :meth:`MutablePolygonalSurface.roots`.

    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface
        sage: S = MutableOrientedSimilaritySurface(QQ)

        sage: from flatsurf import polygons
        sage: S.add_polygon(polygons.square())
        0
        sage: S.add_polygon(polygons.square())
        1

        sage: from flatsurf.geometry.surface import RootedComponents_MutablePolygonalSurface
        sage: components = RootedComponents_MutablePolygonalSurface(S)

    """

    def __init__(self, surface):
        self._surface = surface

    def __getitem__(self, root):
        r"""
        Return the labels of the connected component rooted at the label
        ``root``.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S.add_polygon(polygons.square())
            1
            sage: S.glue((0, 0), (1, 0))

            sage: from flatsurf.geometry.surface import RootedComponents_MutablePolygonalSurface
            sage: components = RootedComponents_MutablePolygonalSurface(S)
            sage: components[0]
            (0, 1)

        """
        return self._surface.component(root)

    def __iter__(self):
        r"""
        Iterate over the keys of this mapping, i.e., the root labels of the
        connected components.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S.add_polygon(polygons.square())
            1
            sage: S.glue((0, 0), (1, 0))

            sage: from flatsurf.geometry.surface import RootedComponents_MutablePolygonalSurface
            sage: components = RootedComponents_MutablePolygonalSurface(S)
            sage: list(components)
            [0]

        """
        # Shortcut enumeration if this is known to be a connected surface.
        connected = "Connected" in self._surface.category().axioms()

        for root in self._surface._roots:
            yield root
            if connected:
                return

        labels = set(self._surface._polygons)
        for root in self._surface._roots:
            for label in self._surface.component(root):
                labels.remove(label)

        while labels:
            root = LabeledView(surface=self._surface, view=labels, finite=True).min()

            yield root
            if connected:
                return
            for label in self._surface.component(root):
                labels.remove(label)

    def __len__(self):
        r"""
        Return the number of connected components of this surface.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)

            sage: from flatsurf import polygons
            sage: S.add_polygon(polygons.square())
            0
            sage: S.add_polygon(polygons.square())
            1
            sage: S.glue((0, 0), (1, 0))

            sage: from flatsurf.geometry.surface import RootedComponents_MutablePolygonalSurface
            sage: components = RootedComponents_MutablePolygonalSurface(S)
            sage: len(components)
            1

        """
        components = 0
        for root in self:
            components += 1
        return components


class LabeledCollection:
    r"""
    Abstract base class for collection of labels as returned by ``labels()``
    methods of surfaces.

    This also serves as a base clas for things such as ``polygons()`` that are
    tied to labels.

    INPUT:

    - ``surface`` -- a polygonal surface, the labels are taken from that
      surface; subclasses might change this to only represent a subset of the
      labels of this surface

    - ``finite`` -- a boolean or ``None`` (default: ``None``); whether this is
      a finite set; if ``None``, it is not known whether the set is finite
      (some operations might not be supported in that case or not terminate if
      the set is actually infinite.)

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: labels = S.labels()

        sage: from flatsurf.geometry.surface import LabeledCollection
        sage: isinstance(labels, LabeledCollection)
        True

    """

    def __init__(self, surface, finite=None):
        if finite is None and surface.is_finite_type():
            finite = True

        self._surface = surface
        self._finite = finite

    def __repr__(self):
        r"""
        Return a printable representation of this set.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.labels()
            (0,)

            sage: S = translation_surfaces.infinite_staircase()
            sage: S.labels()
            (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

        """
        from itertools import islice

        items = list(islice(self, 17))

        if self._finite is True or len(items) < 17:
            return repr(tuple(self))

        return f"({', '.join(str(x) for x in islice(self, 16))}, …)"

    def __len__(self):
        r"""
        Return the size of this set.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: len(S.labels())
            1

        Python does not allow ``__len__`` to return anything but an integer, so
        we cannot return infinity::

            sage: S = translation_surfaces.infinite_staircase()
            sage: len(S.labels())
            Traceback (most recent call last):
            ...
            NotImplementedError: len() of an infinite set

        """
        if self._finite is False:
            raise TypeError("infinite set has no integer length")

        length = 0
        for x in self:  # pylint: disable=not-an-iterable
            length += 1

        return length

    def __contains__(self, x):
        r"""
        Return whether ``x`` is contained in this set.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: labels = S.labels()
            sage: 0 in labels
            True
            sage: 1 in labels
            False

        """
        for item in self:  # pylint: disable=not-an-iterable
            if x == item:
                return True

        return False


class LabeledView(LabeledCollection):
    r"""
    A set of labels (or something resembling labels such as ``polygons()``)
    backed by another collection ``view``.

    INPUT:

    - ``surface`` -- a polygonal surface, the labels in ``view`` are labels of
      that surface

    - ``view`` -- a collection that all queries are going to be redirected to.
      Note that ``labels()`` guarantees that iteration over labels happens in a
      breadth-first-search so iteration over ``view`` must follow that same
      order. However, subclasses can remove this requirement by overriding
      :meth:`__iter__`.

    - ``finite`` -- a boolean or ``None`` (default: ``None``); whether this is
      a finite set; if ``None``, it is not known whether the set is finite
      (some operations might not be supported in that case or not terminate if
      the set is actually infinite.)

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.t_fractal()
        sage: labels = S.labels()

        sage: from flatsurf.geometry.surface import LabeledView
        sage: isinstance(labels, LabeledView)
        True

    """

    def __init__(self, surface, view, finite=None):
        super().__init__(surface, finite=finite)
        self._view = view

    def __iter__(self):
        return iter(self._view)

    def __contains__(self, x):
        return x in self._view

    def __len__(self):
        return len(self._view)

    def min(self):
        r"""
        Return a minimal item in this set.

        If the items can be compared, this is just the actual ``min`` of the
        items.

        Otherwise, we take the one with minimal ``repr``.

        .. NOTE::

            If the items cannot be compared, and there are clashes in the
            ``repr``, this method will fail.

            Also, if comparison of items is not consistent, then this can
            produce somewhat random output.

            Finally, note with this approach the min of a set is not the always
            the min of the mins of a each partition of that set.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.t_fractal()
            sage: S.labels().min()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot determine minimum of an infinite set

        ::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.cathedral(1,2)
            sage: C.labels().min()
            0

        ::

            sage: labels = list(C.labels())[:3]
            sage: from flatsurf.geometry.surface import LabeledView
            sage: LabeledView(C, labels).min()
            0

        """
        if self._finite is False:
            raise NotImplementedError("cannot determine minimum of an infinite set")

        try:
            return min(self)
        except TypeError:
            reprs = {repr(item): item for item in self}
            if len(reprs) != len(self):
                raise TypeError(
                    "cannot determine minimum of tset without ordering and with non-unique repr()"
                )
            return reprs[min(reprs)]


class ComponentLabels(LabeledCollection):
    r"""
    The labels of a connected component.

    INPUT:

    - ``surface`` -- a polygonal surface

    - ``root`` -- a label of the connected component from which enumeration of
      the component starts.

    - ``finite`` -- a boolean or ``None`` (default: ``None``); whether this is
      a finite component; if ``None``, it is not known whether the component is
      finite (some operations might not be supported in that case or not
      terminate if the component is actually infinite.)

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.t_fractal()
        sage: component = S.component(0)

        sage: from flatsurf.geometry.surface import ComponentLabels
        sage: isinstance(component, ComponentLabels)
        True

    """

    def __init__(self, surface, root, finite=None):
        super().__init__(surface, finite=finite)
        self._root = root

    def __iter__(self):
        r"""
        Return an iterator of this component that enumerates labels starting
        from the root label in a breadth-first-search.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.cathedral(1, 2)
            sage: component = C.component(0)
            sage: list(component)
            [0, 1, 3, 2]

        """
        from collections import deque

        seen = set()
        pending = deque([self._root])

        while pending:
            label = pending.popleft()
            if label in seen:
                continue

            seen.add(label)

            yield label
            for e in range(len(self._surface.polygon(label).vertices())):
                cross = self._surface.opposite_edge(label, e)
                if cross is not None:
                    pending.append(cross[0])


class Labels(LabeledCollection, collections.abc.Set):
    r"""
    The labels of a surface.

    .. NOTE::

        This is a generic implementation that represents the set of labels of a
        surface in a breadth-first iteration starting from the root labels of the
        connected components.

        This implementation makes no assumption on the surface and can be very
        slow to answer, e.g., containment or compute the number of labels in
        the surface (because it needs to iterate over the entire surface.)

        When possible, a faster implementation should be used such as
        :class:`LabelsFromView`.

    EXAMPLES::

        sage: from flatsurf import polygons, similarity_surfaces
        sage: T = polygons.triangle(1, 2, 5)
        sage: S = similarity_surfaces.billiard(T)
        sage: S = S.minimal_cover("translation")

        sage: labels = S.labels()
        sage: labels
        ((0, 1, 0), (1, 1, 0), (1, 0, -1), (1, 1/2*c0, 1/2*c0), (0, 1/2*c0, -1/2*c0), (0, 0, 1), (0, -1/2*c0, -1/2*c0), (0, 0, -1), (0, -1/2*c0, 1/2*c0), (0, 1/2*c0, 1/2*c0),
         (1, 1/2*c0, -1/2*c0), (1, -1/2*c0, -1/2*c0), (1, 0, 1), (1, -1/2*c0, 1/2*c0), (1, -1, 0), (0, -1, 0))

    TESTS::

        sage: from flatsurf.geometry.surface import Labels
        sage: type(labels) == Labels
        True

    """

    def __iter__(self):
        for component in self._surface.components():
            yield from component


class LabelsFromView(Labels, LabeledView):
    r"""
    The labels of a surface backed by another set that can quickly compute the
    length of the labels and decide containment in the set.

    .. NOTE::

        Iteration of the view collection does not have to be in breadth-first
        search order in the surface since this class is picking up the generic
        :meth:`Labels.__iter__`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: C = translation_surfaces.cathedral(1, 2)
        sage: labels = C.labels()

        sage: from flatsurf.geometry.surface import LabelsFromView
        sage: type(labels) == LabelsFromView
        True

    """


class Polygons(LabeledCollection, collections.abc.Collection):
    r"""
    The collection of polygons of a surface.

    The polygons are returned in the same order as labels of the surface are
    returned by :class:`.Labels`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: C = translation_surfaces.cathedral(1, 2)
        sage: polygons = C.polygons()

        sage: from flatsurf.geometry.surface import Polygons
        sage: isinstance(polygons, Polygons)
        True

    """

    def __iter__(self):
        r"""
        Iterate over the polygons in the same order as ``labels()`` does.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.cathedral(1, 2)
            sage: labels = C.labels()
            sage: polygons = C.polygons()

            sage: for entry in zip(labels, polygons):
            ....:     print(entry)
            (0, Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            (1, Polygon(vertices=[(1, 0), (1, -2), (3/2, -5/2), (2, -2), (2, 0), (2, 1), (2, 3), (3/2, 7/2), (1, 3), (1, 1)]))
            (3, Polygon(vertices=[(3, 0), (7/2, -1/2), (11/2, -1/2), (6, 0), (6, 1), (11/2, 3/2), (7/2, 3/2), (3, 1)]))
            (2, Polygon(vertices=[(2, 0), (3, 0), (3, 1), (2, 1)]))

        """
        for label in self._surface.labels():
            yield self._surface.polygon(label)

    def __len__(self):
        r"""
        Return the number of polygons in this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.cathedral(1, 2)
            sage: polygons = C.polygons()
            sage: len(polygons)
            4

        """
        return len(self._surface.labels())


class Polygons_MutableOrientedSimilaritySurface(Polygons):
    r"""
    The collection of polygons of a :class:`MutableOrientedSimilaritySurface`.

    This is a faster version of :class:`Polygons`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: C = translation_surfaces.cathedral(1, 2)
        sage: polygons = C.polygons()

        sage: from flatsurf.geometry.surface import Polygons_MutableOrientedSimilaritySurface
        sage: isinstance(polygons, Polygons_MutableOrientedSimilaritySurface)
        True

    """

    def __init__(self, surface):
        # This hack makes __len__ 20% faster (it saves one attribute lookup.)
        self._polygons = surface._polygons
        super().__init__(surface)

    def __len__(self):
        return len(self._polygons)


class Edges(LabeledCollection, collections.abc.Set):
    r"""
    The set of edges of a surface.

    The set of edges contains of pairs (label, index) with the labels of the
    polygons and the actual edges indexed from 0 in the second component.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: C = translation_surfaces.cathedral(1, 2)
        sage: edges = C.edges()
        sage: edges
        ((0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (2, 0), (2, 1), (2, 2), (2, 3))

    TESTS::

        sage: from flatsurf.geometry.surface import Edges
        sage: isinstance(edges, Edges)
        True

    """

    def __iter__(self):
        for label, polygon in zip(self._surface.labels(), self._surface.polygons()):
            for edge in range(len(polygon.vertices())):
                yield (label, edge)

    def __contains__(self, x):
        label, edge = x
        if label not in self._surface.labels():
            return False

        polygon = self._surface.polygon(label)
        return 0 <= len(polygon.vertices()) < edge


class Gluings(LabeledCollection, collections.abc.Set):
    r"""
    The set of gluings of the surface.

    Each gluing consists of two pairs (label, index) that describe the edges
    being glued.

    Note that each gluing (that is not a self-gluing) is reported twice.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: gluings = S.gluings()
        sage: gluings
        (((0, 0), (0, 2)), ((0, 1), (0, 3)), ((0, 2), (0, 0)), ((0, 3), (0, 1)))

    TESTS::

        sage: from flatsurf.geometry.surface import Gluings
        sage: isinstance(gluings, Gluings)
        True

    """

    def __iter__(self):
        for label, edge in self._surface.edges():
            cross = self._surface.opposite_edge(label, edge)
            if cross is None:
                continue
            yield (label, edge), cross

    def __contains__(self, x):
        x, y = x

        if x not in self._surface.edges():
            return False

        cross = self._surface.opposite_edge(*x)
        if cross is None:
            return False

        return y == cross


# Import deprecated symbols so imports using flatsurf.geometry.surface do not break.
from flatsurf.geometry.surface_legacy import (  # noqa, we import at the bottom of the file to break a circular import  # pylint: disable=wrong-import-position
    Surface,
    Surface_list,
    Surface_dict,
    surface_list_from_polygons_and_gluings,
)
