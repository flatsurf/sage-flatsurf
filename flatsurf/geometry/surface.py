# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
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
    def _refine_category_(self, category):
        r"""
        Refine the category of this surface to a subcategory ``category``.

        We need to override this method from Parent since we need to disable a hashing check that is otherwise enabled when doctesting.

        Since our surfaces are not hashable (since equality of infinite
        surfaces is a delicate matter,) that hash check otherwise fails with a
        NotImplementedError.

        """
        from sage.structure.debug_options import debug

        old_refine_category_hash_check = debug.refine_category_hash_check
        debug.refine_category_hash_check = False
        try:
            super()._refine_category_(category)
        finally:
            debug.refine_category_hash_check = old_refine_category_hash_check


class MutablePolygonalSurface(Surface_base):
    def __init__(self, base, category=None):
        from sage.all import ZZ

        self._next_label = ZZ(0)
        self._roots = ()

        self._polygons = {}

        self._mutable = True

        super().__init__(base, category=category)

    def _test_refined_category(self, **options):
        if self._mutable:
            return

        super()._test_refined_category(**options)

    def add_polygon(self, polygon, *, label=None):
        if not self._mutable:
            raise Exception

        if label is None:
            while self._next_label in self._polygons:
                self._next_label += 1
            label = self._next_label

        if label in self._polygons:
            raise ValueError  # must remove first

        self._polygons[label] = polygon
        return label

    def add_polygons(self, polygons):
        import warnings

        warnings.warn(
            "add_polygons() has been deprecated and will be removed in a future version of sage-flatsurf; use labels = [add_polygon(p) for p in polygons] instead"
        )

        return [self.add_polygon(p) for p in polygons]

    def set_default_graphical_surface(self, graphical_surface):
        raise NotImplementedError(
            "set_default_graphical_surface() has been removed from this version of sage-flatsurf. If you want to change the default plotting of a surface create a subclass and override graphical_surface() instead"
        )

    def remove_polygon(self, label):
        if not self._mutable:
            raise Exception

        self._polygons.pop(label)
        self._roots = tuple(root for root in self._roots if root != label)

    def _components(self):
        return RootedComponents_MutablePolygonalSurface(self)

    def roots(self):
        return LabeledView(self, self._components().keys(), finite=True)

    def components(self):
        return LabeledView(self, self._components().values(), finite=True)

    def polygon(self, label):
        return self._polygons[label]

    def set_immutable(self):
        if self._mutable:
            self.set_roots(self.roots())
            self._mutable = False

        self._refine_category_(self.refined_category())

    def is_finite_type(self):
        return True

    def is_mutable(self):
        return self._mutable

    def __eq__(self, other):
        if not isinstance(other, MutablePolygonalSurface):
            return False

        if self._polygons != other._polygons:
            return False

        # Note that the order of the root labels matters since it changes the order of iteration in labels()
        if self._roots != other._roots:
            return False

        if self._mutable != other._mutable:
            return False

        return True

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        if self._mutable:
            raise TypeError("cannot hash a mutable surface")

        return hash((tuple(self.labels()), tuple(self.polygons()), self._roots))

    def _repr_(self):
        if not self.is_finite_type():
            return "Surface built from infinitely many polygons"

        if len(self.labels()) == 0:
            return "Empty Surface"

        return f"{self._describe_surface()} built from {self._describe_polygons()}"

    def _describe_surface(self):
        return "Surface"

    def _describe_polygons(self):
        polygons = [
            (-p.erase_marked_vertices().num_edges(), p.describe_polygon())
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
        if not self._mutable:
            raise Exception

        root = [label for label in self.labels() if label == root]
        if not root:
            raise ValueError("root must be a label in the surface")
        assert len(root) == 1
        root = root[0]

        component = [component for component in self.components() if root in component]

        self._roots = tuple(r for r in self._roots if r not in component) + (root,)

    def set_roots(self, roots):
        if not self._mutable:
            raise Exception

        roots = [[label for label in self.labels() if label == root] for root in roots]

        if any(len(root) == 0 for root in roots):
            raise ValueError("roots must be existing labels in the surface")

        assert all(len(root) == 1 for root in roots)

        roots = tuple(root[0] for root in roots)

        for component in self.components():
            if len([root for root in roots if root in component]) > 1:
                raise ValueError("there must be at most one root for each connected component")

        self._roots = tuple(roots)

    def change_base_label(self, label):
        import warnings

        warnings.warn(
            "change_base_label() has been deprecated and will be removed in a future version of sage-flatsurf; use set_root() instead"
        )

        self.set_root(label)

    @cached_method
    def labels(self):
        return LabelsView(self, self._polygons.keys(), finite=True)

    @cached_method
    def polygons(self):
        return Polygons_MutableOrientedSimilaritySurface(self)


class OrientedSimilaritySurface(Surface_base):
    Element = SurfacePoint

    def __init__(self, base, category=None):
        from sage.categories.all import Rings

        if base not in Rings():
            raise TypeError

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

    def _an_element_(self):
        r"""
        Return a point of this surface.

        EXAMPLES::

            sage: from flatsurf import polygons, MutableOrientedSimilaritySurface

            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))

            sage: S.an_element()
            Point (1/2, 1/2) of polygon 0

        """
        label = next(iter(self.labels()))
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
            for edge in range(polygon.num_edges()):
                if self.opposite_edge(label, edge) != other.opposite_edge(label, edge):
                    return False

        return True


class MutableOrientedSimilaritySurface(
    OrientedSimilaritySurface, MutablePolygonalSurface
):
    def __init__(self, base, category=None):
        self._gluings = {}

        from flatsurf.geometry.categories import SimilaritySurfaces

        if category is None:
            category = SimilaritySurfaces().Oriented().FiniteType()

        category &= SimilaritySurfaces().Oriented().FiniteType()

        super().__init__(base, category=category)

    @classmethod
    def from_surface(cls, surface, category=None):
        if not surface.is_finite_type():
            raise TypeError
        self = MutableOrientedSimilaritySurface(
            surface.base_ring(), category=category or surface.category()
        )
        for label in surface.labels():
            self.add_polygon(surface.polygon(label), label=label)

        for label in surface.labels():
            for edge in range(surface.polygon(label).num_edges()):
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
        label = super().add_polygon(polygon, label=label)
        assert label not in self._gluings
        self._gluings[label] = [None] * polygon.num_edges()
        return label

    def remove_polygon(self, label):
        self._unglue_polygon(label)
        self._gluings.pop(label)

        super().remove_polygon(label)

    def unglue(self, label, edge):
        if not self._mutable:
            raise Exception

        cross = self._gluings[label][edge]
        if cross is not None:
            self._gluings[cross[0]][cross[1]] = None
        self._gluings[label][edge] = None

    def _unglue_polygon(self, label):
        for edge, cross in enumerate(self._gluings[label]):
            if cross is None:
                continue
            cross_label, cross_edge = cross
            self._gluings[cross_label][cross_edge] = None
        self._gluings[label] = [None] * self.polygon(label).num_edges()

    def glue(self, x, y):
        # TODO: Update roots.
        if not self._mutable:
            raise Exception

        if x[0] not in self._polygons:
            raise ValueError

        if y[0] not in self._polygons:
            raise ValueError

        self.unglue(*x)
        self.unglue(*y)

        self._gluings[x[0]][x[1]] = y
        self._gluings[y[0]][y[1]] = x

    def set_edge_pairing(self, label0, edge0, label1, edge1):
        r"""
        TESTS::

            sage: from flatsurf import polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.set_edge_pairing(0, 0, 0, 2)
            doctest:warning
            ...
            UserWarning: set_edge_pairing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead
            sage: S.set_edge_pairing(0, 1, 0, 3)

        """
        import warnings

        warnings.warn(
            "set_edge_pairing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead"
        )
        return self.glue((label0, edge0), (label1, edge1))

    def change_edge_gluing(self, label0, edge0, label1, edge1):
        r"""
        TESTS::

            sage: from flatsurf import polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))

        """
        import warnings

        warnings.warn(
            "change_edge_gluing(label0, edge0, label1, edge1) has been deprecated and will be removed in a future version of sage-flatsurf; use glue((label0, edge0), (label1, edge1)) instead"
        )
        return self.glue((label0, edge0), (label1, edge1))

    def change_polygon_gluings(self, label, gluings):
        r"""
        TESTS::

            sage: from flatsurf import polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))

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
        import warnings

        warnings.warn(
            "change_polygon() has been deprecated and will be removed in a future version of sage-flatsurf; use replace_polygon() or remove_polygon() and add_polygon() instead"
        )

        if not self._mutable:
            raise Exception

        # Note that this obscure feature. If the number of edges is unchanged, we keep the gluings, otherwise we trash them all.
        if polygon.num_edges() != self.polygon(label).num_edges():
            self._unglue_polygon(label)
            self._gluings[label] = [None] * polygon.num_edges()

        self._polygons[label] = polygon

        if gluing_list is not None:
            for i, cross in enumerate(gluing_list):
                self.glue((label, i), cross)

    def replace_polygon(self, label, polygon):
        r"""
        Replace the polygon ``label`` with ``polygon`` while keeping its
        gluings intact.

        INPUT:

        - ``label`` -- an element of :meth:`labels`

        - ``polygon`` -- a Euclidean polygon

        EXAMPLES::

            sage: from flatsurf import polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))

            sage: S.replace_polygon(0, polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)]))

        The replacement of a polygon must have the same number of sides::

            sage: S.replace_polygon(0, polygon(vertices=[(0, 0), (2, 0), (2, 2)]))
            Traceback (most recent call last):
            ...
            ValueError: polygon must be a quadrilateral

        To replace the polygon without keeping its glueings, remove the polygon
        first and then add a new one::

            sage: S.remove_polygon(0)
            sage: S.add_polygon(polygon(vertices=[(0, 0), (2, 0), (2, 2)]), label=0)
            0

        """
        old = self.polygon(label)

        if old.num_edges() != polygon.num_edges():
            from flatsurf.geometry.polygon import Polygon

            article, singular, plural = Polygon._describe_polygon(old.num_edges())
            raise ValueError(f"polygon must be {article} {singular}")

        self._polygons[label] = polygon

    def opposite_edge(self, label, edge=None):
        if edge is None:
            label, edge = label
        return self._gluings[label][edge]

    def __eq__(self, other):
        if not isinstance(other, MutableOrientedSimilaritySurface):
            return False

        if not super().__eq__(other):
            return False

        if self._gluings != other._gluings:
            return False

        return True

    def __hash__(self):
        if self._mutable:
            raise TypeError("cannot hash a mutable surface")

        return hash((super().__hash__(), tuple(self._gluings)))


class BaseRingChangedSurface(OrientedSimilaritySurface):
    def __init__(self, surface, ring, category=None):
        self._reference = surface
        super().__init__(ring, category=category or surface.category())

    def is_mutable(self):
        return False

    def roots(self):
        return self._reference.roots()

    def polygon(self, label):
        return self._reference.polygon(label).change_ring(self.base_ring())

    def opposite_edge(self, label, edge):
        return self._reference.opposite_edge(label, edge)


class RootedComponents_MutablePolygonalSurface(collections.abc.Mapping):
    def __init__(self, surface):
        self._surface = surface

    def __getitem__(self, root):
        return self._surface.component(root)

    def __iter__(self):
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
            try:
                root = min(labels)
            except TypeError:
                if len(set((repr(label) for label in labels))) != len(labels):
                    raise TypeError("cannot determine root label consistently in this surface since the labels cannot be ordered and their repr conversion to string are not unique")
                root = min(labels, key=lambda label: repr(label))

            yield root
            if connected:
                return
            for label in self._surface.component(root):
                labels.remove(label)

    def __len__(self):
        components = 0
        for root in self:
            components += 1
        return components


class LabeledCollection:
    def __init__(self, surface, finite=None):
        if finite is None and surface.is_finite_type():
            finite = True

        self._surface = surface
        self._finite = finite

    def __repr__(self):
        from itertools import islice
        items = list(islice(self, 17))

        if self._finite is True or len(items) < 17:
            return repr(tuple(self))

        return f"({', '.join(str(x) for x in islice(self, 16))}, …)"

    def __len__(self):
        if self._finite is False:
            raise TypeError("infinite set has no integer length")

        length = 0
        for x in self:
            length += 1

        return length

    def __contains__(self, x):
        for label in self:
            if x == label:
                return True

        return False


class LabeledView(LabeledCollection):
    def __init__(self, surface, view, finite=None):
        super().__init__(surface, finite=finite)
        self._view = view

    def __iter__(self):
        return iter(self._view)

    def __contains__(self, x):
        return x in self._view

    def __len__(self):
        return len(self._view)


class ComponentLabels(LabeledCollection):
    def __init__(self, surface, root, finite=None):
        super().__init__(surface, finite=finite)
        self._root = root

    def __iter__(self):
        from collections import deque

        seen = set()
        pending = deque([self._root])

        while pending:
            label = pending.popleft()
            if label in seen:
                continue

            seen.add(label)

            yield label
            for e in range(self._surface.polygon(label).num_edges()):
                cross = self._surface.opposite_edge(label, e)
                if cross is not None:
                    pending.append(cross[0])


class Labels(LabeledCollection, collections.abc.Set):
    def __iter__(self):
        for component in self._surface.components():
            for label in component:
                yield label


class LabelsView(Labels, LabeledView):
    pass


class Polygons(LabeledCollection, collections.abc.Collection):
    def __iter__(self):
        for label in self._surface.labels():
            yield self._surface.polygon(label)

    def __len__(self):
        return len(self._surface.labels())


class Polygons_MutableOrientedSimilaritySurface(Polygons):
    def __init__(self, surface):
        # This hack makes __len__ 20% faster (it saves one attribute lookup.)
        self._polygons = surface._polygons
        super().__init__(surface)

    def __len__(self):
        return len(self._polygons)


class Edges(LabeledCollection, collections.abc.Set):
    def __iter__(self):
        for label, polygon in zip(self._surface.labels(), self._surface.polygons()):
            for edge in range(polygon.num_edges()):
                yield (label, edge)

    def __contains__(self, x):
        label, edge = x
        if label not in self._surface.labels():
            return False

        polygon = self._surface.polygon(label)
        return 0 <= polygon.num_edges() < edge


class Gluings(LabeledCollection, collections.abc.Set):
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
from flatsurf.geometry.surface_legacy import (
    Surface,
    Surface_list,
    Surface_dict,
    surface_list_from_polygons_and_gluings,
    LabelComparator,
)
