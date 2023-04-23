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
    # TODO: Lets us add polygons.
    def __init__(self, base, category=None):
        from sage.all import ZZ
        self._next_label = ZZ(0)
        self._base_label = None

        self._polygons = {}

        self._mutable = True

        super().__init__(base, category=category)

    def add_polygon(self, polygon, label=None):
        if not self._mutable:
            raise Exception

        if label is None:
            while self._next_label in self._polygons:
                self._next_label += 1
            label = self._next_label

        if label in self._polygons:
            raise ValueError  # must remove first

        if self._base_label is None:
            self._base_label = label

        self._polygons[label] = polygon

    def base_label(self):
        if self._base_label is None:
            raise NotImplementedError

        return self._base_label

    def polygon(self, label):
        return self._polygons[label]

    def set_immutable(self):
        if self._mutable:
            self._mutable = False

        self._refine_category_(self.refined_category())

    def is_finite(self):
        return True

    def is_mutable(self):
        return self._mutable

    def __eq__(self, other):
        if not isinstance(other, MutablePolygonalSurface):
            return False

        if self._polygons != other._polygons:
            return False

        if self._base_label != other._base_label:
            return False

        if self._mutable != other._mutable:
            return False

        return True

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        if self._mutable:
            raise Exception
        # TODO
        return 0

    def _repr_(self):
        if not self.is_finite():
            return "Surface built from infinitely many polygons"
        if self.num_polygons() == 1:
            return "Surface built from 1 polygon"

        return "Surface built from {} polygons".format(self.num_polygons())


class OrientedSimilaritySurface(Surface_base):
    Element = SurfacePoint

    def __init__(self, base, category=None):
        from flatsurf.geometry.categories import SimilaritySurfaces
        if category is None:
            category = SimilaritySurfaces().Oriented()

        category &= SimilaritySurfaces().Oriented()

        super().__init__(base, category=category)

    def _an_element_(self):
        r"""
        Return a point of this surface.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: from flatsurf.geometry.surface import Surface_list

            sage: S = Surface_list(QQ)
            sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
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

        if self.num_polygons() == 0:
            return other.num_polygons() == 0
        if other.num_polygons() == 0:
            return False

        if self.num_polygons() != other.num_polygons():
            return False

        if self.base_label() != other.base_label():
            return False
        if self.polygon(self.base_label()) != other.polygon(self.base_label()):
            return False

        if not self.is_finite():
            raise NotImplementedError("cannot compare these infinite surfaces yet")

        for label, polygon in self.label_polygon_iterator():
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


class MutableOrientedSimilaritySurface(OrientedSimilaritySurface, MutablePolygonalSurface):
    def __init__(self, base, category=None):
        self._gluings = {}

        super().__init__(base, category=category)

    def glue(self, x, y):
        if not self._mutable:
            raise Exception

        labels = x[0], y[0]
        edges = x[1], y[1]

        for label in labels:
            if label not in self._polygons:
                raise ValueError

            if label not in self._gluings:
                self._gluings[label] = [None] * self.polygon(label).num_edges()

        for label, edge in zip(labels, edges):
            if self._gluings[label][edge] is not None:
                cross_label, cross_edge = self._gluings[label][edge]
                self._gluings[cross_label][cross_edge] = None
                self._gluings[label][edge] = None

        self._gluings[labels[0]][edges[0]] = (labels[1], edges[1])
        self._gluings[labels[1]][edges[1]] = (labels[0], edges[0])

    def opposite_edge(self, label, edge):
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
            raise Exception
        # TODO
        return 0


class Labels(collections.abc.Set):
    def __init__(self, surface, labels=None, len=None):
        self._surface = surface
        self._labels = labels
        self._len = len

    def __len__(self):
        if self._len is None:
            if not self._surface.is_finite():
                from sage.all import infinity
                self._len = infinity
            elif self._labels is not None:
                self._len = len(self._labels)

            self._len = 0
            for label in self:
                self._len += 1

        return self._len

    def __contains__(self, x):
        if self._labels is not None:
            return x in self._labels

        for label in self:
            if x == label:
                return True

        return False

    def __iter__(self):
        from collections import deque

        seen = set()
        pending = deque()

        pending.append(self._surface.base_label())

        while pending:
            label = pending.popleft()
            if label in seen:
                continue

            seen.add(label)

            yield label
            for e in range(self._surface.polygon(label).num_edges()):
                pending.append(self._surface.opposite_edge(label, e)[0])

    def __repr__(self):
        if self._surface.is_finite():
            return repr(tuple(self))

        from itertools import islice
        return f"({', '.join(str(x) for x in islice(self, 16))}, …)"

# Import deprecated symbols so imports using flatsurf.geometry.surface do not break.
from flatsurf.geometry.surface_legacy import Surface, Surface_list, Surface_dict, surface_list_from_polygons_and_gluings, LabelComparator, BaseRingChangedSurface
