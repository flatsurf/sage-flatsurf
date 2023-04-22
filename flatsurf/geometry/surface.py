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


class Surface_base:
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
    pass


class OrientedSimilaritySurface(Surface_base, Parent):
    Element = SurfacePoint

    def __init__(self, base, category=None):
        Surface_base.__init__(self)

        from flatsurf.geometry.categories import SimilaritySurfaces
        if category is None:
            category = SimilaritySurfaces().Oriented()

        category &= SimilaritySurfaces().Oriented()

        Parent.__init__(self, base, category=category)

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


class MutableOrientedSimilaritySurface(OrientedSimilaritySurface):
    # TODO: Lets us add polygons and change edge gluings.
    pass


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


# Import deprecated symbols so imports using flatsurf.geometry.surface do not break.
from flatsurf.geometry.surface_legacy import Surface, Surface_list, Surface_dict, surface_list_from_polygons_and_gluings, LabelComparator, BaseRingChangedSurface
