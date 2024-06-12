# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016 W. Patrick Hooper
#                      2022 Vincent Delecroix
#                      2023 Julian Rüth
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
from sage.rings.rational_field import QQ
from flatsurf.geometry.surface import OrientedSimilaritySurface


class AbstractOrigami(OrientedSimilaritySurface):
    r"""
    Abstract base class for (connected) origamis.
    """

    def __init__(self, domain, root=None, base_label=None, category=None):
        self._domain = domain

        if base_label is not None:
            import warnings

            warnings.warn(
                "base_label has been deprecated as a keyword argument for AbstractOrigami and will be removed in a future version of sage-flatsurf; use root instead"
            )
            root = base_label
            base_label = None

        if root is None:
            root = domain.an_element()
        self._root = root

        from flatsurf.geometry.categories import TranslationSurfaces

        if category is None:
            category = TranslationSurfaces()

        category &= TranslationSurfaces().WithoutBoundary().Connected()

        finite = domain.is_finite()
        if finite:
            category &= category.FiniteType()
        else:
            category &= category.InfiniteType()

        from flatsurf.geometry.polygon import polygons

        self._square = polygons.square()

        super().__init__(QQ, category=category)

    def roots(self):
        return (self._root,)

    def labels(self):
        from flatsurf.geometry.surface import LabelsFromView

        return LabelsFromView(self, self._domain)

    def is_mutable(self):
        return False

    def up(self, label):
        raise NotImplementedError

    def down(self, label):
        raise NotImplementedError

    def right(self, label):
        raise NotImplementedError

    def left(self, label):
        raise NotImplementedError

    def _repr_(self):
        return "Some AbstractOrigami"

    def polygon_labels(self):
        return self._domain

    def polygon(self, lab):
        if lab not in self._domain:
            # Updated to print a possibly useful error message
            raise ValueError("Label " + str(lab) + " is not in the domain")

        return self._square

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e == 0:
            return self.down(p), 2
        if e == 1:
            return self.right(p), 3
        if e == 2:
            return self.up(p), 0
        if e == 3:
            return self.left(p), 1
        raise ValueError


class Origami(AbstractOrigami):
    def __init__(
        self,
        r,
        u,
        rr=None,
        uu=None,
        domain=None,
        root=None,
        base_label=None,
        category=None,
    ):
        if domain is None:
            domain = r.parent().domain()

        self._r = r
        self._u = u
        if rr is None:
            rr = ~r
        else:
            for a in domain.some_elements():
                if r(rr(a)) != a:
                    raise ValueError(f"r o rr is not identity on {a}")
                if rr(r(a)) != a:
                    raise ValueError(f"rr o r is not identity on {a}")
        if uu is None:
            uu = ~u
        else:
            for a in domain.some_elements():
                if u(uu(a)) != a:
                    raise ValueError(f"u o uu is not identity on {a}")
                if uu(u(a)) != a:
                    raise ValueError(f"uu o u is not identity on {a}")

        self._perms = [uu, r, u, rr]  # down,right,up,left
        AbstractOrigami.__init__(
            self, domain, root=root, base_label=base_label, category=category
        )

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError(
                "Polygon label p=" + str(p) + " is not in domain=" + str(self._domain)
            )
        if e < 0 or e > 3:
            raise ValueError("Edge value e=" + str(e) + " does not satisfy 0<=e<4.")
        return self._perms[e](p), (e + 2) % 4

    def up(self, label):
        return self.opposite_edge(label, 2)[0]

    def down(self, label):
        return self.opposite_edge(label, 0)[0]

    def right(self, label):
        return self.opposite_edge(label, 1)[0]

    def left(self, label):
        return self.opposite_edge(label, 3)[0]

    def _repr_(self):
        return "Origami defined by r={} and u={}".format(self._r, self._u)

    def __eq__(self, other):
        if not isinstance(other, Origami):
            return False

        return (
            self._perms == other._perms
            and self._domain is other._domain
            and self.roots() == other.roots()
        )

    def __hash__(self):
        return hash((Origami, tuple(self._perms), self._domain, self.roots()))
