r"""
Construction of Chamanara's surfaces which depend on a parameter alpha less than one.
See the paper "Affine automorphism groups of surfaces of infinite type" in which the surface
is called $X_\alpha$.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: s = translation_surfaces.chamanara(1/2)
    sage: s.plot()
    ...Graphics object consisting of 129 graphics primitives

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2013-2019 W. Patrick Hooper
#                      2013-2019 Vincent Delecroix
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

from .surface import Surface
from .half_dilation_surface import HalfDilationSurface
from sage.rings.integer_ring import ZZ


def ChamanaraPolygon(alpha):
    from sage.categories.fields import Fields

    field = alpha.parent()
    if field not in Fields():
        ValueError("The value of alpha must lie in a field.")
    if alpha <= 0 or alpha >= 1:
        ValueError("The value of alpha must be between zero and one.")
    # The value of x is $\sum_{n=0}^\infty \alpha^n$.
    x = 1 / (1 - alpha)
    from .polygon import polygons

    return polygons((1, 0), (-x, x), (0, -1), (x - 1, 1 - x))


class ChamanaraSurface(Surface):
    r"""
    The ChamanaraSurface $X_{\alpha}$.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import ChamanaraSurface
        sage: ChamanaraSurface(1/2)
        Chamanara surface with parameter 1/2
    """

    def __init__(self, alpha):
        self._p = ChamanaraPolygon(alpha)

        field = alpha.parent()
        if not field.is_field():
            field = field.fraction_field()

        self.rename("Chamanara surface with parameter {}".format(alpha))

        super().__init__(field, ZZ(0), finite=False, mutable=False)

    def polygon(self, lab):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C.polygon('a')
            Traceback (most recent call last):
            ...
            ValueError: invalid label 'a'
        """
        if lab not in ZZ:
            raise ValueError("invalid label {!r}".format(lab))
        return self._p

    def opposite_edge(self, p, e):
        if e == 0 or e == 2:
            return 1 - p, e
        elif e == 1:
            if p < 0:
                return p + 1, 3
            elif p > 1:
                return p - 1, 3
            else:
                # p==0 or p==1
                return 1 - p, 1
        else:
            # e==3
            if p <= 0:
                return p - 1, 1
            else:
                # p>=1
                return p + 1, 1

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C == C
            True
            sage: D = translation_surfaces.chamanara(1/3)
            sage: C == D
            False

        """
        if isinstance(other, ChamanaraSurface):
            return (
                self._p == other._p
                and self._base_ring == other._base_ring
                and self._base_label == other._base_label
            )

        return super().__eq__(other)


def chamanara_half_dilation_surface(alpha, n=8):
    r"""
    Return Chamanara's surface thought of as a Half Dilation surface.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import chamanara_half_dilation_surface
        sage: s = chamanara_half_dilation_surface(1/2)
        sage: TestSuite(s).run()
    """
    s = HalfDilationSurface(ChamanaraSurface(alpha))
    adjacencies = [(0, 1)]
    for i in range(n):
        adjacencies.append((-i, 3))
        adjacencies.append((i + 1, 3))
    s.graphical_surface(adjacencies=adjacencies)
    return s


def chamanara_surface(alpha, n=8):
    r"""
    Return Chamanara's surface thought of as a translation surface.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: s = chamanara_surface(1/2)
        sage: TestSuite(s).run()
    """
    s = chamanara_half_dilation_surface(alpha).minimal_cover(cover_type="translation")
    label = s.base_label()
    adjacencies = [(label, 1)]
    for i in range(n):
        adjacencies.append((label, 3))
        label = s.opposite_edge(label, 3)[0]
    label = s.base_label()
    label = s.opposite_edge(label, 1)[0]
    for i in range(n):
        adjacencies.append((label, 3))
        label = s.opposite_edge(label, 3)[0]
    s.graphical_surface(adjacencies=adjacencies)
    return s
