r"""
Chamanara's surfaces which depend on a parameter `\alpha` less than one.

REFERENCES:

- Defined as `X_\alpha` in Chamanara, Reza, "Affine automorphism groups of
  surfaces of infinite type", City University of New York, 2002.

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

EXAMPLES:

.. jupyter-execute::

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

from flatsurf.geometry.surface import OrientedSimilaritySurface
from flatsurf.geometry.minimal_cover import MinimalTranslationCover
from flatsurf.geometry.lazy import LazyRelabeledSurface
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
    from flatsurf import Polygon

    return Polygon(edges=[(1, 0), (-x, x), (0, -1), (x - 1, 1 - x)])


class ChamanaraSurface(OrientedSimilaritySurface):
    r"""
    The Chamanara surface $X_{\alpha}$.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import ChamanaraSurface
        sage: S = ChamanaraSurface(1/2); S
        Chamanara surface with parameter 1/2

    TESTS::

        sage: TestSuite(S).run()
    """

    def __init__(self, alpha):
        self._p = ChamanaraPolygon(alpha)

        field = alpha.parent()
        if not field.is_field():
            field = field.fraction_field()

        self.rename("Chamanara surface with parameter {}".format(alpha))

        from flatsurf.geometry.categories import DilationSurfaces

        super().__init__(
            field,
            category=DilationSurfaces()
            .Oriented()
            .InfiniteType()
            .Compact()
            .WithoutBoundary()
            .Connected()
            .Rational(),
        )

    def is_dilation_surface(self, positive=False):
        r"""
        Return whether this surface is a dilation surface, overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_dilation_surface`.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.is_dilation_surface(positive=True)
            False
            sage: S.is_dilation_surface(positive=False)
            True

        """
        return not positive

    def is_cone_surface(self):
        r"""
        Return whether this surface is a cone surfaces, overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_cone_surface`.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.is_cone_surface()
            False

        """
        return False

    def is_translation_surface(self, positive=True):
        r"""
        Return whether this surfaces is a (half-)translation surface, overrides
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_translation_surface`.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.is_translation_surface(positive=True)
            False
            sage: S.is_translation_surface(positive=False)
            False

        """
        return False

    def labels(self):
        r"""
        Return the labels used to identify the polygons that make up this
        surface, i.e., the integers.

        This overrides
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.labels()
            (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

        """
        from flatsurf.geometry.surface import LabelsFromView

        return LabelsFromView(self, ZZ, finite=False)

    def roots(self):
        r"""
        Return a label in each connected component of this surface.

        This overrides
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.roots()
            (0,)

        """
        return (ZZ(0),)

    def is_mutable(self):
        r"""
        Return whether this surface is mutable which it is not.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: S = ChamanaraSurface(1/2)
            sage: S.is_mutable()
            False

        """
        return False

    def polygon(self, lab):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C.polygon(0)
            Polygon(vertices=[(0, 0), (1, 0), (-1, 2), (-1, 1)])

        """
        if lab not in ZZ:
            raise KeyError(lab)

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

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with equality testing.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import ChamanaraSurface
            sage: hash(ChamanaraSurface(1/2)) == hash(ChamanaraSurface(1/2))
            True

        """
        return hash((self._p, self.base_ring()))

    def graphical_surface(self, **kwds):
        adjacencies = [(0, 1)]
        for i in range(8):
            adjacencies.append((-i, 3))
            adjacencies.append((i + 1, 3))
        return super().graphical_surface(adjacencies=adjacencies, **kwds)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C == C
            True
            sage: D = translation_surfaces.chamanara(1/3)
            sage: C == D
            False

        """
        if not isinstance(other, ChamanaraSurface):
            return False

        return self._p == other._p and self.base_ring() == other.base_ring()


def chamanara_half_dilation_surface(alpha, n=None):
    r"""
    Return Chamanara's surface thought of as a Half Dilation surface.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import chamanara_half_dilation_surface
        sage: s = chamanara_half_dilation_surface(1/2)
        sage: TestSuite(s).run()
    """
    if n is not None:
        import warnings

        warnings.warn(
            "the n keyword of chamanara_half_dilation_surface() is not supported anymore; it will be removed in a future version of sage-flatsurf"
        )

    return ChamanaraSurface(alpha)


class ChamanaraTranslationSurface(LazyRelabeledSurface):
    def __init__(self, alpha):
        super().__init__(MinimalTranslationCover(ChamanaraSurface(alpha)))

    def graphical_surface(self, **kwds):
        label = self.root()
        adjacencies = [(label, 1)]
        for i in range(8):
            adjacencies.append((label, 3))
            label = self.opposite_edge(label, 3)[0]
        label = self.root()
        label = self.opposite_edge(label, 1)[0]
        for i in range(8):
            adjacencies.append((label, 3))
            label = self.opposite_edge(label, 3)[0]
        return super().graphical_surface(adjacencies=adjacencies, **kwds)


def chamanara_surface(alpha, n=None):
    r"""
    Return Chamanara's surface thought of as a translation surface.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: s = chamanara_surface(1/2)
        sage: TestSuite(s).run()
    """
    if n is not None:
        import warnings

        warnings.warn(
            "the n keyword of chamanara_half_dilation_surface() is not supported anymore; it will be removed in a future version of sage-flatsurf"
        )

    return ChamanaraTranslationSurface(alpha)
