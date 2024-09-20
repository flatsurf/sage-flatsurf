# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
#                           2023 Sam Freedman
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
from sage.rings.all import ZZ, QQ, RIF, AA, NumberField, polygen
from sage.modules.all import VectorSpace, vector
from sage.structure.coerce import py_scalar_parent
from sage.structure.element import get_coercion_model, parent
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from flatsurf.geometry.polygon import (
    polygons,
    EuclideanPolygon,
    Polygon,
)

from flatsurf.geometry.surface import (
    OrientedSimilaritySurface,
    MutableOrientedSimilaritySurface,
)
from flatsurf.geometry.origami import Origami


ZZ_1 = ZZ(1)
ZZ_2 = ZZ(2)


def flipper_nf_to_sage(K, name="a"):
    r"""
    Convert a flipper number field into a Sage number field

    .. NOTE::

        Currently, the code is not careful at all with root isolation.

    EXAMPLES::

        sage: import flipper  # optional - flipper  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: import realalg  # optional - flipper
        sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_to_sage
        sage: K = realalg.RealNumberField([-2r] + [0r]*5 + [1r])   # optional - flipper
        sage: K_sage = flipper_nf_to_sage(K)                       # optional - flipper
        sage: K_sage                                               # optional - flipper
        Number Field in a with defining polynomial x^6 - 2 with a = 1.122462048309373?
        sage: AA(K_sage.gen())                                     # optional - flipper
        1.122462048309373?
    """
    r = K.lmbda.interval()
    lower = r.lower * ZZ(10) ** (-r.precision)
    upper = r.upper * ZZ(10) ** (-r.precision)

    p = QQ["x"](K.coefficients)
    s = AA.polynomial_root(p, RIF(lower, upper))
    return NumberField(p, name, embedding=s)


def flipper_nf_element_to_sage(x, K=None):
    r"""
    Convert a flipper number field element into Sage

    EXAMPLES::

        sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_element_to_sage
        sage: import flipper                               # optional - flipper
        sage: T = flipper.load('SB_6')                     # optional - flipper
        sage: h = T.mapping_class('s_0S_1S_2s_3s_4s_3S_5') # optional - flipper
        sage: flipper_nf_element_to_sage(h.dilatation())   # optional - flipper
        a
        sage: AA(_)                                        # optional - flipper
        6.45052513748511?
    """
    if K is None:
        K = flipper_nf_to_sage(x.field)
    coeffs = list(map(QQ, x.coefficients))
    coeffs.extend([0] * (K.degree() - len(coeffs)))
    return K(coeffs)


class EInfinitySurface(OrientedSimilaritySurface):
    r"""
    The surface based on the `E_\infinity` graph.

    The biparite graph is shown below, with edges numbered::

          0   1   2  -2   3  -3   4  -4
        *---o---*---o---*---o---*---o---*...
                |
                |-1
                o

    Here, black vertices are colored ``*``, and white ``o``.
    Black nodes represent vertical cylinders and white nodes
    represent horizontal cylinders.
    """

    def __init__(self, lambda_squared=None, field=None):
        if lambda_squared is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

            R = PolynomialRing(ZZ, "x")
            x = R.gen()
            field = NumberField(
                x**3 - ZZ(5) * x**2 + ZZ(4) * x - ZZ(1), "r", embedding=AA(ZZ(4))
            )
            self._lambda_squared = field.gen()
        else:
            if field is None:
                self._lambda_squared = lambda_squared
                field = lambda_squared.parent()
            else:
                self._lambda_squared = field(lambda_squared)
        from flatsurf.geometry.categories import TranslationSurfaces

        super().__init__(
            field,
            category=TranslationSurfaces().InfiniteType().Connected().WithoutBoundary(),
        )

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space, i.e.,
        return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.e_infinity_surface()
            sage: S.is_compact()
            False

        """
        return False

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.e_infinity_surface()
            sage: S.is_mutable()
            False

        """
        return False

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.e_infinity_surface()
            sage: S.roots()
            (0,)

        """
        return (ZZ(0),)

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.e_infinity_surface()
            sage: S
            EInfinitySurface(r)

        """
        return f"EInfinitySurface({repr(self._lambda_squared)})"

    @cached_method
    def get_white(self, n):
        r"""Get the weight of the white endpoint of edge n."""
        if n == 0 or n == 1:
            return self._lambda_squared
        if n == -1:
            return self._lambda_squared - 1
        if n == 2:
            return 1 - 3 * self._lambda_squared + self._lambda_squared**2
        if n > 2:
            x = self.get_white(n - 1)
            y = self.get_black(n)
            return self._lambda_squared * y - x
        return self.get_white(-n)

    @cached_method
    def get_black(self, n):
        r"""Get the weight of the black endpoint of edge n."""
        if n == 0:
            return self.base_ring().one()
        if n == 1 or n == -1 or n == 2:
            return self._lambda_squared - 1
        if n > 2:
            x = self.get_black(n - 1)
            y = self.get_white(n - 1)
            return y - x
        return self.get_black(1 - n)

    def polygon(self, lab):
        r"""
        Return the polygon labeled by ``lab``.
        """
        if lab not in self.labels():
            raise ValueError(f"{lab=} not a valid label")
        return polygons.rectangle(2 * self.get_black(lab), self.get_white(lab))

    def labels(self):
        r"""
        Return the labels of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.e_infinity_surface()
            sage: S.labels()
            (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

        """
        from flatsurf.geometry.surface import LabelsFromView

        return LabelsFromView(self, ZZ, finite=False)

    def opposite_edge(self, p, e):
        r"""
        Return the pair ``(pp,ee)`` to which the edge ``(p,e)`` is glued to.
        """
        if p == 0:
            if e == 0:
                return (0, 2)
            if e == 1:
                return (1, 3)
            if e == 2:
                return (0, 0)
            if e == 3:
                return (1, 1)
        if p == 1:
            if e == 0:
                return (-1, 2)
            if e == 1:
                return (0, 3)
            if e == 2:
                return (2, 0)
            if e == 3:
                return (0, 1)
        if p == -1:
            if e == 0:
                return (2, 2)
            if e == 1:
                return (-1, 3)
            if e == 2:
                return (1, 0)
            if e == 3:
                return (-1, 1)
        if p == 2:
            if e == 0:
                return (1, 2)
            if e == 1:
                return (-2, 3)
            if e == 2:
                return (-1, 0)
            if e == 3:
                return (-2, 1)
        if p > 2:
            if e % 2:
                return -p, (e + 2) % 4
            else:
                return 1 - p, (e + 2) % 4
        else:
            if e % 2:
                return -p, (e + 2) % 4
            else:
                return 1 - p, (e + 2) % 4

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: hash(translation_surfaces.e_infinity_surface()) == hash(translation_surfaces.e_infinity_surface())
            True

        """
        return hash((self.base_ring(), self._lambda_squared))

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.e_infinity_surface()
            sage: S == S
            True

        """
        if not isinstance(other, EInfinitySurface):
            return False

        return (
            self._lambda_squared == other._lambda_squared
            and self.base_ring() == other.base_ring()
        )


class TFractalSurface(OrientedSimilaritySurface):
    r"""
    The TFractal surface.

    The TFractal surface is a translation surface of finite area built from
    infinitely many polygons. The basic building block is the following polygon::

         w/r    w     w/r
        +---+------+---+
        | 1 |   2  | 3 | h2
        +---+------+---+
            |   0  | h1
            +------+
                w

    where ``w``, ``h1``, ``h2``, ``r`` are some positive numbers. Default values
    are ``w=h1=h2=1`` and ``r=2``.

    .. TODO::

        In that surface, the linear flow can be computed more efficiently using
        only one affine interval exchange transformation with 5 intervals. But
        the underlying geometric construction is not a covering.

        Warning: we can not play at the same time with tuples and element of a
        cartesian product (see Sage trac ticket #19555)
    """

    def __init__(self, w=ZZ_1, r=ZZ_2, h1=ZZ_1, h2=ZZ_1):
        from sage.combinat.words.words import Words

        field = Sequence([w, r, h1, h2]).universe()
        if not field.is_field():
            field = field.fraction_field()

        from flatsurf.geometry.categories import TranslationSurfaces

        super().__init__(
            field,
            category=TranslationSurfaces()
            .InfiniteType()
            .WithoutBoundary()
            .Compact()
            .Connected(),
        )

        self._w = field(w)
        self._r = field(r)
        self._h1 = field(h1)
        self._h2 = field(h2)
        self._words = Words("LR", finite=True, infinite=False)
        self._wL = self._words("L")
        self._wR = self._words("R")

        self._root = (self._words(""), 0)

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.t_fractal()
            sage: S.roots()
            ((word: , 0),)

        """
        return (self._root,)

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.t_fractal()
            sage: S.is_mutable()
            False

        """
        return False

    def _repr_(self):
        return "The T-fractal surface with parameters w={}, r={}, h1={}, h2={}".format(
            self._w,
            self._r,
            self._h1,
            self._h2,
        )

    @cached_method
    def labels(self):
        r"""
        Return the labels of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.t_fractal()
            sage: S.labels()
            ((word: , 0), (word: , 2), (word: , 3), (word: , 1), (word: R, 2), (word: R, 0), (word: L, 2), (word: L, 0), (word: R, 3), (word: R, 1), (word: L, 3), (word: L, 1), (word: RR, 2), (word: RR, 0), (word: RL, 2), (word: RL, 0), …)

        """
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        from sage.categories.cartesian_product import cartesian_product

        labels = cartesian_product([self._words, FiniteEnumeratedSet([0, 1, 2, 3])])

        from flatsurf.geometry.surface import LabelsFromView

        return LabelsFromView(self, labels, finite=False)

    def opposite_edge(self, p, e):
        r"""
        Labeling of polygons::

             wl,0             wr,0
            +-----+---------+------+
            |     |         |      |
            | w,1 |   w,2   |  w,3 |
            |     |         |      |
            +-----+---------+------+
                  |         |
                  |   w,0   |
                  |         |
                  +---------+
                       w

        and we always have: bot->0, right->1, top->2, left->3

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: T = sfg.tfractal_surface()
            sage: W = T._words
            sage: w = W('LLRLRL')
            sage: T.opposite_edge((w,0),0)
            ((word: LLRLR, 1), 2)
            sage: T.opposite_edge((w,0),1)
            ((word: LLRLRL, 0), 3)
            sage: T.opposite_edge((w,0),2)
            ((word: LLRLRL, 2), 0)
            sage: T.opposite_edge((w,0),3)
            ((word: LLRLRL, 0), 1)
        """
        w, i = p
        w = self._words(w)
        i = int(i)
        e = int(e)

        def f():
            if e == 0:
                return 2
            if e == 1:
                return 3
            if e == 2:
                return 0
            if e == 3:
                return 1

            raise ValueError("e (={!r}) must be either 0,1,2 or 3".format(e))

        def label():
            if i == 0:
                if e == 0:
                    if w.is_empty():
                        return (w, 2)
                    elif w[-1] == "L":
                        return (w[:-1], 1)
                    elif w[-1] == "R":
                        return (w[:-1], 3)
                if e == 1:
                    return (w, 0)
                if e == 2:
                    return (w, 2)
                if e == 3:
                    return (w, 0)
            elif i == 1:
                if e == 0:
                    return (w + self._wL, 2)
                if e == 1:
                    return (w, 2)
                if e == 2:
                    return (w + self._wL, 0)
                if e == 3:
                    return (w, 3)
            elif i == 2:
                if e == 0:
                    return (w, 0)
                if e == 1:
                    return (w, 3)
                if e == 2:
                    if w.is_empty():
                        return (w, 0)
                    if w[-1] == "L":
                        return (w[:-1], 1)
                    if w[-1] == "R":
                        return (w[:-1], 3)
                if e == 3:
                    return (w, 1)
            elif i == 3:
                if e == 0:
                    return (w + self._wR, 2)
                if e == 1:
                    return (w, 1)
                if e == 2:
                    return (w + self._wR, 0)
                if e == 3:
                    return (w, 2)

            raise ValueError("i (={!r}) must be either 0,1,2 or 3".format(i))

        return label(), f()

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``::

             w/r         w/r
            +---+------+---+
            | 1 |  2   | 3 |
            |   |      |   |  h2
            +---+------+---+
                |  0   | h1
                +------+
                w

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: T = sfg.tfractal_surface()
            sage: T.polygon(('L',0))
            Polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: T.polygon(('LRL',0))
            Polygon(vertices=[(0, 0), (1/8, 0), (1/8, 1/8), (0, 1/8)])
        """
        w = self._words(lab[0])
        return (1 / self._r ** w.length()) * self._base_polygon(lab[1])

    @cached_method
    def _base_polygon(self, i):
        if i == 0:
            w = self._w
            h = self._h1
        elif i == 1 or i == 3:
            w = self._w / self._r
            h = self._h2
        elif i == 2:
            w = self._w
            h = self._h2
        else:
            raise ValueError("i must be one of 0, 1, 2, 3")

        return Polygon(
            base_ring=self.base_ring(), edges=[(w, 0), (0, h), (-w, 0), (0, -h)]
        )

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: hash(translation_surfaces.t_fractal()) == hash(translation_surfaces.t_fractal())
            True

        """
        return hash((self._w, self._h1, self._r, self._h2))

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: translation_surfaces.t_fractal() == translation_surfaces.t_fractal()
            True

        """
        if not isinstance(other, TFractalSurface):
            return False

        return (
            self._w == other._w
            and self._h1 == other._h1
            and self._r == other._r
            and self._h2 == other._h2
        )


def tfractal_surface(w=ZZ_1, r=ZZ_2, h1=ZZ_1, h2=ZZ_1):
    return TFractalSurface(w, r, h1, h2)


class SimilaritySurfaceGenerators:
    r"""
    Examples of similarity surfaces.
    """

    @staticmethod
    def example():
        r"""
        Construct a SimilaritySurface from a pair of triangles.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: ex = similarity_surfaces.example()
            sage: ex
            Genus 1 Surface built from 2 isosceles triangles

        TESTS::

            sage: TestSuite(ex).run()
            sage: from flatsurf.geometry.categories import SimilaritySurfaces
            sage: ex in SimilaritySurfaces()
            True

        """
        s = MutableOrientedSimilaritySurface(QQ)

        s.add_polygon(
            Polygon(vertices=[(0, 0), (2, -2), (2, 0)], base_ring=QQ), label=0
        )
        s.add_polygon(Polygon(vertices=[(0, 0), (2, 0), (1, 3)], base_ring=QQ), label=1)
        s.glue((0, 0), (1, 1))
        s.glue((0, 1), (1, 2))
        s.glue((0, 2), (1, 0))
        s.set_immutable()
        return s

    @staticmethod
    def self_glued_polygon(P):
        r"""
        Return the HalfTranslationSurface formed by gluing all edges of P to themselves.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces
            sage: p = Polygon(edges=[(2,0),(-1,3),(-1,-3)])
            sage: s = similarity_surfaces.self_glued_polygon(p)
            sage: s
            Half-Translation Surface in Q_0(-1^4) built from an isosceles triangle
            sage: TestSuite(s).run()

        """
        s = MutableOrientedSimilaritySurface(P.base_ring())
        s.add_polygon(P)
        for i in range(len(P.vertices())):
            s.glue((0, i), (0, i))
        s.set_immutable()
        return s

    @staticmethod
    def billiard(P, rational=None):
        r"""
        Return the ConeSurface associated to the billiard in the polygon ``P``.

        INPUT:

        - ``P`` -- a polygon

        - ``rational`` -- a boolean or ``None`` (default: ``None``) -- whether
          to assume that all the angles of ``P`` are a rational multiple of π.

        EXAMPLES::

            sage: from flatsurf import Polygon, similarity_surfaces

            sage: P = Polygon(vertices=[(0,0), (1,0), (0,1)])
            sage: Q = similarity_surfaces.billiard(P, rational=True)
            doctest:warning
            ...
            UserWarning: the rational keyword argument of billiard() has been deprecated and will be removed in a future version of sage-flatsurf; rationality checking is now faster so this is not needed anymore
            sage: Q
            Genus 0 Rational Cone Surface built from 2 isosceles triangles
            sage: from flatsurf.geometry.categories import ConeSurfaces
            sage: Q in ConeSurfaces().Rational()
            True
            sage: M = Q.minimal_cover(cover_type="translation")
            sage: M
            Minimal Translation Cover of Genus 0 Rational Cone Surface built from 2 isosceles triangles
            sage: TestSuite(M).run()
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: M in TranslationSurfaces()
            True

        A non-convex examples (L-shape polygon)::

            sage: P = Polygon(vertices=[(0,0), (2,0), (2,1), (1,1), (1,2), (0,2)])
            sage: Q = similarity_surfaces.billiard(P)
            sage: TestSuite(Q).run()
            sage: M = Q.minimal_cover(cover_type="translation")
            sage: TestSuite(M).run()
            sage: M.stratum()
            H_2(2, 0^5)

        A quadrilateral from Eskin-McMullen-Mukamel-Wright::

            sage: from flatsurf import Polygon
            sage: P = Polygon(angles=(1, 1, 1, 7))
            sage: S = similarity_surfaces.billiard(P)
            sage: TestSuite(S).run()
            sage: S = S.minimal_cover(cover_type="translation")
            sage: TestSuite(S).run()
            sage: S = S.erase_marked_points() # optional: pyflatsurf
            sage: TestSuite(S).run()
            sage: S, _ = S.normalized_coordinates()
            sage: TestSuite(S).run()

        Unfolding a triangle with non-algebraic lengths::

            sage: from flatsurf import EuclideanPolygonsWithAngles
            sage: E = EuclideanPolygonsWithAngles((3, 3, 5))
            sage: from pyexactreal import ExactReals # optional: pyexactreal  # random output due to cppyy deprecation warnings
            sage: R = ExactReals(E.base_ring()) # optional: pyexactreal
            sage: angles = (3, 3, 5)
            sage: slopes = EuclideanPolygonsWithAngles(*angles).slopes()
            sage: P = Polygon(angles=angles, edges=[R.random_element() * slopes[0]])  # optional: pyexactreal
            sage: S = similarity_surfaces.billiard(P); S # optional: pyexactreal
            Genus 0 Rational Cone Surface built from 2 isosceles triangles
            sage: TestSuite(S).run() # long time (6s), optional: pyexactreal
            sage: from flatsurf.geometry.categories import ConeSurfaces
            sage: S in ConeSurfaces()
            True

        """
        if not isinstance(P, EuclideanPolygon):
            raise TypeError("invalid input")

        if rational is not None:
            import warnings

            warnings.warn(
                "the rational keyword argument of billiard() has been deprecated and will be removed in a future version of sage-flatsurf; rationality checking is now faster so this is not needed anymore"
            )

        from flatsurf.geometry.categories import ConeSurfaces

        category = ConeSurfaces()
        if P.is_rational():
            category = category.Rational()

        V = P.base_ring() ** 2

        if not P.is_convex():
            # triangulate non-convex ones
            base_ring = P.base_ring()
            comb_edges = P.triangulation()
            vertices = P.vertices()
            comb_triangles = SimilaritySurfaceGenerators._billiard_build_faces(
                len(vertices), comb_edges
            )
            triangles = []
            internal_edges = []  # list (p1, e1, p2, e2)
            external_edges = []  # list (p1, e1)
            edge_to_lab = {}
            for num, (i, j, k) in enumerate(comb_triangles):
                triangles.append(
                    Polygon(
                        vertices=[vertices[i], vertices[j], vertices[k]],
                        base_ring=base_ring,
                    )
                )
                edge_to_lab[(i, j)] = (num, 0)
                edge_to_lab[(j, k)] = (num, 1)
                edge_to_lab[(k, i)] = (num, 2)
            for num, (i, j, k) in enumerate(comb_triangles):
                if (j, i) in edge_to_lab:
                    num2, e2 = edge_to_lab[j, i]
                    internal_edges.append((num, 0, num2, e2))
                else:
                    external_edges.append((num, 0))
                if (k, j) in edge_to_lab:
                    num2, e2 = edge_to_lab[k, j]
                    internal_edges.append((num, 1, num2, e2))
                else:
                    external_edges.append((num, 1))
                if (i, k) in edge_to_lab:
                    num2, e2 = edge_to_lab[i, k]
                    internal_edges.append((num, 2, num2, e2))
                else:
                    external_edges.append((num, 1))
            P = triangles
        else:
            internal_edges = []
            external_edges = [(0, i) for i in range(len(P.vertices()))]
            base_ring = P.base_ring()
            P = [P]

        surface = MutableOrientedSimilaritySurface(base_ring, category=category)

        m = len(P)

        for p in P:
            surface.add_polygon(p)
        for p in P:
            surface.add_polygon(
                Polygon(edges=[V((-x, y)) for x, y in reversed(p.edges())])
            )
        for p1, e1, p2, e2 in internal_edges:
            surface.glue((p1, e1), (p2, e2))
            ne1 = len(surface.polygon(p1).vertices())
            ne2 = len(surface.polygon(p2).vertices())
            surface.glue((m + p1, ne1 - e1 - 1), (m + p2, ne2 - e2 - 1))
        for p, e in external_edges:
            ne = len(surface.polygon(p).vertices())
            surface.glue((p, e), (m + p, ne - e - 1))

        surface.set_immutable()

        return surface

    @staticmethod
    def _billiard_build_faces(n, edges):
        r"""
        Given a combinatorial list of pairs ``edges`` forming a cell-decomposition
        of a polygon (with vertices labeled from ``0`` to ``n-1``) return the list
        of cells.

        This is a helper method for :meth:`billiard`.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: SimilaritySurfaceGenerators._billiard_build_faces(4, [(0,2)])
            [[0, 1, 2], [2, 3, 0]]
            sage: SimilaritySurfaceGenerators._billiard_build_faces(4, [(1,3)])
            [[1, 2, 3], [3, 0, 1]]
            sage: SimilaritySurfaceGenerators._billiard_build_faces(5, [(0,2), (0,3)])
            [[0, 1, 2], [3, 4, 0], [0, 2, 3]]
            sage: SimilaritySurfaceGenerators._billiard_build_faces(5, [(0,2)])
            [[0, 1, 2], [2, 3, 4, 0]]
            sage: SimilaritySurfaceGenerators._billiard_build_faces(5, [(1,4)])
            [[1, 2, 3, 4], [4, 0, 1]]
            sage: SimilaritySurfaceGenerators._billiard_build_faces(5, [(1,3),(3,0)])
            [[1, 2, 3], [3, 4, 0], [0, 1, 3]]
        """
        polygons = [list(range(n))]
        for u, v in edges:
            j = None
            for i, p in enumerate(polygons):
                if u in p and v in p:
                    if j is not None:
                        raise RuntimeError
                    j = i
            if j is None:
                raise RuntimeError
            p = polygons[j]
            i0 = p.index(u)
            i1 = p.index(v)
            if i0 > i1:
                i0, i1 = i1, i0
            polygons[j] = p[i0 : i1 + 1]
            polygons.append(p[i1:] + p[: i0 + 1])
        return polygons

    @staticmethod
    def polygon_double(P):
        r"""
        Return the ConeSurface associated to the billiard in the polygon ``P``.
        Differs from billiard(P) only in the graphical display. Here, we display
        the polygons separately.
        """
        from sage.matrix.constructor import matrix

        n = len(P.vertices())
        r = matrix(2, [-1, 0, 0, 1])
        Q = Polygon(edges=[r * v for v in reversed(P.edges())])

        surface = MutableOrientedSimilaritySurface(P.base_ring())
        surface.add_polygon(P, label=0)
        surface.add_polygon(Q, label=1)
        for i in range(n):
            surface.glue((0, i), (1, n - i - 1))
        surface.set_immutable()
        return surface

    @staticmethod
    def right_angle_triangle(w, h):
        r"""
        TESTS::

            sage: from flatsurf import similarity_surfaces
            sage: R = similarity_surfaces.right_angle_triangle(2, 3)
            sage: R
            Genus 0 Cone Surface built from 2 right triangles
            sage: from flatsurf.geometry.categories import ConeSurfaces
            sage: R in ConeSurfaces()
            True
            sage: TestSuite(R).run()
        """
        F = Sequence([w, h]).universe()

        if not F.is_field():
            F = F.fraction_field()
        V = VectorSpace(F, 2)
        s = MutableOrientedSimilaritySurface(F)
        s.add_polygon(
            Polygon(base_ring=F, edges=[V((w, 0)), V((-w, h)), V((0, -h))]), label=0
        )
        s.add_polygon(
            Polygon(base_ring=F, edges=[V((0, h)), V((-w, -h)), V((w, 0))]), label=1
        )
        s.glue((0, 0), (1, 2))
        s.glue((0, 1), (1, 1))
        s.glue((0, 2), (1, 0))
        s.set_immutable()
        return s


similarity_surfaces = SimilaritySurfaceGenerators()


class DilationSurfaceGenerators:
    @staticmethod
    def basic_dilation_torus(a):
        r"""
        Return a dilation torus built from a `1 \times 1` square and a `a
        \times 1` rectangle. Each edge of the square is glued to the opposite
        edge of the rectangle. This results in horizontal edges glued by a
        dilation with a scaling factor of a, and vertical edges being glued by
        translation::

                b       a
              +----+---------+
              | 0  | 1       |
            c |    |         | c
              +----+---------+
                a       b

        EXAMPLES::

            sage: from flatsurf import dilation_surfaces
            sage: ds = dilation_surfaces.basic_dilation_torus(AA(sqrt(2)))
            sage: ds
            Genus 1 Positive Dilation Surface built from a square and a rectangle
            sage: from flatsurf.geometry.categories import DilationSurfaces
            sage: ds in DilationSurfaces().Positive()
            True
            sage: TestSuite(ds).run()

        """
        s = MutableOrientedSimilaritySurface(a.parent().fraction_field())
        s.add_polygon(
            Polygon(base_ring=s.base_ring(), edges=[(0, 1), (-1, 0), (0, -1), (1, 0)]),
            label=0,
        )
        s.add_polygon(
            Polygon(base_ring=s.base_ring(), edges=[(0, 1), (-a, 0), (0, -1), (a, 0)]),
            label=1,
        )
        # label 1
        s.glue((0, 0), (1, 2))
        s.glue((0, 1), (1, 3))
        s.glue((0, 2), (1, 0))
        s.glue((0, 3), (1, 1))
        s.set_roots([0])
        s.set_immutable()
        return s

    @staticmethod
    def genus_two_square(a, b, c, d):
        r"""
        A genus two dilation surface is returned.

        The unit square is made into an octagon by marking a point on
        each of its edges. Then opposite sides of this octagon are
        glued together by translation. (Since we currently require strictly
        convex polygons, we subdivide the square into a hexagon and two
        triangles as depicted below.) The parameters ``a``, ``b``, ``c``, and
        ``d`` should be real numbers strictly between zero and one. These
        represent the lengths of an edge of the resulting octagon, as below::

                    c
              +--+-------+
            d |2/        |
              |/         |
              +    0     +
              |         /|
              |        /1| b
              +-------+--+
                 a

        The other edges will have length `1-a`, `1-b`, `1-c`, and `1-d`.
        Dilations used to glue edges will be by factors `c/a`, `d/b`,
        `(1-c)/(1-a)` and `(1-d)/(1-b)`.

        EXAMPLES::

            sage: from flatsurf import dilation_surfaces
            sage: ds = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
            sage: ds
            Genus 2 Positive Dilation Surface built from 2 right triangles and a hexagon
            sage: from flatsurf.geometry.categories import DilationSurfaces
            sage: ds in DilationSurfaces().Positive()
            True
            sage: TestSuite(ds).run()
        """
        field = Sequence([a, b, c, d]).universe().fraction_field()
        s = MutableOrientedSimilaritySurface(QQ)

        hexagon = Polygon(
            edges=[(a, 0), (1 - a, b), (0, 1 - b), (-c, 0), (c - 1, -d), (0, d - 1)],
            base_ring=field,
        )
        s.add_polygon(hexagon, label=0)
        s.set_roots([0])
        triangle1 = Polygon(base_ring=field, edges=[(1 - a, 0), (0, b), (a - 1, -b)])
        s.add_polygon(triangle1, label=1)
        triangle2 = Polygon(base_ring=field, edges=[(1 - c, d), (c - 1, 0), (0, -d)])
        s.add_polygon(triangle2, label=2)
        s.glue((0, 0), (0, 3))
        s.glue((0, 2), (0, 5))
        s.glue((0, 1), (1, 2))
        s.glue((0, 4), (2, 0))
        s.glue((1, 0), (2, 1))
        s.glue((1, 1), (2, 2))
        s.set_immutable()
        return s


dilation_surfaces = DilationSurfaceGenerators()


class HalfTranslationSurfaceGenerators:
    # TODO: ideally, we should be able to construct a non-convex polygon and make the construction
    # below as a special case of billiard unfolding.
    @staticmethod
    def step_billiard(w, h):
        r"""
        Return a (finite) step billiard associated to the given widths ``w`` and heights ``h``.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1,1,1,1], [1,1/2,1/3,1/5])
            sage: S
            StepBilliard(w=[1, 1, 1, 1], h=[1, 1/2, 1/3, 1/5])
            sage: from flatsurf.geometry.categories import DilationSurfaces
            sage: S in DilationSurfaces()
            True
            sage: TestSuite(S).run()
        """
        n = len(h)
        if len(w) != n:
            raise ValueError
        if n < 2:
            raise ValueError("w and h must have length at least 2")
        H = sum(h)
        W = sum(w)

        R = Sequence(w + h).universe()
        base_ring = R.fraction_field()

        P = []
        Prev = []
        x = 0
        y = H
        for i in range(n - 1):
            P.append(
                Polygon(
                    vertices=[
                        (x, 0),
                        (x + w[i], 0),
                        (x + w[i], y - h[i]),
                        (x + w[i], y),
                        (x, y),
                    ],
                    base_ring=base_ring,
                )
            )
            x += w[i]
            y -= h[i]
        assert x == W - w[-1]
        assert y == h[-1]
        P.append(
            Polygon(
                vertices=[(x, 0), (x + w[-1], 0), (x + w[-1], y), (x, y)],
                base_ring=base_ring,
            )
        )

        Prev = [
            Polygon(
                vertices=[(x, -y) for x, y in reversed(p.vertices())],
                base_ring=base_ring,
            )
            for p in P
        ]

        S = MutableOrientedSimilaritySurface(base_ring)
        S.rename(
            "StepBilliard(w=[{}], h=[{}])".format(
                ", ".join(map(str, w)), ", ".join(map(str, h))
            )
        )
        for p in P:
            S.add_polygon(p)  # get labels 0, ..., n-1
        for p in Prev:
            S.add_polygon(p)  # get labels n, n+1, ..., 2n-1

        # reflection gluings
        # (gluings between the polygon and its reflection)
        S.glue((0, 4), (n, 4))
        S.glue((n - 1, 0), (2 * n - 1, 2))
        S.glue((n - 1, 1), (2 * n - 1, 1))
        S.glue((n - 1, 2), (2 * n - 1, 0))
        for i in range(n - 1):
            # glue((polygon1, edge1), (polygon2, edge2))
            S.glue((i, 0), (n + i, 3))
            S.glue((i, 2), (n + i, 1))
            S.glue((i, 3), (n + i, 0))

        # translation gluings
        S.glue((n - 2, 1), (n - 1, 3))
        S.glue((2 * n - 2, 2), (2 * n - 1, 3))
        for i in range(n - 2):
            S.glue((i, 1), (i + 1, 4))
            S.glue((n + i, 2), (n + i + 1, 4))

        S.set_immutable()
        return S


half_translation_surfaces = HalfTranslationSurfaceGenerators()


class TranslationSurfaceGenerators:
    r"""
    Common and less common translation surfaces.
    """

    @staticmethod
    def square_torus(a=1):
        r"""
        Return flat torus obtained by identification of the opposite sides of a
        square.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: T
            Translation Surface in H_1(0) built from a square
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: T in TranslationSurfaces()
            True
            sage: TestSuite(T).run()

        Rational directions are completely periodic::

            sage: v = T.tangent_vector(0, (1/33, 1/257), (13,17))
            sage: L = v.straight_line_trajectory()
            sage: L.flow(13+17)
            sage: L.is_closed()
            True

        TESTS::

            sage: TestSuite(T).run()
        """
        return TranslationSurfaceGenerators.torus((a, 0), (0, a))

    @staticmethod
    def torus(u, v):
        r"""
        Return the flat torus obtained as the quotient of the plane by the pair
        of vectors ``u`` and ``v``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, AA(2).sqrt()), (AA(3).sqrt(), 3))
            sage: T
            Translation Surface in H_1(0) built from a quadrilateral
            sage: T.polygon(0)
            Polygon(vertices=[(0, 0), (1, 1.414213562373095?), (2.732050807568878?, 4.414213562373095?), (1.732050807568878?, 3)])
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: T in TranslationSurfaces()
            True

        """
        u = vector(u)
        v = vector(v)
        field = Sequence([u, v]).universe().base_ring()
        if isinstance(field, type):
            field = py_scalar_parent(field)
        if not field.is_field():
            field = field.fraction_field()
        s = MutableOrientedSimilaritySurface(field)
        p = Polygon(vertices=[(0, 0), u, u + v, v], base_ring=field)
        s.add_polygon(p)
        s.glue((0, 0), (0, 2))
        s.glue((0, 1), (0, 3))
        s.set_immutable()
        return s

    @staticmethod
    def veech_2n_gon(n):
        r"""
        The regular 2n-gon with opposite sides identified.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.veech_2n_gon(5)
            sage: s
            Translation Surface in H_2(1^2) built from a regular decagon
            sage: TestSuite(s).run()

        """
        p = polygons.regular_ngon(2 * n)
        s = MutableOrientedSimilaritySurface(p.base_ring())
        s.add_polygon(p)
        for i in range(2 * n):
            s.glue((0, i), (0, (i + n) % (2 * n)))
        s.set_immutable()
        return s

    @staticmethod
    def veech_double_n_gon(n):
        r"""
        A pair of regular n-gons with each edge of one identified to an edge of the other to make a translation surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s=translation_surfaces.veech_double_n_gon(5)
            sage: s
            Translation Surface in H_2(2) built from 2 regular pentagons
            sage: TestSuite(s).run()

        """
        from sage.matrix.constructor import Matrix

        p = polygons.regular_ngon(n)
        s = MutableOrientedSimilaritySurface(p.base_ring())
        m = Matrix([[-1, 0], [0, -1]])
        s.add_polygon(p, label=0)
        s.add_polygon(m * p, label=1)
        for i in range(n):
            s.glue((0, i), (1, i))
        s.set_immutable()
        return s

    @staticmethod
    def regular_octagon():
        r"""
        Return the translation surface built from the regular octagon by
        identifying opposite sides.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.regular_octagon()
            sage: T
            Translation Surface in H_2(2) built from a regular octagon
            sage: TestSuite(T).run()
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: T in TranslationSurfaces()
            True
        """
        return translation_surfaces.veech_2n_gon(4)

    @staticmethod
    def mcmullen_genus2_prototype(w, h, t, e, rel=0, base_ring=None):
        r"""
        McMullen prototypes in the stratum H(2).

        These prototype appear at least in McMullen "Teichmüller curves in genus
        two: Discriminant and spin" (2004). The notation from that paper are
        quadruple ``(a, b, c, e)`` which translates in our notation as
        ``w = b``, ``h = c``, ``t = a`` (and ``e = e``).

        The associated discriminant is `D = e^2 + 4 wh`.

        If ``rel`` is a positive parameter (less than w-lambda) the surface belongs
        to the eigenform locus in H(1,1).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from surface_dynamics import Stratum

            sage: prototypes = {
            ....:      5: [(1,1,0,-1)],
            ....:      8: [(1,1,0,-2), (2,1,0,0)],
            ....:      9: [(2,1,0,-1)],
            ....:     12: [(1,2,0,-2), (2,1,0,-2), (3,1,0,0)],
            ....:     13: [(1,1,0,-3), (3,1,0,-1), (3,1,0,1)],
            ....:     16: [(3,1,0,-2), (4,1,0,0)],
            ....:     17: [(1,2,0,-3), (2,1,0,-3), (2,2,0,-1), (2,2,1,-1), (4,1,0,-1), (4,1,0,1)],
            ....:     20: [(1,1,0,-4), (2,2,1,-2), (4,1,0,-2), (4,1,0,2)],
            ....:     21: [(1,3,0,-3), (3,1,0,-3)],
            ....:     24: [(1,2,0,-4), (2,1,0,-4), (3,2,0,0)],
            ....:     25: [(2,2,0,-3), (2,2,1,-3), (3,2,0,-1), (4,1,0,-3)]}

            sage: for D in sorted(prototypes):  # long time (.5s)
            ....:     for w,h,t,e in prototypes[D]:
            ....:          T = translation_surfaces.mcmullen_genus2_prototype(w,h,t,e)
            ....:          assert T.stratum() == Stratum([2], 1)
            ....:          assert (D.is_square() and T.base_ring() is QQ) or (T.base_ring().polynomial().discriminant() == D)

        An example with some relative homology::

            sage: U8 = translation_surfaces.mcmullen_genus2_prototype(2,1,0,0,1/4)    # discriminant 8
            sage: U8
            Translation Surface in H_2(1^2) built from a rectangle and a quadrilateral
            sage: U12 = translation_surfaces.mcmullen_genus2_prototype(3,1,0,0,3/10)   # discriminant 12
            sage: U12
            Translation Surface in H_2(1^2) built from a rectangle and a quadrilateral

            sage: U8.stratum()
            H_2(1^2)
            sage: U8.base_ring().polynomial().discriminant()
            8
            sage: U8.j_invariant()
            (
                      [4 0]
            (0), (0), [0 2]
            )

            sage: U12.stratum()
            H_2(1^2)
            sage: U12.base_ring().polynomial().discriminant()
            12
            sage: U12.j_invariant()
            (
                      [6 0]
            (0), (0), [0 2]
            )
        """
        w = ZZ(w)
        h = ZZ(h)
        t = ZZ(t)
        e = ZZ(e)
        g = w.gcd(h)
        if (
            w <= 0
            or h <= 0
            or t < 0
            or t >= g
            or not g.gcd(t).gcd(e).is_one()
            or e + h >= w
        ):
            raise ValueError("invalid parameters")

        x = polygen(QQ)
        poly = x**2 - e * x - w * h
        if poly.is_irreducible():
            if base_ring is None:
                emb = AA.polynomial_root(poly, RIF(0, w))
                K = NumberField(poly, "l", embedding=emb)
                λ = K.gen()
            else:
                K = base_ring
                roots = poly.roots(K, multiplicities=False)
                if len(roots) != 2:
                    raise ValueError("invalid base ring")
                roots.sort(key=lambda x: x.numerical_approx())
                assert roots[0] < 0 and roots[0] > 0
                λ = roots[1]
        else:
            if base_ring is None:
                K = QQ
            else:
                K = base_ring
            D = e**2 + 4 * w * h
            d = D.sqrt()
            λ = (e + d) / 2

        try:
            rel = K(rel)
        except TypeError:
            K = get_coercion_model().common_parent(K, parent(rel))
            λ = K(λ)
            rel = K(rel)

        # (lambda,lambda) square on top
        # twisted (w,0), (t,h)
        s = MutableOrientedSimilaritySurface(K)
        if rel:
            if rel < 0 or rel > w - λ:
                raise ValueError("invalid rel argument")
            s.add_polygon(
                Polygon(vertices=[(0, 0), (λ, 0), (λ + rel, λ), (rel, λ)], base_ring=K)
            )
            s.add_polygon(
                Polygon(
                    vertices=[
                        (0, 0),
                        (rel, 0),
                        (rel + λ, 0),
                        (w, 0),
                        (w + t, h),
                        (λ + rel + t, h),
                        (t + λ, h),
                        (t, h),
                    ],
                    base_ring=K,
                )
            )
            s.glue((0, 1), (0, 3))
            s.glue((0, 0), (1, 6))
            s.glue((0, 2), (1, 1))
            s.glue((1, 2), (1, 4))
            s.glue((1, 3), (1, 7))
            s.glue((1, 0), (1, 5))
        else:
            s.add_polygon(
                Polygon(vertices=[(0, 0), (λ, 0), (λ, λ), (0, λ)], base_ring=K)
            )
            s.add_polygon(
                Polygon(
                    vertices=[(0, 0), (λ, 0), (w, 0), (w + t, h), (λ + t, h), (t, h)],
                    base_ring=K,
                )
            )
            s.glue((0, 1), (0, 3))
            s.glue((0, 0), (1, 4))
            s.glue((0, 2), (1, 0))
            s.glue((1, 1), (1, 3))
            s.glue((1, 2), (1, 5))
        s.set_immutable()
        return s

    @staticmethod
    def lanneau_nguyen_genus3_prototype(w, h, t, e):
        r"""
        Return the Lanneau--Ngyuen prototype in the stratum H(4) with
        parameters ``w``, ``h``, ``t``, and ``e``.

        The associated discriminant is `e^2 + 8 wh`.

        INPUT:

        - ``w`` -- a positive integer

        - ``h`` -- a positive integer

        - ``t`` -- a non-negative integer strictly less than the gcd of ``w``
          and ``h``

        - ``e`` -- an integer strictly less than ``w - 2h``

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.lanneau_nguyen_genus3_prototype(2,1,0,-1)
            sage: T.stratum()
            H_3(4)

        REFERENCES:

        - Erwan Lanneau, and Duc-Manh Nguyen, "Teichmüller curves generated by
          Weierstrass Prym eigenforms in genus 3 and genus 4", Journal of
          Topology 7.2, 2014, pp.475-522

        """
        from sage.all import gcd

        w = ZZ(w)
        h = ZZ(h)
        t = ZZ(t)
        e = ZZ(e)

        if not w > 0:
            raise ValueError("w must be positive")
        if not h > 0:
            raise ValueError("h must be positive")
        if not e + 2 * h < w:
            raise ValueError("e + 2h < w must hold")
        if not t >= 0:
            raise ValueError("t must be non-negative")
        if not t < gcd(w, h):
            raise ValueError("t must be smaller than the gcd of w and h")
        if not gcd([w, h, t, e]) == 1:
            raise ValueError("w, h, t, e must be coprime")

        K = QQ
        D = K(e**2 + 8 * w * h)

        if not D.is_square():
            from sage.all import QuadraticField

            K = QuadraticField(D)
            D = K(D)

        d = D.sqrt()

        λ = (e + d) / 2

        # (λ, λ) square in the middle
        # twisted (w, 0), (t, h) on top and bottom
        s = MutableOrientedSimilaritySurface(K)

        from flatsurf import Polygon

        s.add_polygon(
            Polygon(
                vertices=[(0, 0), (λ, 0), (w, 0), (w + t, h), (λ + t, h), (t, h)],
                base_ring=K,
            )
        )
        s.add_polygon(Polygon(vertices=[(0, 0), (λ, 0), (λ, λ), (0, λ)], base_ring=K))
        s.add_polygon(
            Polygon(
                vertices=[
                    (0, 0),
                    (w - λ, 0),
                    (w, 0),
                    (w + t, h),
                    (w - λ + t, h),
                    (t, h),
                ],
                base_ring=K,
            )
        )
        s.glue((0, 0), (2, 3))
        s.glue((0, 1), (0, 3))
        s.glue((0, 2), (0, 5))
        s.glue((0, 4), (1, 0))
        s.glue((1, 1), (1, 3))
        s.glue((1, 2), (2, 1))
        s.glue((2, 0), (2, 4))
        s.glue((2, 2), (2, 5))
        s.set_immutable()
        return s

    @staticmethod
    def lanneau_nguyen_genus4_prototype(w, h, t, e):
        r"""
        Return the Lanneau--Ngyuen prototype in the stratum H(6) with
        parameters ``w``, ``h``, ``t``, and ``e``.

        The associated discriminant is `D = e^2 + 4 wh`.

        INPUT:

        - ``w`` -- a positive integer

        - ``h`` -- a positive integer

        - ``t`` -- a non-negative integer strictly smaller than the gcd of
          ``w`` and ``h``

        - ``e`` -- a positive integer strictly smaller than `2w - \sqrt{D}`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.lanneau_nguyen_genus4_prototype(1,1,0,-2)
            sage: T.stratum()
            H_4(6)

        REFERENCES:

        - Erwan Lanneau, and Duc-Manh Nguyen, "Weierstrass Prym eigenforms in
          genus four", Journal of the Institute of Mathematics of Jussieu,
          19.6, 2020, pp.2045-2085

        """
        from sage.all import gcd

        w = ZZ(w)
        h = ZZ(h)
        t = ZZ(t)
        e = ZZ(e)

        if not w > 0:
            raise ValueError("w must be positive")
        if not h > 0:
            raise ValueError("h must be positive")
        if not t >= 0:
            raise ValueError("t must be non-negative")
        if not t < gcd(w, h):
            raise ValueError("t must be smaller than the gcd of w and h")
        if not gcd([w, h, t, e]) == 1:
            raise ValueError("w, h, t, e must be coprime")

        K = QQ
        D = K(e**2 + 4 * w * h)

        if not D.is_square():
            from sage.all import QuadraticField

            K = QuadraticField(D)
            D = K(D)

        d = D.sqrt()

        λ = (e + d) / 2

        if not λ < w:
            raise ValueError("λ must be smaller than w")

        if λ == w / 2:
            raise ValueError("λ and w/2 must be distinct")

        # (λ/2, λ/2) squares on top and bottom
        # twisted (w/2, 0), (t/2, h/2) in the middle
        s = MutableOrientedSimilaritySurface(K)

        from flatsurf import Polygon

        s.add_polygon(
            Polygon(
                vertices=[(0, 0), (λ / 2, 0), (λ / 2, λ / 2), (0, λ / 2)], base_ring=K
            )
        )
        s.add_polygon(
            Polygon(
                vertices=[
                    (0, 0),
                    (w / 2 - λ, 0),
                    (w / 2 - λ / 2, 0),
                    (w / 2, 0),
                    (w / 2 + t / 2, h / 2),
                    (w / 2 - λ + t / 2, h / 2),
                    (t / 2, h / 2),
                ],
                base_ring=K,
            )
        )
        s.add_polygon(
            Polygon(
                vertices=[
                    (0, 0),
                    (λ, 0),
                    (w / 2, 0),
                    (w / 2 + t / 2, h / 2),
                    (λ + t / 2, h / 2),
                    (λ / 2 + t / 2, h / 2),
                    (t / 2, h / 2),
                ],
                base_ring=K,
            )
        )
        s.add_polygon(
            Polygon(
                vertices=[(0, 0), (λ / 2, 0), (λ / 2, λ / 2), (0, λ / 2)], base_ring=K
            )
        )
        s.glue((0, 0), (2, 4))
        s.glue((0, 1), (0, 3))
        s.glue((0, 2), (1, 2))
        s.glue((1, 0), (1, 5))
        s.glue((1, 1), (3, 2))
        s.glue((1, 3), (1, 6))
        s.glue((1, 4), (2, 0))
        s.glue((2, 1), (2, 3))
        s.glue((2, 2), (2, 6))
        s.glue((2, 5), (3, 0))
        s.glue((3, 1), (3, 3))
        s.set_immutable()
        return s

    @staticmethod
    def mcmullen_L(l1, l2, l3, l4):
        r"""
        Return McMullen's L shaped surface with parameters l1, l2, l3, l4.

        Polygon labels and lengths are marked below::

            +-----+
            |     |
            |  1  |l1
            |     |
            |     |    l4
            +-----+---------+
            |     |         |
            |  0  |    2    |l2
            |     |         |
            +-----+---------+
              l3

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.mcmullen_L(1,1,1,1)
            sage: s
            Translation Surface in H_2(2) built from 3 squares
            sage: TestSuite(s).run()

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: L = translation_surfaces.mcmullen_L(1r, 1r, 1r, 1r)
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: L in TranslationSurfaces()
            True
        """
        field = Sequence([l1, l2, l3, l4]).universe()
        if isinstance(field, type):
            field = py_scalar_parent(field)
        if not field.is_field():
            field = field.fraction_field()

        s = MutableOrientedSimilaritySurface(field)
        s.add_polygon(
            Polygon(edges=[(l3, 0), (0, l2), (-l3, 0), (0, -l2)], base_ring=field)
        )
        s.add_polygon(
            Polygon(edges=[(l3, 0), (0, l1), (-l3, 0), (0, -l1)], base_ring=field)
        )
        s.add_polygon(
            Polygon(edges=[(l4, 0), (0, l2), (-l4, 0), (0, -l2)], base_ring=field)
        )
        s.glue((0, 0), (1, 2))
        s.glue((0, 1), (2, 3))
        s.glue((0, 2), (1, 0))
        s.glue((0, 3), (2, 1))
        s.glue((1, 1), (1, 3))
        s.glue((2, 0), (2, 2))
        s.set_immutable()
        return s

    @staticmethod
    def ward(n):
        r"""
        Return the surface formed by gluing a regular 2n-gon to two regular n-gons.
        These surfaces have Veech's lattice property due to work of Ward.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.ward(3)
            sage: s
            Translation Surface in H_1(0^3) built from 2 equilateral triangles and a regular hexagon
            sage: TestSuite(s).run()
            sage: s = translation_surfaces.ward(7)
            sage: s
            Translation Surface in H_6(10) built from 2 regular heptagons and a regular tetradecagon
            sage: TestSuite(s).run()
        """
        if n < 3:
            raise ValueError
        o = ZZ_2 * polygons.regular_ngon(2 * n)
        p1 = Polygon(edges=[o.edge((2 * i + n) % (2 * n)) for i in range(n)])
        p2 = Polygon(edges=[o.edge((2 * i + n + 1) % (2 * n)) for i in range(n)])
        s = MutableOrientedSimilaritySurface(o.base_ring())
        s.add_polygon(o)
        s.add_polygon(p1)
        s.add_polygon(p2)
        for i in range(n):
            s.glue((1, i), (0, 2 * i))
            s.glue((2, i), (0, 2 * i + 1))
        s.set_immutable()
        return s

    @staticmethod
    def octagon_and_squares():
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: os = translation_surfaces.octagon_and_squares()
            sage: os
            Translation Surface in H_3(4) built from 2 squares and a regular octagon
            sage: TestSuite(os).run()
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: os in TranslationSurfaces()
            True
        """
        return translation_surfaces.ward(4)

    @staticmethod
    def cathedral(a, b):
        r"""
        Return the cathedral surface with parameters ``a`` and ``b``.

        For any parameter ``a`` and ``b``, the cathedral surface belongs to the
        so-called Gothic locus described in McMullen, Mukamel, Wright "Cubic
        curves and totally geodesic subvarieties of moduli space" (2017)::

                     1
                   <--->

                    /\           2a
                   /  \      +------+
               a  b|   | a  /        \
             +----+    +---+          +
             |    |    |   |          |
            1| P0 |P1  |P2 |  P3      |
             |    |    |   |          |
             +----+    +---+          +
                 b|    |    \        /
                   \  /      +------+
                    \/

        If a and b satisfies

        .. MATH::

            a = x + y \sqrt(d) \qquad b = -3x -3/2 + 3y \sqrt(d)

        for some rational x,y and d >= 0 then it is a Teichmüller curve.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.cathedral(1,2)
            sage: C
            Translation Surface in H_4(2^3) built from 2 squares, a hexagon with 4 marked vertices and an octagon
            sage: TestSuite(C).run()

            sage: from pyexactreal import ExactReals # optional: pyexactreal
            sage: K = QuadraticField(5, embedding=AA(5).sqrt())
            sage: R = ExactReals(K) # optional: pyexactreal
            sage: C = translation_surfaces.cathedral(K.gen(), R.random_element([0.1, 0.2])) # optional: pyexactreal
            sage: C  # optional: pyexactreal
            Translation Surface in H_4(2^3) built from 2 rectangles, a hexagon with 4 marked vertices and an octagon
            sage: C.stratum() # optional: pyexactreal
            H_4(2^3)
            sage: TestSuite(C).run() # long time (6s), optional: pyexactreal
        """
        ring = Sequence([a, b]).universe()
        if isinstance(ring, type):
            ring = py_scalar_parent(ring)
        if not ring.has_coerce_map_from(QQ):
            ring = ring.fraction_field()
        a = ring(a)
        b = ring(b)
        s = MutableOrientedSimilaritySurface(ring)
        half = QQ((1, 2))
        p0 = Polygon(base_ring=ring, vertices=[(0, 0), (a, 0), (a, 1), (0, 1)])
        p1 = Polygon(
            base_ring=ring,
            vertices=[
                (a, 0),
                (a, -b),
                (a + half, -b - half),
                (a + 1, -b),
                (a + 1, 0),
                (a + 1, 1),
                (a + 1, b + 1),
                (a + half, b + 1 + half),
                (a, b + 1),
                (a, 1),
            ],
        )
        p2 = Polygon(
            base_ring=ring,
            vertices=[(a + 1, 0), (2 * a + 1, 0), (2 * a + 1, 1), (a + 1, 1)],
        )
        p3 = Polygon(
            base_ring=ring,
            vertices=[
                (2 * a + 1, 0),
                (2 * a + 1 + half, -half),
                (4 * a + 1 + half, -half),
                (4 * a + 2, 0),
                (4 * a + 2, 1),
                (4 * a + 1 + half, 1 + half),
                (2 * a + 1 + half, 1 + half),
                (2 * a + 1, 1),
            ],
        )
        s.add_polygon(p0)
        s.add_polygon(p1)
        s.add_polygon(p2)
        s.add_polygon(p3)
        s.glue((0, 0), (0, 2))
        s.glue((0, 1), (1, 9))
        s.glue((0, 3), (3, 3))
        s.glue((1, 0), (1, 3))
        s.glue((1, 1), (3, 4))
        s.glue((1, 2), (3, 6))
        s.glue((1, 4), (2, 3))
        s.glue((1, 5), (1, 8))
        s.glue((1, 6), (3, 0))
        s.glue((1, 7), (3, 2))
        s.glue((2, 0), (2, 2))
        s.glue((2, 1), (3, 7))
        s.glue((3, 1), (3, 5))
        s.set_immutable()
        return s

    @staticmethod
    def arnoux_yoccoz(genus):
        r"""
        Construct the Arnoux-Yoccoz surface of genus 3 or greater.

        This presentation of the surface follows Section 2.3 of
        Joshua P. Bowman's paper "The Complete Family of Arnoux-Yoccoz
        Surfaces."

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.arnoux_yoccoz(4)
            sage: s
            Translation Surface in H_4(3^2) built from 16 triangles
            sage: TestSuite(s).run()
            sage: s.is_delaunay_decomposed()
            True
            sage: s = s.canonicalize()
            sage: s
            Translation Surface in H_4(3^2) built from 16 triangles
            sage: field=s.base_ring()
            sage: a = field.gen()
            sage: from sage.matrix.constructor import Matrix
            sage: m = Matrix([[a,0],[0,~a]])
            sage: ss = m*s
            sage: ss = ss.canonicalize()
            sage: s.cmp(ss) == 0
            True

        The Arnoux-Yoccoz pseudo-Anosov are known to have (minimal) invariant
        foliations with SAF=0::

            sage: S3 = translation_surfaces.arnoux_yoccoz(3)
            sage: Jxx, Jyy, Jxy = S3.j_invariant()
            sage: Jxx.is_zero() and Jyy.is_zero()
            True
            sage: Jxy
            [ 0  2  0]
            [ 2 -2  0]
            [ 0  0  2]

            sage: S4 = translation_surfaces.arnoux_yoccoz(4)
            sage: Jxx, Jyy, Jxy = S4.j_invariant()
            sage: Jxx.is_zero() and Jyy.is_zero()
            True
            sage: Jxy
            [ 0  2  0  0]
            [ 2 -2  0  0]
            [ 0  0  2  2]
            [ 0  0  2  0]
        """
        g = ZZ(genus)
        if g < 3:
            raise ValueError
        x = polygen(AA)
        p = sum([x**i for i in range(1, g + 1)]) - 1
        cp = AA.common_polynomial(p)
        alpha_AA = AA.polynomial_root(cp, RIF(1 / 2, 1))
        field = NumberField(alpha_AA.minpoly(), "alpha", embedding=alpha_AA)
        a = field.gen()
        V = VectorSpace(field, 2)
        p = [None for i in range(g + 1)]
        q = [None for i in range(g + 1)]
        p[0] = V(((1 - a**g) / 2, a**2 / (1 - a)))
        q[0] = V((-(a**g) / 2, a))
        p[1] = V((-(a ** (g - 1) + a**g) / 2, (a - a**2 + a**3) / (1 - a)))
        p[g] = V((1 + (a - a**g) / 2, (3 * a - 1 - a**2) / (1 - a)))
        for i in range(2, g):
            p[i] = V(((a - a**i) / (1 - a), a / (1 - a)))
        for i in range(1, g + 1):
            q[i] = V(
                (
                    (2 * a - a**i - a ** (i + 1)) / (2 * (1 - a)),
                    (a - a ** (g - i + 2)) / (1 - a),
                )
            )
        s = MutableOrientedSimilaritySurface(field)
        T = [None] * (2 * g + 1)
        Tp = [None] * (2 * g + 1)
        from sage.matrix.constructor import Matrix

        m = Matrix([[1, 0], [0, -1]])
        for i in range(1, g + 1):
            # T_i is (P_0,Q_i,Q_{i-1})
            T[i] = s.add_polygon(
                Polygon(
                    base_ring=field,
                    edges=[q[i] - p[0], q[i - 1] - q[i], p[0] - q[i - 1]],
                )
            )
            # T_{g+i} is (P_i,Q_{i-1},Q_{i})
            T[g + i] = s.add_polygon(
                Polygon(
                    base_ring=field,
                    edges=[q[i - 1] - p[i], q[i] - q[i - 1], p[i] - q[i]],
                )
            )
            # T'_i is (P'_0,Q'_{i-1},Q'_i)
            Tp[i] = s.add_polygon(m * s.polygon(T[i]))
            # T'_{g+i} is (P'_i,Q'_i, Q'_{i-1})
            Tp[g + i] = s.add_polygon(m * s.polygon(T[g + i]))
        for i in range(1, g):
            s.glue((T[i], 0), (T[i + 1], 2))
            s.glue((Tp[i], 2), (Tp[i + 1], 0))
        for i in range(1, g + 1):
            s.glue((T[i], 1), (T[g + i], 1))
            s.glue((Tp[i], 1), (Tp[g + i], 1))
        # P 0 Q 0 is paired with P' 0 Q' 0, ...
        s.glue((T[1], 2), (Tp[g], 2))
        s.glue((Tp[1], 0), (T[g], 0))
        # P1Q1 is paired with P'_g Q_{g-1}
        s.glue((T[g + 1], 2), (Tp[2 * g], 2))
        s.glue((Tp[g + 1], 0), (T[2 * g], 0))
        # P1Q0 is paired with P_{g-1} Q_{g-1}
        s.glue((T[g + 1], 0), (T[2 * g - 1], 2))
        s.glue((Tp[g + 1], 2), (Tp[2 * g - 1], 0))
        # PgQg is paired with Q1P2
        s.glue((T[2 * g], 2), (T[g + 2], 0))
        s.glue((Tp[2 * g], 0), (Tp[g + 2], 2))
        for i in range(2, g - 1):
            # PiQi is paired with Q'_i P'_{i+1}
            s.glue((T[g + i], 2), (Tp[g + i + 1], 2))
            s.glue((Tp[g + i], 0), (T[g + i + 1], 0))
        s.set_immutable()
        return s

    @staticmethod
    def from_flipper(h):
        r"""
        Build a (half-)translation surface from a flipper pseudo-Anosov.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: import flipper                             # optional - flipper

        A torus example::

            sage: t1 = (0r,1r,2r)                            # optional - flipper
            sage: t2 = (~0r,~1r,~2r)                         # optional - flipper
            sage: T = flipper.create_triangulation([t1,t2])  # optional - flipper
            sage: L1 = T.lamination([1r,0r,1r])              # optional - flipper
            sage: L2 = T.lamination([0r,1r,1r])              # optional - flipper
            sage: h1 = L1.encode_twist()                     # optional - flipper
            sage: h2 = L2.encode_twist()                     # optional - flipper
            sage: h = h1*h2^(-1r)                            # optional - flipper
            sage: f = h.flat_structure()                     # optional - flipper
            sage: ts = translation_surfaces.from_flipper(h)  # optional - flipper
            sage: ts                                         # optional - flipper; computation of the stratum fails here, see #227
            Half-Translation Surface built from 2 triangles
            sage: TestSuite(ts).run()                        # optional - flipper
            sage: from flatsurf.geometry.categories import HalfTranslationSurfaces  # optional: flipper
            sage: ts in HalfTranslationSurfaces()  # optional: flipper
            True

        A non-orientable example::

            sage: T = flipper.load('SB_4')                   # optional - flipper
            sage: h = T.mapping_class('s_0S_1s_2S_3s_1S_2')  # optional - flipper
            sage: h.is_pseudo_anosov()                       # optional - flipper
            True
            sage: S = translation_surfaces.from_flipper(h)   # optional - flipper
            sage: TestSuite(S).run()                         # optional - flipper
            sage: len(S.polygons())                          # optional - flipper
            4
            sage: from flatsurf.geometry.similarity_surface_generators import flipper_nf_element_to_sage
            sage: a = flipper_nf_element_to_sage(h.dilatation())  # optional - flipper

        """
        f = h.flat_structure()

        x = next(iter(f.edge_vectors.values())).x
        K = flipper_nf_to_sage(x.field)
        V = VectorSpace(K, 2)
        edge_vectors = {
            i: V(
                (flipper_nf_element_to_sage(e.x, K), flipper_nf_element_to_sage(e.y, K))
            )
            for i, e in f.edge_vectors.items()
        }

        to_polygon_number = {
            k: (i, j) for i, t in enumerate(f.triangulation) for j, k in enumerate(t)
        }

        from flatsurf import MutableOrientedSimilaritySurface

        S = MutableOrientedSimilaritySurface(K)

        for i, t in enumerate(f.triangulation):
            try:
                poly = Polygon(base_ring=K, edges=[edge_vectors[i] for i in tuple(t)])
            except ValueError:
                raise ValueError(
                    "t = {}, edges = {}".format(
                        t, [edge_vectors[i].n(digits=6) for i in t]
                    )
                )

            S.add_polygon(poly)

        for i, t in enumerate(f.triangulation):
            for j, k in enumerate(t):
                S.glue((i, j), to_polygon_number[~k])

        S.set_immutable()

        return S

    @staticmethod
    def origami(r, u, rr=None, uu=None, domain=None):
        r"""
        Return the origami defined by the permutations ``r`` and ``u``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: o
            Origami defined by r=(1,2) and u=(1,3)
            sage: o.stratum()
            H_2(2)
            sage: TestSuite(o).run()

        """
        return Origami(r, u, rr, uu, domain)

    @staticmethod
    def infinite_staircase():
        r"""
        Return the infinite staircase built as an origami.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: S = translation_surfaces.infinite_staircase()
            sage: S
            The infinite staircase
            sage: TestSuite(S).run()
        """
        return TranslationSurfaceGenerators._InfiniteStaircase()

    class _InfiniteStaircase(Origami):
        def __init__(self):
            super().__init__(
                self._vertical,
                self._horizontal,
                self._vertical,
                self._horizontal,
                domain=ZZ,
                root=ZZ(0),
            )

        def is_compact(self):
            return False

        def _vertical(self, x):
            if x % 2:
                return x + 1
            return x - 1

        def _horizontal(self, x):
            if x % 2:
                return x - 1
            return x + 1

        def _position_function(self, n):
            from flatsurf.geometry.similarity import SimilarityGroup

            SG = SimilarityGroup(QQ)
            if n % 2 == 0:
                return SG((n // 2, n // 2))
            else:
                return SG((n // 2, n // 2 + 1))

        def __repr__(self):
            return "The infinite staircase"

        def __hash__(self):
            return 1337

        def __eq__(self, other):
            r"""
            Return whether this surface is indistinguishable from ``other``.

            See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
            on this notion of equality.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S == S
                True

            """
            return isinstance(other, TranslationSurfaceGenerators._InfiniteStaircase)

        def graphical_surface(self, *args, **kwargs):
            default_position_function = kwargs.pop(
                "default_position_function", self._position_function
            )
            graphical_surface = super().graphical_surface(
                *args, default_position_function=default_position_function, **kwargs
            )
            graphical_surface.make_all_visible(limit=10)
            return graphical_surface

        def is_triangulated(self, limit=None):
            r"""
            Return whether this surface is triangulated, which it is not.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.is_triangulated()
                False

            """
            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

            return False

        def is_delaunay_triangulated(self, limit=None):
            r"""
            Return whether this surface is Delaunay triangulated, which it is not.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.is_delaunay_triangulated()
                False

            """
            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for is_delaunay_triangulated() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

            return False

        def is_delaunay_decomposed(self, limit=None):
            r"""
            Return whether this surface is made from Delaunay cells.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()
                sage: S.is_delaunay_decomposed()
                True

            """
            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for is_delaunay_decomposed() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

            return True

    @staticmethod
    def t_fractal(w=ZZ_1, r=ZZ_2, h1=ZZ_1, h2=ZZ_1):
        r"""
        Return the T-fractal with parameters ``w``, ``r``, ``h1``, ``h2``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: tf = translation_surfaces.t_fractal()  # long time (.4s)
            sage: tf  # long time (see above)
            The T-fractal surface with parameters w=1, r=2, h1=1, h2=1
            sage: TestSuite(tf).run()  # long time (see above)
        """
        return tfractal_surface(w, r, h1, h2)

    @staticmethod
    def e_infinity_surface(lambda_squared=None, field=None):
        r"""
        The translation surface based on the `E_\infinity` graph.

        The biparite graph is shown below, with edges numbered::

              0   1   2  -2   3  -3   4  -4
            *---o---*---o---*---o---*---o---*...
                    |
                    |-1
                    o

        Here, black vertices are colored ``*``, and white ``o``.
        Black nodes represent vertical cylinders and white nodes
        represent horizontal cylinders.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.e_infinity_surface()
            sage: TestSuite(s).run()  # long time (1s)
        """
        return EInfinitySurface(lambda_squared, field)

    @staticmethod
    def chamanara(alpha):
        r"""
        Return the Chamanara surface with parameter ``alpha``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C
            Minimal Translation Cover of Chamanara surface with parameter 1/2

        TESTS::

            sage: TestSuite(C).run()
            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: C in TranslationSurfaces()
            True

        """
        from .chamanara import chamanara_surface

        return chamanara_surface(alpha)


translation_surfaces = TranslationSurfaceGenerators()
