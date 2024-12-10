r"""
The category of half-translation surfaces.

A half-translation surface is a surface built by gluing Euclidean polygons. The
sides of the polygons can be glued by translations or half-translations
(rotations by angle π).

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES:

We glue all the sides of a square to themselves. Since each gluing is just a
rotation by π, this is a half-translation surface::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = similarity_surfaces.self_glued_polygon(P)
    sage: S.set_immutable()

    sage: C = S.category()

    sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
    sage: C.is_subcategory(HalfTranslationSurfaces())
    True

"""

# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2013-2019 Vincent Delecroix
#                      2013-2019 W. Patrick Hooper
#                      2023-2024 Julian Rüth
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import (
    SurfaceCategory,
    SurfaceCategoryWithAxiom,
)
from sage.misc.lazy_import import LazyImport
from sage.all import QQ, AA


class HalfTranslationSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations and half-translations (rotations by angle π).

    EXAMPLES::

        sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
        sage: HalfTranslationSurfaces()
        Category of half translation surfaces

    """

    def super_categories(self):
        r"""
        Return the categories that a half-translation surface is always a
        member of.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
            sage: HalfTranslationSurfaces().super_categories()
            [Category of dilation surfaces, Category of rational cone surfaces]

        """
        from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
        from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces

        return [DilationSurfaces(), ConeSurfaces().Rational()]

    # Declare that the "positive" half-translation surfaces are called
    # "translation surfaces".
    Positive = LazyImport(
        "flatsurf.geometry.categories.translation_surfaces", "TranslationSurfaces"
    )

    class ParentMethods:
        r"""
        Provides methods available to all half-translation surfaces.

        If you want to add functionality for such surfaces you most likely want
        to put it here.
        """

        def is_translation_surface(self, positive=True):
            r"""
            Return whether this surface is a (half-)translation surface.

            This overrides
            :meth:`.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_translation_surface`.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                sage: H = B.minimal_cover(cover_type="half-translation")

                sage: H.is_translation_surface(positive=False)
                True
                sage: H.is_translation_surface(positive=True)
                False

            """
            if not positive:
                return True

            # If this is not explicitly a translation surface, we have to
            # decide with the generic checks whether it is a positive
            # half-translation surface.
            return super(  # pylint: disable=bad-super-call
                HalfTranslationSurfaces().parent_class, self
            ).is_translation_surface(positive=positive)

    class Orientable(SurfaceCategoryWithAxiom):
        r"""
        The category of orientable half-translation surfaces.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
            sage: H = B.minimal_cover(cover_type="half-translation")

            sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
            sage: H in HalfTranslationSurfaces().Orientable()
            True

        """

        class WithoutBoundary(SurfaceCategoryWithAxiom):
            r"""
            The category of orientable half-translation surfaces without boundary.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                sage: H = B.minimal_cover(cover_type="half-translation")

                sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
                sage: H in HalfTranslationSurfaces().Orientable().WithoutBoundary()
                True

            """

            class ParentMethods:
                r"""
                Provides methods available to all orientable half-translation
                surfaces.

                If you want to add functionality for such surfaces you most likely
                want to put it here.
                """

                def stratum(self):
                    r"""
                    EXAMPLES::

                        sage: from flatsurf import polygons, similarity_surfaces
                        sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                        sage: H = B.minimal_cover(cover_type="half-translation")
                        sage: H.stratum()
                        Q_1(3, -1^3)

                    TESTS:

                    Verify that the stratum is correct for surfaces with self-glued edges::

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                        sage: S = similarity_surfaces.self_glued_polygon(P)
                        sage: S.stratum()
                        Q_0(0, -1^4)

                    """
                    angles = self.angles()

                    for a, b in self.gluings():
                        if a == b:
                            angles.append(QQ(1 / 2))

                    if all(x.denominator() == 1 for x in angles):
                        raise NotImplementedError

                    from sage.all import ZZ
                    from surface_dynamics import Stratum

                    return Stratum(
                        sorted([ZZ(2 * a - 2) for a in angles], reverse=True), 2
                    )

    class Oriented(SurfaceCategoryWithAxiom):
        r"""
        The category of oriented half-translation surfaces, i.e., orientable
        half-translation surfaces which can be oriented in a way compatible
        with the embedding of their polygons in the real plane.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
            sage: H = B.minimal_cover(cover_type="half-translation")

            sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
            sage: H in HalfTranslationSurfaces().Oriented()
            True

        """

        class ParentMethods:
            r"""
            Provides methods available to all oriented half-translation
            surfaces.

            If you want to add functionality for such surfaces you most likely
            want to put it here.
            """

            def holonomy_field(self):
                r"""
                Return the relative holonomy field of this translation or half-translation surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces, polygons, similarity_surfaces

                    sage: S = translation_surfaces.veech_2n_gon(5)
                    sage: S.holonomy_field()
                    Number Field in a0 with defining polynomial x^2 - x - 1 with a0 = ...
                    sage: S.base_ring()
                    Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.175570504584947?

                    sage: T = translation_surfaces.torus((1, AA(2).sqrt()), (AA(3).sqrt(), 3))
                    sage: T.holonomy_field()
                    Rational Field

                    sage: T = polygons.triangle(1,6,11)
                    sage: S = similarity_surfaces.billiard(T)
                    sage: S = S.minimal_cover("translation")
                    sage: S.base_ring()
                    Number Field in c with defining polynomial x^6 - 6*x^4 + 9*x^2 - 3 with c = 1.969615506024417?
                    sage: S.holonomy_field()
                    Number Field in c0 with defining polynomial x^3 - 3*x - 1 with c0 = 1.879385241571817?
                """
                return self.normalized_coordinates()[0].base_ring()

            def _test_half_translation_surface(self, **options):
                r"""
                Verify that this is a half-translation surface.

                EXAMPLES::

                    sage: from flatsurf import Polygon, similarity_surfaces
                    sage: P = Polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S._test_half_translation_surface()

                """
                tester = self._tester(**options)

                limit = None

                if not self.is_finite_type():
                    from flatsurf import MutableOrientedSimilaritySurface

                    self = MutableOrientedSimilaritySurface.from_surface(
                        self, labels=self.labels()[:32]
                    )

                from flatsurf.geometry.categories import TranslationSurfaces

                tester.assertTrue(
                    TranslationSurfaces.ParentMethods._is_translation_surface(
                        self, positive=False, limit=limit
                    )
                )

        class FiniteType(SurfaceCategoryWithAxiom):
            r"""
            The category of oriented half-translation surfaces built from
            finitely many polygons.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                sage: H = B.minimal_cover(cover_type="half-translation")

                sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
                sage: H in HalfTranslationSurfaces().Oriented().FiniteType()
                True

            """

            class ParentMethods:
                r"""
                Provides methods available to all oriented half-translation
                surfaces built from finitely many polygons.

                If you want to add functionality for such surfaces you most
                likely want to put it here.
                """

                def normalized_coordinates(self):
                    r"""
                    Return a pair ``(new_surface, matrix)`` where ``new_surface`` is defined over the
                    holonomy field and ``matrix`` is the transition matrix that maps this surface to
                    ``new_surface``.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces, polygons, similarity_surfaces

                        sage: S = translation_surfaces.veech_2n_gon(5)
                        sage: U, mat = S.normalized_coordinates()
                        sage: U.base_ring()
                        Number Field in a0 with defining polynomial x^2 - x - 1 with a0 = ...
                        sage: mat
                        [             0 -2/5*a^3 + 2*a]
                        [            -1 -3/5*a^3 + 2*a]

                        sage: T = translation_surfaces.torus((1, AA(2).sqrt()), (AA(3).sqrt(), 3))
                        sage: U, mat = T.normalized_coordinates()
                        sage: U.base_ring()
                        Rational Field
                        sage: U.holonomy_field()
                        Rational Field
                        sage: mat
                        [-2.568914100752347?  1.816496580927726?]
                        [-5.449489742783178?  3.146264369941973?]
                        sage: TestSuite(U).run()

                        sage: T = polygons.triangle(1,6,11)
                        sage: S = similarity_surfaces.billiard(T)
                        sage: S = S.minimal_cover("translation")
                        sage: U, _ = S.normalized_coordinates()
                        sage: U.base_ring()
                        Number Field in c0 with defining polynomial x^3 - 3*x - 1 with c0 = 1.879385241571817?
                        sage: U.holonomy_field() == U.base_ring()
                        True
                        sage: S.base_ring()
                        Number Field in c with defining polynomial x^6 - 6*x^4 + 9*x^2 - 3 with c = 1.969615506024417?
                        sage: TestSuite(U).run()

                        sage: from flatsurf import EuclideanPolygonsWithAngles
                        sage: polygons = EuclideanPolygonsWithAngles((1, 3, 1, 1))
                        sage: p = polygons.an_element()
                        sage: B = similarity_surfaces.billiard(p)
                        sage: B.minimal_cover("translation")
                        Minimal Translation Cover of Genus 0 Rational Cone Surface built from 2 equilateral triangles
                        sage: S = B.minimal_cover("translation")
                        sage: S, _ = S.normalized_coordinates()
                        sage: S
                        Translation Surface in H_1(0^6) built from 6 right triangles

                    TESTS:

                    Verify that #89 has been resolved::

                        sage: from pyexactreal import ExactReals  # optional: pyexactreal  # random output due to pkg_resources deprecation warnings
                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.square_torus()
                        sage: S = S.change_ring(ExactReals())  # optional: pyexactreal
                        sage: S.normalized_coordinates()  # optional: pyexactreal
                        Traceback (most recent call last):
                        ...
                        NotImplementedError: base ring must be a field to normalize coordinates of the surface


                    """
                    from sage.all import matrix

                    if self.base_ring() is QQ:
                        return (self, matrix(QQ, 2, 2, 1))

                    from sage.categories.all import Fields

                    if self.base_ring() not in Fields():
                        raise NotImplementedError(
                            "base ring must be a field to normalize coordinates of the surface"
                        )

                    lab = next(iter(self.labels()))
                    p = self.polygon(lab)
                    u = p.edge(1)
                    v = -p.edge(0)
                    i = 1
                    from flatsurf.geometry.euclidean import ccw

                    while ccw(u, v) == 0:
                        i += 1
                        u = p.edge(i)
                        v = -p.edge(i - 1)
                    M = matrix(2, [u, v]).transpose().inverse()
                    assert M.det() > 0
                    hols = []
                    for lab in self.labels():
                        p = self.polygon(lab)
                        for e in range(len(p.vertices())):
                            w = M * p.edge(e)
                            hols.append(w[0])
                            hols.append(w[1])
                    if self.base_ring() is AA:
                        from flatsurf.geometry.subfield import (
                            number_field_elements_from_algebraics,
                        )

                        K, new_hols = number_field_elements_from_algebraics(hols)
                    else:
                        from flatsurf.geometry.subfield import subfield_from_elements

                        K, new_hols, _ = subfield_from_elements(self.base_ring(), hols)

                    from flatsurf.geometry.polygon import Polygon
                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    S = MutableOrientedSimilaritySurface(K)
                    relabelling = {}
                    k = 0
                    for lab in self.labels():
                        m = len(self.polygon(lab).vertices())
                        relabelling[lab] = S.add_polygon(
                            Polygon(
                                edges=[
                                    (new_hols[k + 2 * i], new_hols[k + 2 * i + 1])
                                    for i in range(m)
                                ],
                                base_ring=K,
                            )
                        )
                        k += 2 * m

                    for (p1, e1), (p2, e2) in self.gluings():
                        S.glue((relabelling[p1], e1), (relabelling[p2], e2))

                    S._refine_category_(self.category())
                    return S, M
