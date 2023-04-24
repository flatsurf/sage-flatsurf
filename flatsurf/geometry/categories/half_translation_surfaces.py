r"""
The category of half-translation surfaces.

A half-translation surface is a surface built by gluing Euclidean polygons. The
sides of the polygons can be glued with translations or half-translations
(translation followed by a rotation of angle π.)

EXAMPLES:

We glue all the sides of a square to themselves. Since each gluing is just a
rotation of π, this is a half-translation surface::

    sage: from flatsurf import polygons, similarity_surfaces
    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import SurfaceCategory, SurfaceCategoryWithAxiom
from sage.misc.lazy_import import LazyImport
from sage.all import QQ, AA


class HalfTranslationSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations and half-translations (translations followed by rotations
    among an angle π.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
        sage: HalfTranslationSurfaces()
        Category of half translation surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
        from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces
        return [DilationSurfaces(), ConeSurfaces().Rational()]

    Positive = LazyImport('flatsurf.geometry.categories.translation_surfaces', 'TranslationSurfaces')

    class Orientable(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def stratum(self):
                r"""
                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                    sage: H = B.minimal_cover(cover_type="half-translation")
                    sage: H.stratum()
                    Q_1(3, -1^3)
                """
                angles = self.angles()
                if all(x.denominator() == 1 for x in angles):
                    raise NotImplementedError
                from surface_dynamics import QuadraticStratum

                return QuadraticStratum(*[2 * a - 2 for a in angles])

    class Oriented(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def angles(self, numerical=False, return_adjacent_edges=False):
                r"""
                Return the set of angles around the vertices of the surface.

                These are given as multiple of `2 \pi`.

                EXAMPLES::

                    sage: import flatsurf.geometry.similarity_surface_generators as sfg
                    sage: sfg.translation_surfaces.regular_octagon().angles()
                    [3]
                    sage: S = sfg.translation_surfaces.veech_2n_gon(5)
                    sage: S.angles()
                    [2, 2]
                    sage: S.angles(numerical=True)
                    [2.0, 2.0]
                    sage: S.angles(return_adjacent_edges=True) # random output
                    [(2, [(0, 1), (0, 5), (0, 9), (0, 3), (0, 7)]),
                     (2, [(0, 0), (0, 4), (0, 8), (0, 2), (0, 6)])]
                    sage: S.angles(numerical=True, return_adjacent_edges=True) # random output
                    [(2.0, [(0, 1), (0, 5), (0, 9), (0, 3), (0, 7)]),
                     (2.0, [(0, 0), (0, 4), (0, 8), (0, 2), (0, 6)])]

                    sage: sfg.translation_surfaces.veech_2n_gon(6).angles()
                    [5]
                    sage: sfg.translation_surfaces.veech_double_n_gon(5).angles()
                    [3]
                    sage: sfg.translation_surfaces.cathedral(1, 1).angles()
                    [3, 3, 3]

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: B = similarity_surfaces.billiard(polygons.triangle(1, 2, 5))
                    sage: H = B.minimal_cover(cover_type="half-translation")
                    sage: S = B.minimal_cover(cover_type="translation")
                    sage: H.angles()
                    [1/2, 5/2, 1/2, 1/2]
                    sage: S.angles()
                    [1, 5, 1, 1]

                    sage: H.angles(return_adjacent_edges=True)
                     [(1/2, [...]), (5/2, [...]), (1/2, [...]), (1/2, [...])]
                    sage: S.angles(return_adjacent_edges=True)
                     [(1, [...]), (5, [...]), (1, [...]), (1, [...])]
                """
                if not self.is_finite():
                    raise NotImplementedError("the set of edges is infinite!")

                edges = set(self.edge_iterator())
                angles = []

                if return_adjacent_edges:
                    while edges:
                        # Note that iteration order here is different for different
                        # versions of Python. Therefore, the output in the doctest
                        # above is random.
                        pair = p, e = next(iter(edges))
                        ve = self.polygon(p).edge(e)
                        angle = 0
                        adjacent_edges = []
                        while pair in edges:
                            adjacent_edges.append(pair)
                            edges.remove(pair)
                            f = (e - 1) % self.polygon(p).num_edges()
                            ve = self.polygon(p).edge(e)
                            vf = -self.polygon(p).edge(f)
                            ppair = pp, ff = self.opposite_edge(p, f)
                            angle += (
                                (ve[0] > 0 and vf[0] <= 0)
                                or (ve[0] < 0 and vf[0] >= 0)
                                or (ve[0] == vf[0] == 0)
                            )
                            pair, p, e = ppair, pp, ff
                        if numerical:
                            angles.append((float(angle) / 2, adjacent_edges))
                        else:
                            angles.append((QQ((angle, 2)), adjacent_edges))
                else:
                    while edges:
                        pair = p, e = next(iter(edges))
                        angle = 0
                        while pair in edges:
                            edges.remove(pair)
                            f = (e - 1) % self.polygon(p).num_edges()
                            ve = self.polygon(p).edge(e)
                            vf = -self.polygon(p).edge(f)
                            ppair = pp, ff = self.opposite_edge(p, f)
                            angle += (
                                (ve[0] > 0 and vf[0] <= 0)
                                or (ve[0] < 0 and vf[0] >= 0)
                                or (ve[0] == vf[0] == 0)
                            )
                            pair, p, e = ppair, pp, ff
                        if numerical:
                            angles.append(float(angle) / 2)
                        else:
                            angles.append(QQ((angle, 2)))

                return angles

            def holonomy_field(self):
                r"""
                Return the relative holonomy field of this translation or half-translation surface.

                EXAMPLES::

                    sage: from flatsurf import *

                    sage: S = translation_surfaces.veech_2n_gon(5)
                    sage: S.holonomy_field()
                    Number Field in a0 with defining polynomial x^2 - x - 1 with a0 = 1.618033988749895?
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

        class FiniteType(SurfaceCategoryWithAxiom):
            class ParentMethods:
                def normalized_coordinates(self):
                    r"""
                    Return a pair ``(new_surface, matrix)`` where ``new_surface`` is defined over the
                    holonomy field and ``matrix`` is the transition matrix that maps this surface to
                    ``new_surface``.

                    EXAMPLES::

                        sage: from flatsurf import *

                        sage: S = translation_surfaces.veech_2n_gon(5)
                        sage: U, mat = S.normalized_coordinates()
                        sage: U.base_ring()
                        Number Field in a0 with defining polynomial x^2 - x - 1 with a0 = 1.618033988749895?
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

                        sage: from flatsurf import EquiangularPolygons
                        sage: E = EquiangularPolygons(1, 3, 1, 1)
                        sage: r1, r2 = [r.vector() for r in E.lengths_polytope().rays()]
                        sage: p = E(r1 + r2)
                        sage: B = similarity_surfaces.billiard(p)
                        sage: B.minimal_cover("translation")
                        MinimalTranslationCover(Surface built from 2 polygons)
                        sage: S = B.minimal_cover("translation")
                        sage: S, _ = S.normalized_coordinates()
                        sage: S
                        Surface built from 6 polygons
                    """
                    from sage.all import matrix
                    if self.base_ring() is QQ:
                        return (self, matrix(QQ, 2, 2, 1))

                    lab = next(self.label_iterator())
                    p = self.polygon(lab)
                    u = p.edge(1)
                    v = -p.edge(0)
                    i = 1
                    from flatsurf.geometry.polygon import wedge_product
                    while wedge_product(u, v) == 0:
                        i += 1
                        u = p.edge(i)
                        v = -p.edge(i - 1)
                    M = matrix(2, [u, v]).transpose().inverse()
                    assert M.det() > 0
                    hols = []
                    for lab in self.label_iterator():
                        p = self.polygon(lab)
                        for e in range(p.num_edges()):
                            w = M * p.edge(e)
                            hols.append(w[0])
                            hols.append(w[1])
                    if self.base_ring() is AA:
                        from flatsurf.geometry.subfield import number_field_elements_from_algebraics

                        K, new_hols = number_field_elements_from_algebraics(hols)
                    else:
                        from flatsurf.geometry.subfield import subfield_from_elements

                        K, new_hols, _ = subfield_from_elements(self.base_ring(), hols)

                    from flatsurf.geometry.polygon import ConvexPolygons
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                    S = MutableOrientedSimilaritySurface(K)
                    C = ConvexPolygons(K)
                    relabelling = {}
                    k = 0
                    for lab in self.label_iterator():
                        m = self.polygon(lab).num_edges()
                        relabelling[lab] = S.add_polygon(
                            C(
                                edges=[
                                    (new_hols[k + 2 * i], new_hols[k + 2 * i + 1]) for i in range(m)
                                ]
                            )
                        )
                        k += 2 * m

                    for (p1, e1), (p2, e2) in self.edge_iterator(gluings=True):
                        S.set_edge_pairing(relabelling[p1], e1, relabelling[p2], e2)

                    S._refine_category_(self.category())
                    return S, M