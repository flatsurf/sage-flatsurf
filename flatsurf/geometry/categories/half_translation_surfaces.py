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

    sage: from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
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

from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.misc.lazy_import import LazyImport
from sage.all import QQ


class HalfTranslationSurfaces(Category):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations and half-translations (translations followed by rotations
    among an angle π.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
        sage: HalfTranslationSurfaces()
        Category of half translation surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
        from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces
        return [DilationSurfaces(), ConeSurfaces().Rational()]

    Positive = LazyImport('flatsurf.geometry.categories.translation_surfaces', 'TranslationSurfaces')

    class Orientable(CategoryWithAxiom):
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

    class Oriented(CategoryWithAxiom):
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
