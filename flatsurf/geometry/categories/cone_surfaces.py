r"""
The category of cone surfaces.

A cone surface is a surface that can be built by gluing Euclidean polygons
along their edges such that the matrix describing monodromy along a closed path
is an isometry; that matrix is given by multiplying the individual matrices
that describe how to transition between pairs of glued edges, see
:meth:`edge_matrix`.

In sage-flatsurf, we restrict cone surfaces slightly by requiring that a cone
surface is given by polygons such that each edge matrix is an isometry.

EXAMPLES:

We glue the sides of a square with a rotation of π/2. Since each gluing is just
a rotation, this is a cone surface::

    sage: from flatsurf import polygon, MutableOrientedSimilaritySurface
    sage: P = polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = MutableOrientedSimilaritySurface(QQ)
    sage: S.add_polygon(P, label=0)
    0
    sage: S.glue((0, 0), (0, 1))
    sage: S.glue((0, 2), (0, 3))
    sage: S.set_immutable()

    sage: C = S.category()

    sage: from flatsurf.geometry.categories import ConeSurfaces
    sage: C.is_subcategory(ConeSurfaces())
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


class ConeSurfaces(Category):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    isometries on the edges.

    See :mod:`flatsurf.geometry.categories.cone_surfaces` and
    :meth:`is_cone_surface` on how this differs slightly from the customary
    definition of a cone surface.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import ConeSurfaces
        sage: ConeSurfaces()
        Category of cone surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.similarity_surfaces import SimilaritySurfaces
        return [SimilaritySurfaces()]

    class ParentMethods:
        @staticmethod
        def _is_cone_surface(surface, limit=None):
            r"""
            Return whether ``surface`` is a cone surface by checking how its
            polygons are glued.

            INPUT:

            - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
              the first ``limit`` polygons are checked

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.infinite_staircase()

                sage: from flatsurf.geometry.categories import ConeSurfaces
                sage: ConeSurfaces.ParentMethods._is_cone_surface(S, limit=8)
                True

            ::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons((2, 0),(-1, 3),(-1, -3))
                sage: S = similarity_surfaces.self_glued_polygon(P)

                sage: ConeSurfaces.ParentMethods._is_cone_surface(S)
                True

            """
            if 'Oriented' not in surface.category().axioms():
                raise NotImplementedError

            labels = surface.labels()

            if limit is not None:
                from itertools import islice
                labels = islice(labels, limit)

            for label in labels:
                for edge in range(surface.polygon(label).num_edges()):
                    cross = surface.opposite_edge(label, edge)

                    if cross is None:
                        continue

                    # We do not call self.edge_matrix() since the surface might
                    # have overriden this (just returning the identity matrix e.g.)
                    # and we want to deduce the matrix from the attached polygon
                    # edges instead.
                    from flatsurf.geometry.categories import SimilaritySurfaces
                    matrix = SimilaritySurfaces.Oriented.ParentMethods.edge_matrix.f(surface, label, edge)

                    if matrix * matrix.transpose() != 1:
                        return False

            return True

    class FiniteType(CategoryWithAxiom):
        class ParentMethods:
            def area(self):
                r"""
                Return the area of this surface.
                """
                return sum(p.area() for p in self.polygons())

    class Oriented(CategoryWithAxiom):
        class ParentMethods:
            def _test_cone_surface(self, **options):
                r"""
                Verify that this is a cone surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S.set_immutable()
                    sage: S._test_cone_surface()

                """
                tester = self._tester(**options)

                limit = None

                if not self.is_finite():
                    limit = 32

                tester.assertTrue(ConeSurfaces.ParentMethods._is_cone_surface(self, limit=limit))

        class FiniteType(CategoryWithAxiom):
            class ParentMethods:
                def angles(self, numerical=False, return_adjacent_edges=False):
                    r"""
                    Return the set of angles around the vertices of the surface.

                    EXAMPLES::

                        sage: from flatsurf import polygons, similarity_surfaces
                        sage: T = polygons.triangle(3, 4, 5)
                        sage: S = similarity_surfaces.billiard(T)
                        sage: S.angles()
                        [1/3, 1/4, 5/12]
                        sage: S.angles(numerical=True)   # abs tol 1e-14
                        [0.333333333333333, 0.250000000000000, 0.416666666666667]

                        sage: S.angles(return_adjacent_edges=True)
                        [(1/3, [(0, 1), (1, 2)]), (1/4, [(0, 0), (1, 0)]), (5/12, [(1, 1), (0, 2)])]
                    """
                    edges = [pair for pair in self.edges()]
                    edges = set(edges)
                    angles = []

                    if return_adjacent_edges:
                        while edges:
                            p, e = edges.pop()
                            adjacent_edges = [(p, e)]
                            angle = self.polygon(p).angle(e, numerical=numerical)
                            pp, ee = self.opposite_edge(p, (e - 1) % self.polygon(p).num_edges())
                            while pp != p or ee != e:
                                edges.remove((pp, ee))
                                adjacent_edges.append((pp, ee))
                                angle += self.polygon(pp).angle(ee, numerical=numerical)
                                pp, ee = self.opposite_edge(
                                    pp, (ee - 1) % self.polygon(pp).num_edges()
                                )
                            angles.append((angle, adjacent_edges))
                    else:
                        while edges:
                            p, e = edges.pop()
                            angle = self.polygon(p).angle(e, numerical=numerical)
                            pp, ee = self.opposite_edge(p, (e - 1) % self.polygon(p).num_edges())
                            while pp != p or ee != e:
                                edges.remove((pp, ee))
                                angle += self.polygon(pp).angle(ee, numerical=numerical)
                                pp, ee = self.opposite_edge(
                                    pp, (ee - 1) % self.polygon(pp).num_edges()
                                )
                            angles.append(angle)

                    return angles
