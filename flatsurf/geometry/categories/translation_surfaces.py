r"""
The category of translation surfaces.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from Euclidean polygons whose glued edges can be transformed into each
other with translations.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.square_torus()
    sage: C = S.category()

    sage: from flatsurf.geometry.categories import TranslationSurfaces
    sage: C.is_subcategory(TranslationSurfaces())
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

from flatsurf.geometry.categories.surface_category import SurfaceCategoryWithAxiom
from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces


class TranslationSurfaces(SurfaceCategoryWithAxiom):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import TranslationSurfaces
        sage: TranslationSurfaces()
        Category of translation surfaces

    """
    _base_category_class_and_axiom = (HalfTranslationSurfaces, 'Positive')

    # TODO: They satisfy the axiom of no self-gluing and are oriented.

    def extra_super_categories(self):
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces
        return (PolygonalSurfaces().Oriented(),)

    class ParentMethods:
        def minimal_translation_cover(self):
            return self

        def edge_matrix(self, p, e=None):
            if e is None:
                p, e = p
            if e < 0 or e >= self.polygon(p).num_edges():
                raise ValueError
            from sage.all import identity_matrix
            return identity_matrix(self.base_ring(), 2)

        def standardize_polygons(self, in_place=False):
            r"""
            Replaces each polygon with a new polygon which differs by
            translation and reindexing. The new polygon will have the property
            that vertex zero is the origin, and all vertices lie either in the
            upper half plane, or on the x-axis with non-negative x-coordinate.

            This is done to the current surface if in_place=True. A mutable
            copy is created and returned if in_place=False (as default).

            EXAMPLES::

                sage: from flatsurf import *
                sage: s=translation_surfaces.veech_double_n_gon(4)
                sage: s.polygon(1)
                Polygon: (0, 0), (-1, 0), (-1, -1), (0, -1)
                sage: [s.opposite_edge(0,i) for i in range(4)]
                [(1, 0), (1, 1), (1, 2), (1, 3)]
                sage: ss=s.standardize_polygons()
                sage: ss.polygon(1)
                Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
                sage: [ss.opposite_edge(0,i) for i in range(4)]
                [(1, 2), (1, 3), (1, 0), (1, 1)]
                sage: TestSuite(ss).run()

            Make sure first vertex is sent to origin::

                sage: from flatsurf import *
                sage: P = ConvexPolygons(QQ)
                sage: p = P(vertices = ([(1,1),(2,1),(2,2),(1,2)]))
                sage: s = Surface_list(QQ)
                sage: s.add_polygon(p)
                0
                sage: s.change_polygon_gluings(0, [(0,2),(0,3),(0,0),(0,1)])
                sage: s.change_base_label(0)
                sage: s.set_immutable()
                sage: s.standardize_polygons().polygon(0)
                Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            """
            if self.is_finite():
                if in_place:
                    if not self.is_mutable():
                        raise ValueError(
                            "An in_place call for standardize_polygons can only be done for a mutable surface."
                        )
                    s = self
                else:
                    s = self.copy(mutable=True)
                cv = {}  # dictionary for non-zero canonical vertices
                for label, polygon in s.label_iterator(polygons=True):
                    best = 0
                    best_pt = polygon.vertex(best)
                    for v in range(1, polygon.num_edges()):
                        pt = polygon.vertex(v)
                        if (pt[1] < best_pt[1]) or (
                            pt[1] == best_pt[1] and pt[0] < best_pt[0]
                        ):
                            best = v
                            best_pt = pt
                    # We replace the polygon if the best vertex is not the zero vertex, or
                    # if the coordinates of the best vertex differs from the origin.
                    if not (best == 0 and best_pt.is_zero()):
                        cv[label] = best
                for label, v in cv.items():
                    s.set_vertex_zero(label, v, in_place=True)
                return s
            else:
                if in_place:
                    raise NotImplementedError(
                        "in place standardization only available for finite surfaces"
                    )

                from flatsurf.geometry.similarity_surface import SimilaritySurface
                from flatsurf.geometry.translation_surface import LazyStandardizedPolygonSurface
                return SimilaritySurface(LazyStandardizedPolygonSurface(self))

        def cmp(self, s2, limit=None):
            r"""
            Compare two surfaces. This is an ordering returning -1, 0, or 1.

            The surfaces will be considered equal if and only if there is a translation automorphism
            respecting the polygons and the base_labels.

            If the two surfaces are infinite, we just examine the first limit polygons.
            """
            if self.is_finite():
                if s2.is_finite():
                    if limit is not None:
                        raise ValueError("limit only enabled for finite surfaces")

                    # print("comparing number of polygons")
                    sign = self.num_polygons() - s2.num_polygons()
                    if sign > 0:
                        return 1
                    if sign < 0:
                        return -1
                    # print("comparing polygons")
                    lw1 = self.walker()
                    lw2 = s2.walker()
                    for p1, p2 in zip(lw1.polygon_iterator(), lw2.polygon_iterator()):
                        # Uses Polygon.cmp:
                        ret = p1.cmp(p2)
                        if ret != 0:
                            return ret
                    # Polygons are identical. Compare edge gluings.
                    # print("comparing edge gluings")
                    for pair1, pair2 in zip(lw1.edge_iterator(), lw2.edge_iterator()):
                        l1, e1 = self.opposite_edge(pair1)
                        l2, e2 = s2.opposite_edge(pair2)
                        num1 = lw1.label_to_number(l1)
                        num2 = lw2.label_to_number(l2)
                        ret = (num1 > num2) - (num1 < num2)
                        if ret:
                            return ret
                        ret = (e1 > e2) - (e1 < e2)
                        if ret:
                            return ret
                    return 0
                else:
                    # s1 is finite but s2 is infinite.
                    return -1
            else:
                if s2.is_finite():
                    # s1 is infinite but s2 is finite.
                    return 1
                else:
                    # both surfaces are infinite.
                    lw1 = self.walker()
                    lw2 = s2.walker()
                    count = 0
                    for (l1, p1), (l2, p2) in zip(
                        lw1.label_polygon_iterator(), lw2.label_polygon_iterator()
                    ):
                        # Uses Polygon.cmp:
                        ret = p1.cmp(p2)
                        if ret != 0:
                            print("Polygons differ")
                            return ret
                        # If here the number of edges should be equal.
                        for e in range(p1.num_edges()):
                            ll1, ee1 = self.opposite_edge(l1, e)
                            ll2, ee2 = s2.opposite_edge(l2, e)
                            num1 = lw1.label_to_number(ll1, search=True, limit=limit)
                            num2 = lw2.label_to_number(ll2, search=True, limit=limit)
                            ret = (num1 > num2) - (num1 < num2)
                            if ret:
                                return ret
                            ret = (ee1 > ee2) - (ee1 < ee2)
                            if ret:
                                return ret
                        if count >= limit:
                            break
                        count += 1
                    return 0

        def stratum(self):
            r"""
            Return the stratum this surface belongs to.

            This uses the package ``surface-dynamics``
            (see http://www.labri.fr/perso/vdelecro/flatsurf_sage.html)

            EXAMPLES::

                sage: import flatsurf.geometry.similarity_surface_generators as sfg
                sage: sfg.translation_surfaces.octagon_and_squares().stratum()
                H_3(4)
            """
            from surface_dynamics import AbelianStratum
            from sage.rings.integer_ring import ZZ

            return AbelianStratum([ZZ(a - 1) for a in self.angles()])

        def canonicalize_mapping(self):
            r"""
            Return a SurfaceMapping canonicalizing this translation surface.
            """
            from flatsurf.geometry.mappings import (
                canonicalize_translation_surface_mapping,
            )

            return canonicalize_translation_surface_mapping(self)

        def canonicalize(self, in_place=False):
            r"""
            Return a canonical version of this translation surface.

            EXAMPLES:

            We will check if an element lies in the Veech group::

                sage: from flatsurf import *
                sage: s = translation_surfaces.octagon_and_squares()
                sage: s
                Surface built from 3 polygons
                sage: from flatsurf.geometry.categories import TranslationSurfaces
                sage: s in TranslationSurfaces()
                True
                sage: a = s.base_ring().gen()
                sage: mat = Matrix([[1,2+a],[0,1]])
                sage: s1 = s.canonicalize()
                sage: s1.set_immutable()
                sage: s2 = (mat*s).canonicalize()
                sage: s2.set_immutable()
                sage: s1.cmp(s2) == 0
                True
                sage: hash(s1) == hash(s2)
                True
            """
            if in_place:
                if not self.is_mutable():
                    raise ValueError(
                        "canonicalize with in_place=True is only defined for mutable translation surfaces."
                    )
                s = self
            else:
                s = self.copy(mutable=True)
            if not s.is_finite():
                raise ValueError(
                    "canonicalize is only defined for finite translation surfaces."
                )
            s.delaunay_decomposition(in_place=True)
            s.standardize_polygons(in_place=True)
            ss = s.copy(mutable=True)
            labels = {label for label in s.label_iterator()}
            labels.remove(s.base_label())
            for label in labels:
                ss.change_base_label(label)
                if ss.cmp(s) > 0:
                    s.change_base_label(label)
            # We now have the base_label correct.
            # We will use the label walker to generate the canonical labeling of polygons.
            w = s.walker()
            w.find_all_labels()
            s.relabel(w.label_dictionary(), in_place=True)
            s.set_immutable()
            return s

        def rel_deformation(self, deformation, local=False, limit=100):
            r"""
            Perform a rel deformation of the surface and return the result.

            This algorithm currently assumes that all polygons affected by this deformation are
            triangles. That should be fixable in the future.

            INPUT:

            - ``deformation`` (dictionary) - A dictionary mapping singularities of
              the surface to deformation vectors (in some 2-dimensional vector
              space). The rel deformation being done will move the singularities
              (relative to each other) linearly to the provided vector for each
              vertex. If a singularity is not included in the dictionary then the
              vector will be treated as zero.

            - ``local`` - (boolean) - If true, the algorithm attempts to deform all
              the triangles making up the surface without destroying any of them.
              So, the area of the triangle must be positive along the full interval
              of time of the deformation.  If false, then the deformation must have
              a particular form: all vectors for the deformation must be parallel.
              In this case we achieve the deformation with the help of the SL(2,R)
              action and Delaunay triangulations.

            - ``limit`` (integer) - Restricts the length of the size of SL(2,R)
              deformations considered. The algorithm should be roughly worst time
              linear in limit.

            .. TODO::

                - Support arbitrary rel deformations.
                - Remove the requirement that triangles be used.

            EXAMPLES::

                sage: from flatsurf import *
                sage: s = translation_surfaces.arnoux_yoccoz(4)
                sage: field = s.base_ring()
                sage: a = field.gen()
                sage: V = VectorSpace(field,2)
                sage: deformation1 = {s.singularity(0,0):V((1,0))}
                doctest:warning
                ...
                UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
                sage: s1 = s.rel_deformation(deformation1).canonicalize()
                sage: deformation2 = {s.singularity(0,0):V((a,0))}
                sage: s2 = s.rel_deformation(deformation2).canonicalize()
                sage: m = Matrix([[a,0],[0,~a]])
                sage: s2.cmp((m*s1).canonicalize())
                0
            """
            s = self
            # Find a common field
            field = s.base_ring()
            for singularity, v in deformation.items():
                if v.parent().base_field() != field:
                    from sage.structure.element import get_coercion_model

                    cm = get_coercion_model()
                    field = cm.common_parent(field, v.parent().base_field())
            from sage.modules.free_module import VectorSpace

            vector_space = VectorSpace(field, 2)

            from collections import defaultdict

            vertex_deformation = defaultdict(
                vector_space.zero
            )  # dictionary associating the vertices.
            deformed_labels = set()  # list of polygon labels being deformed.

            for singularity, vect in deformation.items():
                for label, coordinates in singularity.representatives():
                    v = self.polygon(label).get_point_position(coordinates).get_vertex()
                    vertex_deformation[(label, v)] = vect
                    deformed_labels.add(label)
                    assert s.polygon(label).num_edges() == 3

            from flatsurf.geometry.polygon import wedge_product, ConvexPolygons

            if local:
                ss = s.copy(mutable=True, new_field=field)
                us = ss.underlying_surface()

                P = ConvexPolygons(field)
                for label in deformed_labels:
                    polygon = s.polygon(label)
                    a0 = vector_space(polygon.vertex(1))
                    b0 = vector_space(polygon.vertex(2))
                    v0 = vector_space(vertex_deformation[(label, 0)])
                    v1 = vector_space(vertex_deformation[(label, 1)])
                    v2 = vector_space(vertex_deformation[(label, 2)])
                    a1 = v1 - v0
                    b1 = v2 - v0
                    # We deform by changing the triangle so that its vertices 1 and 2 have the form
                    # a0+t*a1 and b0+t*b1
                    # respectively. We are deforming from t=0 to t=1.
                    # We worry that the triangle degenerates along the way.
                    # The area of the deforming triangle has the form
                    # A0 + A1*t + A2*t^2.
                    A0 = wedge_product(a0, b0)
                    A1 = wedge_product(a0, b1) + wedge_product(a1, b0)
                    A2 = wedge_product(a1, b1)
                    if A2 != field.zero():
                        # Critical point of area function
                        c = A1 / (-2 * A2)
                        if field.zero() < c and c < field.one():
                            if A0 + A1 * c + A2 * c**2 <= field.zero():
                                raise ValueError(
                                    "Triangle with label %r degenerates at critical point before endpoint"
                                    % label
                                )
                    if A0 + A1 + A2 <= field.zero():
                        raise ValueError(
                            "Triangle with label %r degenerates at or before endpoint"
                            % label
                        )
                    # Triangle does not degenerate.
                    us.change_polygon(
                        label, P(vertices=[vector_space.zero(), a0 + a1, b0 + b1])
                    )
                return ss

            else:  # Non local deformation
                # We can only do this deformation if all the rel vector are parallel.
                # Check for this.
                nonzero = None
                for singularity, vect in deformation.items():
                    vvect = vector_space(vect)
                    if vvect != vector_space.zero():
                        if nonzero is None:
                            nonzero = vvect
                        else:
                            assert (
                                wedge_product(nonzero, vvect) == 0
                            ), "In non-local deformation all deformation vectos must be parallel"
                assert nonzero is not None, "Deformation appears to be trivial."
                from sage.matrix.constructor import Matrix

                m = Matrix([[nonzero[0], -nonzero[1]], [nonzero[1], nonzero[0]]])
                mi = ~m
                g = Matrix([[1, 0], [0, 2]], ring=field)
                prod = m * g * mi
                ss = None
                k = 0
                while True:
                    if ss is None:
                        ss = s.copy(mutable=True, new_field=field)
                    else:
                        # In place matrix deformation
                        ss.apply_matrix(prod)
                    ss.delaunay_triangulation(direction=nonzero, in_place=True)
                    deformation2 = {}
                    for singularity, vect in deformation.items():
                        found_start = None
                        for label, coordinates in singularity.representatives():
                            v = (
                                s.polygon(label)
                                .get_point_position(coordinates)
                                .get_vertex()
                            )
                            if (
                                wedge_product(s.polygon(label).edge(v), nonzero) >= 0
                                and wedge_product(
                                    nonzero, -s.polygon(label).edge((v + 2) % 3)
                                )
                                > 0
                            ):
                                found_start = (label, v)
                                found = None
                                for vv in range(3):
                                    if (
                                        wedge_product(ss.polygon(label).edge(vv), nonzero)
                                        >= 0
                                        and wedge_product(
                                            nonzero, -ss.polygon(label).edge((vv + 2) % 3)
                                        )
                                        > 0
                                    ):
                                        found = vv
                                        deformation2[ss.singularity(label, vv)] = vect
                                        break
                                assert found is not None
                                break
                        assert found_start is not None
                    try:
                        sss = ss.rel_deformation(deformation2, local=True)
                    except ValueError:
                        k += 1
                        if limit is not None and k >= limit:
                            raise Exception("exceeded limit iterations")
                        continue

                    sss.apply_matrix(mi * g ** (-k) * m)
                    sss.delaunay_triangulation(direction=nonzero, in_place=True)
                    return sss

        def j_invariant(self):
            r"""
            Return the Kenyon-Smillie J-invariant of this translation surface.

            It is assumed that the coordinates are defined over a number field.

            EXAMPLES::

                sage: from flatsurf import *
                sage: O = translation_surfaces.regular_octagon()
                sage: O.j_invariant()
                (
                          [2 2]
                (0), (0), [2 1]
                )
            """
            it = self.label_iterator()
            lab = next(it)
            P = self.polygon(lab)
            Jxx, Jyy, Jxy = P.j_invariant()
            for lab in it:
                xx, yy, xy = self.polygon(lab).j_invariant()
                Jxx += xx
                Jyy += yy
                Jxy += xy
            return (Jxx, Jyy, Jxy)

        def erase_marked_points(self):
            r"""
            Return an isometric or similar surface with a minimal number of regular
            vertices of angle 2π.

            EXAMPLES::

                sage: import flatsurf

                sage: G = SymmetricGroup(4)
                sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
                sage: S.stratum()
                H_2(2, 0)
                sage: S.erase_marked_points().stratum() # optional: pyflatsurf  # long time (1s)  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
                H_2(2)

                sage: for (a,b,c) in [(1,4,11), (1,4,15), (3,4,13)]: # long time (10s), optional: pyflatsurf
                ....:     T = flatsurf.polygons.triangle(a,b,c)
                ....:     S = flatsurf.similarity_surfaces.billiard(T)
                ....:     S = S.minimal_cover("translation")
                ....:     print(S.erase_marked_points().stratum())
                H_6(10)
                H_6(2^5)
                H_8(12, 2)

            If the surface had no marked points then it is returned unchanged by this
            function::

                sage: O = flatsurf.translation_surfaces.regular_octagon()
                sage: O.erase_marked_points() is O
                True

            TESTS:

            Verify that https://github.com/flatsurf/flatsurf/issues/263 has been resolved::

                sage: from flatsurf import EquiangularPolygons, similarity_surfaces
                sage: E = EquiangularPolygons((10, 8, 3, 1, 1, 1))
                sage: P = E((1, 1, 2, 4), normalized=True)
                sage: B = similarity_surfaces.billiard(P, rational=True)
                sage: S = B.minimal_cover(cover_type="translation")
                sage: S = S.erase_marked_points() # long time (3s), optional: pyflatsurf

            ::

                sage: from flatsurf import EquiangularPolygons, similarity_surfaces
                sage: E = EquiangularPolygons((10, 7, 2, 2, 2, 1))
                sage: P = E((1, 1, 2, 3), normalized=True)
                sage: B = similarity_surfaces.billiard(P, rational=True)
                sage: S_mp = B.minimal_cover(cover_type="translation")
                sage: S = S_mp.erase_marked_points() # long time (3s), optional: pyflatsurf

            """
            if all(a != 1 for a in self.angles()):
                # no 2π angle
                return self
            from flatsurf.geometry.pyflatsurf_conversion import from_pyflatsurf, to_pyflatsurf

            S = to_pyflatsurf(self)
            S.delaunay()
            S = S.eliminateMarkedPoints().surface()
            S.delaunay()
            return from_pyflatsurf(S)
