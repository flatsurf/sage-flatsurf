r"""
The category of similarity surfaces.

This provides shared functionality for all surfaces in sage-flatsurf that are
built from Euclidean polygons that are glued by similarities, i.e., identified
edges can be transformed into each other by application of rotation and
homothety (dilation) and translation.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: C = MutableOrientedSimilaritySurface(QQ).category()

    sage: from flatsurf.geometry.categories import SimilaritySurfaces
    sage: C.is_subcategory(SimilaritySurfaces())
    True

The easiest way to construct a similarity surface is to use the pre-built
constructions from
:class:`flatsurf.geometry.similarity_surface_generators.SimilaritySurfaceGenerators`::

    sage: from flatsurf import polygon, similarity_surfaces
    sage: P = polygon(vertices=[(0,0), (2,0), (1,4), (0,5)])
    sage: similarity_surfaces.self_glued_polygon(P)
    Half-Translation Surface in Q_0(0, -1^4) built from a quadrilateral

The second way is to build a surface (using e.g.
:class:`flatsurf.geometry.surface.MutableOrientedSimilaritySurface`)::

    sage: P = polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = MutableOrientedSimilaritySurface(QQ)
    sage: S.add_polygon(P)
    0
    sage: S.add_polygon(2*P)
    1
    sage: S.add_polygon(3*P)
    2
    sage: S.glue((0, 1), (1, 3))
    sage: S.glue((0, 0), (2, 2))
    sage: S.glue((0, 2), (2, 0))
    sage: S.glue((0, 3), (1, 1))
    sage: S.glue((1, 2), (2, 1))
    sage: S.glue((1, 0), (2, 3))
    sage: S
    Surface built from 3 squares

To perform a sanity check on the obtained surface, you can run its test
suite::

    sage: TestSuite(S).run()

In the following example, we attempt to build a broken surface but edges get
unglued automatically unglued::

    sage: S.glue((0, 0), (0, 3))
    sage: S.glue((0, 1), (0, 3))
    sage: S.glue((0, 2), (0, 3))

    sage: S.gluings()
    (((0, 2), (0, 3)), ((0, 3), (0, 2)))

    sage: TestSuite(S).run()

Here, we build a surface with boundary::

    sage: P = polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
    sage: S = MutableOrientedSimilaritySurface(QQ)
    sage: S.add_polygon(P)
    0
    sage: TestSuite(S).run()

"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
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
# ****************************************************************************

from flatsurf.geometry.categories.surface_category import SurfaceCategory, SurfaceCategoryWithAxiom
from sage.categories.category_with_axiom import all_axioms
from sage.misc.cachefunc import cached_method
from sage.all import ZZ, QQ, AA


class SimilaritySurfaces(SurfaceCategory):
    r"""
    The category of surfaces built from polygons with edges identified by
    similarities.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import SimilaritySurfaces
        sage: SimilaritySurfaces()
        Category of similarity surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.real_projective_polygonal_surfaces import RealProjectivePolygonalSurfaces
        return [RealProjectivePolygonalSurfaces()]

    class ParentMethods:
        def refined_category(self):
            r"""
            Return the smallest subcategory that this surface is in by consulting
            how its edges are glued.

            The result of this method can be fed to ``_refine_category_`` to
            change the category of the surface (and enable functionality
            specific to the smaller classes of surfaces.)


            .. NOTE::

                If a surface cannot implement the various ``is_`` methods used in
                the implementation of this method (i.e., if any of them throws a
                ``NotImplementedError``,) then this method ``refined_category``
                must be overridden to skip that check. We don't want to actively
                catch a ``NotImplementedError`` and instead encourage authors
                to explicitly select the category their surfaces lives in.

            EXAMPLES::

                sage: from flatsurf import MutableOrientedSimilaritySurface
                sage: S = MutableOrientedSimilaritySurface(QQ)

                sage: from flatsurf import polygons
                sage: S.add_polygon(polygons.square(), label=0)
                0
                sage: S.refined_category()
                Category of connected with boundary finite type translation surfaces

                sage: S.glue((0, 0), (0, 2))
                sage: S.glue((0, 1), (0, 3))
                sage: S.refined_category()
                Category of connected without boundary finite type translation surfaces

            """
            from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces
            category = PolygonalSurfaces.ParentMethods.refined_category(self)

            if self.is_cone_surface():
                from flatsurf.geometry.categories.cone_surfaces import ConeSurfaces
                category &= ConeSurfaces()

            if self.is_dilation_surface():
                from flatsurf.geometry.categories.dilation_surfaces import DilationSurfaces
                category &= DilationSurfaces()

                if self.is_dilation_surface(positive=True):
                    category &= DilationSurfaces().Positive()

                    if self.is_translation_surface():
                        from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
                        category &= TranslationSurfaces()
                elif self.is_translation_surface(positive=False):
                    from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
                    category &= HalfTranslationSurfaces()

            if "Rational" not in category.axioms():
                if self.is_rational_surface():
                    category = category.Rational()

            return category

        def is_cone_surface(self):
            r"""
            Return whether this surface is a cone surface, i.e., glued edges
            can be transformed into each other with isometries.

            .. NOTE::

                This method is used to determine whether this surface is in the
                category of :class:`ConeSurfaces`. Surfaces can override this
                method to perform specialized logic, see the note in
                :mod:`flatsurf.geometry.categories` for performance considerations.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_cone_surface()
                True

            """
            if self.is_translation_surface():
                return True

            raise NotImplementedError("surface does not implement is_cone_surface()")

        def is_dilation_surface(self, positive=False):
            r"""
            Return whether this surface is a dilation surface, i.e., whether
            glued edges can be transformed into each other by translation
            followed by a dilation (multiplication by a diagonal matrix.)

            .. NOTE::

                This method is used to determine whether this surface is in the
                category of :class:`DilationSurfaces` or
                :class:`DilationSurfaces.Positive`. Surfaces can override this
                method to perform specialized logic, see the note in
                :mod:`flatsurf.geometry.categories` for performance
                considerations.

            INPUT:

            - ``positive`` -- a boolean (default: ``False``); whether the
              entries of the diagonal matrix must be positive or are allowed to
              be negative.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_dilation_surface()
                True
                sage: S.is_dilation_surface(positive=True)
                False

            """
            if self.is_translation_surface(positive=positive):
                return True

            raise NotImplementedError("surface does not implement is_dilation_surface()")

        def is_translation_surface(self, positive=True):
            r"""
            Return whether this surface is a translation surface, i.e., glued
            edges can be transformed into each other by translations.

            This method must be implemented if this surface is a dilation surface.

            .. NOTE::

                This method is used to determine whether this surface is in the
                category of :class:`HalfTranslationSurfaces` or
                :class:`TranslationSurfaces`. Surfaces can override this method
                to perform specialized logic, see the note in
                :mod:`flatsurf.geometry.categories` for performance
                considerations.

            INPUT:

            - ``positive`` -- a boolean (default: ``True``); whether the
              transformation must be a translation or is allowed to be a
              half-translation, i.e., a translation followed by a reflection in
              a point (equivalently, a rotation by π.)

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_translation_surface()
                False
                sage: S.is_translation_surface(False)
                True

            ::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.is_translation_surface()
                True

            """
            raise NotImplementedError("surface does not implement is_translation_surface()")

        def is_rational_surface(self):
            r"""
            Return whether this surface is a rational surface, i.e., all the
            rotational part of all gluings is a rational multiple of π.

            .. NOTE::

                This method is used to determine whether this surface satisfies
                the :class:`Rational` axiom. Surfaces can override this method
                to perform specialized logic, see the note in
                :mod:`flatsurf.geometry.categories` for performance
                considerations.

            EXAMPLES::

                sage: from flatsurf import polygons, similarity_surfaces
                sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                sage: S = similarity_surfaces.self_glued_polygon(P)
                sage: S.is_rational_surface()
                True

            """
            if self.is_translation_surface():
                return True

            raise NotImplementedError("surface does not implement is_rational_surface()")

        def _mul_(self, matrix, switch_sides=True):
            r"""
            EXAMPLES::

                sage: from flatsurf import *
                sage: s = translation_surfaces.infinite_staircase()
                sage: s
                The infinite staircase
                sage: m=Matrix([[1,2],[0,1]])
                sage: s2=m*s
                sage: TestSuite(s2).run()
                sage: s2.polygon(0)
                polygon(vertices=[(0, 0), (1, 0), (3, 1), (2, 1)])

            Testing multiplication by a matrix with negative determinant::

                sage: from flatsurf import *
                sage: ds1 = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
                sage: ds1.polygon(0)
                polygon(vertices=[(0, 0), (1/2, 0), (1, 1/3), (1, 1), (3/4, 1), (0, 4/5)])
                sage: m = matrix(QQ, [[0, 1], [1, 0]]) # maps (x,y) to (y, x)
                sage: ds2 = m*ds1
                sage: ds2.polygon(0)
                polygon(vertices=[(0, 0), (4/5, 0), (1, 3/4), (1, 1), (1/3, 1), (0, 1/2)])
            """
            if not switch_sides:
                raise NotImplementedError

            from sage.structure.element import is_Matrix
            if not is_Matrix(matrix):
                raise NotImplementedError("Only implemented for matrices.")
            if not matrix.dimensions != (2, 2):
                raise NotImplementedError("Only implemented for 2x2 matrices.")
            from flatsurf.geometry.half_dilation_surface import GL2RImageSurface
            return GL2RImageSurface(self, matrix)

    class Oriented(SurfaceCategoryWithAxiom):
        class ParentMethods:
            @cached_method
            def edge_matrix(self, p, e=None):
                r"""
                Returns the 2x2 matrix representing a similarity which when applied to the polygon with label `p`
                makes it so the edge `e` can be glued to its opposite edge by translation.

                If `e` is not provided, then `p` should be a pair consisting of a polygon label and an edge.

                EXAMPLES::

                    sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                    sage: s = SimilaritySurfaceGenerators.example()
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (2, -2), (2, 0)])
                    sage: s.polygon(1)
                    polygon(vertices=[(0, 0), (2, 0), (1, 3)])
                    sage: s.opposite_edge(0,0)
                    (1, 1)
                    sage: m = s.edge_matrix(0, 0)
                    sage: m
                    [   1  1/2]
                    [-1/2    1]
                    sage: m * vector((2,-2)) == -vector((-1, 3))
                    True

                """
                if e is None:
                    import warnings

                    warnings.warn("edge_matrix will now only take two arguments")
                    p, e = p
                u = self.polygon(p).edge(e)
                pp, ee = self.opposite_edge(p, e)
                v = self.polygon(pp).edge(ee)

                # note the orientation, it is -v and not v
                from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
                from sage.matrix.matrix_space import MatrixSpace

                return similarity_from_vectors(u, -v, MatrixSpace(self.base_ring(), 2))

            def _an_element_(self):
                label = next(iter(self.labels()))
                polygon = self.polygon(label)

                # We use a point that can be constructed without problems on an
                # infinite surface.
                if polygon.is_convex():
                    coordinates = polygon.centroid()
                else:
                    # Sometimes, this is not implemented because it requires the edge
                    # transformation to be known, so we prefer the centroid.
                    coordinates = polygon.edge(0) / 2
                return self(label, coordinates)

            @cached_method
            def _matrix_space(self):
                from sage.matrix.matrix_space import MatrixSpace

                return MatrixSpace(self.base_ring(), 2)

            def underlying_surface(self):
                r"""
                Return this surface.

                EXAMPLES::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.underlying_surface() is S
                    doctest:warning
                    ...
                    UserWarning: underlying_surface() has been deprecated and will be removed in a future version of sage-flatsurf; this function has no effect anymore since there is no distinction between a surface and its underlying surface anymore
                    True

                """
                import warnings
                warnings.warn("underlying_surface() has been deprecated and will be removed in a future version of sage-flatsurf; this function has no effect anymore since there is no distinction between a surface and its underlying surface anymore")

                return self

            def num_singularities(self):
                r"""
                EXAMPLES::

                    sage: from flatsurf import *

                    sage: translation_surfaces.regular_octagon().num_singularities()
                    1

                    sage: S = SymmetricGroup(4)
                    sage: r = S('(1,2)(3,4)')
                    sage: u = S('(2,3)')
                    sage: translation_surfaces.origami(r,u).num_singularities()
                    2

                    sage: S = SymmetricGroup(8)
                    sage: r = S('(1,2,3,4,5,6,7,8)')
                    sage: u = S('(1,8,5,4)(2,3)(6,7)')
                    sage: translation_surfaces.origami(r,u).num_singularities()
                    4
                """
                if not self.is_finite_type():
                    raise ValueError("the method only work for finite surfaces")

                # NOTE:
                # the very same code is implemented in the method angles (translation
                # surfaces). we should factor out the code
                edges = set(
                    (p, e)
                    for p in self.labels()
                    for e in range(self.polygon(p).num_edges())
                )

                n = ZZ(0)
                while edges:
                    p, e = edges.pop()
                    n += 1
                    ee = (e - 1) % self.polygon(p).num_edges()
                    p, e = self.opposite_edge(p, ee)
                    while (p, e) in edges:
                        edges.remove((p, e))
                        ee = (e - 1) % self.polygon(p).num_edges()
                        p, e = self.opposite_edge(p, ee)
                return n

            def edge_transformation(self, p, e):
                r"""
                Return the similarity bringing the provided edge to the opposite edge.

                EXAMPLES::

                    sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
                    sage: s = SimilaritySurfaceGenerators.example()
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (2, -2), (2, 0)])
                    sage: s.polygon(1)
                    polygon(vertices=[(0, 0), (2, 0), (1, 3)])
                    sage: s.opposite_edge(0,0)
                    (1, 1)
                    sage: g = s.edge_transformation(0,0)
                    sage: g((0,0))
                    (1, 3)
                    sage: g((2,-2))
                    (2, 0)
                """
                from flatsurf.geometry.similarity import SimilarityGroup
                G = SimilarityGroup(self.base_ring())
                q = self.polygon(p)
                a = q.vertex(e)
                b = q.vertex(e + 1)
                # This is the similarity carrying the origin to a and (1,0) to b:
                g = G(b[0] - a[0], b[1] - a[1], a[0], a[1])

                pp, ee = self.opposite_edge(p, e)
                qq = self.polygon(pp)
                # Be careful here: opposite vertices are identified
                aa = qq.vertex(ee + 1)
                bb = qq.vertex(ee)
                # This is the similarity carrying the origin to aa and (1,0) to bb:
                gg = G(bb[0] - aa[0], bb[1] - aa[1], aa[0], aa[1])

                # This is the similarity carrying (a,b) to (aa,bb):
                return gg / g

            def set_vertex_zero(self, label, v, in_place=False):
                r"""
                Applies a combinatorial rotation to the polygon with the provided label.

                This makes what is currently vertex v of this polygon vertex 0. In other words,
                what is currently vertex (or edge) e will now become vertex (e-v)%n where
                n is the number of sides of the polygon.

                For the updated polygons, the polygons will be translated so that vertex
                0 is the origin.

                EXAMPLES:

                Example with polygon glued to another polygon::

                    sage: from flatsurf import *
                    sage: s = translation_surfaces.veech_double_n_gon(4)
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                    sage: [s.opposite_edge(0,i) for i in range(4)]
                    [(1, 0), (1, 1), (1, 2), (1, 3)]
                    sage: ss = s.set_vertex_zero(0,1)
                    sage: ss.polygon(0)
                    polygon(vertices=[(0, 0), (0, 1), (-1, 1), (-1, 0)])
                    sage: [ss.opposite_edge(0,i) for i in range(4)]
                    [(1, 1), (1, 2), (1, 3), (1, 0)]
                    sage: TestSuite(ss).run()

                Example with polygon glued to self::

                    sage: s = translation_surfaces.veech_2n_gon(2)
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                    sage: [s.opposite_edge(0,i) for i in range(4)]
                    [(0, 2), (0, 3), (0, 0), (0, 1)]
                    sage: ss = s.set_vertex_zero(0,3)
                    sage: ss.polygon(0)
                    polygon(vertices=[(0, 0), (0, -1), (1, -1), (1, 0)])
                    sage: [ss.opposite_edge(0,i) for i in range(4)]
                    [(0, 2), (0, 3), (0, 0), (0, 1)]
                    sage: TestSuite(ss).run()
                """
                if in_place:
                    us = self
                    if not us.is_mutable():
                        raise ValueError(
                            "set_vertex_zero can only be done in_place for a mutable surface."
                        )
                    p = us.polygon(label)
                    n = p.num_edges()
                    if not (0 <= v < n):
                        raise ValueError
                    glue = []
                    from flatsurf.geometry.polygon import ConvexPolygons
                    P = ConvexPolygons(us.base_ring())
                    pp = P(edges=[p.edge((i + v) % n) for i in range(n)])

                    for i in range(n):
                        e = (v + i) % n
                        ll, ee = us.opposite_edge(label, e)
                        if ll == label:
                            ee = (ee + n - v) % n
                        glue.append((ll, ee))

                    us.remove_polygon(label)
                    us.add_polygon(pp, label=label)
                    for e, cross in enumerate(glue):
                        us.glue((label, e), cross)
                    return self
                else:
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    return MutableOrientedSimilaritySurface.from_surface(self).set_vertex_zero(label, v, in_place=True)

            def _label_comparator(self):
                r"""
                Return a LabelComparator, which provides a fixed total ordering on the polygon labels.
                """
                try:
                    return self._lc
                except AttributeError:
                    from flatsurf.geometry.surface import LabelComparator
                    self._lc = LabelComparator()
                    return self._lc

            def relabel(self, relabeling_map, in_place=False):
                r"""
                Attempt to relabel the polygons according to a relabeling_map, which takes as input
                a current label and outputs a new label for the same polygon. The method returns a pair
                (surface,success) where surface is the relabeled surface, and success is a boolean value
                indicating the success of the operation. The operation will fail if the implementation of the
                underlying surface does not support labels used in the image of the relabeling map. In this case,
                other (arbitrary) labels will be used to replace the labels of the surface, and the resulting
                surface should still be okay.

                Currently, the relabeling_map must be a dictionary.

                If in_place is True then the relabeling is done to the current surface, otherwise a
                mutable copy is made before relabeling.

                ToDo:
                  - Allow relabeling_map to be a function rather than just a dictionary.
                    This will allow it to work for infinite surfaces.

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: s=translation_surfaces.veech_double_n_gon(5)
                    sage: ss,valid=s.relabel({0:1,1:2})
                    sage: valid
                    True
                    sage: ss.base_label()
                    1
                    sage: ss.opposite_edge(1,0)
                    (2, 0)
                    sage: len(ss.polygons())
                    2
                    sage: TestSuite(ss).run()
                """
                if in_place:
                    us = self
                    if not us.is_mutable():
                        raise ValueError(
                            "Your surface is not mutable, so can not be relabeled in place."
                        )
                    if not isinstance(relabeling_map, dict):
                        raise NotImplementedError(
                            "Currently relabeling is only implemented via a dictionary."
                        )
                    domain = set()
                    codomain = set()
                    data = {}
                    for l1, l2 in relabeling_map.items():
                        p = us.polygon(l1)
                        glue = []
                        for e in range(p.num_edges()):
                            ll, ee = us.opposite_edge(l1, e)
                            try:
                                lll = relabeling_map[ll]
                            except KeyError:
                                lll = ll
                            glue.append((lll, ee))
                        data[l2] = (p, glue)
                        domain.add(l1)
                        codomain.add(l2)
                    if len(domain) != len(codomain):
                        raise ValueError(
                            "The relabeling_map must be injective. Received "
                            + str(relabeling_map)
                        )
                    changed_labels = domain.intersection(codomain)
                    added_labels = codomain.difference(domain)
                    removed_labels = domain.difference(codomain)
                    # Pass to add_polygons
                    relabel_errors = {}
                    for l2 in added_labels:
                        p, glue = data[l2]
                        l3 = us.add_polygon(p, label=l2)
                        if not l2 == l3:
                            # This means the label l2 could not be added for some reason.
                            # Perhaps the implementation does not support this type of label.
                            # Or perhaps there is already a polygon with this label.
                            relabel_errors[l2] = l3
                    # Pass to change polygons
                    for l2 in changed_labels:
                        p, glue = data[l2]
                        us.remove_polygon(l2)
                        us.add_polygon(p, label=l2)
                        us.replace_polygon(l2, p)
                    # Deal with the base_label
                    base_label = us.base_label()
                    if base_label in relabeling_map:
                        base_label = relabeling_map[base_label]
                        if base_label in relabel_errors:
                            base_label = relabel_errors[base_label]
                        us.set_base_label(base_label)
                    # Pass to remove polygons:
                    for l1 in removed_labels:
                        us.remove_polygon(l1)
                    # Pass to update the edge gluings
                    if len(relabel_errors) == 0:
                        # No problems. Update the gluings.
                        for l2 in codomain:
                            p, glue = data[l2]
                            for e, cross in enumerate(glue):
                                us.glue((l2, e), cross)
                    else:
                        # Use the gluings provided by relabel_errors when necessary
                        for l2 in codomain:
                            p, glue = data[l2]
                            for e in range(p.num_edges()):
                                ll, ee = glue[e]
                                try:
                                    # First try the error dictionary
                                    us.glue((l2, e), (relabel_errors[ll], ee))
                                except KeyError:
                                    us.glue((l2, e), (ll, ee))
                    return self, len(relabel_errors) == 0
                else:
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    return MutableOrientedSimilaritySurface.from_surface(self).relabel(relabeling_map, in_place=True)

            def copy(
                self,
                relabel=False,
                mutable=False,
                lazy=None,
                new_field=None,
                optimal_number_field=False,
            ):
                r"""
                Returns a copy of this surface. The method takes several flags to modify how the copy is taken.

                If relabel is True, then instead of returning an exact copy, it returns a copy indexed by the
                non-negative integers. This uses the Surface_list implementation. If relabel is False (default),
                then we return an exact copy. The returned surface uses the Surface_dict implementation.

                The mutability flag returns if the resulting surface should be mutable or not. By default, the
                resulting surface will not be mutable.

                If lazy is True, then the surface is copied by reference. This is the only type of copy
                possible for infinite surfaces. The parameter defaults to False for finite surfaces, and
                defaults to True for infinite surfaces.

                The new_field parameter can be used to place the vertices in a larger field than the basefield
                for the original surface.

                The optimal_number_field option can be used to find a best NumberField containing the
                (necessarily finite) surface.

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: ss=translation_surfaces.ward(3)
                    sage: ss.is_mutable()
                    False
                    sage: s=ss.copy(mutable=True)
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead.
                    sage: s.is_mutable()
                    True
                    sage: TestSuite(s).run()
                    sage: s == ss
                    False

                    sage: # Changing the base field
                    sage: from flatsurf import *
                    sage: s=translation_surfaces.veech_double_n_gon(5)
                    sage: ss=s.copy(mutable=False,new_field=AA)
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead. Use set_immutable() to make the resulting surface immutable. Use change_ring() to change the field over which the surface is defined.
                    sage: TestSuite(ss).run()
                    sage: ss.base_ring()
                    Algebraic Real Field

                    sage: # Optimization of number field
                    sage: from flatsurf import *
                    sage: s = translation_surfaces.arnoux_yoccoz(3)
                    sage: ss = s.copy(new_field=AA).copy(optimal_number_field=True)
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead. Use set_immutable() to make the resulting surface immutable. Use change_ring() to change the field over which the surface is defined.
                    doctest:warning
                    ...
                    UserWarning: copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead. Use set_immutable() to make the resulting surface immutable. There is currently no replacement for optimal number field. If you are relying on this features, let the authors of sage-flatsurf know and we will try to make it available again.
                    sage: TestSuite(ss).run()
                    sage: ss.base_ring().discriminant()
                    -44
                """
                message = "copy() has been deprecated and will be removed from a future version of sage-flatsurf; for surfaces of finite type use MutableOrientedSimilaritySurface.from_surface() instead."

                if not mutable:
                    message += " Use set_immutable() to make the resulting surface immutable."

                if relabel:
                    message += " Use relabel({old: new for (new, old) in enumerate(surface.labels())}) for integer labels."

                if not self.is_finite_type():
                    message += " However, there is no immediate replacement for lazy copying of infinite surfaces. Have a look at the implementation of flatsurf.geometry.delaunay.LazyMutableSurface and adapt it to your needs."

                if new_field is not None:
                    message += " Use change_ring() to change the field over which the surface is defined."

                if optimal_number_field:
                    message += " There is currently no replacement for optimal number field. If you are relying on this features, let the authors of sage-flatsurf know and we will try to make it available again."

                import warnings
                warnings.warn(message)

                category = self.category()
                s = None  # This will be the surface we copy. (Likely we will set s=self below.)
                if new_field is not None and optimal_number_field:
                    raise ValueError(
                        "You can not set a new_field and also set optimal_number_field=True."
                    )
                if optimal_number_field is True:
                    if not self.is_finite_type():
                        raise NotImplementedError(
                            "can only optimize_number_field for a finite surface"
                        )
                    if lazy:
                        raise NotImplementedError(
                            "lazy copying is unavailable when optimize_number_field=True"
                        )
                    coordinates_AA = []
                    for label, p in zip(self.labels(), self.polygons()):
                        for e in p.edges():
                            coordinates_AA.append(AA(e[0]))
                            coordinates_AA.append(AA(e[1]))
                    from sage.rings.qqbar import number_field_elements_from_algebraics

                    field, coordinates_NF, hom = number_field_elements_from_algebraics(
                        coordinates_AA, minimal=True
                    )
                    if field is QQ:
                        new_field = QQ
                        # We pretend new_field = QQ was passed as a parameter.
                        # It will now get picked up by the "if new_field is not None:" line below.
                    else:
                        # Unfortunately field doesn't come with an real embedding (which is given by hom!)
                        # So, we make a copy of the field, and add the embedding.
                        from sage.all import NumberField
                        field2 = NumberField(
                            field.polynomial(), name="a", embedding=hom(field.gen())
                        )
                        # The following converts from field to field2:
                        hom2 = field.hom(im_gens=[field2.gen()])

                        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                        ss = MutableOrientedSimilaritySurface(field2)
                        index = 0
                        from flatsurf.geometry.polygon import ConvexPolygons
                        P = ConvexPolygons(field2)
                        for label, p in zip(self.labels(), self.polygons()):
                            new_edges = []
                            for i in range(p.num_edges()):
                                new_edges.append(
                                    (
                                        hom2(coordinates_NF[index]),
                                        hom2(coordinates_NF[index + 1]),
                                    )
                                )
                                index += 2
                            pp = P(edges=new_edges)
                            ss.add_polygon(pp, label=label)
                        ss.set_base_label(self.base_label())
                        for (l1, e1), (l2, e2) in self.gluings():
                            ss.glue((l1, e1), (l2, e2))
                        s = ss
                        if not relabel:
                            if not mutable:
                                s.set_immutable()
                            return s
                        # Otherwise we are supposed to relabel. We will make a relabeled copy of s below.
                if new_field is not None:
                    from flatsurf.geometry.surface import BaseRingChangedSurface

                    s = BaseRingChangedSurface(self, new_field)
                if s is None:
                    s = self
                if s.is_finite_type():
                    if relabel:
                        from flatsurf.geometry.surface import Surface_list
                        return Surface_list(surface=s, copy=not lazy, mutable=mutable, category=category, deprecation_warning=False)
                    else:
                        from flatsurf.geometry.surface import Surface_dict
                        return Surface_dict(surface=s, copy=not lazy, mutable=mutable, category=category, deprecation_warning=False)
                else:
                    if lazy is False:
                        raise ValueError("Only lazy copying available for infinite surfaces.")
                    if self.is_mutable():
                        raise ValueError(
                            "An infinite surface can only be copied if it is immutable."
                        )
                    if relabel:
                        from flatsurf.geometry.surface import Surface_list
                        return Surface_list(surface=s, copy=False, mutable=mutable, category=category, deprecation_warning=False)
                    else:
                        from flatsurf.geometry.surface import Surface_dict
                        return Surface_dict(surface=s, copy=False, mutable=mutable, category=category, deprecation_warning=False)

            def change_ring(self, ring):
                from flatsurf.geometry.surface import BaseRingChangedSurface
                return BaseRingChangedSurface(self, ring)

            def triangle_flip(self, l1, e1, in_place=False, test=False, direction=None):
                r"""
                Flips the diagonal of the quadrilateral formed by two triangles
                glued together along the provided edge (l1,e1). This can be broken
                into two steps: join along the edge to form a convex quadilateral,
                then cut along the other diagonal. Raises a ValueError if this
                quadrilateral would be non-convex. In this case no changes to the
                surface are made.

                The direction parameter defaults to (0,1). This is used to decide how
                the triangles being glued in are labeled. Let p1 be the triangle
                associated to label l1, and p2 be the triangle associated to l2
                but moved by a similarity to share the edge (l1,e1). Each triangle
                has a exactly one separatrix leaving a vertex which travels in the
                provided direction or its opposite. (For edges we only count as sepatrices
                traveling counter-clockwise around the triangle.) This holds for p1
                and p2 and the separatrices must point in opposite directions.

                The above description gives two new triangles t1 and t2 which must be
                glued in (obtained by flipping the diagonal of the quadrilateral).
                Up to swapping t1 and t2 we can assume the separatrix in t1 in the
                provided direction (or its opposite) points in the same direction as
                that of p1. Further up to cyclic permutation of vertex labels we can
                assume that the separatrices in p1 and t1 start at the vertex with the
                same index (an element of {0,1,2}). The same can be done for p2 and t2.
                We apply the label l1 to t1 and the label l2 to t2. This precisely
                determines how t1 and t2 should be used to replace p1 and p2.

                INPUT:

                - ``l1`` - label of polygon

                - ``e1`` - (integer) edge of the polygon

                - ``in_place`` (boolean) - If True do the flip to the current surface
                  which must be mutable. In this case the updated surface will be
                  returned.  Otherwise a mutable copy is made and then an edge is
                  flipped, which is then returned.

                - ``test`` (boolean) - If True we don't actually flip, and we return
                  True or False depending on whether or not the flip would be
                  successful.

                - ``direction`` (2-dimensional vector) - Defaults to (0,1). The choice
                  of this vector determines how the newly added triangles are labeled.

                EXAMPLES::

                    sage: from flatsurf import similarity_surfaces, MutableOrientedSimilaritySurface, polygon

                    sage: s = similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (1, 0), (0, 1)])
                    sage: s.triangle_flip(0, 0, test=True)
                    False
                    sage: s.triangle_flip(0, 1, test=True)
                    True
                    sage: s.triangle_flip(0, 2, test=True)
                    False

                    sage: s = similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
                    sage: s = MutableOrientedSimilaritySurface.from_surface(s)
                    sage: s.triangle_flip(0,0,in_place=True)
                    Traceback (most recent call last):
                    ...
                    ValueError: Gluing triangles along this edge yields a non-convex quadrilateral.
                    sage: s.triangle_flip(0,1,in_place=True)
                    Genus 0 Rational Cone Surface built from 2 isosceles triangles
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (1, 1), (0, 1)])
                    sage: s.polygon(1)
                    polygon(vertices=[(0, 0), (-1, -1), (0, -1)])
                    sage: s.gluings()
                    (((0, 0), (1, 0)), ((0, 1), (0, 2)), ((0, 2), (0, 1)), ((1, 0), (0, 0)), ((1, 1), (1, 2)), ((1, 2), (1, 1)))
                    sage: s.triangle_flip(0,2,in_place=True)
                    Traceback (most recent call last):
                    ...
                    ValueError: Gluing triangles along this edge yields a non-convex quadrilateral.

                    sage: p = polygon(edges=[(2,0),(-1,3),(-1,-3)])
                    sage: s = similarity_surfaces.self_glued_polygon(p)
                    sage: s = MutableOrientedSimilaritySurface.from_surface(s)
                    sage: s.triangle_flip(0,1,in_place=True)
                    Half-Translation Surface in Q_0(-1^4) built from a triangle

                    sage: from flatsurf.geometry.categories import DilationSurfaces
                    sage: s in DilationSurfaces()
                    True
                    sage: s.labels()
                    (0,)
                    sage: s.polygons()
                    (polygon(vertices=[(0, 0), (-3, -3), (-1, -3)]),)
                    sage: s.gluings()
                    (((0, 0), (0, 0)), ((0, 1), (0, 1)), ((0, 2), (0, 2)))
                    sage: TestSuite(s).run()

                """
                if test:
                    # Just test if the flip would be successful
                    p1 = self.polygon(l1)
                    if not p1.num_edges() == 3:
                        return False
                    l2, e2 = self.opposite_edge(l1, e1)
                    p2 = self.polygon(l2)
                    if not p2.num_edges() == 3:
                        return False
                    sim = self.edge_transformation(l2, e2)
                    hol = sim(p2.vertex((e2 + 2) % 3) - p1.vertex((e1 + 2) % 3))
                    from flatsurf.geometry.polygon import wedge_product

                    return (
                        wedge_product(p1.edge((e1 + 2) % 3), hol) > 0
                        and wedge_product(p1.edge((e1 + 1) % 3), hol) > 0
                    )

                if in_place:
                    s = self
                    if not s.is_mutable():
                        raise ValueError("surface must be mutable for in place triangle_flip")
                else:
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    s = MutableOrientedSimilaritySurface.from_surface(self)

                p1 = s.polygon(l1)
                if not p1.num_edges() == 3:
                    raise ValueError("The polygon with the provided label is not a triangle.")
                l2, e2 = s.opposite_edge(l1, e1)

                sim = s.edge_transformation(l2, e2)
                p2 = s.polygon(l2)
                if not p2.num_edges() == 3:
                    raise ValueError(
                        "The polygon opposite the provided edge is not a triangle."
                    )
                P = p1.parent()
                p2 = P(vertices=[sim(v) for v in p2.vertices()])

                if direction is None:
                    direction = s.vector_space()((0, 1))
                # Get vertices corresponding to separatices in the provided direction.
                v1 = p1.find_separatrix(direction=direction)[0]
                v2 = p2.find_separatrix(direction=direction)[0]
                # Our quadrilateral has vertices labeled:
                # * 0=p1.vertex(e1+1)=p2.vertex(e2)
                # * 1=p1.vertex(e1+2)
                # * 2=p1.vertex(e1)=p2.vertex(e2+1)
                # * 3=p2.vertex(e2+2)
                # Record the corresponding vertices of this quadrilateral.
                q1 = (3 + v1 - e1 - 1) % 3
                q2 = (2 + (3 + v2 - e2 - 1) % 3) % 4

                new_diagonal = p2.vertex((e2 + 2) % 3) - p1.vertex((e1 + 2) % 3)
                # This list will store the new triangles which are being glued in.
                # (Unfortunately, they may not be cyclically labeled in the correct way.)
                new_triangle = []
                try:
                    new_triangle.append(
                        P(edges=[p1.edge((e1 + 2) % 3), p2.edge((e2 + 1) % 3), -new_diagonal])
                    )
                    new_triangle.append(
                        P(edges=[p2.edge((e2 + 2) % 3), p1.edge((e1 + 1) % 3), new_diagonal])
                    )
                    # The above triangles would be glued along edge 2 to form the diagonal of the quadrilateral being removed.
                except ValueError:
                    raise ValueError(
                        "Gluing triangles along this edge yields a non-convex quadrilateral."
                    )

                # Find the separatrices of the two new triangles, and in particular which way they point.
                new_sep = []
                new_sep.append(new_triangle[0].find_separatrix(direction=direction)[0])
                new_sep.append(new_triangle[1].find_separatrix(direction=direction)[0])
                # The quadrilateral vertices corresponding to these separatrices are
                # new_sep[0]+1 and (new_sep[1]+3)%4 respectively.

                # i=0 if the new_triangle[0] should be labeled l1 and new_triangle[1] should be labeled l2.
                # i=1 indicates the opposite labeling.
                if new_sep[0] + 1 == q1:
                    assert (new_sep[1] + 3) % 4 == q2, (
                        "Bug: new_sep[1]=" + str(new_sep[1]) + " and q2=" + str(q2)
                    )
                    i = 0
                else:
                    assert (new_sep[1] + 3) % 4 == q1
                    assert new_sep[0] + 1 == q2
                    i = 1

                # These quantities represent the cyclic relabeling of triangles needed.
                cycle1 = (new_sep[i] - v1 + 3) % 3
                cycle2 = (new_sep[1 - i] - v2 + 3) % 3

                # This will be the new triangle with label l1:
                tri1 = P(
                    edges=[
                        new_triangle[i].edge(cycle1),
                        new_triangle[i].edge((cycle1 + 1) % 3),
                        new_triangle[i].edge((cycle1 + 2) % 3),
                    ]
                )
                # This will be the new triangle with label l2:
                tri2 = P(
                    edges=[
                        new_triangle[1 - i].edge(cycle2),
                        new_triangle[1 - i].edge((cycle2 + 1) % 3),
                        new_triangle[1 - i].edge((cycle2 + 2) % 3),
                    ]
                )
                # In the above, edge 2-cycle1 of tri1 would be glued to edge 2-cycle2 of tri2
                diagonal_glue_e1 = 2 - cycle1
                diagonal_glue_e2 = 2 - cycle2

                assert p1.find_separatrix(direction=direction) == tri1.find_separatrix(
                    direction=direction
                )
                assert p2.find_separatrix(direction=direction) == tri2.find_separatrix(
                    direction=direction
                )

                # Two opposite edges will not change their labels (label,edge) under our regluing operation.
                # The other two opposite ones will change and in fact they change labels.
                # The following finds them (there are two cases).
                # At the end of the if statement, the following will be true:
                # * new_glue_e1 and new_glue_e2 will be the edges of the new triangle with label l1 and l2 which need regluing.
                # * old_e1 and old_e2 will be the corresponding edges of the old triangles.
                # (Note that labels are swapped between the pair. The appending 1 or 2 refers to the label used for the triangle.)
                if p1.edge(v1) == tri1.edge(v1):
                    # We don't have to worry about changing gluings on edge v1 of the triangles with label l1
                    # We do have to worry about the following edge:
                    new_glue_e1 = (
                        3 - diagonal_glue_e1 - v1
                    )  # returns the edge which is neither diagonal_glue_e1 nor v1.
                    # This corresponded to the following old edge:
                    old_e1 = 3 - e1 - v1  # Again this finds the edge which is neither e1 nor v1
                else:
                    temp = (v1 + 2) % 3
                    assert p1.edge(temp) == tri1.edge(temp)
                    # We don't have to worry about changing gluings on edge (v1+2)%3 of the triangles with label l1
                    # We do have to worry about the following edge:
                    new_glue_e1 = (
                        3 - diagonal_glue_e1 - temp
                    )  # returns the edge which is neither diagonal_glue_e1 nor temp.
                    # This corresponded to the following old edge:
                    old_e1 = (
                        3 - e1 - temp
                    )  # Again this finds the edge which is neither e1 nor temp
                if p2.edge(v2) == tri2.edge(v2):
                    # We don't have to worry about changing gluings on edge v2 of the triangles with label l2
                    # We do have to worry about the following edge:
                    new_glue_e2 = (
                        3 - diagonal_glue_e2 - v2
                    )  # returns the edge which is neither diagonal_glue_e2 nor v2.
                    # This corresponded to the following old edge:
                    old_e2 = 3 - e2 - v2  # Again this finds the edge which is neither e2 nor v2
                else:
                    temp = (v2 + 2) % 3
                    assert p2.edge(temp) == tri2.edge(temp)
                    # We don't have to worry about changing gluings on edge (v2+2)%3 of the triangles with label l2
                    # We do have to worry about the following edge:
                    new_glue_e2 = (
                        3 - diagonal_glue_e2 - temp
                    )  # returns the edge which is neither diagonal_glue_e2 nor temp.
                    # This corresponded to the following old edge:
                    old_e2 = (
                        3 - e2 - temp
                    )  # Again this finds the edge which is neither e2 nor temp

                # remember the old gluings.
                old_opposite1 = s.opposite_edge(l1, old_e1)
                old_opposite2 = s.opposite_edge(l2, old_e2)

                us = s

                # Replace the triangles.
                us.replace_polygon(l1, tri1)
                us.replace_polygon(l2, tri2)
                # Glue along the new diagonal of the quadrilateral
                us.glue((l1, diagonal_glue_e1), (l2, diagonal_glue_e2))
                # Now we deal with that pair of opposite edges of the quadrilateral that need regluing.
                # There are some special cases:
                if old_opposite1 == (l2, old_e2):
                    # These opposite edges were glued to each other.
                    # Do the same in the new surface:
                    us.glue((l1, new_glue_e1), (l2, new_glue_e2))
                else:
                    if old_opposite1 == (l1, old_e1):
                        # That edge was "self-glued".
                        us.glue((l2, new_glue_e2), (l2, new_glue_e2))
                    else:
                        # The edge (l1,old_e1) was glued in a standard way.
                        # That edge now corresponds to (l2,new_glue_e2):
                        us.glue(
                            (l2, new_glue_e2), (old_opposite1[0], old_opposite1[1])
                        )
                    if old_opposite2 == (l2, old_e2):
                        # That edge was "self-glued".
                        us.glue((l1, new_glue_e1), (l1, new_glue_e1))
                    else:
                        # The edge (l2,old_e2) was glued in a standard way.
                        # That edge now corresponds to (l1,new_glue_e1):
                        us.glue(
                            (l1, new_glue_e1), (old_opposite2[0], old_opposite2[1])
                        )
                return s

            def join_polygons(self, p1, e1, test=False, in_place=False):
                r"""
                Join polygons across the provided edge (p1,e1). By default,
                it returns the surface obtained by joining the two polygons
                together. It raises a ValueError if gluing the two polygons
                together results in a non-convex polygon. This is done to the
                current surface if in_place is True, and otherwise a mutable
                copy is made and then modified.

                If test is True then instead of changing the surface, it just
                checks to see if the change would be successful and returns
                True if successful or False if not.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
                    sage: ss = translation_surfaces.ward(3)
                    sage: s = MutableOrientedSimilaritySurface.from_surface(ss)
                    sage: s.join_polygons(0,0, in_place=True)
                    Translation Surface in H_1(0^3) built from an equilateral triangle and a pentagon with 2 marked vertices
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (0, 2*a), (-1, a)])
                    sage: s.join_polygons(0,4, in_place=True)
                    Translation Surface in H_1(0^3) built from a rhombus
                    sage: s.polygon(0)
                    polygon(vertices=[(0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (1, 3*a), (0, 2*a), (-1, a)])

                TESTS::

                    sage: from flatsurf.geometry.categories import TranslationSurfaces
                    sage: s in TranslationSurfaces()
                    True

                """
                poly1 = self.polygon(p1)
                p2, e2 = self.opposite_edge(p1, e1)
                poly2 = self.polygon(p2)
                if p1 == p2:
                    if test:
                        return False
                    else:
                        raise ValueError("Can't glue polygon to itself.")
                t = self.edge_transformation(p2, e2)
                dt = t.derivative()
                vs = []
                edge_map = {}  # Store the pairs for the old edges.
                for i in range(e1):
                    edge_map[len(vs)] = (p1, i)
                    vs.append(poly1.edge(i))
                ne = poly2.num_edges()
                for i in range(1, ne):
                    ee = (e2 + i) % ne
                    edge_map[len(vs)] = (p2, ee)
                    vs.append(dt * poly2.edge(ee))
                for i in range(e1 + 1, poly1.num_edges()):
                    edge_map[len(vs)] = (p1, i)
                    vs.append(poly1.edge(i))

                try:
                    from flatsurf.geometry.polygon import ConvexPolygons
                    new_polygon = ConvexPolygons(self.base_ring())(vs)
                except (ValueError, TypeError):
                    if test:
                        return False
                    else:
                        raise ValueError(
                            "Joining polygons along this edge results in a non-convex polygon."
                        )

                if test:
                    # Gluing would be successful
                    return True

                # Now no longer testing. Do the gluing.
                if in_place:
                    ss = self
                else:
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    ss = MutableOrientedSimilaritySurface.from_surfaces(self)
                s = ss

                inv_edge_map = {}
                for key, value in edge_map.items():
                    inv_edge_map[value] = (p1, key)

                glue_list = []
                for i in range(len(vs)):
                    p3, e3 = edge_map[i]
                    p4, e4 = self.opposite_edge(p3, e3)
                    if p4 == p1 or p4 == p2:
                        glue_list.append(inv_edge_map[(p4, e4)])
                    else:
                        glue_list.append((p4, e4))

                if s.base_label() == p2:
                    s.set_base_label(p1)

                s.remove_polygon(p2)

                s.remove_polygon(p1)
                s.add_polygon(new_polygon, label=p1)
                for e, cross in enumerate(glue_list):
                    s.glue((p1, e), cross)

                return ss

            def subdivide_polygon(self, p, v1, v2, test=False, new_label=None):
                r"""
                Cut the polygon with label p along the diagonal joining vertex
                v1 to vertex v2. This cuts p into two polygons, one will keep the same
                label. The other will get a new label, which can be provided
                via new_label. Otherwise a default new label will be provided.
                If test=False, then the surface will be changed (in place). If
                test=True, then it just checks to see if the change would be successful

                The convention is that the resulting subdivided polygon which has an oriented
                edge going from the original vertex v1 to vertex v2 will keep the label p.
                The other polygon will get a new label.

                The change will be done in place.
                """
                poly = self.polygon(p)
                ne = poly.num_edges()
                if v1 < 0 or v2 < 0 or v1 >= ne or v2 >= ne:
                    if test:
                        return False
                    else:
                        raise ValueError("Provided vertices out of bounds.")
                if abs(v1 - v2) <= 1 or abs(v1 - v2) >= ne - 1:
                    if test:
                        return False
                    else:
                        raise ValueError("Provided diagonal is not actually a diagonal.")

                if v2 < v1:
                    v2 = v2 + ne

                newedges1 = [poly.vertex(v2) - poly.vertex(v1)]
                for i in range(v2, v1 + ne):
                    newedges1.append(poly.edge(i))
                from flatsurf.geometry.polygon import ConvexPolygons
                newpoly1 = ConvexPolygons(self.base_ring())(newedges1)

                newedges2 = [poly.vertex(v1) - poly.vertex(v2)]
                for i in range(v1, v2):
                    newedges2.append(poly.edge(i))
                newpoly2 = ConvexPolygons(self.base_ring())(newedges2)

                # Store the old gluings
                old_gluings = {(p, i): self.opposite_edge(p, i) for i in range(ne)}

                # Update the polygon with label p, add a new polygon.
                self.remove_polygon(p)
                self.add_polygon(newpoly1, label=p)
                if new_label is None:
                    new_label = self.add_polygon(newpoly2)
                else:
                    new_label = self.add_polygon(newpoly2, label=new_label)
                # This gluing is the diagonal we used.
                self.glue((p, 0), (new_label, 0))

                # Setup conversion from old to new labels.
                old_to_new_labels = {}
                for i in range(v1, v2):
                    old_to_new_labels[(p, i % ne)] = (new_label, i - v1 + 1)
                for i in range(v2, ne + v1):
                    old_to_new_labels[(p, i % ne)] = (p, i - v2 + 1)

                for e in range(1, newpoly1.num_edges()):
                    pair = old_gluings[(p, (v2 + e - 1) % ne)]
                    if pair in old_to_new_labels:
                        pair = old_to_new_labels[pair]
                    self.glue((p, e), (pair[0], pair[1]))

                for e in range(1, newpoly2.num_edges()):
                    pair = old_gluings[(p, (v1 + e - 1) % ne)]
                    if pair in old_to_new_labels:
                        pair = old_to_new_labels[pair]
                    self.glue((new_label, e), (pair[0], pair[1]))

            def singularity(self, label, v, limit=None):
                r"""
                Represents the Singularity associated to the v-th vertex of the polygon
                with label ``label``.

                If the surface is infinite, the limit can be set. In this case the
                construction of the singularity is successful if the sequence of
                vertices hit by passing through edges closes up in limit or less steps.

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: s = translation_surfaces.square_torus()
                    sage: pc = s.minimal_cover(cover_type="planar")
                    sage: pc.singularity(pc.base_label(),0)
                    doctest:warning
                    ...
                    UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
                    Vertex 0 of polygon (0, (x, y) |-> (x, y))
                    sage: pc.singularity(pc.base_label(),0,limit=1)
                    Traceback (most recent call last):
                    ...
                    ValueError: number of edges at singularity exceeds limit

                """
                from flatsurf.geometry.surface_objects import Singularity
                return Singularity(self, label, v, limit)

            def point(self, label, point, ring=None, limit=None):
                r"""
                Return a point in this surface.

                INPUT:

                - ``label`` - label of the polygon

                - ``point`` - coordinates of the point inside the polygon

                - ``ring`` (optional) - a ring for the coordinates

                - ``limit`` (optional) - undocumented (only relevant if the point
                  corresponds to a singularity in an infinite surface)

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: s = translation_surfaces.square_torus()
                    sage: pc = s.minimal_cover(cover_type="planar")
                    sage: pc.point(pc.base_label(),(0,0))
                    Vertex 0 of polygon (0, (x, y) |-> (x, y))
                    sage: z = pc.point(pc.base_label(),(sqrt(2)-1,sqrt(3)-1),ring=AA)
                    doctest:warning
                    ...
                    UserWarning: the ring parameter is deprecated and will be removed in a future version of sage-flatsurf; define the surface over a larger ring instead so that this points' coordinates live in the base ring
                    sage: next(iter(z.coordinates(next(iter(z.labels()))))).parent()
                    Vector space of dimension 2 over Algebraic Real Field

                """
                return self(label, point, limit=limit, ring=ring)

            def surface_point(self, *args, **kwargs):
                import warnings
                warnings.warn("surface_point() has been deprecated and will be removed in a future version of sage-flatsurf; use point() instead")

                return self.point(*args, **kwargs)

            def minimal_cover(self, cover_type="translation"):
                r"""
                Return the minimal translation or half-translation cover of the surface.

                Cover type may be either "translation", "half-translation" or "planar".

                The minimal planar cover of a surface S is the smallest cover C so that
                the developing map from the universal cover U to the plane induces a
                well defined map from C to the plane. This is an infinite translation
                surface that is naturally a branched cover of the plane.

                EXAMPLES::

                    sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
                    sage: s = MutableOrientedSimilaritySurface(QQ)
                    sage: square = polygons.square(field=QQ)
                    sage: s.add_polygon(square)
                    0
                    sage: s.glue((0,0), (0,1))
                    sage: s.glue((0,2) ,(0,3))
                    sage: cs = s
                    sage: ts = cs.minimal_cover(cover_type="translation")
                    sage: ts
                    Minimal Translation Cover of Rational Cone Surface built from a square
                    sage: from flatsurf.geometry.categories import TranslationSurfaces
                    sage: ts in TranslationSurfaces()
                    True
                    sage: hts = cs.minimal_cover(cover_type="half-translation")
                    sage: hts
                    Minimal Half-Translation Cover of Rational Cone Surface built from a square
                    sage: from flatsurf.geometry.categories import HalfTranslationSurfaces
                    sage: hts in HalfTranslationSurfaces()
                    True
                    sage: TestSuite(hts).run()
                    sage: ps = cs.minimal_cover(cover_type="planar")
                    sage: ps
                    Minimal Planar Cover of Rational Cone Surface built from a square
                    sage: ps in TranslationSurfaces()
                    True
                    sage: TestSuite(ps).run()

                    sage: from flatsurf import similarity_surfaces
                    sage: S = similarity_surfaces.example()
                    sage: T = S.minimal_cover(cover_type="translation")
                    sage: T
                    Minimal Translation Cover of Genus 1 Surface built from 2 isosceles triangles
                    sage: T in TranslationSurfaces()
                    True
                    sage: T.polygon(T.base_label())
                    polygon(vertices=[(0, 0), (2, -2), (2, 0)])

                """
                if cover_type == "translation":
                    from flatsurf.geometry.minimal_cover import MinimalTranslationCover

                    return MinimalTranslationCover(self)
                if cover_type == "half-translation":
                    from flatsurf.geometry.minimal_cover import MinimalHalfTranslationCover

                    return MinimalHalfTranslationCover(self)
                if cover_type == "planar":
                    from flatsurf.geometry.minimal_cover import MinimalPlanarCover

                    return MinimalPlanarCover(self)
                raise ValueError("Provided cover_type is not supported.")

            def vector_space(self):
                r"""
                Return the vector space in which self naturally embeds.
                """
                from sage.modules.free_module import VectorSpace

                return VectorSpace(self.base_ring(), 2)

            def fundamental_group(self, base_label=None):
                r"""
                Return the fundamental group of this surface.
                """
                if not self.is_finite_type():
                    raise ValueError("the method only work for finite surfaces")
                if base_label is None:
                    base_label = self.base_label()
                from flatsurf.geometry.fundamental_group import FundamentalGroup

                return FundamentalGroup(self, base_label)

            def tangent_bundle(self, ring=None):
                r"""
                Return the tangent bundle

                INPUT:

                - ``ring`` -- an optional field (defaults to the coordinate field of the
                  surface)
                """
                if ring is None:
                    ring = self.base_ring()

                try:
                    return self._tangent_bundle_cache[ring]
                except AttributeError:
                    self._tangent_bundle_cache = {}
                except KeyError:
                    pass

                from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle

                self._tangent_bundle_cache[ring] = SimilaritySurfaceTangentBundle(self, ring)
                return self._tangent_bundle_cache[ring]

            def tangent_vector(self, lab, p, v, ring=None):
                r"""
                Return a tangent vector.

                INPUT:

                - ``lab`` -- label of a polygon

                - ``p`` -- coordinates of a point in the polygon

                - ``v`` -- coordinates of a vector in R^2

                EXAMPLES::

                    sage: from flatsurf.geometry.chamanara import chamanara_surface
                    sage: S = chamanara_surface(1/2)
                    sage: S.tangent_vector(S.base_label(), (1/2,1/2), (1,1))
                    SimilaritySurfaceTangentVector in polygon (1, -1, 0) based at (1/2, -3/2) with vector (1, 1)
                    sage: K.<sqrt2> = QuadraticField(2)
                    sage: S.tangent_vector(S.base_label(), (1/2,1/2), (1,sqrt2), ring=K)
                    SimilaritySurfaceTangentVector in polygon (1, -1, 0) based at (1/2, -3/2) with vector (1, sqrt2)
                """
                from sage.all import vector
                p = vector(p)
                v = vector(v)

                if p.parent().dimension() != 2 or v.parent().dimension() != 2:
                    raise ValueError("p (={!r}) and v (={!v}) should have two coordinates")

                if ring is None:
                    ring = self.base_ring()
                    try:
                        return self.tangent_bundle(ring)(lab, p, v)
                    except TypeError:
                        raise TypeError(
                            "Use the ring=??? option to construct tangent vectors in other field different from the base_ring()."
                        )
                else:
                    return self.tangent_bundle(ring)(lab, p, v)

            def reposition_polygons(self, in_place=False, relabel=None):
                r"""
                We choose a maximal tree in the dual graph of the decomposition into
                polygons, and ensure that the gluings between two polygons joined by
                an edge in this tree is by translation.

                This guarantees that the group generated by the edge identifications is
                minimal among representations of the surface. In particular, if for instance
                you have a translation surface which is anot representable as a translation
                surface (because polygons are presented with rotations) then after this
                change it will be representable as a translation surface.
                """
                if relabel is not None:
                    if relabel:
                        raise NotImplementedError("the relabel keyword has been removed from reposition_polygon; use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead")
                    else:
                        import warnings
                        warnings.warn("the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to reposition_polygons()")
                if not self.is_finite_type():
                    raise NotImplementedError("Only implemented for finite surfaces.")
                if in_place:
                    if not self.is_mutable():
                        raise ValueError(
                            "reposition_polygons in_place is only available "
                            + "for mutable surfaces."
                        )
                    s = self
                else:
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    s = MutableOrientedSimilaritySurface.from_surface(self)
                labels = list(s.labels())
                from flatsurf.geometry.similarity import SimilarityGroup

                S = SimilarityGroup(self.base_ring())
                identity = S.one()
                it = iter(labels)
                label = next(it)
                changes = {label: identity}
                for label in it:
                    polygon = self.polygon(label)
                    adjacencies = {edge: self.opposite_edge(label, edge)[0] for edge in range(polygon.num_edges())}
                    edge = min(adjacencies, key=lambda edge: labels.index(adjacencies[edge]))
                    label2, edge2 = s.opposite_edge(label, edge)
                    changes[label] = changes[label2] * s.edge_transformation(label, edge)
                it = iter(labels)
                # Skip the base label:
                label = next(it)
                for label in it:
                    p = s.polygon(label)
                    p = changes[label].derivative() * p
                    s.replace_polygon(label, p)
                return s

            def triangulation_mapping(self):
                r"""
                Return a ``SurfaceMapping`` triangulating the surface
                or ``None`` if the surface is already triangulated.
                """
                from flatsurf.geometry.mappings import triangulation_mapping

                return triangulation_mapping(self)

            def triangulate(self, in_place=False, label=None, relabel=False):
                r"""
                Return a triangulated version of this surface. (This may be mutable
                or not depending on the input.)

                If label=None (as default) all polygons are triangulated. Otherwise,
                label should be a polygon label. In this case, just this polygon
                is split into triangles.

                This is done in place if in_place is True (defaults to False).

                If we are not doing triangulation in_place, then we must make a copy.
                This can be a relabeled copy (indexed by the non-negative ints)
                or a label preserving copy. The copy is relabeled if relabel=True
                (default False).

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: s=translation_surfaces.mcmullen_L(1,1,1,1)
                    sage: ss=s.triangulate()
                    sage: gs=ss.graphical_surface()
                    sage: gs.make_all_visible()
                    sage: gs
                    Graphical representation of Translation Surface in H_2(2) built from 6 isosceles triangles

                A non-strictly convex example that caused trouble:

                    sage: from flatsurf import *
                    sage: s=similarity_surfaces.self_glued_polygon(polygons(edges=[(1,1),(-3,-1),(1,0),(1,0)]))
                    sage: s=s.triangulate()
                    sage: s.polygon(0).num_edges()
                    3
                """
                if label is None:
                    # We triangulate the whole surface
                    if self.is_finite_type():
                        # Store the current labels.
                        labels = [label for label in self.labels()]
                        if in_place:
                            s = self
                        else:
                            from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                            s = MutableOrientedSimilaritySurface.from_surface(self)
                        # Subdivide each polygon in turn.
                        for label in labels:
                            s = s.triangulate(in_place=True, label=label)
                        return s
                    else:
                        if in_place:
                            raise ValueError(
                                "You can't triangulate an infinite surface in place."
                            )
                        from flatsurf.geometry.delaunay import LazyTriangulatedSurface

                        return LazyTriangulatedSurface(self)
                else:
                    poly = self.polygon(label)
                    n = poly.num_edges()
                    if n > 3:
                        if in_place:
                            s = self
                        else:
                            from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                            s = MutableOrientedSimilaritySurface.from_surface(self)
                    else:
                        # This polygon is already a triangle.
                        return self
                    from flatsurf.geometry.polygon import wedge_product

                    for i in range(n - 3):
                        poly = s.polygon(label)
                        n = poly.num_edges()
                        for i in range(n):
                            e1 = poly.edge(i)
                            e2 = poly.edge((i + 1) % n)
                            if wedge_product(e1, e2) != 0:
                                # This is in case the polygon is a triangle with subdivided edge.
                                e3 = poly.edge((i + 2) % n)
                                if wedge_product(e1 + e2, e3) != 0:
                                    s.subdivide_polygon(label, i, (i + 2) % n)
                                    break
                    return s
                raise RuntimeError("Failed to return anything!")

            def _edge_needs_flip(self, p1, e1):
                p2, e2 = self.opposite_edge(p1, e1)
                poly1 = self.polygon(p1)
                poly2 = self.polygon(p2)
                if poly1.num_edges() != 3 or poly2.num_edges() != 3:
                    raise ValueError("Edge must be adjacent to two triangles.")
                from flatsurf.geometry.matrix_2x2 import similarity_from_vectors

                sim1 = similarity_from_vectors(poly1.edge(e1 + 2), -poly1.edge(e1 + 1))
                sim2 = similarity_from_vectors(poly2.edge(e2 + 2), -poly2.edge(e2 + 1))
                sim = sim1 * sim2
                return sim[1][0] < 0

            def _edge_needs_join(self, p1, e1):
                p2, e2 = self.opposite_edge(p1, e1)
                poly1 = self.polygon(p1)
                poly2 = self.polygon(p2)
                from flatsurf.geometry.matrix_2x2 import similarity_from_vectors

                sim1 = similarity_from_vectors(
                    poly1.vertex(e1) - poly1.vertex(e1 + 2), -poly1.edge(e1 + 1)
                )
                sim2 = similarity_from_vectors(
                    poly2.vertex(e2) - poly2.vertex(e2 + 2), -poly2.edge(e2 + 1)
                )
                sim = sim1 * sim2

                return sim[1][0] == 0

            def delaunay_single_flip(self):
                r"""
                Does a single in place flip of a triangulated mutable surface.
                """
                if not self.is_finite_type():
                    raise NotImplementedError("Not implemented for infinite surfaces.")
                lc = self._label_comparator()
                for (l1, e1), (l2, e2) in self.gluings():
                    if (lc.lt(l1, l2) or (l1 == l2 and e1 <= e2)) and self._edge_needs_flip(
                        l1, e1
                    ):
                        self.triangle_flip(l1, e1, in_place=True)
                        return True
                return False

            def is_delaunay_triangulated(self, limit=None):
                r"""
                Return if the surface is triangulated and the triangulation is Delaunay.
                If limit is set, then it checks this only limit many edges.
                Limit must be set for infinite surfaces.
                """
                if limit is None:
                    if not self.is_finite_type():
                        raise NotImplementedError("A limit must be set for infinite surfaces.")
                    limit = self.num_edges()
                count = 0
                for (l1, e1), (l2, e2) in self.gluings():
                    if count >= limit:
                        break
                    count = count + 1
                    if self.polygon(l1).num_edges() != 3:
                        print("Polygon with label " + str(l1) + " is not a triangle.")
                        return False
                    if self.polygon(l2).num_edges() != 3:
                        print("Polygon with label " + str(l2) + " is not a triangle.")
                        return False
                    if self._edge_needs_flip(l1, e1):
                        print("Edge " + str((l1, e1)) + " needs to be flipped.")
                        print("This edge is glued to " + str((l2, e2)) + ".")
                        return False
                return True

            def is_delaunay_decomposed(self, limit=None):
                r"""
                Return if the decomposition of the surface into polygons is Delaunay.
                If limit is set, then it checks this only limit many polygons.
                Limit must be set for infinite surfaces.
                """
                if limit is None:
                    if not self.is_finite_type():
                        raise NotImplementedError("A limit must be set for infinite surfaces.")
                    limit = len(self.polygons())
                for l1, p1 in zip(self.labels(), self.polygons()):
                    try:
                        c1 = p1.circumscribing_circle()
                    except ValueError:
                        # p1 is not circumscribed
                        return False
                    for e1 in range(p1.num_edges()):
                        c2 = self.edge_transformation(l1, e1) * c1
                        l2, e2 = self.opposite_edge(l1, e1)
                        if c2.point_position(self.polygon(l2).vertex(e2 + 2)) != -1:
                            # The circumscribed circle developed into the adjacent polygon
                            # contains a vertex in its interior or boundary.
                            return False
                    return True

            def delaunay_triangulation(
                self,
                triangulated=False,
                in_place=False,
                limit=None,
                direction=None,
                relabel=None,
            ):
                r"""
                Returns a Delaunay triangulation of a surface, or make some
                triangle flips to get closer to the Delaunay decomposition.

                INPUT:

                - ``triangulated`` (boolean) - If true, the algorithm assumes the
                  surface is already triangulated. It does this without verification.

                - ``in_place`` (boolean) - If true, the triangulating and the
                  triangle flips are done in place.  Otherwise, a mutable copy of the
                  surface is made.

                - ``limit`` (None or Integer) - If None, this will return a
                  Delaunay triangulation. If limit is an integer 1 or larger, then at
                  most limit many diagonal flips will be done.

                - ``direction`` (None or Vector) - with two entries in the base field
                    Used to determine labels when a pair of triangles is flipped. Each triangle
                    has a unique separatrix which points in the provided direction or its
                    negation. As such a vector determines a sign for each triangle.
                    A pair of adjacent triangles have opposite signs. Labels are chosen
                    so that this sign is preserved (as a function of labels).

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: from flatsurf.geometry.delaunay import *

                    sage: m = matrix([[2,1],[1,1]])
                    sage: s = m*translation_surfaces.infinite_staircase()
                    sage: ss = s.delaunay_triangulation()
                    sage: ss.base_label()
                    (0, (0, 1, 2))
                    sage: ss.polygon((0, (0, 1, 2)))
                    polygon(vertices=[(0, 0), (1, 0), (1, 1)])
                    sage: TestSuite(ss).run()
                    sage: ss.is_delaunay_triangulated(limit=10)
                    True
                """
                if relabel is not None:
                    if relabel:
                        raise NotImplementedError("the relabel keyword has been removed from delaunay_triangulation(); use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead")
                    else:
                        import warnings
                        warnings.warn("the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_triangulation()")

                if not self.is_finite_type() and limit is None:
                    if in_place:
                        raise ValueError(
                            "in_place delaunay triangulation is not possible for infinite surfaces unless a limit is set."
                        )
                    if self.is_mutable():
                        raise ValueError(
                            "delaunay_triangulation only works on infinite "
                            + "surfaces if they are immutable or if a limit is set."
                        )
                    from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface

                    return LazyDelaunayTriangulatedSurface(self, direction=direction, category=self.category())
                if in_place and not self.is_mutable():
                    raise ValueError(
                        "in_place delaunay_triangulation only defined for mutable surfaces"
                    )
                if triangulated:
                    if in_place:
                        s = self
                    else:
                        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                        s = MutableOrientedSimilaritySurface.from_surface(self)
                else:
                    if in_place:
                        s = self
                        self.triangulate(in_place=True)
                    else:
                        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                        s = MutableOrientedSimilaritySurface.from_surface(self)
                        s.triangulate(in_place=True)
                loop = True
                if direction is None:
                    base_ring = self.base_ring()
                    direction = self.vector_space()((base_ring.zero(), base_ring.one()))

                if direction.is_zero():
                    raise ValueError

                if s.is_finite_type() and limit is None:
                    from collections import deque

                    unchecked_labels = deque(s.labels())
                    checked_labels = set()
                    while unchecked_labels:
                        label = unchecked_labels.popleft()
                        flipped = False
                        for edge in range(3):
                            if s._edge_needs_flip(label, edge):
                                # Record the current opposite edge:
                                label2, edge2 = s.opposite_edge(label, edge)
                                # Perform the flip.
                                s.triangle_flip(label, edge, in_place=True, direction=direction)
                                # Move the opposite polygon to the list of labels we need to check.
                                if label2 != label:
                                    try:
                                        checked_labels.remove(label2)
                                        unchecked_labels.append(label2)
                                    except KeyError:
                                        # Occurs if label2 is not in checked_labels
                                        pass
                                flipped = True
                                break
                        if flipped:
                            unchecked_labels.append(label)
                        else:
                            checked_labels.add(label)
                    return s
                else:
                    # Old method for infinite surfaces, or limits.
                    count = 0
                    lc = self._label_comparator()
                    while loop:
                        loop = False
                        for (l1, e1), (l2, e2) in s.gluings():
                            if (
                                lc.lt(l1, l2) or (l1 == l2 and e1 <= e2)
                            ) and s._edge_needs_flip(l1, e1):
                                s.triangle_flip(l1, e1, in_place=True, direction=direction)
                                count += 1
                                if limit is not None and count >= limit:
                                    return s
                                loop = True
                                break
                    return s

            def delaunay_single_join(self):
                if not self.is_finite_type():
                    raise NotImplementedError("Not implemented for infinite surfaces.")
                lc = self._label_comparator()
                for (l1, e1), (l2, e2) in self.gluings():
                    if (lc.lt(l1, l2) or (l1 == l2 and e1 <= e2)) and self._edge_needs_join(
                        l1, e1
                    ):
                        self.join_polygons(l1, e1, in_place=True)
                        return True
                return False

            def delaunay_decomposition(
                self,
                triangulated=False,
                delaunay_triangulated=False,
                in_place=False,
                direction=None,
                relabel=None,
            ):
                r"""
                Return the Delaunay Decomposition of this surface.

                INPUT:

                - ``triangulated`` (boolean) - If true, the algorithm assumes the
                  surface is already triangulated. It does this without verification.

                - ``delaunay_triangulated`` (boolean) - If true, the algorithm assumes
                  the surface is already delaunay_triangulated. It does this without
                  verification.

                - ``in_place`` (boolean) - If true, the triangulating and the triangle
                  flips are done in place. Otherwise, a mutable copy of the surface is
                  made.

                - ``direction`` - (None or Vector with two entries in the base field) -
                  Used to determine labels when a pair of triangles is flipped. Each triangle
                  has a unique separatrix which points in the provided direction or its
                  negation. As such a vector determines a sign for each triangle.
                  A pair of adjacent triangles have opposite signs. Labels are chosen
                  so that this sign is preserved (as a function of labels).

                EXAMPLES::

                    sage: from flatsurf import *
                    sage: s0 = translation_surfaces.octagon_and_squares()
                    sage: a = s0.base_ring().gens()[0]
                    sage: m = Matrix([[1,2+a],[0,1]])
                    sage: s = m*s0
                    sage: s = s.triangulate()
                    sage: ss = s.delaunay_decomposition(triangulated=True)
                    sage: len(ss.polygons())
                    3

                    sage: p = polygons((4,0),(-2,1),(-2,-1))
                    sage: s0 = similarity_surfaces.self_glued_polygon(p)
                    sage: s = s0.delaunay_decomposition()
                    sage: TestSuite(s).run()

                    sage: m = matrix([[2,1],[1,1]])
                    sage: s = m*translation_surfaces.infinite_staircase()
                    sage: ss = s.delaunay_decomposition()
                    sage: ss.base_label()
                    (0, (0, 1, 2))
                    sage: ss.polygon(ss.base_label())
                    polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
                    sage: TestSuite(ss).run()
                    sage: ss.is_delaunay_decomposed(limit=10)
                    True
                """
                if relabel is not None:
                    if relabel:
                        raise NotImplementedError("the relabel keyword has been removed from delaunay_decomposition(); use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead")
                    else:
                        import warnings
                        warnings.warn("the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to delaunay_decomposition()")
                if not self.is_finite_type():
                    if in_place:
                        raise ValueError(
                            "in_place delaunay_decomposition is not possible for infinite surfaces."
                        )
                    if self.is_mutable():
                        raise ValueError(
                            "delaunay_decomposition only works on infinite "
                            + "surfaces if they are immutable."
                        )
                    from flatsurf.geometry.delaunay import LazyDelaunaySurface

                    return LazyDelaunaySurface(self, direction=direction, category=self.category())
                if in_place:
                    s = self
                else:
                    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                    s = MutableOrientedSimilaritySurface.from_surface(self)
                if not delaunay_triangulated:
                    s.delaunay_triangulation(
                        triangulated=triangulated, in_place=True, direction=direction
                    )
                # Now s is Delaunay Triangulated
                loop = True
                lc = self._label_comparator()
                while loop:
                    loop = False
                    for (l1, e1), (l2, e2) in s.gluings():
                        if (lc.lt(l1, l2) or (l1 == l2 and e1 <= e2)) and s._edge_needs_join(
                            l1, e1
                        ):
                            s.join_polygons(l1, e1, in_place=True)
                            loop = True
                            break
                return s

            def saddle_connections(
                self,
                squared_length_bound,
                initial_label=None,
                initial_vertex=None,
                sc_list=None,
                check=False,
            ):
                r"""
                Returns a list of saddle connections on the surface whose length squared is less than or equal to squared_length_bound.
                The length of a saddle connection is measured using holonomy from polygon in which the trajectory starts.

                If initial_label and initial_vertex are not provided, we return all saddle connections satisfying the bound condition.

                If initial_label and initial_vertex are provided, it only provides saddle connections emanating from the corresponding
                vertex of a polygon. If only initial_label is provided, the added saddle connections will only emanate from the
                corresponding polygon.

                If sc_list is provided the found saddle connections are appended to this list and the resulting list is returned.

                If check==True it uses the checks in the SaddleConnection class to sanity check our results.

                EXAMPLES::
                    sage: from flatsurf import *
                    sage: s = translation_surfaces.square_torus()
                    sage: sc_list = s.saddle_connections(13, check=True)
                    sage: len(sc_list)
                    32
                """
                if squared_length_bound <= 0:
                    raise ValueError

                if sc_list is None:
                    sc_list = []
                if initial_label is None:
                    if not self.is_finite_type():
                        raise NotImplementedError
                    if initial_vertex is not None:
                        raise ValueError(
                            "when initial_label is not provided, then initial_vertex must not be provided either"
                        )
                    for label in self.labels():
                        self.saddle_connections(
                            squared_length_bound, initial_label=label, sc_list=sc_list
                        )
                    return sc_list
                if initial_vertex is None:
                    for vertex in range(self.polygon(initial_label).num_edges()):
                        self.saddle_connections(
                            squared_length_bound,
                            initial_label=initial_label,
                            initial_vertex=vertex,
                            sc_list=sc_list,
                        )
                    return sc_list

                # Now we have a specified initial_label and initial_vertex
                from flatsurf.geometry.similarity import SimilarityGroup
                SG = SimilarityGroup(self.base_ring())
                start_data = (initial_label, initial_vertex)
                from flatsurf.geometry.circle import Circle
                circle = Circle(
                    self.vector_space().zero(), squared_length_bound, base_ring=self.base_ring()
                )
                p = self.polygon(initial_label)
                v = p.vertex(initial_vertex)
                last_sim = SG(-v[0], -v[1])

                # First check the edge eminating rightward from the start_vertex.
                e = p.edge(initial_vertex)
                if e[0] ** 2 + e[1] ** 2 <= squared_length_bound:
                    from flatsurf.geometry.surface_objects import SaddleConnection
                    sc_list.append(SaddleConnection(self, start_data, e))

                # Represents the bounds of the beam of trajectories we are sending out.
                wedge = (
                    last_sim(p.vertex((initial_vertex + 1) % p.num_edges())),
                    last_sim(p.vertex((initial_vertex + p.num_edges() - 1) % p.num_edges())),
                )

                # This will collect the data we need for a depth first search.
                chain = [
                    (
                        last_sim,
                        initial_label,
                        wedge,
                        [
                            (initial_vertex + p.num_edges() - i) % p.num_edges()
                            for i in range(2, p.num_edges())
                        ],
                    )
                ]

                while len(chain) > 0:
                    # Should verts really be edges?
                    sim, label, wedge, verts = chain[-1]
                    if len(verts) == 0:
                        chain.pop()
                        continue
                    vert = verts.pop()
                    p = self.polygon(label)
                    # First check the vertex
                    vert_position = sim(p.vertex(vert))
                    from flatsurf.geometry.polygon import wedge_product
                    if (
                        wedge_product(wedge[0], vert_position) > 0
                        and wedge_product(vert_position, wedge[1]) > 0
                        and vert_position[0] ** 2 + vert_position[1] ** 2
                        <= squared_length_bound
                    ):
                        sc_list.append(
                            SaddleConnection(
                                self,
                                start_data,
                                vert_position,
                                end_data=(label, vert),
                                end_direction=~sim.derivative() * -vert_position,
                                holonomy=vert_position,
                                end_holonomy=~sim.derivative() * -vert_position,
                                check=check,
                            )
                        )
                    # Now check if we should develop across the edge
                    vert_position2 = sim(p.vertex((vert + 1) % p.num_edges()))
                    if (
                        wedge_product(vert_position, vert_position2) > 0
                        and wedge_product(wedge[0], vert_position2) > 0
                        and wedge_product(vert_position, wedge[1]) > 0
                        and circle.line_segment_position(vert_position, vert_position2) == 1
                    ):
                        if wedge_product(wedge[0], vert_position) > 0:
                            # First in new_wedge should be vert_position
                            if wedge_product(vert_position2, wedge[1]) > 0:
                                new_wedge = (vert_position, vert_position2)
                            else:
                                new_wedge = (vert_position, wedge[1])
                        else:
                            if wedge_product(vert_position2, wedge[1]) > 0:
                                new_wedge = (wedge[0], vert_position2)
                            else:
                                new_wedge = wedge
                        new_label, new_edge = self.opposite_edge(label, vert)
                        new_sim = sim * ~self.edge_transformation(label, vert)
                        p = self.polygon(new_label)
                        chain.append(
                            (
                                new_sim,
                                new_label,
                                new_wedge,
                                [
                                    (new_edge + p.num_edges() - i) % p.num_edges()
                                    for i in range(1, p.num_edges())
                                ],
                            )
                        )
                return sc_list

            def ramified_cover(self, d, data):
                r"""
                Make a ramified cover of this surface.

                INPUT:

                - ``d`` - integer (the degree of the cover)

                - ``data`` - a dictionary which to a pair ``(label, edge_num)`` associates a permutation
                  of {1,...,d}

                EXAMPLES:

                The L-shape origami::

                    sage: import flatsurf
                    sage: T = flatsurf.translation_surfaces.square_torus()
                    sage: T.ramified_cover(3, {(0,0): '(1,2)', (0,1): '(1,3)'})
                    Translation Surface in H_2(2) built from 3 squares
                    sage: O = T.ramified_cover(3, {(0,0): '(1,2)', (0,1): '(1,3)'})
                    sage: O.stratum()
                    H_2(2)

                TESTS::

                    sage: import flatsurf
                    sage: T = flatsurf.translation_surfaces.square_torus()
                    sage: T.ramified_cover(3, {(0,0): '(1,2)', (0,2): '(1,3)'})
                    Traceback (most recent call last):
                    ...
                    ValueError: inconsistent covering data
                """
                from sage.groups.perm_gps.permgroup_named import SymmetricGroup

                G = SymmetricGroup(d)
                for k in data:
                    data[k] = G(data[k])
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                cover = MutableOrientedSimilaritySurface(self.base_ring())
                edges = set(self.edges())
                cover_labels = {}
                for i in range(1, d + 1):
                    for lab in self.labels():
                        cover_labels[(lab, i)] = cover.add_polygon(self.polygon(lab))
                while edges:
                    lab, e = elab = edges.pop()
                    llab, ee = eelab = self.opposite_edge(lab, e)
                    edges.remove(eelab)
                    if elab in data:
                        if eelab in data:
                            if not (data[elab] * data[eelab]).is_one():
                                raise ValueError("inconsistent covering data")
                        s = data[elab]
                    elif eelab in data:
                        s = ~data[eelab]
                    else:
                        s = G.one()

                    for i in range(1, d + 1):
                        p0 = cover_labels[(lab, i)]
                        p1 = cover_labels[(lab, s(i))]
                        cover.glue((p0, e), (p1, ee))
                cover.set_immutable()
                return cover

            def subdivide(self):
                r"""
                Return a copy of this surface whose polygons have been partitioned into
                smaller triangles with
                :meth:`.polygon.ConvexPolygon.subdivide`.

                EXAMPLES:

                A surface consisting of a single triangle::

                    sage: from flatsurf import MutableOrientedSimilaritySurface, polygon

                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(polygon(edges=[(1, 0), (0, 1), (-1, -1)]), label="Δ")
                    'Δ'

                Subdivision of this surface yields a surface with three triangles::

                    sage: T = S.subdivide()
                    sage: T.labels()
                    (('Δ', 0), ('Δ', 1), ('Δ', 2))

                Note that the new labels are old labels plus an index. We verify that
                the triangles are glued correctly::

                    sage: list(T.gluings())
                    [((('Δ', 0), 1), (('Δ', 1), 2)),
                     ((('Δ', 0), 2), (('Δ', 2), 1)),
                     ((('Δ', 1), 1), (('Δ', 2), 2)),
                     ((('Δ', 1), 2), (('Δ', 0), 1)),
                     ((('Δ', 2), 1), (('Δ', 0), 2)),
                     ((('Δ', 2), 2), (('Δ', 1), 1))]

                If we add another polygon to the original surface and glue things, we
                can see how existing gluings are preserved when subdividing::

                    sage: S.add_polygon(polygon(edges=[(1, 0), (0, 1), (-1, 0), (0, -1)]), label='□')
                    '□'

                    sage: S.glue(("Δ", 0), ("□", 2))
                    sage: S.glue(("□", 1), ("□", 3))

                    sage: T = S.subdivide()

                    sage: T.labels()
                    (('Δ', 0), ('□', 2), ('Δ', 1), ('Δ', 2), ('□', 3), ('□', 1), ('□', 0))
                    sage: list(sorted(T.gluings()))
                    [((('Δ', 0), 0), (('□', 2), 0)),
                     ((('Δ', 0), 1), (('Δ', 1), 2)),
                     ((('Δ', 0), 2), (('Δ', 2), 1)),
                     ((('Δ', 1), 1), (('Δ', 2), 2)),
                     ((('Δ', 1), 2), (('Δ', 0), 1)),
                     ((('Δ', 2), 1), (('Δ', 0), 2)),
                     ((('Δ', 2), 2), (('Δ', 1), 1)),
                     ((('□', 0), 1), (('□', 1), 2)),
                     ((('□', 0), 2), (('□', 3), 1)),
                     ((('□', 1), 0), (('□', 3), 0)),
                     ((('□', 1), 1), (('□', 2), 2)),
                     ((('□', 1), 2), (('□', 0), 1)),
                     ((('□', 2), 0), (('Δ', 0), 0)),
                     ((('□', 2), 1), (('□', 3), 2)),
                     ((('□', 2), 2), (('□', 1), 1)),
                     ((('□', 3), 0), (('□', 1), 0)),
                     ((('□', 3), 1), (('□', 0), 2)),
                     ((('□', 3), 2), (('□', 2), 1))]

                """
                labels = list(self.labels())
                polygons = [self.polygon(label) for label in labels]

                subdivisions = [p.subdivide() for p in polygons]

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                surface = MutableOrientedSimilaritySurface(self.base())

                # Add subdivided polygons
                for s, subdivision in enumerate(subdivisions):
                    label = labels[s]
                    for p, polygon in enumerate(subdivision):
                        surface.add_polygon(polygon, label=(label, p))

                surface.set_base_label((self._base_label, 0))

                # Add gluings between subdivided polygons
                for s, subdivision in enumerate(subdivisions):
                    label = labels[s]
                    for p in range(len(subdivision)):
                        surface.glue(
                            ((label, p), 1), ((label, (p + 1) % len(subdivision)), 2)
                        )

                        # Add gluing from original surface
                        opposite = self.opposite_edge(label, p)
                        if opposite is not None:
                            surface.glue(((label, p), 0), (opposite, 0))

                return surface

            def subdivide_edges(self, parts=2):
                r"""
                Return a copy of this surface whose edges have been split into
                ``parts`` equal pieces each.

                INPUT:

                - ``parts`` -- a positive integer (default: 2)

                EXAMPLES:

                A surface consisting of a single triangle::

                    sage: from flatsurf import MutableOrientedSimilaritySurface
                    sage: from flatsurf.geometry.polygon import Polygon, ConvexPolygons

                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: P = ConvexPolygons(QQ)
                    sage: S.add_polygon(P([(1, 0), (0, 1), (-1, -1)]), label="Δ")
                    'Δ'

                Subdividing this triangle yields a triangle with marked points along
                the edges::

                    sage: T = S.subdivide_edges()

                If we add another polygon to the original surface and glue them, we
                can see how existing gluings are preserved when subdividing::

                    sage: S.add_polygon(P([(1, 0), (0, 1), (-1, 0), (0, -1)]), label='□')
                    '□'

                    sage: S.glue(("Δ", 0), ("□", 2))
                    sage: S.glue(("□", 1), ("□", 3))

                    sage: T = S.subdivide_edges()
                    sage: list(sorted(T.gluings()))
                    [(('Δ', 0), ('□', 5)),
                     (('Δ', 1), ('□', 4)),
                     (('□', 2), ('□', 7)),
                     (('□', 3), ('□', 6)),
                     (('□', 4), ('Δ', 1)),
                     (('□', 5), ('Δ', 0)),
                     (('□', 6), ('□', 3)),
                     (('□', 7), ('□', 2))]

                """
                labels = list(self.labels())
                polygons = [self.polygon(label) for label in labels]

                subdivideds = [p.subdivide_edges(parts=parts) for p in polygons]

                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                surface = MutableOrientedSimilaritySurface(self.base())

                # Add subdivided polygons
                for s, subdivided in enumerate(subdivideds):
                    surface.add_polygon(subdivided, label=labels[s])

                surface.set_base_label(self._base_label)

                # Reestablish gluings between polygons
                for label, polygon, subdivided in zip(labels, polygons, subdivideds):
                    for e in range(polygon.num_edges()):
                        opposite = self.opposite_edge(label, e)
                        if opposite is not None:
                            for p in range(parts):
                                surface.glue(
                                    (label, e * parts + p),
                                    (opposite[0], opposite[1] * parts + (parts - p - 1)),
                                )

                return surface

    class Rational(SurfaceCategoryWithAxiom):
        r"""
        The axiom satisfied by similarity surfaces with rational monodromy,
        i.e., a walk around a point leads to a turn by a rational multiple of
        2π.
        """
        class ParentMethods:
            @staticmethod
            def _is_rational_surface(surface, limit=None):
                r"""
                Return whether the gluings of this surface lead to a rational
                surface, i.e., whether the rotational part of all gluings is a
                rational multiple of π.

                This is a helper method for
                :meth:`flatsurf.geometry.categories.similarity_surfaces.ParentMethods.is_rational_surface`.

                INPUT:

                - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
                  the first ``limit`` polygons are checked

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.infinite_staircase()

                    sage: from flatsurf.geometry.categories import SimilaritySurfaces
                    sage: SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(S, limit=8)
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
                        # have overridden this (just returning the identity matrix e.g.)
                        # and we want to deduce the matrix from the attached polygon
                        # edges instead.
                        matrix = SimilaritySurfaces.Oriented.ParentMethods.edge_matrix.f(surface, label, edge)

                        a = AA(matrix[0, 0])
                        b = AA(matrix[1, 0])
                        q = (a**2 + b**2).sqrt()

                        from flatsurf.geometry.matrix_2x2 import is_cosine_sine_of_rational
                        if not is_cosine_sine_of_rational(a / q, b / q):
                            return False

                return True

            def is_rational_surface(self):
                return True

            def _test_rational_surface(self, **options):
                r"""
                Verify that this is a rational similarity surface.

                EXAMPLES::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S._test_rational_surface()

                """
                tester = self._tester(**options)

                limit = None

                if not self.is_finite_type():
                    limit = 32

                tester.assertTrue(SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(self, limit=limit))

    class FiniteType(SurfaceCategoryWithAxiom):
        class ParentMethods:
            def is_cone_surface(self):
                r"""
                Return whether this finite type surface is a cone surface,
                i.e., glued edges can be transformed into each other with
                isometries.

                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_cone_surface()
                    True

                """
                from flatsurf.geometry.categories import ConeSurfaces
                return ConeSurfaces.ParentMethods._is_cone_surface(self)

            def is_dilation_surface(self, positive=False):
                r"""
                Return whether this finite type surface is a dilation surface,
                i.e., whether glued edges can be transformed into each other by
                translation followed by a dilation (multiplication by a
                diagonal matrix.)

                INPUT:

                - ``positive`` -- a boolean (default: ``False``); whether the
                  entries of the diagonal matrix must be positive or are allowed to
                  be negative.

                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_dilation_surface()
                    True
                    sage: S.is_dilation_surface(positive=True)
                    False

                """
                from flatsurf.geometry.categories import DilationSurfaces
                return DilationSurfaces.ParentMethods._is_dilation_surface(self, positive=positive)

            def is_translation_surface(self, positive=True):
                r"""
                Return whether this finite type surface is a translation
                surface, i.e., glued edges can be transformed into each other
                by translations.

                INPUT:

                - ``positive`` -- a boolean (default: ``True``); whether the
                  transformation must be a translation or is allowed to be a
                  half-translation, i.e., a translation followed by a reflection in
                  a point (equivalently, a rotation by π.)

                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_translation_surface()
                    False
                    sage: S.is_translation_surface(False)
                    True

                ::

                    sage: from flatsurf import translation_surfaces
                    sage: S = translation_surfaces.square_torus()
                    sage: S.is_translation_surface()
                    True

                """
                from flatsurf.geometry.categories import TranslationSurfaces
                return TranslationSurfaces.ParentMethods._is_translation_surface(self, positive=positive)

            def is_rational_surface(self):
                r"""
                Return whether this finite type surface is a rational surface,
                i.e., all the rotational part of all gluings is a rational
                multiple of π.

                EXAMPLES::

                    sage: from flatsurf import polygons, similarity_surfaces
                    sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
                    sage: S = similarity_surfaces.self_glued_polygon(P)
                    sage: S.is_rational_surface()
                    True

                """
                return SimilaritySurfaces.Rational.ParentMethods._is_rational_surface(self)

    class SubcategoryMethods:
        def Rational(self):
            return self._with_axiom("Rational")


all_axioms += ("Rational",)
