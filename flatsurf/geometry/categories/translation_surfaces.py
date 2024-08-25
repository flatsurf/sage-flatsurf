r"""
The category of translation surfaces.

This module provides shared functionality for all surfaces in sage-flatsurf
that are built from Euclidean polygons whose glued edges can be transformed
into each other with translations.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

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

from flatsurf.geometry.categories.surface_category import SurfaceCategoryWithAxiom
from flatsurf.geometry.categories.half_translation_surfaces import (
    HalfTranslationSurfaces,
)
from sage.misc.cachefunc import cached_method


class TranslationSurfaces(SurfaceCategoryWithAxiom):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons with
    translations.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import TranslationSurfaces
        sage: TranslationSurfaces()
        Category of translation surfaces

    """
    # The category of translation surfaces is identical to the category of
    # half-translation surfaces with the positive axiom.
    _base_category_class_and_axiom = (HalfTranslationSurfaces, "Positive")

    def extra_super_categories(self):
        r"""
        Return the other categories that a translation surface is automatically
        a member of (apart from being a positive half-translation surface, its
        orientation is compatible with the orientation of the polygons in the
        real plane, so it's "oriented.")

        EXAMPLES::

            sage: from flatsurf.geometry.categories import TranslationSurfaces
            sage: C = TranslationSurfaces()
            sage: C.extra_super_categories()
            (Category of oriented polygonal surfaces,)

        """
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces

        return (PolygonalSurfaces().Oriented(),)

    class ParentMethods:
        r"""
        Provides methods available to all translation surfaces in
        sage-flatsurf.

        If you want to add functionality for such surfaces you most likely want
        to put it here.
        """

        def is_translation_surface(self, positive=True):
            r"""
            Return whether this surface is a translation surface, i.e., return
            ``True``.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.is_translation_surface(positive=True)
                True
                sage: S.is_translation_surface(positive=False)
                True

            """
            return True

        @staticmethod
        def _is_translation_surface(surface, positive=True, limit=None):
            r"""
            Return whether ``surface`` is a translation surface by checking how its
            polygons are glued.

            This is a helper method for
            :meth:`flatsurf.geometry.categories.similarity_surfaces.ParentMethods.is_translation_surface.

            INPUT:

            - ``surface`` -- an oriented similarity surface

            - ``positive`` -- a boolean (default: ``True``); whether the
              transformation must be a translation or is allowed to be a
              half-translation, i.e., a translation followed by a reflection in
              a point (equivalently, a rotation by π.)

            - ``limit`` -- an integer or ``None`` (default: ``None``); if set, only
              the first ``limit`` polygons are checked

            EXAMPLES::

                sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
                sage: S = translation_surfaces.infinite_staircase()

                sage: from flatsurf.geometry.categories import TranslationSurfaces
                sage: TranslationSurfaces.ParentMethods._is_translation_surface(S, limit=8)
                doctest:warning
                ...
                UserWarning: limit has been deprecated as a keyword argument for _is_translation_surface() and will be removed from a future version of sage-flatsurf; ...
                True
                sage: TranslationSurfaces.ParentMethods._is_translation_surface(MutableOrientedSimilaritySurface.from_surface(S, labels=S.labels()[:8]))
                True

            ::

                sage: from flatsurf import Polygon, similarity_surfaces
                sage: P = Polygon(edges=[(2, 0),(-1, 3),(-1, -3)])
                sage: S = similarity_surfaces.self_glued_polygon(P)

                sage: TranslationSurfaces.ParentMethods._is_translation_surface(S)
                False
                sage: TranslationSurfaces.ParentMethods._is_translation_surface(S, positive=False)
                True

            """
            if limit is not None:
                import warnings

                warnings.warn(
                    "limit has been deprecated as a keyword argument for _is_translation_surface() and will be removed from a future version of sage-flatsurf; "
                    "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
                )

            if "Oriented" not in surface.category().axioms():
                raise NotImplementedError(
                    "cannot decide whether a non-oriented surface is a translation surface yet"
                )

            labels = surface.labels()

            if limit is not None:
                from itertools import islice

                labels = islice(labels, limit)

            checked = set()

            for label in labels:
                for edge in range(len(surface.polygon(label).vertices())):
                    cross = surface.opposite_edge(label, edge)

                    if cross is None:
                        continue

                    if cross in checked:
                        continue

                    checked.add((label, edge))

                    # We do not call self.edge_matrix() since the surface might
                    # have overridden this (just returning the identity matrix e.g.)
                    # and we want to deduce the matrix from the attached polygon
                    # edges instead.
                    from flatsurf.geometry.categories import SimilaritySurfaces

                    matrix = SimilaritySurfaces.Oriented.ParentMethods.edge_matrix.f(  # pylint: disable=no-member
                        surface, label, edge
                    )

                    if not matrix.is_diagonal():
                        return False

                    if matrix[0][0] == 1 and matrix[1][1] == 1:
                        continue

                    if matrix[0][0] == -1 and matrix[1][1] == -1:
                        if not positive:
                            continue

                    return False

            return True

        def minimal_translation_cover(self):
            r"""
            Return the minimal cover of this surface that makes this surface a
            translation surface, i.e., return this surface itself.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.minimal_translation_cover() is S
                True

            """
            return self

        def edge_matrix(self, p, e=None):
            r"""
            Returns the 2x2 matrix representing a similarity which when
            applied to the polygon with label `p` makes it so the edge `e`
            can be glued to its opposite edge by translation.

            Since this is a translation surface, this is just the identity.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.edge_matrix(0, 0)
                [1 0]
                [0 1]

            """
            if e is None:
                import warnings

                warnings.warn(
                    "passing only a single tuple argument to edge_matrix() has been deprecated and will be deprecated in a future version of sage-flatsurf; pass the label and edge index as separate arguments instead"
                )
                p, e = p

            if e < 0 or e >= len(self.polygon(p).vertices()):
                raise ValueError("invalid edge index for this polygon")

            from sage.all import identity_matrix

            return identity_matrix(self.base_ring(), 2)

        def _test_translation_surface(self, **options):
            r"""
            Verify that this is a translation surface.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S._test_translation_surface()

            """
            tester = self._tester(**options)

            limit = None

            if not self.is_finite_type():
                from flatsurf import MutableOrientedSimilaritySurface

                self = MutableOrientedSimilaritySurface.from_surface(
                    self, labels=self.labels()[:32]
                )

            tester.assertTrue(
                TranslationSurfaces.ParentMethods._is_translation_surface(
                    self, limit=limit
                )
            )

    class FiniteType(SurfaceCategoryWithAxiom):
        r"""
        The category of translation surfaces built from finitely many polygons.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.octagon_and_squares()
            sage: s.category()
            Category of connected without boundary finite type translation surfaces

        """

        class ParentMethods:
            r"""
            Provides methods available to all translation surfaces built from
            finitely many polygons.

            If you want to add functionality to such surfaces you most likely
            want to put it here.
            """

            # TODO: Make sure this is cached for immutable surfaces.
            def pyflatsurf(self):
                r"""
                Return an isomorphism to a surface backed by libflatsurf.

                EXAMPLES::

                    sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface

                    sage: S = MutableOrientedSimilaritySurface(QQ)
                    sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1)]), label=0)
                    0
                    sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 1), (0, 1)]), label=1)
                    1

                    sage: S.glue((0, 0), (1, 1))
                    sage: S.glue((0, 1), (1, 2))
                    sage: S.glue((0, 2), (1, 0))

                    sage: S.set_immutable()

                    sage: S.pyflatsurf().codomain()  # optional: pyflatsurf
                    FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

                """
                from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

                return Surface_pyflatsurf._from_flatsurf(self)

            def saddle_connections(
                self,
                squared_length_bound=None,
                initial_label=None,
                initial_vertex=None,
                algorithm=None,
            ):
                if squared_length_bound is not None and squared_length_bound < 0:
                    raise ValueError("length bound must be non-negative")

                if algorithm is None:
                    from flatsurf.features import pyflatsurf_feature

                    if (
                        pyflatsurf_feature.is_present()
                        and pyflatsurf_feature.is_saddle_connection_enumeration_functional()
                    ):
                        algorithm = "pyflatsurf"
                    else:
                        algorithm = "generic"

                if algorithm == "pyflatsurf":
                    from flatsurf.features import pyflatsurf_feature

                    if (
                        not flatsurf_feature.is_saddle_connection_enumeration_functional()
                    ):
                        import warnings

                        warnings.warn(
                            "enumerating saddle connections is broken in your version of pyflatsurf, namely saddle connections are enumerated in the wrong order and consequently some might be missing when enumerating with a length bound; upgrade to pyflatsurf >=3.14.1 to resolve this warning."
                        )
                    return self._saddle_connections_pyflatsurf(
                        squared_length_bound, initial_label, initial_vertex
                    )

                from flatsurf.geometry.categories.similarity_surfaces import (
                    SimilaritySurfaces,
                )

                return SimilaritySurfaces.Oriented.ParentMethods.saddle_connections(
                    self, squared_length_bound, initial_label, initial_vertex
                )

            def _saddle_connections_pyflatsurf(
                self, squared_length_bound, initial_label, initial_vertex
            ):
                pyflatsurf_conversion = self.pyflatsurf()

                connections = (
                    pyflatsurf_conversion.codomain()._flat_triangulation.connections()
                )
                connections = connections.byLength()

                if initial_label is not None:
                    raise NotImplementedError
                    if initial_vertex is not None:
                        raise NotImplementedError

                for connection in connections:
                    from flatsurf.geometry.pyflatsurf.saddle_connection import (
                        SaddleConnection_pyflatsurf,
                    )

                    connection = SaddleConnection_pyflatsurf(
                        connection, pyflatsurf_conversion.codomain()
                    )
                    connection = pyflatsurf_conversion.section()(connection)
                    if squared_length_bound is not None:
                        holonomy = connection.holonomy()
                        # TODO: Use dot_product everywhere.
                        if holonomy.dot_product(holonomy) > squared_length_bound:
                            break
                    yield connection

        class WithoutBoundary(SurfaceCategoryWithAxiom):
            r"""
            The category of translation surfaces without boundary built from
            finitely many polygons.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: s = translation_surfaces.octagon_and_squares()
                sage: s.category()
                Category of connected without boundary finite type translation surfaces

            """

            class ParentMethods:
                r"""
                Provides methods available to all translation surfaces without
                boundary that are built from finitely many polygons.

                If you want to add functionality for such surfaces you most likely
                want to put it here.
                """

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

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: s = translation_surfaces.octagon_and_squares()
                        sage: s.canonicalize_mapping()
                        doctest:warning
                        ...
                        UserWarning: canonicalize_mapping() has been deprecated and will be removed in a future version of sage-flatsurf; use canonicalize() instead
                        Composite morphism:
                          From: Translation Surface in H_3(4) built from 2 squares and a regular octagon
                          To:   Translation Surface in H_3(4) built from 2 squares and a regular octagon
                          Defn:   Delaunay Decomposition morphism:
                                  From: Translation Surface in H_3(4) built from 2 squares and a regular octagon
                                  To:   Translation Surface in H_3(4) built from 2 squares and a regular octagon
                                then
                                  Polygon Standardization morphism:
                                  From: Translation Surface in H_3(4) built from 2 squares and a regular octagon
                                  To:   Translation Surface in H_3(4) built from 2 squares and a regular octagon
                                then
                                  Relabeling morphism:
                                  From: Translation Surface in H_3(4) built from 2 squares and a regular octagon
                                  To:   Translation Surface in H_3(4) built from 2 squares and a regular octagon

                    """
                    import warnings

                    warnings.warn(
                        "canonicalize_mapping() has been deprecated and will be removed in a future version of sage-flatsurf; use canonicalize() instead"
                    )

                    return self.canonicalize()

                def canonicalize(self, in_place=None):
                    r"""
                    Return a canonical version of this translation surface.

                    EXAMPLES:

                    We will check if an element lies in the Veech group::

                        sage: from flatsurf import translation_surfaces
                        sage: s = translation_surfaces.octagon_and_squares()
                        sage: s
                        Translation Surface in H_3(4) built from 2 squares and a regular octagon
                        sage: from flatsurf.geometry.categories import TranslationSurfaces
                        sage: s in TranslationSurfaces()
                        True
                        sage: a = s.base_ring().gen()
                        sage: mat = matrix([[1, 2 + a], [0, 1]])
                        sage: s1 = s.canonicalize().codomain()
                        sage: s1.set_immutable()
                        sage: s2 = (mat * s).canonicalize().codomain()
                        sage: s2.set_immutable()
                        sage: s1.cmp(s2) == 0
                        True
                        sage: hash(s1) == hash(s2)
                        True

                    """
                    if in_place is not None:
                        if in_place:
                            raise NotImplementedError(
                                "calling canonicalize(in_place=True) is not supported anymore"
                            )

                        import warnings

                        warnings.warn(
                            "the in_place keyword of canonicalize() has been deprecated and will be removed in a future version of sage-flatsurf"
                        )

                    delaunay_decomposition = self.delaunay_decomposition()

                    standardization = (
                        delaunay_decomposition.codomain().standardize_polygons()
                    )

                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    s = MutableOrientedSimilaritySurface.from_surface(
                        standardization.codomain()
                    )

                    ss = MutableOrientedSimilaritySurface.from_surface(
                        standardization.codomain()
                    )

                    for label in ss.labels():
                        ss.set_roots([label])
                        if ss.cmp(s) > 0:
                            s.set_roots([label])

                    # We have determined the root label such that this surface
                    # is minimal. Now we relabel all the polygons so that they
                    # are natural numbers in the order of the walk on the
                    # surface.
                    relabeling = s.relabel(
                        {label: i for (i, label) in enumerate(s.labels())},
                        in_place=True,
                    )
                    s.set_immutable()
                    return s

                def j_invariant(self):
                    r"""
                    Return the Kenyon-Smillie J-invariant of this translation surface.

                    It is assumed that the coordinates are defined over a number field.

                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: O = translation_surfaces.regular_octagon()
                        sage: O.j_invariant()
                        (
                                  [2 2]
                        (0), (0), [2 1]
                        )
                    """
                    it = iter(self.labels())
                    lab = next(it)
                    P = self.polygon(lab)
                    Jxx, Jyy, Jxy = P.j_invariant()
                    for lab in it:
                        xx, yy, xy = self.polygon(lab).j_invariant()
                        Jxx += xx
                        Jyy += yy
                        Jxy += xy
                    return (Jxx, Jyy, Jxy)

                @cached_method
                def erase_marked_points(self):
                    r"""
                    Return an isomorphism to a surface with a minimal regular vertices
                    of angle 2π.

                    ALGORITHM:

                    We use the erasure of marked points implemented in pyflatsurf. For
                    this we triangulate, then we Delaunay triangulate, then erase
                    marked points with pyflatsurf, and then Delaunay triangulate again.

                    EXAMPLES::

                        sage: import flatsurf

                        sage: G = SymmetricGroup(4)
                        sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
                        sage: S.stratum()
                        H_2(2, 0)
                        sage: S.erase_marked_points().codomain().stratum() # optional: pyflatsurf  # long time (1s)
                        H_2(2)

                        sage: for (a,b,c) in [(1,4,11), (1,4,15), (3,4,13)]: # long time (10s), optional: pyflatsurf
                        ....:     T = flatsurf.polygons.triangle(a,b,c)
                        ....:     S = flatsurf.similarity_surfaces.billiard(T)
                        ....:     S = S.minimal_cover("translation")
                        ....:     print(S.erase_marked_points().codomain().stratum())
                        H_6(10)
                        H_6(2^5)
                        H_8(12, 2)

                    If the surface had no marked points then it is returned unchanged by this
                    function::

                        sage: O = flatsurf.translation_surfaces.regular_octagon()
                        sage: O.erase_marked_points().codomain() is O
                        True

                    This method produces a morphism from the surface with marked points
                    to the surface without marked points::

                        sage: G = SymmetricGroup(4)
                        sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
                        sage: erasure = S.erase_marked_points()  # optional: pyflatsurf  # long time (1s)
                        sage: marked_point = S(1, 1); marked_point  # optional: pyflatsurf  # long time (from above)
                        Vertex 0 of polygon 2
                        sage: erasure(marked_point)  # optional: pyflatsurf  # long time (from above)
                        Point (-1, 0) of polygon (-3, -7, -2)
                        sage: unmarked_point = S(1, 0); unmarked_point  # optional: pyflatsurf  # long time (from above)
                        Vertex 0 of polygon 1
                        sage: erasure(unmarked_point)  # optional: pyflatsurf  # long time (from above)
                        Vertex 0 of polygon (-3, -7, -2)

                    TESTS:

                    Verify that https://github.com/flatsurf/flatsurf/issues/263 has been resolved::

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(angles=(10, 8, 3, 1, 1, 1), lengths=(1, 1, 2, 4))
                        sage: B = similarity_surfaces.billiard(P)
                        sage: S = B.minimal_cover(cover_type="translation")
                        sage: S = S.erase_marked_points().codomain() # long time (3s), optional: pyflatsurf

                    ::

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(angles=(10, 7, 2, 2, 2, 1), lengths=(1, 1, 2, 3))
                        sage: B = similarity_surfaces.billiard(P)
                        sage: S_mp = B.minimal_cover(cover_type="translation")
                        sage: S = S_mp.erase_marked_points().codomain() # long time (3s), optional: pyflatsurf

                    """
                    if all(a != 1 for a in self.angles()):
                        # no 2π angle
                        from flatsurf.geometry.morphism import IdentityMorphism

                        return IdentityMorphism._create_morphism(self)

                    # Triangulate the surface: to_pyflatsurf maps self to a
                    # triangulated libflatsurf surface (later called delaunay0_domain.)
                    to_pyflatsurf = self.pyflatsurf()

                    # Delaunay triangulate: delaunay0 maps delaunay0_domain to delaunay0_codomain.
                    # Since the flips of delaunay() are performed in-place, we create
                    # the mapping using the Tracked[Deformation] feature of pyflatsurf.
                    delaunay0_codomain = (
                        to_pyflatsurf.codomain()._flat_triangulation.clone()
                    )

                    from pyflatsurf import flatsurf

                    delaunay0 = flatsurf.Tracked(
                        delaunay0_codomain.combinatorial(),
                        flatsurf.Deformation[type(delaunay0_codomain)](
                            delaunay0_codomain.clone()
                        ),
                    )

                    delaunay0_codomain.delaunay()

                    # Erase marked points: elimination maps delaunay0_codomain to a
                    # surface without marked points, later called delaunay1_domain.
                    elimination = delaunay0_codomain.eliminateMarkedPoints()

                    # Delaunay triangulate again: delaunay1 maps delaunay1_domain to delaunay1_codomain.
                    # Again, we use a Tracked[Deformation] to create this mapping.
                    delaunay1_codomain = elimination.codomain().clone()

                    delaunay1 = flatsurf.Tracked(
                        delaunay1_codomain.combinatorial(),
                        flatsurf.Deformation[type(delaunay1_codomain)](
                            delaunay1_codomain.clone()
                        ),
                    )

                    delaunay1_codomain.delaunay()

                    # Bring the surface back into sage-flatsurf: from_pyflatsurf maps a
                    # sage-flatsurf surface to delaunay1_codomain.
                    from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

                    codomain_pyflatsurf = Surface_pyflatsurf(delaunay1_codomain)

                    from flatsurf.geometry.pyflatsurf.morphism import (
                        Morphism_from_Deformation,
                    )

                    pyflatsurf_morphism = Morphism_from_Deformation._create_morphism(
                        to_pyflatsurf.codomain(),
                        codomain_pyflatsurf,
                        delaunay1.value() * elimination * delaunay0.value(),
                    )

                    from flatsurf.geometry.pyflatsurf_conversion import (
                        FlatTriangulationConversion,
                    )

                    from_pyflatsurf = FlatTriangulationConversion.from_pyflatsurf(
                        delaunay1_codomain
                    )

                    from flatsurf.geometry.pyflatsurf.morphism import (
                        Morphism_from_pyflatsurf,
                    )

                    from_pyflatsurf = Morphism_from_pyflatsurf._create_morphism(
                        codomain_pyflatsurf, from_pyflatsurf.domain(), from_pyflatsurf
                    )

                    return from_pyflatsurf * pyflatsurf_morphism * to_pyflatsurf

                def rel_deformation(self, deformation, local=None, limit=None):
                    r"""
                    Return a deformed surface obtained by shifting the vertices
                    by ``deformation``.

                    INPUT:

                    - ``deformation`` -- a dict which maps the vertices of this
                      surfaces to vectors. The rel deformation will move each
                      vertex by that amount (relative to the others); any
                      vertex not present in the dict will be treated as a
                      deformation by the zero vector.


                    EXAMPLES::

                        sage: from flatsurf import translation_surfaces
                        sage: S = translation_surfaces.arnoux_yoccoz(4)
                        sage: S1 = S.rel_deformation({S(0, 0): (1, 0)}).canonicalize()  # optional: pyflatsurf

                        sage: a = S.base_ring().gen()
                        sage: S2 = S.rel_deformation({S(0, 0): (a, 0)}).canonicalize()  # optional: pyflatsurf

                        sage: M = matrix([[a, 0], [0, ~a]])
                        sage: S2.cmp((M*S1).canonicalize())  # optional: pyflatsurf
                        0

                    """
                    if local is not None:
                        import warnings

                        warnings.warn(
                            "the local keyword has been removed from rel_deformation() without a replacement; do not use it anymore"
                        )

                    if limit is not None:
                        import warnings

                        warnings.warn(
                            "the limit keyword has been removed from rel_deformation() without a replacement; do not use it anymore"
                        )

                    if not self.is_triangulated():
                        raise NotImplementedError(
                            "only triangulated surfaces can be rel-deformed"
                        )

                    from flatsurf.geometry.pyflatsurf_conversion import (
                        FlatTriangulationConversion,
                    )

                    conversion = FlatTriangulationConversion.to_pyflatsurf(self)

                    from sage.all import vector

                    deformation = {
                        vertex: vector(
                            self.base_ring(), deformation.get(vertex, (0, 0))
                        )
                        for vertex in self.vertices()
                    }

                    # pylint: disable=not-callable
                    deformation = {
                        edge: deformation[self(*self.opposite_edge(*edge))]
                        - deformation[self(*edge)]
                        for edge in self.edges()
                    }
                    # pylint: enable=not-callable

                    vector_space_conversion = conversion.vector_space_conversion()
                    deformation = [
                        vector_space_conversion(
                            deformation[conversion.section(edge.positive())]
                        )
                        for edge in conversion.codomain().edges()
                    ]

                    deformed = (conversion.codomain() + deformation).codomain()

                    return FlatTriangulationConversion.from_pyflatsurf(
                        deformed
                    ).domain()
