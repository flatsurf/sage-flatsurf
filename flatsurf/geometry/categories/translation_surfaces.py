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
from flatsurf.cache import cached_surface_method


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
              a point (equivalently, a rotation by π).

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

        def canonicalize_mapping(self):
            r"""
            Return a SurfaceMapping canonicalizing this translation surface.
            """
            from flatsurf.geometry.mappings import (
                canonicalize_translation_surface_mapping,
            )

            return canonicalize_translation_surface_mapping(self)

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

            @cached_surface_method
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

                    sage: T = S.pyflatsurf().codomain()  # optional: pyflatsurf  # random output due to cppyy deprecation warnings
                    sage: T  # optional: pyflatsurf
                    Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

                """
                from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

                return Surface_pyflatsurf._from_flatsurf(self)

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
                    from surface_dynamics import Stratum
                    from sage.rings.integer_ring import ZZ

                    return Stratum([ZZ(a - 1) for a in self.angles()], 1)

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
                    if in_place is not None:
                        if in_place:
                            raise NotImplementedError(
                                "calling canonicalize(in_place=True) is not supported anymore"
                            )

                        import warnings

                        warnings.warn(
                            "the in_place keyword of canonicalize() has been deprecated and will be removed in a future version of sage-flatsurf"
                        )

                    s = self.delaunay_decompose().codomain().standardize_polygons()

                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    s = MutableOrientedSimilaritySurface.from_surface(s)

                    from flatsurf.geometry.surface import (
                        MutableOrientedSimilaritySurface,
                    )

                    ss = MutableOrientedSimilaritySurface.from_surface(s)

                    for label in ss.labels():
                        ss.set_roots([label])
                        if ss.cmp(s) > 0:
                            s.set_roots([label])

                    # We have chosen the root label such that this surface is minimal.
                    # Now we relabel all the polygons so that they are natural
                    # numbers in the order of the walk on the surface.
                    labels = {label: i for (i, label) in enumerate(s.labels())}
                    s.relabel(labels, in_place=True)
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

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(angles=(10, 8, 3, 1, 1, 1), lengths=(1, 1, 2, 4))
                        sage: B = similarity_surfaces.billiard(P)
                        sage: S = B.minimal_cover(cover_type="translation")
                        sage: S = S.erase_marked_points() # long time (3s), optional: pyflatsurf

                    ::

                        sage: from flatsurf import Polygon, similarity_surfaces
                        sage: P = Polygon(angles=(10, 7, 2, 2, 2, 1), lengths=(1, 1, 2, 3))
                        sage: B = similarity_surfaces.billiard(P)
                        sage: S_mp = B.minimal_cover(cover_type="translation")
                        sage: S = S_mp.erase_marked_points() # long time (3s), optional: pyflatsurf

                    """
                    if all(a != 1 for a in self.angles()):
                        # no 2π angle
                        return self
                    from flatsurf.geometry.pyflatsurf.conversion import (
                        from_pyflatsurf,
                        to_pyflatsurf,
                    )

                    S = to_pyflatsurf(self)
                    S.delaunay()
                    S = S.eliminateMarkedPoints().surface()
                    S.delaunay()
                    return from_pyflatsurf(S)

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

                    from flatsurf.geometry.pyflatsurf.conversion import (
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
