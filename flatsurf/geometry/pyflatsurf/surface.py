r"""
Surfaces backed by pyflatsurf.

There should be no need to create such surfaces directly, not even for authors
of sage-flatsurf. Instead just call
:meth:`~flatsurf.geometry.categories.translation_surfaces.TranslationSurfaces.FiniteType.ParentMethods.pyflatsurf`
on a surface which returns such a surface (or rather, a morphism to such a
surface).

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.square_torus()

    sage: T = S.pyflatsurf().codomain()  # optional: pyflatsurf  # random output due to cppyy deprecation warnings
    sage: T  # optional: pyflatsurf
    Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

Ideally, there should be no need to use the underlying ``FlatTriangulation``
directly, but it can be accessed with
:meth:`Surface_pyflatsurf.flat_triangulation`::

    sage: T.flat_triangulation()  # optional: pyflatsurf
    FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2024 Julian Rüth
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

from sage.misc.cachefunc import cached_method


class Surface_pyflatsurf(OrientedSimilaritySurface):
    r"""
    A translation surface backed by pyflatsurf.

    Most surfaces in sage-flatsurf, such as the
    :class:`~flatsurf.geometry.surface.MutableOrientedSimilaritySurface` are
    implemented in Python instead.

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

        sage: T = S.pyflatsurf().codomain()  # optional: pyflatsurf

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

        sage: isinstance(T, Surface_pyflatsurf)  # optional: pyflatsurf
        True
        sage: TestSuite(T).run()  # optional: pyflatsurf

    Disconnected surfaces can be import from pyflatsurf::

        sage: S = MutableOrientedSimilaritySurface(QQ)
        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
        0
        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]))
        1

        sage: S.glue((0, 0), (0, 2))
        sage: S.glue((0, 1), (0, 3))
        sage: S.glue((1, 0), (1, 2))
        sage: S.glue((1, 1), (1, 3))

        sage: S.set_immutable()

        sage: T = S.pyflatsurf().codomain()  # optional: pyflatsurf

        sage: T.is_connected()  # optional: pyflatsurf
        False

        sage: T.flat_triangulation()  # optional: pyflatsurf
        FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2)(4, -6, 5, -4, 6, -5), faces = (1, 2, 3)(-1, -2, -3)(4, 5, 6)(-4, -5, -6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1), 4: (1, 0), 5: (0, 1), 6: (-1, -1)}

        sage: TestSuite(T).run()  # optional: pyflatsurf

    """

    def __init__(self, flat_triangulation):
        self._flat_triangulation = flat_triangulation

        from flatsurf.geometry.pyflatsurf.conversion import RingConversion

        self._ring_conversion = RingConversion.from_pyflatsurf_from_flat_triangulation(
            flat_triangulation
        )

        from flatsurf.geometry.pyflatsurf.conversion import VectorSpaceConversion

        self._vector_space_conversion = VectorSpaceConversion.to_pyflatsurf(
            self._ring_conversion.domain() ** 2, ring_conversion=self._ring_conversion
        )

        from flatsurf.geometry.categories import TranslationSurfaces

        category = TranslationSurfaces().FiniteType()
        if flat_triangulation.hasBoundary():
            category = category.WithBoundary()
        else:
            category = category.WithoutBoundary()

        super().__init__(base=self._ring_conversion.domain(), category=category)

        if self.is_connected():
            self._refine_category_(category.Connected())

    def flat_triangulation(self):
        r"""
        Return the pyflatsurf ``FlatTriangulation`` object underlying this
        surface.

        .. WARNING::

            This surface is supposed to be immutable. However, the
            ``FlatTriangulation`` is not immutable. If you make any changes to
            it, things might break in strange ways. If you want to make
            modifications to the ``FlatTriangulation`` be sure to work on a
            ``clone()`` of it instead.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.pyflatsurf().codomain().flat_triangulation()  # optional: pyflatsurf
            FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

        """
        return self._flat_triangulation

    def is_mutable(self):
        r"""
        Return whether this surface is mutable.

        .. WARNING::

            This surface is supposed to be immutable. However, the underlying
            ``FlatTriangulation`` is not immutable. If you make any changes to
            it, things might break in strange ways. If you want to make
            modifications to the ``FlatTriangulation`` be sure to work on a
            ``clone()`` of it instead. See :meth:`flat_triangulation`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.pyflatsurf().codomain().is_mutable()  # optional: pyflatsurf
            False

        """
        return False

    @cached_method
    def roots(self):
        r"""
        Return some root labels for this surface, one for each connected
        component of the surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.pyflatsurf().codomain()  # optional: pyflatsurf

            sage: list(T.labels())  # optional: pyflatsurf
            [(1, 2, 3), (-3, -1, -2)]

            sage: T.roots()  # optional: pyflatsurf
            ((1, 2, 3),)

        """
        roots = []
        seen = set()

        for face in self.flat_triangulation().faces():
            root = self._normalize_label(face)
            if root not in seen:
                roots.append(root)
                for root in self.component(root):
                    seen.add(root)

        return tuple(roots)

    def pyflatsurf(self):
        r"""
        Return a version of this surface that is backed by pyflatsurf. Since
        this surface is already backed by pyflatsurf, this returns a trivial
        morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.pyflatsurf().codomain()  # optional: pyflatsurf

            sage: T.pyflatsurf()  # optional: pyflatsurf
            Identity endomorphism of Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

        """
        from flatsurf.geometry.morphism import IdentityMorphism

        return IdentityMorphism._create_morphism(self)

    @classmethod
    def _from_flatsurf(cls, surface):
        r"""
        Return an isomorphism to a :class:`Surface_pyflatsurf` built from
        ``surface``, i.e., represent ``surface`` in pyflatsurf wrapped as a
        :class:`Surface` for sage-flatsurf.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: Surface_pyflatsurf._from_flatsurf(S)  # optional: pyflatsurf
            pyflatsurf conversion morphism:
              From: Triangulation of Translation Surface in H_1(0) built from a square
              To:   Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

        """
        if isinstance(surface, Surface_pyflatsurf):
            return surface

        if not surface.is_triangulated():
            triangulation = surface.triangulate()
            to_pyflatsurf = cls._from_flatsurf(triangulation.codomain())
            return to_pyflatsurf * triangulation

        from flatsurf.geometry.pyflatsurf.conversion import FlatTriangulationConversion

        to_pyflatsurf = FlatTriangulationConversion.to_pyflatsurf(surface)

        surface_pyflatsurf = Surface_pyflatsurf(to_pyflatsurf.codomain())

        from flatsurf.geometry.pyflatsurf.morphism import Morphism_to_pyflatsurf

        return Morphism_to_pyflatsurf._create_morphism(
            surface, surface_pyflatsurf, to_pyflatsurf
        )

    def __repr__(self):
        r"""
        Return a printable representation of this surface, namely, print this
        surface essentially as pyflatsurf would.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: S.pyflatsurf().codomain()  # optional: pyflatsurf
            Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

        """
        return f"Surface backed by {self._flat_triangulation!r}"

    def apply_matrix(self, m, in_place=None):
        r"""
        Overrides the generic
        :meth:`~flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.ParentMethods.apply_matrix` for this
        pyflatsurf backed surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().pyflatsurf().codomain()  # optional: pyflatsurf

            sage: S.apply_matrix(matrix([[1, 2], [0, 1]])).codomain()  # optional: pyflatsurf
            Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (2, 1), 3: (-3, -1)}

        """
        if in_place:
            return super().apply_matrix(m, in_place=in_place)

        from sage.all import matrix

        m = matrix(m, base_ring=self.base_ring())
        m = [self._ring_conversion(x) for x in m.list()]

        deformation = self._flat_triangulation.applyMatrix(*m)
        codomain = Surface_pyflatsurf(deformation.codomain())

        from flatsurf.geometry.pyflatsurf.morphism import Morphism_Deformation

        return Morphism_Deformation._create_morphism(self, codomain, deformation)

    @classmethod
    def _normalize_label(cls, label):
        r"""
        Return a normalized version of the ``label`` of this surface.

        Labels are formed by the triple of the attached half edges in order.
        They are normalized such that the minimal half edge is first.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: Surface_pyflatsurf._normalize_label((1, 2, 3))
            (1, 2, 3)
            sage: Surface_pyflatsurf._normalize_label((1, -2, 3))
            (-2, 3, 1)
            sage: Surface_pyflatsurf._normalize_label((1, -2, -3))
            (-3, 1, -2)

        """
        label = tuple(edge.id() if hasattr(edge, "id") else int(edge) for edge in label)

        shift = label.index(min(*label))

        label = label[shift:] + label[:shift]
        return label

    def polygon(self, label):
        r"""
        Return the polygon with ``label`` in this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()
            sage: S = S.pyflatsurf().codomain()  # optional: pyflatsurf

            sage: S.polygon((1, 2, 3))  # optional: pyflatsurf
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])

        """
        label = Surface_pyflatsurf._normalize_label(label)

        from pyflatsurf import flatsurf

        half_edges = (flatsurf.HalfEdge(half_edge) for half_edge in label)
        vectors = [
            self._flat_triangulation.fromHalfEdge(half_edge) for half_edge in half_edges
        ]
        vectors = [self._vector_space_conversion.section(vector) for vector in vectors]

        from flatsurf.geometry.polygon import Polygon

        return Polygon(edges=vectors)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon and edge that is across from ``edge`` of ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()
            sage: S = S.pyflatsurf().codomain()  # optional: pyflatsurf

            sage: S.opposite_edge((1, 2, 3), 0)  # optional: pyflatsurf
            ((-3, -1, -2), 1)
            sage: S.opposite_edge((1, 2, 3), 1)  # optional: pyflatsurf
            ((-3, -1, -2), 2)
            sage: S.opposite_edge((1, 2, 3), 2)  # optional: pyflatsurf
            ((-3, -1, -2), 0)

        """
        label = Surface_pyflatsurf._normalize_label(label)

        from pyflatsurf import flatsurf

        half_edge = flatsurf.HalfEdge(label[edge])

        opposite_half_edge = -half_edge

        opposite_label = self._flat_triangulation.face(opposite_half_edge)
        opposite_label = Surface_pyflatsurf._normalize_label(opposite_label)

        return opposite_label, opposite_label.index(opposite_half_edge.id())

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().pyflatsurf().codomain()  # optional: pyflatsurf
            sage: T = translation_surfaces.square_torus().pyflatsurf().codomain()  # optional: pyflatsurf
            sage: S == T  # optional: pyflatsurf
            True

        """
        if not isinstance(other, Surface_pyflatsurf):
            return False

        return self._flat_triangulation == other._flat_triangulation

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().pyflatsurf().codomain()  # optional: pyflatsurf
            sage: T = translation_surfaces.square_torus().pyflatsurf().codomain()  # optional: pyflatsurf
            sage: hash(S) == hash(T)  # optional: pyflatsurf
            True

        """
        return hash(self._flat_triangulation)
