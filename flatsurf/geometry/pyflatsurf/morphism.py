r"""
Morphisms involving pylatsurf backed surfaces.

This module extends :mod:`flatsurf.geometry.morphism` with morphisms that rely
on the C++/Python library ``pyflatsurf``.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.veech_double_n_gon(5)
    sage: to_pyflatsurf = S.pyflatsurf()  # optional: pyflatsurf  # random output due to cppyy deprecation warnings
    sage: to_pyflatsurf  # optional: pyflatsurf
    Composite morphism:
      From: Translation Surface in H_2(2) built from 2 regular pentagons
      To:   Surface backed by FlatTriangulationCombinatorial(...) with vectors ...
      Defn: Triangulation morphism:
              ...
            then pyflatsurf conversion morphism:
              ...

    sage: to_pyflatsurf.codomain().flat_triangulation()  # optional: pyflatsurf
    FlatTriangulationCombinatorial(...) with vectors ...

.. SEEALSO:

    :mod:`flatsurf.geometry.pyflatsurf.conversion` for much of the underlying
    machinery.

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
from flatsurf.geometry.morphism import SurfaceMorphism


class Morphism_to_pyflatsurf(SurfaceMorphism):
    r"""
    A trivial isomorphism from a sage-flatsurf translation surface to a
    pyflatsurf backed translation surface.

    You should not create such morphisms directly but rely on the caching
    provided by
    :meth:`~flatsurf.geometry.categories.translation_surfaces.TranslationSurfaces.FiniteType.ParentMethods.pyflatsurf`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
        sage: to_pyflatsurf = S.pyflatsurf()  # optional: pyflatsurf

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.morphism import Morphism_to_pyflatsurf
        sage: isinstance(to_pyflatsurf, Morphism_to_pyflatsurf)  # optional: pyflatsurf
        True

        sage: TestSuite(to_pyflatsurf).run()  # optional: pyflatsurf

    """

    def __init__(self, parent, pyflatsurf_conversion):
        super().__init__(parent)
        self._pyflatsurf_conversion = pyflatsurf_conversion

    def _image_edge(self, label, edge):
        r"""
        Helper method for :meth:`_image_homology_edge` and others.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: to_pyflatsurf = S.pyflatsurf()  # optional: pyflatsurf
            sage: to_pyflatsurf._image_edge((0, 0), 0)  # optional: pyflatsurf
            ((1, 2, 3), 0)

        """
        half_edge = self._pyflatsurf_conversion((label, edge))
        face = tuple(self.codomain().flat_triangulation().face(half_edge))
        label = type(self.codomain())._normalize_label(face)
        edge = label.index(half_edge)
        return (label, edge)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: to_pyflatsurf = S.pyflatsurf()  # optional: pyflatsurf
            sage: to_pyflatsurf._image_homology_edge((0, 0), 0, codomain=to_pyflatsurf.codomain().homology())  # optional: pyflatsurf
            B[((1, 2, 3), 0)]

        """
        return codomain(self._image_edge(label, edge))

    def section(self):
        r"""
        Return the inverse of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: to_pyflatsurf = S.pyflatsurf()  # optional: pyflatsurf
            sage: to_pyflatsurf.section()  # optional: pyflatsurf
            pyflatsurf reconversion morphism:
              From: Surface backed by FlatTriangulationCombinatorial(...) with vectors ...
              To:   Triangulation of Translation Surface in H_2(2) built from 2 regular pentagons

        """
        return Morphism_from_pyflatsurf._create_morphism(
            self.codomain(), self.domain(), self._pyflatsurf_conversion
        )

    def _repr_type(self):
        r"""
        Helper method for printing this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: S.pyflatsurf()  # optional: pyflatsurf
            pyflatsurf conversion morphism:
              From: Triangulation of Translation Surface in H_2(2) built from 2 regular pentagons
              To:   Surface backed by FlatTriangulationCombinatorial(...) with vectors ...

        """
        return "pyflatsurf conversion"

    def _test_section_point(self, **options):
        r"""
        Do not verify that :meth:`_section_point` has been implemented
        correctly.

        Surfaces backed by pyflatsurf cannot represent points yet.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: to_pyflatsurf = S.pyflatsurf()  # optional: pyflatsurf
            sage: to_pyflatsurf._test_section_point()  # optional: pyflatsurf

        """
        try:
            super()._test_section_point(**options)
        except NotImplementedError:
            # This test is known to fail until #211 adds mapping of points through pyflatsurf morphisms
            return

        raise Exception("this test is expected to fail until #211 is merged")

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: S.pyflatsurf() == S.pyflatsurf()  # optional: pyflatsurf
            True

        """
        if not isinstance(other, Morphism_to_pyflatsurf):
            return False

        return (
            self.parent() == other.parent()
            and self._pyflatsurf_conversion == other._pyflatsurf_conversion
        )

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: hash(S.pyflatsurf()) == hash(S.pyflatsurf())  # optional: pyflatsurf
            True

        """
        return hash(self.parent())


class Morphism_from_pyflatsurf(SurfaceMorphism):
    r"""
    A trivial isomorphism from a pyflatsurf backed translation surface to a
    sage-flatsurf translation surface.

    You should not create such morphisms directly but only create them as the
    :meth:`Morphism_to_pyflatsurf.section` of another morphism.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
        sage: from_pyflatsurf = S.pyflatsurf().section()  # optional: pyflatsurf

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.morphism import Morphism_from_pyflatsurf
        sage: isinstance(from_pyflatsurf, Morphism_from_pyflatsurf)  # optional: pyflatsurf
        True

        sage: TestSuite(from_pyflatsurf).run()  # optional: pyflatsurf

    """

    def __init__(self, parent, pyflatsurf_conversion):
        super().__init__(parent)
        self._pyflatsurf_conversion = pyflatsurf_conversion

    def _image_half_edge(self, half_edge):
        r"""
        Helper method for :meth:`_image_homology_edge` and others.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: from_pyflatsurf = S.pyflatsurf().section()  # optional: pyflatsurf
            sage: from_pyflatsurf._image_half_edge(1R)  # optional: pyflatsurf
            ((1, 2, 3), 0)

        """
        face = tuple(self.domain().flat_triangulation().face(half_edge))
        label = type(self.domain())._normalize_label(face)
        edge = label.index(half_edge)
        return (label, edge)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: from_pyflatsurf = S.pyflatsurf().section()  # optional: pyflatsurf
            sage: from_pyflatsurf._image_homology_edge((1, 2, 3), 0, codomain=from_pyflatsurf.codomain().homology())  # optional: pyflatsurf
            B[((0, 0), 0)]

        """
        half_edge = label[edge]

        import pyflatsurf

        half_edge = pyflatsurf.flatsurf.HalfEdge(int(half_edge))

        return codomain(self._pyflatsurf_conversion._preimage_half_edge(half_edge))

    def _repr_type(self):
        r"""
        Helper method for printing this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: S.pyflatsurf().section()  # optional: pyflatsurf
            pyflatsurf reconversion morphism:
              From: Surface backed by FlatTriangulationCombinatorial(...) with vectors ...
              To:   Triangulation of Translation Surface in H_2(2) built from 2 regular pentagons

        """
        return "pyflatsurf reconversion"

    def _test_section_point(self, **options):
        r"""
        Do not verify that :meth:`_section_point` has been implemented
        correctly.

        Surfaces backed by pyflatsurf cannot represent points yet.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: from_pyflatsurf = S.pyflatsurf().section()  # optional: pyflatsurf
            sage: from_pyflatsurf._test_section_point()  # optional: pyflatsurf

        """
        try:
            super()._test_section_point(**options)
        except NotImplementedError:
            # This test is known to fail until #211 adds mapping of points through pyflatsurf morphisms
            return

        raise Exception("this test is expected to fail until #211 is merged")

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: S.pyflatsurf().section() == S.pyflatsurf().section()  # optional: pyflatsurf
            True

        """
        if not isinstance(other, Morphism_from_pyflatsurf):
            return False

        return (
            self.parent() == other.parent()
            and self._pyflatsurf_conversion == other._pyflatsurf_conversion
        )

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: hash(S.pyflatsurf().section()) == hash(S.pyflatsurf().section())  # optional: pyflatsurf
            True

        """
        return hash(self.parent())


class Morphism_Deformation(SurfaceMorphism):
    r"""
    A morphism of
    :class:`~flatsurf.geometry.pyflatsurf.surface.Surface_pyflatsurf` surfaces
    that is backed by a libflatsurf ``Deformation``.

    These morphisms are usually hidden deep inside the machinery of some
    complex morphism constructions.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: T = S.relabel({0: 1})
        sage: isomorphism = S.delaunay_decompose(codomain=T)  # optional: pyflatsurf

        sage: deformation = isomorphism._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf

        sage: deformation  # optional: pyflatsurf
        pyflatsurf deformation morphism:
          From: Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}
          To:   Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}
          Defn: FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)} → ...

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.morphism import Morphism_Deformation
        sage: isinstance(deformation, Morphism_Deformation)  # optional: pyflatsurf
        True

        sage: TestSuite(deformation).run()  # optional: pyflatsurf

    """

    def __init__(self, parent, deformation):
        super().__init__(parent)
        self._deformation = deformation

    def _repr_type(self):
        r"""
        Helper method for the printing of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: isomorphism = S.delaunay_decompose(codomain=S)  # optional: pyflatsurf

            sage: deformation = isomorphism._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf

            sage: deformation  # optional: pyflatsurf
            pyflatsurf deformation morphism:
              From: ...
              To: ...
              Defn: ...

        """
        return "pyflatsurf deformation"

    def _repr_defn(self):
        r"""
        Helper method for the printing of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: isomorphism = S.delaunay_decompose(codomain=S)  # optional: pyflatsurf

            sage: deformation = isomorphism._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf

            sage: deformation  # optional: pyflatsurf
            pyflatsurf deformation morphism:
              From: ...
              To: ...
              Defn: FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)} → ...

        """
        return repr(self._deformation)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: isomorphism = S.delaunay_decompose(codomain=S)  # optional: pyflatsurf

            sage: deformation = isomorphism._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf

            sage: H = deformation.domain().homology()  # optional: pyflatsurf
            sage: H.hom(deformation).matrix()  # optional: pyflatsurf
            [1 0]
            [0 1]

        """
        half_edge = label[edge]

        from pyflatsurf import flatsurf

        saddle_connection = flatsurf.SaddleConnection[
            type(self.domain().flat_triangulation())
        ](self.domain().flat_triangulation(), flatsurf.HalfEdge(half_edge))
        path = flatsurf.Path[type(self.domain().flat_triangulation())](
            saddle_connection
        )

        path = self._deformation(path)

        if not path:
            raise NotImplementedError(
                "cannot map edge through this deformation in pyflatsurf yet"
            )

        path = path.value()

        image = codomain.zero()
        for step in path:
            chain = step.chain()
            for edge, coefficient in chain:
                from flatsurf.geometry.pyflatsurf.conversion import RingConversion

                coefficient = RingConversion.from_pyflatsurf_from_elements(
                    [coefficient]
                ).section(coefficient)

                half_edge = edge.positive()

                face = tuple(self.codomain().flat_triangulation().face(half_edge))
                label = type(self.codomain())._normalize_label(face)
                edge = label.index(half_edge)

                image += coefficient * codomain((label, edge))

        return image

    def _test_section_point(self, **options):
        r"""
        Do not verify that :meth:`_section_point` has been implemented
        correctly.

        Surfaces backed by pyflatsurf cannot represent points yet.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.relabel({0: 1})
            sage: isomorphism = S.delaunay_decompose(codomain=T)  # optional: pyflatsurf
            sage: deformation = isomorphism._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf

            sage: deformation._test_section_point()  # optional: pyflatsurf

        """
        try:
            super()._test_section_point(**options)
        except NotImplementedError:
            # This test is known to fail until #211 adds mapping of points through pyflatsurf morphisms
            return

        raise Exception("this test is expected to fail until #211 is merged")

    def _test_pickling(self, **options):
        r"""
        Do not test that this morphism can be serialized and deserialized.

        This cannot work because libflatsurf deformations do not implement
        cereal serialization yet.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.relabel({0: 1})
            sage: isomorphism = S.delaunay_decompose(codomain=T)  # optional: pyflatsurf
            sage: deformation = isomorphism._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf

            sage: deformation._test_pickling()  # optional: pyflatsurf

        """
        try:
            super()._test_pickling(**options)
        except Exception:
            return

        raise Exception("libflatsurf is not expected to implement serialization yet")

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.relabel({0: 1})
            sage: S.delaunay_decompose(codomain=T)._factorization()._factorization()._morphisms[2] == S.delaunay_decompose(codomain=T)._factorization()._factorization()._morphisms[2]  # optional: pyflatsurf
            Traceback (most recent call last):
            ...
            NotImplementedError: deformations do not implement the == operator yet

        """
        if not isinstance(other, Morphism_Deformation):
            return False

        if self.parent() != other.parent():
            return False

        if self is other:
            return True

        # This is not implemented in Deformation in libflatsurf
        raise NotImplementedError("deformations do not implement the == operator yet")

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.relabel({0: 1})
            sage: hash(S.delaunay_decompose(codomain=T)._factorization()._factorization()._morphisms[2]) == hash(S.delaunay_decompose(codomain=T)._factorization()._factorization()._morphisms[2])  # optional: pyflatsurf
            True

        """
        # In libflatsurf there is no implementation for hashing of Deformation yet
        return hash((self.parent(), repr(self)))
