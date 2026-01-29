r"""
Absolute and relative (simplicial) homology of surfaces.

EXAMPLES:

The absolute homology of the regular octagon::

    sage: from flatsurf import translation_surfaces, SimplicialHomology
    sage: S = translation_surfaces.regular_octagon()
    sage: H = SimplicialHomology(S)

A basis of homology, with generators written as (sums of) oriented edges::

    sage: H.gens()
    ([(0, 1)], [(0, 2)], [(0, 3)], [(0, 0)])

The absolute homology of the unfolding of the (3, 4, 13) triangle::

    sage: from flatsurf import Polygon, similarity_surfaces
    sage: P = Polygon(angles=[3, 4, 13])
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")
    sage: S.genus()
    8
    sage: H = SimplicialHomology(S)
    sage: len(H.gens())
    16

Relative homology, relative to the singularities of the surface::

    sage: S = S.erase_marked_points()  # optional: pyflatsurf  # random output due to deprecation warnings
    sage: H1 = SimplicialHomology(S, relative=S.vertices())  # optional: pyflatsurf
    sage: len(H1.gens())  # optional: pyflatsurf
    17

We can also form relative `H_0` and `H_2`, though they are not overly
interesting of course::

    sage: H0 = SimplicialHomology(S, relative=S.vertices(), k=0)
    sage: len(H0.gens())
    0

    sage: H2 = SimplicialHomology(S, relative=S.vertices(), k=2)
    sage: len(H2.gens())
    1

We create the homology class corresponding to the core curve of a cylinder (the
interface here is terrible at the moment, see
https://github.com/flatsurf/sage-flatsurf/issues/166)::

    sage: from flatsurf import Polygon, similarity_surfaces, SimplicialHomology, GL2ROrbitClosure

    sage: P = Polygon(angles=[3, 4, 13])
    sage: S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation").triangulate().codomain()

    sage: from flatsurf.geometry.pyflatsurf.conversion import FlatTriangulationConversion  # optional: pyflatsurf
    sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)  # optional: pyflatsurf
    sage: T = conversion.codomain()  # optional: pyflatsurf

    sage: D = S._decomposition((13, 37))  # optional: pyflatsurf
    sage: cylinder = D.cylinders()[0]  # optional: pyflatsurf

    sage: H = SimplicialHomology(S)
    sage: C = H.chain_module()
    sage: core = H(sum(int(str(chain[edge])) * C(conversion.section(edge.positive())) for segment in cylinder.right() for chain in [segment.saddleConnection().chain()] for edge in T.edges()))  # optional: pyflatsurf
    sage: core  # optional: pyflatsurf  # random output, the chosen generators vary between operating systems
    972725347814111665129717*[((0, -1/2*c0, -1/2*c0^2 + 3/2), 2)] + 587352809047576581321682*[((0, -1/2*c0^2 + 1, -1/2*c0^3 + 3/2*c0), 2)] + 60771110563809382932401*[((0, -1/2*c0^2 + 1, 1/2*c0^3 - 3/2*c0), 2)] ...

"""

######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2026 Julian Rüth
#                           2023 Julien Boulanger
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
######################################################################
from typing import List, Tuple

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.morphism import Morphism
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.morphism import MorphismSpace


class SimplicialHomologyClass(Element):
    r"""
    An element of a homology group.

    INPUT:

    - ``parent`` -- a :class:`SimplicialHomology`

    - ``coefficients`` -- a vector of coefficients in the base ring, one for
      each element of the homology's generators.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: S = translation_surfaces.regular_octagon()
        sage: H0 = SimplicialHomology(S, k=0)
        sage: g0 = H0.gens()[0]
        sage: g0
        [Vertex 0 of polygon 0]

        sage: H1 = SimplicialHomology(S, k=1)
        sage: g1 = H1.gens()[0]
        sage: g1
        [(0, 1)]

        sage: H2 = SimplicialHomology(S, k=2)
        sage: g2 = H2.gens()[0]
        sage: g2
        [0]

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyClass
        sage: isinstance(g0, SimplicialHomologyClass)
        True

        sage: isinstance(g1, SimplicialHomologyClass)
        True

        sage: isinstance(g2, SimplicialHomologyClass)
        True

    """

    def __init__(self, parent, coefficients):
        super().__init__(parent)

        assert len(coefficients) == parent.ngens()
        assert not coefficients.is_mutable()

        self._coefficients = coefficients

    def algebraic_intersection(self, other):
        r"""
        Return the algebraic intersection of this class of a closed curve with
        ``other``.

        INPUT:

        - ``other`` - a :class:`SimplicialHomologyClass` in the same homology

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()
            sage: H = SimplicialHomology(S)

            sage: H((0, 0)).algebraic_intersection(H((0, 1)))
            1

            sage: a = H((0, 1))
            sage: b = 3 * H((0, 0)) + 5 * H((0, 2)) - 2 * H((0, 4))
            sage: a.algebraic_intersection(b)
            0

            sage: a = 2 * H((0, 0)) + H((0, 1)) + 3 * H((0, 2)) + H((0, 3)) + H((0, 4)) + H((0, 5)) + H((0, 7))
            sage: b = H((0, 0)) + 2 * H((0, 1)) + H((0, 2)) + H((0, 3)) + 2 * H((0, 4)) + 3 * H((0, 5)) + 4 * H((0, 6)) + 3 * H((0, 7))
            sage: a.algebraic_intersection(b)
            -6

            sage: S = translation_surfaces.cathedral(1, 4)
            sage: H = SimplicialHomology(S)
            sage: C = H.chain_module()
            sage: a = H((0, 3))
            sage: b = H((2, 1))
            sage: a.algebraic_intersection(b)
            0

            sage: a = H((0, 3))
            sage: b = H(C((3, 4)) + 3 * C((0, 3)) + 2 * C((0, 0)) - C((1, 7)) + 7 * C((2, 1)) - 2 * C((2, 2)))
            sage: a.algebraic_intersection(b)
            2

        """
        if not self.parent().is_absolute():
            raise NotImplementedError(
                "algebraic intersection only available for absolute homology classes"
            )

        other = self.parent()(other)

        if self.parent().degree() != 1:
            raise NotImplementedError(
                "algebraic intersections only available for homology in degree 1"
            )

        intersection = 0

        multiplicities = dict(self.chain())
        other_multiplicities = dict(other.chain())

        for vertex in self.parent().surface().vertices():
            counter = 0
            other_counter = 0

            for edge in [edge for (edge, _) in vertex.edges_ccw()[::2]]:
                opposite_edge = self.parent().surface().opposite_edge(*edge)

                counter += multiplicities.get(edge, 0)
                intersection += counter * other_multiplicities.get(edge, 0)
                intersection -= counter * other_multiplicities.get(opposite_edge, 0)

                counter -= multiplicities.get(opposite_edge, 0)
                other_counter += other_multiplicities.get(edge, 0)
                other_counter -= other_multiplicities.get(opposite_edge, 0)

            if counter:
                raise TypeError("homology class does not correspond to a closed curve")
            if other_counter:
                raise ValueError("homology class does not correspond to a closed curve")

        return intersection

    def _acted_upon_(self, c, self_on_left=None):
        r"""
        Return this homology class scaled by ``c``.

        INPUT:

        - ``c`` -- an element of the base ring of scalars

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: 3 * H.gens()[0]
            3*[(0, 1)]
            sage: H.gens()[0] * 0
            0

        """
        del self_on_left  # parameter intentionally ignored, the side does not matter
        return self.parent()(c * self._coefficients)

    def coefficients(self):
        r"""
        Return the coefficients of this element in terms of the generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0].coefficients()
            (1, 0)
            sage: H.gens()[1].coefficients()
            (0, 1)

        """
        return self._coefficients

    def holonomy(self):
        r"""
        Return the holonomy vector of this class.

        OUTPUT:

        A two-dimensional vector over the compositum of the coefficient ring of
        homology and the base ring of the translation surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0].holonomy()
            (0, 1)

            sage: H = H.change(relative=T.vertices())
            sage: H.gens()[0].holonomy()
            (0, 1)

        """
        from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
        if not self.surface() in TranslationSurfaces():  # pyright: ignore[reportCallIssue]
            raise NotImplementedError("cannot compute holonomies for non-translation surfaces yet")

        return self._coefficients * self.parent()._holonomy()

    def __eq__(self, other):
        r"""
        Return whether this class is indistinguishable from ``other``.

        .. NOTE::

            We cannot implement ``_richcmp_`` since we have to handle elements
            with differing parents here. Also, we cannot implement
            ``__richcmp__`` easily here since we are not in Cython.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0] == H.gens()[0]
            True
            sage: H.gens()[0] == H.gens()[1]
            False

        Since surfaces are not unique parents, this treats classes on equal
        surfaces as being equal::

            sage: S = translation_surfaces.square_torus()
            sage: T == S
            True
            sage: T is S
            False

        ::

            sage: h = T.homology().gens()[0]
            sage: g = S.homology().gens()[0]

            sage: g == h
            True
        
        ::

            sage: H.gens()[0] != H.gens()[0]
            False
            sage: H.gens()[0] != H.gens()[1]
            True

        """
        if self is other:
            return True

        if not isinstance(other, SimplicialHomologyClass):
            return False

        if self.parent() != other.parent():
            return False

        return self.coefficients() == other.coefficients()

    __ne__ = object.__ne__

    def __hash__(self):
        r"""
        Return a hash value of this class that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: hash(H.gens()[0]) == hash(H.gens()[0])
            True

        """
        return hash(self.coefficients())

    def _repr_(self):
        r"""
        Return a printable representation of this class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0]
            [(0, 1)]

        """
        return repr(self.chain())

    @cached_method
    def chain(self):
        r"""
        Return a lift of this element to the
        :meth:`SimplicialHomologyGroup.chain_module`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a.chain()
            [(0, 1)]

        We can use the chain representation to write a homology class as
        simplices, i.e., edges, with multiplicities::

            sage: coeffs = (a - b).chain().monomial_coefficients()
            sage: coeffs  # random output due to random ordering of edges
            {(0, 1): 1, (0, 0): -1}

        From this representation, we determine the holonomy vector that a chain
        encodes on a translation surface::

            sage: sum(c * T.polygon(label).edge(edge) for ((label, edge), c) in coeffs.items())
            (-1, 1)

        """
        return self.parent().chain_module()(self._coefficients * self.parent()._chain())

    def coefficient(self, gen):
        r"""
        Return the multiplicity of this class at a generator of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a.coefficient(a)
            1
            sage: a.coefficient(b)
            0

        TESTS::

            sage: a.coefficient(a + b)
            Traceback (most recent call last):
            ...
            ValueError: gen must be a generator not [(0, 0)] + [(0, 1)]

        """
        coefficients = gen.coefficients()
        indexes = [i for (i, c) in enumerate(coefficients) if c]

        if len(indexes) != 1 or coefficients[indexes[0]] != 1:
            raise ValueError(f"gen must be a generator not {gen}")

        index = indexes[0]

        return self.coefficients()[index]

    def _add_(self, other):
        r"""
        Return the formal sum of this class and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a + b
            [(0, 0)] + [(0, 1)]

        """
        return self.parent()(self._coefficients + other._coefficients)

    def _sub_(self, other):
        r"""
        Return the formal difference of this class and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a - b
            -[(0, 0)] + [(0, 1)]

        """
        return self.parent()(self._coefficients - other._coefficients)

    def _neg_(self):
        r"""
        Return the negative of this homology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a + b
            [(0, 0)] + [(0, 1)]
            sage: -(a + b)
            -[(0, 0)] - [(0, 1)]

        """
        return self.parent()(-self._coefficients)

    def surface(self):
        r"""
        Return the surface on which this homology class in defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: h = H.gens()[0]
            sage: h.surface()
            Translation Surface in H_1(0) built from a square

        """
        return self.parent().surface()

    def __bool__(self):
        r"""
        Return whether this class is non-trivial.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: h = H.gens()[0]
            sage: bool(h)
            True
            sage: bool(h-h)
            False

        """
        return bool(self._coefficients)

    def __setstate__(self, state):
        r"""
        Helper method to restore this element from a pickle.

        TESTS:

        Verify that we can unpickle old pickles that used to track the chain instead of the coefficient vector::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.cathedral(1, 4)
            sage: h = T.homology().gens()[-1]; h
            -[(1, 0)] - [(1, 6)] + [(2, 0)]

            sage: g = loads(
            ....:           b'x\x9c\x95Xw|\x14E\x14N\x83\xc0\x12@\x90"\x16,\x88\x1c\x96(`\xef'
            ....:           b"\x12\xaa\x91\x13'\xa8g\x89\xeb\xe5n\x92\xdd\xe4n7ov/&\xea\xd9oW\xec\x05{\xaf"
            ....:           b'\xd8{\xc5\x8e\xbd\xf7\xde\xc5^\xb0\xf7\xfa\xe6\xcd\xde\xee]H~\xea'
            ....:           b'?\xbb\xf3\xe6\xbd\xf7\xcd\x9b7\xdf\xd4\xc3\xaaRN\xb2\x8d\xd7;\xae'
            ....:           b'\xc8\xa5\xdc\x9c\xe0\xf5\xe9\x1e+\x995Sz*\x93t\x1c\xad\\\x82\x8aXbZEEE\x93'
            ....:           b'\x99\xed\xcc\x98)3\x99\x99cg\xed\x8c\xdd\xd63[\xd8\xb9N\xfd\x00\xd35\xf4'
            ....:           b'T\xd2\xe5m\xb6\xe8\xa9\xe7\x19\x9e\xe5\x96\xab|\xa12\xd5\x9aI\xbaNN\xb4\xd6'
            ....:           b'\xb7q;\xcb]\xb40\x02omE\xc0\x06\xd5`UJ\xd7[rf\xc65-]\xd7\xda\xb8\x9bt]\xa1'
            ....:           b'Au\xef\xa8s\x96\t9\xae\x0b\xde)\xb8\x83\x8d&]\xd3\xb6\xb4\x9c%x:'
            ....:           b'\x97\xe2\x1a\xd4(\x87 6\x93;\xf5Y;\x9d\xcbpG\x9b\x17\xfca\x802\x11'
            ....:           b'\xa6\xd5\xe6\xd4\x9b\x16\x1ar\xa1KI\x9b\xab\x04&\xcb0p2\x83\xda\x02\x0c\xca'
            ....:           b'\xc3\xe0\xc4 \xccE\xdat:\x93n\xca\x00m\xa1\xe3\xc3\x10\x06u\x89\xa1X]'
            ....:           b'\xde\xf9\xa1\x1e\x0cc0\xdc\x83\x95\xe2\xf1\xb8\x0b#\x18\x8c\x9c|8\xacl`>'
            ....:           b"'\xfe\xd7|\xc2\xa8\xff\x97A\x82\xd0`\xb4Qm\xd4\x18\x03\x8c\x153P\x04\xd6"
            ....:           b'v\xb6M\xab\xa1(\xc0\x98X\x1f\xb6i\x9e\xe6\x1d\xa6\x95\xd6\xd3v6iZ\x8e'
            ....:           b'6#\xa8\x98\x11\xc80vr\x1eV\xf1a\x1c\x83U\xfb\xf0\xe7\xb9T\xc6L\xf3'
            ....:           b'\xa4\x15\x02\xcc,\xd6\x84\x08\xab!\xc2\xea>\xac\xc1`|\x1f\x08\x96'
            ....:           b'\xcd]\x83\x0b\x13!h\x88\xb4xX\xc1H\x865\xd1\x7f-\x1f\xd6f\xb0\x8e\xf2'
            ....:           b'\xce\x9aN\n!2\x19M~t\xcc\x98a\xa75\x98`L\xe8\x03\xdf\xe1\xae#s\xad5a'
            ....:           b'A\x83u\x11l\xa2\x0f\xeb1\x98\x94\x18\x82C\xa4\xc6"\xd9m\xdaY\x88%4'
            ....:           b'\xac\x99i\xe5\xb2\\ @\x1a&\xfb\xb0>\x83\r\x8c\x181b\xae\xd5jZ\xa6\xcbaC\x1f'
            ....:           b'6bP\xdfGcr\xf0pV!sR\x92\x82$5)\x0166&\x15`\x93<L\xf1a*\x83i.lZ\x80'
            ....:           b'\xcd\xf2\xb0\xb9\x0f[0\xd8\xb2\x00[\xe5ak\x83\xc8\xb6\r\x83m\x13u\xd8bgRD\\'
            ....:           b'\xdb\xce\x83\xed\x19\xec\xe0\xc1\x8e\xf1\xb81\xda\x85\x9d\x18LG\xb25'
            ....:           b'\xe4a\x86\x8aO\x97\x14\xc2\xb6`\xa6\xe4_=\xd6\xcc\xcb\xb9\xc9\x96'
            ....:           b'\x0c\xdf\x15cC\xaa\xa7\x91Kf&)L\xb7\xa7IY\xf6"\xe2\xac>\x88\x18`j\xff\x06'
            ....:           b'\xa5\xc1l"\xe4\x18c\x82\xa1\x06bE\xac\x92D\x19\xc9L\xab\xee\x8a\xa4\xe5d'
            ....:           b'hJ\x17cw\xb49\xa8Y\x10)\x9a\x8a\xf50\x07Gn\xae\x0f;3h\x0c\x06d\xbe\xed'
            ....:           b'\x98\xae\xd9\xc5a\x17\x1f\xe61\x88\x1bj\x00g\xd1(-\xe8\xe9\xe4'
            ....:           b'\xb0\xab\x0f\xf3\x19\xec\x86\x8a\xe1\xa8\xd8\x13{k\xe7\xdc\xe9v\xceJ'
            ....:           b"'\xb1\xbf\xcc\x87&\x06\x0b\x82\x88\xff5\xe6\x94m\xf1(\xce\x06\x94"
            ....:           b'\xa2\xe0v\xc7\xe0\xf6\xf0aO\x06\x89\x15\xa3\xd8\xcb\x87\xbd\x19\xec\x13D]L!'
            ....:           b'\xec\xebC3\x83\xfd\xfa\tN\xf7a\x7f\x06I\xd4\x0eF-\xb6f\xf1\x94\xf4j'
            ....:           b'\xf1!\xc5 \xed\x01/@k\x1e\xda|0\x18\x98\x06\xb2\xa3\x9dA\x87\x07\x19d'
            ....:           b'\xc7l\x17\xb2\x0c,d\x87\x9d\x87\xce\x80\x1dm\x99\x9c\x9cO\x00y\x10\xb1>\xd6'
            ....:           b"D-\x9b\xec\xe0zQ\x00'Q\x89^\x9b\x80[\x80\x1c\x83\xaef8 \xd6X\xd1X\xe5A"
            ....:           b'wce\xe3`\x0fzP\xaa\xf0\xe0\xc0\xc6\xea\xc6j\x0f\x0e\xe2\x86r\x98\x02\x07'
            ....:           b'\x17 \xcf\xe0\x90f84\x86\x96\xa8;\x0cMj<8\x1c\x7f\x03=8\x02+\xd1\xef\xc8\xc6'
            ....:           b'*\xa9;\n\xa5A\x1e\x14P\x87\x95\x1e\xfe\xb0\x05\x1f+\x07xp4\xb6P\xe9'
            ....:           b'\xc1\xc2"\xf4T8\xa6\x00\xc728\xae\x19\x8e\x8f\xa1;Z\x9e\x80\x0e\xb5\x1e\x9c'
            ....:           b'\x88\x12\xba\x9f\x84~\xd8\xd0\xc9E\x87ipJ\x01Ne\xb0\xa8\x19N\x93\xb1`'
            ....:           b'\xeb\xa7\xa3\x03B\x9f\x81\x12\xfa\x9d\x89-`\x10g\xa1\x84\r\x9d\x8d'
            ....:           b':\xfc\x9d\x83\x12B\x9f\x8b\x98(\x9d\xc7sj\x9d\xb0x\xb7\xabg\x92-<\x03'
            ....:           b'\xe7\x1b\x8b\x12\x03e\x9d\xb0m\xd7\x81\x0b\x8c\xae\x02\\H\xc3\xa4w'
            ....:           b'\xda\x99\x9e6\xdbr\xe0\xa2<\\\x1c3\xba\xe4L\x1c/\x17\x95\xe2\xaa8_\x19'
            ....:           b'\xf4\x9ay\x97\xf4\xc1\xbc\x00J\xeb\xed\xaa\xc1\xa58\xd3\xfe\x03_\xa3'
            ....:           b'\xb5\xb9\x18\xd5\nX\xc8\xdc\xcbJ\xa9 h\xce%3z\xab\xc93i\x8d\x05\xe2'
            ....:           b',\x92\xe0r\xdc"\x17\x17\xe0\x8a<\\\xe9\xc3U\x0c\xaeFr\x0e,ns\x1c'
            ....:           b'\xae\xf1\xe1Z\x06\xd7\x05\x95\xc8\xd8.\xde\r\xd7\xfbp\x03\x83\x1b%?obp'
            ....:           b'\xb3\x07\xb7 ?/u\xe1V\x06\xb7!?o\xcf\xc3\x1d\xb1D\x95L]\x17\xdc\x19\xd02\xd8'
            ....:           b'\xc2\xeb\xbb\x90\xf16\xee\xd5\xc5\xa0\xd2\xdcr8n\xfe\x9df\xaa#\xc3\xf5\xae)'
            ....:           b'\x1a,\x89\xf5>1\xe0\x84t\xe5~\xd7\xc6-.\xd7\xe2@\xd6\x8bn\x1a\xdc\xd5'
            ....:           b'\xafO\xc6\xb6;p{n\xcb\xd8-\xc9\x8c\x06w\xab\x89,8W\x87\t\xb8\xa7'
            ....:           b'\x00\xf72\xb8\xafQ\x93\xe4\xb9?f,n\xacZ\x88\xdb\xfe\x03yX\xea\xc2\x83'
            ....:           b'\x0c\x1ej\x86\x87c}\xe5S\xcd\xadP\x82G\x0c\x9cV\x8f2x\xcc\xa0\xd2\xe3'
            ....:           b'\x0c\x9e\xe0\x88\xe5\xc2\x93\x0c\x9e2\x96\xc4\x0cDz:\x86J\x9cL\xcf0xV'
            ....:           b'\x99=\xc7\xe0ye\xf6\x02\x83\x17\x03\xb3\x97\x02\xb3\x97\x19\xbc\xa2J\xaf2xM'
            ....:           b'\x99\xbd\xce\xe0\x8d\xc0\xec\xcd\x98\xc2x\x8b\xc1\xdb\xca\xec\x1d\x06\xef*'
            ....:           b'\xb3\xf7\x18\xbc\xef\xc2\x07\x89\xd1r\x18\xc2\x03\x0e\xb2\x83\x12\x84#\x00'
            ....:           b'\xcbJ6\xde\xd6 \x8f\xc5\x84\xce\x0b\xb6\xdf\x0f\x13\x93\xfa\xf3o\x15vV/?'
            ....:           b"9}d\xdcn\xdc\xe6\xc3\xc7\x0c>Q3'\x9c\x0b\x9f\x1a7&\x06\xc8\x9a\x96"
            ....:           b'\xa4\xc3\xe13c\xb1\x9akx^\xe5\x0e|\x1e\xa7\xc5R/R\x80\x0b\x07S\n_4V\xd2'
            ....:           b"X\xe9<\xdb\xc2\xd3i\xcc=|\x19O\x8c\xea'\x1e\xf8\xca\xf8D)S\x92\xa4"
            ....:           b"\xc2\r\x8e\x10\xd4\x06,\x8f'\xc6\x95:\xcau\\\xef\xc4\x98u\xb5"
            ....:           b'\x1d\xc3\xd7\x0b\x95\xb3T\xe0\xb9\x0cA\xb8\xc0\xbd4c:.|\xd3\x0c\xdf\x96k'
            ....:           b'\x91Zr\x7f#\xedw\xcd\xf0\xbd\xcar\xe4\xab\x02 \xf5\x0f\xcd\xf0c\xae\xc58'
            ....:           b'\xc4\x903\xe3\xa7<\xfc\x1c3\xee\x8c\x05\xe3\xf7K0\xcc\xbf2\xf8M\x8d\xe4\xef'
            ....:           b"\x0c\xfeP\xe3\xf7'\x83\xbf\x02\xb3\xbf\x95\x99\xa8\xc0\xf5\x8f\x89J\xfc"
            ....:           b'\x1a\x8f\xd0\xf4\xdahSQ\x85RAT\x93\xa6\x06\xbf\xe4,\x06\x90<PZ\x12'
            ....:           b'\x84\xa8\xc5"\xa2$\xaa\xe5\n\xba\xf1T1\x88\xdc\x06\x93\x99\xa6\x00k$\xe0V'
            ....:           b'\xa8\x1bB\xba:\xd2\r\r!\x87\x91<<\x82\\IA\x1a\xc7\x14\xc4\x08\xd2\x8d,\x0ble'
            ....:           b'B\x19E\x9a\xd1!\xca\x18\x92\xc7F(\xabD(\xe3H\xb7*\xa1`.\xc4j$\xaf\x1e\xfa'
            ....:           b'\xaeA\xf2\xf8\xc8w\xcd\xc8w-\xd2\xad\xad|1U\xeb\x90<!\xf4]\x97\xe4'
            ....:           b"\x89\x91\xefz\x91\xef$\xd2\xc5T\xf4r\x8b\xd9LL\xa6\xe0\xd7'\xc5\x06!\xc8"
            ....:           b'\x86$o\x14\x81\xd4\xf7\xca\xea\xc6\xe4\xb6\t\x99MQxR\xd5\x82\xaa\xa9\xa4\x9a'
            ....:           b'F\xaaMC\xc4\xcdH\xde<B\xdc"\x08\x0b\xbb\xb0%\xe9\xb6"\x98v\x15\xd0\xd6T\xb5M'
            ....:           b'\xe8\xbe-\xc9\xdbE\xee\xdbG\xee;\x90n\xc70#;\x91<=\xf4m y\x06~]1SZ-3>4>2'
            ....:           b'~\xc2\t,f\x91n\xb6\xac\xc5\x89k\xe0t5>\x8f\x1b8\x1f\x8d/\xe3\xc6W\xed\xa4X'
            ....:           b'\x1e7\xbe^h|\xd3,\xe6H\xe9\xbbf1W\xfe\x7fh\x16;\xe3\x1f\t\x7f\x9c$\xbchD!/v'
            ....:           b'\xa1\xb0\x8a\xb4\x17\xf3\xa2\xdc\xc7\xa9\xa5]\xc31\x9fO\xf2na\x94'
            ....:           b'\x8c\xe4\xa6\xa8\x87\x0b\x02\xdfS\nbw\xd2\xed\x11\xfa\xeeIr"\xf4\xdd\x8b\xe4'
            ....:           b'\xbd#\xdf}"\xdf}I\xd7\x1cfg?\x92\xf5\xd0w\x7f\x92\x93\x91oK\x14s\x8at\xe9'
            ....:           b'\xd0\x97\x93\xdc\x1a\xfa\xb6\x91lPf\xcd0\xb3\xed2\x152\xb9\xed\xa4\xee\xe8/'
            ....:           b'\xb9\x1de\xc9\xcd\x04\xc9\xcd\x06\xc9\xb5Tr\x17QrmJngyr!\xea\xa4\xa0\x96\x9c'
            ....:           b'0A.\xc9\xb90\xd0.\x92\x0f\x88:\xd9]\xc2\xe7-\x90\xb4=\xc4\xba\x03\xc9'
            ....:           b'\xec\xa0\x92Ub\n\xea\x0e&]\x9et\x87\x84\x90\x87\x92|X\x04yx\t\xa4'
            ....:           b'\x9c\x07G\x90\xdb\x91dvT/\xc8\x02\xe9<\xd2\xf9!\xe4\xd1$/\x8c \x8f\t \xe5T'
            ....:           b'\xdd\\\x1cKN\xc7\x91\xd1\xf1aWO \xf9\xc4\x10\xe4$\x92O\x8e@NQ \xed\xca\xfd'
            ....:           b'TR/\n\x87\xf44\x92O\x0f\xdd\xcf \xf9\xcc\xc8\xfd\xac^\xdd:\x9b`\xce'
            ....:           b'!\xb3s\xa3\x99/\x17\x85\xf3Hu>\xa9.\x08\x11/$\xf9\xa2\x08\xf1\xe2^'
            ....:           b'\xb9\xbf\x84\xdc.%\xb3\xcb\xca\x11/\'\xd5bR]\x11"^I\xf2U\x11\xe2'
            ....:           b'\xd5\x11\x13\xae!\xdd\xb5a\xff\xae#\xf9\xfa\xd0\xf7\x06\x92o$\xca'
            ....:           b'\xde\x14Q\xd6\x0e({3\xa9o\xe9\x8f\xb2\xb7\x94Q\xf6\xd6\x80\xb2'
            ....:           b'\xb7\x05\x94\xbd\x9d(\x9bS\xf7\x95\xac\xbaq\x8a;\xb0ra\xb0\xe4'
            ....:           b'\xd8\x86\xe5\x8b;\xa9\x89%\xf8\xa5\x93\x01\x9d\xc9\x1dq\x17\xcae\xcf'
            ....:           b'\x03)\x83\xb7\xe6\xac\x94\xd6 KiuJ\x99\xaf\x8e.\xe2n\xd9\xa0\xdd.}<q\x0f\xe1'
            ....:           b'\xdd+\xf1F\xa8\x93\x08:\xe8\xba^|\x83\x11\xf7\xf5\x07\x1d\xb7\xd5'
            ....:           b'\x19\x08\xcf\x1b3\xcc\x94\xab\x89\xfb\xd1r2\x13\x0fH,:\xd5\x84\x18K\xb1\xaa'
            ....:           b'=hv)5\xfb 5\xfb\x904\x95\xfd-\x9e\xd1\xc5\xc3%\x96\x0f\x93\xe5#d\xf9\xa8\xb4'
            ....:           b'\x1c+A\x93-x\x92\xc1sE\xf9\x91J<&}\xaa\r\xd3\x18\xea\x89\xc7\xc9\xe5\tY3'
            ....:           b'T\xdeA\xb6\xfd\xdf\xaf\x01\xe5\xcf|\xe2I\x99\x83~\x9f\x07t\xbb'
            ....:           b'\xa5\x1d\x0f\xeb\x8e\x16@\xcd\xb7\xf1\x1e\xa9\x89\xa7dW\x9e\xa0><\x8d_\xf9F'
            ....:           b'&\x9e\xa1\xc8\x9e\x95\x9d\x19]\x9a\xed\xa4\xd5\xa6\x0et.\xef\x16\xcf'
            ....:           b'I\xc7 \x99\xcfK\xd3\xa1\xd4\xefR\x93\x17J\xd2\xf4\x025\xf1"!\xbf\xa4\xa8'
            ....:           b'g\xae@\xbd%e\xd4{9\xa0\xde+\x01\xf5^%\xea\xa9\x9bH\x87xM\xd6:\x92\xfe\xaf'
            ....:           b'\x13\xe8\x1ba\x0cx\xd8km5S2\x85\x8exS\x9a\xd5\xaa\x91\x16\\>X`po\xc9'
            ....:           b'L\x95\xben\xe2\xe9\xf7@ys\xc1\x84\xbc\x8d\xbaf\xf1\x0eM\xcaw\t\xf9=\x89LGZ'
            ....:           b'\xbaDe\x8bo\x9c\xe2}\xa9\xa8EEp\x91\x11\x1f\xc8\x0e\xd3\x18\x0f\xc71\x1d'
            ....:           b'I\x0fD)#\x89M\xa8\xfb\x92X\x16f\xa4\xa1}\x19e\xe4Cj\xe2\xa3"\x1b\x1d\xf5\x96'
            ....:           b'\xc8\x1d\xf1q\x89\xe9\xc7d\xfa\t\x99~*MkT\xa3\x8e\xf8\xac\xc4\xea3\xb2'
            ....:           b'\xfa\x9c\xac\xbe\x90\xf5KK\x94\x8a\xd0_\x92\xf2+\t1\x8cRE\xc1\xa5l\xd9'
            ....:           b'\xafn\xb1\xbc\xc4~9\xd9\x7fM\xf6\xdfD\xf6\x8a\t\xe8f\xb5q\xf1m\t\x05'
            ....:           b'\xbe+\xce\xf5@\xf7}\t\xd6\xf7\x84\xf5\x03a\xfd\xa8\xc6\xbe\x0e\xc7'
            ....:           b'\xbe\xb6t\xecq\xedXf4\x18\xd3}\xf1\x13\xd9\xfd\\F\x85_\x02*\xfc\x1a'
            ....:           b'P\xe17\xa2B^\xfc^lU\xf5D\xfc\x11-\x03\xc1\x05\x95&\x0c^SL\xf9\x88\x92\xe6'
            ....:           b'\xdd<]\x9c\x92Zx\xdf\x0ck\xc4\x9f\x12\xbf\xf8fhg[L+\xe9\xd6\xb7\xe2\xed2\x18'
            ....:           b'@\xad!\xa8\xc4\xeb;\xde\xb9\xc3k\xa7&\xfe"\x9e\x15}\xe5\x8bf\xbdz\x8a'
            ....:           b'\xd4y\xf8X\xa9K\x82\xa9G\xa7\xe8\x05\xb3I\x92\xeeoZ\xdc\xbb\x1a+=Q'
            ....:           b'\x81\xfc\xa6\xe7\rQ\x89%\xac\xc3\xd4U\xa9\xba\x01\x9e\xa8\x96%\xf9\x06"jT'
            ....:           b"\x1d\x96\x06(\xbbjO\x0cT\xdaZO\xd4\xca\x92|q\x11\x83\x94]\x8d'\x06\xab\x12j5"
            ....:           b'UB\xed\x10U\x1a\xe8\x89\xba*\xb9]\x0c\xad\x92\xbc\x1fV%\x8f"\xc3\xf1\xeb'
            ....:           b'\x8b\x95\xaa\xe4p\x8c\xa8\x92C8\x92\xeaW\xc6\xafz\xb0*.Db\x94\x0cA'
            ....:           b'\xbew\xd4\x05\x8fa\xf2\xa9l\xba\xcc\xba\x18M(c\x08e\xac4\x8b\xd1\x12\xae'
            ....:           b'\xd20\xc3\xc4\xbc;t\xe3\x16\xab\x90\xe182\\\xb5*\x18\xd7N\xc1[\xcdn'
            ....:           b'\xb1\x9a\x941G\x15\xd3\xc5\xea\xb2(/\x9ft\xdb\x14k\xa0\x18\xcf\xf9b<\xf9'
            ....:           b'\xadI\x01\xae%\x03l\x97\x1dW\x0b\xc4\xda\xa4[\x07\xbf\xedC\xa8R\xddj\xa6'
            ....:           b'\x88\t\xd4\xdbuI=Q\xaa\xeb\xca\xd4\xeb\x91z\x12\xa9c\xf8\xcdyb2\t'
            ....:           b'\xeb\xe3\xd7\xf1\xc4\x06\xf8k\xa9\xff\x07]\xc1`\xec')
            sage: g
            -[(1, 0)] - [(1, 6)] + [(2, 0)]

        """
        parent = state[0]
        if "_chain" in state[1]:
            if "coefficients" in state[1]:
                del state[1]["coefficients"]

            C = parent.chain_module()
            chain = sum((c * C(gen) for (gen, c) in list(state[1]["_chain"])), start=C.zero())

            state[1]["_coefficients"] = parent(chain)._coefficients

            del state[1]["_chain"]

        super().__setstate__(state)

    ### def lift(self, v):
    ###     # TODO: This should return a 1-chain instead so it's just a lift of a homology class. And thus should live on cohomology classes.
    ###     r"""
    ###     Given a vector in the "spanning set basis" return a vector on the full basis of
    ###     edges.

    ###     The vectors are returned as columns in a matrix. Each column
    ###     corresponds to one edge in the underlying :meth:`_flat_triangulation`
    ###     ordered as returned by its ``.edges()``.

    ###     EXAMPLES::

    ###         sage: from flatsurf import polygons, translation_surfaces, similarity_surfaces
    ###         sage: from flatsurf import GL2ROrbitClosure  # optional: pyflatsurf

    ###         sage: S = translation_surfaces.mcmullen_genus2_prototype(4,2,1,1,0)
    ###         sage: O = GL2ROrbitClosure(S)  # optional: pyflatsurf
    ###         sage: u0,u1 = O.tangent_space_basis()  # optional: pyflatsurf
    ###         sage: v0 = O.lift(u0)  # optional: pyflatsurf
    ###         sage: v1 = O.lift(u1)  # optional: pyflatsurf
    ###         sage: span([v0, v1])  # optional: pyflatsurf
    ###         Vector space of degree 9 and dimension 2 over Number Field in l with defining polynomial x^2 - x - 8 with l = 3.372281323269015?
    ###         Basis matrix:
    ###         [            1             0            -1   1/8*l + 7/8  -1/8*l + 1/8            -1   5/8*l - 5/8  -1/2*l + 3/2 -5/8*l + 13/8]
    ###         [            0             1            -1   1/4*l - 1/4  -1/4*l + 1/4             0   1/4*l - 1/4             0  -1/4*l + 1/4]

    ###     This can be used to deform the surface::

    ###         sage: T = polygons.triangle(3,4,13)
    ###         sage: S = similarity_surfaces.billiard(T)
    ###         sage: S = S.minimal_cover("translation").erase_marked_points() # long time (3s, #122), optional: pyflatsurf
    ###         sage: O = GL2ROrbitClosure(S)  # long time (above), optional: pyflatsurf
    ###         sage: for slope in S.slopes(bound=4): # long time (2s, #124), optional: pyflatsurf
    ###         ....:     d = S._decomposition(slope, limit=20)
    ###         ....:     O.update_tangent_space_from_flow_decomposition(d)
    ###         ....:     if O.dimension() == 4:
    ###         ....:         break
    ###         sage: d1,d2,d3,d4 = [O.lift(b) for b in O.tangent_space_basis()]  # long time (above), optional: pyflatsurf
    ###         sage: dreal = d1/132 + d2/227 + d3/1280 - d4/13201  # long time (above), optional: pyflatsurf
    ###         sage: dimag = d1/141 - d2/233 + d4/1230 + d4/14250  # long time (above), optional: pyflatsurf
    ###         sage: d = [O._vector_space_conversion()((x,y)) for x,y in zip(dreal,dimag)]  # long time (above), optional: pyflatsurf
    ###         sage: S2 = O._flat_triangulation() + d  # long time (6s), optional: pyflatsurf

    ###         sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf  # optional: pyflatsurf
    ###         sage: S2 = Surface_pyflatsurf(S2.surface())  # long time (above), optional: pyflatsurf
    ###         sage: O2 = GL2ROrbitClosure(S2)  # long time (above), optional: pyflatsurf
    ###         sage: for slope in S2.slopes(bound=1):  # long time (25s, #124), optional: pyflatsurf
    ###         ....:     d = S2._decomposition(slope, limit=20)
    ###         ....:     O2.update_tangent_space_from_flow_decomposition(d)

    ###     TESTS:

    ###     Verify that this also works with exact-real coefficients::

    ###         sage: from flatsurf import Polygon, EuclideanPolygonsWithAngles
    ###         sage: from pyexactreal import ExactReals  # optional: pyexactreal  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib

    ###         sage: E = EuclideanPolygonsWithAngles((1, 5, 5, 5))
    ###         sage: R = ExactReals(E.base_ring())  # optional: pyexactreal
    ###         sage: slopes = E.slopes()
    ###         sage: T = Polygon(angles=(1, 5, 5, 5), edges=[slopes[0], R.random_element(1/4) * slopes[1]])  # optional: pyexactreal
    ###         sage: S = similarity_surfaces.billiard(T)  # optional: pyexactreal
    ###         sage: S = S.minimal_cover(cover_type="translation")  # optional: pyexactreal
    ###         sage: O = GL2ROrbitClosure(S)  # optional: pyflatsurf, optional: pyexactreal
    ###         sage: d1, d2, d3, d4 = [O.lift(b) for b in O.tangent_space_basis()]  # optional: pyflatsurf, optional: pyexactreal

    ###     """
    ###     # given the values on the spanning edges we reconstruct the unique vector that
    ###     # vanishes on the boundary
    ###     bdry = self.boundaries()
    ###     n = self._flat_triangulation().edges().size()
    ###     k = len(self.spanning_set)
    ###     assert k + len(bdry) == n + 1
    ###     A = matrix(QQ, n + 1, n)
    ###     for i, e in enumerate(self.spanning_set):
    ###         A[i, e.index()] = 1
    ###     for i, b in enumerate(bdry):
    ###         A[k + i, :] = b

    ###     u = vector(self._surface.base_ring(), n + 1)
    ###     u[:k] = v


    ###     from pyexactreal.exact_reals import ExactReals
    ###     if isinstance(u.base_ring(), ExactReals):
    ###         u = u.change_ring(u.base_ring().base_ring())

    ###     return A.solve_right(u)


class SimplicialHomologyGroup(Parent):
    r"""
    The ``k``-th simplicial homology group of the ``surface`` with
    ``coefficients``.

    .. NOTE:

        This method should not be called directly since it leads to problems
        with pickling and uniqueness. Instead use :meth:`SimplicialHomology` or
        :meth:`homology` on a surface.

    INPUT:

    - ``surface`` -- a finite type surface without boundary

    - ``k`` -- an integer

    - ``coefficients`` -- a ring

    - ``relative`` -- a subset of points of the ``surface``

    - ``implementation`` -- a string; the algorithm used to compute the
      homology, only ``"generic"`` is supported at the moment which uses the
      generic homology machinery of SageMath.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialHomology, MutableOrientedSimilaritySurface
        sage: T = translation_surfaces.square_torus()

    Surfaces must be immutable to compute their homology::

        sage: T = MutableOrientedSimilaritySurface.from_surface(T)
        sage: SimplicialHomology(T)
        Traceback (most recent call last):
        ...
        ValueError: surface must be immutable to compute homology

    ::

        sage: T.set_immutable()
        sage: SimplicialHomology(T)
        H₁(Translation Surface in H_1(0) built from a square)

    TESTS::

        sage: T = translation_surfaces.square_torus()
        sage: H = SimplicialHomology(T, implementation="generic")
        sage: TestSuite(H).run()

    Verify that #310 has been resolved::

        sage: S = translation_surfaces.mcmullen_L(1,1,1,1)
        sage: H = S.homology(coefficients=QQ)
        sage: H.gens()
        ([(0, 1)], [(0, 0)], [(1, 1)], [(2, 0)])

    """

    Element = SimplicialHomologyClass

    def __init__(self, surface, k, coefficients, relative, implementation, category):
        Parent.__init__(self, base=coefficients, category=category)

        if surface.is_mutable():
            raise TypeError("surface must be immutable")

        from sage.all import ZZ

        if k not in ZZ:
            raise TypeError("k must be an integer")

        from sage.categories.all import Rings

        if coefficients not in Rings():  # pyright: ignore[reportCallIssue]
            raise TypeError("coefficients must be a ring")

        if relative:
            for point in relative:
                if point not in surface.vertices():
                    raise NotImplementedError(
                        "can only compute homology relative to a subset of the vertices"
                    )

        if implementation in ["generic", "spanning-set"]:
            if not surface.is_finite_type():
                raise NotImplementedError(
                    "homology only implemented for surfaces with finitely many polygons"
                )

            if surface.is_with_boundary():
                raise NotImplementedError(
                    "homology only implemented for surfaces without boundary"
                )
        else:
            raise NotImplementedError(
                "cannot compute homology with this implementation yet"
            )

        self._surface = surface
        self._k = k
        self._coefficients = coefficients
        self._relative = relative
        self._implementation = implementation

    def _an_element_(self):
        r"""
        Return a typical homology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: H = T.homology()
            sage: H.an_element()
            [(0, 0)] + [(0, 1)]

        """
        return sum(self.gens(), start=self.zero())

    def is_absolute(self):
        r"""
        Return whether this is absolute homology (and not relative to some set
        of points).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.is_absolute()
            True

        """
        return not self._relative

    def some_elements(self):
        r"""
        Return some typical homology classes (for testing).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.some_elements()
            [0, [(0, 1)], [(0, 0)]]

        """
        return [self.zero()] + list(self.gens())

    def surface(self):
        r"""
        Return the surface of which this is the homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.surface() == T
            True

        """
        return self._surface

    @cached_method
    def chain_module(self):
        r"""
        Return the free module of simplicial chains of the surface i.e., formal
        sums of simplicies, e.g., formal sums of edges of the surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.chain_module()
            C₁(Translation Surface in H_1(0) built from a square)

        """
        return SimplicialChains(surface=self.surface(), coefficients=self._coefficients, k=self._k, relative=self._relative)

    @cached_method(key=lambda _, k, relative: (k, None if relative is None else frozenset(relative)))
    def change(self, k=None, relative=None):
        r"""
        Return a variant of this homology.

        INPUT:

        - ``k`` -- if set, return this homology but in degree ``k``.
        - ``relative`` -- if set, return this homology but relative to
          ``relative`` instead; set to an empty tuple for absolute homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H
            H₁(Translation Surface in H_1(0) built from a square)

            sage: H.change(k=0)
            H₀(Translation Surface in H_1(0) built from a square)

            sage: H.change(relative=T.vertices())
            H₁(Translation Surface in H_1(0) built from a square, {Vertex 0 of polygon 0})

        """
        return SimplicialHomology(
            surface=self._surface,
            k=self._k if k is None else k,
            coefficients=self._coefficients,
            relative=self._relative if relative is None else relative,
            implementation=self._implementation,
            category=self.category(),
        )

    def zero(self) -> SimplicialHomologyClass:
        r"""
        Return the zero homology class.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.zero()
            0

        """
        return self(self.chain_module().zero())

    @cached_method
    def _homology(self):
        r"""
        Return the free module isomorphic to homology, a lift from that module
        to the free module isomorphic to the chain module, and a left inverse,
        (i.e., the map from cycles to homology.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._homology()
            (Ambient free module of rank 2 over the principal ideal domain Integer Ring,
             Generic endomorphism of Ambient free module of rank 2 over the principal ideal domain Integer Ring,
             Generic endomorphism of Ambient free module of rank 2 over the principal ideal domain Integer Ring)

        ::

            sage: T = translation_surfaces.cathedral(1, 3)
            sage: H = T.homology()
            sage: H._homology()
            (Ambient free module of rank 8 over the principal ideal domain Integer Ring,
             Generic morphism:
               From: Ambient free module of rank 8 over the principal ideal domain Integer Ring
               To:   Ambient free module of rank 13 over the principal ideal domain Integer Ring,
             Generic morphism:
               From: Ambient free module of rank 13 over the principal ideal domain Integer Ring
               To:   Ambient free module of rank 8 over the principal ideal domain Integer Ring)

            sage: H = T.homology(relative=T.vertices())
            sage: H._homology()
            (Ambient free module of rank 10 over the principal ideal domain Integer Ring,
             Generic morphism:
               From: Ambient free module of rank 10 over the principal ideal domain Integer Ring
               To:   Ambient free module of rank 13 over the principal ideal domain Integer Ring,
             Generic morphism:
               From: Ambient free module of rank 13 over the principal ideal domain Integer Ring
               To:   Ambient free module of rank 10 over the principal ideal domain Integer Ring)

        """
        C = self.chain_module()

        # We compute the spaces of cycles and boundaries over the integers.
        # (The relations are all integer anyway.)
        boundaries = C.change(k=self._k + 1)._boundary().transpose().image()

        generators = self._homology_generators()

        assert all(g.parent() is C for g in generators), "all generators of homology must be chains"
        assert all(not g.boundary() for g in generators), "all generators of homology must be cycles"

        # The spaces of chains and homology in terms of distinguished generators.
        # We are going to return these spaces with maps between them essentially.
        free_chains = C._chains()
        free_homology = self.base_ring() ** len(generators)

        # Construct the map lifting our _homology_generators(), i.e., free_homology -> free_chains.
        to_chain = free_homology.module_morphism(on_basis=lambda i: free_chains(generators[i].coefficients()), codomain=free_chains)

        # Construct the map reducing cycles to homology; we produce a matrix
        # that we then use to represent each cycle as <generators> +
        # boundaries.
        from sage.all import matrix
        T = matrix(tuple(g.coefficients() for g in generators) + boundaries.gens()).transpose()

        assert not C._boundary() * T, "<generators> + boundary must span the space of cycles"

        to_homology = free_chains.module_morphism(function=lambda x: free_homology(T.solve_right(x)[:len(generators)]), codomain=free_homology)

        assert all(to_homology(to_chain(gen)) == gen for gen in free_homology.gens())

        return free_homology, to_chain, to_homology

    def _homology_generators(self):
        r"""
        Return a set of generators of homology as chains that we want to use in our computations.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._homology_generators()
            ([(1, 0)], [(0, 1)])

        ::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T, implementation="spanning-set")
            sage: H._homology_generators()

        """
        if self._implementation == "generic":
            # Take the homology generator that we get from SageMath when naïvely taking the quotient cycles/boundaries.
            C = self.chain_module()

            cycles = C._boundary().right_kernel()
            boundaries = C.change(k=self._k + 1)._boundary().transpose().image()
            homology = cycles.quotient(boundaries)

            return tuple(C(C._chains()(gen.lift())) for gen in homology.gens())

        raise NotImplementedError("cannot compute generators of homology for this implementation yet")

    def _test_homology(self, **options):
        r"""
        Test that :meth:`_homology` computes homology correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._test_homology()

        """
        tester = self._tester(**options)

        homology, to_chain, to_homology = self._homology()
        chains = self.chain_module()._chains()

        tester.assertEqual(homology, to_homology.codomain())
        tester.assertEqual(homology, to_chain.domain())
        tester.assertEqual(chains, to_homology.domain())
        tester.assertEqual(chains, to_chain.codomain())

        for gen in homology.gens():
            tester.assertEqual(to_homology(to_chain(gen)), gen)

    @cached_method
    def _chain(self):
        r"""
        Helper method for :meth:`SimplicialHomologyClass.chain` that records
        the coefficients of a lift to a chain of the basis elements as a
        matrix.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._chain()
            [1 0]
            [0 1]

        """
        homology, to_chain, _ = self._homology()

        from sage.all import matrix
        return matrix(to_chain(gen) for gen in homology.gens())

    @cached_method
    def _holonomy(self):
        r"""
        Helper method for meth:`SimplicialHomologyClass.holonomy` that records
        the holonomies of the basis element as a matrix.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._holonomy()
            [0 1]
            [1 0]

        """
        if not self.surface().is_translation_surface():
            raise NotImplementedError("holonomy can only be computed in translation surfaces")

        if not self.surface().is_finite_type():
            raise NotImplementedError("holonomy can only be represented as a matrix in surfaces of finite type")

        if not self._k == 1:
            raise NotImplementedError("holonomy only implemented for 1-dimensional homology")

        from sage.all import matrix
        return matrix(h.chain().holonomy() for h in self.gens())

    def _repr_(self):
        r"""
        Return a printable representation of this homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H
            H₁(Translation Surface in H_1(0) built from a square)

        """
        k = self._k
        if k == 0:
            k = "₀"
        elif k == 1:
            k = "₁"
        elif k == 2:
            k = "₂"
        else:
            k = f"_{k}"

        H_k = f"H{k}"

        X = repr(self.surface())
        if not self.is_absolute():
            X = f"{X}, {set(self._relative)}"

        from sage.all import ZZ

        if self._coefficients is not ZZ:
            sep = ";"
            X = f"{X}{sep} {self._coefficients}"

        return f"{H_k}({X})"

    def _element_constructor_(self, x, check=True):
        r"""
        Return ``x`` as an element of this homology.

        TESTS::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)

        ::

            sage: H(0)
            0
            sage: H(None)
            0

        ::

            sage: H((0, 0))
            [(0, 0)]
            sage: H((0, 2))
            -[(0, 0)]

        ::

            sage: H = SimplicialHomology(T, 0)
            sage: H(H.chain_module().gens()[0])
            [Vertex 0 of polygon 0]

            sage: H = SimplicialHomology(T, 1)
            sage: H(H.chain_module().gens()[0])
            [(0, 1)]

            sage: H = SimplicialHomology(T, 2)
            sage: H(H.chain_module().gens()[0])
            [0]

        """
        homology, _, to_homology = self._homology()

        # If x is the zero element in its parent, we return the zero homology class.
        if x == 0 or x is None:
            return self.element_class(self, homology.zero())

        # We allow cycles to be specified directly from surface data.
        try:
            x = self.chain_module()(x)
        except NotImplementedError:
            pass

        # We turn a chain into a quotient element (if it is a cycle.)
        if x.parent() is self.chain_module():
            if check and x.boundary():
                raise ValueError("chain is not a cycle so it has no representation in this homology")
            x = to_homology(x.coefficients())

        # We interpret a coefficient vector as a homology element.
        if x.parent() is homology:
            if x.is_mutable():
                x = x.parent()(x)
                x.set_immutable()
            return self.element_class(self, x)

        # If nothing else worked, we look for a special _homology_ method that
        # can provide a custom cast to homology.
        try:
            hom_method = x._homology_
        except AttributeError:
            pass
        else:
            return hom_method(self)

        raise NotImplementedError("cannot convert this element to a homology class yet")

    @cached_method
    def gens(self) -> Tuple[SimplicialHomologyClass, ...]:
        r"""
        Return generators of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        ::

            sage: H = SimplicialHomology(T)
            sage: H.gens()
            ([(0, 1)], [(0, 0)])

        ::

            sage: H = SimplicialHomology(T, 0)
            sage: H.gens()
            ([Vertex 0 of polygon 0],)

        ::

            sage: H = SimplicialHomology(T, 2)
            sage: H.gens()
            ([0],)

        """
        if self._k < 0 or self._k > 2:
            return ()

        homology, _, _ = self._homology()
        return tuple(self(g) for g in homology.gens())

    @cached_method
    def ngens(self):
        r"""
        Return the Betti number of this homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

            sage: H = T.homology()
            sage: H.ngens()
            2

        """
        homology, _, _ = self._homology()
        return homology.ngens()

    def _test_ngens(self, **options):
        r"""
        Validate the Betti number of this homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

            sage: H = T.homology()
            sage: H._test_ngens()

        """
        tester = self._tester(**options)

        expected = 2 * self.surface().genus()
        if self._relative:
            expected += len(self._relative) - 1

        tester.assertEqual(expected, self.ngens())

    def degree(self):
        r"""
        Return the degree `k` for this homology `H_k`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

            sage: H = SimplicialHomology(T)
            sage: H.degree()
            1

        """
        return self._k

    def symplectic_basis(self) -> List[SimplicialHomologyClass]:
        r"""
        Return a symplectic basis of generators of this homology group.

        A symplectic basis is a basis of the form `(e_1, \ldots, e_n, f_1,
        \ldots, f_n)` such that for the
        :meth:`SimplicialHomologyClass.algebraic_intersection`, `e_i
        \cdot e_j = f_i \cdot f_j = 0` and `e_i \cdot f_j = \delta_{ij}` for
        all `i` and `j`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        ::

            sage: H = SimplicialHomology(T)
            sage: H.symplectic_basis()
            [[(0, 0)], [(0, 1)]]

        """
        from sage.all import matrix

        E = matrix(
            self.base_ring(),
            [[g.algebraic_intersection(h) for h in self.gens()] for g in self.gens()],
        )

        from sage.categories.all import Fields

        if self.base_ring() in Fields:
            from sage.matrix.symplectic_basis import symplectic_basis_over_field

            F, C = symplectic_basis_over_field(E)
        else:
            from sage.matrix.symplectic_basis import symplectic_basis_over_ZZ

            F, C = symplectic_basis_over_ZZ(E)

        if any(entry not in [-1, 0, 1] for row in F for entry in row):
            raise NotImplementedError(
                "cannot determine symplectic basis for this homology group over this ring yet"
            )

        return [
            sum(c * g for (c, g) in zip(row, self.gens())) for row in C
        ]  # pyright: ignore

    def _test_symplectic_basis(self, **options):
        r"""
        Verify that :meth:`symplectic_basis` has been implemented correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H._test_symplectic_basis()

        """
        tester = self._tester(**options)

        basis = self.symplectic_basis()
        n = len(basis)

        tester.assertEqual(len(self.gens()), n)
        tester.assertEqual(n % 2, 0)

        A = basis[: n // 2]
        B = basis[n // 2 :]

        for i, a in enumerate(A):
            for j, b in enumerate(B):
                tester.assertEqual(a.algebraic_intersection(b), i == j)

        for a in A:
            for aa in A:
                tester.assertEqual(a.algebraic_intersection(aa), 0)

        for b in B:
            for bb in B:
                tester.assertEqual(b.algebraic_intersection(bb), 0)

    def hom(self, f, codomain=None):
        r"""
        Return the homomorphism of homology induced by ``f``.

        INPUT:

        - ``f`` -- a morphism of surfaces or a matrix

        - ``codomain`` -- the simplicial homology this morphism maps into

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: g.matrix()  # optional: pyflatsurf
            [1 0]
            [2 1]

            sage: H.gens()
            ([(0, 1)], [(0, 0)])
            sage: [g(h) for h in H.gens()]  # optional: pyflatsurf
            [2*[(0, 0)] + [(0, 1)], [(0, 0)]]

        """
        from flatsurf.geometry.veech_group import SurfaceMorphism
        from sage.matrix.matrix0 import Matrix
        from sage.all import Hom

        if isinstance(f, SurfaceMorphism):
            if codomain is None:
                # TODO: This is wrong.
                codomain = f.codomain().homology()

            if codomain.surface() is not f.codomain():
                raise ValueError("codomain must be codomain of morphism or None")

            if f.domain() is self.surface():
                parent = Hom(self, codomain)
                return parent.__make_element_class__(
                    SimplicialHomologyMorphism_induced
                )(parent, f)
        elif isinstance(f, Matrix):
            if codomain is None:
                if f.is_square():
                    codomain = self

            if codomain is None:
                raise NotImplementedError("cannot deduce codomain from this matrix")

            if f.ncols() != self.ngens():
                raise ValueError(
                    "matrix must have one column for each generator of homology"
                )

            if f.nrows() != codomain.ngens():
                raise ValueError(
                    "matrix must have one row for each generator of the codomain"
                )

            parent = Hom(self, codomain)
            return parent.__make_element_class__(SimplicialHomologyMorphism_matrix)(
                parent, f
            )

        raise NotImplementedError(
            "cannot create a morphism in homology from this data yet"
        )

    def _Hom_(self, Y, category=None):
        r"""
        Return the space of morphisms from this homology to ``Y``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: H = S.homology()

            sage: End(H)
            Endomorphisms of H₁(Translation Surface in H_1(0) built from a square)
        """
        if isinstance(Y, SimplicialHomologyGroup):
            return SimplicialHomologyMorphismSpace(self, Y, category=category)

        return super()._Hom_(Y, category=category)

    def __eq__(self, other):
        r"""
        Return whether this homology is indistinguishable from ``other``.

        .. NOTE::

            We cannot rely on the builtin `==` by ``id`` since we need to
            detect homologies over equal but distinct surfaces to be equal. See
            :meth:`homology` for ideas on how to fix this.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)

            sage: T = translation_surfaces.square_torus()
            sage: HH = SimplicialHomology(T)

            sage: H == HH
            True

        """
        if not isinstance(other, SimplicialHomologyGroup):
            return False

        return (
            self._surface == other._surface
            and self._coefficients == other._coefficients
            and self._relative == other._relative
            and self._implementation == other._implementation
            and self.category() == other.category()
        )

    def __hash__(self):
        r"""
        Return a hash value for this homology that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)

            sage: T = translation_surfaces.square_torus()
            sage: HH = SimplicialHomology(T)

            sage: hash(H) == hash(HH)
            True

        """
        return hash(
            (
                self._surface,
                self._coefficients,
                self._relative,
                self._implementation,
                self.category(),
            )
        )


class SimplicialChain(Element):
    r"""
    A chain of simplexes, i.e., a formal sum.

    INPUT:

    - ``parent`` -- a :class:`SimplicialChainModule`

    - ``coefficients`` -- a vector of coefficients in the base ring, one for
      each simplex in this dimension

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: C0 = S.chains(k=0)
        sage: g0 = C0.gens()[0]
        sage: g0
        [Vertex 0 of polygon 0]

        sage: C1 = S.chains(k=1)
        sage: g1 = C1.gens()[0]
        sage: g1
        [(0, 1)]

        sage: C2 = S.chains(k=2)
        sage: g2 = C2.gens()[0]
        sage: g2
        [0]

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialChain
        sage: isinstance(g0, SimplicialChain)
        True

        sage: isinstance(g1, SimplicialChain)
        True

        sage: isinstance(g2, SimplicialChain)
        True

    """
    def __init__(self, parent, coefficients):
        super().__init__(parent)

        assert len(coefficients) == parent.ngens()
        assert not coefficients.is_mutable()

        self._coefficients = coefficients

    def _acted_upon_(self, c, self_on_left=None):
        r"""
        Return this chain scaled by ``c``.

        INPUT:

        - ``c`` -- an element of the base ring of scalars

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: 3 * C.gens()[0]
            3*[(0, 1)]
            sage: C.gens()[0] * 0
            0

        """
        del self_on_left  # parameter intentionally ignored, the side does not matter
        return self.parent()(c * self._coefficients)

    def coefficients(self):
        r"""
        Return the coefficients of this element in terms of the generator simplices.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C.gens()[0].coefficients()
            (1, 0)
            sage: C.gens()[1].coefficients()
            (0, 1)

        """
        return self._coefficients

    def holonomy(self):
        r"""
        Return the holonomy vector of this chain.

        OUTPUT:

        A two-dimensional vector over the compositum of the coefficient ring of
        the chains and the base ring of the translation surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C.gens()[0].holonomy()
            (0, 1)

        """
        from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
        if not self.surface() in TranslationSurfaces():  # pyright: ignore[reportCallIssue]
            raise NotImplementedError("cannot compute holonomies for non-translation surfaces yet")

        return sum((c * self.surface().polygon(label).edge(edge) for ((label, edge), c) in self), start=(self.base_ring()**2).zero())

    def __eq__(self, other):
        r"""
        Return whether this chain is indistinguishable from ``other``*.

        .. NOTE::

            We cannot implement ``_richcmp_`` since we have to handle elements
            with differing parents here. Also, we cannot implement
            ``__richcmp__`` easily here since we are not in Cython.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialChains
            sage: T = translation_surfaces.square_torus()
            sage: C = SimplicialChains(T)
            sage: C.gens()[0] == C.gens()[0]
            True

        Since surfaces are not unique parents, this treats chains on equal
        surfaces as being equal::

            sage: S = translation_surfaces.square_torus()
            sage: T == S
            True
            sage: T is S
            False

        ::

            sage: h = T.chains().gens()[0]
            sage: g = S.chains().gens()[0]

            sage: g == h
            True
        
        ::

            sage: C.gens()[0] != C.gens()[0]
            False
            sage: C.gens()[0] != C.gens()[1]
            True

        """
        if self is other:
            return True

        if not isinstance(other, SimplicialChain):
            return False

        if self.parent() != other.parent():
            return False

        return self.coefficients() == other.coefficients()

    __ne__ = object.__ne__

    def __hash__(self):
        r"""
        Return a hash value of this cchain that is compatible with
        :meth:`_richcmp_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: hash(C.gens()[0]) == hash(C.gens()[0])
            True

        """
        return hash(self.coefficients())

    def _repr_(self):
        r"""
        Return a printable representation of this chain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: H = SimplicialHomology(T)
            sage: H.gens()[0]
            [(0, 1)]

        """
        from sage.all import CombinatorialFreeModule
        R = CombinatorialFreeModule(self.base_ring(), self.parent().simplices(), prefix="")
        return str(R.from_vector(self.coefficients()))

    def coefficient(self, gen):
        r"""
        Return the multiplicity of this class at a generator simplex.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: a, b = C.gens()
            sage: a.coefficient(a)
            1
            sage: a.coefficient(b)
            0

        TESTS::

            sage: a.coefficient(a + b)
            Traceback (most recent call last):
            ...
            ValueError: gen must be a generator not [(0, 0)] + [(0, 1)]

        """
        coefficients = gen.coefficients()
        indexes = [i for (i, c) in enumerate(coefficients) if c]

        if len(indexes) != 1 or coefficients[indexes[0]] != 1:
            raise ValueError(f"gen must be a generator not {gen}")

        index = indexes[0]

        return self.coefficients()[index]

    def _add_(self, other):
        r"""
        Return the formal sum of this chain and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: a, b = C.gens()
            sage: a + b
            [(0, 0)] + [(0, 1)]

        """
        return self.parent()(self._coefficients + other._coefficients)

    def _sub_(self, other):
        r"""
        Return the formal difference of this chain and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: a, b = C.gens()
            sage: a - b
            -[(0, 0)] + [(0, 1)]

        """
        return self.parent()(self._coefficients - other._coefficients)

    def _neg_(self):
        r"""
        Return the negative of this chain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: a, b = C.gens()
            sage: a + b
            [(0, 0)] + [(0, 1)]
            sage: -(a + b)
            -[(0, 0)] - [(0, 1)]

        """
        return self.parent()(-self._coefficients)

    def surface(self):
        r"""
        Return the surface on which this chain in defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: h = C.gens()[0]
            sage: h.surface()
            Translation Surface in H_1(0) built from a square

        """
        return self.parent().surface()

    def __bool__(self):
        r"""
        Return whether this chain is non-trivial.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: h = C.gens()[0]
            sage: bool(h)
            True
            sage: bool(h-h)
            False

        """
        return bool(self._coefficients)

    def boundary(self):
        r"""
        Return the boundary of this chain as a chain in lower degree.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

        ::

            sage: C = T.chains(k=0)
            sage: c = C.an_element(); c
            [Vertex 0 of polygon 0]
            sage: c.boundary()
            0

        ::

            sage: C = T.chains()
            sage: c = C.an_element(); c
            [(0, 0)] + [(0, 1)]
            sage: c.boundary()
            0

        ::

            sage: C = T.chains(k=2)
            sage: c = C.an_element(); c
            [0]
            sage: c.boundary()
            0

        """
        return self.parent().change(k=self.parent().degree() - 1)(self.parent()._boundary() * self.coefficients())

    def is_cycle(self):
        r"""
        Return whether this chain is a cycle.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.cathedral(1, 3)

        ::

            sage: C = T.chains(k=0)
            sage: c = C.an_element(); c
            [Vertex 0 of polygon 0] + [Vertex 1 of polygon 1] + [Vertex 0 of polygon 1]
            sage: c.is_cycle()
            True

        ::

            sage: C = T.chains()
            sage: c = C.an_element(); c
            [(0, 0)] + [(0, 1)] + [(0, 3)] + [(1, 0)] + [(1, 1)] + [(1, 2)] + [(1, 4)] + [(1, 5)] + [(1, 6)] + [(1, 7)] + [(2, 0)] + [(3, 1)] + [(3, 7)]
            sage: c.is_cycle()
            False

        ::

            sage: C = T.chains(k=2)
            sage: c = C.an_element(); c
            [0] + [1] + [2] + [3]
            sage: c.is_cycle()
            True

        """
        return not self.boundary()

    def is_boundary(self):
        r"""
        Return whether this chain is a boundary.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.cathedral(1, 3)

        ::

            sage: C = T.chains(k=0)
            sage: c = C.an_element(); c
            [Vertex 0 of polygon 0] + [Vertex 1 of polygon 1] + [Vertex 0 of polygon 1]
            sage: c.is_boundary()
            False
            sage: c.boundary().is_boundary()
            True

        ::

            sage: C = T.chains()
            sage: c = C.an_element(); c
            [(0, 0)] + [(0, 1)] + [(0, 3)] + [(1, 0)] + [(1, 1)] + [(1, 2)] + [(1, 4)] + [(1, 5)] + [(1, 6)] + [(1, 7)] + [(2, 0)] + [(3, 1)] + [(3, 7)]
            sage: c.is_boundary()
            False
            sage: c.boundary().is_boundary()
            True

        ::

            sage: C = T.chains(k=2)
            sage: c = C.an_element(); c
            [0] + [1] + [2] + [3]
            sage: c.is_boundary()
            False
            sage: c.boundary().is_boundary()
            True

        """
        return self.coefficients() in self.parent().change(self.parent().degree() + 1)._boundary().transpose().image()

    def __iter__(self):
        r"""
        Return an iterator over the generators together with their coefficients
        if they are non-zero.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

        ::

            sage: C = T.chains(k=0)
            sage: c = C.an_element(); c
            [Vertex 0 of polygon 0]
            sage: list(c)
            [(Vertex 0 of polygon 0, 1)]

        ::

            sage: C = T.chains()
            sage: c = C.an_element(); c
            [(0, 0)] + [(0, 1)]
            sage: list(c)
            [((0, 1), 1), ((0, 0), 1)]

        ::

            sage: C = T.chains(k=2)
            sage: c = C.an_element(); c
            [0]
            sage: list(c)
            [(0, 1)]

        """
        for simplex, coefficient in zip(self.parent().simplices(), self.coefficients()):
            if coefficient:
                yield simplex, coefficient

    def monomial_coefficients(self):
        r"""
        Return the coefficients at the basis elements as a ``dict``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

        ::

            sage: C = T.chains(k=0)
            sage: c = C.an_element(); c
            [Vertex 0 of polygon 0]
            sage: c.monomial_coefficients()
            {Vertex 0 of polygon 0: 1}

        ::

            sage: C = T.chains()
            sage: c = C.an_element(); c
            [(0, 0)] + [(0, 1)]
            sage: c.monomial_coefficients()
            {(0, 0): 1, (0, 1): 1}

        ::

            sage: C = T.chains(k=2)
            sage: c = C.an_element(); c
            [0]
            sage: c.monomial_coefficients()
            {0: 1}

        """
        return dict(self)


class SimplicialChainModule(Parent):
    r"""
    The free module of formal sums of ``k``-simplices.

    INPUT:

    - ``surface`` -- a finite type surface without boundary

    - ``k`` -- an integer

    - ``coefficients`` -- a ring

    - ``relative`` -- a set of vertices

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, SimplicialChains
        sage: T = translation_surfaces.square_torus()

        sage: SimplicialChains(T)
        C₁(Translation Surface in H_1(0) built from a square)

    TESTS::

        sage: C = SimplicialChains(T)
        sage: TestSuite(C).run()

    """

    Element = SimplicialChain

    def __init__(self, surface, k, coefficients, relative, category):
        Parent.__init__(self, base=coefficients, category=category)

        if surface.is_mutable():
            raise TypeError("surface must be immutable")

        from sage.all import ZZ

        if k not in ZZ:
            raise TypeError("k must be an integer")

        from sage.categories.all import Rings

        if coefficients not in Rings():  # pyright: ignore[reportCallIssue]
            raise TypeError("coefficients must be a ring")

        if relative:
            for point in relative:
                if point not in surface.vertices():
                    raise NotImplementedError(
                        "can only compute chains relative to a subset of the vertices"
                    )

        self._surface = surface
        self._k = k
        self._coefficients = coefficients
        self._relative = relative

    def _an_element_(self):
        r"""
        Return a typical chain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C.an_element()
            [(0, 0)] + [(0, 1)]

        """
        return sum(self.gens(), start=self.zero())

    def some_elements(self):
        r"""
        Return some typical chains (for testing).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C.some_elements()
            [0, [(0, 1)], [(0, 0)]]

        """
        return [self.zero()] + list(self.gens())

    def surface(self):
        r"""
        Return the surface over which these chains are defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C.surface() == T
            True

        """
        return self._surface

    @cached_method
    def _chains(self):
        r"""
        Return a free module that can be used to implement the underlying chain
        machinery.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

        ::

            sage: C = T.chains()
            sage: C._chains()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring

        """
        return self.base_ring() ** self.ngens()

    @cached_method
    def gens(self) -> Tuple[SimplicialChain, ...]:
        r"""
        Return generators of the free chain module.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

        ::

            sage: C = T.chains()
            sage: C.gens()
            ([(0, 1)], [(0, 0)])

        ::

            sage: C = T.chains(k=0)
            sage: C.gens()
            ([Vertex 0 of polygon 0],)

        ::

            sage: C = T.chains(k=2)
            sage: C.gens()
            ([0],)

        """
        return tuple(self(g) for g in self._chains().gens())

    def ngens(self):
        r"""
        Return the number of simplexes that generate this free module.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C.ngens()
            2

        ::

            sage: C.change(k=3).ngens()
            0

        """
        return len(self.simplices())

    @cached_method
    def simplices(self):
        r"""
        Return the simplices that form the generators of this module.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()

        In dimension 1, this is the set of edges::

            sage: C = T.chains()
            sage: C.simplices()
            ((0, 1), (0, 0))

        In dimension 0, this is the set of vertices::

            sage: C = T.chains(k=0)
            sage: C.simplices()
            (Vertex 0 of polygon 0,)

        In dimension 2, this is the set of polygons::

            sage: C = T.chains(k=2)
            sage: C.simplices()
            (0,)

        In all other dimensions, there are no simplices::

            sage: C = T.chains(k=12)
            sage: C.simplices()
            ()

        """
        if self._k == 0:
            return tuple(
                vertex
                for vertex in self._surface.vertices()
                if vertex not in self._relative
            )
        if self._k == 1:
            simplices = set()
            for edge in self._surface.edges():
                if self._surface.opposite_edge(*edge) not in simplices:
                    simplices.add(edge)
            return tuple(simplices)
        if self._k == 2:
            return tuple(self._surface.labels())

        return ()

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as an element of this chain module.

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()

        ::

            sage: C(0)
            0
            sage: C(None)
            0

        ::

            sage: C((0, 0))
            [(0, 0)]
            sage: C((0, 2))
            -[(0, 0)]

        Note that we do not "allow" generating chains from polygon labels
        (since they are very likely to clash with other constructions)::

            sage: C(0)
            0

        """
        chains = self._chains()

        if x == 0 or x is None:
            return self.element_class(self, chains.zero())

        # We allow chains to be specified directly from surface data.
        if self._k == 0:
            if isinstance(x, tuple) and len(x) == 2:
                x = self.surface().point(*x)

            from flatsurf.geometry.surface_objects import SurfacePoint
            if isinstance(x, SurfacePoint) and x in self.surface().vertices():
                if x in self._relative:
                    return self.element_class(self, chains.zero())
                x = chains.gen(self.simplices().index(x))
        if self._k == 1:
            if isinstance(x, tuple) and len(x) == 2:
                if x in self.simplices():
                    x = chains.gen(self.simplices().index(x))
                elif self.surface().opposite_edge(*x) in self.simplices():
                    x = -chains.gen(self.simplices().index(self.surface().opposite_edge(*x)))
                    x.set_immutable()

        # We allow elements to be given as vectors in the underlying free implementation.
        if x.parent() is chains:
            if x.is_mutable():
                x = x.parent()(x)
                x.set_immutable()
            return self.element_class(self, x)

        # If nothing else worked, we look for a special _chain_ method that
        # can provide a custom cast to homology.
        try:
            chain_method = x._chain_
        except AttributeError:
            pass
        else:
            return chain_method(self)

        raise NotImplementedError("cannot convert this element to a chain yet")

    def degree(self):
        r"""
        Return the degree `k` for this chain module `C_k`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()

            sage: C = T.chains()
            sage: C.degree()
            1

        """
        return self._k

    def change(self, k=None, relative=None):
        r"""
        Return a variant of this chain module.

        INPUT:

        - ``k`` -- if set, return this module but in degree ``k``.
        - ``relative`` -- if set, return this module but as the quotient
          relative to ``relative`` instead; set to an empty tuple for the full
          chain module.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C
            C₁(Translation Surface in H_1(0) built from a square)

            sage: C.change(k=0)
            C₀(Translation Surface in H_1(0) built from a square)

            sage: C.change(relative=T.vertices())
            C₁(Translation Surface in H_1(0) built from a square)/C₁({Vertex 0 of polygon 0})

        """
        return SimplicialChains(
            surface=self._surface,
            k=self._k if k is None else k,
            coefficients=self._coefficients,
            relative=self._relative if relative is None else relative,
            category=self.category(),
        )

    @cached_method
    def _boundary(self):
        r"""
        Return a matrix that represents the boundary on simplices.

        The matrix returned is such that multiplying it with a coefficient
        vector of a chain yields the coefficient vector of the boundary chain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.square_torus()
            sage: C = T.chains()
            sage: C._boundary()
            [0 0]

            sage: C.change(k=2)._boundary()
            [0]
            [0]

            sage: C.change(k=0)._boundary()
            []

            sage: C.change(relative=T.vertices())._boundary()
            []

        """
        boundary = []

        C = self.change(k=self._k-1)

        for gen in self.simplices():
            if self._k == 0:
                # All vertices map to 0
                boundary.append([])
            elif self._k == 1:
                # Edges map to the difference of their end points
                boundary.append(
                    (C(self.surface().opposite_edge(*gen)) - C(gen)).coefficients()
                )
            elif self._k == 2:
                # Faces map to the sum of their edges
                boundary.append(
                    (sum([C((gen, edge)) for edge in range(len(self.surface().polygon(gen).edges()))], start=C.zero())).coefficients())
            else:
                assert False, "there can only be simplices in degrees 0, 1, and 2"

        from sage.all import matrix, ZZ
        return matrix(boundary, base_ring=ZZ).transpose()

    def is_absolute(self):
        return not self._relative

    def _repr_(self):
        r"""
        Return a printable representation of this chain module.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.square_torus()
            sage: T.chains()
            C₁(Translation Surface in H_1(0) built from a square)
            sage: T.chains(relative=T.vertices())
            C₁(Translation Surface in H_1(0) built from a square)/C₁({Vertex 0 of polygon 0})

        """
        k = self._k
        if k == 0:
            k = "₀"
        elif k == 1:
            k = "₁"
        elif k == 2:
            k = "₂"
        else:
            k = f"_{k}"

        C_k = f"C{k}"

        X = repr(self.surface())

        from sage.all import ZZ
        if self._coefficients is not ZZ:
            sep = ";"
            X = f"{X}{sep} {self._coefficients}"

        if not self.is_absolute():
            A = f"{set(self._relative)}"
            if self._coefficients is not ZZ:
                sep = ";"
                A = f"{A}{sep} {self._coefficients}"

            return f"{C_k}({X})/{C_k}({A})"

        return f"{C_k}({X})"

    def hom(self, f, codomain=None):
        r"""
        Return the homomorphism of chains induced by ``f``.

        INPUT:

        - ``f`` -- a morphism of surfaces or a matrix

        - ``codomain`` -- the chains this morphism maps into

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)

            sage: g.matrix()  # optional: pyflatsurf
            [0 0]
            [0 1]
            [1 0]

            sage: C.gens()
            ([(0, 1)], [(0, 0)])
            sage: [g(h) for h in C.gens()]  # optional: pyflatsurf
            [[((0, 0), 1)], [((0, 0), 0)]]

        """
        from flatsurf.geometry.veech_group import SurfaceMorphism
        from sage.matrix.matrix0 import Matrix
        from sage.all import Hom

        if isinstance(f, SurfaceMorphism):
            if codomain is None:
                # TODO: This is wrong.
                codomain = f.codomain().chains()

            if codomain.surface() is not f.codomain():
                raise ValueError("codomain must be codomain of morphism or None")

            if f.domain() is self.surface():
                parent = Hom(self, codomain)
                return parent.__make_element_class__(
                    SimplicialChainMorphism_induced
                )(parent, f)
        elif isinstance(f, Matrix):
            if codomain is None:
                if f.is_square():
                    codomain = self

            if codomain is None:
                raise NotImplementedError("cannot deduce codomain from this matrix")

            if f.ncols() != self.ngens():
                raise ValueError(
                    "matrix must have one column for each generator of homology"
                )

            if f.nrows() != codomain.ngens():
                raise ValueError(
                    "matrix must have one row for each generator of the codomain"
                )

            parent = Hom(self, codomain)
            return parent.__make_element_class__(SimplicialChainMorphism_matrix)(
                parent, f
            )

        raise NotImplementedError(
            "cannot create a morphism of chains from this data yet"
        )

    def _Hom_(self, Y, category=None):
        r"""
        Return the space of morphisms from these chains to ``Y``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: C = S.chains()

            sage: End(C)
            Endomorphisms of C₁(Translation Surface in H_1(0) built from a square)
        """
        if isinstance(Y, SimplicialChainModule):
            return SimplicialChainMorphismSpace(self, Y, category=category)

        return super()._Hom_(Y, category=category)

    def __eq__(self, other):
        r"""
        Return whether these chains are indistinguishable from ``other``.

        .. NOTE::

            We cannot rely on the builtin `==` by ``id`` since we need to
            detect chains over equal but distinct surfaces to be equal. See
            :meth:`homology` for ideas on how to fix this.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialChains
            sage: T = translation_surfaces.square_torus()
            sage: C = SimplicialChains(T)

            sage: T = translation_surfaces.square_torus()
            sage: CC = SimplicialChains(T)

            sage: C == CC
            True

        """
        if not isinstance(other, SimplicialChainModule):
            return False

        return (
            self._surface == other._surface
            and self._coefficients == other._coefficients
            and self._relative == other._relative
            and self.category() == other.category()
        )

    def __hash__(self):
        r"""
        Return a hash value for these chains that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialChains
            sage: T = translation_surfaces.square_torus()
            sage: C = SimplicialChains(T)

            sage: T = translation_surfaces.square_torus()
            sage: CC = SimplicialChains(T)

            sage: hash(C) == hash(CC)
            True

        """
        return hash(
            (
                self._surface,
                self._coefficients,
                self._relative,
                self.category(),
            )
        )


class SimplicialHomologyMorphismSpace(MorphismSpace):
    r"""
    The space of homomorphisms from the homology ``domain`` to ``codomain``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: M = matrix([[1, 0], [0, 1]])
        sage: f = A.derivative().section()(M, check=False)

        sage: from flatsurf import SimplicialHomology
        sage: H = SimplicialHomology(S)
        sage: g = H.hom(f)
        sage: G = g.parent()

    Since these are homomorphisms in homology, they preserve the linear
    structure of the homology::

        sage: g.category()
        Category of endsets of modules over Integer Ring

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphismSpace
        sage: isinstance(G, SimplicialHomologyMorphismSpace)
        True

        sage: TestSuite(G).run()

    """

    def __init__(self, domain, codomain, category=None):
        from sage.all import Hom

        super().__init__(
            domain,
            codomain,
            category=category or Hom(domain, codomain).homset_category(),
        )

    def _an_element_(self):
        r"""
        Return some homomorphism in homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S.homology()).an_element()
            Generic endomorphism of H₁(Translation Surface in H_1(0) built from a square)
              Defn: [1 0]
                    [0 1]

            sage: Hom(S.homology(), T.homology()).an_element()
            Generic morphism:
              From: H₁(Translation Surface in H_1(0) built from a square)
              To:   H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [0 0]
                    [0 0]
                    [0 0]
                    [0 0]

        """
        if self.domain() is self.codomain():
            return self.identity()
        return self.zero()

    def identity(self):
        r"""
        Return the identity homomorphism in this space (if it exists).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S.homology()).identity()
            Generic endomorphism of H₁(Translation Surface in H_1(0) built from a square)
              Defn: [1 0]
                    [0 1]

        """
        if self.is_endomorphism_set():
            from sage.all import identity_matrix

            matrix = identity_matrix(
                self.codomain().base_ring(), self.domain().ngens(), sparse=True
            )
            return self.__make_element_class__(SimplicialHomologyMorphism_matrix)(
                self, matrix
            )
        return super().identity()

    def zero(self):
        r"""
        Return the zero homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: End(S.homology()).zero()
            Generic endomorphism of H₁(Translation Surface in H_1(0) built from a square)
              Defn: [0 0]
                    [0 0]

        """
        from sage.all import zero_matrix

        return self.domain().hom(
            zero_matrix(
                self.codomain().base_ring(),
                nrows=self.codomain().ngens(),
                ncols=self.domain().ngens(),
                sparse=True,
            ),
            codomain=self.codomain(),
        )

    def base_ring(self):
        r"""
        Return the ring over which these homomorphisms are defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: End(S.homology()).base_ring()
            Integer Ring

        """
        if self.domain().base_ring() is self.codomain().base_ring():
            return self.domain().base_ring()

        return super().base_ring()

    def __reduce__(self):
        r"""
        Return a picklable version of this space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: H = End(S.homology()).base_ring()
            sage: loads(dumps(H)) == H
            True

        """
        return SimplicialHomologyMorphismSpace, (
            self.domain(),
            self.codomain(),
            self.homset_category(),
        )

    def __repr__(self):
        r"""
        Return a printable representation of this space of homomorphisms.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: H = End(S.homology())
            sage: H
            Endomorphisms of H₁(Translation Surface in H_1(0) built from a square)

        """
        if self.domain() == self.codomain():
            return f"Endomorphisms of {self.domain()!r}"
        return f"Homomorphisms from {self.domain()!r} to {self.codomain()!r}"


class SimplicialChainMorphismSpace(MorphismSpace):
    r"""
    The space of homomorphisms from the chains ``domain`` to ``codomain``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.triangulate()

        sage: from flatsurf import SimplicialChains
        sage: C = SimplicialChains(S)
        sage: g = C.hom(f)
        sage: G = g.parent()

    Since these are homomorphisms of chains, they preserve the linear
    structure::

        sage: g.category()
        Category of homsets of modules over Integer Ring

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialChainMorphismSpace
        sage: isinstance(G, SimplicialChainMorphismSpace)
        True

        sage: TestSuite(G).run()

    """

    def __init__(self, domain, codomain, category=None):
        from sage.all import Hom

        super().__init__(
            domain,
            codomain,
            category=category or Hom(domain, codomain).homset_category(),
        )

    def _an_element_(self):
        r"""
        Return some homomorphism of chains.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S.chains()).an_element()
            Generic endomorphism of C₁(Translation Surface in H_1(0) built from a square)
              Defn: [1 0]
                    [0 1]

            sage: Hom(S.chains(), T.chains()).an_element()
            Generic morphism:
              From: C₁(Translation Surface in H_1(0) built from a square)
              To:   C₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [0 0]
                    [0 0]
                    [0 0]
                    [0 0]
                    [0 0]
                    [0 0]

        """
        if self.domain() is self.codomain():
            return self.identity()
        return self.zero()

    def identity(self):
        r"""
        Return the identity homomorphism in this space (if it exists).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S.chains()).identity()
            Generic endomorphism of C₁(Translation Surface in H_1(0) built from a square)
              Defn: [1 0]
                    [0 1]

        """
        if self.is_endomorphism_set():
            from sage.all import identity_matrix

            matrix = identity_matrix(
                self.codomain().base_ring(), self.domain().ngens(), sparse=True
            )
            return self.__make_element_class__(SimplicialChainMorphism_matrix)(
                self, matrix
            )
        return super().identity()

    def zero(self):
        r"""
        Return the zero homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: End(S.chains()).zero()
            Generic endomorphism of C₁(Translation Surface in H_1(0) built from a square)
              Defn: [0 0]
                    [0 0]

        """
        from sage.all import zero_matrix

        return self.domain().hom(
            zero_matrix(
                self.codomain().base_ring(),
                nrows=self.codomain().ngens(),
                ncols=self.domain().ngens(),
                sparse=True,
            ),
            codomain=self.codomain(),
        )

    def base_ring(self):
        r"""
        Return the ring over which these homomorphisms are defined.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: End(S.chains()).base_ring()
            Integer Ring

        """
        if self.domain().base_ring() is self.codomain().base_ring():
            return self.domain().base_ring()

        return super().base_ring()

    def __reduce__(self):
        r"""
        Return a picklable version of this space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: H = End(S.chains())
            sage: loads(dumps(H)) == H
            True

        """
        return SimplicialChainMorphismSpace, (
            self.domain(),
            self.codomain(),
            self.homset_category(),
        )

    def __repr__(self):
        r"""
        Return a printable representation of this space of homomorphisms.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: H = End(S.chains())
            sage: H
            Endomorphisms of C₁(Translation Surface in H_1(0) built from a square)

        """
        if self.domain() == self.codomain():
            return f"Endomorphisms of {self.domain()!r}"
        return f"Homomorphisms from {self.domain()!r} to {self.codomain()!r}"


class SimplicialHomologyMorphism_base(Morphism):
    r"""
    Base class for all homomorphisms in homology.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: M = matrix([[1, 0], [0, 1]])
        sage: f = A.derivative().section()(M, check=False)

        sage: from flatsurf import SimplicialHomology
        sage: H = SimplicialHomology(S)
        sage: g = H.hom(f)

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphism_base
        sage: isinstance(g, SimplicialHomologyMorphism_base)
        True

        sage: TestSuite(g).run()  # optional: pyflatsurf

    """

    @cached_method
    def matrix(self):
        r"""
        Return the matrix describing this homomorphism on the generators of
        homology (as a multiplication from the left).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: H.gens()
            ([(0, 1)], [(0, 0)], [(1, 1)], [(2, 0)])
            sage: g = H.hom(f)

            sage: g.matrix()  # optional: pyflatsurf
            [1 0 0 0]
            [1 1 2 0]
            [0 0 1 0]
            [1 0 0 1]

        """
        from sage.all import matrix

        return matrix(
            [list(self(gen).coefficients()) for gen in self.domain().gens()]
        ).transpose()

    def _add_(self, other):
        r"""
        Return the pointwise sum of this morphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = End(H).one()

            sage: g + h  # optional: pyflatsurf
            Generic endomorphism of H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [2 0 0 0]
                    [1 2 2 0]
                    [0 0 2 0]
                    [1 0 0 2]

        """
        return self.domain().hom(
            self.matrix() + other.matrix(), codomain=self.codomain()
        )

    def _acted_upon_(self, x, self_on_left):
        r"""
        Return the morphism given by pointwise multiplying with ``x``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: (2**1234567 * g).matrix().trace() == 2**1234569  # optional: pyflatsurf
            True

        """
        return self.domain().hom(x * self.matrix(), codomain=self.codomain())

    def _neg_(self):
        r"""
        Return the pointwise negative of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: -g  # optional: pyflatsurf
            Generic endomorphism of H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [-1  0  0  0]
                    [-1 -1 -2  0]
                    [ 0  0 -1  0]
                    [-1  0  0 -1]

        """
        return self.domain().hom(-self.matrix(), codomain=self.codomain())

    def _composition(self, other):
        r"""
        Return the composition of this homomorphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 2)
            sage: U = translation_surfaces.mcmullen_L(1, 1, 1, 3)

            sage: f = S.homology().hom(2 * identity_matrix(4), codomain=T.homology())
            sage: g = T.homology().hom(3 * identity_matrix(4), codomain=U.homology())

            sage: g * f
            Generic morphism:
              From: H₁(Translation Surface in H_2(2) built from 3 squares)
              To:   H₁(Translation Surface in H_2(2) built from 2 squares and a rectangle)
              Defn: [6 0 0 0]
                    [0 6 0 0]
                    [0 0 6 0]
                    [0 0 0 6]

        """
        return other.domain().hom(
            self.matrix() * other.matrix(), codomain=self.codomain()
        )

    def __bool__(self):
        r"""
        Return whether this is not the homommorphism that is zero everywhere.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: bool(g)  # optional: pyflatsurf
            True

            sage: bool(g.parent().zero())
            False

        """
        return bool(self.matrix())


class SimplicialHomologyMorphism_matrix(SimplicialHomologyMorphism_base):
    r"""
    A homomorphism of homology that is given by a matrix that describes the
    homomorphism on the generators.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
        sage: T = translation_surfaces.square_torus()

        sage: f = S.homology().hom(matrix([[1, 2, 3, 4], [5, 6, 7, 8]]), codomain=T.homology())
        sage: f
        Generic morphism:
          From: H₁(Translation Surface in H_2(2) built from 3 squares)
          To:   H₁(Translation Surface in H_1(0) built from a square)
          Defn: [1 2 3 4]
                [5 6 7 8]

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphism_matrix
        sage: isinstance(f, SimplicialHomologyMorphism_matrix)
        True

        sage: TestSuite(f).run()

    """

    def __init__(self, parent, matrix):
        super().__init__(parent)

        if matrix.is_mutable():
            from sage.all import matrix as copy

            matrix = copy(matrix)
            matrix.set_immutable()

        self._matrix = matrix

    def _call_(self, g):
        r"""
        Return the image of ``g`` under this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = translation_surfaces.square_torus()

            sage: f = S.homology().hom(matrix([[1, 2, 3, 4], [5, 6, 7, 8]]), codomain=T.homology())
            sage: [f(gen) for gen in S.homology().gens()]
            [5*[(0, 0)] + [(0, 1)],
             6*[(0, 0)] + 2*[(0, 1)],
             7*[(0, 0)] + 3*[(0, 1)],
             8*[(0, 0)] + 4*[(0, 1)]]

        """
        homology, _, _ = self.codomain()._homology()

        return self.codomain()(homology(self._matrix * g.coefficients()))

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        .. NOTE::

            We cannot implement ``_richcmp_`` since we have to handle elements
            with differing parents here. Also, we cannot implement
            ``__richcmp__`` easily here since we are not in Cython.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = H.hom(g.matrix())  # optional: pyflatsurf
            sage: h == h  # optional: pyflatsurf
            True

        Note that this determines whether two morphisms are indistinguishable,
        not whether they are pointwise the same::

            sage: h == g  # optional: pyflatsurf
            False

        """
        if self is other:
            return True
        if not isinstance(other, SimplicialHomologyMorphism_matrix):
            return False

        return self.parent() == other.parent() and self._matrix == other._matrix

    __ne__ = object.__ne__

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = End(S.homology()).one()
            sage: g = End(S.homology()).one()

            sage: hash(f) == hash(g)
            True

        """
        return hash((self.parent(), self._matrix))

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = End(S.homology()).one()
            sage: f
            Generic endomorphism of H₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 1]

        """
        return repr(self._matrix)


class SimplicialHomologyMorphism_induced(SimplicialHomologyMorphism_base):
    r"""
    A homomorphism of homology induced by a morphism of surfaces.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: A = S.affine_automorphism_group()
        sage: M = matrix([[1, 2], [0, 1]])
        sage: f = A.derivative().section()(M, check=False)

        sage: from flatsurf import SimplicialHomology
        sage: H = SimplicialHomology(S)
        sage: g = H.hom(f)

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialHomologyMorphism_induced
        sage: isinstance(g, SimplicialHomologyMorphism_induced)
        True

        sage: TestSuite(g).run()  # optional: pyflatsurf

    """

    def __init__(self, parent, morphism):
        super().__init__(parent)

        self._morphism = morphism

    def _call_(self, x):
        r"""
        Return the image of the homology class ``x`` under this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)

            sage: H.gens()
            ([(0, 1)], [(0, 0)])
            sage: [g(h) for h in H.gens()]  # optional: pyflatsurf
            [2*[(0, 0)] + [(0, 1)], [(0, 0)]]

        """
        return self._morphism._image_homology(x, codomain=self.codomain())

    def _repr_type(self):
        r"""
        Helper method for :meth:`_repr_` to produce a printable representation
        of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: g._repr_type()
            'Induced'

        """
        return "Induced"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_` to produce a printable representation
        of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: print(g._repr_defn())
            Induced by Affine endomorphism of Translation Surface in H_1(0) built from a square
              Defn: Lift of linear action given by
                    [1 2]
                    [0 1]

        """
        return f"Induced by {self._morphism!r}"

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        .. NOTE::

            We cannot implement ``_richcmp_`` since we have to handle elements
            with differing parents here. Also, we cannot implement
            ``__richcmp__`` easily here since we are not in Cython.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = H.hom(f)

            sage: g == h
            True

        Note that this does not compare homomorphisms pointwise::

            sage: h = H.hom(g.matrix())  # optional: pyflatsurf
            sage: g == h  # optional: pyflatsurf
            False

        """
        if self is other:
            return True

        if not isinstance(other, SimplicialHomologyMorphism_induced):
            return False

        return self.parent() == other.parent() and self._morphism == other._morphism

    __ne__ = object.__ne__

    def __hash__(self):
        r"""
        Return a hash value for this homomorphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: A = S.affine_automorphism_group()
            sage: M = matrix([[1, 2], [0, 1]])
            sage: f = A.derivative().section()(M, check=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: g = H.hom(f)
            sage: h = H.hom(f)
            sage: hash(g) == hash(h)
            True

        """
        return hash((self.parent(), self._morphism))


class SimplicialChainMorphism_base(Morphism):
    r"""
    Base class for all homomorphisms of chains.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.triangulate()

        sage: from flatsurf import SimplicialChains
        sage: C = SimplicialChains(S)
        sage: g = C.hom(f)

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialChainMorphism_base
        sage: isinstance(g, SimplicialChainMorphism_base)
        True

        sage: TestSuite(g).run()

    """

    @cached_method
    def matrix(self):
        r"""
        Return the matrix describing this homomorphism on the generators of
        chains (as a multiplication from the left).

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: C.gens()
            ([(0, 1)], [(0, 0)], [(1, 1)], [(0, 3)], [(2, 0)], [(0, 2)])
            sage: g = C.hom(f)

            sage: g.matrix()
            [ 0  0  0  0 -1  0]
            [ 0  0 -1  0  0  0]
            [ 0  0  0  0  0  1]
            [ 0  0  0  0  0  0]
            [ 1  0  0  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  1  0  0]
            [ 0  0  0  0  0  0]
            [ 0  1  0  0  0  0]

        """
        from sage.all import matrix

        return matrix(
            [list(self(gen).coefficients()) for gen in self.domain().gens()]
        ).transpose()

    def _add_(self, other):
        r"""
        Return the pointwise sum of this morphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)

            sage: g + g
            Generic morphism:
              From: C₁(Translation Surface in H_2(2) built from 3 squares)
              To:   C₁(Triangulation of Translation Surface in H_2(2) built from 3 squares)
              Defn: [ 0  0  0  0 -2  0]
                    [ 0  0 -2  0  0  0]
                    [ 0  0  0  0  0  2]
                    [ 0  0  0  0  0  0]
                    [ 2  0  0  0  0  0]
                    [ 0  0  0  0  0  0]
                    [ 0  0  0  2  0  0]
                    [ 0  0  0  0  0  0]
                    [ 0  2  0  0  0  0]

        """
        return self.domain().hom(
            self.matrix() + other.matrix(), codomain=self.codomain()
        )

    def _acted_upon_(self, x, self_on_left):
        r"""
        Return the morphism given by pointwise multiplying with ``x``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = S.triangulate()
            sage: g = f.section()

            sage: C = S.chains()
            sage: h = C.hom(g*f)

            sage: (2**1234567 * h).matrix().trace() == 6 * 2**1234567
            True

        """
        return self.domain().hom(x * self.matrix(), codomain=self.codomain())

    def _neg_(self):
        r"""
        Return the pointwise negative of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)

            sage: -g
            Generic morphism:
              From: C₁(Translation Surface in H_2(2) built from 3 squares)
              To:   C₁(Triangulation of Translation Surface in H_2(2) built from 3 squares)
              Defn: [ 0  0  0  0  1  0]
                    [ 0  0  1  0  0  0]
                    [ 0  0  0  0  0 -1]
                    [ 0  0  0  0  0  0]
                    [-1  0  0  0  0  0]
                    [ 0  0  0  0  0  0]
                    [ 0  0  0 -1  0  0]
                    [ 0  0  0  0  0  0]
                    [ 0 -1  0  0  0  0]

        """
        return self.domain().hom(-self.matrix(), codomain=self.codomain())

    def _composition(self, other):
        r"""
        Return the composition of this homomorphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = translation_surfaces.mcmullen_L(1, 1, 1, 2)
            sage: U = translation_surfaces.mcmullen_L(1, 1, 1, 3)

            sage: f = S.chains().hom(2 * identity_matrix(6), codomain=T.chains())
            sage: g = T.chains().hom(3 * identity_matrix(6), codomain=U.chains())

            sage: g * f
            Generic morphism:
              From: C₁(Translation Surface in H_2(2) built from 3 squares)
              To:   C₁(Translation Surface in H_2(2) built from 2 squares and a rectangle)
              Defn: [6 0 0 0 0 0]
                    [0 6 0 0 0 0]
                    [0 0 6 0 0 0]
                    [0 0 0 6 0 0]
                    [0 0 0 0 6 0]
                    [0 0 0 0 0 6]

        """
        return other.domain().hom(
            self.matrix() * other.matrix(), codomain=self.codomain()
        )

    def __bool__(self):
        r"""
        Return whether this is not the homommorphism that is zero everywhere.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)

            sage: bool(g)
            True

            sage: bool(g.parent().zero())
            False

        """
        return bool(self.matrix())


class SimplicialChainMorphism_matrix(SimplicialChainMorphism_base):
    r"""
    A homomorphism of chains that is given by a matrix that describes the
    homomorphism on the generators.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
        sage: T = translation_surfaces.square_torus()

        sage: f = S.chains().hom(matrix([[1, 2, 3, 4, 5, 6], [7, 8, 9, 0, 1, 2]]), codomain=T.chains())
        sage: f
        Generic morphism:
          From: C₁(Translation Surface in H_2(2) built from 3 squares)
          To:   C₁(Translation Surface in H_1(0) built from a square)
          Defn: [1 2 3 4 5 6]
                [7 8 9 0 1 2]

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialChainMorphism_matrix
        sage: isinstance(f, SimplicialChainMorphism_matrix)
        True

        sage: TestSuite(f).run()

    """

    def __init__(self, parent, matrix):
        super().__init__(parent)

        if matrix.is_mutable():
            from sage.all import matrix as copy

            matrix = copy(matrix)
            matrix.set_immutable()

        self._matrix = matrix

    def _call_(self, g):
        r"""
        Return the image of ``g`` under this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: T = translation_surfaces.square_torus()

            sage: f = S.chains().hom(matrix([[1, 2, 3, 4, 5, 6], [7, 8, 9, 0, 1, 2]]), codomain=T.chains())
            sage: [f(gen) for gen in S.chains().gens()]
            [7*[(0, 0)] + [(0, 1)],
             8*[(0, 0)] + 2*[(0, 1)],
             9*[(0, 0)] + 3*[(0, 1)],
             4*[(0, 1)],
             [(0, 0)] + 5*[(0, 1)],
             2*[(0, 0)] + 6*[(0, 1)]]

        """
        chains = self.codomain()._chains()

        return self.codomain()(chains(self._matrix * g.coefficients()))

    def __eq__(self, other):
        r"""
        Return whethir this morphism is indistinguishable from ``other``.

        .. NOTE::

            We cannot implement ``_richcmp_`` since we have to handle elements
            with differing parents here. Also, we cannot implement
            ``__richcmp__`` easily here since we are not in Cython.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)
            sage: h = C.hom(g.matrix(), codomain=g.codomain())
            sage: h == h
            True

        Note that this determines whether two morphisms are indistinguishable,
        not whether they are pointwise the same::

            sage: h == g
            False

        """
        if self is other:
            return True

        if not isinstance(other, SimplicialChainMorphism_matrix):
            return False

        return self.parent() == other.parent() and self._matrix == other._matrix

    __ne__ = object.__ne__

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = End(S.chains()).one()
            sage: g = End(S.chains()).one()

            sage: hash(f) == hash(g)
            True

        """
        return hash((self.parent(), self._matrix))

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: f = End(S.chains()).one()
            sage: f
            Generic endomorphism of C₁(Translation Surface in H_2(2) built from 3 squares)
              Defn: [1 0 0 0 0 0]
                    [0 1 0 0 0 0]
                    [0 0 1 0 0 0]
                    [0 0 0 1 0 0]
                    [0 0 0 0 1 0]
                    [0 0 0 0 0 1]

        """
        return repr(self._matrix)


class SimplicialChainMorphism_induced(SimplicialChainMorphism_base):
    r"""
    A homomorphism of chains induced by a morphism of surfaces.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.triangulate()

        sage: from flatsurf import SimplicialChains
        sage: C = SimplicialChains(S)
        sage: g = C.hom(f)

    TESTS::

        sage: from flatsurf.geometry.homology import SimplicialChainMorphism_induced
        sage: isinstance(g, SimplicialChainMorphism_induced)
        True

        sage: TestSuite(g).run()

    """

    def __init__(self, parent, morphism):
        super().__init__(parent)

        self._morphism = morphism

    def _call_(self, x):
        r"""
        Return the image of the chain ``x`` under this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)

            sage: C.gens()
            ([(0, 1)], [(0, 0)])
            sage: [g(h) for h in C.gens()]
            [[((0, 0), 1)], [((0, 0), 0)]]

        """
        return self._morphism._image_chain(x, codomain=self.codomain())

    def _repr_type(self):
        r"""
        Helper method for :meth:`_repr_` to produce a printable representation
        of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)
            sage: g._repr_type()
            'Induced'

        """
        return "Induced"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`_repr_` to produce a printable representation
        of this homomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)
            sage: print(g._repr_defn())
            Induced by Triangulation morphism:
              From: Translation Surface in H_1(0) built from a square
              To:   Triangulation of Translation Surface in H_1(0) built from a square

        """
        return f"Induced by {self._morphism!r}"

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        .. NOTE::

            We cannot implement ``_richcmp_`` since we have to handle elements
            with differing parents here. Also, we cannot implement
            ``__richcmp__`` easily here since we are not in Cython.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)
            sage: h = C.hom(f)

            sage: g == h
            True

        Note that this does not compare homomorphisms pointwise::

            sage: h = C.hom(g.matrix(), codomain=g.codomain())
            sage: g == h
            False

        """
        if self is other:
            return True

        if not isinstance(other, SimplicialChainMorphism_induced):
            return False

        return self.parent() == other.parent() and self._morphism == other._morphism

    __ne__ = object.__ne__

    def __hash__(self):
        r"""
        Return a hash value for this homomorphism that is compatible with
        :meth:`__ne__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()

            sage: from flatsurf import SimplicialChains
            sage: C = SimplicialChains(S)
            sage: g = C.hom(f)
            sage: h = C.hom(f)
            sage: hash(g) == hash(h)
            True

        """
        return hash((self.parent(), self._morphism))


def SimplicialHomology(
    surface,
    k=1,
    coefficients=None,
    relative=None,
    implementation="generic",
    category=None,
):
    r"""
    Return the ``k``-th simplicial homology group of ``surface``.

    INPUT:

    - ``surface`` -- a surface

    - ``k`` -- an integer (default: ``1``)

    - ``coefficients`` -- a ring (default: the integer ring);
      consider the homology with coefficients in this ring

    - ``relative`` -- a set (default: the empty set); if non-empty,
      then relative homology with respect to this set is
      constructed.

    - ``implementation`` -- a string (default: ``"generic"``); the
      algorithm used to compute the homology groups. Currently only
      ``"generic"`` is supported, i.e., the groups are computed
      with the generic homology machinery from SageMath.

    - ``category`` -- a category; if not specified, a category for
      the homology group is chosen automatically depending on
      ``coefficients``.

    TESTS:

    Homology is unique and cached::

        sage: from flatsurf import translation_surfaces, SimplicialHomology
        sage: T = translation_surfaces.square_torus()
        sage: SimplicialHomology(T) is SimplicialHomology(T)
        True

    """
    return surface.homology(k, coefficients, relative, implementation, category)


def SimplicialChains(
    surface,
    k=1,
    coefficients=None,
    relative=None,
    category=None,
):
    r"""
    Return the ``k``-th simplicial chain group of ``surface``.

    INPUT:

    - ``surface`` -- a surface

    - ``k`` -- an integer (default: ``1``)

    - ``coefficients`` -- a ring (default: the integer ring);
      consider the homology with coefficients in this ring

    - ``relative`` -- a set (default: the empty set); if non-empty,
      then the chains modulo chains in this set are constructed.

    - ``category`` -- a category; if not specified, a category for
      the homology group is chosen automatically depending on
      ``coefficients``.

    TESTS:

    Chains are unique and cached::

        sage: from flatsurf import translation_surfaces, SimplicialChains
        sage: T = translation_surfaces.square_torus()
        sage: SimplicialChains(T) is SimplicialChains(T)
        True

    """
    return surface.chains(k, coefficients, relative, category)
