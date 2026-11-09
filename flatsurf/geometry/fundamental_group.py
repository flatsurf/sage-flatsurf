# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
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

from sage.misc.cachefunc import cached_method

from sage.structure.element import parent

from sage.categories.groups import Groups
from sage.groups.group import Group
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation


def intersection(i0, j0, i1, j1):
    r"""
    Intersection inside a polygon.

    In case of equality we consider that i0 < i1 < j1 < j0

    INPUT:

    - ``i0``, ``j0`` -- start and end of the first curve

    - ``i1``, ``j1`` -- start and end of the second curve

    TESTS::

        sage: from flatsurf.geometry.fundamental_group import intersection
        sage: intersection(3,0,3,2)
        1
        sage: intersection(0,1,1,2)
        1
        sage: intersection(0,2,3,1)
        -1
        sage: intersection(0,3,3,2)
        0
        sage: intersection(1,1,1,1)
        0
        sage: intersection(3,2,3,2)
        0
        sage: intersection(3,2,0,3)
        0
        sage: intersection(0,3,0,3)
        0
        sage: intersection(0,2,0,2)
        0
    """
    if i0 <= j0:  # that is i0 < j0 or i0 == j0
        return (i0 <= i1 <= j0) - (i0 <= j1 <= j0)
    else:
        return (j0 < j1 < i0) - (j0 < i1 < i0)


class Path(MultiplicativeGroupElement):
    r"""
    A closed walk in the dual graph of a surface, an element of a
    :class:`FundamentalGroup`.

    The walk is stored as the sequence of its letters: the oriented edges
    ``(label, edge)`` through which it enters its successive polygons. It is
    freely reduced at construction.

    INPUT:

    - ``parent`` -- a :class:`FundamentalGroup`

    - ``walk`` -- a sequence of pairs ``(label, edge)``

    - ``check`` -- a boolean (default: ``True``); whether to check that
      ``walk`` is a closed walk at the base polygon of ``parent``

    - ``reduced`` -- a boolean (default: ``False``); whether ``walk`` is
      known to be freely reduced

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: t = translation_surfaces.square_torus()
        sage: F = t.fundamental_group(punctured_vertices=True)
        sage: a, b = F.gens()
        sage: (a * b)._walk
        ((0, 2), (0, 3))

    TESTS::

        sage: F([(0, 2), (0, 0)]).is_one()
        True

        sage: S = translation_surfaces.veech_double_n_gon(5)
        sage: F = S.fundamental_group(punctured_vertices=True)
        sage: F([(0, 1)])
        Traceback (most recent call last):
        ...
        ValueError: the walk does not start at the base polygon 0

    A letter that is not a pair is refused::

        sage: from flatsurf.geometry.fundamental_group import Path
        sage: Path(F, [(1, 0), None])
        Traceback (most recent call last):
        ...
        ValueError: the letter None is not a pair (label, edge)
    """

    # activating the following somehow break the discovery of the Python _mul_
    # method below...
    #    __slots__ = ['_base_label', '_walk']

    def __init__(self, parent, walk, check=True, reduced=False):
        self._base_label = parent._b
        letters = []
        for letter in walk:
            try:
                letters.append(tuple(letter))
            except TypeError:
                if not check:
                    raise
                raise ValueError(
                    "the letter {} is not a pair (label, edge)".format(letter)
                ) from None
        self._walk = tuple(letters)
        MultiplicativeGroupElement.__init__(self, parent)

        if check:
            self._check()
        if not reduced:
            self._reduce()

    def _reduce(self):
        r"""
        Freely reduce the walk of this element.

        A letter followed by its partner, that is, the walk crossing an edge
        and immediately crossing it back, cancels.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group(punctured_vertices=True)
            sage: F([0, 1, 3, 2]).is_one()
            True
        """
        surface = self.parent()._s
        stack = []
        for letter in self._walk:
            if stack and surface.opposite_edge(*letter) == stack[-1]:
                stack.pop()
            else:
                stack.append(letter)
        self._walk = tuple(stack)

    @property
    def _polys(self):
        r"""
        Return the labels of the polygons visited by this walk, starting
        and ending at the base polygon.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a, b = t.fundamental_group(punctured_vertices=True).gens()
            sage: (a * b)._polys
            (0, 0, 0)
        """
        return (self._base_label,) + tuple(label for label, _ in self._walk)

    @property
    def _edges(self):
        r"""
        Return the edges through which this walk leaves its successive
        polygons.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a, b = t.fundamental_group(punctured_vertices=True).gens()
            sage: (a * b)._edges
            (0, 1)
        """
        surface = self.parent()._s
        return tuple(surface.opposite_edge(*letter)[1] for letter in self._walk)

    @property
    def _edges_rev(self):
        r"""
        Return the edges through which this walk enters its successive
        polygons.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a, b = t.fundamental_group(punctured_vertices=True).gens()
            sage: (a * b)._edges_rev
            (2, 3)
        """
        return tuple(edge for _, edge in self._walk)

    def _poly_cross_dict(self):
        r"""
        Return, for each polygon, the pairs ``(entry, exit)`` of the edges
        through which this walk enters and leaves it.

        The empty walk of the identity crosses no polygon.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group(punctured_vertices=True)
            sage: a, b = F.gens()
            sage: a._poly_cross_dict()
            {0: [(2, 0)]}
            sage: F.one()._poly_cross_dict()
            {0: []}
        """
        d = {p: [] for p in self.parent()._s.labels()}
        if not self._walk:
            return d
        d[self._polys[0]].append((self._edges_rev[-1], self._edges[0]))
        for i in range(1, len(self._polys) - 1):
            p = self._polys[i]
            e0 = self._edges_rev[i - 1]
            e1 = self._edges[i]
            d[p].append((e0, e1))
        return d

    def __hash__(self):
        r"""
        Return a hash value for this element.

        Only the elements of a group punctured at every vertex are hashable:
        otherwise distinct reduced walks may be equal in the group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a, b = t.fundamental_group(punctured_vertices=True).gens()
            sage: hash(a) == hash(a * b * ~b)
            True

        TESTS::

            sage: a, b = t.fundamental_group().gens()
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError: unhashable element of a fundamental group not punctured at every vertex
        """
        if not self.parent()._all_punctured:
            raise TypeError(
                "unhashable element of a fundamental group not punctured at every vertex"
            )
        return hash((self._base_label, self._walk))

    def __eq__(self, other):
        r"""
        Return whether this element is equal to ``other``.

        Two elements with the same reduced walk are equal. When every vertex
        is punctured, the group is free and elements with distinct reduced
        walks differ; otherwise deciding the equality needs a homotopy test,
        which is not implemented.

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group(punctured_vertices=True)
            sage: a,b = F.gens()
            sage: a == b
            False
            sage: a*b == b*a
            False
            sage: a*b == a*b
            True

        Without punctures, only the equality of reduced walks is decided::

            sage: F = t.fundamental_group()
            sage: a, b = F.gens()
            sage: a*b == a*b
            True
            sage: F([0, 2]) == F.one()
            True
            sage: a*b == b*a
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot decide the equality of elements with distinct reduced walks in a fundamental group not punctured at every vertex

        Elements of different groups are never equal::

            sage: a == t.fundamental_group(punctured_vertices=True).gens()[0]
            False
        """
        if parent(self) is not parent(other):
            return False
        if self._walk == other._walk:
            return True
        if self.parent()._all_punctured:
            return False
        raise NotImplementedError(
            "cannot decide the equality of elements with distinct reduced walks in a fundamental group not punctured at every vertex"
        )

    def __ne__(self, other):
        r"""
        Return whether this element is not equal to ``other``.

        See :meth:`__eq__` for the cases that are decided.

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group(punctured_vertices=True)
            sage: a,b = F.gens()
            sage: a != b
            True
            sage: a*b != b*a
            True
            sage: a*b != a*b
            False

        Without punctures, only the equality of reduced walks is decided::

            sage: F = t.fundamental_group()
            sage: a, b = F.gens()
            sage: a*b != a*b
            False
            sage: a*b != b*a
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot decide the equality of elements with distinct reduced walks in a fundamental group not punctured at every vertex
        """
        return not self.__eq__(other)

    def _check(self):
        r"""
        Check that the walk of this element is a closed walk at the base
        polygon.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5)
            sage: F = S.fundamental_group(punctured_vertices=True)
            sage: F([(1, 0), (0, 1)])
            0 --0-- 1 --1-- 0
            sage: F([(1, 0), (1, 1)])
            Traceback (most recent call last):
            ...
            ValueError: the letter (1, 1) does not start at the polygon 1 where the walk is
            sage: F([(1, 0)])
            Traceback (most recent call last):
            ...
            ValueError: the walk does not end at the base polygon 0
        """
        surface = self.parent()._s
        current = self._base_label
        for position, letter in enumerate(self._walk):
            if len(letter) != 2:
                raise ValueError(
                    "the letter {} is not a pair (label, edge)".format(letter)
                )
            partner = surface.opposite_edge(*letter)
            if partner is None:
                raise ValueError("the edge {} is not glued".format(letter))
            if partner[0] != current:
                if position == 0:
                    raise ValueError(
                        "the walk does not start at the base polygon {}".format(
                            self._base_label
                        )
                    )
                raise ValueError(
                    "the letter {} does not start at the polygon {} where the walk is".format(
                        letter, current
                    )
                )
            current = letter[0]
        if current != self._base_label:
            raise ValueError(
                "the walk does not end at the base polygon {}".format(self._base_label)
            )

    def is_one(self):
        r"""
        Return whether this element is the identity, that is, whether its
        reduced walk is empty.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a, b = t.fundamental_group().gens()
            sage: (a * ~a).is_one()
            True
            sage: a.is_one()
            False
        """
        return not self._walk

    def _repr_(self):
        return "".join(
            "{} --{}-- ".format(p, e) for p, e in zip(self._polys, self._edges)
        ) + "{}".format(self._polys[-1])

    def _mul_(self, other):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a,b = t.fundamental_group(punctured_vertices=True).gens()
            sage: a*b
            0 --0-- 0 --1-- 0
        """
        P = self.parent()
        return P.element_class(P, self._walk + other._walk, check=False)

    def __invert__(self):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: o = translation_surfaces.octagon_and_squares()
            sage: F = o.fundamental_group(punctured_vertices=True)
            sage: a1,a2,a3,a4,a5,a6 = F.gens()
            sage: (a1 * a2 * ~a2 * ~a1).is_one()
            True
            sage: (a1 * ~a2 * a2 * a1) == a1 * a1
            True
        """
        P = self.parent()
        return P.element_class(
            P, P._inverse_walk(self._walk), check=False, reduced=True
        )

    def intersection(self, other):
        r"""
        The intersection number of this element with ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: a,b = t.fundamental_group(punctured_vertices=True).gens()
            sage: a.intersection(b)
            1
            sage: b.intersection(a)
            -1

        The identity intersects every element zero times::

            sage: F = a.parent()
            sage: F.one().intersection(a), a.intersection(F.one())
            (0, 0)

        This is an Abelian invariant::

            sage: x1 = a*b*b*~a*~a*b*b*a*~b*~b
            sage: x2 = a*b*a*b
            sage: x3 = ~b*~b*a
            sage: (x1*x2).intersection(x3) == x1.intersection(x3) + x2.intersection(x3)
            True
            sage: (x2*x1*~x2).intersection(x2*x3*~x2) == x1.intersection(x3)
            True

        A little bit more involved example::

            sage: S = SymmetricGroup(4)
            sage: r = S('(1,2)(3,4)')
            sage: u = S('(2,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: F = o.fundamental_group(punctured_vertices=True)
            sage: x = F([1,2,2,3])
            sage: x.intersection(x)
            0

            sage: a = F.gens()
            sage: m = matrix(ZZ, 5, lambda i,j: a[i].intersection(a[j]))
            sage: m
            [ 0  1  0  0  0]
            [-1  0 -1  0  0]
            [ 0  1  0  1  0]
            [ 0  0 -1  0 -1]
            [ 0  0  0  1  0]

        A slightly more involved example::

            sage: S = SymmetricGroup(8)
            sage: r = S('(1,2,3,4,5,6,7,8)')
            sage: u = S('(1,8,5,4)(2,3)(6,7)')
            sage: o = translation_surfaces.origami(r,u)
            sage: a = o.fundamental_group(punctured_vertices=True).gens()
            sage: m = matrix(ZZ, 9, lambda i,j: a[i].intersection(a[j]))
            sage: m
            [ 0 -1  1 -1  1 -1 -1 -1  1]
            [ 1  0  1  0  0 -1 -1 -1  1]
            [-1 -1  0  0  0  0  0  0  0]
            [ 1  0  0  0 -1  0  0  0  0]
            [-1  0  0  1  0  0  0  0  0]
            [ 1  1  0  0  0  0  0  0  0]
            [ 1  1  0  0  0  0  0  1 -1]
            [ 1  1  0  0  0  0 -1  0 -1]
            [-1 -1  0  0  0  0  1  1  0]
        """
        si = self._poly_cross_dict()
        oi = other._poly_cross_dict()
        n = 0
        for p in self.parent()._s.labels():
            for e0, e1 in si[p]:
                for f0, f1 in oi[p]:
                    n += intersection(e0, e1, f0, f1)
        return n


class FundamentalGroup(UniqueRepresentation, Group):
    r"""
    The fundamental group of a surface built from polygons, punctured at
    some of its vertices.

    The elements are closed walks in the dual graph of the surface, based at
    a polygon: the vertices of the dual graph are the polygons, its edges
    the gluings, and its faces the vertices of the surface. These walks
    avoid the vertices of the surface; up to free reduction they form the
    fundamental group of the surface punctured at all of its vertices, a
    free group on :meth:`gens`. Filling a vertex adds its loop (see
    :meth:`vertex_loops`) as a relation.

    A walk is the sequence of its letters: the oriented edges ``(label,
    edge)`` through which it enters its successive polygons. An oriented
    edge runs counterclockwise around its polygon, and its partner is the
    edge it is glued to. A walk starts at the polygon of the partner of its
    first letter. The product ``u * v`` is the walk of ``u`` followed by the
    walk of ``v``; the inverse reverses the walk and replaces each letter by
    its partner. A letter followed by its partner cancels.

    When every vertex is punctured, two elements are equal if and only if
    their reduced walks agree, and the elements are hashable. Otherwise,
    ``==`` answers ``True`` when the reduced walks agree and raises a
    ``NotImplementedError`` when they differ, since deciding equality
    then needs a homotopy test; the elements are not hashable.

    INPUT:

    - ``surface`` -- an oriented surface built from finitely many polygons

    - ``base`` -- the label of the base polygon

    - ``punctured_vertices`` -- a boolean or a list of vertices of
      ``surface`` (default: ``False``); the vertices at which the surface is
      punctured, ``True`` standing for all of them

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: t = translation_surfaces.square_torus()
        sage: F = t.fundamental_group(punctured_vertices=True)
        sage: F
        Fundamental group of Translation Surface in H_1(0) built from a square punctured at all of its vertices based at polygon 0
        sage: TestSuite(F).run()

    Without punctures, the equality of elements is only partially decided.
    The default test suite passes because ``some_elements()`` is a single
    element here, which is only compared with itself::

        sage: F = t.fundamental_group()
        sage: F
        Fundamental group of Translation Surface in H_1(0) built from a square based at polygon 0
        sage: list(F.some_elements())
        [0 --0-- 0 --1-- 0]
        sage: TestSuite(F).run()

    On the generators, the three hooks below compare distinct reduced walks,
    whose equality raises a ``NotImplementedError``, and are skipped::

        sage: TestSuite(F).run(elements=list(F.gens()),
        ....:     skip=['_test_elements_eq_symmetric', '_test_elements_eq_transitive',
        ....:           '_test_elements_neq'])

    TESTS::

        sage: from flatsurf.geometry.fundamental_group import FundamentalGroup
        sage: FundamentalGroup(t, 0, True) is FundamentalGroup(t, 0, list(t.vertices()))
        True
        sage: FundamentalGroup(t, 0, [t(0, 0)])
        Fundamental group of Translation Surface in H_1(0) built from a square punctured at all of its vertices based at polygon 0
        sage: FundamentalGroup(t, 0, [0])
        Traceback (most recent call last):
        ...
        ValueError: 0 is not a vertex of the surface
        sage: t.fundamental_group(base_label=5)
        Traceback (most recent call last):
        ...
        ValueError: 5 is not a label of a polygon of the surface
    """

    Element = Path

    @staticmethod
    def __classcall__(cls, surface, base, punctured_vertices=False):
        r"""
        Normalize the arguments so that equal groups are identical.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.fundamental_group import FundamentalGroup
            sage: t = translation_surfaces.square_torus()
            sage: FundamentalGroup(t, 0) is FundamentalGroup(t, 0, [])
            True
        """
        if not surface.is_finite_type():
            raise ValueError("the method only work for finite surfaces")
        if base not in surface.labels():
            raise ValueError(
                "{} is not a label of a polygon of the surface".format(base)
            )
        vertices = surface.vertices()
        if punctured_vertices is True:
            punctured_vertices = frozenset(vertices)
        elif punctured_vertices is False:
            punctured_vertices = frozenset()
        else:
            punctured_vertices = frozenset(punctured_vertices)
            for vertex in punctured_vertices:
                if vertex not in vertices:
                    raise ValueError("{} is not a vertex of the surface".format(vertex))
        return super().__classcall__(cls, surface, base, punctured_vertices)

    def __init__(self, surface, base, punctured_vertices):
        self._s = surface
        self._b = base
        self._punctured_vertices = punctured_vertices
        self._all_punctured = len(punctured_vertices) == len(surface.vertices())

        Group.__init__(self, category=Groups().Infinite())

    def _element_constructor_(self, *args):
        r"""
        Return the element given by a walk from the base polygon.

        The walk is a sequence whose entries are either the index of the
        edge through which the walk leaves its current polygon, or a letter
        ``(label, edge)``, the edge through which it enters the next polygon.

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: S = SymmetricGroup(4)
            sage: r = S('(1,2)(3,4)')
            sage: u = S('(2,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: F = o.fundamental_group(punctured_vertices=True)
            sage: F([1,1])
            1 --1-- 2 --1-- 1
            sage: F([1,2,2,3])
            1 --1-- 2 --2-- 3 --2-- 2 --3-- 1
            sage: F([(2, 3), (1, 3)]) == F([1, 1])
            True

            sage: F([1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: the walk does not end at the base polygon 1

        An element of another fundamental group is not converted, and a
        step that is neither an edge index nor a pair is refused::

            sage: G = o.fundamental_group(base_label=2, punctured_vertices=True)
            sage: F(G.gens()[0])
            Traceback (most recent call last):
            ...
            ValueError: the element 2 --3-- 1 --0-- 1 --1-- 2 belongs to another fundamental group
            sage: F([(2, 3), None])
            Traceback (most recent call last):
            ...
            ValueError: the letter None is not a pair (label, edge)
        """
        if len(args) == 1 and isinstance(args[0], Path):
            raise ValueError(
                "the element {} belongs to another fundamental group".format(args[0])
            )
        if len(args) == 1 and isinstance(args[0], (tuple, list)):
            args = args[0]
        s = self._s
        current = self._b
        walk = []
        for step in args:
            if isinstance(step, (tuple, list)):
                letter = tuple(step)
            else:
                try:
                    step = int(step)
                except (TypeError, ValueError):
                    raise ValueError(
                        "the letter {} is not a pair (label, edge)".format(step)
                    ) from None
                step = step % len(s.polygon(current).vertices())
                letter = s.opposite_edge(current, step)
                if letter is None:
                    raise ValueError("the edge {} is not glued".format((current, step)))
            walk.append(letter)
            current = letter[0]
        return self.element_class(self, walk)

    def _repr_(self):
        r"""
        Return a printable representation of this group.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5)
            sage: S.fundamental_group(punctured_vertices=[S(0, 0)])
            Fundamental group of Translation Surface in H_2(2) built from 2 regular pentagons punctured at all of its vertices based at polygon 0
        """
        punctures = ""
        if self._all_punctured:
            punctures = " punctured at all of its vertices"
        elif self._punctured_vertices:
            punctures = " punctured at {} of its vertices".format(
                len(self._punctured_vertices)
            )
        return "Fundamental group of {}{} based at polygon {}".format(
            self._s, punctures, self._b
        )

    @cached_method
    def one(self):
        r"""
        Return the identity element, the empty walk at the base polygon.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: t.fundamental_group().one()
            0
        """
        return self.element_class(self, (), check=False, reduced=True)

    def _inverse_walk(self, walk):
        r"""
        Return the reverse of ``walk`` with each letter replaced by its
        partner.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group()
            sage: F._inverse_walk(((0, 2), (0, 3)))
            ((0, 1), (0, 0))
        """
        return tuple(self._s.opposite_edge(*letter) for letter in reversed(walk))

    @cached_method
    def _spanning_tree(self):
        r"""
        Return a spanning tree of the dual graph rooted at the base polygon,
        and the gluings off it.

        OUTPUT:

        A pair ``(tree, cotree)``. The dictionary ``tree`` maps the base
        label to ``(None, None, None)`` and every other label ``p`` to
        ``(q, e, f)``: the walk from the base reaches ``p`` from its parent
        ``q``, leaving ``q`` through ``f`` and entering ``p`` through ``e``.
        The list ``cotree`` holds a tuple ``(p1, e1, p2, e2)`` for each
        gluing of the edge ``e1`` of ``p1`` with the edge ``e2`` of ``p2``
        that is not in the tree, with ``(p1, e1) <= (p2, e2)``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: t.fundamental_group()._spanning_tree()
            ({0: (None, None, None)}, [(0, 1, 0, 3), (0, 0, 0, 2)])
        """
        p = self._b
        s = self._s
        tree = {}  # a tree whose root is base_label
        cotree = []

        tree[p] = (None, None, None)

        wait = []  # list of edges of the dual graph, ie p1 -- (e1,e2) --> p2
        for e in range(len(s.polygon(p).vertices())):
            opposite = s.opposite_edge(p, e)
            if opposite is not None:
                wait.append(opposite + (p, e))
        while wait:
            p1, e1, p2, e2 = wait.pop()
            assert p2 in tree
            if p1 in tree:  # new cycle?
                if (p1, e1) > (p2, e2):
                    continue
                cotree.append((p1, e1, p2, e2))

            else:  # new branch
                tree[p1] = (p2, e1, e2)
                for e in range(len(s.polygon(p1).vertices())):
                    if e != e1:
                        opposite = s.opposite_edge(p1, e)
                        if opposite is not None:
                            wait.append(opposite + (p1, e))

        return tree, cotree

    def _tree_walk(self, label):
        r"""
        Return the walk in the spanning tree from the base polygon to the
        polygon ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.veech_double_n_gon(5)
            sage: S.fundamental_group()._tree_walk(1)
            ((1, 4),)
        """
        tree = self._spanning_tree()[0]
        walk = []
        parent_label, edge, _ = tree[label]
        while parent_label is not None:
            walk.append((label, edge))
            label = parent_label
            parent_label, edge, _ = tree[label]
        return tuple(reversed(walk))

    @cached_method
    def gens(self):
        r"""
        Return a generating set: one generator for each gluing off a
        spanning tree of the dual graph.

        With every vertex punctured, the group is free on these generators.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = SymmetricGroup(8)
            sage: r = S('(1,2,3,4,5,6,7,8)')
            sage: u = S('(1,8,5,4)(2,3)(6,7)')
            sage: o = translation_surfaces.origami(r,u)
            sage: len(o.fundamental_group(punctured_vertices=True).gens())
            9

        Unglued edges are not crossed by any generator::

            sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)]))
            0
            sage: S.glue((0, 1), (0, 3))
            sage: S.set_immutable()
            sage: S.fundamental_group().gens()
            (0 --1-- 0,)
        """
        tree, cotree = self._spanning_tree()
        basis = []
        for p1, e1, p2, e2 in cotree:
            walk = (
                self._tree_walk(p1)
                + ((p2, e2),)
                + self._inverse_walk(self._tree_walk(p2))
            )
            basis.append(self.element_class(self, walk, check=False))

        basis.sort(key=lambda g: (g._polys, g._edges, g._edges_rev))
        return tuple(basis)

    def vertex_loops(self):
        r"""
        Return the loops around the vertices of the surface.

        For a vertex ``v``, let ``(q_0, j_0)`` be its first corner in the
        order of the edges of the surface, ``v`` being the start of the edge
        ``j_0`` of ``q_0``. The loop around ``v`` has one letter per corner,
        counterclockwise around ``v``: ``(q_0, j_0), (q_1, j_1), ...,
        (q_{k-1}, j_{k-1})``, where ``(q_{i+1}, j_{i+1})`` is the edge glued
        to the edge ``j_i - 1`` of ``q_i``. It is a walk at ``q_{k-1}``,
        which is conjugated into the group by the walk in the spanning tree
        from the base polygon to ``q_{k-1}``.

        A vertex on the boundary of the surface has no loop and is left out.

        OUTPUT:

        A dictionary from the vertices of the surface (as returned by
        ``vertices()``) to elements of this group.

        EXAMPLES:

        On the square torus, the loop around the vertex is a commutator of
        the generators::

            sage: from flatsurf import translation_surfaces
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group(punctured_vertices=True)
            sage: a, b = F.gens()
            sage: F.vertex_loops()
            {Vertex 0 of polygon 0: 0 --2-- 0 --3-- 0 --0-- 0 --1-- 0}
            sage: F.vertex_loops()[t(0, 0)] == ~a * ~b * a * b
            True

        On the square glued by rotations of angle π/2, the loops around the
        three vertices multiply to the identity::

            sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
            sage: P = Polygon(vertices=[(0,0), (1,0), (1,1), (0,1)])
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(P, label=0)
            0
            sage: S.glue((0, 0), (0, 1))
            sage: S.glue((0, 2), (0, 3))
            sage: S.set_immutable()
            sage: F = S.fundamental_group(punctured_vertices=True)
            sage: x, y = F.gens()
            sage: loops = F.vertex_loops()
            sage: loops[S(0, 1)] == x
            True
            sage: loops[S(0, 3)] == y
            True
            sage: loops[S(0, 0)] == ~x * ~y
            True
            sage: (loops[S(0, 3)] * loops[S(0, 1)] * loops[S(0, 0)]).is_one()
            True
        """
        s = self._s
        loops = {}
        seen = set()
        for q0, j0 in s.edges():
            if (q0, j0) in seen:
                continue
            letters = [(q0, j0)]
            seen.add((q0, j0))
            boundary = s.opposite_edge(q0, j0) is None
            while not boundary:
                q, j = letters[-1]
                letter = s.opposite_edge(q, (j - 1) % len(s.polygon(q).vertices()))
                if letter is None:
                    boundary = True
                elif letter == (q0, j0):
                    break
                else:
                    letters.append(letter)
                    seen.add(letter)
            if boundary:
                continue
            # pylint: disable-next=not-callable
            vertex = s(q0, j0)
            t = self._tree_walk(letters[-1][0])
            walk = t + tuple(letters) + self._inverse_walk(t)
            loops[vertex] = self.element_class(self, walk, check=False)
        return loops
