# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2016-2019 Vincent Delecroix
#                     2016-2019 W. Patrick Hooper
#                          2023 Julian Rüth
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
from sage.groups.group import Group
from sage.categories.groups import Groups
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import MultiplicativeGroupElement
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.rings.integer_ring import ZZ

from flatsurf.geometry.origami import AbstractOrigami

_Q = QuaternionAlgebra(-1, -1)
_i, _j, _k = _Q.gens()


class MegaWollmilchsauGroupElement(MultiplicativeGroupElement):
    @staticmethod
    def quat_to_tuple(r):
        r"""Convert an element in the quaternion algebra to a quadruple"""
        if isinstance(r, int):
            return (r, 0, 0, 0)
        else:
            return (r[0], r[1], r[2], r[3])

    @staticmethod
    def wedge(r1, r2):
        r"""Wedge two quaterions. Returns an integer."""
        x = MegaWollmilchsauGroupElement.quat_to_tuple(r1)
        y = MegaWollmilchsauGroupElement.quat_to_tuple(r2)
        return -x[0] * y[3] + x[1] * y[2] - x[2] * y[1] + x[3] * y[0]

    def __init__(self, parent, i, r, q):
        if parent is None:
            raise ValueError("The parent must be provided")
        # I should assert that the element lives in the domain of the group.
        if i not in ZZ:
            raise ValueError
        if r not in _Q:
            raise ValueError
        # Actually q should be in {1,-1,-i,i,j,-j,k,-k}. I'm not testing for that.
        if q not in _Q:
            raise ValueError
        # There is one more condition. The group doesn't have full image...
        self._i = i
        self._r = r
        self._q = q
        self._parent = parent
        MultiplicativeGroupElement.__init__(self, parent)

    def _repr_(self):
        return "[" + str(self._i) + ", " + str(self._r) + ", " + str(self._q) + "]"

    def _cmp_(self, other):
        return (
            (self._i > other._i - self._i < other._i)
            or (self._r > other._r - self._r < other._r)
            or (self._q > other._q - self._q < other._q)
        )

    def _mul_(self, m):
        return MegaWollmilchsauGroupElement(
            self._parent,
            self._i
            + m._i
            + MegaWollmilchsauGroupElement.wedge(self._r, self._q * m._r),
            self._r + self._q * m._r,
            self._q * m._q,
        )

    def __invert__(self):
        q1 = ~self._q
        r1 = -(q1 * self._r)
        i1 = -(self._i + MegaWollmilchsauGroupElement.wedge(r1, q1 * self._r))
        return MegaWollmilchsauGroupElement(self._parent, i1, r1, q1)

    def _div_(self, m):
        return self._mul_(m.__invert__())

    def __hash__(self):
        return (
            67 * hash(self._i)
            + 23 * hash(MegaWollmilchsauGroupElement.quat_to_tuple(self._r))
            - 17 * hash(MegaWollmilchsauGroupElement.quat_to_tuple(self._q))
        )


class MegaWollmilchsauGroup(UniqueRepresentation, Group):
    Element = MegaWollmilchsauGroupElement

    def _element_constructor_(self, *args, **kwds):
        if len(args) != 1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        return self.element_class(self, x, **kwds)

    def __init__(self):
        Group.__init__(self, category=Groups().Infinite())

    def _repr_(self):
        return "MegaWollmilchsauGroup"

    def a(self):
        return self.element_class(self, 0, 1, _i)

    def b(self):
        return self.element_class(self, 0, 1, _j)

    def one(self):
        return self.element_class(self, 0, 0, 1)

    def gens(self):
        return (self.a(), self.b())

    def is_abelian(self):
        return False

    def _an_element_(self):
        return self.a()

    def some_elements(self):
        return [self.a(), self.b()]

    def _test_relations(self, **options):
        tester = self._tester(**options)
        a, b = self.gens()
        e = self.one()
        tester.assertEqual(a**4, e)
        tester.assertEqual(b**4, e)
        tester.assertEqual((a * b) ** 4, e)
        tester.assertEqual((a / b) ** 4, e)
        tester.assertEqual((a * a * b) ** 4, e)
        tester.assertEqual((a * a / b) ** 4, e)
        tester.assertNotEqual((a * b / a / b) ** 2, e)


class MegaWollmilchsau(AbstractOrigami):
    def __init__(self):
        self._G = self._domain = MegaWollmilchsauGroup()
        self._a, self._b = self._G.gens()
        self._ai = ~self._a
        self._bi = ~self._b

    def up(self, label):
        return self._b * label

    def down(self, label):
        return self._bi * label

    def right(self, label):
        return self._a * label

    def left(self, label):
        return self._ai * label

    def _repr_(self):
        return "MegaWollmilchsau Origami"
