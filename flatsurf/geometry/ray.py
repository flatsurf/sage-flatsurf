r"""
Geometry with rays in the Euclidean plane.

EXAMPLES::

    sage: from flatsurf.geometry.ray import Rays
    sage: R = Rays(QQ)
    sage: R((1, 0))
    Ray towards (1, 0)

    sage: R((1, 0)) == R((2, 0))
    True

    sage: R((1, 0)) == R((-2, 0))
    False

"""
######################################################################
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
######################################################################
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method


class Ray(Element):
    r"""
    A ray in the Euclidean plane.

    EXAMPLES::

        sage: from flatsurf.geometry.ray import Rays
        sage: R = Rays(QQ)
        sage: r = R((1, 0)); r
        Ray towards (1, 0)

        sage: R((0, 0))
        Traceback (most recent call last):
        ...
        ValueError: direction must not be the zero vector

    TESTS::

        sage: from flatsurf.geometry.ray import Ray
        sage: isinstance(r, Ray)
        True
        sage: TestSuite(r).run()

    """

    def __init__(self, parent, direction):
        super().__init__(parent)

        direction = parent.ambient_space()(direction)
        direction.set_immutable()

        if direction.is_zero():
            raise ValueError("direction must not be the zero vector")

        self._direction = direction

    def vector(self):
        r"""
        Return a vector in this ray.

        EXAMPLES::

            sage: from flatsurf.geometry.ray import Rays
            sage: R = Rays(QQ)
            sage: r = R((1, 0)); r
            Ray towards (1, 0)
            sage: r.vector()
            (1, 0)

        """
        return self._direction

    def _neg_(self):
        return self.parent()(-self._direction)

    def _repr_(self):
        return f"Ray towards {self._direction}"

    def _richcmp_(self, other, op):
        r"""
        Return how this ray compares to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.ray import Rays
            sage: R = Rays(QQ)

            sage: R((0, 1)) == R((0, 2))
            True
            sage: R((1, 0)) == R((2, 0))
            True
            sage: R((1, 1)) == R((2, 2))
            True
            sage: R((1, 1)) == R((2, 1))
            False
            sage: R((0, 1)) == R((0, -2))
            False
            sage: R((1, 1)) == R((-2, -2))
            False

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            from sage.all import sgn

            return (
                self._direction[0] * other._direction[1]
                == other._direction[0] * self._direction[1]
                and sgn(self._direction[0]) == sgn(other._direction[0])
                and sgn(self._direction[1]) == sgn(other._direction[1])
            )


class Rays(UniqueRepresentation, Parent):
    r"""
    The space of rays from the origin in the Euclidean plane.

    EXAMPLES::

        sage: from flatsurf.geometry.ray import Rays
        sage: R = Rays(QQ)
        sage: R
        Rays in Vector space of dimension 2 over Rational Field

    TESTS::

        sage: isinstance(R, Rays)
        True
        sage: TestSuite(R).run()

    """
    Element = Ray

    def __init__(self, base_ring, category=None):
        from sage.categories.all import Sets, Rings

        if base_ring not in Rings():
            raise TypeError("base ring must be a ring")

        super().__init__(base=base_ring, category=category or Sets())

    def _an_element_(self):
        return self((1, 0))

    def some_element(self):
        return [self(v) for v in self.ambient_space().some_elements() if v]

    @cached_method
    def ambient_space(self):
        r"""
        Return the ambient Euclidean space containing these rays.

        EXAMPLES::

            sage: from flatsurf.geometry.ray import Rays
            sage: Rays(QQ).ambient_space()
            Vector space of dimension 2 over Rational Field

        """
        return self.base_ring() ** 2

    def _repr_(self):
        return f"Rays in {self.ambient_space()}"
