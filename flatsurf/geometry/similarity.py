# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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
# *********************************************************************
from sage.misc.cachefunc import cached_method

from sage.structure.element import MultiplicativeGroupElement, parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.groups import Groups

from sage.all import Rings

from sage.modules.free_module_element import vector
from sage.groups.group import Group
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import FreeModuleElement

from sage.structure.element import Matrix

ZZ_0 = Integer(0)
ZZ_1 = Integer(1)
ZZ_m1 = -ZZ_1


class Similarity(MultiplicativeGroupElement):
    r"""
    Class for a similarity of the plane with possible reflection.

    Construct the similarity (x,y) mapsto (ax-by+s,bx+ay+t) if sign=1,
    and (ax+by+s,bx-ay+t) if sign=-1
    """

    def __init__(self, p, a, b, s, t, sign):
        r"""
        Construct the similarity (x,y) mapsto (ax-by+s,bx+ay+t) if sign=1,
        and (ax+by+s,bx-ay+t) if sign=-1
        """
        if p is None:
            raise ValueError("The parent must be provided")

        if parent(a) is not p.base_ring():
            raise ValueError("wrong parent for a")
        if parent(b) is not p.base_ring():
            raise ValueError("wrong parent for b")
        if parent(s) is not p.base_ring():
            raise ValueError("wrong parent for s")
        if parent(t) is not p.base_ring():
            raise ValueError("wrong parent for t")
        if parent(sign) is not ZZ or not sign.is_unit():
            raise ValueError("sign must be either 1 or -1.")

        self._a = a
        self._b = b
        self._s = s
        self._t = t
        self._sign = sign

        MultiplicativeGroupElement.__init__(self, p)

    def sign(self):
        return self._sign

    def is_translation(self):
        r"""
        Return whether this element is a translation.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)).is_translation()
            True
            sage: S((1,0,3,-1/2)).is_translation()
            True
            sage: S((0,1,0,0)).is_translation()
            False
        """
        return self._sign.is_one() and self._a.is_one() and self._b.is_zero()

    def is_half_translation(self):
        r"""
        Return whether this element is a half translation.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)).is_half_translation()
            True
            sage: S((-1, 0, 0, 2)).is_half_translation()
            True
            sage: S((0,1,0,0)).is_half_translation()
            False
        """
        return (
            self._sign.is_one()
            and (self._a.is_one() or ((-self._a).is_one()))
            and self._b.is_zero()
        )

    def is_orientable(self):
        return self._sign.is_one()

    def is_rotation(self):
        r"""
        Check whether this element is a rotation

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)).is_rotation()
            False
            sage: S((0,-1,0,0)).is_rotation()
            True
            sage: S.one().is_rotation()
            True
        """
        return self.is_one() or (self.det().is_one() and not self.is_translation())

    def is_isometry(self):
        r"""
        Check whether this element is an isometry

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S.one().is_isometry()
            True
            sage: S((0,1,0,0)).is_isometry()
            True
            sage: S((0,1,0,0,-1)).is_isometry()
            True
            sage: S((1,1,0,0)).is_isometry()
            False
            sage: S((3,-1/2)).is_isometry()
            True
        """
        det = self.det()
        return det.is_one() or (-det).is_one()

    def det(self):
        r"""
        Return the determinant of this element
        """
        return self._sign * (self._a * self._a + self._b * self._b)

    def _mul_(self, right):
        r"""
        Composition

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)) * S((3,-5)) == S((4,-3))
            True

            sage: from itertools import product
            sage: a1 = S((0,2,0,0,1))
            sage: a2 = S((1,0,0,0,-1))
            sage: a3 = S((1,1,0,0))
            sage: a4 = S((1,0,-1,1))
            sage: a5 = S((2,-1,3/5,2/3,-1))
            sage: for g1,g2,g3 in product([a1,a2,a3,a4,a5], repeat=3):
            ....:     assert g1.matrix()*g2.matrix() == (g1*g2).matrix()
            ....:     assert (g1*g2).matrix()*g3.matrix() == (g1*g2*g3).matrix()
        """
        a = self._a * right._a - self._sign * self._b * right._b
        b = self._b * right._a + self._sign * self._a * right._b
        s = self._a * right._s - self._sign * self._b * right._t + self._s
        t = self._b * right._s + self._sign * self._a * right._t + self._t
        sign = self._sign * right._sign
        P = self.parent()
        return P.element_class(P, a, b, s, t, sign)

    def __invert__(self):
        r"""
        Invert a similarity.

        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: from itertools import product
            sage: for a in [S((0,2,0,0,1)), S((1,0,0,0,-1)), S((1,1,0,0)),
            ....:           S((1,0,-1,1)), S((2,-1,3/5,2/3,-1))]:
            ....:     assert (a*~a).is_one() and (~a*a).is_one()
        """
        P = self.parent()
        sign = self._sign
        det = self.det()
        a = sign * self._a / det
        b = -self._b / det
        return P.element_class(
            P,
            a,
            b,
            -a * self._s + sign * b * self._t,
            -b * self._s - sign * a * self._t,
            sign,
        )

    def _div_(self, right):
        det = right.det()

        inv_a = right._sign * right._a
        inv_b = -right._b
        inv_s = -right._sign * right._a * right._s - right._sign * right._b * right._t
        inv_t = right._b * right._s - right._a * right._t

        a = (self._a * inv_a - self._sign * self._b * inv_b) / det
        b = (self._b * inv_a + self._sign * self._a * inv_b) / det
        s = (self._a * inv_s - self._sign * self._b * inv_t) / det + self._s
        t = (self._b * inv_s + self._sign * self._a * inv_t) / det + self._t

        return self.parent().element_class(
            self.parent(),
            self.base_ring()(a),
            self.base_ring()(b),
            self.base_ring()(s),
            self.base_ring()(t),
            self._sign * right._sign,
        )

    def __hash__(self):
        return (
            73 * hash(self._a)
            - 19 * hash(self._b)
            + 13 * hash(self._s)
            + 53 * hash(self._t)
            + 67 * hash(self._sign)
        )

    def __call__(self, w, ring=None):
        r"""
        Return the image of ``w`` under the similarity. Here ``w`` may be a
        convex polygon or a vector (or something that can be indexed in the
        same way as a vector). If a ring is provided, the objects returned will
        be defined over this ring.

        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(AA)
            sage: a = S((1,-1,AA(2).sqrt(),0))
            sage: a((1,2))
            (4.414213562373095?, 1)
            sage: a.matrix()*vector((1,2,1))
            (4.414213562373095?, 1, 1)

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: SG = SimilarityGroup(QQ)
            sage: from flatsurf import Polygon
            sage: p = Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
            sage: g = SG.an_element()**2
            sage: g
            (x, y) |-> (25*x + 4, 25*y + 10)
            sage: g(p)
            Polygon(vertices=[(4, 10), (29, 10), (29, 35), (4, 35)])
            sage: g(p, ring=AA).category()
            Category of convex simple euclidean polygons over Algebraic Real Field

        """
        if ring is not None and ring not in Rings():
            raise TypeError("ring must be a ring")

        from flatsurf.geometry.polygon import EuclideanPolygon

        if isinstance(w, EuclideanPolygon) and w.is_convex():
            if ring is None:
                ring = self.parent().base_ring()

            from flatsurf import Polygon

            try:
                return Polygon(vertices=[self(v) for v in w.vertices()], base_ring=ring)
            except ValueError:
                if not self._sign.is_one():
                    raise ValueError("Similarity must be orientation preserving.")

                # Not sure why this would happen:
                raise

        if ring is None:
            if self._sign.is_one():
                return vector(
                    [
                        self._a * w[0] - self._b * w[1] + self._s,
                        self._b * w[0] + self._a * w[1] + self._t,
                    ]
                )
            else:
                return vector(
                    [
                        self._a * w[0] + self._b * w[1] + self._s,
                        self._b * w[0] - self._a * w[1] + self._t,
                    ]
                )
        else:
            if self._sign.is_one():
                return vector(
                    ring,
                    [
                        self._a * w[0] - self._b * w[1] + self._s,
                        self._b * w[0] + self._a * w[1] + self._t,
                    ],
                )
            else:
                return vector(
                    ring,
                    [
                        self._a * w[0] + self._b * w[1] + self._s,
                        self._b * w[0] - self._a * w[1] + self._t,
                    ],
                )

    def _repr_(self):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S.one()
            (x, y) |-> (x, y)
            sage: S((1,-2/3))
            (x, y) |-> (x + 1, y - 2/3)
            sage: S((-1,0,2/3,3))
            (x, y) |-> (-x + 2/3, -y + 3)
            sage: S((-1,0,2/3,3,-1))
            (x, y) |-> (-x + 2/3, y + 3)
        """
        R = self.parent().base_ring()["x", "y"]
        x, y = R.gens()
        return "(x, y) |-> ({}, {})".format(
            self._a * x - self._sign * self._b * y + self._s,
            self._b * x + self._sign * self._a * y + self._t,
        )

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,0)) == S((1,0))
            True
            sage: S((1,0)) == S((0,1))
            False
            sage: S((1,0,0,0)) == S((0,1,0,0))
            False
            sage: S((1,0,0,0,1)) == S((1,0,0,0,-1))
            False
        """
        if other is None:
            return False
        if isinstance(other, int):
            return False
        if self.parent() != other.parent():
            return False
        return (
            self._a == other._a
            and self._b == other._b
            and self._s == other._s
            and self._t == other._t
            and self._sign == other._sign
        )

    def __ne__(self, other):
        return not (self == other)

    def matrix(self):
        r"""
        Return the 3x3 matrix representative of this element

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)

            sage: S((1,-2/3,1,1,-1)).matrix()
            [   1 -2/3    1]
            [-2/3   -1    1]
            [   0    0    1]
        """
        P = self.parent()
        M = P._matrix_space_3x3()
        z = P._ring.zero()
        o = P._ring.one()
        return M(
            [
                self._a,
                -self._sign * self._b,
                self._s,
                self._b,
                +self._sign * self._a,
                self._t,
                z,
                z,
                o,
            ]
        )

    def derivative(self):
        r"""
        Return the 2x2 matrix corresponding to the derivative of this element

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)

            sage: S((1,-2/3,1,1,-1)).derivative()
            [   1 -2/3]
            [-2/3   -1]
        """
        M = self.parent()._matrix_space_2x2()
        return M([self._a, -self._sign * self._b, self._b, self._sign * self._a])


class SimilarityGroup(UniqueRepresentation, Group):
    r"""
    The group of possibly orientation reversing similarities in the plane.

    This is the group generated by rotations, translations and dilations.
    """

    Element = Similarity

    def __init__(self, base_ring):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: TestSuite(SimilarityGroup(QQ)).run()
            sage: TestSuite(SimilarityGroup(AA)).run()
        """
        self._ring = base_ring
        Group.__init__(self, category=Groups().Infinite())

    @cached_method
    def _matrix_space_2x2(self):
        from sage.matrix.matrix_space import MatrixSpace

        return MatrixSpace(self._ring, 2)

    @cached_method
    def _matrix_space_3x3(self):
        from sage.matrix.matrix_space import MatrixSpace

        return MatrixSpace(self._ring, 3)

    @cached_method
    def _vector_space(self):
        from sage.modules.free_module import VectorSpace

        return VectorSpace(self._ring, 2)

    def _element_constructor_(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,1))  # translation
            (x, y) |-> (x + 1, y + 1)

            sage: V = QQ^2
            sage: S(V((1,-1)))
            (x, y) |-> (x + 1, y - 1)

            sage: S(vector((1,1)))
            (x, y) |-> (x + 1, y + 1)
        """
        if len(args) == 1:
            x = args[0]
        else:
            x = args

        a = self._ring.one()
        b = s = t = self._ring.zero()
        sign = ZZ_1

        # TODO: 2x2 and 3x3 matrix input

        if isinstance(x, (tuple, list)):
            if len(x) == 2:
                s, t = map(self._ring, x)
            elif len(x) == 4:
                a, b, s, t = map(self._ring, x)
            elif len(x) == 5:
                a, b, s, t = map(self._ring, x[:4])
                sign = ZZ(x[4])
            else:
                raise ValueError(
                    "can not construct a similarity from a list of length {}".format(
                        len(x)
                    )
                )
        elif isinstance(x, Matrix):
            #   a -sb
            #   b sa
            if x.nrows() == x.ncols() == 2:
                a, c, b, d = x.list()
                if a == d and b == -c:
                    sign = ZZ_1
                elif a == -d and b == c:
                    sign = ZZ_m1
                else:
                    raise ValueError("not a similarity matrix")
            elif x.nrows() == x.ncols() == 3:
                raise NotImplementedError
            else:
                raise ValueError("invalid dimension for matrix input")
        elif isinstance(x, FreeModuleElement):
            if len(x) == 2:
                if x.base_ring() is self._ring:
                    s, t = x
                else:
                    s, t = map(self._ring, x)
            else:
                raise ValueError("invalid dimension for vector input")
        else:
            p = parent(x)
            if self._ring.has_coerce_map_from(p):
                a = self._ring(x)
            else:
                raise ValueError(
                    "element in {} cannot be used to create element in {}".format(
                        p, self
                    )
                )

        if (a * a + b * b).is_zero():
            raise ValueError("not invertible")

        return self.element_class(self, a, b, s, t, sign)

    def _coerce_map_from_(self, S):
        if self._ring.has_coerce_map_from(S):
            return True
        if isinstance(S, SimilarityGroup):
            return self._ring.has_coerce_map_from(S._ring)

    def _repr_(self):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: SimilarityGroup(QQ)
            Similarity group over Rational Field
        """
        return "Similarity group over {}".format(self._ring)

    def one(self):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: SimilarityGroup(QQ).one()
            (x, y) |-> (x, y)
            sage: SimilarityGroup(QQ).one().is_one()
            True
        """
        return self.element_class(
            self,
            self._ring.one(),  # a
            self._ring.zero(),  # b
            self._ring.zero(),  # s
            self._ring.zero(),  # t
            ZZ_1,
        )  # sign

    def _an_element_(self):
        r"""
        Return a typical element of this group.

        EXAMPLES:

            sage: from flatsurf.geometry.similarity import SimilarityGroup

            sage: SimilarityGroup(QQ)._an_element_()
            (x, y) |-> (3*x + 4*y + 2, 4*x - 3*y - 1)

            sage: SimilarityGroup(QQ).an_element()
            (x, y) |-> (3*x + 4*y + 2, 4*x - 3*y - 1)

        .. SEEALSO::

            :meth:`sage.structure.parent.Parent.an_element` which relies on
            this method and should be called instead

        """
        return self(3, 4, 2, -1, -1)

    def is_abelian(self):
        return False

    def base_ring(self):
        return self._ring


def similarity_from_vectors(u, v, matrix_space=None):
    r"""
    Return the unique similarity matrix that maps ``u`` to ``v``.

    EXAMPLES::

        sage: from flatsurf.geometry.similarity import similarity_from_vectors

        sage: V = VectorSpace(QQ,2)
        sage: u = V((1,0))
        sage: v = V((0,1))
        sage: m = similarity_from_vectors(u,v); m
        [ 0 -1]
        [ 1  0]
        sage: m*u == v
        True

        sage: u = V((2,1))
        sage: v = V((1,-2))
        sage: m = similarity_from_vectors(u,v); m
        [ 0  1]
        [-1  0]
        sage: m * u == v
        True

    An example built from the Pythagorean triple 3^2 + 4^2 = 5^2::

        sage: u2 = V((5,0))
        sage: v2 = V((3,4))
        sage: m = similarity_from_vectors(u2,v2); m
        [ 3/5 -4/5]
        [ 4/5  3/5]
        sage: m * u2 == v2
        True

    Some test over number fields::

        sage: K.<sqrt2> = NumberField(x^2-2, embedding=1.4142)
        sage: V = VectorSpace(K,2)
        sage: u = V((sqrt2,0))
        sage: v = V((1, 1))
        sage: m = similarity_from_vectors(u,v); m
        [ 1/2*sqrt2 -1/2*sqrt2]
        [ 1/2*sqrt2  1/2*sqrt2]
        sage: m*u == v
        True

        sage: m = similarity_from_vectors(u, 2*v); m
        [ sqrt2 -sqrt2]
        [ sqrt2  sqrt2]
        sage: m*u == 2*v
        True

    """
    if u.parent() is not v.parent():
        raise ValueError

    if matrix_space is None:
        from sage.matrix.matrix_space import MatrixSpace

        matrix_space = MatrixSpace(u.base_ring(), 2)

    if u == v:
        return matrix_space.one()

    sqnorm_u = u[0] * u[0] + u[1] * u[1]
    cos_uv = (u[0] * v[0] + u[1] * v[1]) / sqnorm_u
    sin_uv = (u[0] * v[1] - u[1] * v[0]) / sqnorm_u

    m = matrix_space([cos_uv, -sin_uv, sin_uv, cos_uv])
    m.set_immutable()
    return m
