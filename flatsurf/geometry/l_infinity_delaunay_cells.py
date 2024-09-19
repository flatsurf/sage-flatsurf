r"""
Cells for the L-infinity Delaunay triangulation.

Each cell of the L-infinity Delaunay triangulation can be identified with a
marked triangulation. The marking corresponds to the position of the horizontal
and vertical separatrices. Each triangle hence get one of the following types:
bottom-left, bottom-right, top-left, top-right.
"""

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
from sage.misc.cachefunc import cached_method

# the types of edges
V_NONE = 0  # start vertex has no horizontal/vertical separatrix
V_LEFT = 1  # horizontal separatrix going right
V_RIGHT = 2  # horizontal separatrix going left
V_BOT = 3  # vertical separatrix going up
V_TOP = 4  # vertical separatrix going down


def sign_and_norm_conditions(dim, i, s):
    r"""
    Inequalities:

    if `s=+1`, encode the inequalities `+1 >= x_i >= 0`
    if `s=-1`, encode the inequalities `-1 <= x_i <= 0`

    EXAMPLES::

        sage: from flatsurf.geometry.l_infinity_delaunay_cells import sign_and_norm_conditions

        sage: sorted(Polyhedron(ieqs=sign_and_norm_conditions(1, 0, 1)).vertices_list())
        [[0], [1]]
        sage: sorted(Polyhedron(ieqs=sign_and_norm_conditions(1, 0, -1)).vertices_list())
        [[-1], [0]]

        sage: ieqs = []
        sage: ieqs.extend(sign_and_norm_conditions(2, 0, 1))
        sage: ieqs.extend(sign_and_norm_conditions(2, 1, -1))
        sage: sorted(Polyhedron(ieqs=ieqs).vertices_list())
        [[0, -1], [0, 0], [1, -1], [1, 0]]
    """
    l_sign = [0] * (dim + 1)
    l_sign[i + 1] = s

    l_norm = [0] * (dim + 1)
    l_norm[i + 1] = -s
    l_norm[0] = 1

    return (l_sign, l_norm)


def opposite_condition(dim, i, j):
    r"""
    encode the equality `x_i = -x_j`

    EXAMPLES::

        sage: from flatsurf.geometry.l_infinity_delaunay_cells import \
        ....:     opposite_condition, sign_and_norm_conditions

        sage: eq1 = opposite_condition(2, 0, 1)
        sage: eq2 = opposite_condition(2, 1, 0)

        sage: ieqs1 = sign_and_norm_conditions(2, 0, 1)
        sage: ieqs2 = sign_and_norm_conditions(2, 1, -1)

        sage: sorted(Polyhedron(eqns=[eq1], ieqs=ieqs1).vertices_list())
        [[0, 0], [1, -1]]
        sage: sorted(Polyhedron(eqns=[eq1,eq2], ieqs=ieqs1).vertices_list())
        [[0, 0], [1, -1]]
        sage: sorted(Polyhedron(eqns=[eq1,eq2], ieqs=ieqs1+ieqs2).vertices_list())
        [[0, 0], [1, -1]]
    """
    linear = [0] * (dim + 1)
    linear[i + 1] = 1
    linear[j + 1] = 1

    return linear


def bottom_top_delaunay_condition(dim, p1, e1, p2, e2):
    r"""
    Delaunay condition for bottom-top pairs of triangles
    """
    # re(e2p1) <= im(e1) + im(e2m1)
    e2m1 = (e2 + 2) % 3
    e2p1 = (e2 + 1) % 3

    im_e2m1 = 2 * (3 * p2 + e2m1) + 1
    re_e2p1 = 2 * (3 * p2 + e2p1)
    im_e1 = 2 * (3 * p1 + e1) + 1

    linear = [0] * (dim + 1)
    linear[im_e1 + 1] = 1
    linear[im_e2m1 + 1] = 1
    linear[re_e2p1 + 1] = -1

    return linear


def right_left_delaunay_condition(dim, p1, e1, p2, e2):
    r"""
    Delaunay condition for right-left pairs of triangles
    """
    # im(e2p1) <= re(e2) + re(e1m1)
    e1m1 = (e1 + 2) % 3
    e2p1 = (e2 + 1) % 3

    im_e2p1 = 3 * (p2 + e2p1) + 1
    re_e2 = 3 * (p2 + e2)
    re_e1m1 = 3 * (p1 + e1m1)

    linear = [0] * (dim + 1)
    linear[re_e2 + 1] = 1
    linear[re_e1m1 + 1] = 1
    linear[im_e2p1 + 1] = -1

    return linear


class LInfinityMarkedTriangulation:
    r"""
    EXAMPLES::

        sage: from flatsurf.geometry.l_infinity_delaunay_cells import \
        ....:     V_NONE, V_BOT, V_TOP, V_RIGHT, V_LEFT, LInfinityMarkedTriangulation

        sage: gluings = {(0,0):(1,2), (1,2):(0,0), (0,1):(1,0), (1,0):(0,1),
        ....:            (0,2):(1,1), (1,1):(0,2)}
        sage: types = [(V_BOT, V_NONE, V_RIGHT), (V_NONE, V_LEFT, V_TOP)]
        sage: T = LInfinityMarkedTriangulation(2, gluings, types)
    """

    def __init__(self, num_faces, edge_identifications, edge_types, check=True):
        from sage.rings.integer_ring import ZZ

        self._n = ZZ(num_faces)
        self._edge_identifications = edge_identifications
        self._edge_types = edge_types
        if check:
            self._check()

    def _check(self):
        if self._n % 2:
            raise ValueError("the number of faces must be even")

        if sorted(self._edge_identifications.keys()) != [
            (i, j) for i in range(self._n) for j in range(3)
        ]:
            raise ValueError("should be a triangulation")

        if (
            not isinstance(self._edge_types, list)
            or len(self._edge_types) != self.num_faces()
        ):
            raise ValueError("edge_types invalid")

        for i in range(self.num_faces()):
            if len(self._edge_types[i]) != 3:
                raise ValueError("edge_types invalid")
            for j in range(3):
                if self._edge_types[i][j] not in [
                    V_NONE,
                    V_LEFT,
                    V_RIGHT,
                    V_BOT,
                    V_TOP,
                ]:
                    raise ValueError("types[{}] = {} invalid", i, self._edge_types[i])

        seen = [False] * self._n
        for p in range(self._n):
            if seen[p]:
                continue
            sh = sum(
                self._edge_types[p][r] == V_LEFT or self._edge_types[p][r] == V_RIGHT
                for r in (0, 1, 2)
            )
            sv = sum(
                self._edge_types[p][r] == V_BOT or self._edge_types[p][r] == V_TOP
                for r in (0, 1, 2)
            )
            if sh != 1 or sv != 1:
                raise ValueError(
                    "triangle {} has invalid types {}".format(p, self._edge_types[p])
                )

    def num_faces(self):
        return self._n

    def num_edges(self):
        return 3 * self._n // 2

    def opposite_edge(self, p, e):
        return self._edge_identifications[(p, e)]

    def __repr__(self):
        return "Marked triangulation made of {} triangles".format(self.num_faces())

    def bottom_top_pairs(self):
        r"""
        Return a list ``(p1,e1,p2,e2)``.

        EXAMPLES::

            sage: from flatsurf.geometry.l_infinity_delaunay_cells import \
            ....:     V_NONE, V_BOT, V_TOP, V_RIGHT, V_LEFT, LInfinityMarkedTriangulation

            sage: gluings = {(0,0):(1,2), (1,2):(0,0), (0,1):(1,0), (1,0):(0,1),
            ....:            (0,2):(1,1), (1,1):(0,2)}
            sage: types = [(V_BOT, V_NONE, V_RIGHT), (V_NONE, V_LEFT, V_TOP)]
            sage: T = LInfinityMarkedTriangulation(2, gluings, types)
            sage: T.bottom_top_pairs()
            [(0, 0, 1, 2)]
        """
        pairs = []
        for p1 in range(self._n):
            for e1 in range(3):
                if self._edge_types[p1][e1] == V_BOT:
                    e1p1 = (e1 + 1) % 3
                    p2, e2 = self.opposite_edge(p1, e1p1)
                    e2m1 = (e2 - 1) % 3
                    if self._edge_types[p2][e2m1] == V_TOP:
                        pairs.append((p1, e1, p2, e2m1))
        return pairs

    def right_left_pairs(self):
        r"""
        Return a list ``(p1,e1,p2,e2)``

        EXAMPLES::

            sage: from flatsurf.geometry.l_infinity_delaunay_cells import \
            ....:     V_NONE, V_BOT, V_TOP, V_RIGHT, V_LEFT, LInfinityMarkedTriangulation

            sage: gluings = {(0,0):(1,2), (1,2):(0,0), (0,1):(1,0), (1,0):(0,1),
            ....:            (0,2):(1,1), (1,1):(0,2)}
            sage: types = [(V_BOT, V_NONE, V_LEFT), (V_NONE, V_RIGHT, V_TOP)]
            sage: T = LInfinityMarkedTriangulation(2, gluings, types)
            sage: T.right_left_pairs()
            [(1, 1, 0, 2)]
        """
        pairs = []
        for p1 in range(self._n):
            for e1 in range(3):
                if self._edge_types[p1][e1] == V_RIGHT:
                    e1p1 = (e1 + 1) % 3
                    p2, e2 = self.opposite_edge(p1, e1p1)
                    e2m1 = (e2 - 1) % 3
                    if self._edge_types[p2][e2m1] == V_LEFT:
                        pairs.append((p1, e1, p2, e2m1))
        return pairs

    @cached_method
    def polytope(self):
        r"""
        Each edge correspond to a vector in RR^2 (identified to CC)

        We assign the following coordinates

        (p,e) -> real part at 2*(3*p + e) and imag part at 2*(3*p + e) + 1

        The return polyhedron is compact as we fix each side to be of L-infinity
        length less than 1.
        """
        dim = 4 * self.num_edges()
        eqns = []
        ieqs = []

        signs = [None] * dim

        # edges should sum up to zero
        for p in range(self._n):
            linear = [0] * (dim + 1)
            linear[6 * p + 1] = 1
            linear[6 * p + 3] = 1
            linear[6 * p + 5] = 1
            eqns.append(linear)

            linear = [0] * (dim + 1)
            linear[6 * p + 2] = 1
            linear[6 * p + 4] = 1
            linear[6 * p + 6] = 1
            eqns.append(linear)

        # opposite edges are opposite vectors
        for p1 in range(self._n):
            for e1 in range(3):
                p2, e2 = self.opposite_edge(p1, e1)
                re1 = 2 * (3 * p1 + e1)
                im1 = 2 * (3 * p1 + e1) + 1
                re2 = 2 * (3 * p2 + e2)
                im2 = 2 * (3 * p2 + e2) + 1
                if re1 < re2:
                    eqns.append(opposite_condition(dim, re1, re2))
                    eqns.append(opposite_condition(dim, im1, im2))

        # Compute the signs depending on edge types
        for p in range(self._n):
            for e1, e2 in ((2, 0), (0, 1), (1, 2)):
                t = self._edge_types[p][e2]
                re1 = 2 * (3 * p + e1)
                im1 = 2 * (3 * p + e1) + 1
                re2 = 2 * (3 * p + e2)
                im2 = 2 * (3 * p + e2) + 1
                if t == V_BOT:
                    signs[re1] = signs[re2] = +1
                    signs[im1] = -1
                    signs[im2] = +1
                elif t == V_TOP:
                    signs[re1] = signs[re2] = -1
                    signs[im1] = +1
                    signs[im2] = -1
                elif t == V_RIGHT:
                    signs[re1] = +1
                    signs[re2] = -1
                    signs[im1] = signs[im2] = +1
                elif t == V_LEFT:
                    signs[re1] = -1
                    signs[re2] = +1
                    signs[im1] = signs[im2] = -1

        # adding sign conditions
        for i in range(dim):
            ieqs.extend(sign_and_norm_conditions(dim, i, signs[i]))

        # Delaunay conditions
        for p1, e1, p2, e2 in self.bottom_top_pairs():
            ieqs.append(bottom_top_delaunay_condition(dim, p1, e1, p2, e2))
        for p1, e1, p2, e2 in self.right_left_pairs():
            ieqs.append(right_left_delaunay_condition(dim, p1, e1, p2, e2))

        #        return eqns, ieqs

        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.rings.rational_field import QQ

        return Polyhedron(ieqs=ieqs, eqns=eqns, base_ring=QQ)

    def barycenter(self):
        r"""
        Return the translation surface in the barycenter of this polytope.

        EXAMPLES::

            sage: from flatsurf.geometry.l_infinity_delaunay_cells import \
            ....:     V_NONE, V_BOT, V_TOP, V_RIGHT, V_LEFT, LInfinityMarkedTriangulation

            sage: gluings = {(0,0):(1,2), (1,2):(0,0), (0,1):(1,0), (1,0):(0,1),
            ....:            (0,2):(1,1), (1,1):(0,2)}
            sage: types = [(V_BOT, V_NONE, V_LEFT), (V_NONE, V_RIGHT, V_TOP)]
            sage: T = LInfinityMarkedTriangulation(2, gluings, types)
            sage: S = T.barycenter()
            sage: S.polygon(0)
            Polygon(vertices=[(0, 0), (3/7, 13/21), (-3/7, 11/21)])
            sage: S.polygon(1)
            Polygon(vertices=[(0, 0), (6/7, 2/21), (3/7, 13/21)])
        """
        verts = [v.vector() for v in self.polytope().vertices()]
        b = sum(verts) / len(verts)

        from flatsurf import Polygon
        from sage.rings.rational_field import QQ

        from flatsurf import MutableOrientedSimilaritySurface

        barycenter = MutableOrientedSimilaritySurface(QQ)

        for p in range(self._n):
            e1 = (b[6 * p], b[6 * p + 1])
            e2 = (b[6 * p + 2], b[6 * p + 3])
            e3 = (b[6 * p + 4], b[6 * p + 5])
            barycenter.add_polygon(Polygon(edges=[e1, e2, e3], base_ring=QQ))

        for gluing in self._edge_identifications.items():
            barycenter.glue(*gluing)

        return barycenter
