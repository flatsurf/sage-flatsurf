r"""
Thurston-Veech construction

A Thurston-Veech surface is a translation surface whose horizontal and vertical
directions are parabolic, that is stabilized by a multitwist. Combinatorially,
such surface decomposes into rectangles with horizontal and vertical sides.

We encode such surface by labelling the rectangles by `\{1, 2, \ldots, n\}`.
Then, the gluings is described by two permutations `h` and `v` that describe
respectively horizontal gluings (going right) and vertical gluings (going up).
The horizontal and vertical cylinders are then respectively the cycles
of `h` and `v`.

As shown by Thurston and Veech, to fix the flat structure of the surface one
only has to specify the topological data of the twists, that is their
multiplicities in each cylinder.

REFERENCES:

- Pascal Hubert, Erwan Lanneau
  "Veech groups without parabolic elements"
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2020 Vincent Delecroix
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

from sage.all import ZZ, QQ, AA, matrix, diagonal_matrix, NumberField
from sage.rings.qqbar import do_polred

from surface_dynamics.flat_surfaces.origamis.origami import Origami
from surface_dynamics.misc.permutation import perm_dense_cycles

from .polygon import ConvexPolygons
from .surface import Surface_list
from .translation_surface import TranslationSurface

class ThurstonVeech:
    def __init__(self, hp, vp):
        r"""
        INPUT:

        - ``hp`` - permutation describing the horizontal gluing of the rectangles

        - ``vp`` - permutation describing the vertical gluing of the rectangles

        EXAMPLES::

            sage: from flatsurf.geometry.thurston_veech import ThurstonVeech
            sage: TV = ThurstonVeech('(1,2)(3,4)', '(2,3)')
            sage: TV
            ThurstonVeech("(1,2)(3,4)", "(1)(2,3)(4)")
            sage: TV.stratum()
            H_2(1^2)

            sage: S = TV([1,2], [3,1,1])
            sage: S
            TranslationSurface built from 4 polygons
            sage: S.stratum()
            H_2(1^2)
            sage: S.base_ring()
            Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?

            sage: S = TV([1,1], [2,1,2])
            sage: S.base_ring()
            Rational Field
        """
        o = self._o = Origami(hp, vp)
        n = o.nb_squares()
        hcycles, hsizes = perm_dense_cycles(o.r_tuple(), n)
        vcycles, vsizes = perm_dense_cycles(o.u_tuple(), n)
        self._hcycles = hcycles
        self._vcycles = vcycles
        self._num_hcyls = len(hsizes)
        self._num_vcyls = len(vsizes)
        E = self._E = matrix(ZZ, self._num_hcyls, self._num_vcyls)
        for i in range(n):
            E[hcycles[i], vcycles[i]] += 1
        E.set_immutable()

    def __repr__(self):
        return "ThurstonVeech(\"{}\", \"{}\")".format(self._o.r().cycle_string(singletons=True),
                                              self._o.u().cycle_string(singletons=True))
    
    def stratum(self):
        return self._o.stratum()

    def stratum_component(self):
        return self._o.stratum_component()

    def cylinder_intersection_matrix(self):
        return self._E

    def __call__(self, hmult, vmult):
        r"""
        INPUT:

        - ``hmult`` -- multiplicities of the horizontal twists

        - ``vmult`` -- multiplicities of the vertical twists
        """
        if len(hmult) != self._num_hcyls or len(vmult) != self._num_vcyls:
            raise ValueError("invalid input lengths")

        E = self._E
        H = E * diagonal_matrix(vmult)
        V = E.transpose() * diagonal_matrix(hmult)

        if self._num_hcyls < self._num_vcyls:
            F = H * V
        else:
            F = V * H
        p = F.charpoly()
        assert F.nrows() == F.ncols() == min([self._num_hcyls, self._num_vcyls])

        pf = max(p.roots(AA, False))
        mp = pf.minpoly()
        if mp.degree() == 1:
            K = QQ
            pf = QQ(pf)
        else:
            fwd, bck, q = do_polred(pf.minpoly())
            im_gen = fwd(pf)
            K = NumberField(q, 'a', embedding=im_gen)
            pf = bck(K.gen())

        # Compute widths of the cylinders via Perron-Frobenius
        if self._num_hcyls < self._num_vcyls:
            hcirc = (F - pf).right_kernel_matrix()
            assert hcirc.nrows() == 1
            assert hcirc.ncols() == self._num_hcyls
            hcirc = hcirc[0]
            assert all(x > 0 for x in hcirc)
            vcirc = V * hcirc
            c = 1
            d = pf
        else:
            vcirc = (F - pf).right_kernel_matrix()
            assert vcirc.nrows() == 1
            assert vcirc.ncols() == self._num_vcyls
            vcirc = vcirc[0]
            assert all(x > 0 for x in vcirc)
            hcirc = H * vcirc
            d = 1
            c = pf

        # Solve linear systems to get heights
        h = [hcirc[i] * hmult[i] / c for i in range(self._num_hcyls)]
        v = [vcirc[i] * vmult[i] / d for i in range(self._num_vcyls)]

        C = ConvexPolygons(K)
        P = []
        for i in range(self._o.nb_squares()):
            hi = h[self._hcycles[i]]
            vi = v[self._vcycles[i]]
            P.append(C(edges=[(vi,0),(0,hi),(-vi,0),(0,-hi)]))

        surface = Surface_list(base_ring=K)
        for p in P:
            surface.add_polygon(p)
        r = self._o.r_tuple()
        u = self._o.u_tuple()
        for i in range(self._o.nb_squares()):
            surface.set_edge_pairing(i, 1, r[i], 3)
            surface.set_edge_pairing(i, 0, u[i], 2)
        surface.set_immutable()
        return TranslationSurface(surface)
