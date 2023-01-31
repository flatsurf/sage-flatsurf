######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
#                      2022 Sam Freedman
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

from flatsurf.geometry.thurston_veech import ThurstonVeech

class PeriodicPoints():
    def __init__(self, T):
        # we will need hp, vp, hmult, vmult
        self._surface = T
        pass

    def list(self):
        r'''Return the periodic points'''

        G = self._candidate_graph()
        G = self._prune_candidate_graph(G)
        return G.vertices()

    def _horizontal_constraint(self, height):
        r'''Constrain y coordinate of point with rational height lemma

        EXAMPLE::

            sage: import flatsurf as fs
            sage: X = fs.translation_surfaces.mcmullen_genus2_prototype(1, 1, 0, -1)
            sage: P = PeriodicPoints(X)
            sage: R = X.base_ring()
            sage: P._horizontal_constraint(R(R.gen()))
            ([0 0 1 0], (0))
        '''

        A = (~height).matrix().submatrix(1)
        A = A.parent().zero().augment(A)
        return A, A.column_space().zero()

    def _vertical_constraint(self, width):
        r'''Constrain x coordinate of point with rational height lemma

        EXAMPLE::

            sage: import flatsurf as fs
            sage: X = fs.translation_surfaces.mcmullen_genus2_prototype(1, 1, 0, -1)
            sage: P = PeriodicPoints(X)
            sage: R = X.base_ring()
            sage: P._vertical_constraint(R(R.gen()))
            ([1 0 0 0], (0))

        '''
        A = (~width).matrix().submatrix(1)
        A = A.augment(A.parent().zero())
        return A, A.column_space().zero()

    def _connecting_constraint_horizontal(self, R, distance):
        return self._connecting_constraint_generic(self._surface.horizontal_twist()[0][1], distance, R.width())

    def _connecting_constraint_vertical(self, R, distance):
        return self._connecting_constraint_generic(self._surface.vertical_twist()[0][1], distance, R.height())

    def _connecting_constraint_generic(self, twist, distance, side_length):
        A_t = twist.matrix()
        A = (~side_length).matrix() * A_t.parent().one().augment(A_t)
        v = (~side_length).matrix() * distance.vector()
        return A.submatrix(1), v[1:]

    def _constraint_segments_in_region(self, R):
        r'''Return segments containing all periodic points in region ``R`` '''
        segments = []

        width, height = R.width(), R.height()
        horizontal_cylinder = R.horizontal_cylinder()
        vertical_cylinder = R.vertcial_cylinder()

        horizontal_constraint =  self._horizontal_constraint(height)
        vertical_constraint = self._vertical_constraint(width)

        from sage.all import QQ
        if (horizontal_cylinder.circumference() / width) not in QQ:
            for R1, wrap in horizontal_cylinder.regions(start=R, multiplicity=horizontal_cylinder.multiplicity()):
                connecting_constraint = self._connecting_constraint(R1, wrap)
                segments.append(self._solve_constraints(horizontal_constraint, vertical_constraint, connecting_constraint))

        else:
            assert (vertical_cylinder.circumference() / height) not in QQ
            for R1 in vertical_cylinder.regions(start=R, multiplicity=vertical_cylinder.multiplicity()):
                connecting_constraint = self._connecting_equation(R, R1)
                segments.append(self._solve_constraints(horizontal_constraint, vertical_constraint, connecting_constraint))

    def _constraint_segments(self):
        r'''Return set of segments containing all periodic points'''
        return {l
                for R in self._surface.regions()
                for l in self._constraint_segments_in_region(R)}

    def _candidate_points(self):
        r'''Return a (potentially large) set of candidate periodic points'''

        candidate_points = set()
        S = self._constraint_segments()

        while True:
            g = self._surface.veech_group.random_element()
            S1 = {l.apply_matrix(g) for l in S}
            if {l.slope() for l in S}.isdisjoint({l1.slope() for l1 in S1}):
                break

        for l in S:
            for l1 in S1:
                candidate_points.add(l.intersection(l1))

        return candidate_points

    def _candidate_graph(self):
        r'''Return graph with:
        vertices being either a set S of potential periodic points or ``None``, and
        p <--> q iff there is a generator of the Veech Group sending p to q,
        where p <--> ``None`` when a generator g sends p outside of S
        '''
        S = self._candidate_points()
        S.add(None)

        G = Graph()

        for p in S:
            for gen in self._surface.veech_group.gens():
                q = p.apply_matrix(gen)
                if q in S:
                    G.add_edge(p, q)
                else:
                    G.add_edge(p, None)

        return G

    def _prune_candidate_graph(self, G):
        r'''Given candidate graph, return the subgraph of periodic points'''

        # remove the connected component of ``None``
        raise NotImplementedError

