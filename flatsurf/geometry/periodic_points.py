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

from thurston_veech import ThurstonVeech

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
            for R1 in horizontal_cylinder.regions(start=R, multiplicity=horizontal_cylinder.multiplicity()):
                connecting_constraint = self._connecting_equation(R, R1)
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


    '''
    given a rectangle (H cap V), write down the constraints in that rectangle
    assuming that c(H)/h(V) is irrational
    we'd need:
    all the rectangles in the horizontal cylinders
    global horizontal multitwist
    multiplicities of the horizontal cylinder

    y / ht(H) in Q
    x / ht(V) in Q
    (x + ay - d) / ht(V') in Q,
    where d = sum ht(V_k) adding from V --> V_k
    and T = [1 a | 0 1]
    Note: the left-hand sides are elements of K[x, y]

    pi_j : (K = Q^n) --> Q
    pi_j (y / ht(H)) = 0 for all 1 <= j <= n
    pi_j (x / ht) ...
    pi_j (...)

    each constraint --> k1 x + k2 y + k3 in Q
    and now project k1, k2, k3 to get coefficient for each number field basis elt.

    let alpha_1, .. alpha_n be a Q-basis for K
    then
    x = x^i alpha_i
    k = a^i alpha_i
    (a^i alpha_i) * (x^j alpha_j) + (b^i alpha_i) * (y^j alpha_j) + (c^i alpha_i) in Q

    y / ht(H) in Q

    e.g. (y^1 + y^2 rt2) / (a^1 + a^2 rt2) in Q
    [(y1 + y2 rt2) * (a1 - a2rt2)] / (a1**2 - 2 a2**2)
    [a1 y1 - a2y1rt2 + y2 a1 rt2 - 2 a2y2] / (a1**2 - 2 a2**2) in Q
    (a1 y1 - 2 a2 y2)/(a1**2 - 2 a2**2) + rt2 (y2 a1 - a2 y1)/(a1**2 - 2 a2**2) in Q
    for each element of the basis that is not rational, get equations
    rt2 equation is ==> (a1 y2 - a2 y1)/(a1**2 - 2 a2**2) = 0
    i.e. [a1/(a1**2 - 2 a2**2)] y2 - [a2/(a1**2 - 2 a2**2)] y1 = 0

    now solve the resulting system of equations
    '''




