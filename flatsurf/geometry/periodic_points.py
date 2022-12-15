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
        pass

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

    


