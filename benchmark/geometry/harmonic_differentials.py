# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
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
# ********************************************************************
import flatsurf


def time_harmonic_differential(surface):
    if surface == "TORUS":
        surface = flatsurf.translation_surfaces.torus((1, 0), (0, 1))
    elif surface == "3413":
        E = flatsurf.EquiangularPolygons(3, 4, 13)
        P = E.an_element()
        surface = flatsurf.similarity_surfaces.billiard(P, rational=True).minimal_cover(cover_type="translation")
    else:
        raise NotImplementedError

    surface = surface.delaunay_triangulation()
    surface.set_immutable()
    Ω = flatsurf.HarmonicDifferentials(surface)
    a = flatsurf.SimplicialHomology(surface).gens()[0]
    H = flatsurf.SimplicialCohomology(surface)
    Ω(H({a: 1}), check=False)


time_harmonic_differential.params = ([
    "TORUS",
    "3413",
])
