r"""
Wrapper for flow decompositions of surfaces that are backed by libflatsurf.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.regular_octagon()

    sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
    sage: decomposition = S.pyflatsurf().codomain().flow_decomposition(direction)  # optional: pyflatsurf

    sage: decomposition  # optional: pyflatsurf
    Flow decomposition of Surface backed by FlatTriangulationCombinatorial(vertices = (1, -3, -8, 4, -5, -6, 9, 8, -7, 6, 2, -1, -4, -9, 7, 3, -2, 5), faces = (1, 2, 3)(-1, 5, 4)(-2, 6, -5)(-3, 7, 8)(-4, -8, 9)(-6, -7, -9)) with vectors {1: (1, 0), 2: ((1/2*a ~ 0.70710678), (1/2*a ~ 0.70710678)), 3: ((-1/2*a-1 ~ -1.7071068), (-1/2*a ~ -0.70710678)), 4: (1, (a+1 ~ 2.4142136)), 5: (0, (-a-1 ~ -2.4142136)), 6: ((1/2*a ~ 0.70710678), (-1/2*a-1 ~ -1.7071068)), 7: (0, 1), 8: ((-1/2*a-1 ~ -1.7071068), (-1/2*a-1 ~ -1.7071068)), 9: ((-1/2*a ~ -0.70710678), (1/2*a ~ 0.70710678))} into 1 undetermined component

"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2026 Julian Rüth
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
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.flow_decomposition import FlowDecomposition_base, FlowComponent_base


class FlowDecomposition_pyflatsurf(FlowDecomposition_base):
    r"""
    Flow decomposition backed by a ``FlowDecomposition`` in libflatsurf.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()

        sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
        sage: decomposition = S.pyflatsurf().codomain().flow_decomposition(direction)  # optional: pyflatsurf

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.flow_decomposition import FlowDecomposition_pyflatsurf
        sage: isinstance(decomposition, FlowDecomposition_pyflatsurf)
        True

    """
    def __init__(self, surface, direction):
        vector_space_conversion = surface.vector_space_conversion()

        if type(direction) is not vector_space_conversion.codomain():
            direction = vector_space_conversion(vector_space_conversion.domain()(direction))

        self._surface = surface
        self._direction = direction

    @cached_method
    def _decomposition(self):
        from flatsurf.features import pyflatsurf_feature

        pyflatsurf_feature.require()
        import pyflatsurf

        return pyflatsurf.flatsurf.makeFlowDecomposition(
            self._surface.flat_triangulation(), self._direction
        )

    def decompose(self, limit=None):
        if limit == 0:
            return

        if limit is None:
            limit = -1

        self._decomposition().decompose(int(limit))

    def components(self):
        # TODO: Track components so we can trash them when they are modified.
        return [FlowComponent_pyflatsurf(self, component) for component in self._decomposition().components()]


class FlowComponent_pyflatsurf(FlowComponent_base):
    def __init__(self, parent: FlowDecomposition_pyflatsurf, component):
        super().__init__(parent)

        self._component = component

    @staticmethod
    def _from_tribool(tribool):
        if tribool == True:
            return True
        if tribool == False:
            return False

        from sage.misc.unknown import Unknown
        return Unknown

    def is_cylinder(self):
        return self._from_tribool(self._component.cylinder())

    def is_minimal(self):
        return self._from_tribool(self._component.withoutPeriodicTrajectory())

    def is_undetermined(self):
        from sage.misc.unknown import Unknown
        return self._from_tribool(self._component.cylinder()) == Unknown
