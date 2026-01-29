r"""
Flow decompositions of surfaces into cylinders and minimal components.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.regular_octagon()

    sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
    sage: decomposition = S.flow_decomposition(direction)  # optional: pyflatsurf

    sage: decomposition  # optional: pyflatsurf
    Flow decomposition of Translation Surface in H_2(2) built from a regular octagon into 1 undetermined component

    sage: decomposition.decompose()  # optional: pyflatsurf
    sage: decomposition  # optional: pyflatsurf
    Flow decomposition of Translation Surface in H_2(2) built from a regular octagon into 2 cylinders

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

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.unknown import UnknownClass


class FlowDecomposition_base(SageObject):
    r"""
    Base class for flow decompositions of a surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()

        sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
        sage: decomposition = S.flow_decomposition(direction)  # optional: pyflatsurf

    TESTS::

        sage: from flatsurf.geometry.flow_decomposition import FlowDecomposition_base
        sage: isinstance(decomposition, FlowDecomposition_base)  # optional: pyflatsurf
        True

    """
    def __init__(self, surface, direction):
        from flatsurf.geometry.categories.translation_surfaces import TranslationSurfaces
        if surface not in TranslationSurfaces():
            raise NotImplementedError("surface must be a translation surface")

        if not direction or direction not in surface.base_ring() ** 2:
            raise TypeError("direction must a direction in the surface")


        self._surface = surface
        self._direction = direction

    def surface(self):
        return self._surface

    def direction(self):
        return self._direction

    def _repr_(self):
        r"""
        Return a printable representation of this decomposition.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
            sage: S.flow_decomposition(direction)  # optional: pyflatsurf
            Flow decomposition of Translation Surface in H_2(2) built from a regular octagon into 1 undetermined component
            
        """
        components = []

        cylinders = len(self.cylinders())
        if cylinders:
            components.append(f"{cylinders} cylinder{'s' if cylinders > 1 else ''}")

        minimals = len(self.minimal_components())
        if minimals:
            components.append(f"{minimals} mininmal component{'s' if minimals > 1 else ''}")

        undetermineds = len(self.undetermined_components())
        if undetermineds:
            components.append(f"{undetermineds} undetermined component{'s' if undetermineds > 1 else ''}")

        return f"Flow decomposition of {self._surface} into {'and'.join(components)}"

    def cylinders(self):
        r"""
        Return the cylinders components in this decomposition.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
            sage: decomposition = S.flow_decomposition(direction)  # optional: pyflatsurf
            sage: decomposition.decompose()

            sage: decomposition.cylinders()
            [Cylinder component, Cylinder component]

        """
        return [c for c in self.components() if c.is_cylinder() == True]

    def minimal_components(self):
        r"""
        Return the minimal components in this decomposition.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
            sage: decomposition = S.flow_decomposition(direction)  # optional: pyflatsurf
            sage: decomposition.decompose()

            sage: decomposition.minimal_components()
            []

        """
        return [c for c in self.components() if c.is_minimal() == True]

    def undetermined_components(self):
        r"""
        Return the undetermined components in this decomposition.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
            sage: decomposition = S.flow_decomposition(direction)  # optional: pyflatsurf

            sage: decomposition.undetermined_components()
            [Undetermined component]

            sage: decomposition.decompose()
            sage: decomposition.undetermined_components()
            []

        """
        return [c for c in self.components() if c.is_undetermined()]

    def components(self) -> list["FlowComponent_base"]:
        raise NotImplementedError

    def decompose(self, limit=None):
        raise NotImplementedError


class FlowDecomposition(FlowDecomposition_base):
    r"""
    A flow decomposition on a native sage-flatsurf surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()

        sage: direction = next(iter(S.slopes()))  # optional: pyflatsurf
        sage: decomposition = S.flow_decomposition(direction)  # optional: pyflatsurf

    TESTS::

        sage: from flatsurf.geometry.flow_decomposition import FlowDecomposition
        sage: isinstance(decomposition, FlowDecomposition)  # optional: pyflatsurf
        True

    """
    def __init__(self, surface, direction):
        super().__init__(surface, direction)


    @cached_method
    def _backend(self):
        to_backend = self.surface().pyflatsurf()
        return to_backend, to_backend.codomain().flow_decomposition(self._direction)

    def components(self) -> list["FlowComponent_base"]:
        to_backend, backend = self._backend()
        return [FlowComponent_mapped(self, component, to_backend.section()) for component in backend.components()]

    def decompose(self, limit=None):
        _, backend = self._backend()
        backend.decompose(limit=limit)


class FlowComponent_base(SageObject):
    def __init__(self, parent: FlowDecomposition_base):
        pass

    def is_cylinder(self) -> bool | UnknownClass:
        raise NotImplementedError

    def is_minimal(self) -> bool | UnknownClass:
        raise NotImplementedError

    def is_undetermined(self) -> bool:
        raise NotImplementedError

    def area(self):
        raise NotImplementedError

    def _repr_(self):
        if self.is_cylinder() == True:
            return "Cylinder component"
        if self.is_minimal() == True:
            return "Minimal component"
        return "Undetermined component"

    ## TODO: Provide this somehow. E.g., by returning the cylinder circumference as a chain.
    ## def cylinder_circumference(self, component, A, sc_index, proj):
    ##     r"""
    ##     Return the circumference of the cylinder ``component`` in the homology
    ##     of the underlying surface.

    ##     INPUT:

    ##     - ``component`` -- a cylinder

    ##     - ``A``, ``sc_index``, ``proj`` -- the output of
    ##       ``flow_decomposition_kontsevich_zorich_cocycle``

    ##     EXAMPLES::

    ##         sage: from flatsurf import translation_surfaces
    ##         sage: from flatsurf import GL2ROrbitClosure  # optional: pyflatsurf

    ##         sage: S = translation_surfaces.veech_double_n_gon(5)
    ##         sage: O = GL2ROrbitClosure(S)  # optional: pyflatsurf
    ##         sage: dec = next(iter(S._decomposition(slope) for slope in S.slopes(bound=1)))  # optional: pyflatsurf
    ##         sage: c0, c1 = dec.components() # optional: pyflatsurf
    ##         sage: kz = O.flow_decomposition_kontsevich_zorich_cocycle(dec) # optional: pyflatsurf
    ##         sage: O.cylinder_circumference(c0, *kz) # optional: pyflatsurf
    ##         (1, 1, 1, -1)
    ##         sage: O.cylinder_circumference(c1, *kz) # optional: pyflatsurf
    ##         (0, 0, 1, 0)

    ##     """
    ##     # TODO: Change the returned value to be an actual relative homology class.
    ##     # TODO: Do we really need the weird input from kontsevich_zorich_cocycle?
    ##     if (
    ##         component.cylinder() != True
    ##     ):  # we are comparing to a boost tribool so this cannot be replaced by "is not True"  # noqa
    ##         raise ValueError

    ##     perimeters = list(component.perimeter())
    ##     per = perimeters[0]
    ##     assert not per.vertical()
    ##     sc = per.saddleConnection()
    ##     i = sc_index[sc]
    ##     if i < 0:
    ##         s = -1
    ##         i = -i - 1
    ##     else:
    ##         s = 1
    ##     v = s * proj.column(i)
    ##     circumference = -A.solve_right(v)

    ##     # check
    ##     hol = self.holonomy_dual(circumference)
    ##     holbis = self._vector_space_conversion().section(component.circumferenceHolonomy())
    ##     assert hol == holbis, (hol, holbis)

    ##     return circumference
    ##
    ## def cylinder_circumferences(self, decomposition):
    ##     # TODO: Return relative homology classes. (Though they actually lift to elements of absolute homology which we could assert maybe.)
    ##     # TODO: This is just a base change applied to the relative homology induced by decomposition. Maybe we could model it there trivially.
    ##     kz = self.flow_decomposition_kontsevich_zorich_cocycle(decomposition)

    ##     vcyls = []

    ##     for component in decomposition.components():
    ##         if (
    ##             component.cylinder() == False
    ##         ):  # we are comparing to a boost tribool so this cannot be replaced by "is False"  # noqa
    ##             continue
    ##         elif (
    ##             component.cylinder() == True
    ##         ):  # we are comparing to a boost tribool so this cannot be replaced with "is True"  # noqa
    ##             vcyls.append(self.cylinder_circumference(component, *kz))

    ##         else:
    ##             return []

    ##     return vcyls

    ## TODO: Expose this.
    ## def cylinder_module(self, cylinder):
    ##     r"""
    ##     Return the modulus of ``cylinder``, i.e., its width and height.

    ##     EXAMPLES::

    ##         # TODO: Add example.

    ##     """
    ##     section  = self._vector_space_conversion().ring_conversion().section
    ##     width = section(cylinder.width())
    ##     height = section(cylinder.vertical().project(cylinder.circumferenceHolonomy()))

    ##     return width, height


class FlowComponent_mapped(FlowComponent_base):
    def __init__(self, parent: FlowDecomposition_base, component: FlowComponent_base, morphism):
        super().__init__(parent)

        self._component = component

    def is_cylinder(self):
        return self._component.is_cylinder()

    def is_minimal(self):
        return self._component.is_minimal()

    def is_undetermined(self):
        return self._component.is_undetermined()
