from flatsurf.geometry.flow_decomposition import FlowDecomposition_base, FlowComponent_base


def tribool_to_bool_or_none(tribool):
    if tribool:
        return True

    import cppyy
    if cppyy.gbl.boost.logic.indeterminate(tribool):
        return None

    return False


class FlowDecomposition_pyflatsurf(FlowDecomposition_base):
    def __init__(self, flow_decomposition, surface):
        super().__init__(surface)

        self._flow_decomposition = flow_decomposition
        self._components = None

    def decompose(self, limit=-1):
        if limit != 0:
            self.invalidate_components()
            self._flow_decomposition.decompose(int(limit))

    def is_parabolic(self):
        return tribool_to_bool_or_none(self._flow_decomposition.parabolic())

    def invalidate_components(self):
        for component in self._components or []:
            component._invalidate()

        self._components = None

    def components(self):
        if self._components is None:
            self._components = [FlowComponent_pyflatsurf(component, self) for component in self._flow_decomposition.components()]

        return self._components

    def has_cylinder(self):
        return tribool_to_bool_or_none(self._flow_decomposition.hasCylinder())

    def is_completely_periodic(self):
        return tribool_to_bool_or_none(self._flow_decomposition.completelyPeriodic())


class FlowComponent_pyflatsurf(FlowComponent_base):
    def __init__(self, flow_component, flow_decomposition):
        super().__init__(flow_decomposition)

        self.__flow_component = flow_component

    def _invalidate(self):
        self.__flow_component = None

    def _flow_component(self):
        if self.__flow_component is None:
            raise NotImplementedError("component of flow decomposition has been modified externally since it was created")
        return self.__flow_component

    def is_cylinder(self):
        return tribool_to_bool_or_none(self._flow_component().cylinder())
