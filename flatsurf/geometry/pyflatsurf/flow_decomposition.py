from flatsurf.geometry.flow_decomposition import FlowDecomposition_base


class FlowDecomposition_pyflatsurf(FlowDecomposition_base):
    def __init__(self, flow_decomposition, surface):
        super().__init__(surface)

        self._flow_decomposition = flow_decomposition
        self._surface = surface

    def surface(self):
        return self._surface
