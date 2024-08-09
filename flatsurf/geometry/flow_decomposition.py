from sage.all import SageObject


class FlowDecomposition_base(SageObject):
    def __init__(self, surface):
        self._surface = surface

    def surface(self):
        return self._surface


class FlowComponent_base(SageObject):
    def __init__(self, flow_decomposition):
        self._flow_decomposition = flow_decomposition

    def flow_decomposition(self):
        return self._flow_decomposition
