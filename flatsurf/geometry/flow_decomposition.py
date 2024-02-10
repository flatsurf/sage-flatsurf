from sage.all import SageObject


class FlowDecomposition_base(SageObject):
    def __init__(self, surface):
        self._surface = surface
