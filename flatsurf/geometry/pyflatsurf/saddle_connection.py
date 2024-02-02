from flatsurf.geometry.saddle_connection import SaddleConnection_base


class SaddleConnection_pyflatsurf(SaddleConnection_base):
    def __init__(self, connection, surface):
        self._connection = connection
        self._surface = surface

    def surface(self):
        return self._surface
