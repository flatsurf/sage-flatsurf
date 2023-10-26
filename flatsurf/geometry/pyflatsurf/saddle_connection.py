class SaddleConnection_pyflatsurf:
    def __init__(self, connection, surface):
        self._connection = connection
        self._surface = surface

    def surface(self):
        return self._surface
