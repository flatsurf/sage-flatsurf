from flatsurf.geometry.surface_objects import SurfacePoint_base


class SurfacePoint_pyflatsurf(SurfacePoint_base):
    def __init__(self, point, surface):
        self._point = point

        super().__init__(surface)
