from flatsurf.geometry.cone_surface import (
    ConeSurface_generic, 
    ConeSurface_polygons_and_gluings,
    ConeSurface_wrapper)

from flatsurf.geometry.surface import SurfaceType, convert_to_type

class RationalConeSurface_generic(ConeSurface_generic):
    r"""
    A Euclidean cone surface whose cone angles are all rational multiples of pi.
    """
    def surface_type(self):
        return SurfaceType.RATIONAL_CONE

class RationalConeSurface_polygons_and_gluings(
        ConeSurface_polygons_and_gluings,
        RationalConeSurface_generic):
    pass

class RationalConeSurface_wrapper(
        ConeSurface_wrapper,
        RationalConeSurface_generic):
    pass

def convert_to_rational_cone_surface(surface):
    r"""
    Returns a cone surface version of the provided surface.
    """
    if surface.is_finite():
        return RationalConeSurface_polygons_and_gluings(surface)
    else:
        return RationalConeSurface_wrapper(surface)

