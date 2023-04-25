r"""
sage-flatsurf: Sagemath module for similitude surfaces
"""
from .version import version as __version__

from .geometry.polygon import polygons, EquiangularPolygons, Polygons, ConvexPolygons

from .geometry.similarity_surface_generators import (
    similarity_surfaces,
    dilation_surfaces,
    half_translation_surfaces,
    translation_surfaces,
)

from .geometry.surface import Surface_list, Surface_dict, MutableOrientedSimilaritySurface

from .geometry.gl2r_orbit_closure import GL2ROrbitClosure

from .geometry.hyperbolic import HyperbolicPlane
