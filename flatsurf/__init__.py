r"""
sage-flatsurf: Sagemath module for similitude surfaces
"""
from __future__ import absolute_import, print_function

from .version import version as __version__

from .geometry.polygon import polygons, EquiangularPolygons, Polygons, ConvexPolygons

from .geometry.similarity_surface_generators import (similarity_surfaces,
        dilation_surfaces, half_translation_surfaces, translation_surfaces)

from .geometry.surface import Surface_list, Surface_dict

# The various surface types
from .geometry.similarity_surface import SimilaritySurface
from .geometry.half_dilation_surface import HalfDilationSurface
from .geometry.dilation_surface import DilationSurface
from .geometry.cone_surface import ConeSurface
from .geometry.rational_cone_surface import RationalConeSurface
from .geometry.half_translation_surface import HalfTranslationSurface
from .geometry.translation_surface import TranslationSurface

try:
    import pyflatsurf
except ImportError:
    pass
else:
    from .geometry.gl2r_orbit_closure import GL2ROrbitClosure

del absolute_import, print_function
