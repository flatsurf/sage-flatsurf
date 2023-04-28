r"""
sage-flatsurf: Sagemath module for similitude surfaces
"""
from flatsurf.version import version as __version__

from flatsurf.geometry.polygon import polygon, polygons, EquiangularPolygons, Polygons, ConvexPolygons

from flatsurf.geometry.similarity_surface_generators import (
    similarity_surfaces,
    dilation_surfaces,
    half_translation_surfaces,
    translation_surfaces,
)

from flatsurf.geometry.surface import Surface_list, Surface_dict, MutableOrientedSimilaritySurface

from flatsurf.geometry.gl2r_orbit_closure import GL2ROrbitClosure

from flatsurf.geometry.hyperbolic import HyperbolicPlane
