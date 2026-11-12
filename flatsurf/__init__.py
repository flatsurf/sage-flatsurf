r"""
Flat Surfaces in SageMath
"""

from flatsurf.version import version as __version__

from flatsurf.geometry.polygon import (
    Polygon,
    polygons,
    EquiangularPolygons,
    EuclideanPolygonsWithAngles,
    EuclideanPolygons as Polygons,
    ConvexPolygons,
)

from flatsurf.geometry.similarity_surface_generators import (
    similarity_surfaces,
    dilation_surfaces,
    half_translation_surfaces,
    translation_surfaces,
)

from flatsurf.geometry.hyperbolic_surface_generators import (
    hyperbolic_polygons,
    hyperbolic_surfaces,
    fuchsian_groups,
)

from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

from flatsurf.geometry.gl2r_orbit_closure import GL2ROrbitClosure

from flatsurf.geometry.hyperbolic import HyperbolicPlane
from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon
from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
from flatsurf.geometry.signature import SurfaceSignature

from flatsurf.geometry.homology import SimplicialHomology
from flatsurf.geometry.cohomology import SimplicialCohomology

from flatsurf.geometry.veech_group import AffineAutomorphismGroup, VeechGroup

from flatsurf.geometry.surface_legacy import (
    Surface_list,
    Surface_dict,
    SimilaritySurface,
    HalfDilationSurface,
    DilationSurface,
    ConeSurface,
    RationalConeSurface,
    HalfTranslationSurface,
    TranslationSurface,
)
