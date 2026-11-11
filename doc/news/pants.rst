**Added:**

* Added ``hyperbolic_surfaces.pants()``, the hyperbolic pair of pants with
  given cuffs, by the traces of their holonomy, and cone points of angle less
  than ``2π``, by their orders, so that ``pants(**data).signature()`` gives
  ``data`` back. The front polygon is moved off the origin when this makes
  its coordinates rational, e.g. the pants with three cuffs of trace ``21/4``
  is defined over ``QQ``.
* Added ``hyperbolic_surfaces.regular_octagon()``, the hyperbolic surface of
  genus 2 built from the regular octagon of angle ``π/4``, with opposite
  edges glued or with the identification of Katok.

