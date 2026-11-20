**Added:**

* Added ``hyperbolic_surfaces.pants()``, the hyperbolic pair of pants with
  given cuffs, by the traces of their holonomy, and cone points of angle less
  than ``2π``, by their orders, so that ``pants(**data).signature()`` gives
  ``data`` back. The pants lies in the hyperbolic plane of a quadratic form
  of its own, over the field of the ``t_k + 1`` of its cuffs of traces
  ``t_k``, of ``4 cos^2`` of the angles of its cone points and of the product
  ``p`` of their cosines and hyperbolic cosines, e.g. ``QQ`` for three cuffs
  of trace ``21/4``; ``standard=True`` carries it to the plane of the form
  ``x^2 + y^2 - z^2``.
* Added ``hyperbolic_surfaces.regular_octagon()``, the hyperbolic surface of
  genus 2 built from the regular octagon of angle ``π/4``, with opposite
  edges glued or with the identification of Katok.

