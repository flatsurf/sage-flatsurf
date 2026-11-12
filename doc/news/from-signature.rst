**Added:**

* Added ``hyperbolic_surfaces.from_signature()``, a compact hyperbolic
  surface, possibly with geodesic boundary and cone points of angle less than
  ``2π``, of a given ``SurfaceSignature``, so that
  ``from_signature(s).signature() == s``. It is glued without twist from
  pairs of pants along a chain, the cone points of angle ``π`` paired into
  folds, every interior cuff of trace ``21/4`` unless ``cuff_traces`` is
  given. The base ring is the number field of the coordinates and side
  pairings, e.g. ``QQ`` for genus 2 and for the orbifolds ``(1; 2)``,
  ``(0; 2, 2, 2, 2, 2)`` and ``(0; 2, 2, 2, 3)``.

