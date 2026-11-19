**Added:**

* Added ``apply_transformation()`` to Euclidean and hyperbolic polygons, which
  moves a polygon by a similarity, which may reverse the orientation, or by an
  orientation preserving hyperbolic isometry (or anything the corresponding
  group converts), keeping the base ring of the polygon.
* Added ``apply_similarity()`` to Euclidean polygons, which moves a polygon by
  a similarity that may reverse the orientation, the counterpart of
  ``apply_isometry()`` for hyperbolic polygons.
* Added the product of a hyperbolic isometry and a hyperbolic polygon,
  ``h * P`` and ``P * h``.

**Fixed:**

* Fixed calling a similarity on a non-convex polygon, which raised an
  ``AttributeError``.
* Fixed the conversion of a similarity between similarity groups over
  different rings, which was always refused.
* Fixed the error raised by a hyperbolic isometry group when asked to convert
  an element it cannot convert.

