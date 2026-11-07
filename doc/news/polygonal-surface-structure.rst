**Added:**

* Added ``is_right_triangle()`` and ``is_isosceles_triangle()`` to Euclidean polygons.

**Changed:**

* Changed ``MutablePolygonalSurface`` to store the gluings of its polygons: ``glue()``, ``unglue()``, ``replace_polygon()`` and ``opposite_edge()`` are now available on it and no longer specific to ``MutableOrientedSimilaritySurface``.
