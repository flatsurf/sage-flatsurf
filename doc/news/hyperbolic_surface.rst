**Added:**

* Added ``is_right_triangle()`` and ``is_isosceles_triangle()`` to polygons.

* Added categories for hyperbolic surfaces, namely ``HyperbolicPolygonalSurfaces`` and ``HyperbolicIsometrySurfaces``.

* Added ``angle()`` method to points on a surface to return the total angle at a point.

* Added ``adjacencies()`` method to polygons that reports the edges adjacent to the vertices in order. This is necessary for hyperbolic polygons of infinite area where some vertices might not be adjacent to two edges.

**Changed:**

* Changed surface plotting quite a bit to make the same code work for hyperbolic surfaces and similarity surfaces. Edge labels can now only be configured differently for adjacent edges, non-adjacent edges, self-glued edges, and boundary edges. The new defaults there should make it much easier to interactively build a surface.

**Fixed:**

* Fixed merging of hyperbolic sets of vertices.

* Fixed ``vertices()`` of hyperbolic polygons with marked vertices in some edge cases.

* Fixed ``angle(numerical=True)`` of polygons which are now consistently returned as elements in RR.

* Fixed implementation of ``an_element()`` for testing to also create points on non-similarity surfaces.

* Fixed vectors returned by methods. The returned vectors are now immutable in more places.

* Fixed computation of pole of hyperbolic geodesics so it works without extending the base ring.

* Fixed handling of points at vertices of surfaces with boundary.

* Fixed complicated printing of hyperbolic segments. Now, segments are printed as ``start → end`` instead of giving the intersection of the half planes that define the segment.

**Deprecated:**

* Deprecated ``angles()`` method of surfaces.

**Removed:**

* Removed ``plot_polygon()`` on surfaces. All plotting is now in graphical surfaces. (This function can largely be emulated by hiding all but one polygon in a graphical surface before plotting.)
