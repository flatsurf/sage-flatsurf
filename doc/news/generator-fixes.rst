**Added:**

* Added ``HyperbolicPolygon.edge_cosh_squared()``, the square of the
  hyperbolic cosine of the length of an edge, an element of the base ring of
  the polygon where the hyperbolic cosine itself often does not lie.
* Exported ``HyperbolicPolygon`` from the ``flatsurf`` package.

**Changed:**

* The hyperbolic generators with a ``base_ring`` argument
  (``hyperbolic_polygons.triangle()``, ``regular_ngon()``,
  ``right_angled_hexagon()``, ``right_angled_pentagon()``,
  ``right_angled_quadrilateral()``, ``hyperbolic_surfaces.pants()``,
  ``regular_octagon()``, ``from_signature()`` and
  ``fuchsian_groups.triangle()``) raise a ``ValueError`` when ``base_ring``
  is not exact.
* The hyperbolic polygon constructors raise a ``ValueError`` on an angle that
  is not a rational number, a float included.
* ``hyperbolic_polygons.regular_ngon()`` raises a ``ValueError`` instead of a
  ``TypeError`` when the number of vertices is not an integer.
* ``SurfaceSignature`` raises a ``ValueError`` when only some of the boundary
  traces are ``None``.

