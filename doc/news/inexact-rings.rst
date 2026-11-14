**Changed:**

* ``Polygon()`` and ``change_ring()`` of a Euclidean polygon raise a
  ``ValueError`` when the base ring is not exact, whether it is given or
  computed from float coordinates, as the hyperbolic polygon constructors do;
  so ``polygons.square(base_ring=RDF)`` is refused too.

