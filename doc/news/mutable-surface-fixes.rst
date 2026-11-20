**Fixed:**

* A flat surface no longer accepts a hyperbolic polygon in ``add_polygon()``
  and ``replace_polygon()``; a ``TypeError`` is raised.
* An empty flat surface is no longer equal to an empty hyperbolic surface over
  the same ring.

