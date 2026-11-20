**Added:**

* Added ``similarity_surfaces.pillowcase()``, which glues a Euclidean polygon
  to its mirror image in the line of one of its edges, along all of its edges
  or along some of them, the others forming the boundary. A non-convex
  polygon is not triangulated, and the surface keeps the base ring of the
  polygon.
* Added ``flatsurf.hyperbolic_surfaces`` with ``pillowcase()``, the same
  construction for a hyperbolic polygon, e.g. a pair of pants from a
  right-angled hexagon glued along alternate edges.

**Deprecated:**

* Deprecated ``similarity_surfaces.polygon_double()``; use ``pillowcase()``
  instead, which reflects the polygon in the line of an edge rather than in
  a horizontal line, and numbers the vertices of the mirror polygon
  differently for an even number of vertices.

