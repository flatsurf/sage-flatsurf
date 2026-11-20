**Added:**

* Added ``ramified_cover()`` for hyperbolic surfaces: it is now available
  to all oriented surfaces of finite type built from flat or hyperbolic
  polygons, and the cover of a hyperbolic surface stays in the hyperbolic
  plane of the surface, with its quadratic form.

**Fixed:**

* Fixed ``ramified_cover()`` on surfaces of several polygons: an edge of a
  sheet was glued to the sheet of the polygon of that edge instead of the
  sheet of the opposite polygon, which gave a wrong, often disconnected,
  cover.
* Fixed ``ramified_cover()`` on surfaces with self-glued edges (a
  ``KeyError``; the permutation of such an edge must now be an involution)
  and with unglued edges (a ``TypeError``; they now stay unglued in every
  sheet). The dictionary of permutations is no longer changed in place.

