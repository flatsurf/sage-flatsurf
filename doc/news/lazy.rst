**Added:**

* <news item>

**Changed:**

* Renamed ``flatsurf.geometry.delaunay`` to ``flatsurf.geometry.lazy``. (If for some reason you used this module directly, you need to update your imports.)
* Moved ``GL2RImageSurface`` from ``flatsurf.geometry.half_dilation_surface`` to ``flatsurf.geometry.lazy``. (If for some reason you used this class directly, you need to update your imports.)

**Deprecated:**

* <news item>

**Removed:**

* Removed the ``ring`` keyword argument from ``GL2RImageSurface``, the ring is now always the common parent of the surface base ring and the matrix base ring. (Use ``change_ring`` if you want the surface to live over another ring.)

* Removed the option ``in_place=True`` from ``delaunay_triangulation()`` and ``delaunay_decomposition()``. There is no asymptotic runtime advantage in performing this operation in place (and we need to maintain two very different copies of the same functionality.)

* Removed the option ``triangulated=True`` from ``delaunay_triangulation()`` and ``delaunay_decomposition()`` since there is no significant runtime advantage in practice.

* Removed the option ``delaunay_triangulated=True`` from ``delaunay_decomposition()`` since there is no significant runtime advantage in practice.

* Removed the ``direction`` keyword from ``delaunay_triangulation()`` and ``delaunay_decomposition()``. The claimed behavior about separatrices does not actually hold so the behavior here was somewhat random. Instead, we now always turn each edge counterclockwise when flipping.

**Fixed:**

* Fixed comparison of infinite sets of labels. Labels can now be compared and hashed in some very limited cases.

**Performance:**

* <news item>

