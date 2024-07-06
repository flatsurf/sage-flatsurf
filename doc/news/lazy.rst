**Added:**

* Added ``triangulate()`` method to Euclidean polygons.

* Added ``__getitem__()`` to all ``labels()`` so that accessing ``surface.labels()[i]`` is possible.

* Added ``is_triangulated()``, ``is_delaunay_triangulated()``, and ``is_delaunay_decomposed()`` for the infinite staircase.

**Changed:**

* Renamed ``flatsurf.geometry.delaunay`` to ``flatsurf.geometry.lazy``. (If for some reason you used this module directly, you need to update your imports.)

* Moved ``GL2RImageSurface`` from ``flatsurf.geometry.half_dilation_surface`` to ``flatsurf.geometry.lazy``. (If for some reason you used this class directly, you need to update your imports.)

* Changed ``labels()`` and ``polygons()`` not to inherit from ``collections.abc.Set`` anymore. These are not just sets because their order matter, in particular, their order is compatible. However, ``edges()`` and ``gluings()`` are still ``collections.abc.Set``.

* Changed ``relabel()`` on surfaces to default to relabeling to integer labels. Also the keyword parameter ``relabeling_map`` is now called ``relabeling``.

**Deprecated:**

* Deprecated triangulation with ``triangulate(in_place=True)``. There is no performance advantage in such a triangulation and it complicates future work on morphisms.

* Deprecated the ``limit`` keywords for ``is_triangulated``, ``is_delaunay_triangulated``, ``_is_*_surface``, ``cmp``. Querying infinite surfaces for properties up to a limit is often not very useful and at the same time the question can be answered trivially with knowledge of the surfaces. Also, this was implemented very inconsistently.

* Deprecated the ``singularity_limit`` keywords for ``to_surface``. There is no replacement planned for this feature.

**Removed:**

* Removed the ``ring`` keyword argument from ``GL2RImageSurface`` and ``GL2RMapping``, the ring is now always the common parent of the surface base ring and the matrix base ring. (Use ``change_ring`` if you want the surface to live over another ring.)

* Removed the option ``in_place=True`` from ``delaunay_triangulation()`` and ``delaunay_decomposition()``. There is no asymptotic runtime advantage in performing this operation in place (and we need to maintain two very different copies of the same functionality.)

* Removed the option ``triangulated=True`` from ``delaunay_triangulation()`` and ``delaunay_decomposition()`` since there is no significant runtime advantage in practice.

* Removed the option ``delaunay_triangulated=True`` from ``delaunay_decomposition()`` since there is no significant runtime advantage in practice.

* Removed the ``direction`` keyword from ``delaunay_triangulation()`` and ``delaunay_decomposition()``. The claimed behavior about separatrices does not actually hold so the behavior here was somewhat random. Instead, we now always turn each edge counterclockwise when flipping.

* Removed ``rel_deformation()``. It should be replaced by the ``operator+`` that is implemented in libflatsurf and performs much better. If you need this method, let us know and we'll add something back quickly.

**Fixed:**

* Fixed comparison of infinite sets of labels. Labels can now be compared and hashed in some very limited cases.

**Performance:**

* <news item>
