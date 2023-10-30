**Added:**

* Added ``bidict`` (providing dict encoding a bijection) as a dependency.

* Added ``apply_matrix`` to all oriented similarity surfaces. (We apply the matrix to all polygons and keep the gluings intact. Probably not the most meaningful operation for non-dilation surfaces but it can be useful while building surfaces.)

**Changed:**

* Changed return type of ``flow_to_exit()``. It now returns just the point where the flow exits a polygon and not a description of the point anymore. **This is a breaking change.**

* Changed ``billiard()`` to not triangulate non-convex polygons before creating the billiard. To restore the old behavior call ``triangulate()`` explicitly on the returned surface. Since surfaces built from non-convex polygons are quite limited, this might be a breaking change for some.

**Deprecated:**

* Deprecated ``polygon_double()`` since it is now identical with ``billiard()``.

* Deprecated ``triangulation_mapping()`` since ``triangulate()`` now returns a morphism.

* Deprecated ``canonicalize_mapping()`` (and moved it to the correct category) since ``canonicalize()`` now returns a morphism.

**Removed:**

* Removed ``flatsurf.geometry.mapping`` module. It has been replaced by ``flatsurf.geometry.morphism`` that can handle more general situations.

* Removed the ``mapping`` keyword from ``apply_matrix()`` of dilation surfaces. The method now always returns a morphism.

* Removed the previously deprecated ``rational`` keyword from ``billiard()``.

* Removed ability to ``apply_matrix(in_place=True)`` with matrix with negative determinant. If you rely on this feature for some reason, you can use ``MutableOrientedSimilaritySurface.from_surface(apply_matrix(in_place=False))`` instead.

**Fixed:**

* Fixed ``flow_to_exit()`` to also work for polygons that are not strictly convex.

**Performance:**

* <news item>
