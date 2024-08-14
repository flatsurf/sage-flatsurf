**Added:**

* Added ``VeechGroup`` and ``AffineAutomorphismGroup`` (also accessible as ``veech_gruop()`` and ``affine_automorphism_group()`` on surfaces.)

* Added ``apply_matrix()`` to all oriented similarity surfaces. (We apply the matrix to all polygons and keep the gluings intact. Probably not the most meaningful operation for non-dilation surfaces but it can be useful while building surfaces.)

* Added ``delaunay_triangulate()`` to similarity surfaces which returns a morphisms to a Delaunay triangulation of the surface.

* Added ``delaunay_decompose()`` to similarity surfaces which returns a morphisms to a Delaunay cell decomposition of the surface. This method always has an optional parameter ``codomain`` which can be an existing cell decomposition, this effectively exposes the ``isomorphism()`` function of a ``FlatTriangulation`` in libflatsurf.

* Added ``ngens()`` to homology to get the Betti number of homology.

* Added ``hom()`` for homology to produce the morphism in homology induced by a morphism of a surface or the morphism in homology given by a matrix.

* Added ``pyflatsurf()`` to translation surfaces to produce an isomorphic surface that is backed by libflatsurf (and exposes the underlying libflatsurf object in its ``flat_triangulation()`` method.)

**Changed:**

* Changed ``apply_matrix()``, it now returns a morphism to the deformed surface. To recover the old behavior, use ``apply_matrix().codomain()``.

* Changed ``triangle_flip()`` to now always turn the diagonal counterclockwise from the perspetive of ``label``.

**Deprecated:**

* Deprecated ``circumscribing_circle()`` of a polygon in favor of ``circumscribed_circle()``.

* Deprecated the ``test`` keyword of ``triangle_flip`` since there is already ``is_convex(strict=True)``.

* Deprecated ``delaunay_triangulation()`` in favor of ``delaunay_triangulate()``.

* Deprecated ``delaunay_decomposition()`` in favor of ``delaunay_decompose()``.

**Removed:**

* Removed the ``direction`` keyword of ``triangle_flip``, ``veering_triangulation``, ``l_infinity_delaunay_triangulation``, ``delaunay_triangulation``, and ``delaunay_decomposition`` since it was not clear what it actually did in general, instead the diagonal is now always turned counterclockwise.

**Fixed:**

* Fixed ``is_convex()`` for edges of surfaces. It does now allow this check for self-glued polygons.

* Fixed ``==`` and hashing for circles. Circles can now be used in dicts and sets.

**Performance:**

* <news item>
