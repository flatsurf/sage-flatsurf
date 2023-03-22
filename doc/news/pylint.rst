**Removed:**

* Removed ``MatrixListDeformedSurface`` and ``MatrixListDeformedSurfaceMapping`` that have not been functional for many years.

**Fixed:**

* Fixed some corner cases in hyperbolic isometry detection.
* Fixed some cases of ``finitely_generated_matrix_group.matrix_multiplicative_order()``.
* Fixed checks in half translation surfaces with an infinite number of polygons.
* Fixed comparison operators of relative homology classes.
* Fixed error reporting in some corner cases of tangent bundles.

**Deprecated:**

* Deprecated (and disabled) the ``cache`` keyword argument of `SaddleConnection.trajectory()`. Trajectories are now always cached.
