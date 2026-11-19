**Added:**

* Added ``is_reflection()`` and ``is_glide_reflection()`` to similarities of the plane, to recognize orientation-reversing isometries.

* Added ``reflection()`` in ``flatsurf.geometry.similarity``, the affine matrix of the reflection in the line through two points, and the keyword ``orientation_preserving`` of ``similarity_from_vectors()`` to get the orientation-reversing similarity that maps one vector to another.

* Added orientation-reversing isometries (reflections and glide reflections) to ``HyperbolicIsometryGroup``, which is now the group ``O(Q)⁺`` instead of ``SO(Q)⁺``, with ``det()``, ``is_orientable()``, ``is_reflection()``, ``is_glide_reflection()`` and ``axis()`` of hyperbolic isometries, the geodesic preserved by a hyperbolic isometry, a glide reflection or a reflection.

* Added ``SimilarityGroup.subgroup()``, and ``is_orientable()`` and ``orientation_preserving_subgroup()`` to the finitely generated subgroups of ``SimilarityGroup`` and ``HyperbolicIsometryGroup``, which share the class ``FinitelyGeneratedSubgroup``; ``invariant_trace_field()`` of a subgroup of hyperbolic isometries with a reversing generator is read on its subgroup of orientation preserving elements.

* Added ``change_ring()`` to ``SimilarityGroup`` and ``HyperbolicIsometryGroup``.

**Changed:**

* Changed calling an orientation-reversing similarity on a polygon to return its image, pointed at the image of vertex 0 with its vertices counterclockwise, as ``matrix * P`` already did, instead of raising a ``ValueError``.

* Changed ``SimilarityGroup`` to accept 3×3 affine matrices, which raised a ``NotImplementedError``, and to convert the entries of matrix input into its base ring.
