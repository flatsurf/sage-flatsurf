**Added:**

* Added ``is_reflection()`` and ``is_glide_reflection()`` to similarities of the plane, to recognize orientation-reversing isometries.

* Added ``reflection()`` in ``flatsurf.geometry.similarity``, the affine matrix of the reflection in the line through two points, and the keyword ``orientation_preserving`` of ``similarity_from_vectors()`` to get the orientation-reversing similarity that maps one vector to another.

**Changed:**

* Changed calling an orientation-reversing similarity on a polygon to return its image, pointed at the image of vertex 0 with its vertices counterclockwise, as ``matrix * P`` already did, instead of raising a ``ValueError``.

* Changed ``SimilarityGroup`` to accept 3×3 affine matrices, which raised a ``NotImplementedError``, and to convert the entries of matrix input into its base ring.
