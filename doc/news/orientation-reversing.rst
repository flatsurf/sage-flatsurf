**Added:**

* Added ``is_reflection()`` and ``is_glide_reflection()`` to similarities of the plane, to recognize orientation-reversing isometries.

**Changed:**

* Changed calling an orientation-reversing similarity on a polygon to return its image, pointed at the image of vertex 0 with its vertices counterclockwise, as ``matrix * P`` already did, instead of raising a ``ValueError``.
