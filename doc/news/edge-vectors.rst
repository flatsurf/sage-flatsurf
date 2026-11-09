**Added:**

* Added ``edge_vector()`` and ``edge_vectors()`` to Euclidean polygons, the vectors going from one vertex to the next.

* Added ``edge_segment()`` and ``edge_segments()`` to Euclidean polygons, the oriented segments going from one vertex to the next, and ``EuclideanOrientedSegment`` to represent them.

**Deprecated:**

* Deprecated ``edge()`` of Euclidean polygons; use ``edge_vector()`` instead.

* Deprecated ``edges()`` of Euclidean polygons; use ``edge_vectors()`` instead.

