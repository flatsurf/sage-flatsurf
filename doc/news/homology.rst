**Added:**

* Added ``homology()`` and ``cohomology()`` for surfaces; also exposed as ``SimplicialHomology()`` and ``SimplicialCohomology()``.

**Removed:**

* Removed the old ``flatsurf.geometry.relative_homology`` without prior deprecation (since it was not exposed publicly anywhere.) The new implementation should cover all the relevant features.
