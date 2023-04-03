**Fixed:**

* `Singularity`, `SurfacePoint`, `SaddleConnection` do not throw an exception anymore when compared to a different kind of object or to an object defined on another surface. (Instead, they now return to be non-equal in such cases.)
