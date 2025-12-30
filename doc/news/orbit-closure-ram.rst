**Added:**

* <news item>

**Changed:**

* Changed return types of methods of GL2ROrbitClosure to always return SageMath objects, e.g., instead of return e-antic number field elements, they now return SageMath number field elements. This change might break calling code if it is too reliant on the exact types returned.

**Deprecated:**

* Deprecated construction of GL2ROrbitClosure from pyflatsurf flat triangulations, instead all orbit closures should be created from actual sage-flatsurf surfaces.

**Removed:**

* <news item>

**Fixed:**

* Fixed error when calling ``FlatTriangulationConversion.vector_space_conversion`` for some exact-real surfaces.

**Performance:**

* Improved computation of orbit closures of large surfaces.
