**Added:**

* Added ``slopes()`` and ``_decomposition()`` to all translation surfaces; these were originally only implemented on ``GL2ROrbitClosure``. Note that ``_decomposition()`` does not return a proper sage-flatsurf object yet but a libflatsurf ``FlowDecomposition`` defined over the corresponding flat triangulation.
* Added ``vector_space_conversion()`` and ``ring_conversion()`` to pyflatsurf backed surfaces to get direct access to the underlying map from SageMath to libflatsurf objects.

**Changed:**

* Changed return types of methods of GL2ROrbitClosure to always return SageMath objects, e.g., instead of return e-antic number field elements, they now return SageMath number field elements. This change might break calling code if it is too reliant on the exact types returned.

**Deprecated:**

* Deprecated construction of GL2ROrbitClosure from pyflatsurf flat triangulations, instead all orbit closures should be created from actual sage-flatsurf surfaces.
* Deprecated flow decomposition machinery on GL2ROrbitClosure, i.e., ``decomposition()``, ``decompositions()``, ``decompositions_depth_first()``, and ``decompositions_breadth_first()``.

**Removed:**

* <news item>

**Fixed:**

* Fixed error when calling ``FlatTriangulationConversion.vector_space_conversion`` for some exact-real surfaces.
* Fixed ``pyflatsurf()`` for pyflatsurf backed surfaces.

**Performance:**

* Improved computation of orbit closures of large surfaces.
