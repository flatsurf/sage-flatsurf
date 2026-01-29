**Added:**

* Added ``slopes()`` and ``flow_decomposition()`` to all translation surfaces; these were originally only implemented on ``GL2ROrbitClosure``.
* Added ``chains()`` to all translation surfaces and added support for mapping chains similar to the existing machinery on homology classes.
* Added ``holonomy()`` to homology classes to replace ``holonomy()`` on ``GL2ROrbitClosure``.
* Added ``vector_space_conversion()`` and ``ring_conversion()`` to pyflatsurf backed surfaces to get direct access to the underlying map from SageMath to libflatsurf objects.

**Changed:**

* Changed return types of methods of GL2ROrbitClosure to always return SageMath objects, e.g., instead of returning e-antic number field elements, they now return SageMath number field elements. This change might break calling code if it is too reliant on the exact types returned.

**Deprecated:**

* Deprecated construction of GL2ROrbitClosure from pyflatsurf flat triangulations, instead all orbit closures should be created from actual sage-flatsurf surfaces.
* Deprecated flow decomposition machinery on GL2ROrbitClosure, i.e., ``decomposition()``, ``decompositions()``, ``decompositions_depth_first()``, and ``decompositions_breadth_first()``.
* Deprecated (co)homology machinery on GL2ROrbitClosure, i.e., ``holonomy()``, ...

**Removed:**

* <news item>

**Fixed:**

* Fixed error when calling ``FlatTriangulationConversion.vector_space_conversion`` for some exact-real surfaces.
* Fixed ``pyflatsurf()`` for pyflatsurf backed surfaces.
* Fixed homology classes to check whether they are actually cycles and not just an arbitrary chain. (This is a side effect of using an homology basis to represent homology classes. As a side effect, one cannot build homology classes by doing ``H((label, edge)) + H((label', edge'))`` if the summands are not cycles. Instead, one needs to build this from chains, e.g., ``C = H.chain; H(C((label, edge)) + +C((label', edge')))``.
* Fixed multiplication of homomorphisms in homology by scalars.

**Performance:**

* Improved computation of orbit closures of large surfaces.
* Improved internal representation of homology classes to store coefficients in a basis of homology instead of a basis of the chain module.
