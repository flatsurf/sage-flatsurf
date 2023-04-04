**Added:**

* Added some convenience functions for ``SurfacePoint`` such as ``representatives()``, ``representative()``. Also ported over the methods from ``Singularity`` to ``SurfacePoint`` (though most are marked deprecated.)

**Changed:**

* The parameter ``limit`` for `SurfacePoint` and `Singularity` is now optional. If not given for infinite surfaces, the search will keep going until the object has been constructed.

* Changed how points print. Instead of `Surface point with n coordinate representations`, points now print as one of these coordinate representations.

**Deprecated:**

* Deprecated `Singularity` and most of its methods in favor of the more generic `SurfacePoint`.

**Performance:**

* Some operations of ``Singularity`` and ``SurfacePoint`` used to be linear and are now constant time and vice versa. This change should not be noticeable on most surfaces. Things could be easily sped up if this is a problem for some applications.
