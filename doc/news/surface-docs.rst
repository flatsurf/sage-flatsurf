**Changed:**

* The `copy` parameter of `Surface_list.__init__()` and
  `__Surface_dict.__init__()` now defaults to `surface.is_mutable()`. Before
  the default was `True`. However, in principle this should not break any
  existing code but only change the runtime slightly in some cases.

* The `mutable` parameter of `Surface_list.__init__()` and
  `Surface_dict.__init__()` now defaults to `True`. Before its default was
  `False` in many cases. This change might break some existing code. If it
  does, one needs to either explicitly set `mutable=True` in this invocation or
  call `surface.set_immutable()`.

**Fixed:**

* Fixed some issues in documentation of Surface classes and simplified some of their implementation.

* Fixed typos that lead to runtime errors in rare cases.
