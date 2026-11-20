**Changed:**

* Renamed the category ``ConeSurfaces`` to ``FlatConeSurfaces``; its repr is
  now "Category of flat cone surfaces".
* ``ConeSurfaces`` is now the category common to flat and hyperbolic cone
  surfaces. Code that used ``ConeSurfaces`` for flat surfaces should use
  ``FlatConeSurfaces`` instead.
* The path ``flatsurf.geometry.categories.cone_surfaces.ConeSurfaces`` is a
  deprecated alias of ``FlatConeSurfaces``, so that stored pickles of surfaces
  made by earlier releases load; the common category lives in
  ``flatsurf.geometry.categories.common_cone_surfaces``.

