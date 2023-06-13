**Added:**

* Added ``point()`` on surfaces which contains the features of both, ``singularity()`` and ``surface_point()``. (When passed a vertex id, it returns the point corresponding to that "singularity", when passed coordinates, it creates the point from the coordinates.)

* Added a hierarchy of surface categories, ``TopologicalSurfaces``, ``PolygonalSurfaces``, ``SimilaritySurfaces``, ``DilationSurfaces``, ``ConeSurfaces``, ``TranslationSurfaces``, ``HyperbolicSurfaces``, and ``RealProjectivePolygonalSurfaces`` together with axioms ``Connected``, ``Orientable``, ``Oriented``, ``WithBoundary``, ``WithoutBoundary``, ``FiniteType``, ``InfiniteType``, ``Positive``. These categories replace the existing hierarchy of ``SimilaritySurface``, ``DilationSurface``, ``HalfDilationSurface``, ``TranslationSurface``, ``HalfTranslationSurface``, ``RationalConeSurface``, ``RationalSimilaritySurface``, ``ConeSurface``. They serve essentially the same purpose, providing a place where to put functionality that applies to a certain kind of surface. However, this allows for more granularity, e.g., if some computation is only possible for translation surfaces of finite types, the method now lives in ``TranslationSurfaces.FiniteType`` and simply won't be available to surfaces that are of infinite type.

* Added a hierarchy of polygon categories, ``Polygons``, ``HyperbolicPolygons``, ``RealProjectivePolygons`` and ``RealProjectivePolygonsWithAngles`` together with axioms ``Convex``, ``Rational``, and ``Simple``. These replace the existing differentiation between ``Polygon``, ``ConvexPolygon`` and ``EquiangularPolygons``.

* Added a ``labels()`` method to surfaces. This replaces the deprecated ``label_iterator()`` method. The object returned by ``labels()`` can (often) be efficiently queried for its ``len`` and decide containment.

* Added a ``polygons()`` method to surfaces. This essentially replaces the deprecated ``label_polygon_iterator()`` and ``label_iterator(polygons=True)`` method. The object returned by ``polygons()`` can be efficiently queried for its ``len`` so this also replaces the deprecated ``num_polygons``.

* Added a ``roots()`` method to surfaces that returns labels from which iteration by ``labels()`` should start exploration of the connected components of a surface. There is also a ``root()`` method that returns the only such label on a connected surface.

* Added a ``is_finite_type()`` method to surfaces. This replaces the deprecated ``is_finite()``.

* Added an ``is_compact()`` method to surfaces that returns whether the surface is compact as a topological space.

* Added an ``is_connected()`` method to surfaces that returns whether the surface is connected as a topological space. (Before such surfaces were not well supported.)

* Added an ``is_with_boundary()`` method that returns whether a surfaces has polygons with unglued edges. (Before such surfaces were considered to be invalid now they are supported to some limited extent.)

* Added an ``euler_characteristic()`` method for surfaces of finite type.

* Added a ``change_ring()`` method to all surfaces to create a copy of the surface with polygons defined over a different base ring.

* Added a ``vertices()`` method to all surfaces built from polygons that returns the set of equivalence classes of vertices of those polygons.

* Added points of polygons as explicit objects (so that polygons become proper SageMath parents.) These points currently do not have many features.

* Added a ``describe_polygon()`` method to all polygons to create nicer textual representation of polygons for error messages.

* Added a ``is_degenerate()`` method to detect polygons with zero area or marked vertices.

* Added a ``marked_vertices`` keyword to ``vertices()`` of polygons to control whether vertices with angle π are included in the output.

* Added a ``erase_marked_vertices()`` method to polygons to produce a copy of the polygon with any vertices with angle π.

* Added ``is_equilateral()`` and ``is_equiangular()`` methods to polygons.

**Changed:**

* Changed the notion of "dilation surface" in some places. What was previously called a "half-dilation surface" is now called a "dilation surface", and what was previously called a "dilation surface" is now called a "positive dilation surface". (Existing code should not be affected by this but the documentation has been updated and new functions use this naming.) Note that the notions of a half-translation surface and a translation surface have not changed (though internally a translation surface is just a positive half-translation surface.)

* Changed the mutability of the surface returned by ``polyhedron_to_cone_surface``; the returned surface is now immutable.

* Changed the mutability of methods taking an ``in_place`` keyword argument. These methods now consistentlny return an immutable surface when ``in_place`` is ``False``.

* Changed ``Polygon`` in ``flatsurf.geometry.polygon``; it has been renamed to ``EuclideanPolygon``. ``Polygon()`` is now a factory function that creates an euclidean polygon, compatible with the way that ``polygons()`` used to create such a polygon.

* Changed the structure of hyperbolic sets. Points are now SageMath elements of their parent, the hyperbolic plane. Other sets are now (facade) parents. This change allows us to bring hyperbolic convex polygons as parents into the category of polygons. So, hyperbolic polygons and Euclidean polygons are on a similar footing and we can (eventually) build surfaces from both of them in the same way.

* Renamed ``flatsurf.geometry.matrix_2x2`` to ``flatsurf.geometry.euclidean`` and moved more geometry helpers there.

* Changed placement and naming of some geometry helper functions. For example, unified ``is_parallel`` and ``is_same_direction``, renamed ``is_opposite_direction`` to ``is_anti_parallel`` (and simplified implementation), renamed ``wedge_product`` to ``ccw`` and renamed ``wedge`` to ``wedge_product``. Most of these now live in ``flatsurf.geometry.euclidean``.

* Changed the label parameter of ``add_polygon()``. It is now required to be a keyword argument.

**Deprecated:**

* Deprecated the ``walker()`` method on surfaces. The ``labels()`` are now always guaranteed to be iterated in a canonical order (starting from the ``roots()``, a breadth-first search is performed.)

* Deprecated the ``base_label()`` method on surfaces. The ``root()`` and ``roots()`` methods serve the same purpose but have clearer semantics for disconnected surfaces.

* Deprecated the ``num_polygons()`` method on surfaces. For finite type surfaces, ``len(polygons())`` serves the same purpose (and is sufficiently fast.)

* Deprecated the ``label_iterator()``, ``edge_iterator()``, ``label_polygon_iterator()``, ``edge_gluing_iterator()`` and ``polygon_iterator()`` methods on surfaces. The ``_iterator`` suffix has always been confusing to some. Also, returning an iterator has limitations, e.g., containment and length cannot be queried easily.

* Deprecated the ``is_finite()`` method on surfaces. Since surfaces are now in the category of sets, ``is_finite()`` is also understood to answer whether the surface is a finite set of points. Eventually, ``is_finite()`` will change to return False for all non-empty surfaces.

* Deprecated the ``field()`` method on polygons since its semantics were a bit confusing. (Does it return the fraction field of the base ring or complain if the base ring is not a field?)

* Deprecated calling the object returned by ``Polygons()`` and ``ConvexPolygons()`` since they do not transfer well to the category framework (subcategories do not inherit ``__call__`` and ``Polygon()`` seems to be a more convenient alternative anyway.)

* Deprecated ``convexity()`` and ``strict_convexity()`` from ``EquiangularPolygons`` in favor of ``is_convex(strict=False)`` that is identical to the convexity method on polygons.

* Deprecated ``module()`` and ``vector_space()`` on polygons and ``EquiangularPolygons`` in favor of ``base_ring()**2``.

* Deprecated ``num_singularities()`` in favor of ``vertices()`` since the count does not distinguish between singularities and marked points.

* Deprecated the ``translation`` keyword argument of ``vertices()`` of a polygon; ``translate().vertices()`` seems to be the more straightforward approach.

* Deprecated ``is_strictly_convex()`` for polygons; it has been replaced with a ``strict`` keyword for ``is_convex()``.

* Deprecated ``num_edges()`` for polygons; it is essentially equivalent to ``len(vertices())`` (and "There should be one-- and preferably only one --obvious way to do it.")

* Deprecated implicitly iterating over the vertices of a polygon; this is problematic since a polygon is now the parent of its infinitely many points (and iterating over vertices() is easier to understand and equivalent anyway.)

* Deprecated ``add_polygons()`` for mutable surfaces; there is no benefit over adding polygons in a loop with ``add_polygon()``.

* Deprecated ``change_base_label()`` on surfaces; it has been replaced by ``set_root()`` and ``set_roots()`` to also support disonnected surfaces.

* Deprecated ``set_edge_pairing()`` and ``change_edge_gluing()`` for similarity surfaces; they have been replaced by ``glue()``.

* Deprecated ``change_polygon_gluings()``; it has a confusing syntax (and semantics) and using ``glue()`` in a loop does the same.

* Deprecated ``change_polygon()`` for surfaces; it had strange side effects in some cases; ``replace_polygon()`` should be easier to use.

**Removed:**

* Removed the ``cached`` parameter from ``.graphical_surface()`` of surfaces. The graphical surface returned is now never cached. If you want to customize the graphical surface returned you need to subclass the surface and add custom logic explicitly. (The "caching" that used to happen here made immutable surfaces in fact mutable and led to problems with serialization and equality testing in the past.)

* Removed ``flatsurf.geometry.xml``. It did not correctly serialize all kinds of surfaces and most likely nobody has been using it. If you relied on this functionality please let us know so we can bring it back in some form.

* Removed the ``limit`` keyword from ``delaunay_triangulation()``.

* Removed undocumented and untested method ``delaunay_single_join()`` from surfaces.

* Removed the ``_label_comparator`` from surfaces since it did not produce a consistent ordering on different architectures. There is now a ``min`` on labels, e.g., on a ``LabeledView`` which just uses the builtin ``min`` when it works and otherwise compares the ``repr`` of the labels. This approach also has problems (see documentation) but at least it is not platform dependent on the most common inputs such as polygons labeled by strings or numbers.

* Removed the untested ``standardize_polygons()`` for infinite type surfaces.

* Removed the possibility to ``canonicalize()`` a translation surface in-place. (This is a very expensive operation anyway and there does not seem to be a benefit to do this operation in-place.)

* Removed the ``n`` keyword argument in ``chamanara_surface(alpha, n)``. This keyword only affected plotting. It is ignored now and will be an error in a future version of sage-flatsurf.

* Removed the ``relabel`` argument in ``LazyTriangulatedSurface``, ``LazyDelaunayTriangulatedSurface``, and ``LazyDelaunaySurface``.

* Removed unused and untested ``translation_surface_cmp()`` from ``flatsurf.geometry.mappings``.

* Removed ``set_default_graphical_surface()``; if we allow this, we need to add the graphical surface to equality checks and hashing which is confusing.

**Fixed:**

* Fixed conversion of surfaces from flipper to sage-flatsurf.

* Fixed ``genus()`` for surfaces with self-glued edges.

* Fixed ``stratum()`` for half-translation surfaces with self-glued edges.

* Fixed (the deprecated) ``num_polygons()`` for disconnected surfaces.

* Fixed (most of the deprecated) ``*_iterator()`` methods for disconnected surfaces.

* Fixed ``subdivide_polygon(test=True)`` which sometimes returned ``None``.

* Fixed ``is_delaunay_triangulated()`` which now does not print to stdout anymore.

* Fixed ``is_delaunay_decomposed()`` which now checks not only the first polygon to decide whether a surface is Delaunay decomposed.

* Fixed ``standardize_polygons()`` which is now available on all similarity surfaces and not only on translation surfaces.

* Fixed ``LazyTriangulatedSurface``, ``LazyDelaunayTriangulatedSurface``, and ``LazyDelaunaySurface``. It is not necessary to walk the labels of such a surface before accessing the structure of the surface anymore.

* Fixed ``change_ring()`` for hyperbolic polygons. We do not forget about marked points when changing the base ring anymore.

* Fixed ``is_strictly_convex()`` for non-convex polygons.

**Performance:**

* Improved performance of computations related to angles (e.g., ``polygons.triangle(26, 48, 75)`` takes 50ms instead of 6s now,  asking that polygon for its angles is immediate now, ``%timeit similarity_surfaces.billiard(polygons.triangle(2, 13, 26)).minimal_cover("translation")`` takes 200ms instead of 15s now.)
