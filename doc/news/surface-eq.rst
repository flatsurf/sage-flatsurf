**Fixed:**

* Fixed comparison of `Surface` for equality in some cases. Equality comparison does not throw exceptions in as many cases anymore in particular surfaces can now be compared to non surfaces and thus be put into sets and dicts (they will return that they are not equal to a non-surface.)
