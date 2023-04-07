**Changed:**

* Changed == of surfaces in some limited cases. Equality now means more strictly that surfaces are indistinguishable. For example, a mutable surface is now always distinct from an immutable surface.

**Fixed:**

* Fixed pickling for most infinite surfaces.
* Fixed missing base label in mutable surfaces. When the first polygon is added to the surface, it is set as the base label now.
