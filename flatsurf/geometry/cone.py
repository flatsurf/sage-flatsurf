from sage.all import UniqueRepresentation, Parent
from sage.structure.element import Element
from sage.misc.cachefunc import cached_method


class Cone(Element):
    # This is an open cone.
    def __init__(self, parent, start, end):
        super().__init__(parent)
        self._start = parent.rays()(start)
        self._end = parent.rays()(end)

    def is_empty(self):
        return self._start == self._end

    def is_convex(self):
        from flatsurf.geometry.euclidean import ccw
        return ccw(self._start.vector(), self._end.vector()) >= 0

    # TODO: Rename to contains_cone() [arguments are reversed!]
    def is_subset(self, other):
        r"""
        Return whether this cone is contained in ``other``.
        """
        if not isinstance(other, Cone):
            raise NotImplementedError

        if other.parent() is not self.parent():
            raise NotImplementedError

        if self.is_empty():
            return True

        if other.is_empty():
            return False

        from flatsurf.geometry.euclidean import ccw
        start_ccw = ccw(self._start.vector(), other._start.vector())
        end_ccw = ccw(self._end.vector(), other._end.vector())

        # Check whether the interior of this cone is contained in other.
        if self.is_convex():
            if other.is_convex():
                if start_ccw > 0:
                    return False
                if end_ccw < 0:
                    return False
            else:
                raise NotImplementedError
        else:
            if other.is_convex():
                return False

            raise NotImplementedError

        return True

    def contains_ray(self, ray):
        if ray.parent() is not self.parent().rays():
            raise NotImplementedError

        if self.is_empty():
            return False

        from flatsurf.geometry.euclidean import ccw
        ccw_from_start = ccw(self._start.vector(), ray.vector())
        ccw_to_end = ccw(ray.vector(), self._end.vector())

        if self.is_convex():
            if ccw_from_start <= 0:
                return False

            if ccw_to_end <= 0:
                return False

            return True
        else:
            raise NotImplementedError

    def sorted_rays(self, rays):
        class Key:
            def __init__(self, cone, ray):
                self._cone = cone
                self._ray = ray

            def __lt__(self, rhs):
                if not self._cone.is_convex():
                    raise NotImplementedError
                from flatsurf.geometry.euclidean import ccw
                return ccw(self._ray.vector(), rhs._ray.vector()) > 0

        rays = sorted(rays, key=lambda ray: Key(self, ray))

        from itertools import groupby
        return [ray for ray, _ in groupby(rays)]

    def a_ray(self):
        if self.is_empty():
            raise TypeError

        if self.is_convex():
            return self.parent().rays()((self._start.vector() + self._end.vector()) / 2)

        raise NotImplementedError

    def start(self):
        return self._start

    def end(self):
        return self._end

    def contains_point(self, p):
        if self.is_empty():
            return False
        if p.is_zero():
            return True
        return self.contains_ray(self.parent().rays()(p))

    def _repr_(self):
        if self.is_empty():
            return "Empty cone"
        return f"Open cone between {self._start} and {self._end}"


class Cones(UniqueRepresentation, Parent):
    Element = Cone

    def __init__(self, base_ring, category=None):
        from sage.categories.all import Sets
        super().__init__(base_ring, category=category or Sets())

    @cached_method
    def rays(self):
        from flatsurf.geometry.ray import Rays
        return Rays(self.base_ring())
