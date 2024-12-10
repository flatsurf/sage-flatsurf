r"""
Adapts :mod:`sage.misc.cachefunc` for sage-flatsurf.

Some methods on surfaces in sage-flatsurf benefit a lot from caching. However,
we cannot simply put ``@cached_method`` on these surfaces since they are often
defined in a class that is used by both, mutable and immutable surfaces.
Caching for mutable surfaces would be incorrect, or we would have to manually
clear caches upon mutation.

Therefore, we add a decorator ``@cached_surface_method`` that only caches if
the surface is immutable.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface, translation_surfaces
    sage: S = MutableOrientedSimilaritySurface.from_surface(translation_surfaces.square_torus())
    sage: S.edge_matrix(0, 0) is S.edge_matrix(0, 0)
    False

When we call
:meth:`flatsurf.geometry.surface.MutablePolygonalSurface.set_immutable`,
caching is enabled for this method::

    sage: S.set_immutable()
    sage: S.edge_matrix(0, 0) is S.edge_matrix(0, 0)
    True

"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2024 Julian Rüth
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# *********************************************************************

from sage.misc.cachefunc import (
    CachedMethod,
    CachedMethodCaller,
    CachedMethodCallerNoArgs,
)
from sage.misc.decorators import decorator_keywords


@decorator_keywords
def cached_surface_method(f, name=None, key=None, do_pickle=None):
    r"""
    Make a surface method cached as soon as the surface is immutable.

    This is a drop-in replacement for SageMath's
    ``sage.misc.cachefunc.cached_method`` but it does not provide any caching
    for mutable surfaces.
    """
    from sage.misc.cachefunc import CachedMethod

    fname = name or f.__name__
    if fname.startswith("__"):
        # Copy the special names and wrap CachedSpecialMethod to make this work.
        raise NotImplementedError("cannot cache special methods of surfaces yet")

    return CachedSurfaceMethod(f, name=name, key=key, do_pickle=do_pickle)


class CachedSurfaceMethodCaller(CachedMethodCaller):
    r"""
    Adapts the ``sage.misc.cachefunc.CachedMethodCaller`` from SageMath to not
    do any caching if the defining surface is mutable.
    """

    def __init__(self, cached_caller, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._cached_caller = cached_caller

    def __call__(self, *args, **kwargs):
        r"""
        Invoke the underlying method with the given arguments.

        The underlying ``sage.misc.cachefunc.CachedMethodCaller`` is written in
        Cython and highly optimized. This method is much slower (by about 250ns
        per call). However, once the surface is immutable, this method is
        replaced with the underlying caching one, so the overhead is only paid
        once.
        """
        if self._instance.is_mutable():
            return self._instance_call(*args, **kwargs)

        setattr(self._instance, self.__name__, self._cached_caller)

        return self._cached_caller.__call__(*args, **kwargs)


class CachedSurfaceMethodCallerNoArgs(CachedMethodCallerNoArgs):
    def __init__(self, cached_caller, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._cached_caller = cached_caller

    def __call__(self):
        r"""
        Invoke the underlying method with the given arguments.

        The underlying ``sage.misc.cachefunc.CachedMethodCallerNoArgs`` is written in
        Cython and highly optimized. This method is much slower (by about 250ns
        per call). However, once the surface is immutable, this method is
        replaced with the underlying caching one, so the overhead is only paid
        once.
        """
        if self._instance.is_mutable():
            return self._instance_call()

        setattr(self._instance, self.__name__, self._cached_caller)

        return self._cached_caller.__call__()


class CachedSurfaceMethod(CachedMethod):
    r"""
    Customizes a method in a class so that it conditionally enables caching
    just like SageMath's ``sage.misc.cachefunc.CachedMethod`` does.
    """

    def __init__(self, f, name, key, do_pickle):
        super().__init__(f, name, key, do_pickle)
        self._f = f
        self._key = key
        self._do_pickle = do_pickle

    def __get__(self, inst, cls):
        try:
            caller = super().__get__(inst, cls)

            if inst is not None and inst.is_mutable():
                # If the surface is mutable, we replace the default cached caller
                # with a caller that considers mutability (since the surface might
                # later become immutable).
                if isinstance(caller, CachedMethodCallerNoArgs):
                    caller = CachedSurfaceMethodCallerNoArgs(
                        caller,
                        inst,
                        self._f,
                        name=caller.__name__,
                        do_pickle=self._do_pickle,
                    )
                else:
                    caller = CachedSurfaceMethodCaller(
                        caller,
                        self,
                        inst,
                        cache=self._get_instance_cache(inst),
                        name=caller.__name__,
                        key=self._key,
                        do_pickle=self._do_pickle,
                    )

            if inst is not None:
                setattr(inst, caller.__name__, caller)

            return caller
        except Exception:
            import traceback

            traceback.print_exc()
