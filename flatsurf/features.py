r"""
Tests for optional packages used by sage-flatsurf.
"""

# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2021-2024 Julian Rüth
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
# ####################################################################

from sage.features import PythonModule

cppyy_feature = PythonModule(
    "cppyy", url="https://cppyy.readthedocs.io/en/latest/installation.html"
)


class PyeanticModule(PythonModule):
    def __init__(self):
        super().__init__(
            "pyeantic", url="https://github.com/flatsurf/e-antic/#install-with-conda"
        )

    @staticmethod
    def fix_unwrap_intrusive_ptr():
        r"""
        Unconditionally backports a fix from e-antic 2.0.1, see
        https://github.com/flatsurf/e-antic/pull/260.
        """
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import pyeantic
            import cppyy

            def unwrap_intrusive_ptr(K):
                if isinstance(K, pyeantic.eantic.renf_class):
                    K = cppyy.gbl.boost.intrusive_ptr["const eantic::renf_class"](K)

                if not isinstance(
                    K, cppyy.gbl.boost.intrusive_ptr["const eantic::renf_class"]
                ):
                    raise TypeError("argument must be an intrusive_ptr to a renf_class")

                wrapped = K.get()
                # pylint: disable=unused-private-member
                wrapped.__lifeline = K
                # pylint: enable=unused-private-member

                return wrapped

            import pyeantic.cppyy_eantic

            pyeantic.cppyy_eantic.unwrap_intrusive_ptr = unwrap_intrusive_ptr

            pyeantic.eantic.renf = lambda *args: unwrap_intrusive_ptr(
                pyeantic.eantic.renf_class.make(*args)
            )


pyeantic_feature = PyeanticModule()

pyexactreal_feature = PythonModule(
    "pyexactreal", url="https://github.com/flatsurf/exact-real/#install-with-conda"
)

pyflatsurf_feature = PythonModule(
    "pyflatsurf", url="https://github.com/flatsurf/flatsurf/#install-with-conda"
)

gmpxxyy_feature = PythonModule(
    "gmpxxyy", url="https://github.com/flatsurf/flatsurf/#install-with-conda"
)
