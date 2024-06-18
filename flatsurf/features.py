r"""
Tests for optional packages used by sage-flatsurf.
"""
# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2021-2023 Julian Rüth
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

pyeantic_feature = PythonModule(
    "pyeantic", url="https://github.com/flatsurf/e-antic/#install-with-conda"
)

pyexactreal_feature = PythonModule(
    "pyexactreal", url="https://github.com/flatsurf/exact-real/#install-with-conda"
)


class PyflatsurfModule(PythonModule):
    def __init__(self):
        super().__init__(
            "pyflatsurf", url="https://github.com/flatsurf/flatsurf/#install-with-conda"
        )

    def is_saddle_connection_enumeration_functional(self):
        # TODO: Check whether version is >=3.14.1
        return False


pyflatsurf_feature = PyflatsurfModule()
