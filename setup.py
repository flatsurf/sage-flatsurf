# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2017-2019 Vincent Delecroix
#                     2021-2023 Julian Rüth
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
# ****************************************************************************
from setuptools import setup

with open("README.md") as f:
    long_description = f.read()

setup(
    name="sage_flatsurf",
    author="Vincent Delecroix, W. Patrick Hooper, and Julian Rüth",
    author_email="contact@flatsurf.org",
    description="Flat surfaces in SageMath",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version='0.5.0',
    url="https://github.com/flatsurf/sage-flatsurf",
    license="GNU General Public License, version 2",
    packages=[
        "flatsurf",
        "flatsurf.geometry",
        "flatsurf.geometry.categories",
        "flatsurf.graphical",
    ],
    install_requires=["surface-dynamics"],
    setup_requires=["wheel"],
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    keywords="surfaces, dynamics, geometry, flat surfaces, Abelian differentials, quadratic differentials, Riemann surfaces",
)
