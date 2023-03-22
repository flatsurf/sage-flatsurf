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
    version="0.4.7",
    url="https://github.com/flatsurf/sage-flatsurf",
    license="GNU General Public License, version 2",
    packages=["flatsurf", "flatsurf.geometry", "flatsurf.graphical"],
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
