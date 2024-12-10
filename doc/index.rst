sage-flatsurf: Flat Surfaces in SageMath
========================================

sage-flatsurf is a Python package for working with flat surfaces in `SageMath
<https://sagemath.org>`_.

We aim for sage-flatsurf to support the investigation of geometric, algebraic
and dynamical questions related to flat surfaces. By flat surface we mean a
surface modeled on the plane with monodromy given by similarities of the plane,
though current efforts are focused on `translation surfaces
<https://en.wikipedia.org/wiki/Translation_surface>`_ and `half-translation
surfaces
<https://en.wikipedia.org/wiki/Translation_surface#Half-translation_surfaces>`_.

Installation
============

The preferred way to install software should be to use your package manager
(e.g. ``apt`` on Debian or Ubuntu, ``pacman`` on Arch Linux, ``brew`` on macOS,
etc). However, as of this writing, sage-flatsurf has not been picked up by `any
of the major distributions yet <https://repology.org/project/python:sage-flatsurf/packages>`_.

We therefore recommend to install sage-flatsurf using our installers from the
`Releases Page <https://github.com/flatsurf/sage-flatsurf/releases>`_.

Detailed installation instructions:

* :ref:`Install with the pixi tarball <installation-tarball>` for Linux and macOS (recommended)
* :ref:`Install with the installer <installation-installer>` for Windows (recommended)
* :ref:`Install with Conda <installation-conda>`
* :ref:`Install into an existing SageMath source build <installation-sagemath>`
* :ref:`Install with pip <installation-pip>`

If you are planning to work on the sage-flatsurf source code, you can also
build sage-flatsurf from source. For this, please have a look at our
:ref:`Developer's Guide <developers-guide>`.

.. toctree::
   :hidden:

   install
   developer

A Tour of sage-flatsurf
=======================

Demos of some of the capabilities of sage-flatsurf:

.. toctree::
   :maxdepth: 1

   examples/tour
   examples/defining_surfaces
   examples/graphics_configuration
   examples/linear_action_and_delaunay
   examples/rel_deformations
   examples/saddle_connections
   examples/apisa_wright
   examples/siegel_veech
   examples/straight_line_flow
   examples/warwick-2017
   examples/boshernitzan_conjecture

Module Reference
================

The sage-flatsurf source code is split into two packages: ``graphical`` which
contains plotting logic, and ``geometry`` containing everything else. The links
below lead to the documentation for the modules that make up these packages
with more usage examples and the source code for all the classes and functions
implemented in sage-flatsurf. (Note that you can also access this documentation
from an interactive SageMath session with |help|_.

.. |help| replace:: ``?`` and ``??``
.. _help: https://ipython.readthedocs.io/en/stable/interactive/python-ipython-diff.html#accessing-help

.. toctree::
   :maxdepth: 1

   cache
   features
   geometry
   geometry/categories
   geometry/pyflatsurf
   graphical
