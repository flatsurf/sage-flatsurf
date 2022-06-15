sage-flatsurf: Flat Surfaces in SageMath
========================================

sage-flatsurf is a Python package for working with flat surfaces in `SageMath
<https://sagemath.org>`.

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
(e.g. `apt-get` on Debian or Ubuntu, `pacman` on Arch Linux, `brew` on MacOS,
etc). However, as of this writing sage-flatsurf has not been picked up by `any
of the major distributions yet <https://repology.org/project/sage-flatsurf/packages>`.

We therefore recommend to install sage-flatsurf with the `Mamba
<https://github.com/mamba-org/mamba>` package manager.

Detailed installation instructions:

* :ref:`Install with Mamba <installation-mamba>` and :ref:`Upgrade an existing
  installation with Mamba <upgrade-mamba>`
* :ref:`Install into an existing SageMath source build <installation-sagemath>`
* :ref:`Install with pip <installation-pip>`

If you are planning to work on the sage-flatsurf source code, you can also
build sage-flatsurf from source. For this, please have a look at our
:ref:`Developer's Guide <developers-guide>`.

A Tour of sage-flatsurf
=======================

The following presents some of the capabilities of sage-flatsurf.

.. toctree::
   :maxdepth: 1

   examples/defining_surfaces
   examples/linear_action_and_delaunay
   examples/rel_deformations
   examples/saddle_connections
   examples/apisa_wright
   examples/siegel_veech
   examples/straight_line_flow
   examples/warwick-2017

These examples can also be explored interactively by clicking |binder|_. The interactive session might take a moment to start. Once ready, press Shift + Enter to execute the cells, starting from the top.

.. |binder| image:: https://mybinder.org/badge_logo.svg
.. _binder: https://mybinder.org/v2/gh/flatsurf/sage-flatsurf/0.4.6?filepath=doc%2Fexamples

Reference Manual
================

.. toctree::
   :maxdepth: 2

   geometry
   graphical
