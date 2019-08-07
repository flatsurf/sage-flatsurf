sage-flatsurf
=============

This is a software package for working with flat surfaces in
`SageMath`_. The flatsurf documentation
is available at http://www.labri.fr/perso/vdelecro/flatsurf/latest/

We aim for this module to support the investigation of geometric, algebraic and
dynamical questions related to flat surfaces. By flat surface we mean a surface
modeled on the plane with monodromy given by similarities of the plane, though
current efforts are focused on translation surfaces and half-translation
surfaces.

Currently, the program can generate images of flat surfaces, plot straight-line
trajectories, deform surfaces through the SL(2,R) action, and compute Delaunay
decompositions. `SageMath`_ is used to perform exact arithmetic.

This software is open source released under GPL v2 (see the COPYING file). We
welcome any help to improve the package and especially to broaden the package's
mathematical abilities.

The package is currently in active development. If you would like assistance
in using the package, please contact the authors.

There is also a related `surface-dynamics`_.
The code in this repository is currently independent of the package but
the aim is to get them merged.

Installing the dependency
-------------------------

Our software depends on the `surface-dynamics`_.
The module is distributed on PyPI. To install it, you just need to run the
following command::

    $ sage -pip install surface_dynamics [--user]

The --user option is optional and allows to install the module in your user
space (and does not require administrator rights).

Installing the module
---------------------

sage-flatsurf is available on PyPI at https://pypi.org/project/sage-flatsurf/. To install the
released version of sage-flatsurf, run the following command::

    $ sage -pip install sage-flatsurf [--user] [--upgrade]

To install the development version of sage-flatsurf, run::

    $ sage -pip install git+https://github.com/videlec/sage-flatsurf [--user] [--upgrade]

The options `--user` and `--upgrade` are optional. The option `--user` make
the installation in your home directory instead of the Sage sources. The
option `--upgrade` allows you to upgrade if the package is already installed.

This might fail because [git](https://git-scm.com/) is not installed on your computer
(this is for example if you run Sage in Windows). In this situation you have two options.
Either you install git. Or you download this project from the "Clone or download" drop
menu above (you should get a zip file). Then you need to run the command::

    $ sage -pip install TARBALL_NAME [--user] [--upgrade]

where `TARBALL_NAME` has to be replaced with the full path to your tarball. If you
run windows, it should be a cygwin path and will looks something like
`/cygdrive/c/Users/You/Downloads/sage-flatsurf-master.zip`.

Then you should be able to use the following within sage::

    sage: import flatsurf.geometry.similarity_surface_generators as sfg
    sage: T = sfg.translation_surfaces.regular_octagon()
    sage: T
    Translation surface built from 1 polygon
    sage: T.stratum()
    H_2(2)

To uninstall the package, you can do `$ sage -pip uninstall flatsurf`.

Run the tests
-------------
::
    $ sage -t --force-lib ARG

where `ARG` is either a directory or file. In particular, to test all the
files in the module just do::

    $ sage -t --force-lib flatsurf

Tests on the master branch are automatically run through `Travis-CI <https://travis-ci.org/videlec/sage-flatsurf?branch=master>`_.

Primary Contributors
--------------------

* Vincent Delecroix (Bordeaux)
* Pat Hooper (City College of New York and CUNY Graduate Center)

We welcome others to contribute.

Acknowledgements
----------------

* This software project was created during a thematic semester at `ICERM <https://icerm.brown.edu>`_.
* Hooper's contribution to the project has been supported by the National
  Science Foundation under Grant Number DMS-1500965. Any opinions, findings,
  and conclusions or recommendations expressed in this material are those of
  the authors and do not necessarily reflect the views of the National
  Science Foundation.
* Delecroix' contribution to the project has been supported by OpenDreamKit,
  Horizon 2020 European Research Infrastructures project #676541.

.. _SageMath: http://sagemath.org
.. _surface-dynamics:
