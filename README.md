sage-flatsurf
=============

This is a software package (Sage module) for working with flat surfaces in 
[Sage](http://sagemath.org). 

We aim for this module to support the investigation of geometric, algebraic and 
dynamical questions related to flat surfaces. By flat surface we mean a surface
modeled on the plane with monodromy given by similarities of the plane, though
current efforts are focused on translation surfaces and half-translation 
surfaces.

Currently, the program can generate images of flat surfaces, plot straight-line
trajectories, deform surfaces through the SL(2,R) action, and compute Delaunay
decompositions. [Sage](http://sagemath.org) is used to perform exact arithmetic.

This software is open source (see the file LICENSE). We welcome any help to 
improve the package and especially to broaden the package's mathematical 
abilities.

The package is currently in active development. If you would like assistance
in using the package, please contact the authors.

There is also a related [flatsurf package](http://www.labri.fr/perso/vdelecro/programming.html).
The code in this repository is currently independent of the package but
the aim is to get them merged.

Installing the module
---------------------

Get the file in this git repository. Then run

    $ sage -pip install git+https://github.com/videlec/sage-flatsurf [--user] [--upgrade]

The options `--user` and `--upgrade` are optional. The option `--user` make
the installation in your home directory instead of the Sage sources. The
option `--upgrade` allows you to upgrade if the package is already installed.

Then you should be able to use the following within sage

    sage: import flatsurf.geometry.similarity_surface_generators as sfg
    sage: T = sfg.translation_surfaces.regular_octagon()
	sage: T
	Translation surface built from 1 polygon
	sage: T.stratum()
	H(2)

To uninstall the package, you can do `$ sage -pip uninstall flatsurf`.

Run the tests
-------------
:
    $ sage -t --force-lib ARG

where `ARG` is either a directory or file. In particular, to test all the
files in the module just do

    $ sage -t --force-lib flatsurf

Tests on the master branch are automatically run through [Travis-CI](https://travis-ci.org/videlec/sage-flatsurf?branch=master).

Primary Contributors
--------------------
* Vincent Delecroix (Bordeaux)
* Pat Hooper (City College of New York and CUNY Graduate Center)

We welcome others to contribute.

Acknowledgements
----------------
* This software project was created during a thematic semester at [ICERM](https://icerm.brown.edu).
* Hooper's contribution to the project has been supported by the National 
Science Foundation under Grant Number DMS-1500965. Any opinions, findings, 
and conclusions or recommendations expressed in this material are those of 
the authors and do not necessarily reflect the views of the National 
Science Foundation.

