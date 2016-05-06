sage-flatsurf
=============

This is a github repository that was started during a thematic semester
at [ICERM](https://icerm.brown.edu/home/index.php) by Vincent Delecroix
and Pat Hooper. It contains some code to interactively draw translation
surfaces from [Sage](http://sagemath.org).

It currently contains a lot of bugs... if you intend to use it, please
contact the authors.

There is also a related [flatsurf package](http://www.labri.fr/perso/vdelecro/programming.html).
The code in this repository is currently independent of the package but
the aim is to get them merged.

Installing the module
---------------------

Get the file in this git repository. Then run

    $ sage -python setup.py install

Then you should be able to use the following within sage

    sage: import flatsurf.geometry.similarity_surface_generators as sfg
    sage: T = sfg.translation_surfaces.regular_octagon()
	sage: T
	Translation surface built from 1 polygon
	sage: T.stratum()
	H(2)

Run the tests
-------------

    $ sage -t --force-lib ARG

where `ARG` is either a directory or file. In particular, to test all the
files in the module just do

    $ sage -t --force-lib flatsurf
