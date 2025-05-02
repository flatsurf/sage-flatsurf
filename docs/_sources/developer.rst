.. _developers-guide:

The Developer's Guide for sage-flatsurf
=======================================

sage-flatsurf is a pure Python module that depends on SageMath. We try to
follow the `programming conventions of SageMath
<https://doc.sagemath.org/html/en/developer/coding_basics.html>`_.

Contributions to sage-flatsurf are always welcome. If you want to contribute,
don't hesitate to `reach out to us <https://flatsurf.github.io>`_, create an
`issue <https://github.com/flatsurf/sage-flatsurf/issues>`_, or contribute a
`pull request <https://github.com/flatsurf/sage-flatsurf/pulls>`_.

We recommend to work on sage-flatsurf with a pixi environment which guarantees
that you are using dependencies that are known to work.

Once you intalled `pixi <https://pixi.sh>`_, you can enter a shell with your
version of sage-flatsurf installed by typing::

  pixi shell -e dev

You can then start SageMath normally, by typing ``sage``.

You can also run SageMath directly with::

  pixi run sage

To run the test suite, you can use ``sage -tp flatsurf`` and ``pytest``
directly or just run::

  pixi run -e dev doctest
  pixi run -e dev doctest-long  # to run the tests with the --long flag
  pixi run -e dev pytest

Or to just run all tests and doctests::

  pixi run test

To check your code for style errors::

  pixi run lint

And to preview the documentation::

  pixi run doc

This should cover the very basics of development but there are certainly lots
of things that we missed here, so don't hesitate to `contact us
<https://flatsurf.github.io>`_ if anything does not work out right away :)
