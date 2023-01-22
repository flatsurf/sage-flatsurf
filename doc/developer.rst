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

We recommend to work on sage-flatsurf in a conda environment since currently
that is the easiest way to make sure that all optional dependencies of
sage-flatsurf are satisfied.

To setup a conda environment, it's maybe best to first install sage-flatsurf
with mamba, as explained in our `installation guide <installation-mamba>`_.

If that worked, you should create a separate environment to work on
sage-flatsurf::

        mamba env create -n sage-flatsurf-build -f https://github.com/flatsurf/sage-flatsurf/raw/master/environment.yml

This environment has all the dependencies of sage-flatsurf installed, but not
sage-flatsurf itself::

        mamba activate sage-flatsurf-build

`Clone <https://swcarpentry.github.io/git-novice/>`_ the sage-flatsurf
repository::

        git clone https://github.com/flatsurf/sage-flatsurf.git

You can now install an editable version of sage-flatsurf inside the ``sage-flatsurf-build`` environment, so any changes you make to ``sage-flatsurf/`` are going to be available immediately in this environment::

        pip install -e ./sage-flatsurf

You can now run our doctests, and run our test suite::

        sage -tp --initial --optional=sage,flipper,eantic,exactreal,pyflatsurf flatsurf doc
        pytest -n auto

Note that you can use ``mamba upgrade -n sage-flatsurf-build --all`` to update all of sage-flatsurf's dependencies. You can also recreate the environment to make sure that it's identical to the one that is used in our automated tests::

        conda deactivate
        mamba uninstall -n sage-flatsurf-build --all
        mamba env create -n sage-flatsurf-build -f ./sage-flatsurf/environment.yml
        mamba activate sage-flatsurf-build

This should cover the very basics of development but there are certainly lots
of things that we missed here, so don't hesitate to `contact us
<https://flatsurf.github.io>`_ if anything does not work out right away :)
