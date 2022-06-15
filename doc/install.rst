Installing sage-flatsurf
========================

Currently, sage-flatsurf is not widely available in package managers. We
therefore recommend to use mamba to install sage-flatsurf. This should work
nicely, if you are on Linux or macOS. If you are using Windows, you could try
to install sage-flatsurf with mamba in Windows Subsystem for Linux (WSL). If
you are having trouble with this or are on another operating system, please
`contact us <https://flatsurf.github.io>`_. We're thrilled about any new user of
sage-flatsurf and we're very happy to help and in the process improve these
installation instructions.

.. _installation-mamba:

Install with Mamba
------------------

Please first install `mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`_.

Here is some vocabulary, in case you are interested. Feel free to skip this
paragraph. mambaforge, is a program that installs both conda and mamba on your
system. `conda <https://conda.io>`_ is essentially a package manager for
Windows, macOS, and Linux that installs software without requiring superuser
permissions. `mamba <https://github.com/mamba-org/mamba>`_ is a fast version of
conda. mambaforge installs conda and mamba preconfigured to use packages from
conda-forge. `conda-forge <https://conda-forge.org>`_ is a community effort to
package lots of software and make it easily installable with conda and mamba.
The entire flatsurf stack, `in particular sage-flatsurf
<https://github.com/conda-forge/sage-flatsurf-feedstock/>`_, is available from
conda-forge. conda and mamba install packages in *conda environments*,
independent sets of software. You could for example have different environments
for different versions of sage-flatsurf and SageMath. To use software from an
environment, it needs to be *activated*. You can only activate a single
environment at a time.

Once you installed mambaforge, open a new command line terminal. The ``mamba``
command should now be available. The ``base`` environment is already activate,
it contains the mamba package manager. (If your terminal always prints
``(base)`` and you find that annoying, run ``conda config --set
auto_activate_base false`` to get rid of that prompt.)

You can now create a ``flatsurf`` environment (feel free to pick any other
name) with the entire flatsurf stack (`SageMath <https://sagemath.org>`_,
sage-flatsurf, `libflatsurf <https://github.com/flatsurf/flatsurf>`_,
`ipyvue-flatsurf <https://github.com/flatsurf/ipyvue-flatsurf>`_ widgets, ...)
installed::

        mamba env create -n flatsurf -f https://github.com/flatsurf/sage-flatsurf/raw/master/flatsurf.yml

Immediately after a new version of sage-flatsurf is released, this command
might fail because not all packages have been built for your platform yet. When
this is the case, you can also install a previous version of sage-flatsurf;
replace the ``0.5.0`` below with any recent `version of sage-flatsurf
<https://github.com/flatsurf/sage-flatsurf/releases>`_::

        mamba env create -n flatsurf -f https://raw.githubusercontent.com/flatsurf/sage-flatsurf/0.5.0/flatsurf.yml

You can now activate the flatsurf environment and start SageMath::

        conda activate flatsurf
        sage

If you prefer to work in a graphical environment, you can also launch Jupyter::

        conda activate flatsurf
        jupyter notebook

.. _upgrade-mamba:

Upgrade with Mamba
------------------

If you followed the above guide, it should be easy to upgrade to the latest
release of sage-flatsurf. Before you do that, it's usually a good idea to make
sure that you are using the latest conda and mamba::

        mamba upgrade -n base --all

To update sage-flatsurf and all the other packages in the flatsurf environment itself::

        mamba -n flatsurf upgrade --all

If after an upgrade something does not work correctly, please `let us know
<https://flatsurf.github.io>`_; we probably missed to declare some dependency.

You can also reinstall the flatsurf environment when this happens. First,
uninstall the flatsurf environment::

        mamba uninstall -n flatsurf --all

Then reinstall as above::

        mamba env create -n flatsurf -f https://github.com/flatsurf/sage-flatsurf/raw/master/flatsurf.yml

This does not delete or otherwise modify any of your notebooks or SageMath archives.

.. _installation-sagemath:

Install into SageMath
---------------------

If you are using a `source build of SageMath
<https://doc.sagemath.org/html/en/installation/source.html>`_ or if you
downloaded a SageMath binary, you can install sage-flatsurf into SageMath. Note
that this does not install all the optional dependencies of sage-flatsurf so
some computations might fail in this setup::

        sage -pip install sage-flatsurf

.. _installation-pip:

Install from PyPI
-----------------

You can install sage-flatsurf from `PyPI
<https://pypi.org/project/sage-flatsurf/>`_ if you installed sagelib as a
Python package. Again, this does not come with the optional dependencies of
sage-flatsurf, so some computations might fail in this setup::

        pip install --user sage-flatsurf
