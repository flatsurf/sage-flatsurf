Installing sage-flatsurf
========================

There are several different ways to install sage-flatsurf on your machine. You can

* extract `our pixi tarball <#installation-tarball>`_ for Linux and macOS (recommended),
* use `our executable installer <#installation-installer>`_ on Windows (recommended),
* create a `conda environment <#installation-conda>`_ on Linux and macOS,
* install sage-flatsurf `into an existing source build of SageMath <#installation-sagemath>`_,
* or `pip install <#installation-pip>`_ sage-flatsurf.

If you are having trouble with this or are on another operating system, please
`contact us <https://flatsurf.github.io>`_. We're thrilled about any new user
of sage-flatsurf and we're very happy to help and in the process improve these
installation instructions.

.. _installation-tarball:

Install with the pixi tarball
-----------------------------

Open a terminal and run the following command::

  curl -fsSL https://github.com/flatsurf/sage-flatsurf/releases/download/0.7.3/sage-flatsurf-0.7.3.unix.tar.gz | tar zxf -

This will download the latest pixi tarball from our `Releases Page
<https://github.com/flatsurf/sage-flatsurf/releases/>`_ and extract it into a
subdirectory of the directory where you opened the terminal.

Due to conda limitations [1]_, this only works if the current directory does not
contain any spaces. Run ``pwd`` to double check that your absolute path is free
of whitespace.

The entire installation of sage-flatsurf and its dependencies are going to
reside in this subdirectory. This is not making any changes to your system. If
you later change your mind, you can safely delete that directory, or move or
rename it.

You can now use sage-flatsurf in a terminal or through Jupyter notebooks. To
launch sage-flatsurf in the terminal, run::

  ./sage-flatsurf-0.7.3/sage

To launch a browser with Jupyter Lab instead, run::

  ./sage-flatsurf-0.7.3/jupyterlab

The first time you run either of these commands, the installer downloads a copy
of SageMath and some other dependencies (in total this is going to use about
7GB of disk space). In particular on macOS this can take a very long time on
the first launch (probably due to security measures in the operating system),
please be patient.

When a new version of sage-flatsurf is released and you want to upgrade, just
download the latest tarball and extract it elsewhere. The two installations do
not interfere with each other and you can just delete the old version if you do
not need it anymore.

There is also a "nightly" build of sage-flatsurf that contains the latest
development version. To download it, go to `installer workflow runs
<https://github.com/flatsurf/sage-flatsurf/actions/workflows/installer.yml>`_
and click on the latest run, then from the Artifacts section download the
Unix or Windows installer. (This requires being signed in at GitHub.)

.. _installation-installer:

Install with the Windows Installer
----------------------------------

Download the latest .exe file from the `Releases Page
<https://github.com/flatsurf/sage-flatsurf/releases/>`_ and execute it
normally. The installer is going to configure your system to be able to run
`Windows Subsystem for Linux
<https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_. It is going to
ask for administrator privileges to install the necessary dependencies on your
system and a reboot might be required during the installation process.

Upon first launch of sage-flatsurf from your Start menu, a virtual machine
running Ubuntu is created in your ``%LOCALAPPDATA%``. This machine won't
interfere with any other WSL machines you might have set up on your system.

This Ubuntu system is then going to download all the dependencies of
sage-flatsurf (about 7GB once extracted). The first time, this might take some
time, please be patient.

To remove sage-flatsurf again, please use the provided "Uninstall Virtual
Machine" link from your Start menu, or just delete the sage-flatsurf directory
from ``%LOCALAPPDATA%``.

When a new version of sage-flatsurf is released and you want to upgrade, just
download the latest installer and run it. The two versions of sage-flatsurf do
not interfere with each other. If you do not need the old version of
sage-flatsurf anymore, just uninstall it.

.. _installation-conda:

Install with Conda
------------------

Almost the entire flatsurf stack is available at `conda-forge
<https://conda-forge.org>`_.

If you already have `miniforge
<https://github.com/conda-forge/miniforge>`_ installed, you can
create an environment with the entire flatsurf stack by running::

  conda create -n flatsurf sage-flatsurf pyflatsurf pyexactreal sage pip

Some optional bits of the flatsurf stack are only available on PyPI, to install
them as well run::

  conda activate flatsurf
  pip install ipyvue-flatsurf flipper realalg veerer

.. _installation-sagemath:

Install into SageMath
---------------------

If you are using a `source build of SageMath
<https://doc.sagemath.org/html/en/installation/source.html>`_ or if you
downloaded a SageMath binary, you can install sage-flatsurf into SageMath. Note
that this does not install all the optional dependencies of sage-flatsurf so
some computations might fail in this setup::

        sage -pip install sage-flatsurf

To uninstall sage-flatsurf again later::

        sage -pip uninstall sage-flatsurf

.. _installation-pip:

Install from PyPI
-----------------

You can install sage-flatsurf from `PyPI
<https://pypi.org/project/sage-flatsurf/>`_ if you installed sagelib as a
Python package. Again, this does not come with the optional dependencies of
sage-flatsurf, so some computations might fail in this setup::

        pip install --user sage-flatsurf

To uninstall sage-flatsurf again later::

        pip uninstall --user sage-flatsurf


.. [1] conda (which pixi uses internally) rewrites the shebangs in executable
   scripts. Shebangs must not contain spaces and conda presently has no way to
   work around this limitation.
