# sage-flatsurf installer for Linux & macOS

Create the installer by running `./build-installer.sh sage-flatsurf-0.5.2`.

This packages a directory with:

* a copy of flatsurf
* launcher scripts: ./sage ./shell ./jupyterlab

The launchers download pixi if it's not present in the directory and then
launch things through pixi.

## Installation

```sh
tar zxf sage-flatsurf-*.unix.tar.gz
cd sage-flatsurf-*
./sage or ./jupyterlab
```

## Limitations

This only works for the platforms for which sage-flatsurf and its dependencies
are available, i.e., at the time of this writing, Linux and macOS on x86\_64.

## Graphical Installer for Linux

Linux users do not expect programs to have graphical installers. Rather they
expect things to be present in their distribution directly. Eventually we
should bring sage-flatsurf downstream to the distributions.

For a graphical macOS installer, see the README in ../macos.
