# sage-flatsurf installer for Linux & macOS

Create the installer by running `./build-installer.sh sage-flatsurf-0.5.2`.

This packages a directory with:

* a copy of flatsurf
* launcher scripts: ./sage ./shell ./jupyterlab

The launchers download pixi if it's not present in the directory and then
launches things through pixi.

## Installation

```sh
tar zxf sage-flatsurf-0.5.2.tar.gz
cd sage-flatsurf-0.5.2
./sage or ./jupyterlab
```

## Limitations

This only works for the platforms for which sage-flatsurf and its dependencies
are available, i.e., at the time of this writing, Linux and macOS on x86\_64.
