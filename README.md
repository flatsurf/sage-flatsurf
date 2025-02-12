<p align="center">
    <img alt="logo" src="https://github.com/flatsurf/sage-flatsurf/raw/master/doc/static/logo.svg?sanitize=true">
</p>

<h1><p align="center">sage-flatsurf</p></h1>

sage-flatsurf is a Python package for working with flat surfaces in
[SageMath](https://sagemath.org).

We aim for sage-flatsurf to support the investigation of geometric, algebraic
and dynamical questions related to flat surfaces. By flat surface we mean a
surface modeled on the plane with monodromy given by similarities of the plane,
though current efforts are focused on [translation
surfaces](https://en.wikipedia.org/wiki/Translation_surface) and
[half-translation
surfaces](https://en.wikipedia.org/wiki/Translation_surface#Half-translation_surfaces).

Take the [Tour of flatsurf](https://flatsurf.github.io/sage-flatsurf/examples/tour)
to see some of the capabilities of sage-flatsurf.

sage-flatsurf is free software, released under the [GPL v2 (or later)](./COPYING).

We welcome any help to improve sage-flatsurf. If you would like to help, have
ideas for improvements, or if you need any assistance in using sage-flatsurf,
please don't hesitate to [contact us](https://flatsurf.github.io#contact).

## Installation

If you are on **Linux or macOS**, download the latest `.unix.tar.gz` file from our
[Releases page](https://github.com/flatsurf/sage-flatsurf/releases).

Extract it anywhere (make sure there are no spaces in the directory name) and
run `./sage` or `./jupyterlab`.

```sh
tar zxf sage-flatsurf-0.7.3.unix.tar.gz
./sage-flatsurf-0.7.3/jupyterlab  # or
./sage-flatsurf-0.7.3/sage
```

If you are on **Windows**, download the latest `.exe` installer from our [Releases
page](https://github.com/flatsurf/sage-flatsurf/releases).

Please also consult [our
documentation](https://flatsurf.github.io/sage-flatsurf/#installation) for
other options and more detailed instructions.

## Developing sage-flatsurf

We recommend you install [pixi](https://pixi.sh) to provide all the
dependencies for sage-flatsurf. Once installed, `git clone` this repository and
then

```sh
pixi run sage  # to run SageMath with your version of sage-flatsurf installed
pixi run test  # to run the test suite
pixi run lint  # to check for errors and formatting issues
```

Please consult our [Developer's
Guide](https://flatsurf.github.io/sage-flatsurf/developer.html) for more
details.

## Contributors

The main authors and current maintainers of sage-flatsurf are:

* Vincent Delecroix (Bordeaux)
* W. Patrick Hooper (City College of New York and CUNY Graduate Center)
* Julian Rüth

We welcome others to [contribute](https://flatsurf.github.io#contact).

## How to Cite This Project

If you have used this project, please cite us as described [on our
zenodo website](https://zenodo.org/badge/latestdoi/13970050).

## Acknowledgements

* sage-flatsurf was started during a thematic semester at
  [ICERM](https://icerm.brown.edu).
* Vincent Delecroix's contribution to the project has been supported by
  OpenDreamKit, Horizon 2020 European Research Infrastructures project #676541.
* W. Patrick Hooper's contribution to the project has been supported by the National
  Science Foundation under Grant Number DMS-1500965. Any opinions, findings,
  and conclusions or recommendations expressed in this material are those of
  the authors and do not necessarily reflect the views of the National Science
  Foundation.
* Julian Rüth's contributions to this project have been supported by the Simons
  Foundation Investigator grant of Alex Eskin.
