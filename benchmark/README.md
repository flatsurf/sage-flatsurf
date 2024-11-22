# Airspeed Velocity Benchmarks

The benchmarks in this directory can be run with `asv`, e.g., `asv run
--dry-run --python=same --quick` to run all benchmarks once.

The benchmarks also run automatically for each Pull Request and timings can be
found in the "Benchmarks" section of the Pull Request checks.

Whenever a change is merged into our `master` branch, the benchmark results are
collected in [this public asv
page](https://flatsurf.github.io/sage-flatsurf/asv). Because these benchmarks
run on shared GitHub runners, there is a lot of noise in the output. However,
regressions should become clearly visible over time.
