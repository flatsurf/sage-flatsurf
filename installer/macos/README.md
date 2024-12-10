# Graphical macOS Installer

Currently, there is no graphical macOS installer. The generic unix installer
works though.

From some experimentation, it seems that it would be relatively easy to build a
graphical installer. The easiest would be to create one `.app` directory per
entrypoint (sage, shell, jupyterlab) and embed the Unix installer tarball in
each of them as a `Resources/`. In `macOS`, we can then put a launcher script
that unpacks the tarball to `$HOME/sage-flatsurf-whatever` and runs the
launcher script from the Unix installer with `open -a Console`. Finally, these
three apps can be trivially packaged up into a .pkg with pkgbuild. To make this
useful, they need to be signed and notarized (not sure if Apple will notarize
something like that).

This is quite untypical for a macOS app (such as the 3-manifolds version of
SageMath) which usually ships all the binaries and installs everything into
/Applications (and gets everything notarized). However, installing into
user-space (similarly to how we do this in the Windows installer) gives us a
bit more flexibility, e.g., when installing further dependencies.

One obvious improvement to the above would be to have proper status indicators
in the macOS bar, i.e., by not just opening a terminal with the launcher but
wrapping it a bit more nicely so we can keep track of the lifetime of things
better.
