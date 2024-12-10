# sage-flatsurf installer for Windows

This [Inno Setup](https://jrsoftware.org/isinfo.php) installer performs the
following steps:

* Enable WSL2 (and perform a reboot if needed).
* Copy the ../unix installer to Program Files.
* Create shortcuts for the entrypoints of the ../unix installer.
* Update the WSL kernel.

The launcher shortcuts are then going to perform the following steps:

* Create a virtual machine running Ubuntu for this sage-flatsurf version with
  Ubuntu (unless already present).
* Extract the ../unix installer into this Ubuntu.
* Launch the entrypoint.

Note that this approach is a bit odd. Normally, programs would install everything into Program Files and only put a little configuration into user space.
