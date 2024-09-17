This installer performs the following steps:

* Enable WSL2 (and perform a reboot if needed.)
* Install ArchLinux into WSL2 with https://github.com/yuk7/wsldl
* Copy the ../unix installer into the users home directory in ArchLinux
* Create shortcuts to invoke the launcher scripts of the unix installer
* or, uninstall which just deletes the directory where ArchLinux lives.
* Mirror the rootfs so we do not rely on that random tarball.
