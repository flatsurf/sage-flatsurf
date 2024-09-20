#!/bin/bash

VERSION=`git describe --tags --abbrev=0`
POST=`git rev-list --count $(git describe --tags --abbrev=0)..HEAD`
if [ $POST != "0" ]; then
  VERSION=$VERSION.post$POST
fi
