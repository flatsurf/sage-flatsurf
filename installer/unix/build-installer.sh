#!/bin/bash
###################################################################
# 
###################################################################
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ $# -eq 0 ]; then
  NAME=`git describe --tags --abbrev=0`
  POST=`git rev-list --count $(git describe --tags --abbrev=0)..HEAD`
  if [ $POST != "0" ]; then
    NAME=$NAME.post$POST
  fi
elif [ $# -eq 1]; then
  NAME=$1
else
  echo "Usage: $0 [NAME]"
  exit 1
fi

TARBALL="$NAME.tar.gz"

echo "Generating $TARBALL in $SCRIPT_DIR"

( cd "$SCRIPT_DIR" && tar --transform "s?^?$NAME/?" -chzf "$TARBALL" sage-flatsurf LICENSE sage shell jupyterlab .ensure-pixi.sh .pixi-install.sh )
