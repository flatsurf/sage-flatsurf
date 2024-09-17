#!/bin/bash
###################################################################
# 
###################################################################
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ $# -ne 1 ]; then
  echo "Usage: $0 name"
  exit 1
fi

TARBALL="$1.tar.gz"

echo "Generating $TARBALL in $SCRIPT_DIR"

( cd "$SCRIPT_DIR" && tar --transform "s?^?$1/?" -czf "$TARBALL" pixi.lock pixi.toml sage shell jupyterlab .ensure-pixi.sh .pixi-install.sh )
