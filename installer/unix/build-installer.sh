#!/bin/bash
###################################################################
# 
###################################################################
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ $# -eq 0 ]; then
  NAME=sage-flatsurf-`git describe --tags --abbrev=0`
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

TARBALL="$NAME.pixi.tar.gz"

echo "Generating $TARBALL in $SCRIPT_DIR"

tmp=$(mktemp -d)

trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/$NAME"

git clone ../../ "$tmp/$NAME/sage-flatsurf"

rm -rf "$tmp/$NAME/sage-flatsurf/.git"

cp LICENSE sage shell jupyterlab .ensure-pixi.sh "$tmp/$NAME"
curl -fsSL https://pixi.sh/install.sh > "$tmp/$NAME/.pixi-install.sh"

( cd "$tmp" && tar czf "$TARBALL" "$NAME")

mv "$tmp/$TARBALL" .
