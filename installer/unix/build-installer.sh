#!/bin/bash
###########################################################################
# Creates a tarball containing a pixi environment and launcher scripts for
# sage-flatsurf.
###########################################################################
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ $# -eq 0 ]; then
  source "$SCRIPT_DIR/version.sh"
  NAME=sage-flatsurf-$VERSION
elif [ $# -eq 1]; then
  NAME=$1
else
  echo "Usage: $0 [NAME]"
  exit 1
fi

TARBALL="$NAME.unix.tar.gz"

echo "Generating $TARBALL in $SCRIPT_DIR"

tmp=$(mktemp -d)

trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/$NAME"

git clone ${SCRIPT_DIR}/../../ "$tmp/$NAME/sage-flatsurf"

rm -rf "$tmp/$NAME/sage-flatsurf/.git"

( cd "$SCRIPT_DIR" && cp LICENSE sage shell jupyterlab .ensure-pixi.sh "$tmp/$NAME" )
curl -fsSL https://pixi.sh/install.sh > "$tmp/$NAME/.pixi-install.sh"
chmod +x "$tmp/$NAME/.pixi-install.sh"

( cd "$tmp" && tar czf "$TARBALL" "$NAME")

mv "$tmp/$TARBALL" .
