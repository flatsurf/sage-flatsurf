#!/bin/bash
###############################################################################
# This script downloads pixi if it's not present here yet by running the official pixi installer script.
###############################################################################
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PIXI="${SCRIPT_DIR}/.pixi-home/bin/pixi"

if [[ ! -x "$PIXI" ]]; then
  # We force a version of pixi that we used in our CI, see test.yml
  PIXI_VERSION="v0.39.5" PIXI_HOME="${SCRIPT_DIR}/.pixi-home" PIXI_NO_PATH_UPDATE="yes" "${SCRIPT_DIR}"/.pixi-install.sh
fi
