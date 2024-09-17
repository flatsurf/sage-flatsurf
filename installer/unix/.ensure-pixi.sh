#!/bin/bash
###############################################################################
# This script downloads pixi if it's not present in the current directory yet #
# essentially, we just run the official pixi script which we copied over into #
# this repository.                                                            #
###############################################################################
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PIXI="${SCRIPT_DIR}/.pixi-home/bin/pixi"

if [[ ! -x "$PIXI" ]]; then
  PIXI_HOME="${SCRIPT_DIR}/.pixi-home" PIXI_NO_PATH_UPDATE="yes" "${SCRIPT_DIR}"/.pixi-install.sh
fi
